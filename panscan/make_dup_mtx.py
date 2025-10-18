import os
import glob
import subprocess
import pandas as pd
import argparse
import shlex
from io import StringIO

def ensure_symlink(target_path, link_path):
    target_abs = os.path.abspath(target_path)
    # If link exists but is wrong, or a regular file/dir exists -> remove it
    if os.path.lexists(link_path):
        try:
            if os.path.islink(link_path):
                if os.readlink(link_path) == target_abs:
                    return  # correct link already there
            # Anything else: remove and recreate
            os.unlink(link_path)
        except OSError:
            # If it's a directory, etc.
            if os.path.isdir(link_path):
                raise RuntimeError(f"Refusing to replace directory at {link_path}")
            os.remove(link_path)
    os.symlink(target_abs, link_path)

def run_liftoff(assembly_fa, reference_fa, gencode_gff3, out_prefix, threads, log_path):
    # Build argv list to avoid shell quoting issues
    cmd = [
        "liftoff",
        "-p", str(threads),
        "-sc", "0.90",
        "-copies",
        "-g", gencode_gff3,
        "-u", f"{out_prefix}.unmapped.txt",
        "-o", f"{out_prefix}.genedup.gff3",
        "-polish", assembly_fa, reference_fa
    ]
    with open(log_path, "w") as logf:
        subprocess.run(cmd, check=True, stdout=logf, stderr=subprocess.STDOUT)

def process_assemblies(gencode_gff3, fofn_path, ref_fa, threads):
    os.makedirs('copy-num', exist_ok=True)

    assembly_info = []
    with open(fofn_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                assembly_info.append(parts)

    for assembly_path, sample_name, haplotype in assembly_info:
        assembly_dir_name = f"{sample_name}_{haplotype}"
        assembly_dir = os.path.join('copy-num', assembly_dir_name)
        os.makedirs(assembly_dir, exist_ok=True)

        gencode_symlink = os.path.join(assembly_dir, 'gencode.gff3')
        ensure_symlink(gencode_gff3, gencode_symlink)

        out_prefix = os.path.join(assembly_dir, assembly_dir_name)
        log_path = os.path.join(assembly_dir, 'liftoff.log')
        run_liftoff(assembly_path, ref_fa, gencode_symlink, out_prefix, threads, log_path)

    # Reference self-run
    ref_dir = os.path.join('copy-num', 'reference')
    os.makedirs(ref_dir, exist_ok=True)
    gencode_symlink = os.path.join(ref_dir, 'gencode.gff3')
    ensure_symlink(gencode_gff3, gencode_symlink)
    out_prefix = os.path.join(ref_dir, 'reference')
    log_path = os.path.join(ref_dir, 'liftoff.log')
    run_liftoff(ref_fa, ref_fa, gencode_symlink, out_prefix, threads, log_path)

    # --- Merge & matrix ---
    gff3_files = glob.glob('copy-num/**/*.genedup.gff3', recursive=True)
    gff_data_list = []
    for file_path in gff3_files:
        sample_name = os.path.basename(file_path).rsplit('.', 2)[0]
        with open(file_path, 'r') as f:
            lines = [line for line in f if not line.startswith('#')]
        gff_data = pd.read_csv(StringIO('\n'.join(lines)), sep='\t', header=None, dtype=str)
        gff_data = gff_data[gff_data[2] == 'gene']
        gff_data['sample'] = sample_name
        gff_data_list.append(gff_data)

    combined = pd.concat(gff_data_list, ignore_index=True)
    attrs = combined[8]
    combined['cn'] = attrs.str.extract(r'extra_copy_number=([^;]+)').astype(float)
    combined['gene_id'] = attrs.str.extract(r'gene_id=([^;]+)').iloc[:,0].str.split('_').str[0]
    combined['gene_name'] = attrs.str.extract(r'gene_name=([^;]+)')
    combined['gene_type'] = attrs.str.extract(r'gene_type=([^;]+)')
    combined.to_csv('all_dup_data_gene.csv', index=False)

    filtered = combined[combined['cn'] >= 1]
    matrix = filtered.pivot_table(index='gene_name', columns='sample', values='cn', aggfunc='max').fillna(0)
    matrix.to_csv('gene-dup-matrix-before-ref-correction.csv')

    # Optional reference subtraction
    if 'reference' in matrix.columns:
        matrix = matrix.sub(matrix['reference'], axis=0).drop(columns=['reference']).clip(lower=0)

    # Now always prune all-zero rows
    matrix = matrix.loc[~(matrix == 0).all(axis=1)]
    matrix.to_csv('gene-dup-matrix-final.csv')

    # Ideogram input
    max_cn = matrix.max(axis=1) * 1000
    gene_coords = {}
    with open(gencode_gff3, 'r') as gf:
        for line in gf:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if cols[2] != 'gene':
                continue
            seqid, start, end = cols[0], int(cols[3]), int(cols[4])
            # safer split
            attr_pairs = [x.split('=', 1) for x in cols[8].split(';') if '=' in x]
            attr = {k: v for k, v in attr_pairs}
            gname = attr.get('gene_name')
            if gname:
                if not seqid.startswith('chr'):
                    seqid = 'chr' + seqid
                gene_coords[gname] = (seqid, start, end)

    rows = []
    for gene, val in max_cn.items():
        coords = gene_coords.get(gene)
        if coords:
            seqid, st, en = coords
            rows.append({'Chr': seqid, 'Start': st, 'End': en, 'Value': int(val)})

    pd.DataFrame(rows, columns=['Chr','Start','End','Value']).to_csv('ideo-ip.csv', index=False)

def add_subparser(subparsers):
    parser = subparsers.add_parser(
        "make_dup_mtx",
        help="Generate a duplication matrix based on gene duplication data."
    )
    parser.add_argument("--gencode_gff3", required=True, help="Path to the GENCODE GFF3 file.")
    parser.add_argument("--fofn", required=True, help="Path to the FOFN (assemblies list).")
    parser.add_argument("--ref_fa", required=True, help="Path to the reference FASTA.")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads.")
    parser.set_defaults(func=run_make_dup_mtx)

def run_make_dup_mtx(args):
    """Run the make-dup-mtx process with the provided arguments."""
    process_assemblies(args.gencode_gff3, args.fofn, args.ref_fa, args.threads)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a duplication matrix based on gene duplication data.")
    parser.add_argument("--gencode_gff3", required=True)
    parser.add_argument("--fofn", required=True)
    parser.add_argument("--ref_fa", required=True)
    parser.add_argument("--threads", type=int, required=True)
    args = parser.parse_args()
    process_assemblies(args.gencode_gff3, args.fofn, args.ref_fa, args.threads)
