import os
import glob
import subprocess
import pandas as pd
import sys

def process_assemblies(gencode_gff3, fofn_path, ref_fa, threads):
    # Step 1: Set up the main folder and directories
    os.makedirs('copy-num', exist_ok=True)

    # Read assembly paths from the file of files (fofn_path)
    assembly_info = []
    with open(fofn_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 3:
                assembly_info.append(parts)

    # Step 2: Process each assembly sequentially
    for assembly_path, sample_name, haplotype in assembly_info:
        assembly_dir_name = f"{sample_name}_{haplotype}"
        assembly_dir = os.path.join('copy-num', assembly_dir_name)
        os.makedirs(assembly_dir, exist_ok=True)
        gencode_symlink = os.path.join(assembly_dir, 'gencode.gff3')
        if not os.path.islink(gencode_symlink):
            os.symlink(gencode_gff3, gencode_symlink)
        liftoff_cmd = (
            f"liftoff -p {threads} -sc 0.90 -copies -g {gencode_symlink} "
            f"-u {os.path.join(assembly_dir, assembly_dir_name)}.unmapped.txt "
            f"-o {os.path.join(assembly_dir, assembly_dir_name)}.genedup.gff3 "
            f"-polish {assembly_path} {ref_fa} "
            f"> {os.path.join(assembly_dir, 'liftoff.log')} 2>&1"
        )
        subprocess.run(liftoff_cmd, shell=True, check=True)

    # Also process the reference genome
    ref_dir = os.path.join('copy-num', 'reference')
    os.makedirs(ref_dir, exist_ok=True)
    gencode_symlink = os.path.join(ref_dir, 'gencode.gff3')
    if not os.path.islink(gencode_symlink):
        os.symlink(gencode_gff3, gencode_symlink)
    liftoff_cmd = (
        f"liftoff -p {threads} -sc 0.90 -copies -g {gencode_symlink} "
        f"-u {os.path.join(ref_dir, 'reference')}.unmapped.txt "
        f"-o {os.path.join(ref_dir, 'reference')}.genedup.gff3 "
        f"-polish {ref_fa} {ref_fa} "
        f"> {os.path.join(ref_dir, 'liftoff.log')} 2>&1"
    )
    subprocess.run(liftoff_cmd, shell=True, check=True)

    # Modified file reading and merging section
    gff3_files = glob.glob('copy-num/**/*.genedup.gff3', recursive=True)
    gff_data_list = []
    
    for file_path in gff3_files:
        sample_name = os.path.basename(file_path).rsplit('.', 2)[0]
        # Read only lines that don't start with #
        with open(file_path, 'r') as f:
            lines = [line for line in f if not line.startswith('#')]

        # Convert lines to DataFrame
        gff_data = pd.read_csv(
            pd.io.common.StringIO('\n'.join(lines)),
            sep='\t',
            header=None,
            dtype=str  # Read all columns as strings initially
        )
        # Filter for gene entries
        gff_data = gff_data[gff_data[2] == 'gene']
        gff_data['sample'] = sample_name
        gff_data_list.append(gff_data)
        print(f"Processed {file_path}")

    # Combine all the DataFrames
    combined_gff_data = pd.concat(gff_data_list, ignore_index=True)
    
    
    # Rest of your original code
    combined_gff_data['cn'] = combined_gff_data[8].str.extract(r'extra_copy_number=([^;]+)').astype(float)
    combined_gff_data['gene_id'] = combined_gff_data[8].str.extract(r'gene_id=([^;]+)')
    combined_gff_data['gene_id'] = combined_gff_data['gene_id'].str.split('_').str[0]
    combined_gff_data['gene_name'] = combined_gff_data[8].str.extract(r'gene_name=([^;]+)')
    combined_gff_data['gene_type'] = combined_gff_data[8].str.extract(r'gene_type=([^;]+)')
    combined_gff_data.to_csv('all_dup_data_gene.csv', index=False)
    
    filtered_data = combined_gff_data[combined_gff_data['cn'] >= 1]
    matrix = filtered_data.pivot_table(index='gene_name', columns='sample', values='cn', aggfunc='max').fillna(0)
    
    matrix.to_csv('gene-dup-matrix-before-ref-correction.csv')
    # Subtract reference column from all others
    if 'reference' in matrix.columns:
        matrix = matrix.sub(matrix['reference'], axis=0)
        matrix = matrix.drop(columns=['reference'])
        matrix = matrix.clip(lower=0)
        matrix = matrix.loc[~(matrix == 0).all(axis=1)]
    matrix.to_csv('gene-dup-matrix-final.csv')


    max_cn = matrix.max(axis=1) * 1000

    # 2) Parse the GENCODE GFF3 for gene coordinates
    gene_coords = {}
    with open(gencode_gff3, 'r') as gf:
        for line in gf:
            if line.startswith('#'):
                continue
            cols = line.rstrip('\n').split('\t')
            if cols[2] != 'gene':
                continue
            seqid, start, end = cols[0], int(cols[3]), int(cols[4])
            # pull gene_name=XYZ from the 9th column
            attr = dict(item.split('=') for item in cols[8].split(';') if '=' in item)
            gname = attr.get('gene_name')
            if gname:
                # normalize seqid to always start with 'chr'
                if not seqid.startswith('chr'):
                    seqid = 'chr' + seqid
                gene_coords[gname] = (seqid, start, end)

    rows = []
    for gene, val in max_cn.items():
        coords = gene_coords.get(gene)
        if not coords:
            # gene not found in GFF3, skip
            continue
        seqid, st, en = coords
        rows.append({
            'Chr': seqid,
            'Start': st,
            'End': en,
            'Value': int(val)
        })

    # Write out the ideo‐ip.csv
    ideo_df = pd.DataFrame(rows, columns=['Chr','Start','End','Value'])
    ideo_df.to_csv('ideo-ip.csv', index=False)
    print(f"Written {len(ideo_df)} entries to ideo-ip.csv")


def add_subparser(subparsers):
    """Add subparser for the make_dup_mtx command."""
    parser = subparsers.add_parser(
        "make_dup_mtx",
        help="Generate a duplication matrix based on gene duplication data."
    )
    parser.add_argument("--gencode_gff3", required=True, help="Path to the GENCODE GFF3 file.")
    parser.add_argument("--fofn", required=True, help="Path to the file of filenames (FOFN).")
    parser.add_argument("--ref_fa", required=True, help="Path to the reference FASTA file.")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads to use.")
    parser.set_defaults(func=run_make_dup_mtx)

def run_make_dup_mtx(args):
    """Run the make-dup-mtx process with the provided arguments."""
    process_assemblies(args.gencode_gff3, args.fofn, args.ref_fa, args.threads)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate a duplication matrix based on gene duplication data.")
    parser.add_argument("--gencode_gff3", required=True, help="Path to the GENCODE GFF3 file.")
    parser.add_argument("--fofn", required=True, help="Path to the file of filenames (FOFN).")
    parser.add_argument("--ref_fa", required=True, help="Path to the reference FASTA file.")
    parser.add_argument("--threads", type=int, required=True, help="Number of threads to use.")
    
    args = parser.parse_args()
    process_assemblies(args.gencode_gff3, args.fofn_path, args.ref_fa, args.threads)

