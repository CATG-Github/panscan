import argparse
import subprocess
import os

def run_find_uniq_variants(vcf, var_type, pbsv, db, overlap, op, db_path, debug):
    # Convert VCF path to an absolute path.
    vcf = os.path.abspath(vcf)
    if debug:
        print("Absolute VCF path:", vcf)
    
    # Determine the directory of this Python file (package installation directory)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if debug:
        print("Script directory:", script_dir)
    
    # Construct the absolute path to the Perl script (assumed in the "scripts" subdirectory)
    perl_script = os.path.join(script_dir, "scripts", "findUniqVariants.pl")
    perl_script = os.path.abspath(perl_script)
    if debug:
        print("Perl script path:", perl_script)
    if not os.path.isfile(perl_script):
        raise FileNotFoundError(f"Perl script not found: {perl_script}")
    
    # Determine output directory if provided.
    if op:
        op = os.path.abspath(op)
    if debug:
        print("Output directory:", op if op else "(default will be used by Perl script)")
    
    # Build the command list.
    # The Perl script accepts:
    #   -i <input VCF>
    #   -t <variant type>
    #   -pbsv flag (if enabled)
    #   -db <database or ALL>
    #   -overlap <overlap percentage>
    #   -op <output directory> (if provided)
    #   -db-path <database base path> (if provided)
    cmd = [
        "perl", perl_script,
        "-i", vcf,
        "-t", var_type,
        "-db", db,
        "-overlap", str(overlap)
    ]
    if pbsv:
        cmd.append("-pbsv")
    if op:
        cmd.extend(["-op", op])
    if db_path:
        # Pass the provided db-path flag.
        cmd.extend(["-db-path", db_path])
    
    if debug:
        print("Command to execute:", " ".join(cmd))
    
    # Set the working directory to the Perl script's "scripts" folder so that config.yaml and modules are found.
    cwd = os.path.join(script_dir, "scripts")
    env = os.environ.copy()
    # Add the scripts directory to PERL5LIB so that Perl locates its modules.
    env["PERL5LIB"] = cwd + ":" + env.get("PERL5LIB", "")
    if debug:
        print("Setting PERL5LIB to:", env["PERL5LIB"])
    
    subprocess.run(cmd, check=True, cwd=cwd, env=env)

def main(args):
    run_find_uniq_variants(
        args.vcf, args.var_type, args.pbsv, args.db, args.overlap,
        args.op, args.db_path, args.debug
    )

def add_subparser(subparsers):
    parser = subparsers.add_parser("find_uniq_variants",
                                   help="Find unique variants in a VCF file by comparing with databases.")
    parser.add_argument("-i", "--vcf", required=True, help="Path to the input VCF file.")
    parser.add_argument("-t", "--var_type", required=True, choices=["SNP", "INDEL", "SV"],
                        help="Specify the variant type to compare (SNP/INDEL/SV).")
    parser.add_argument("--pbsv", action="store_true",
                        help="Enable pbsv mode for SV analysis.")
    parser.add_argument("-db", "--db", default="ALL",
                        help="""Specify the databases (comma-separated).
For SNP: dbSNP, gnomAD, 1000Genomes, GME.
For INDEL: GNOMAD_INDEL, 1000Genome_INDEL, GME_INDEL.
For SV: 1000Genome_SV_DEL, DGV_SV_DEL.
Default: ALL""")
    parser.add_argument("-overlap", "--overlap", type=int, default=80,
                        help="Specify the percentage of overlap for SV comparison (default: 80).")
    parser.add_argument("--op", help="Output directory")
    parser.add_argument("--db-path", help="Optional base directory for databases. Overrides the paths in config.yaml if provided.")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode to print diagnostic messages.")
    parser.set_defaults(func=main)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find unique variants in a VCF file using the updated Perl script."
    )
    subparsers = parser.add_subparsers()
    add_subparser(subparsers)
    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


