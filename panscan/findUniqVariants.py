import argparse
import subprocess
import os

def run_find_uniq_variants(vcf, var_type, pbsv, db, overlap, db_path, output, debug):
    # Convert the VCF file to an absolute path if not already.
    if not os.path.isabs(vcf):
        vcf = os.path.abspath(vcf)
        if debug:
            print("Absolute VCF path:", vcf)
    
    # If output directory is not provided, use the current working directory
    if not output:
        output = os.getcwd()
    else:
        output = os.path.abspath(output)
    if debug:
        print("Output directory:", output)
    
    # Determine the directory of this Python file (installed package directory)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if debug:
        print("Script directory:", script_dir)
    
    # Construct the absolute path to the Perl script (assumed in "scripts" subdirectory)
    perl_script = os.path.join(script_dir, "scripts", "findUniqVariants.pl")
    perl_script = os.path.abspath(perl_script)
    if debug:
        print("Perl script path:", perl_script)
    
    # Ensure the config file exists in the same directory as the Perl script.
    config_path = os.path.join(os.path.dirname(perl_script), "config.yaml")
    if debug:
        print("Config file path:", config_path)
    if not os.path.isfile(config_path):
        raise FileNotFoundError(f"Config file 'config.yaml' not found in {os.path.dirname(perl_script)}")
    
    # Set working directory to where the Perl script (and config.yaml) reside.
    cwd = os.path.dirname(perl_script)
    if debug:
        print("Working directory set to:", cwd)
    
    # Prepare environment variables:
    # Add the "scripts" directory (cwd) to PERL5LIB so that Perl can locate perlModules/panscan.pm.
    env = os.environ.copy()
    old_perl5lib = env.get("PERL5LIB", "")
    env["PERL5LIB"] = cwd + ":" + old_perl5lib
    if debug:
        print("Setting PERL5LIB to:", env["PERL5LIB"])
    
    # Build the command using the options expected by the Perl script.
    # (Note: We pass single-dash options to match what the Perl GetOptions expects.)
    cmd = [
        "perl", perl_script,
        "-i", vcf,
        "-t", var_type,
        "-db", db,
        "-overlap", str(overlap)
    ]
    
    if pbsv:
        cmd.append("-pbsv")
    
    if db_path:
        cmd.extend(["-db-path", db_path])
    
    # Pass the output directory.
    cmd.extend(["-output", output])
    
    if debug:
        print("Command to execute:", " ".join(cmd))
    
    # Execute the command using the specified working directory and environment.
    subprocess.run(cmd, check=True, cwd=cwd, env=env)

def main(args):
    run_find_uniq_variants(
        args.vcf_file,
        args.var_type,
        args.pbsv,
        args.db,
        args.overlap,
        args.db_path,
        args.op,
        args.debug
    )

def add_subparser(subparsers):
    parser = subparsers.add_parser(
        "find_uniq_variants",
        help="Find unique variants in a VCF file using the updated Perl script."
    )
    parser.add_argument(
        "-i", dest="vcf_file", required=True,
        help="Specify the input VCF file."
    )
    parser.add_argument(
        "-t", dest="var_type", choices=["SNP", "INDEL", "SV"], required=True,
        help="Specify the variant type to compare (SNP/INDEL/SV)."
    )
    parser.add_argument(
        "-pbsv", "--pbsv", action="store_true",
        help="Enable pbsv mode for SV analysis."
    )
    parser.add_argument(
        "-db", "--db", dest="db", default="ALL",
        help="""Specify the databases (comma-separated).
For SNP: dbSNP, gnomAD, 1000Genomes, GME.
For INDEL: gnomAD, 1000Genomes, GME.
For SV: the appropriate SV databases are used.
Default: ALL"""
    )
    parser.add_argument(
        "-overlap", "--overlap", dest="overlap", type=int, default=80,
        help="Specify the percentage of overlap for SV comparison. Default: 80"
    )
    parser.add_argument(
        "-db-path", "--db-path", dest="db_path", default=None,
        help="Specify the path to the database files. If not provided, the Perl script uses the paths from config.yaml."
    )
    parser.add_argument(
        "--op", dest="op", default=None,
        help="Optional output directory. If not provided, the output will be saved to the directory from which the command was called."
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="Enable debug mode to print diagnostic messages."
    )
    parser.set_defaults(func=main)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Find unique variants in a VCF file using the updated Perl script."
    )
    subparsers = parser.add_subparsers()
    add_subparser(subparsers)
    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)
    else:
        parser.print_help()
