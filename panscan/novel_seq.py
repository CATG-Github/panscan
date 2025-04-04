import argparse
import subprocess
import os

def run_novel_seq(input_file, ref_file, threads, dpi, exclude, genome, output, debug, use_pInp, use_pRef):
    # All provided paths are assumed to be absolute by now.
    # Determine the directory of this Python file (the installed package directory)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if debug:
        print("Script directory:", script_dir)
    
    # Construct the absolute path to the Perl script (assumed in the "scripts" subdirectory)
    perl_script = os.path.join(script_dir, "scripts", "findNovelSeq.pl")
    perl_script = os.path.abspath(perl_script)
    if debug:
        print("Perl script path:", perl_script)
    
    if not os.path.isfile(perl_script):
        raise FileNotFoundError(f"Perl script not found: {perl_script}")
    
    # Build the command list.
    cmd = ["perl", perl_script]
    
    # For input: if use_pInp is True, pass -pInp; else, pass -i.
    if use_pInp:
        cmd += ["-pInp", input_file]
    else:
        cmd += ["-i", input_file]
    
    # For reference: if use_pRef is True, pass -pRef; else, pass -r.
    if use_pRef:
        cmd += ["-pRef", ref_file]
    else:
        cmd += ["-r", ref_file]
    
    # Add other flags.
    cmd += ["-t", str(threads),
            "-dpi", str(dpi),
            "-exclude", exclude,
            "-genome", genome]
    
    # Add output directory if provided.
    if output:
        cmd += ["-op", output]
    
    if debug:
        print("Command to execute:", " ".join(cmd))
    
    # Execute the command.
    subprocess.run(cmd, check=True)

def main(args):
    # Check mutual exclusivity for input options.
    if (args.vcf and args.pInp) or (not args.vcf and not args.pInp):
        raise ValueError("Please provide either a raw input VCF (--vcf) OR a preprocessed input VCF (--pInp), but not both.")
    if (args.ref and args.pRef) or (not args.ref and not args.pRef):
        raise ValueError("Please provide either a raw reference VCF (--ref) OR a preprocessed reference VCF (--pRef), but not both.")
    
    # Determine which input to use.
    if args.pInp:
        input_file = os.path.abspath(args.pInp)
        use_pInp = True
    else:
        input_file = os.path.abspath(args.vcf)
        use_pInp = False
    
    if args.pRef:
        ref_file = os.path.abspath(args.pRef)
        use_pRef = True
    else:
        ref_file = os.path.abspath(args.ref)
        use_pRef = False

    # For output directory, convert to absolute if provided.
    output = os.path.abspath(args.op) if args.op else None
    
    run_novel_seq(input_file, ref_file, args.threads, args.dpi,
                  args.exclude, args.genome, output, args.debug,
                  use_pInp, use_pRef)

def add_subparser(subparsers):
    parser = subparsers.add_parser("novel_seq", 
                                   help="Detect novel sequences from two VCF files using the updated Perl script.")
    # For input VCF: either --vcf or --pInp (mutually exclusive)
    parser.add_argument("-i", "--vcf", help="Path to the raw input pangenome VCF file.")
    parser.add_argument("-pInp", "--pInp", help="Path to the preprocessed input VCF file (mutually exclusive with --vcf).")
    # For reference VCF: either --ref or --pRef
    parser.add_argument("-r", "--ref", help="Path to the raw reference pangenome VCF file.")
    parser.add_argument("-pRef", "--pRef", help="Path to the preprocessed reference VCF file (mutually exclusive with --ref).")
    
    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads (default: 1).")
    parser.add_argument("-dpi", "--dpi", type=int, default=600, help="DPI for ideogram figure (default: 600).")
    parser.add_argument("-exclude", "--exclude", default="NA", help="Sample(s) to exclude (default: NA).")
    parser.add_argument("-genome", "--genome", default="HG38", help="Genome version (HG38 or CHM13; default: HG38).")
    parser.add_argument("--op", help="Optional output directory. If not provided, the Perl script uses its default (NovelSeq_Results).")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode to print diagnostic messages.")
    parser.set_defaults(func=main)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Detect novel sequences from two VCF files using the updated Perl script."
    )
    subparsers = parser.add_subparsers()
    add_subparser(subparsers)
    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
