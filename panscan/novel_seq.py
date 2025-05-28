import argparse
import subprocess
import os

def run_novel_seq(infile, ref, threads, dthreads, dpi, pInp, pRef, exclude, genome, output, debug):
    # Convert provided paths to absolute paths.
    if infile:
        infile = os.path.abspath(infile)
    if ref:
        ref = os.path.abspath(ref)
    if pInp:
        pInp = os.path.abspath(pInp)
    if pRef:
        pRef = os.path.abspath(pRef)
    if output:
        output = os.path.abspath(output)
    
    if debug:
        print("Input VCF (raw or preprocessed):", infile if infile else pInp)
        print("Reference VCF (raw or preprocessed):", ref if ref else pRef)
        print("Threads:", threads)
        print("RTG decompose threads:", dthreads)
        print("DPI:", dpi)
        print("Exclude:", exclude)
        print("Genome:", genome)
        print("Output directory:", output if output else "(default NovelSeq_Results)")
    
    # Determine the directory of this Python file (i.e. where the package is installed)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    if debug:
        print("Script directory:", script_dir)
    
    # Construct the absolute path to the Perl script (assumed in "scripts" subdirectory)
    perl_script = os.path.join(script_dir, "scripts", "findNovelSeq.pl")
    perl_script = os.path.abspath(perl_script)
    if debug:
        print("Perl script path:", perl_script)
    if not os.path.isfile(perl_script):
        raise FileNotFoundError(f"Perl script not found: {perl_script}")
    
    # Build the command list.
    cmd = ["perl", perl_script,
           "-t", str(threads),
           "-dt", str(dthreads),
           "-dpi", str(dpi),
           "-exclude", exclude,
           "-genome", genome]
    
    # For input VCF: either raw (--vcf) or preprocessed (--pInp) must be provided.
    if pInp and infile:
        raise ValueError("Provide either a raw input VCF (--vcf) OR a preprocessed input VCF (--pInp), not both.")
    elif pInp:
        cmd.extend(["-pInp", pInp])
    elif infile:
        cmd.extend(["-i", infile])
    else:
        raise ValueError("Please provide an input VCF file (--vcf or --pInp).")
    
    # For reference VCF: either raw (--ref) or preprocessed (--pRef) must be provided.
    if pRef and ref:
        raise ValueError("Provide either a raw reference VCF (--ref) OR a preprocessed reference VCF (--pRef), not both.")
    elif pRef:
        cmd.extend(["-pRef", pRef])
    elif ref:
        cmd.extend(["-r", ref])
    else:
        raise ValueError("Please provide a reference VCF file (--ref or --pRef).")
    
    # Add output directory flag if provided.
    if output:
        cmd.extend(["-op", output])
    
    if debug:
        print("Command to execute:", " ".join(cmd))
    
    # Run the command.
    subprocess.run(cmd, check=True)

def main(args):
    run_novel_seq(args.vcf, args.ref, args.threads, args.dt, args.dpi,
                  args.pInp, args.pRef, args.exclude, args.genome, args.op, args.debug)

def add_subparser(subparsers):
    parser = subparsers.add_parser("novel_seq", help="Generate a novel sequence FASTA from your VCF file by comparing with reference pangenome VCF.")
    # For input VCF: either raw or preprocessed.
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument(
        "-i", "--vcf",
        metavar="VCF",
        help="Raw input pangenome VCF."
    )
    inp.add_argument(
        "-pInp", "--pInp",
        metavar="VCF",
        help="Pre-processed input VCF (mutually exclusive with --vcf)."
    )
    # For reference VCF: either raw or preprocessed.
    ref = parser.add_mutually_exclusive_group(required=True)
    ref.add_argument(
        "-r", "--ref",
        metavar="VCF",
        help="Raw reference pangenome VCF."
    )
    ref.add_argument(
        "-pRef", "--pRef",
        metavar="VCF",
        help="Pre-processed reference VCF (mutually exclusive with --ref)."
    )

    parser.add_argument("-t", "--threads", type=int, default=1, help="Number of threads (default: 1).")
    parser.add_argument("-dt", "--dt", type=int, default=1, help="Number of threads for RTG decompose (default: 1).")
    parser.add_argument("-dpi", "--dpi", type=int, default=600, help="DPI for ideogram figure (default: 600).")
    parser.add_argument("-exclude", "--exclude", default="NA", help="Samples to exclude (default: NA).")
    parser.add_argument("-genome", "--genome", default="HG38", help="Genome version (HG38 or CHM13; default: HG38).")
    parser.add_argument("--op", help="Optional output directory. If not provided, the tool uses its default ('NovelSeq_Results').")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode with diagnostic messages.")
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
