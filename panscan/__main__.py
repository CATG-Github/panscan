import argparse
from . import complex, novel_seq, findUniqVariants, make_dup_mtx, gene_dup

def main():
    parser = argparse.ArgumentParser(description="Panscan tool to analyze complex loci, perform variant analysis, and gene duplication detection.")
    
    subparsers = parser.add_subparsers(dest="command", help="Available commands")
    # Add subcommands
    complex.add_subparser(subparsers)
    novel_seq.add_subparser(subparsers)
    findUniqVariants.add_subparser(subparsers)
    make_dup_mtx.add_subparser(subparsers)
    gene_dup.add_subparser(subparsers)

    args = parser.parse_args()
    if hasattr(args, 'func'):
        args.func(args)  # Executes the function corresponding to the subcommand
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
