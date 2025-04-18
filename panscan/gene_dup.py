import argparse
import glob
import subprocess
import pandas as pd
import sys
import os
import tempfile

def process_and_plot_data(csv_file, hprc_file, cpc_file, cohort_name):
    # Read the input files
    result_df = pd.read_csv(csv_file, index_col=0)
    hprc = pd.read_csv(hprc_file, skipfooter=1, engine='python', index_col=0)
    cpc = pd.read_csv(cpc_file, index_col=0)

    # Count non-zero values (duplicated genes per genome)
    duplicates_count = (result_df != 0).sum().sort_values()

    # Dynamically adjust figure size based on the number of bars
    num_bars = len(duplicates_count)
    fig_width = max(20, num_bars * 5)
    fig_height = 30

    import matplotlib.pyplot as plt
    plt.figure(figsize=(fig_width, fig_height))

    base_font_size = 40
    adjusted_font_size = max(30, min(base_font_size, 2000 / num_bars))
    plt.rcParams.update({'font.size': adjusted_font_size})

    padding_factor = 1.5 

    duplicates_count.plot(kind='bar', color='#1F1F6A', edgecolor='black')

    plt.xlabel(f'Genome by number of duplications ({cohort_name})',
               fontsize=adjusted_font_size, 
               labelpad=adjusted_font_size * padding_factor, weight='bold')
    plt.ylabel('Duplicated genes per genome', fontsize=adjusted_font_size, 
               labelpad=adjusted_font_size * padding_factor, weight='bold')
    plt.xticks(rotation=45, ha='right', fontsize=adjusted_font_size * 0.8)
    plt.yticks(fontsize=adjusted_font_size * 0.8)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig('dup-per-hap.png', dpi=300)
    plt.show()

    common_indices = hprc.index.intersection(result_df.index).intersection(cpc.index)
    venn_data = {
        '100': len(hprc.index.difference(result_df.index).difference(cpc.index)),
        '010': len(result_df.index.difference(hprc.index).difference(cpc.index)),
        '001': len(cpc.index.difference(hprc.index).difference(result_df.index)),
        '110': len(hprc.index.intersection(result_df.index).difference(cpc.index)),
        '101': len(hprc.index.difference(result_df.index).intersection(cpc.index)),
        '011': len(result_df.index.difference(hprc.index).intersection(cpc.index)),
        '111': len(common_indices)
    }
    plt.figure(figsize=(40, 40))
    from matplotlib_venn import venn3
    venn = venn3(subsets=venn_data, set_labels=('HPRC', cohort_name, 'CPC'))
    for text in (venn.set_labels if venn.set_labels else []):
        if text:
            text.set_fontsize(80)
    if venn.subset_labels:
        for text in venn.subset_labels:
            if text:
                text.set_fontsize(60)
    plt.title('Gene Duplication Venn Diagram', weight='bold', fontsize=90, pad=30)
    plt.savefig('venn-comparison.png')
    plt.show()

import os
import subprocess
import pandas as pd
import tempfile
import uuid

def generate_ideogram(ideo_ip, karyotype, ideo_dpi, csv_file , debug=False):

    # Read the ideogram input file (assumed tab-separated, no header)
    df = pd.read_csv(ideo_ip)
    if debug:
        print("Original ideogram input shape:", df.shape)
    
    # Process column 8 to extract copy number and gene information.
    df['cn'] = df['8'].str.extract(r'extra_copy_number=([^;]+)').astype(float)
    df['gene_id'] = df['8'].str.extract(r'gene_id=([^;]+)')
    df['gene_id'] = df['gene_id'].str.split('_').str[0]
    df['gene_name'] = df['8'].str.extract(r'gene_name=([^;]+)')
    df['gene_type'] = df['8'].str.extract(r'gene_type=([^;]+)')
    
    # Compute group size (the number of occurrences for each gene_name)
    df['group_size'] = df.groupby('gene_name')['gene_name'].transform('size')
    
    # Retain rows where:
    # - the gene appears more than once, OR
    # - it appears only once but its copy number is at least 1.
    df = df[(df['group_size'] > 1) | ((df['group_size'] == 1) & (df['cn'] >= 1))]
    df = df.drop(columns=['group_size'])
    df= df.loc[df.groupby('gene_id')['cn'].idxmax()].copy()
    df=df[df['cn']!=0]
    # Now, extract the ideogram input columns:
    # "Chr" from column 0, "Start" from column 3, and "End" from column 4.
    # Add a constant "Value" of 1000.

    # Read the gene duplication table.
    # This table has gene names as the index (or as a column); here we assume gene_name is the index.
    dup_df = pd.read_csv(csv_file, index_col=0)  # change filename as needed

    # If gene names are in the index, reset them to a column for merging.
    dup_df = dup_df.reset_index().rename(columns={'index': 'gene_name'})

    # For each gene, compute the maximum value from the numeric columns.
    # If all columns are numeric, simply:
    dup_df['max_val'] = dup_df.drop(columns='gene_name').max(axis=1)

    # Multiply the maximum value by 1000.
    dup_df['max_val'] = dup_df['max_val'] * 1000

    merged_df = pd.merge(df , dup_df[['gene_name', 'max_val']], on='gene_name', how='inner')

    # Update the ideogram Value column with the computed value.
    merged_df["Value"] = merged_df["max_val"]

    # Optionally drop the temporary max_val column.
    merged_df.drop(columns=["max_val"], inplace=True)


    # Save the processed ideogram CSV to a temporary file.
    ideo_df = pd.DataFrame({
        "Chr": merged_df['chr'],
        "Start": merged_df['3'],
        "End": merged_df['4'],
        "Value": merged_df['Value']
    })
    ideo_df["Chr"] = merged_df["chr"].astype(str).str.replace(r'\D+', '', regex=True)
        


    temp_ideo_file = os.path.join(os.getcwd(), f"temp_ideogram_{uuid.uuid4().hex}.tsv")
    ideo_df.to_csv(temp_ideo_file, index=False, sep='\t')
    if debug:
        print("Processed ideogram input saved to:", temp_ideo_file)
    
    # Create a temporary folder for the R script.
    tmp_r_dir = os.path.join(os.getcwd(), "tmp_rscript")
    os.makedirs(tmp_r_dir, exist_ok=True)
    r_script_filename = f"temp_{uuid.uuid4().hex}.R"
    r_script_path = os.path.join(tmp_r_dir, r_script_filename)
    
    # Build the R script.
    # Note: We assume Rscript is in PATH.
    r_script = f'''setwd("{os.getcwd()}")
library(RIdeogram)
karyo <- read.table("{karyotype}", sep="\\t", header=TRUE, stringsAsFactors=FALSE)
ideo <- read.table("{temp_ideo_file}", sep="\\t", header=TRUE, stringsAsFactors=FALSE)
# Force the 'Value' column to a constant value of 1000.

ideogram(karyotype=karyo, synteny=NULL, label=NULL, label_type=NULL,
         overlaid=ideo, colorset1=c("#FF0000", "#FF0000", "#FF0000"),
         output="GeneDup_Ideogram.svg", width=200, Lx=10000)
svg2tiff("GeneDup_Ideogram.svg", file="GeneDup_Ideogram.tiff", height=5.4, width=7, dpi={ideo_dpi})
file.remove("GeneDup_Ideogram.svg")
if(file.exists("Rplots.pdf")) file.remove("Rplots.pdf")
'''
    with open(r_script_path, 'w') as f:
        f.write(r_script)
    if debug:
        print("Temporary R script written to:", r_script_path)
    
    try:
        subprocess.run(["Rscript", r_script_path], check=True)
        if debug:
            print("Ideogram generated successfully as 'GeneDup_Ideogram.tiff'")
    except subprocess.CalledProcessError as e:
        print("Error generating ideogram:", e)
    finally:
        # Cleanup: remove the temporary R script file and ideogram input file.
        try:
            #os.remove(r_script_path)
            if debug:
                print("Temporary R script file removed.")
        except Exception as e:
            print("Warning: could not remove temporary R script file:", e)
        try:
            #os.remove(temp_ideo_file)
            if debug:
                print("Temporary ideogram input file removed.")
        except Exception as e:
            print("Warning: could not remove temporary ideogram input file:", e)
        try:
            #os.rmdir(tmp_r_dir)
            if debug:
                print("Temporary R script directory removed.")
        except Exception as e:
            if debug:
                print("Temporary R script directory not removed (likely not empty):", e)


def run_gene_dup(args):
    # Process gene duplication and plot results.
    #process_and_plot_data(args.csv_file, args.hprc_file, args.cpc_file, args.cohort_name)

    # Frequency comparisons, etc.—existing functionality omitted for brevity.
    
    # Ideogram generation.
    if args.ideogram:
        # The ideogram input file MUST be provided.
        if not args.ideo_ip:
            sys.exit("Error: When --ideogram is specified, --ideo-ip (ideogram input CSV) must be provided.")
        ideo_input = os.path.abspath(args.ideo_ip)
        # Instead of asking for a karyotype file, use the new --ideo-ref option.
        # It must be either 'hg38' or 'chm13'.
        if not args.ideo_ref:
            sys.exit("Error: When --ideogram is specified, --ideo-ref (hg38 or chm13) must be provided.")
        
        ref_choice = args.ideo_ref.lower()
        # Use the base module directory (where the module is installed) as the source for karyotype files.
        script_dir = os.path.dirname(os.path.abspath(__file__))
        if ref_choice == "hg38":
            karyotype = os.path.join(script_dir, "scripts/HG38_karyotype.bed")
        elif ref_choice == "chm13":
            karyotype = os.path.join(script_dir, "scripts/CHM13_karyotype.bed")
        else:
            sys.exit("Error: --ideo-ref must be either 'hg38' or 'chm13'.")
            
        
        #generate_ideogram(ideo_input, karyotype, args.ideo_dpi , args.csv_file, debug=args.debug)

def add_subparser(subparsers):
    parser = subparsers.add_parser("gene_dup", 
                                   help="Detect gene duplications from a CSV file and optionally generate an ideogram.")
    parser.add_argument("--csv_file", required=True, help="Path to the CSV file for gene duplication detection.")
    parser.add_argument("--hprc_file", default="hprc-matrix.csv", help="Path to the HPRC CSV file.")
    parser.add_argument("--cpc_file", default="cpc-matrix.csv", help="Path to the CPC CSV file.")
    parser.add_argument("--cohort_name", default="Test-Cohort", help="Optional cohort name for plots.")
    parser.add_argument("--top_n", type=int, default=5, help="Number of top genes to plot for comparison (default: 5).")
    # Ideogram options: if ideogram is desired, then both --ideo-ip and --ideo-ref are required.
    parser.add_argument("--ideogram", action="store_true", help="Generate ideogram for your duplication data.")
    parser.add_argument("--ideo-ip", help="(Required with --ideogram). Input file for the ideogram. Produced in --make-dup-mtx step as *gene.csv file.  ")
    parser.add_argument("--ideo-ref", choices=["hg38", "chm13"], help="(Required with --ideogram) Reference for karyotype: 'hg38' or 'chm13'.")
    parser.add_argument("--ideo-dpi", type=int, default=300, help="DPI for the ideogram TIFF output (default: 300).")
    parser.add_argument("--debug", action="store_true", help="Enable debug mode to print diagnostic messages.")
    parser.set_defaults(func=run_gene_dup)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run gene duplication analysis and produce plots")
    subparsers = parser.add_subparsers()
    add_subparser(subparsers)
    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()
