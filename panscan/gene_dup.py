import argparse
import glob
import subprocess
import pandas as pd
import sys
import os
import tempfile
import uuid
import matplotlib.pyplot as plt
import pandas as pd
import math
from matplotlib_venn import venn3

MAX_RENDERER_PIXELS = (2**16 - 256)  # renderer hard cap


def _safe_square_figure(base_dpi=100, target_in=8, max_in=20):
    """
    Return (fig_w, fig_h, dpi) for a square figure that won't exceed the renderer cap.
    """
    w_in = min(target_in, max_in)
    dpi_cap = int(MAX_RENDERER_PIXELS / max(w_in, 1e-6))
    dpi = max(80, min(base_dpi, dpi_cap))
    return w_in, w_in, dpi


def plot_gene_dup_venn(genes_main, label_main, genes_ip1, label_ip1, genes_ip2, label_ip2,
                       *, out_png, out_svg=None):
    """
    genes_* are sets of gene IDs (row index). Labels are strings.
    """
    fig_w, fig_h, dpi = _safe_square_figure(base_dpi=150, target_in=8, max_in=18)

    fig = plt.figure(figsize=(fig_w, fig_h), dpi=dpi)
    v = venn3([genes_main, genes_ip1, genes_ip2],
              set_labels=(label_main, label_ip1, label_ip2))

    # Title
    plt.title("Gene Duplication Venn Diagram", fontsize=18, fontweight="bold", pad=25)

    # Style the set labels (A, B, C labels)
    for label in v.set_labels:
        if label:
            label.set_fontsize(14)
            label.set_fontweight("bold")

    # Style the subset labels (numbers)
    for label in v.subset_labels:
        if label:
            label.set_fontsize(12)
            label.set_fontweight("normal")

    # Make axis labels bold, if they exist
    plt.xlabel("X-axis", fontsize=12, fontweight="bold")  # optional
    plt.ylabel("Y-axis", fontsize=12, fontweight="bold")  # optional

    try:
        plt.tight_layout()
    except Exception:
        pass

    fig.savefig(out_png, dpi=dpi, bbox_inches="tight")
    if out_svg:
        fig.savefig(out_svg, bbox_inches="tight")
    plt.close(fig)



def _safe_bar_figure(n_bars,
                     base_dpi=100,
                     min_width_in=40,
                     max_width_in=300,
                     inches_per_bar=0.25,
                     height_in=30):
    """
    Compute a figsize + dpi that won’t exceed renderer pixel cap.
    """
    desired_w_in = max(min_width_in, min(max_width_in, n_bars * inches_per_bar))
    dpi_cap = int(MAX_RENDERER_PIXELS / max(desired_w_in, 1e-6))
    dpi = max(50, min(base_dpi, dpi_cap))
    return desired_w_in, height_in, dpi

def _thin_xticks(ax, n_bars, target_labels=100):
    """
    Show at most ~target_labels x-tick labels.
    """
    if n_bars <= target_labels:
        return
    step = max(1, math.ceil(n_bars / target_labels))
    for i, label in enumerate(ax.get_xticklabels()):
        label.set_visible((i % step) == 0)

def process_and_plot_data(csv_file, cohort_name):
    # Read input
    result_df = pd.read_csv(csv_file, index_col=0)

    # Count duplicated genes per genome
    duplicates_count = (result_df != 0).sum().sort_values()
    num_bars = len(duplicates_count)

    # Safe figure sizing
    fig_w, fig_h, dpi = _safe_bar_figure(num_bars)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h), dpi=dpi)

    # Font scaling
    base_font_size = 40
    adjusted_font_size = max(18, min(base_font_size, 1600 / max(num_bars, 1)))

    # Plot
    duplicates_count.plot(kind='bar',
                          width=0.2,
                          color='#1F1F6A',
                          edgecolor='black',
                          ax=ax)

    # Labels
    ax.set_xlabel(
        f'Genome by number of duplications ({cohort_name})',
        fontsize=adjusted_font_size,
        labelpad=adjusted_font_size * 1.2,
        weight='bold'
    )
    ax.set_ylabel(
        'Duplicated genes per genome',
        fontsize=adjusted_font_size,
        labelpad=adjusted_font_size * 1.2,
        weight='bold'
    )

    # Tick styling
    for lbl in ax.get_xticklabels():
        lbl.set_rotation(45)
        lbl.set_ha('right')
        lbl.set_fontsize(adjusted_font_size * 0.8)
    _thin_xticks(ax, num_bars)

    for lbl in ax.get_yticklabels():
        lbl.set_fontsize(adjusted_font_size * 0.8)

    # Grid + margins
    ax.grid(axis='y', linestyle='--', alpha=0.7)
    ax.margins(x=0.01)

    # Layout guardrail
    try:
        plt.tight_layout()
    except Exception:
        pass

    # Save outputs
    fig.savefig('dup-per-hap.png', dpi=dpi, bbox_inches='tight')
    if (fig_w * dpi) >= (0.9 * MAX_RENDERER_PIXELS) or num_bars > 2000:
        fig.savefig('dup-per-hap.svg', bbox_inches='tight')

    plt.close(fig)



def frequency_comparison(df1, df2, output_file, cohort_label, data_label, color1, color2, top_n):
    # Compute frequency
    df1_freq = df1.apply(lambda row: (row != 0).sum() / len(row), axis=1).to_frame(name='Frequency_x')
    df2_freq = df2.apply(lambda row: (row != 0).sum() / len(row), axis=1).to_frame(name='Frequency_y')

    # Merge and compute differences
    comparison_df = df1_freq.merge(df2_freq, left_index=True, right_index=True)
    comparison_df['Absolute_Difference'] = comparison_df['Frequency_x'] - comparison_df['Frequency_y']
    comparison_df['Percentage_Difference'] = (
        comparison_df['Absolute_Difference'] /
        comparison_df[['Frequency_x','Frequency_y']].max(axis=1)
    ) * 100
    comparison_df.sort_values('Percentage_Difference', ascending=False) \
                 .to_csv(f'Frequency-comparison-data-{data_label}.csv')

    # Find threshold yielding at least top_n genes
    threshold = 5
    while threshold >= 1:
        filtered = comparison_df[comparison_df['Percentage_Difference'] >= threshold]
        if len(filtered) >= top_n:
            break
        print(f"Lowering threshold to {threshold-1}% for {data_label}")
        threshold -= 1

    top_diff = filtered.sort_values('Percentage_Difference', ascending=False).head(top_n)
    if top_diff.empty:
        print(f"No significant differences for {data_label}")
        return

    import matplotlib.pyplot as plt
    plt.figure(figsize=(40, 5 + top_n * 5))
    plt.barh(top_diff.index, top_diff['Frequency_x'], color=color1, label=cohort_label)
    plt.barh(top_diff.index, -top_diff['Frequency_y'], color=color2, label=data_label)
    plt.xlabel('CNV frequency', labelpad=50, fontsize=70)
    plt.xticks(fontsize=60)
    plt.yticks(fontsize=60)
    plt.title(f'Frequency Comparison: {data_label}', weight='bold', fontsize=80)
    plt.legend(fontsize=60)
    plt.tight_layout(pad=4)
    plt.savefig(output_file)
    plt.show()


def generate_ideogram(ideo_ip, karyotype, ideo_dpi, debug=False):
    df = pd.read_csv(ideo_ip)
    df['Chr'] = df['Chr'].astype(str).str.replace(r'^[cC]hr','', regex=True)

    # Save temporary ideogram file
    temp_ideo = os.path.join(os.getcwd(), f"ideo_{uuid.uuid4().hex}.tsv")
    df.to_csv(temp_ideo, index=False, sep='\t')
    if debug: print("Ideo input:", temp_ideo)

    # Prepare R script
    tmpdir = os.path.join(os.getcwd(),'tmp_rscript')
    os.makedirs(tmpdir, exist_ok=True)
    rpath = os.path.join(tmpdir, f"ideo_{uuid.uuid4().hex}.R")
    rscr = f'''setwd("{os.getcwd()}")
library(RIdeogram)
karyo<-read.table("{karyotype}",sep="\t",header=TRUE,stringsAsFactors=FALSE)
ideo <- read.table("{temp_ideo}",sep="\t",header=TRUE,stringsAsFactors=FALSE)
ideogram(karyotype=karyo,overlaid=ideo,colorset1=c("#FF0000","#FF0000","#FF0000"),output="GeneDup_Ideogram.svg",width=150,Lx=10000)
svg2tiff("GeneDup_Ideogram.svg",file="GeneDup_Ideogram.tiff",height=5.4,width=7,dpi={ideo_dpi})
unlink("GeneDup_Ideogram.svg")
if(file.exists("Rplots.pdf")) unlink("Rplots.pdf")
'''    
    with open(rpath,'w') as f: f.write(rscr)
    if debug: print("R script:",rpath)

    try:
        subprocess.run(["Rscript",rpath],check=True)
    except subprocess.CalledProcessError as e:
        print("Error generating ideogram:",e)
    finally:
        os.remove(rpath)
        os.remove(temp_ideo)
        try: os.rmdir(tmpdir)
        except: pass


def run_gene_dup(args):
    # Primary bar plot
    process_and_plot_data(args.csv_file, args.cohort_name)

    # Read for comparisons
    df_main = pd.read_csv(args.csv_file, index_col=0)
    ip1_df = pd.read_csv(args.ip1, index_col=0)
    ip2_df = pd.read_csv(args.ip2, index_col=0)

    # Frequency comparisons
    frequency_comparison(df_main, ip1_df, f"freq_{args.ip1_label}.png",
                         args.cohort_name, args.ip1_label, "#1F1F6A", "#E63946", args.top_n)
    frequency_comparison(df_main, ip2_df, f"freq_{args.ip2_label}.png",
                         args.cohort_name, args.ip2_label, "#1F1F6A", "#2ECC71", args.top_n)
    
    g_main = set(df_main.index)
    g_ip1  = set(ip1_df.index)
    g_ip2  = set(ip2_df.index)

    venn_png = f"venn_{args.cohort_name}_{args.ip1_label}_{args.ip2_label}.png"
    venn_svg = f"venn_{args.cohort_name}_{args.ip1_label}_{args.ip2_label}.svg"
    plot_gene_dup_venn(g_main, args.cohort_name, g_ip1, args.ip1_label, g_ip2, args.ip2_label,
                       out_png=venn_png, out_svg=venn_svg)

    # Ideogram
    if args.ideogram:
        if not args.ideo_ip or not args.ideo_ref:
            sys.exit("--ideogram requires --ideo-ip and --ideo-ref")
        ref = args.ideo_ref.lower()
        script_dir = os.path.dirname(os.path.abspath(__file__))
        karyo = os.path.join(script_dir, f"scripts/{args.ideo_ref.upper()}_karyotype.bed")
        generate_ideogram(os.path.abspath(args.ideo_ip), karyo, args.ideo_dpi, debug=args.debug)


def add_subparser(subparsers):
    p = subparsers.add_parser("gene_dup", help="Gene duplication analysis + ideogram")
    p.add_argument("--csv_file", required=True)
    p.add_argument("--ip1", required=True, help="First comparison matrix CSV.")
    p.add_argument("--ip2", required=True, help="Second comparison matrix CSV.")
    p.add_argument("--ip1_label", default="HPRC", help="Label for --ip1.")
    p.add_argument("--ip2_label", default="CPC", help="Label for --ip2.")
    p.add_argument("--cohort_name", default="Cohort", help="Name for bar plot.")
    p.add_argument("--top_n", type=int, default=5)
    p.add_argument("--ideogram", action="store_true")
    p.add_argument("--ideo_ip", help="Ideogram input CSV.")
    p.add_argument("--ideo_ref", choices=["hg38","chm13"], help="Reference karyotype.")
    p.add_argument("--ideo_dpi", type=int, default=300)
    p.add_argument("--debug", action="store_true")
    p.set_defaults(func=run_gene_dup)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subs = parser.add_subparsers()
    add_subparser(subs)
    args = parser.parse_args()
    if hasattr(args,'func'): args.func(args)
    else: parser.print_help()
