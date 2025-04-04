import argparse
import pandas as pd  # type: ignore
from pathlib import Path
import subprocess
import os
import random
import re
import json
import glob
from PIL import Image, ImageDraw, ImageFont
import shutil

def extract_info(info_str, key):
    for field in info_str.split(';'):
        if field.startswith(key + '='):
            return field[len(key)+1:]
    return None

def allele_lengths(ref, alt):
    alts = alt.split(',')
    return [len(a) for a in alts]

def get_interesting_alleles(vcf_file, n_10kb, site_size):
    if n_10kb == None:
        n_10kb = 5
    if site_size == None:
        site_size=10000

    interesting_alleles = []
    with open(vcf_file, 'r') as vcf:
        for line in vcf:
            if not line.startswith('#'):  # Exclude header lines
                columns = line.strip().split('\t')
                ref = columns[3]
                alt = columns[4]
                
                lengths = allele_lengths(ref, alt)

                if len(lengths) >= n_10kb and any((l-len(ref)) > site_size for l in lengths):
                    interesting_alleles.append(line.strip())

    return interesting_alleles

def get_df(sites):
    data = []
    for site in sites:
        data.append((site.split()[0], int(site.split()[1])))
    
    return pd.DataFrame(data, columns=['chrom', 'pos'])



def find_complex_sites(vcf_file, n_10kb, site_size):
    interesting_alleles = get_interesting_alleles(vcf_file, n_10kb, site_size)
    locations = [(allele.split()[0], allele.split()[1]) for allele in interesting_alleles]
    regions = [(region[0], int(region[1] )- 10000, int(region[1]) + 10000) for region in locations]

    site_df = get_df(interesting_alleles)

    return site_df

def find_interesting_regions(site_df, window_size):
    if window_size == None:
        window_size = 100000

    regions = []
    for chrom in site_df['chrom'].unique():
        temp_df = site_df.loc[site_df['chrom'] == chrom]
        max_position = temp_df['pos'].max()
        min_position = temp_df['pos'].min()
        
        ptr = min_position
        end = max_position

        while ptr < end:
            sites_in_window = temp_df.loc[(temp_df['pos']>ptr) & (temp_df['pos']<ptr+window_size)]
            if sites_in_window.shape[0] >= 2:
                regions.append((chrom, ptr, ptr+window_size))
                
            ptr += window_size

    return regions




def produce_plottable(*task):
    gfab, gff3_file, region, cutpoints, graph_base, viz_output, connected_output, ref, chrom, start, end, query_region, workdir, gene_alignments, df_all, sep_pattern = task
    
    Path(workdir).mkdir(parents=True, exist_ok=True)
    
    # Modify query region to match gfabase requirements
    modified_query_region = f'{ref}{sep_pattern}{chrom}:{start}-{end}'
    
    connected_command = f'gfabase sub {gfab} -o {connected_output} {modified_query_region} --range --cutpoints {cutpoints} --view --connected'
    viz_command = f'gfabase sub {gfab} -o {viz_output} {modified_query_region} --range --view --cutpoints {cutpoints}'
    
    print(connected_command)
    print(viz_command)
    
    os.system(connected_command)
    os.system(viz_command)
    def run_command(command_string):
        cmd = command_string.split(' ')
        try:
            result = subprocess.run(cmd, check=True, text=True)
            return (cmd, None)
        except subprocess.CalledProcessError as e:
            return (cmd, None, e)
    if not os.path.isfile(viz_output):
        os.system(viz_command)
    
    if not os.path.isfile(connected_output):
        os.system(connected_command)
  
        
        
    def extract_genes_from_region(filename, chrom_filter, start_filter, end_filter):
        genes = []
        with open(filename, 'r') as gff:
            for line in gff:
                if not line.startswith("#"):  # Ignore comment lines
                    parts = line.strip().split("\t")
                    chrom = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    
                    # Filter for the specific region
                    if chrom == chrom_filter and parts[2] == "gene" and start_filter <= start and end_filter >= end:
                        attributes = parts[8].split(";")
                        gene_name = None
                        for attribute in attributes:
                            if attribute.startswith("gene_name="):
                                gene_name = attribute.split("=")[1]
                                break
                        if gene_name:
                            genes.append((chrom,start, end, gene_name))
        return genes

    def extract_protein_coding_genes_from_region(filename, chrom_filter, start_filter, end_filter):
        genes = []
        with open(filename, 'r') as gff:
            for line in gff:
                if not line.startswith("#"):  # Ignore comment lines
                    parts = line.strip().split("\t")
                    chrom = parts[0]
                    start = int(parts[3])
                    end = int(parts[4])
                    
                    # Filter for the specific region
                    if chrom == chrom_filter and parts[2] == "gene" and start_filter <= start and end_filter >= end:
                        attributes = parts[8].split(";")
                        gene_name = None
                        is_protein_coding = False
                        for attribute in attributes:
                            if attribute.startswith("gene_name="):
                                gene_name = attribute.split("=")[1]
                            if attribute == "gene_biotype=protein_coding":
                                is_protein_coding = True

                        if gene_name and is_protein_coding:
                            
                            genes.append((chrom,start, end, gene_name))

        return genes



    print('got_genes')

    all_genes_in_region = extract_genes_from_region(gff3_file, chrom, start, end) 
    pc_genes_in_region= extract_protein_coding_genes_from_region(gff3_file, chrom, start, end)
    all_genes = [items[-1] for items in all_genes_in_region]
    pc_genes = [items[-1] for items in pc_genes_in_region]
    
    print('getting genes')
    def hex_to_rgb(value):
        value = value.lstrip('#')
        length = len(value)
        return tuple(int(value[i:i + length // 3], 16) for i in range(0, length, length // 3))

    def rgb_to_hex(rgb):
        return '#{:02x}{:02x}{:02x}'.format(rgb[0], rgb[1], rgb[2])

    def gradient(start_hex, end_hex, num):
        start_rgb = hex_to_rgb(start_hex)
        end_rgb = hex_to_rgb(end_hex)
        colors = []
        for i in range(num):
            r = int(start_rgb[0] + (i / (num - 1)) * (end_rgb[0] - start_rgb[0]))
            g = int(start_rgb[1] + (i / (num - 1)) * (end_rgb[1] - start_rgb[1]))
            b = int(start_rgb[2] + (i / (num - 1)) * (end_rgb[2] - start_rgb[2]))
            colors.append(rgb_to_hex((r,g,b)))
        return colors
    
    
    def parse_string(s):
        return set(map(int, ''.join(c if c.isdigit() else ' ' for c in s).split()))

    
    
    def get_segement_colors_from_alignment(df_all, section_genes, color_files_dir):
        distinct_colors = [
            "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd", "#8c564b",
            "#e377c2", "#7f7f7f", "#bcbd22", "#17becf", "#393b79", "#637939",
            "#8c6d31", "#843c39", "#7b4173"
        ]
        used_colors = set()
        color_map = {}  # To store assigned colors

        def get_new_color():
            """Assign a new color: first use distinct colors, then randomize."""
            if len(used_colors) < len(distinct_colors):
                new_color = distinct_colors[len(used_colors)]
            else:
                while True:
                    new_color = "#%06x" % random.randint(0, 0xFFFFFF)
                    if new_color not in used_colors:
                        break
            used_colors.add(new_color)
            return new_color

        # Keep the original logic intact
        df_all = df_all[(df_all[0].isin(section_genes)) & (df_all[11] > 0)]
        df_all['index'] = df_all.groupby(0).cumcount() + 1
        df_all['Numbers'] = df_all[5].apply(parse_string)

        to_drop = []
        for i in range(len(df_all)):
            for j in range(i+1, len(df_all)):
                if df_all.iloc[i]['Numbers'].intersection(df_all.iloc[j]['Numbers']):
                    if df_all.iloc[i][1] < df_all.iloc[j][1]:
                        to_drop.append(df_all.iloc[i].name)
                    else:
                        to_drop.append(df_all.iloc[j].name)

        df_all.drop(to_drop, inplace=True)
        df_all.drop(columns=['Numbers'], inplace=True)

        df = df_all[[0, 5]]
        df = df.groupby(0).agg(lambda x: ''.join(x)).reset_index()

        items = []
        for i, row in df.iterrows():
            gene = str(row[0])
            nodes = re.split(r'[<>]', row[5])

            if gene not in color_map:
                color_map[gene] = get_new_color()  # Assign color only if not assigned

            with open(f'{color_files_dir}/gene_colors/{gene}.csv', 'w') as f:
                for node in nodes:
                    if node:
                        f.write(f'{node}, {color_map[gene]}\n')
                        items.append((node, color_map[gene]))

        return pd.DataFrame(items), df_all


    
    
    color_files_dir = f'{workdir}/color_csvs'
    Path(color_files_dir).mkdir(parents=True, exist_ok=True)
    Path(f'{color_files_dir}/gene_colors').mkdir(parents=True, exist_ok=True)
    all_genes, all_genes_all = get_segement_colors_from_alignment(df_all,all_genes, color_files_dir)
    all_genes.to_csv(f'{color_files_dir}/all_colors.csv', header=None, index=None)
    pc_genes, pc_genes_all = get_segement_colors_from_alignment(df_all, pc_genes, color_files_dir)
    pc_genes.to_csv(f'{color_files_dir}/pc_colors.csv', header = None, index=None)
    viz_list = []
    sample_walks = {} 
    with open(viz_output, 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            if line[0] == 'S':
                segment = line.split('\t')[1]
                viz_list.append(segment)
                
    with open(connected_output, 'r') as f:
        lines = f.read().splitlines()
        for line in lines:
            if line[0] == 'W':
                items = line.split('\t')
                nodes = re.split(r'[><]', items[6].strip())
                if items[1] not in sample_walks:
                    sample_walks[items[1]] = []
                
                sample_walks[items[1]].extend([node for node in nodes if node != '' ])

    genes_walks = {}

    for i, row in pc_genes_all.iterrows():
        gene = row[0]
        gene_index = row['index']
        path = re.split(r'[<>]',row[5])
        s_i = f'{gene}_{gene_index}'
        if gene not in genes_walks:
            genes_walks[s_i] = []
            
        genes_walks[s_i].extend(path)
        
    for k,v in sample_walks.items():
        sample_walks[k] = set(viz_list).intersection(set(v))
    for k, v in sample_walks.items():
        colors = gradient('#073f40', '#2b184d', len(v))
        solid_color = '#073f40'  # Define the solid color you want to assign
        with open(f'{color_files_dir}/{k}.walks.csv', 'w') as f:
            for item in v:
                f.write(f'{item},{solid_color}\n')
    
    # Existing sample_haps processing for protein-coding genes
    sample_haps = {}
    hap_directory = f'{workdir}/sample_haps/'
    Path(hap_directory).mkdir(parents=True, exist_ok=True)
    
    for k, v in sample_walks.items():
        if k not in sample_haps:
            sample_haps[k] = []
            
        for kj, vj in genes_walks.items():
            if len(set(v).intersection(set(vj))) > 1:
                sample_haps[k].append(kj)
                
    # Save protein-coding sample haps
    with open(f'{hap_directory}/sample_haps_pc.csv', 'w') as f:
        for i,haps in sample_haps.items():
            f.write(f'{i} : {",".join(haps)}\n')
    
    # Create a new genes_walks dictionary using all_genes_all
    all_genes_walks = {}
    for i, row in all_genes_all.iterrows():
        gene = row[0]
        gene_index = row['index']
        path = re.split(r'[<>]',row[5])
        s_i = f'{gene}_{gene_index}'
        if gene not in all_genes_walks:
            all_genes_walks[s_i] = []
            
        all_genes_walks[s_i].extend(path)
    
    # Generate sample haps for all genes
    sample_haps_all = {}
    for k, v in sample_walks.items():
        if k not in sample_haps_all:
            sample_haps_all[k] = []
            
        for kj, vj in all_genes_walks.items():
            if len(set(v).intersection(set(vj))) > 1:
                sample_haps_all[k].append(kj)
    
    # Save all genes sample haps
    with open(f'{hap_directory}/sample_haps_all.csv', 'w') as f:
        for i,haps in sample_haps_all.items():
            f.write(f'{i} : {",".join(haps)}\n')
            
    haps = open(f'{hap_directory}/sample_haps_pc.csv').read().splitlines()
    haps_dict = {hap.split(':')[0] : tuple(hap.split(':')[-1].split(',')) for hap in haps}

    
    def count_unique_lists(d):
        count_dict = {}
    
    
    
        for key, value in d.items():
            # Convert list to tuple
            t = tuple(value)
            
            if t in count_dict:
                count_dict[t]['count'] += 1
                count_dict[t]['keys'].append(key)
            else:
                count_dict[t] = {'count': 1, 'keys': [key]}
        
        return count_dict
    
    result = count_unique_lists(haps_dict)
    json_ready_result = {str(list(key)): value for key, value in result.items()}

    # Save the result to a JSON file
    with open(f'{hap_directory}/haps_counts_pc.json', 'w') as f:
        json.dump(json_ready_result, f, indent=4)

    return result

def run_panscan_complex(vcf_file, a, n, s, l, sites, sv, gfab, ref_fasta, gaf, sep_pattern, gff3_file, ref_name=None):
    if ref_name is None:
        ref_name = os.path.splitext(os.path.basename(ref_fasta))[0]


    print(f"Running complex with vcf_file={vcf_file}, a={a}, n={n}, s={s}, l={l}, sites={sites}, sv={sv}")
    site_df = find_complex_sites(vcf_file, n, s)
    complex_regions = find_interesting_regions(site_df, l)
 
    pd.DataFrame(complex_regions).to_csv("complex_regions.csv")
    graphs = [gfab]
    print('reading alignments')
    
        
    if not os.path.exists(gaf):
        os.system(f'GraphAligner -g {gfab} -f {gene_seq_fa} -t 2 -a {gaf} -x vg --multimap-score-fraction 0.1')
    
    # %%
    df_arp = pd.read_csv(gaf, sep = '\t', header=None)
    print(df_arp[0].head())

    df_arp[0] = df_arp[0].str.split('_', expand = True)[0].str.split('-', expand=True).iloc[:, 1:].fillna('').apply('-'.join, axis = 1).str.strip('-')
    print(df_arp[0].head())
    # %%
    # df_cpc = pd.read_csv('chm13.genes.cpc.gaf', sep = '\t', header=None)
    # df_cpc[0] = df_cpc[0].str.split('_', expand = True)[0].str.split('-', expand=True).iloc[:, 1:].fillna('').apply('-'.join, axis = 1).str.strip('-')


    # %%
    tasks = []
    for gfab in graphs:
        for chrom,start,end in complex_regions:
        
            # gfab = 'hprc-v2.1-mc-chm13.gfab'
            #gff3_file = "chm13v2.0_RefSeq_Liftoff_v5.1.gff3"
            fasta_file = ref_fasta
 
            region = 'complex'
            cutpoints = 1
            graph_base = os.path.splitext(os.path.basename(gfab))[0]
            
            
            gene_alignments = gaf
          
            df_all = df_arp
            
            # if gfab.split('.')[0] == 'CPC': 
            #     ref = 'CHM13v2'
            #     gene_alignments = 'chm13.genes.cpc.gaf'
            #     sep = '.'
            #     df_all = df_cpc
            
            
            
            query_region = f'{ref_name}{sep_pattern}{chrom}:{start}-{end}'

            region = f'{chrom}{sep_pattern}{start}_{end}'

            workdir = f'all_plottables/{graph_base}.{region}.wd'
            viz_output = f'{workdir}/{region}.{cutpoints}.{graph_base}.gfa'
            connected_output = f'all_walks_gfas/{region}.{cutpoints}.{graph_base}.walks.gfa'
        
            os.system('mkdir -p all_walks_gfas')
        
            tasks.append((gfab, gff3_file, region, cutpoints, graph_base, viz_output, connected_output, ref_name, chrom, int(start), int(end), query_region, workdir, gene_alignments, df_all, sep_pattern))

        # Add this line to generate and potentially save the DataFrame
         

    
    for task in tasks:
        produce_plottable(*task)

def generate_complex_sites_dataframe():
    # Initialize lists to store data
    complex_sites = []
    complex_regions = []
    genes_detected = []
    sample_haps = []
    pc_gene_counts = []
    all_gene_counts = []
    unique_hap_counts = []

    # Iterate through directories in all_plottables
    for complex_site_dir in glob.glob('all_plottables/*'):
        # Extract complex site name
        complex_site = os.path.basename(complex_site_dir)

        # Check for gene color CSV files
        gene_colors_dir = os.path.join(complex_site_dir, 'color_csvs', 'gene_colors')
        gene_color_files = glob.glob(os.path.join(gene_colors_dir, '*.csv'))

        # Determine if complex region
        has_complex_region = len(gene_color_files) > 0

        # Extract gene names
        genes = [os.path.splitext(os.path.basename(f))[0] for f in gene_color_files]
        genes_str = ','.join(genes) if genes else ''

        # Read sample haps for protein coding and all genes
        pc_sample_haps_file = os.path.join(complex_site_dir, 'sample_haps', 'sample_haps_pc.csv')
        all_sample_haps_file = os.path.join(complex_site_dir, 'sample_haps', 'sample_haps_all.csv')

        # Function to parse sample haps and get counts
        def parse_sample_haps(filepath):
            sample_haps_info = []
            gene_counts = {}

            if os.path.exists(filepath):
                with open(filepath, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line:  # Skip empty lines
                            continue

                        try:
                            sample, haps = line.split(' : ')
                            haps_list = haps.split(',')

                            # Count unique genes for each sample
                            uniquegenes = set(hap.split('')[0] for hap in haps_list)
                            gene_counts[sample] = len(unique_genes)

                            sample_haps_info.append(f"{sample}:{haps}")
                        except ValueError:
                            print(f"Skipping malformed line in {filepath}: {line}")

            return ' ; '.join(sample_haps_info), gene_counts        

        # Process sample haps
        pc_sample_haps_str, pc_gene_counts_dict = parse_sample_haps(pc_sample_haps_file)
        all_sample_haps_str, all_gene_counts_dict = parse_sample_haps(all_sample_haps_file)

        # Prepare gene count strings
        pc_gene_count_str = ' ; '.join(f"{sample}:{count}" for sample, count in pc_gene_counts_dict.items())
        all_gene_count_str = ' ; '.join(f"{sample}:{count}" for sample, count in all_gene_counts_dict.items())

        # Optional: count unique haplotypes
        unique_haps_file = os.path.join(complex_site_dir, 'sample_haps', 'haps_counts.json')
        unique_haps_count = 0
        if os.path.exists(unique_haps_file):
            with open(unique_haps_file, 'r') as f:
                haps_data = json.load(f)
                unique_haps_count = len(haps_data)

        # Store information
        complex_sites.append(complex_site)
        complex_regions.append('yes' if has_complex_region else 'no')
        genes_detected.append(genes_str)
        sample_haps.append(pc_sample_haps_str)
        pc_gene_counts.append(pc_gene_count_str)
        all_gene_counts.append(all_gene_count_str)
        unique_hap_counts.append(unique_haps_count)

    # Create DataFrame
    df = pd.DataFrame({
        'Complex Site': complex_sites,
         #'Complex Region': complex_regions,
        'Genes Detected': genes_detected,
        'Total Genes Detected': [len(genes.split(',')) if genes else 0 for genes in genes_detected],
        'Protein Coding Samples': sample_haps,
        'Protein Coding Gene Counts per Sample': pc_gene_counts,
        'All Gene Counts per Sample': all_gene_counts
    })

    #df=df.sort_values(by='Total Genes Detected', ascending=False)
    return df
            
def plot_complex_sites(workdir=None):
    """Generate colored GFA files and Bandage images for complex sites with automatic legends."""
    colored_gfa_dir = 'colored_gfas'
    bandage_img_dir = 'bandage_images'
    legend_img_dir = 'bandage_images_with_legend'
    
    Path(colored_gfa_dir).mkdir(parents=True, exist_ok=True)
    Path(bandage_img_dir).mkdir(parents=True, exist_ok=True)
    Path(legend_img_dir).mkdir(parents=True, exist_ok=True)
    
    # Find all complex sites
    if workdir:
        complex_sites = [workdir]
    else:
        complex_sites = glob.glob('all_plottables/*')
    
    for site_dir in complex_sites:
        site_name = os.path.basename(site_dir)
        color_csv_dir = os.path.join(site_dir, 'color_csvs')
        
        # Only process the GFA files in the working directory
        viz_gfa_files = glob.glob(f"{site_dir}/*.gfa")
        
        # Process each GFA file with appropriate color files
        for gfa_file in viz_gfa_files:
            gfa_base = os.path.basename(gfa_file)
            region_name = gfa_base.replace('.gfa', '')
            
            # Process with different color files
            for color_type in ['all_colors.csv', 'pc_colors.csv']:
                color_file = os.path.join(color_csv_dir, color_type)
                if os.path.exists(color_file):
                    colored_prefix = f"{color_type.split('.')[0]}_"
                    colored_gfa = os.path.join(colored_gfa_dir, f"{colored_prefix}{gfa_base}")
                    
                    # Run color addition script
                    add_color_tags(color_file, gfa_file, colored_gfa)
                    
                    # Generate Bandage image
                    output_image = os.path.join(bandage_img_dir, f"{colored_prefix}{gfa_base.replace('.gfa', '.png')}")
                    bandage_cmd = f"Bandage image {colored_gfa} {output_image} --height 3000 --nodewidth 200 --singlearr --colour custom"
                    
                    print(f"Running: {bandage_cmd}")
                    os.system(bandage_cmd)
                    
                    # Create legend from color file
                    legend_title = f"{site_name} - {'All Genes' if 'all_colors' in color_type else 'Protein Coding Genes'}"
                    legend_output = os.path.join(legend_img_dir, f"{colored_prefix}{gfa_base.replace('.gfa', '.png')}")
                    create_legend_image(color_file, output_image, legend_output, legend_title, is_walk=False)
            
            # Also process sample walks color files
            sample_walks_files = glob.glob(os.path.join(color_csv_dir, '*.walks.csv'))
            for walk_file in sample_walks_files:
                sample_name = os.path.basename(walk_file).replace('.walks.csv', '')
                colored_gfa = os.path.join(colored_gfa_dir, f"{sample_name}_{gfa_base}")
                
                # Add sample-specific colors
                add_color_tags(walk_file, gfa_file, colored_gfa)
                
                # Generate sample-specific Bandage image
                output_image = os.path.join(bandage_img_dir, f"{sample_name}_{gfa_base.replace('.gfa', '.png')}")
                bandage_cmd = f"Bandage image {colored_gfa} {output_image} --height 3000 --nodewidth 200 --singlearr --colour custom"
                
                print(f"Running: {bandage_cmd}")
                os.system(bandage_cmd)
                
                # Create legend for walk
                legend_title = f"{site_name} - Sample {sample_name} Walk"
                legend_output = os.path.join(legend_img_dir, f"{sample_name}_{gfa_base.replace('.gfa', '.png')}")
                create_legend_image(walk_file, output_image, legend_output, legend_title, is_walk=True)
    
    print(f"Bandage images with legends saved to {legend_img_dir}/")


def create_legend_image(color_file, image_path, output_path, title, is_walk=False):
    """
    Create a legend for the given color file and add it to the image.
    
    Args:
        color_file: Path to the CSV file with color information
        image_path: Path to the input image
        output_path: Path to save the output image with legend
        title: Title for the image
        is_walk: Whether this is a walk color file (True) or gene color file (False)
    """
    # Build legend dictionary from the color file
    legend_dict = {}
    gene_to_color = {}
    
    with open(color_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
                
            parts = line.split(',')
            if len(parts) >= 2:
                node_id = parts[0].strip()
                color_hex = parts[1].strip()
                
                if is_walk:
                    # For walks, we'll use a gradient of the same color
                    sample_name = os.path.basename(color_file).replace('.walks.csv', '')
                    if color_hex not in legend_dict:
                        legend_dict[color_hex] = ""#f"Sample {sample_name}"
                else:
                    # For genes, extract gene name from directory
                    gene_dir = os.path.join(os.path.dirname(color_file), 'gene_colors')
                    if os.path.exists(gene_dir):
                        gene_files = glob.glob(os.path.join(gene_dir, '*.csv'))
                        for gene_file in gene_files:
                            gene_name = os.path.basename(gene_file).replace('.csv', '')
                            # Check if this node is in the gene's file
                            with open(gene_file, 'r') as gf:
                                for gline in gf:
                                    if node_id in gline:
                                        parts = gline.split(',')
                                        if len(parts) >= 2:
                                            gene_color = parts[1].strip()
                                            gene_to_color[gene_name] = gene_color
                                            break
    
    # If it's a gene color file, use gene_to_color instead
    if not is_walk:
        legend_dict = {color: gene for gene, color in gene_to_color.items()}
    
    # Limit legend to a reasonable number of entries
    if len(legend_dict) > 20:
        # Take the first 19 entries and add a "..." entry
        keys = list(legend_dict.keys())[:19]
        trimmed_dict = {k: legend_dict[k] for k in keys}
        trimmed_dict["#000000"] = "... and more"
        legend_dict = trimmed_dict
    
    # Open the image
    image = Image.open(image_path)
    draw = ImageDraw.Draw(image)
    
    # Try to load font
    try:
        font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSansMono-Bold.ttf", 48)
        title_font = ImageFont.truetype("/usr/share/fonts/dejavu/DejaVuSansMono-Bold.ttf", 56)
    except:
        # Fall back to default font if the specified font is not available
        font = ImageFont.load_default()
        title_font = ImageFont.load_default()
    
    # Draw title
    W, H = image.size
    title_w, title_h = get_text_dimensions(title, title_font)
    title_x = (W - title_w) / 2
    draw.text((title_x, 20), title, fill="black", font=title_font)
    
    # Draw legend
    box_size = 50
    spacing = 20
    line_height = box_size + spacing
    
    # Position legend in the top-right corner
    x_start = W - 400
    y_start = 100
    
    for i, (color, label) in enumerate(legend_dict.items()):
        x0 = x_start
        y0 = y_start + i * line_height
        x1 = x0 + box_size
        y1 = y0 + box_size
        
        # Draw color box
        draw.rectangle([x0, y0, x1, y1], fill=color, outline="black")
        
        # Draw label
        text_x = x1 + spacing
        text_y = y0
        draw.text((text_x, text_y), label, fill="black", font=font)
    
    # Save the image with legend
    image.save(output_path)
    print(f"Added legend to {output_path}")


def get_text_dimensions(text, font):
    """Get the dimensions of a text using the given font."""
    try:
        # For newer versions of PIL
        return font.getbbox(text)[2:4]
    except AttributeError:
        try:
            # For older versions of PIL
            return font.getsize(text)
        except:
            # Fallback method
            dummy_img = Image.new('RGB', (1, 1))
            dummy_draw = ImageDraw.Draw(dummy_img)
            return dummy_draw.textsize(text, font=font)

# Helper function to add color tags to GFA
def add_color_tags(csv_file, gfa_file, output_file):
    """
    Read the CSV file (node_id, color_code), build a color map,
    then read the GFA, and append CL:z:<color> for matching nodes.
    """
    # 1) Build a dictionary from node_id -> color_code
    color_map = {}
    with open(csv_file, 'r') as f_csv:
        for line in f_csv:
            line = line.strip()
            if not line:
                continue
            # e.g. line might be: "1449372, #4b3710"
            parts = line.split(',')
            node_id = parts[0].strip()
            if len(parts) > 1:
                # second column is the color (e.g. "#4b3710")
                color_code = parts[1].strip()
            else:
                color_code = "green"  # default color if none provided
            color_map[node_id] = color_code
    
    # 2) Go through GFA, append CL:z:<color_code> for matching lines
    with open(gfa_file, 'r') as f_gfa, open(output_file, 'w') as f_out:
        for line in f_gfa:
            # Lines that define segments start with "S\t" per GFA 1.0
            if line.startswith('S\t'):
                columns = line.rstrip('\n').split('\t')
                node_id = columns[1]
                # If node_id is in color_map, append the color tag
                if node_id in color_map:
                    color = color_map[node_id]
                    # Append color tag: CL:z:<color>
                    columns.append(f"CL:Z:{color}")
                # Rebuild the line
                line = "\t".join(columns) + "\n"
            # Write (possibly updated) line to output
            f_out.write(line)

def main(args):
    run_panscan_complex(
        args.vcf_file, 
        args.a, 
        args.n, 
        args.s, 
        args.l, 
        args.sites, 
        args.sv, 
        args.gfab_file, 
        args.ref_fasta, 
        args.gaf_file, 
        args.sep_pattern, 
        args.gff3,
        ref_name=args.ref_name
    )
    
    # Generate summary DataFrame
    complex_sites_df = generate_complex_sites_dataframe()
    complex_sites_df.to_csv('complex_sites_summary.csv', index=False)
    
    # Generate Bandage plots if requested
    if args.plot_complex:
        plot_complex_sites()
    

    



def add_subparser(subparsers):
    parser = subparsers.add_parser("complex", help="Analyze complex loci from a VCF file.")
    parser.add_argument("vcf_file", help="Path to the VCF file generated from Minigraph-cactus pipeline.")
    parser.add_argument("gfab_file", help="Path to the gfab file generated from the GFA file")
    parser.add_argument("--ref_fasta", help="Path to the reference file" )
    parser.add_argument("--gaf_file", help="Path to the gene alignments to the graph in GAF format")
    parser.add_argument("--sep_pattern", help="Separator for the sample and the sequence as present in GAF file (eg. '#')")
    parser.add_argument('--gff3', help='Path to the gff3 file')
    parser.add_argument("-a", type=int, default=5, help="Number of alleles to define a complex site (default: 5).")
    parser.add_argument("-n", type=int, default=1, help="Number of 10kb variants to define a complex site (default: 1).")
    parser.add_argument("-s", type=int, default=10000, help="Minimum size of a variant to define a complex site (default: 10000).")
    parser.add_argument("--regions", action="store_true", help="List complex regions.")
    parser.add_argument("-l", type=int, default=100000, help="Length of the region to define a complex region (default: 100000).")
    parser.add_argument("--sites", type=int, default=1, help="Number of complex sites in a region to define it as complex (default: 1).")
    parser.add_argument("--sv", type=int, default=1, help="Number of secondary SVs in a region to define it as complex (default: 1).")
    parser.add_argument("--ref_name", default=None, help="Reference name to use in region queries as present in GAF file (eg. CHM13)(default: None)")
    parser.add_argument("--plot_complex", action="store_true", help="Generate Bandage plots for complex regions")
    parser.set_defaults(func=main)





# %%
