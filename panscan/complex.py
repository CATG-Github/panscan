import os
os.environ.setdefault("XDG_RUNTIME_DIR", "/tmp")
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

import re
import csv
import glob
import sys
import json
import hashlib
import colorsys
import argparse
import subprocess
import multiprocessing
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
import time
import threading

# ------------------------------
# Utility
# ------------------------------

import gzip
import shutil
import queue


def _decompress_if_gz(path: str, dest_dir: str = ".") -> str:
    """
    If path ends with .gz, decompress to dest_dir and return new path.
    Otherwise return path unchanged. Skips decompression if output already exists.
    """
    if not path.endswith(".gz"):
        return path
    out_path = os.path.join(dest_dir, os.path.basename(path)[:-3])
    if os.path.exists(out_path):
        print(f"[decompress] using existing {out_path}")
        return out_path
    print(f"[decompress] {path} → {out_path}")
    with gzip.open(path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return out_path


def _open_maybe_gz(path: str, mode: str = "r"):
    """Open a file that may or may not be gzipped."""
    if path.endswith(".gz"):
        return gzip.open(path, mode + "t" if "b" not in mode else mode)
    return open(path, mode)


def up_to_date(out_path, *inputs):
    if not os.path.exists(out_path):
        return False
    t_out = os.path.getmtime(out_path)
    return all(os.path.exists(p) and os.path.getmtime(p) <= t_out for p in inputs)


def _normalize_node_id(s: str) -> str:
    return s.strip()


# ------------------------------
# Gene index cache (fast extraction)
# ------------------------------
_GENE_INDEX = None  # { chrom_alias: [(start,end,name,is_pc), ...] }
_GENE_ALIASES = {}  # canonical alias map


# --- BEGIN: legacy helpers ported (trimmed) ---


def allele_lengths(ref, alt):
    return [len(a) for a in alt.split(',')]

def get_interesting_alleles(vcf_file, n_10kb, site_size):
    if n_10kb is None: n_10kb = 5
    if site_size is None: site_size = 10000
    out = []
    with _open_maybe_gz(vcf_file, 'r') as vcf:
        for line in vcf:
            if line.startswith('#'): 
                continue
            cols = line.rstrip('\n').split('\t')
            ref, alt = cols[3], cols[4]
            lens = allele_lengths(ref, alt)
            if len(lens) >= n_10kb and any((L - len(ref)) > site_size for L in lens):
                out.append(line.strip())
    return out

def _sites_df(rows):
    data = [(row.split()[0], int(row.split()[1])) for row in rows]
    return pd.DataFrame(data, columns=['chrom', 'pos'])

def find_complex_sites(vcf_file, n_10kb, site_size):
    alleles = get_interesting_alleles(vcf_file, n_10kb, site_size)
    return _sites_df(alleles)

def find_interesting_regions(site_df, window_size, min_sites=1):
    if window_size is None: window_size = 100000
    regions = []
    for chrom in site_df['chrom'].unique():
        tmp = site_df.loc[site_df['chrom'] == chrom]
        if tmp.empty:
            continue
        mn, mx = tmp['pos'].min(), tmp['pos'].max()
        ptr = mn
        while ptr <= mx:
            count = tmp.loc[(tmp['pos'] >= ptr) & (tmp['pos'] < ptr + window_size)].shape[0]
            if count >= min_sites:
                regions.append((chrom, ptr, ptr + window_size))
            ptr += window_size
    return regions

# --- The pipeline entry that complex-fast expects to call ---
def run_panscan_complex(
    vcf_file, a, n, s, l, sites, sv, gfab, ref_fasta, gaf_file,
    sep_pattern, gff3_file, ref_name=None, regions=False, no_walks=False,
    max_workers=0, min_free_gb=20.0, **_unused
):
    # default ref name from FASTA if user didn't pass one
    if ref_name is None and ref_fasta:
        ref_name = os.path.splitext(os.path.basename(ref_fasta))[0]
    if ref_name is None:
        ref_name = "CHM13"

    # decompress inputs if gz
    vcf_file = _decompress_if_gz(vcf_file)
    gfab     = _decompress_if_gz(gfab)
    if gaf_file:
        gaf_file = _decompress_if_gz(gaf_file)
    if gff3_file:
        gff3_file = _decompress_if_gz(gff3_file)

    print(f"[complex-fast] scanning complex sites: vcf={vcf_file} a={a} n={n} s={s} l={l} sites={sites} sv={sv}")
    site_df = find_complex_sites(vcf_file, n, s)
    print(f"[debug] complex sites found: {len(site_df)}")
    print(site_df)
    complex_regions = find_interesting_regions(site_df, l, min_sites=sites)
    print(f"[debug] complex regions found: {len(complex_regions)}")
    print(complex_regions)
    pd.DataFrame(complex_regions).to_csv("complex_regions.csv", index=False)

    # Read provided GAF (do NOT try to regenerate here)
    if not gaf_file or not os.path.exists(gaf_file):
        raise FileNotFoundError(f"--gaf_file not found: {gaf_file}")
    
    # normalize gene/sample column to match legacy logic
    def normalize_gaf_query(q: str) -> str:
        q = str(q).strip()
        # drop coordinates suffix
        q = q.split("::", 1)[0]
        # handle optional "gene-" prefix
        if q.startswith("gene-"):
            q = q[len("gene-"):]
        # some pipelines add a trailing underscore before "::"
        if q[-1]=="_":
            q = q.rstrip("_")
        # if anything still has underscores, keep the first token (common in your good format)
        # q = q.split("_", 1)[0]
        return q

    df_arp = pd.read_csv(gaf_file, sep="\t", header=None)
    df_arp[0] = df_arp[0].map(normalize_gaf_query)

    tasks = []
    graph_base = os.path.splitext(os.path.basename(gfab))[0]
    os.makedirs('all_walks_gfas', exist_ok=True)

    for (chrom, start, end) in complex_regions:
        region = f"{chrom}{sep_pattern}{start}_{end}"
        workdir = f"all_plottables/{graph_base}.{region}.wd"
        cutpoints = 1
        viz_output = f"{workdir}/{region}.{cutpoints}.{graph_base}.gfa"
        connected_output = f"all_walks_gfas/{region}.{cutpoints}.{graph_base}.walks.gfa"

        tasks.append((
            gfab, gff3_file, region, cutpoints, graph_base, viz_output, connected_output,
            ref_name, chrom, int(start), int(end),
            f"{ref_name}{sep_pattern}{chrom}:{start}-{end}",  # query_region (for logging)
            workdir, gaf_file, df_arp, sep_pattern, no_walks
        ))

    total_workers  = max_workers or max(1, multiprocessing.cpu_count())
    # auto split: 90% produce, 10% plot (min 1 each)
    plot_workers    = max(1, total_workers // 10)
    produce_workers = max(1, total_workers - plot_workers)
    print(f"[complex] {len(tasks)} region(s) | produce_workers={produce_workers} plot_workers={plot_workers}")

    memgate = _PlotMemoryGate(min_free_gb=min_free_gb, check_sec=2.0, label="produce_plottable")

    # queue connecting producers to consumer plot pool
    plot_queue = queue.Queue()
    _SENTINEL = object()

    def _run_task(t):
        memgate.wait_for_memory()
        produce_plottable(*t)
        # signal consumer: this site is ready to plot
        # extract site workdir from task tuple for plot_complex_sites
        plot_queue.put(t)

    def _plot_consumer():
        """Consume completed regions and plot as they arrive."""
        while True:
            item = plot_queue.get()
            if item is _SENTINEL:
                plot_queue.task_done()
                break
            try:
                # unpack enough to identify the workdir
                (t_gfab, t_gff3, t_region, t_cutpoints, t_graph_base,
                 t_viz_output, t_connected_output, t_ref, t_chrom,
                 t_start, t_end, t_query_region, t_workdir,
                 t_gene_alignments, t_df_all, t_sep_pattern, t_no_walks) = item
                plot_complex_sites(
                    workdir=t_workdir,
                    no_walks=t_no_walks,
                    no_images=False,
                    max_workers=plot_workers,
                    plot_min_free_gb=min_free_gb,
                )
            except Exception as e:
                print(f"[warn] plot_complex_sites failed: {e}")
            finally:
                plot_queue.task_done()

    # start plot consumer threads
    consumer_threads = []
    for _ in range(plot_workers):
        t = threading.Thread(target=_plot_consumer, daemon=True)
        t.start()
        consumer_threads.append(t)

    # run producers
    with ThreadPoolExecutor(max_workers=produce_workers) as ex:
        futs = [ex.submit(_run_task, t) for t in tasks]
        for fut in as_completed(futs):
            try:
                fut.result()
            except Exception as e:
                print(f"[error] produce_plottable failed: {e}")

    # signal consumers to stop
    for _ in range(plot_workers):
        plot_queue.put(_SENTINEL)
    for t in consumer_threads:
        t.join()

# --- END: port ---

def _normalize_args(args):
    # allow either gfab or gfab_file
    if hasattr(args, "gfab_file") and not hasattr(args, "gfab"):
        args.gfab = args.gfab_file
    if hasattr(args, "gfab") and not hasattr(args, "gfab_file"):
        args.gfab_file = args.gfab

    # unify threads/max_workers aliases
    if getattr(args, "threads", None) and not getattr(args, "max_workers", 0):
        args.max_workers = args.threads
    if not hasattr(args, "sep_pattern"):
        args.sep_pattern = "#0#"
    return args

def _build_gene_index(gff3_path):
    global _GENE_INDEX, _GENE_ALIASES
    _GENE_INDEX, _GENE_ALIASES = {}, {}

    biotype_keys = ("gene_biotype", "biotype", "gene_type")
    name_keys = ("gene_name", "Name", "gene", "ID")

    with open(gff3_path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            chrom = parts[0]
            try:
                s, e = int(parts[3]), int(parts[4])
            except Exception:
                continue
            # parse attributes
            attr = {}
            for kv in parts[8].split(";"):
                if "=" in kv:
                    k, v = kv.split("=", 1)
                    attr[k] = v

            # name
            gname = None
            for k in name_keys:
                if k in attr and attr[k]:
                    gname = attr[k]
                    break
            if not gname:
                continue

            # biotype
            is_pc = False
            for k in biotype_keys:
                if k in attr and str(attr[k]).lower() == "protein_coding":
                    is_pc = True
                    break

            # alias both chrX and X
            aliases = {chrom}
            if chrom.startswith("chr"):
                aliases.add(chrom[3:])
            else:
                aliases.add("chr" + chrom)

            for al in aliases:
                _GENE_INDEX.setdefault(al, []).append((s, e, gname, is_pc))
                _GENE_ALIASES[al] = aliases

    # sort by start
    for al in _GENE_INDEX:
        _GENE_INDEX[al].sort(key=lambda x: x[0])


def extract_genes_from_region_fast(chrom, start, end, protein_coding=None):
    if _GENE_INDEX is None:
        return []
    out = set()
    chroms = _GENE_ALIASES.get(chrom, {chrom})
    for c in chroms:
        lst = _GENE_INDEX.get(c, [])
        if not lst:
            continue
        # linear scan with early break (sorted by start)
        for s, e, g, pc in lst:
            if s > end:
                break
            if e < start:
                continue
            if protein_coding is True and not pc:
                continue
            if protein_coding is False and pc:
                continue
            out.add(g)
    return sorted(out)


# ------------------------------
# Coloring helpers
# ------------------------------
PALETTE = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#ff9896", "#bcbd22", "#17becf",
    "#393b79", "#637939", "#8c6d31", "#843c39", "#7b4173",
    "#dbdb8d", "#9edae5", "#9edae5", "#ff9896", "#98df8a",
]
GOLDEN_ANGLE = 0.61803398875


def _hash_idx(g):
    return int(hashlib.md5(g.encode("utf-8")).hexdigest()[:8], 16) % len(PALETTE)


def _fallback_color(k):
    h = (k * GOLDEN_ANGLE) % 1.0
    l, s = 0.50, 0.65
    r, g, b = colorsys.hls_to_rgb(h, l, s)
    return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))


def get_segement_colors_from_alignment(df_all, section_genes, color_files_dir, node_filter=None):
    """
    Build per-gene node->color CSVs in color_files_dir/gene_colors and
    return a DataFrame of [(node, color)] rows.
    Colors are visually distinct, collision-free per site, stable across All/PC.
    """
    section_genes = [str(g) for g in sorted(set(section_genes))]
    if not section_genes:
        return pd.DataFrame([], columns=[0, 1]), df_all

    df_all = df_all.copy()
    df_all = df_all[(df_all[0].astype(str).isin(section_genes)) & (df_all[11] > 0)]
    df = (
        df_all[[0, 5]]
        .groupby(0, sort=True)[5]
        .apply("".join)
        .reset_index()
    )

    # assign colors deterministically, avoid reuse within site
    assigned, taken = {}, set()
    for g in section_genes:
        idx = _hash_idx(g)
        for step in range(len(PALETTE)):
            cand = PALETTE[(idx + step) % len(PALETTE)]
            if cand not in taken:
                assigned[g] = cand
                taken.add(cand)
                break
        else:
            k = len(assigned)
            col = _fallback_color(k)
            while col in taken:
                k += 1
                col = _fallback_color(k)
            assigned[g] = col
            taken.add(col)

    os.makedirs(os.path.join(color_files_dir, "gene_colors"), exist_ok=True)

    items = []
    for _, row in df.iterrows():
        gene = str(row[0])
        col = assigned[gene]
        nodes = [n for n in re.split(r"[<>]", row[5]) if n]
        nodes = [_normalize_node_id(n) for n in nodes]
        if node_filter is not None:
            nodes = [n for n in nodes if n in node_filter]
        with open(os.path.join(color_files_dir, "gene_colors", f"{gene}.csv"), "w") as f:
            for nid in nodes:
                f.write(f"{nid}, {col}\n")
                items.append((nid, col))

    return pd.DataFrame(items), df_all


# ------------------------------
# Tagging + legend extraction
# ------------------------------

def add_color_tags(csv_file, gfa_file, output_file):
    # node -> color
    color_map = {}
    with open(csv_file, "r") as f_csv:
        for line in f_csv:
            line = line.strip()
            if not line:
                continue
            parts = [p.strip() for p in line.split(",")]
            node_id = parts[0]
            color_code = parts[1] if len(parts) > 1 else "#00ff00"
            if color_code and not color_code.startswith("#"):
                color_code = "#" + color_code
            color_map[node_id] = color_code

    # node -> gene (from per-gene CSVs)
    node2gene = {}
    gene_dir = os.path.join(os.path.dirname(csv_file), "gene_colors")
    if os.path.isdir(gene_dir):
        for gene_csv in glob.glob(os.path.join(gene_dir, "*.csv")):
            gene_name = os.path.splitext(os.path.basename(gene_csv))[0]
            try:
                with open(gene_csv, "r") as f:
                    for ln in f:
                        ln = ln.strip()
                        if not ln:
                            continue
                        parts = [p.strip() for p in ln.split(",")]
                        nid = parts[0]
                        if nid and gene_name and nid not in node2gene:
                            node2gene[nid] = gene_name
            except OSError:
                continue

    added_cl = added_gn = 0
    with open(gfa_file, "r") as fin, open(output_file, "w") as fout:
        for line in fin:
            if not line.startswith("S\t"):
                fout.write(line)
                continue

            cols = line.rstrip("\n").split("\t")
            node_id = cols[1]
            tags = cols[3:] if len(cols) > 3 else []

            has_cl = any(t.startswith("CL:Z:") for t in tags)
            has_gn = any(t.startswith("GN:Z:") for t in tags)

            if node_id in color_map and not has_cl:
                tags.append(f"CL:Z:{color_map[node_id]}")
                added_cl += 1

            if node_id in color_map and not has_gn:
                gene = node2gene.get(node_id, f"node_{node_id}")
                tags.append(f"GN:Z:{gene}")
                added_gn += 1

            line = "\t".join(cols[:3] + tags) + "\n"
            fout.write(line)

    print(
        f"[add_color_tags] Wrote {output_file} | +CL={added_cl}, +GN={added_gn}, nodes_colored={len(color_map)}"
    )


def extract_gene_colors_from_gfa(gfa_path):
    gene2color = {}
    if not os.path.exists(gfa_path):
        return gene2color
    with open(gfa_path) as f:
        for line in f:
            if not line.startswith("S\t"):
                continue
            parts = line.rstrip("\n").split("\t")
            tags = parts[3:] if len(parts) > 3 else []
            gene = None
            color = None
            for t in tags:
                if t.startswith("GN:Z:"):
                    g = t.split(":", 2)[-1]
                    if g and not g.startswith("node_"):
                        gene = g
                if t.startswith("CL:Z:"):
                    color = t.split(":", 2)[-1]
            if gene and color:
                gene2color[gene] = color
    return gene2color


# ------------------------------
# Walks (per-sample)
# ------------------------------

def _build_sample_walk_csvs(walks_gfa, color_csv_dir, all_colors_csv, viz_nodes=None):
    if not os.path.exists(walks_gfa):
        return 0

    node2col = {}
    if os.path.exists(all_colors_csv):
        with open(all_colors_csv) as f:
            for ln in f:
                ln = ln.strip()
                if not ln:
                    continue
                parts = [p.strip() for p in ln.split(",")[:2]
                        ]  # node,color CSV
                if len(parts) == 2:
                    nid, col = parts
                    if col and not col.startswith("#"):
                        col = "#" + col
                    node2col[nid] = col

    made = 0
    # key = (sample, hap) to keep haplotypes separate — fixes collapse bug
    sample_hap_nodes = {}

    with open(walks_gfa) as f:
        for ln in f:
            if not ln.startswith("W\t"):
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            sample_raw = parts[1].strip()
            hap        = parts[2].strip()
            path       = parts[6]
            # sanitize sample name but keep hap separate
            sample = re.sub(r'[^A-Za-z0-9._+-]+', '_', sample_raw)
            key = f"{sample}.hap{hap}"
            nodes = [_normalize_node_id(n) for n in re.split(r"[<>]", path) if n]
            sample_hap_nodes.setdefault(key, set()).update(nodes)

    os.makedirs(color_csv_dir, exist_ok=True)
    for hap_key, nodes in sample_hap_nodes.items():
        out_csv = os.path.join(color_csv_dir, f"{hap_key}.walks.csv")
        with open(out_csv, "w") as w:
            for n in nodes:
                if viz_nodes and n not in viz_nodes:
                    continue
                col = node2col.get(n, "#000000")
                w.write(f"{n}, {col}\n")
        made += 1

    print(f"[walks] built {made} per-haplotype walk CSVs in {color_csv_dir}")
    return made


# ------------------------------
# Legend rendering (scalable, top-centered title)
# ------------------------------

def __text_size(draw, text, font):
    try:
        box = draw.textbbox((0, 0), text, font=font)
        return (box[2] - box[0], box[3] - box[1])
    except Exception:
        try:
            box = font.getbbox(text)
            return (box[2] - box[0], box[3] - box[1])
        except Exception:
            return (int(0.6 * len(text) * 16), 24)


def create_legend_image(color_file, image_path, output_path, title, is_walk=False, colored_gfa=None):
    from collections import Counter

    legend = {}
    try:
        with open(color_file, "r", newline="") as cf:
            r = csv.DictReader(cf)
            headers = {h.strip().lower() for h in (r.fieldnames or [])}
            if {"gene", "color"}.issubset(headers):
                for row in r:
                    g = (row.get("gene") or "").strip()
                    c = (row.get("color") or "").strip()
                    if g and c:
                        if not c.startswith("#"):
                            c = "#" + c
                        legend[g] = c
            else:
                raise ValueError("no headered gene,color")
    except Exception:
        #if is_walk:
        #    sample = os.path.splitext(os.path.basename(color_file))[0]
        #    cols = []
        #    with open(color_file) as f:
        #        for ln in f:
        #            ln = ln.strip()
        #            if not ln:
        #                continue
        #            parts = [p.strip() for p in ln.split(",")]
        #            if len(parts) >= 2 and parts[1].startswith("#"):
        #                cols.append(parts[1])
        #    c = Counter(cols).most_common(1)[0][0] if cols else "#000000"
        #    legend = {sample: c}
        if is_walk:
                sample = os.path.splitext(os.path.basename(color_file))[0]
                # load gene2color from the colored GFA's GN:Z tag
                gene2color = {}
                if colored_gfa and os.path.exists(colored_gfa):
                        gene2color = extract_gene_colors_from_gfa(colored_gfa)
                if gene2color:
                        legend = gene2color  # {gene_name: color}
                        legend["(traversed, no gene)"] = "#000000"  # ← ADD HERE
                else:
                # fallback to single dominant color
                        cols = []
                        with open(color_file) as f:
                           for ln in f:
                                ln = ln.strip()
                                if not ln:
                                         continue
                                parts = [p.strip() for p in ln.split(",")]
                                if len(parts) >= 2 and parts[1].startswith("#"):
                                         cols.append(parts[1])
                        c = Counter(cols).most_common(1)[0][0] if cols else "#000000"
                        legend = {sample: c}
        else:
                legend = {}

    base = Image.open(image_path).convert("RGBA")
    W, H = base.width, base.height
    draw = ImageDraw.Draw(base)

    s = max(1.0, H / 2000.0)
    title_pt = int(52 * s)
    label_pt = int(38 * s)
    swatch_px = int(34 * s)
    pad = int(12 * s)

    try:
        font = ImageFont.truetype("DejaVuSans.ttf", label_pt)
        fontB = ImageFont.truetype("DejaVuSans-Bold.ttf", title_pt)
    except Exception:
        font = fontB = ImageFont.load_default()

    # title
    tw, th = __text_size(draw, title, fontB)
    tx = max(0, (W - tw) // 2)
    draw.rectangle([tx - pad, pad, tx + tw + pad, pad + th + pad], fill=(255, 255, 255, 220))
    draw.text((tx, pad + pad // 2), title, fill=(0, 0, 0, 255), font=fontB)

    entries = list(legend.items())
    if entries:
        max_label_w = 0
        for lbl, _c in entries:
            lw, _ = __text_size(draw, str(lbl), font)
            if lw > max_label_w:
                max_label_w = lw
        panel_w = pad + swatch_px + pad + max_label_w + pad
        panel_h = pad + len(entries) * (max(swatch_px, label_pt) + pad) + pad
        x0 = W - panel_w - pad
        y0 = th + 3 * pad
        draw.rectangle([x0, y0, x0 + panel_w, y0 + panel_h], fill=(255, 255, 255, 220), outline=(0, 0, 0, 80))

        y = y0 + pad
        for lbl, c in entries:
            draw.rectangle([x0 + pad, y, x0 + pad + swatch_px, y + swatch_px], fill=c, outline=(0, 0, 0, 180))
            draw.text((x0 + pad + swatch_px + pad, y), str(lbl), fill=(0, 0, 0, 255), font=font)
            y += max(swatch_px, label_pt) + pad

    base.save(output_path)
    print(f"Added legend to {output_path}")


# ------------------------------
# Bandage runner (parallel)
# ------------------------------

#def _auto_workers(user_max):
#    if user_max and user_max > 0:
#        return user_max
#    return max(2, min(6, multiprocessing.cpu_count() // 3))

def _cpu_auto_workers():
    """Safe default when the user doesn't specify a thread budget."""
    return max(2, min(6, multiprocessing.cpu_count() // 3))

def _effective_plot_workers(user_threads: int, plot_threads: int) -> int:
    """
    Decide how many parallel Bandage runs to allow.

    Rules:
      - If --plot_threads is set (>0), use it.
      - Else if user provided a thread budget (e.g. --threads / --max_workers),
        use 1/3 of that budget (min 1).
      - Else fall back to a small CPU-based default.
    """
    if plot_threads and plot_threads > 0:
        return max(1, int(plot_threads))
    if user_threads and user_threads > 0:
        return max(1, int(user_threads) // 3)
    return _cpu_auto_workers()

def _mem_available_bytes() -> int:
    """Return available system memory in bytes (Linux-friendly, no psutil required)."""
    try:
        import psutil  # type: ignore
        return int(psutil.virtual_memory().available)
    except Exception:
        pass

    try:
        with open("/proc/meminfo") as f:
            info = {}
            for ln in f:
                if ":" not in ln:
                    continue
                k, v = ln.split(":", 1)
                info[k.strip()] = v.strip()
        if "MemAvailable" in info:
            kb = int(info["MemAvailable"].split()[0])
            return kb * 1024
    except Exception:
        pass

    return 0


class _PlotMemoryGate:
    """Blocks launching new Bandage processes when free RAM is too low."""

    def __init__(self, min_free_gb: float = 0.0, check_sec: float = 1.0, label: str = "Bandage"):
        self.min_free_bytes = int(max(0.0, float(min_free_gb)) * (1024 ** 3))
        self.check_sec = max(0.1, float(check_sec))
        self.label = label
        self._lock = threading.Lock()
        self._last_log = 0.0

    def wait_for_memory(self):
        if self.min_free_bytes <= 0:
            return
        while True:
            avail = _mem_available_bytes()
            if avail <= 0 or avail >= self.min_free_bytes:
                return
            now = time.time()
            with self._lock:
                if now - self._last_log >= 5.0:
                    self._last_log = now
                    need_gb = self.min_free_bytes / (1024 ** 3)
                    avail_gb = avail / (1024 ** 3)
                    print(f"[memwatch] {self.label}: waiting for RAM (free {avail_gb:.2f} GB < {need_gb:.2f} GB)")
            time.sleep(self.check_sec)
 

def _run_bandage(colored_gfa, out_png):
    env = os.environ.copy()
    env.setdefault("XDG_RUNTIME_DIR", "/tmp")
    env.setdefault("QT_QPA_PLATFORM", "offscreen")
    cmd = [
        "Bandage",
        "image",
        colored_gfa,
        out_png,
        "--height",
        "3000",
        "--nodewidth",
        "200",
        "--singlearr",
        "--colour",
        "custom",
    ]
    subprocess.run(cmd, check=True, env=env)

def _run_bandage_with_memwatch(colored_gfa, out_png, memgate):
    if memgate is not None:
        memgate.wait_for_memory()
    _run_bandage(colored_gfa, out_png)


# ------------------------------
# Pipeline steps
# ------------------------------

def produce_plottable(*task):
    (
        gfab,
        gff3,
        region,
        cutpoints,
        graph_base,
        viz_output,
        connected_output,
        ref,
        chrom,
        start,
        end,
        query_region,
        workdir,
        gene_alignments,
        df_all,
        sep_pattern,
        no_walks,
    ) = task

    Path(workdir).mkdir(parents=True, exist_ok=True)

    modified_query_region = f"{ref}{sep_pattern}{chrom}:{start}-{end}"

    # ── viz GFA (no --connected, fast regional subgraph) ────────────────────
    viz_command = f"gfabase sub {gfab} -o {viz_output} {modified_query_region} --range --view --cutpoints {cutpoints}"
    print(viz_command)
    if not os.path.isfile(viz_output):
        os.system(viz_command)
    if not os.path.isfile(viz_output):
        print(f"[warn] viz GFA missing after extraction: {viz_output}")

    # ── walks GFA: extract full connected → filter W lines → compact ────────
    # tmp file keeps the full connected GFA for reference/debugging
    tmp_connected = connected_output.replace(".walks.gfa", ".walks.tmp.gfa")
    connected_command = (f"gfabase sub {gfab} -o {tmp_connected} "
                         f"{modified_query_region} --range --cutpoints {cutpoints} "
                         f"--view --connected")
    print(connected_command)
    if not os.path.isfile(tmp_connected):
        os.system(connected_command)
    if not os.path.isfile(tmp_connected):
        print(f"[warn] tmp connected GFA missing: {tmp_connected}")

    # build regional node set from viz GFA
    # AFTER
    regional_nodes = set()
    if os.path.exists(tmp_connected):
        with open(tmp_connected) as f:
            for ln in f:
                if ln.startswith("S\t"):
                    regional_nodes.add(ln.split("\t", 2)[1].strip())
    print(f"[debug] regional nodes from connected GFA: {len(regional_nodes)}")
    viz_nodes = set()
    if os.path.exists(viz_output):
        with open(viz_output) as f:
            for ln in f:
                if ln.startswith("S\t"):
                    viz_nodes.add(ln.split("\t", 2)[1].strip())
    print(f"[debug] viz nodes: {len(viz_nodes)}")

    # filter W lines from tmp connected GFA to only those traversing regional nodes
    w_written = 0
    if os.path.exists(tmp_connected) and regional_nodes:
        with open(connected_output, "w") as fout:
            # copy S/E/H lines from viz GFA
            with open(viz_output) as f:
                for ln in f:
                    fout.write(ln)
            # append filtered W lines
            with open(tmp_connected) as f:
                for ln in f:
                    if not ln.startswith("W\t"):
                        continue
                    parts = ln.rstrip("\n").split("\t")
                    if len(parts) < 7:
                        continue
                    walk_nodes = set(n for n in re.split(r"[<>]", parts[6]) if n)
                    if walk_nodes & regional_nodes:
                        fout.write(ln)
                        w_written += 1
        print(f"[debug] compact walks GFA: {w_written} W lines → {connected_output}")
    elif not regional_nodes:
        print(f"[warn] no regional nodes found, skipping W line filtering")

    # Extract genes fast
    all_genes = extract_genes_from_region_fast(chrom, start, end, protein_coding=None)
    pc_genes = extract_genes_from_region_fast(chrom, start, end, protein_coding=True)

    site_name = f"{ref}.{chrom}#{start}_{end}"

    # Debug
    debug_dir = os.path.join(workdir, "debug_lists")
    os.makedirs(debug_dir, exist_ok=True)
    with open(os.path.join(debug_dir, f"{site_name}.all_genes.txt"), "w") as w:
        for g in sorted(set(all_genes)):
            w.write(g + "\n")
    with open(os.path.join(debug_dir, f"{site_name}.pc_genes.txt"), "w") as w:
        for g in sorted(set(pc_genes)):
            w.write(g + "\n")
    print(f"[debug] {site_name}: all={len(all_genes)}, pc={len(pc_genes)} (lists in {debug_dir})")

    # Prepare per-site dirs
    color_csv_dir = os.path.join(workdir, "color_csvs")
    Path(color_csv_dir).mkdir(parents=True, exist_ok=True)
    Path(os.path.join(color_csv_dir, "gene_colors")).mkdir(parents=True, exist_ok=True)

    # ── FIX: node_filter from connected GFA nodes (not viz GFA) ─────────────
    node_filter = set()
    if os.path.exists(connected_output):
        with open(connected_output) as f:
            for ln in f:
                if ln.startswith("S\t"):
                    node_filter.add(ln.split("\t", 2)[1].strip())
    print(f"[debug] node_filter from connected GFA: {len(node_filter)} nodes")

    # Build colors + CSVs
    all_df, _ = get_segement_colors_from_alignment(
        df_all, all_genes, color_csv_dir, node_filter=node_filter
    )
    all_csv = os.path.join(color_csv_dir, "all_colors.csv")
    if not all_df.empty:
        all_df.to_csv(all_csv, header=False, index=False)

    pc_df, _ = get_segement_colors_from_alignment(
        df_all, pc_genes, color_csv_dir, node_filter=node_filter
    )
    pc_csv = os.path.join(color_csv_dir, "pc_colors.csv")
    if not pc_df.empty:
        pc_df.to_csv(pc_csv, header=False, index=False)

    print(f"[debug] wrote color CSVs in {color_csv_dir} (all={len(all_df)} rows, pc={len(pc_df)} rows)")

    # Build per-haplotype walk CSVs (skip when --no-walks)
    if not no_walks:
        _ = _build_sample_walk_csvs(connected_output, color_csv_dir, all_csv, viz_nodes=viz_nodes)
    else:
        print('[debug] --no-walks set: skipping walk CSV generation')


def plot_complex_sites(workdir=None, no_images=False, no_walks=False, max_workers=0, plot_threads=0, plot_min_free_gb: float = 0.0, plot_mem_check_sec: float = 1.0):
    if workdir:
        complex_sites = [workdir]
    else:
        complex_sites = glob.glob("all_plottables/*")

    for site_dir in complex_sites:
        site_name = os.path.basename(site_dir)
        colored_gfa_dir = os.path.join(site_dir, "colored_gfas")
        bandage_img_dir = os.path.join(site_dir, "bandage_images")
        legend_img_dir = os.path.join(site_dir, "bandage_images_with_legend")
        color_csv_dir = os.path.join(site_dir, "color_csvs")

        for d in (colored_gfa_dir, bandage_img_dir, legend_img_dir):
            Path(d).mkdir(parents=True, exist_ok=True)

        viz_gfa_files = glob.glob(f"{site_dir}/*.gfa")
        bandage_jobs = []  # (colored_gfa, output_png, legend_src, title, is_walk)

        for gfa_file in viz_gfa_files:
            gfa_base = os.path.basename(gfa_file)

            for color_type in ["all_colors.csv", "pc_colors.csv"]:
                color_file = os.path.join(color_csv_dir, color_type)
                if not os.path.exists(color_file):
                    continue
                prefix = f"{color_type.split('.')[0]}_"
                colored_gfa = os.path.join(colored_gfa_dir, f"{prefix}{gfa_base}")
                output_png = os.path.join(
                    bandage_img_dir, f"{prefix}{gfa_base.replace('.gfa', '.png')}"
                )
                legend_png = os.path.join(
                    legend_img_dir, f"{prefix}{gfa_base.replace('.gfa', '.png')}"
                )
                legend_csv = os.path.join(legend_img_dir, f"{prefix}{gfa_base}.legend.csv")
                legend_title = (
                    f"{site_name} - {'All Genes' if 'all_colors' in prefix else 'Protein Coding Genes'}"
                )

                if not up_to_date(colored_gfa, gfa_file, color_file):
                    add_color_tags(color_file, gfa_file, colored_gfa)

                gene2color = extract_gene_colors_from_gfa(colored_gfa)
                with open(legend_csv, "w") as w:
                    w.write("gene,color\n")
                    for g, c in sorted(gene2color.items()):
                        w.write(f"{g},{c}\n")
                print(
                    f"[legend] {os.path.basename(colored_gfa)}: {len(gene2color)} genes found for legend"
                )

                if no_images:
                    continue
                if not up_to_date(output_png, colored_gfa):
                    bandage_jobs.append((colored_gfa, output_png, legend_csv, legend_title, False))

        # sample walks — parallel add_color_tags with memory gate
        if not no_walks:
            sample_walks_files = glob.glob(os.path.join(color_csv_dir, "*.walks.csv"))
            walk_tag_jobs = []
            for walk_file in sample_walks_files:
                for gfa_file in viz_gfa_files:
                    gfa_base = os.path.basename(gfa_file)
                    sample = os.path.basename(walk_file).replace(".walks.csv", "")
                    colored_gfa = os.path.join(colored_gfa_dir, f"{sample}_{gfa_base}")
                    output_png = os.path.join(
                        bandage_img_dir, f"{sample}_{gfa_base.replace('.gfa', '.png')}"
                    )
                    if not up_to_date(colored_gfa, gfa_file, walk_file):
                        walk_tag_jobs.append((walk_file, gfa_file, colored_gfa, output_png, sample))

            if walk_tag_jobs:
                walk_memgate = _PlotMemoryGate(
                    min_free_gb=plot_min_free_gb, check_sec=plot_mem_check_sec,
                    label=f"{site_name}-walks"
                )
                W_walk = _effective_plot_workers(max_workers, plot_threads)

                def _tag_walk(job):
                    wf, gf, cg, _, _ = job
                    walk_memgate.wait_for_memory()
                    add_color_tags(wf, gf, cg)

                with ThreadPoolExecutor(max_workers=W_walk) as ex:
                    futs = [ex.submit(_tag_walk, j) for j in walk_tag_jobs]
                    for fut in as_completed(futs):
                        try:
                            fut.result()
                        except Exception as e:
                            print(f"[warn] walk tag failed: {e}")

            # collect bandage jobs after tagging
            for walk_file in sample_walks_files:
                for gfa_file in viz_gfa_files:
                    gfa_base = os.path.basename(gfa_file)
                    sample = os.path.basename(walk_file).replace(".walks.csv", "")
                    colored_gfa = os.path.join(colored_gfa_dir, f"{sample}_{gfa_base}")
                    output_png = os.path.join(
                        bandage_img_dir, f"{sample}_{gfa_base.replace('.gfa', '.png')}"
                    )
                    if no_images:
                        continue
                    if not up_to_date(output_png, colored_gfa):
                        bandage_jobs.append(
                            (colored_gfa, output_png, walk_file, f"{site_name} - Sample {sample} Walk", True)
                        )

        #if not no_images and bandage_jobs:
        #    W = _auto_workers(max_workers)
        #    with ThreadPoolExecutor(max_workers=W) as ex:
        #        futs = [ex.submit(_run_bandage, cg, png) for (cg, png, _, _, _) in bandage_jobs]
        if not no_images and bandage_jobs:
             W = _effective_plot_workers(max_workers, plot_threads)
             memgate = _PlotMemoryGate(min_free_gb=plot_min_free_gb, check_sec=plot_mem_check_sec, label=site_name)
             with ThreadPoolExecutor(max_workers=W) as ex:
                 futs = [ex.submit(_run_bandage_with_memwatch, cg, png, memgate) for (cg, png, _, _, _) in bandage_jobs]
                 for f in as_completed(futs):
                     f.result()
             for (cg, png, legend_src, title, is_walk) in bandage_jobs:
                 legend_out = os.path.join(legend_img_dir, os.path.basename(png))
                 create_legend_image(legend_src, png, legend_out, title, is_walk, colored_gfa=cg)

    print("Plotting done.")


# ------------------------------
# Summary
# ------------------------------

def generate_complex_sites_dataframe():
    complex_sites, genes_detected, pc_samples, pc_gene_counts, all_gene_counts = [], [], [], [], []
    total_genes = []

    def _genes_from_colored_gfa(gfa_path):
        genes = set()
        if not os.path.exists(gfa_path):
            return genes
        with open(gfa_path) as f:
            for line in f:
                if not line.startswith("S\t"):
                    continue
                parts = line.rstrip("\n").split("\t")
                for t in parts[3:] if len(parts) > 3 else []:
                    if t.startswith("GN:Z:"):
                        g = t.split(":", 2)[-1]
                        if g and not g.startswith("node_"):
                            genes.add(g)
                        break
        return genes

    def _genes_by_sample_from_gfas(colored_gfa_dir, gfa_base):
        all_counts, pc_counts = [], []
        pc_set = _genes_from_colored_gfa(os.path.join(colored_gfa_dir, f"pc_colors_{gfa_base}"))
        for gfa in sorted(glob.glob(os.path.join(colored_gfa_dir, f"*{gfa_base}"))):
            base = os.path.basename(gfa)
            if base.startswith(("all_colors_", "pc_colors_")):
                continue
            #sample = base.split("_", 1)[0]
            if base.endswith(gfa_base):
               sample = base[:-(len(gfa_base) + 1)]  # remove the "_" plus gfa_base
            else:
            # fallback (shouldn't happen, but safe)
               sample = base.split("_", 1)[0]
            genes = _genes_from_colored_gfa(gfa)
            if not genes:
                continue
            all_counts.append(f"{sample}:{len(genes)}")
            pc_counts.append(f"{sample}:{len(genes & pc_set)}" if pc_set else f"{sample}:0")
        return all_counts, pc_counts

    for site_dir in glob.glob("all_plottables/*"):
        site = os.path.basename(site_dir)
        colored_gfa_dir = os.path.join(site_dir, "colored_gfas")
        color_csv_dir = os.path.join(site_dir, "color_csvs", "gene_colors")

        if not os.path.exists(colored_gfa_dir):
            continue

        # find base from all_colors_*.gfa
        all_list = [f for f in os.listdir(colored_gfa_dir) if f.startswith("all_colors_") and f.endswith(".gfa")]
        if not all_list:
            # fallback when images/GFAs skipped: derive from gene_colors CSVs
            gene_files = glob.glob(os.path.join(color_csv_dir, "*.csv"))
            genes = [os.path.splitext(os.path.basename(f))[0] for f in gene_files]
            complex_sites.append(site)
            genes_detected.append(",".join(sorted(set(genes))))
            total_genes.append(len(set(genes)))
            pc_samples.append("")
            pc_gene_counts.append("")
            all_gene_counts.append("")
            continue

        gfa_base = all_list[0].replace("all_colors_", "")
        all_genes_set = _genes_from_colored_gfa(os.path.join(colored_gfa_dir, f"all_colors_{gfa_base}"))
        genes_detected.append(",".join(sorted(all_genes_set)))
        total_genes.append(len(all_genes_set))

        all_counts_list, pc_counts_list = _genes_by_sample_from_gfas(colored_gfa_dir, gfa_base)
        pc_samples_str = " ; ".join(sorted({s.split(":")[0] for s in pc_counts_list if not s.endswith(":0")}))
        pc_gene_counts.append(" ; ".join(pc_counts_list))
        all_gene_counts.append(" ; ".join(all_counts_list))
        complex_sites.append(site)
        pc_samples.append(pc_samples_str)

    return pd.DataFrame(
        {
            "Complex Site": complex_sites,
            "Genes Detected": genes_detected,
            "Total Genes Detected": total_genes,
            "Protein Coding Samples": pc_samples,
            "Protein Coding Gene Counts per Sample": pc_gene_counts,
            "All Gene Counts per Sample": all_gene_counts,
        }
    )


# ------------------------------
# CLI driver (example)
# ------------------------------

#def main():
#    p = argparse.ArgumentParser()
#    sub = p.add_subparsers(dest="cmd")

#    sp = sub.add_parser("complex")
#    sp.add_argument("vcf_file")
#    sp.add_argument("gfab")
#    sp.add_argument("--gff3",  required=True)
#    sp.add_argument("--ref_fasta")
#    sp.add_argument("--sep_pattern", default="#0#")
#    sp.add_argument("--ref_name", default="CHM13")
#    sp.add_argument("--plot_complex", action="store_true")
#    sp.add_argument("--no-walks", dest="no_walks", action="store_true", help="Do not plot sample walks when rendering plots")
#    sp.add_argument("--no_images", action="store_true", help=argparse.SUPPRESS)
#    sp.add_argument("--max_workers", type=int, default=0, help="Max parallel Bandage workers (0=auto)")
#    sp.add_argument("--gaf_file", default=None, help="Path to the gene alignments to the graph in GAF format")
#    sp.add_argument("--min_free_gb", type=float, default=20.0, help="Minimum free RAM (GB) before launching parallel workers")
#    sp.add_argument("--threads", "-t", dest="max_workers", type=int, default=0, help="Thread budget")
#    sp.add_argument("-a", type=int, default=5)
#    sp.add_argument("-n", type=int, default=1)
#    sp.add_argument("-s", type=int, default=10000)
#    sp.add_argument("-l", type=int, default=100000)
#    sp.add_argument("--sites", type=int, default=1)
#    sp.add_argument("--sv", type=int, default=1)
#    sp.add_argument("--plot_threads", type=int, default=0, help=("Bandage plotting workers. If set, overrides auto sizing. If not set, plotting uses 1/3 of --threads (min 1) instead of 1/3 o"))
#    sp.add_argument("--plot_min_free_gb", type=float, default=0.0, help=("Memory guard for Bandage: require at least this much free RAM (GB) before launching a plot. 0 disables the guard.") )
    # other args (a, n, s, l, sites, sv, etc.) omitted for brevity; keep your originals
#    sp.set_defaults(func=run_complex)

#    args = p.parse_args()
#    if not hasattr(args, "func"):
#        p.print_help(); sys.exit(1)
#    args.func(args)

def main():
    p = argparse.ArgumentParser()
    sub = p.add_subparsers(dest="cmd")
    add_subparser(sub)
    args = p.parse_args()
    if not hasattr(args, "func"):
        p.print_help()
        sys.exit(1)
    args.func(args)


#def run_complex(args):
#    # Build gene index once
#    _build_gene_index(args.gff3)#

    # Your pipeline to create tasks and call produce_plottable(*) goes here.
    # For brevity, we assume produce_plottable has already been run for each site
    # and per-site directories exist under all_plottables/.

#    if args.plot_complex:
#        plot_complex_sites(no_images=args.no_images, max_workers=args.threads)

#    df = generate_complex_sites_dataframe()
#    df.to_csv("complex_sites_summary.csv", index=False)
#    print("Saved complex_sites_summary.csv")

def run_complex(args):
    args = _normalize_args(args)

    # Backward compatibility: --no_images is deprecated; treat as --no-walks
    if getattr(args, 'no_images', False) and not getattr(args, 'no_walks', False):
        args.no_walks = True
        print('[WARN] --no_images is deprecated; use --no-walks instead.', file=sys.stderr)

    # 1) Build gene index once (fast region queries)
    _build_gene_index(args.gff3)

    # 2) RUN YOUR ORIGINAL PIPELINE to generate per-site outputs
    # If your file still has run_panscan_complex(...), call it here:
    if "run_panscan_complex" in globals():
        run_panscan_complex(
            vcf_file=args.vcf_file,
            gfab=args.gfab,                         # normalized above
            ref_fasta=getattr(args, "ref_fasta", None),
            gaf_file=getattr(args, "gaf_file", None),
            sep_pattern=args.sep_pattern,
            gff3_file=args.gff3,
            a=getattr(args, "a", 5),
            n=getattr(args, "n", 1),
            s=getattr(args, "s", 10000),
            regions=getattr(args, "regions", False),
            l=getattr(args, "l", 100000),
            sites=getattr(args, "sites", 1),
            sv=getattr(args, "sv", 1),
            ref_name=getattr(args, "ref_name", "CHM13"),
            no_walks=getattr(args, 'no_walks', False),
            max_workers=getattr(args, 'max_workers', 0),
            min_free_gb=getattr(args, 'plot_min_free_gb', 20.0),
        )
    else:
        print("[warn] run_panscan_complex(...) not found in this module.")
        print("       Summary/plots require per-site outputs under all_plottables/.")

    # 3) Plot / tag GFAs
    if getattr(args, "plot_complex", False):
        # render images + legends
        #plot_complex_sites(no_images=False, no_walks=getattr(args, 'no_walks', False), max_workers=getattr(args, "max_workers", 0))
        plot_complex_sites(no_images=False, no_walks=getattr(args, "no_walks", False), max_workers=getattr(args, "max_workers", 0), plot_threads=getattr(args, "plot_threads", 0), plot_min_free_gb=getattr(args, "plot_min_free_gb", 0.0), plot_mem_check_sec=getattr(args, "plot_mem_check_sec", 1.0))
    else:
        # default to no-image mode so legends/colored GFAs get produced for summary
        plot_complex_sites(no_images=True, no_walks=getattr(args, 'no_walks', False), max_workers=getattr(args, "max_workers", 0), plot_threads=getattr(args, "plot_threads", 0), plot_min_free_gb=getattr(args, "plot_min_free_gb", 0.0), plot_mem_check_sec=getattr(args, "plot_mem_check_sec", 1.0))

    # 4) Build summary
    df = generate_complex_sites_dataframe()
    df.to_csv("complex_sites_summary.csv", index=False)
    print("Saved complex_sites_summary.csv")




def add_subparser(subparsers):
    parser = subparsers.add_parser(
        "complex",
        help="Analyze complex loci from a VCF file."
    )
    parser.add_argument("vcf_file", help="Path to the VCF file generated from Minigraph-cactus pipeline.")
    parser.add_argument("gfab_file", help="Path to the gfab file generated from the GFA file")
    parser.add_argument("--ref_fasta", help="Path to the reference file")
    parser.add_argument("--gaf_file", help="Path to the gene alignments to the graph in GAF format")
    parser.add_argument("--sep_pattern", default="#0#", help="Separator for the sample and the sequence as present in GAF file (eg. '#0#')")
    parser.add_argument("--gff3", help="Path to the gff3 file")

    parser.add_argument("-a", type=int, default=5, help="Number of alleles to define a complex site (default: 5).")
    parser.add_argument("-n", type=int, default=1, help="Number of 10kb variants to define a complex site (default: 1).")
    parser.add_argument("-s", type=int, default=10000, help="Minimum size of a variant to define a complex site (default: 10000).")
    parser.add_argument("--regions", action="store_true", help="List complex regions.")
    parser.add_argument("-l", type=int, default=100000, help="Length of the region to define a complex region (default: 100000).")
    parser.add_argument("--sites", type=int, default=1, help="Number of complex sites in a region to define it as complex (default: 1).")
    parser.add_argument("--sv", type=int, default=1, help="Number of secondary SVs in a region to define it as complex (default: 1).")
    parser.add_argument("--ref_name", default=None, help="Reference name to use in region queries as present in GAF file (eg. CHM13) (default: None)")

    # Plot control
    parser.add_argument("--plot_complex", action="store_true", help=argparse.SUPPRESS)  
    parser.add_argument("--no-walks", dest="no_walks", action="store_true", help="Do not plot sample walks when rendering plots")

    # legacy alias (deprecated)
    parser.add_argument("--no_images", action="store_true", help=argparse.SUPPRESS)

    # Parallelism (alias both --threads/-t and --max_workers to same destination)
    #parser.add_argument("--threads", "-t", dest="max_workers", type=int, default=0,
    #                    help="Max parallel Bandage workers (0 = auto)")
    parser.add_argument("--threads", "-t", dest="max_workers", type=int, default=0, help=("Thread budget for this command. By default, Bandage plotting uses 1/3 of this value (min 1). " "Use --plot_threads to override plotting concurrency explicitly. (0 = auto)"))
    parser.add_argument("--max_workers", dest="max_workers", type=int, default=0,
                        help=argparse.SUPPRESS)  # keep for backward compat, hidden in help
    parser.add_argument("--plot_threads", type=int, default=0, help=("Bandage plotting workers. If set, overrides auto sizing. " "If not set, plotting uses 1/3 of --threads (min 1) instead of 1/3 of system CPUs."))
    parser.add_argument("--plot_min_free_gb", type=float, default=0.0, help=("Memory guard for Bandage: require at least this much free RAM (GB) before launching a plot. " "0 disables the guard.") )
    parser.add_argument("--min_free_gb", type=float, default=20.0, help="Minimum free RAM (GB) required before launching parallel region workers (default: 20.0).")
    parser.add_argument("--plot_mem_check_sec", type=float, default=1.0, help="How often to re-check free RAM when the memory guard is active (seconds).")
    parser.set_defaults(func=run_complex)



if __name__ == "__main__":
    main()
