# PanScan
PanScan is a computational toolkit for the exploration of **novel sequences, structural variants (SVs), and repetitive regions** in genomic data.  
It supports **gene duplication annotation**, **complex region visualization**, and provides modules for **variant and sequence comparison**.

## 🚀 Quick Start (Containerized Installation)
The easiest way to use **PanScan** is with its **Docker** or **Singularity** image, both preloaded with all dependencies.

### **Option 1: Docker**
Pull the latest version:
```bash
docker pull catg/panscan:latest
```

Run:
```bash
docker run --rm -it catg/panscan:latest panscan
```

### **Option 2: Singularity**
Pull the image:
```bash
singularity pull docker://catg/panscan:latest
```

Run:
```bash
singularity exec panscan_latest.sif panscan
```

You should see:
```
usage: panscan [-h] {complex,preprocess_vcf,novel_seq,find_uniq_variants,make_dup_mtx,gene_dup} ...
```

> 💡 *Using Docker or Singularity removes the need to install Perl, Python, R, or other dependencies manually.*

## 🧰 Installation from Source
If you prefer to install from source:

```bash
git clone https://github.com/CATG-Github/panscan.git
cd panscan
pip install .
```

Run with:
```bash
panscan
```

## 📦 Requirements (for manual installations)
**Programming Languages**
- [Perl](https://www.perl.org/)
- [Python ≥ 3.7](https://www.python.org/)
- [R](https://www.r-project.org/)

**Python Packages**
`argparse, Image, ImageDraw, ImageFont, matplotlib, numpy, os, pandas, pickle, sys`

**Perl Modules**
`Getopt::Long, YAML::XS, Cwd, File::Copy, Exporter`

**Bioinformatics Tools**
- [Bandage](https://rrwick.github.io/Bandage)
- [BCFtools](https://github.com/samtools/bcftools)
- [bgzip](https://www.htslib.org/doc/bgzip.html)
- [cd-hit](https://github.com/weizhongli/cdhit)
- [GFABase](https://github.com/mlin/gfabase)
- [GraphAligner](https://github.com/maickrau/GraphAligner)
- [Liftoff](https://github.com/agshumate/Liftoff)
- [Panscan Zenodo Databases](https://zenodo.org/records/15314528)  
  BED files derived from [dbSNP](https://www.ncbi.nlm.nih.gov/snp), [1000 Genomes](https://www.internationalgenome.org/home), [gnomAD](https://gnomad.broadinstitute.org), [GME](https://illumina.github.io/NirvanaDocumentation/data-sources/gme), and [DGV](https://dgv.tcag.ca/dgv/app/home)
- [rtg-tools](https://github.com/RealTimeGenomics/rtg-tools)
- [tabix](https://www.htslib.org/doc/tabix.html)
- [truvari](https://github.com/ACEnglish/truvari)
- [RIdeogram](https://cran.r-project.org/web/packages/RIdeogram/vignettes/RIdeogram.html)



## 📄 Citation
If you use **PanScan** in your work, please cite the Zenodo DOI:  
👉 [https://zenodo.org/records/15314528](https://zenodo.org/records/15314528)

## Gene-duplication analyses

There are 2 parts to this:
You have to first run 
```
panscan make_dup_mtx --gencode_gff3 gencode.gff3 --fofn assemblies.fofn --ref_fa hg38_ref.fa --threads 64
```
to produce the gene-duplication matrix from all your assmeblies.

**Please Note: The gene duplication modules are only compatible with gencode gff3 files**

A sample assembly fofn has been uploaded with the package, it has to follow that format for panscan to use your assemblies. The format is as follows:

```
/path/to/assembly.fa SAMPLE_ID HAPLOTYPE
/path/to/assembly.fa SAMPLE_ID HAPLOTYPE
/path/to/assembly.fa SAMPLE_ID HAPLOTYPE
```
Then the second command 
```
panscan gene_dup --csv_file gene-dup-matrix.csv --ip1 gene_dup/hprc-matrix.csv --ip2 gene_dup/cpc-matrix.csv 
``` 
takes in your gene duplication matrix, and visualizes the duplications in your data and compares them with the hprc and cpc duplications as well. The plots made are :
 - Duplications per assembly
 - Venn diagram of duplications w.r.t HPRC and CPC
 - Frequency comparison of your duplications w.r.t. HPRC and CPC . (Plots the most distinct ones)

**As the HPRC and CPC duplciation matrices have only been shared in HG38, we do not recommend using the comparison plots if you use another reference**

The ```panscan gene_dup``` command will ask you for paths to HPRC and CPC matrices, they are present in the directory you cloned the repo into. Specifically in ```panscan/gene_dup```

To plot an ideogram showing gene duplications use the flag '''--ideogram''' with the command. It also needs the ideo-ip.csv

## Complex region analyses 
### Complex Sites
The complex sites (not to be interchanged with complex region) can be anlayzed from the vcf file generated from [Minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pipeline. 
The complex sites, as defined in the draft [Arab Pangenome reference](https://www.biorxiv.org/content/10.1101/2024.07.09.602638v1), are the sites with more than 5 haplotypes (alleles) and atleast one 10kb structural variant (SV). 

You can use  ```-a```, ```-n```, ```-s``` to redefine the parameters for a complex site. 

### Complex regions
Complex regions are regions of 100Kb with atleast one complex site and another SV. To list complex regions use the --regions flag with
```panscan complex``` subcommand
```-l``` can be used to define the length of the region
```--sites``` can be used to define how many sites should be present
```--sv``` can be used to define how many secondary SVs should be present

for a region to be considered a complex region.

### End-to-end
Run the command below in full to produce complex regions and haplotype walks for each sample in all the regions. 

```
panscan complex --ref_fasta chm13v2.0.fa --gaf_file chm13_mapped_genes.gaf --sep_pattern '#0#' --gff3 chm13v2.0_RefSeq_Liftoff_v5.1.gff3 -a 5 -n 1 -s 10000 --regions -l 100000 --sites 1 --sv 1 --ref_name CHM13 panscan.vcf panscan.gfab
```

### Plotting
To plot the complex regions (all sample walks and genes present in detected complex sequences) that were produced using the above command, the ```--plot_complex``` option can be used.
```
panscan complex --ref_fasta chm13v2.0.fa --gaf_file chm13_mapped_genes.gaf --sep_pattern '#0#' --gff3 chm13v2.0_RefSeq_Liftoff_v5.1.gff3 -a 5 -n 1 -s 10000 --regions -l 100000 --sites 1 --sv 1 --ref_name CHM13 --plot_complex panscan.vcf panscan.gfab 
```

**The gaf files needed for the complex command should be produced by aligning the gene sequences file to your pangenome.** 
The gene sequence files (and scripts to produce them for other references)  are present in the **complex.tar.gz** file present in the [Panscan Complex Region Files](https://drive.google.com/drive/folders/16O6InjctvIsGSTzroDu2366_wMrTFR3p).



## Pangenome VCF Processing

For all modules in this section you can use the Sample.vcf file provided in the repo for testing

### Novel Variants
This module identifies novel variants (SNPs, InDels, and SVs) in the input Pangenome VCF file by comparing them against public databases like dbSNP, gnomAD, 1000 Genomes, GME, and DGV.

```
panscan find_uniq_variants -i /path/to//panscan/Sample.vcf -t SNP --db ALL --overlap 80 --db-path /path/to/databases --op /path/to/output --debug
```

**This function needs the path for the databases to be provided to it**. These databases are present in the **database.tar.gz** present in the [Panscan Zenodo Database](https://zenodo.org/records/15314528).

### Novel seq
This module identifies novel sequences present in Pangenome VCF file1 by comparing SV insertions with those in Pangenome VCF file2 and reports them in FASTA format. Initially, the input VCF files undergo pre-processing, which involves splitting multi-allelic variants into single-allelic ones and decomposing complex variants into indels and SNPs using the "decompose" program from RTG Tools. After pre-processing, the SV insertions in the VCF files are compared using the "truvari bench" command, identifying novel SV insertions in VCF file1. These novel insertions at the same locus are clustered using the CD-HIT program, and the final novel sequence FASTA file is generated.

```
panscan novel_seq -i /path/to/panscan/Sample.vcf -r /path/to/panscan/Reference.vcf --genome HG38 -t 4 -dt 2 --dpi 600 --op /path/to/my/output --debug
```
As the pre-processing of the HPRC and CPC pangenomes can be very time-consuming, we recommend using the pre-processed VCF files for these references, available on [Zenodo](https://zenodo.org/records/15314528).

```
panscan novel_seq -i VCF -pRef VCF
```

The default number of threads for RTG decompose is set to 1 (recommended). If you choose to increase the thread count, ensure that your system has sufficient memory to handle the workload.
 
```
-dt DT, --dt DT       Number of threads for RTG decompose (default: 1).
```


