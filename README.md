
# Panscan
Pangenome analyses toolkit.

## Requirements
- Python >= 3.7 
- perl
- Matplotlib
- pandas
- [liftoff](https://github.com/agshumate/Liftoff) 
- [cd-hit](https://github.com/weizhongli/cdhit)
- [rtg-tools](https://github.com/RealTimeGenomics/rtg-tools)
- [truvari](https://github.com/ACEnglish/truvari)
- [GFABase](https://github.com/mlin/gfabase)
- [GraphAligner](https://github.com/maickrau/GraphAligner)
- [Panscan Databases](https://drive.google.com/drive/folders/16O6InjctvIsGSTzroDu2366_wMrTFR3p)


## Installation
To install Panscan, run:

```
git clone https://github.com/muddinmbru/panscan.git
cd panscan
pip install .
```

After successful installation run the tool with :

```
panscan
```


## Complex loci analyses 
The complex loci can be anlayzed from the vcf file generated from [Minigraph-cactus](https://github.com/ComparativeGenomicsToolkit/cactus/blob/master/doc/pangenome.md) pipeline. 
The complex sites, as defined in the draft [Arab Pangenome reference](https://www.biorxiv.org/content/10.1101/2024.07.09.602638v1), are the sites with more than 5 haplotypes (alleles) and atleast one 10kb structural variant (SV). To list all of the complex sites within this definition
```
panscan complex --regions -a 5 -n1 -s 10000
```

use ```-a```, ```-n```, ```-s``` to redefine the parameters for a complex site. 

### Complex regions
Complex regions are regions of 100Kb with atleast one complex site and another SV. To list complex regions use the --regions flag with
```panscan complex``` subcommand
```-l``` can be used to define the length of the region
```--sites``` can be used to define how many sites should be present
```--sv``` can be used to define how many secondary SVs should be present

for a region to be considered complex.



### End-to-end
Run the command below in full to produce complex regions and haplotype walks for each sample in all the regions. 

```
panscan complex --ref_fasta chm13v2.0.fa --gaf_file chm13_mapped_genes.gaf --sep_pattern '#0#' --gff3 chm13v2.0_RefSeq_Liftoff_v5.1.gff3 -a 5 -n 1 -s 10000 --regions -l 100000 --sites 1 --sv 1 --ref_name CHM13 panscan.vcf panscan.gfab
```

**The gaf files needed for the complex command should be produced by aligning the gene sequences file to your pangenome.** 
The gene sequence files (and scripts to produce them for other references)  are present in the **complex.tar.gz** file present in the [Panscan Database](https://drive.google.com/drive/folders/16O6InjctvIsGSTzroDu2366_wMrTFR3p).


## Gene-duplication analyses

There are 2 parts to this:
You have to first run 
```
panscan make_dup_mtx gencode.gff3 assemblies.fofn hg38_ref.fa 64
```
to produce the gene-duplication matrix from all your assmeblies.

**Please Note: The gene duplication modules are only compatible with gencode gff3 files**

Then the second command 
```
panscan gene_dup gene-dup-matrix.csv gene_dup/hprc-matrix.csv gene_dup/cpc-matrix.csv
``` 
takes in your gene duplication matrix, and visualizes the duplications in your data and compares them with the hprc and cpc duplications as well. The plots made are :
 - Duplications per assembly
 - Venn diagram of duplications w.r.t HPRC and CPC
 - Frequency comparison of your duplications w.r.t. HPRC and CPC . (Plots the most distinct ones)

**As the HPRC and CPC duplciation matrices have only been shared in HG38, we do not recommend using the comparison plots if you use another reference**

The ```panscan gene_dup``` command will ask you for paths to HPRC and CPC matrices, they are present in the directory you cloned the repo into. Specifically in ```panscan/gene_dup```

## Pangenome VCF Processing

For all modules in this section you can use the Sample.vcf file provided in the repo for testing

### Preprocess Vcf

This module will convert multi-allelic Pangenome VCF records to single-allelic ones. Next, complex indels will be decomposed into SNPs and indels using the RTG tools "decompose" program. Finally, the genotypes of variants at the same locus will be merged to produce the final pre-processed VCF file.

```
panscan preprocess_vcf --i Sample.vcf
```

### Novel seq
This module identifies novel sequences present in Pangenome VCF file1 by comparing SV insertions with those in Pangenome VCF file2 and reports them in FASTA format. Initially, the input VCF files undergo pre-processing, which involves splitting multi-allelic variants into single-allelic ones and decomposing complex variants into indels and SNPs using the "decompose" program from RTG Tools. After pre-processing, the SV insertions in the VCF files are compared using the "truvari bench" command, identifying novel SV insertions in VCF file1. These novel insertions at the same locus are clustered using the CD-HIT program, and the final novel sequence FASTA file is generated.

```
panscan novel_seq APR.vcf CPC-HPRC.vcf 
```

The provided Pangenome VCF files of the APR and CPC-HPRC can be used to be compared with your VCF files as well.

### Novel Variants
This module identifies novel variants (SNPs, InDels, and SVs) in the input Pangenome VCF file by comparing them against public databases like dbSNP, gnomAD, 1000 Genomes, GME, and DGV.

```
panscan find_uniq_variants --i Sample.vcf --t SNP --db ALL --op 80 --output novel_variants  --db-path downloads/database
```

**This function needs the path for the databases to be provided to it**. These databases are present in the **database.tar.gz** present in the [Panscan Database](https://drive.google.com/drive/folders/16O6InjctvIsGSTzroDu2366_wMrTFR3p).
