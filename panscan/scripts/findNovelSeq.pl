#!/usr/bin/perl -w
use strict;

#BEGIN#######################################################
# Version:   1.0
# Filename:   findNovelSeq.pl
# Author: Dr. Bipin Balan
# Date: 15/04/2024
# Purpose: Generate a novel sequence FASTA file from VCF file1 by
#          comparing it with VCF file2.
# External modules used:
# Notes:
#    1.) To execute :
#        perl findNovelSeq.pl sample1.vcf
#    2.) Prerequisite : https://samtools.github.io/bcftools/bcftools.html
#                       https://github.com/RealTimeGenomics/rtg-tools
#			https://www.htslib.org/doc/bgzip.html
#			https://www.htslib.org/doc/tabix.html
#			https://github.com/ACEnglish/truvari
#			https://sites.google.com/view/cd-hit
##############################################################

#!/usr/bin/perl -w
#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use File::Copy;
use YAML::XS 'LoadFile';
use Getopt::Long;
use Cwd qw(getcwd abs_path);
use perlModules::panscan qw(preprocessVCF mergeGT cleanUp print_novel_seq_help extractSVIns extractNovelInsertions clusterNovelIns getSVsampleInfo generateIdeogram runTruvari chrSplitVcf mergeTruvariVcfs);

# Specify the path to the config file relative to the script directory.
my $script_dir = $FindBin::Bin;
my $config_file = "$script_dir/config.yaml";

if (! -e $config_file) {
    die("Config file '$config_file' not found.\n");
}

my $cfg = LoadFile($config_file);
if (!$cfg) {
    die("Failed to read YAML file: $config_file\n");
}

# Retrieve tool paths from config file.
my $truvari   = $cfg->{tools}->{truvari};
my $bgzip     = $cfg->{tools}->{bgzip};
my $tabix     = $cfg->{tools}->{tabix};
my $cdhit_est = $cfg->{tools}->{cdhit_est};
my $Rscript   = $cfg->{tools}->{Rscript};

# Variable declarations.
my $infile1 = undef;
my $infile2 = undef;
my $pInfile1 = undef;
my $pInfile2 = undef;
my $threads = 1;
my $dthreads = 1;
my $excludeSample = 'NA';
my $dpi = 600;
my $genome = 'HG38';
my $help = undef;
my $output = '';   # New output directory flag.

GetOptions(
    'i=s'       => \$infile1,
    'r=s'       => \$infile2,
    't=i'       => \$threads,
    'dt=i'      => \$dthreads,
    'dpi=i'     => \$dpi,
    'pInp=s'    => \$pInfile1,
    'pRef=s'    => \$pInfile2,
    'exclude=s' => \$excludeSample,
    'genome=s'  => \$genome,
    'op=s'      => \$output,  # New output option.
    'help'      => \$help
) or die "Error in command line arguments\n";

# Determine the result folder.
my $resultfolder;
if ($output) {
    $resultfolder = abs_path($output) or die "Cannot resolve output directory: $output\n";
    mkdir $resultfolder unless (-d $resultfolder);
} else {
    $resultfolder = 'NovelSeq_Results';
    mkdir $resultfolder unless (-d $resultfolder);
}

my $karyo = "$cfg->{databases}->{dbpath}/HG38_karyotype.bed";
if ($genome eq 'CHM13') {
    $karyo = "$cfg->{databases}->{dbpath}/CHM13_karyotype.bed";
}

if ($dthreads > 1) {
    print "WARNING !! Please make sure the system/server has enough memory to run RTG decompose !!\n";
}

if ($help) {
    print_novel_seq_help();
    exit;
}

if((!defined($infile1)) && (!defined($pInfile1))) {
    print "\nERROR: Please provide input pangenome VCF file !!\n\n";
    print_novel_seq_help();
    exit;
}

if((!defined($infile2)) && (!defined($pInfile2))) {
    print "\nERROR: Please provide reference pangenome VCF file to compare !!\n\n";
    print_novel_seq_help();
    exit;
}

my $tmpFolder = 'tmp_novel_seq';
if (-d $tmpFolder) {
    cleanUp($tmpFolder);
    mkdir $tmpFolder;
} else {
    mkdir $tmpFolder;
}

# Pre-processing of the VCF files.
my $ppvcf1 = '';
my $ppvcf2 = '';
my %ExcludeSamples = ();
if ($excludeSample ne 'NA') {
    foreach my $eachSamples (split /,/, $excludeSample, -1) {
        $ExcludeSamples{$eachSamples} = 0;
    }
}
if (defined($pInfile1)) {
    my $sample = (split /_/, (split /\./, (split /\//, $pInfile1)[-1], -1)[0], -1)[0];
    my $resultfile = "$tmpFolder/${sample}_preprocessed.vcf";
    $ppvcf1 = mergeGT($pInfile1, $resultfile, \%ExcludeSamples);
} else {
    $ppvcf1 = preprocessVCF($cfg, $infile1, $threads, $tmpFolder, $dthreads, \%ExcludeSamples);
}
if (defined($pInfile2)) {
    my $sample = (split /_/, (split /\./, (split /\//, $pInfile2)[-1], -1)[0], -1)[0];
    my $resultfile = "$tmpFolder/${sample}_preprocessed.vcf";
    $ppvcf2 = mergeGT($pInfile2, $resultfile, \%ExcludeSamples);
} else {
    $ppvcf2 = preprocessVCF($cfg, $infile2, $threads, $tmpFolder, $dthreads, \%ExcludeSamples);
}

# Extract SV Insertions.
my $InsSV1 = extractSVIns($ppvcf1, $tmpFolder, $bgzip, $tabix);
my $NovelSampleFile = getSVsampleInfo($InsSV1, $tmpFolder);
my $InsSV2 = extractSVIns($ppvcf2, $tmpFolder, $bgzip, $tabix);

# Comparison of the VCFs using Truvari bench command.
runTruvari($InsSV1, $InsSV2, $tmpFolder, $threads, $truvari, $bgzip, $tabix);
extractNovelInsertions($tmpFolder, $resultfolder, $NovelSampleFile, $karyo, $Rscript, $dpi, $threads, $cdhit_est);

cleanUp($tmpFolder);
exit;
