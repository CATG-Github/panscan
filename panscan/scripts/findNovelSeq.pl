#!/usr/bin/perl -w
use strict;
use File::Copy;
use lib '.';
use FindBin qw($Bin);
use lib "$Bin";  
use lib "$Bin/perlModules";
use Getopt::Long;

use perlModules::panscan qw(preprocessVCF mergeGT cleanUp print_novel_seq_help extractSVIns extractNovelInsertions clusterNovelIns getSVsampleInfo generateIdeogram runTruvari chrSplitVcf mergeTruvariVcfs);

#BEGIN#######################################################
# Version:   1.0
# Filename:   findNovelSeq.pl
# Author: Dr. Bipin Balan
# Date: 15/04/2024
# Purpose: Generate a novel sequence FASTA file from VCF file1 by
#          comparing it with VCF file2.
# External modules used:
# Notes:
#    To execute :
#        perl findNovelSeq.pl ( provides the complete arguments )
##############################################################
#       Modification Block
# Number #   Author    Change Description.      Date     Version
# 1.    Dr.Bipin Balan Implemented threading   9/04/2025 1.0
#
#########################################################END#


# Retrieve the paths of the tools from the config file
my $truvari = 'truvari';
my $bgzip = 'bgzip';
my $tabix = 'tabix';
my $cdhit_est = 'cdhit_est';
my $Rscript = 'Rscript';

# Variable Declaration #
my $infile1 = undef; # Input vcf file #
my $infile2 = undef; # Input refernce vcf file #
my $pInfile1 = undef; # Input pre-processed vcf file #
my $pInfile2 = undef; # Input pre-processed refernce vcf file #
my $threads = 1; # Number of Threads #
my $dthreads = 1; # Number of Threads for RTG decompose #
my $db_path = 'NA';
my $excludeSample = 'NA';
my $dpi = 600; # dpi for Ideogram figure #
my $genome = 'NA';
my $help = undef;

GetOptions(
        'i=s'   => \$infile1,
        'r=s'   => \$infile2,
        't=i'  => \$threads,
        'dt=i'  => \$dthreads,
        'db_path=s'  => \$db_path,
        'dpi=i'  => \$dpi,
	'pInp=s' => \$pInfile1,
	'pRef=s' => \$pInfile2,
	'exclude=s' => \$excludeSample,
        'help'      => \$help
) or die "Error in command line arguments\n";


# Print warning message #
if($dthreads>1)
{
	print "WARNING !! Please make sure the system/server is having enogh memory to run RTG decompose !!\n";
}


# Print the help message #
if ($help)
{
    print_novel_seq_help();
    exit;
}

if($db_path eq 'NA')
{
        print "ERROR: Please mention the pangenome database path !! \n";
        print_novel_seq_help();
        exit;
}


if( (!(defined($infile1))) && (!(defined($pInfile1))) )
{
        print "\nERROR : Please provide input pangenome vcf file !! \n\n";
        print_novel_seq_help();
        exit;
}

if( (!(defined($infile2))) && (!(defined($pInfile2))) )
{
        print "\nERROR : Please provide refernce pangenome vcf file to compare !! \n\n";
        print_novel_seq_help();
        exit;
}

my $karyo = "$db_path\/karyotype.bed";

my $tmpFolder = 'tmp_novel_seq';
if(-d $tmpFolder)
{
	&cleanUp($tmpFolder);
	mkdir $tmpFolder;
}
else
{
	mkdir $tmpFolder;
}

my $resultfolder = 'NovelSeq_Results';
if(-d $resultfolder)
{
        my $prefix = (split'\_',$resultfolder,-1)[-1];
        my $fNo = 0;
        if($prefix eq 'Results')
        {
                $fNo=1;
        }
        elsif($prefix=~/^\d$/gi)
        {
                $fNo=$prefix+1;
        }
        $resultfolder = "NovelSeq_Results\_$fNo";
        mkdir $resultfolder;
}
else
{
        mkdir $resultfolder;
}

# Pre-processing of the VCF files #
my $ppvcf1 = '';
my $ppvcf2 = '';
my %ExcludeSamples = ();
if($excludeSample ne 'NA')
{
	foreach my $eachSamples(split'\,',$excludeSample,-1)
	{
		$ExcludeSamples{$eachSamples} = 0;
	}
}
if(defined($pInfile1))
{
	my $sample = (split'\_',(split'\.',(split'\/',$pInfile1)[-1],-1)[0],-1)[0];
	my $resultfile = "$tmpFolder\/$sample\_preprocessed.vcf";
	$ppvcf1 = &mergeGT($pInfile1,$resultfile,\%ExcludeSamples);
}
else
{
	$ppvcf1 = &preprocessVCF($infile1,$threads,$tmpFolder,$dthreads,\%ExcludeSamples);
}
if(defined($pInfile2))
{
	my $sample = (split'\_',(split'\.',(split'\/',$pInfile2)[-1],-1)[0],-1)[0];
	my $resultfile = "$tmpFolder\/$sample\_preprocessed.vcf";
	$ppvcf2 = &mergeGT($pInfile2,$resultfile,\%ExcludeSamples);
}
else
{
	$ppvcf2 = &preprocessVCF($infile2,$threads,$tmpFolder,$dthreads,\%ExcludeSamples);
}

# Extract SV Insertions #
my $InsSV1 = &extractSVIns($ppvcf1,$tmpFolder,$bgzip,$tabix);
my $NovelSampleFile = &getSVsampleInfo($InsSV1,$tmpFolder,\%ExcludeSamples);
my $InsSV2 = &extractSVIns($ppvcf2,$tmpFolder,$bgzip,$tabix);

# Comparison of the VCFs using Truvari bench command #
&runTruvari($InsSV1,$InsSV2,$tmpFolder,$threads,$truvari,$bgzip,$tabix);
&extractNovelInsertions($tmpFolder,$resultfolder,$NovelSampleFile,$karyo,$Rscript,$dpi,$threads,$cdhit_est,\%ExcludeSamples);

# Removing temp folders #
&cleanUp($tmpFolder);
exit;
