#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use YAML::XS 'LoadFile';
use perlModules::panscan qw(cleanUp print_help1);
use Cwd;
my $currpath = getcwd();
use lib '$currpath';
use FindBin;
use File::Spec;

#BEGIN#######################################################
# Version:   1.0
# Filename:   preprocessVCF.pl
# Author: Dr. Bipin Balan
# Date: 15/04/2024
# Purpose:   Convert the multi allelic vcf records to single
#            allelic, followed by rtg decompose and then
#            perform genotype merging to generate the pre-
#            processed vcf file.
# External modules used:
# Notes:
#    1.) To execute :
#        perl preprocessVCF.pl sample.vcf
#    2.) Prerequisite : https://samtools.github.io/bcftools/bcftools.html
#                       https://github.com/RealTimeGenomics/rtg-tools
##############################################################

# Specify the path to the config file
my $script_dir = $FindBin::Bin;

# Construct the full path to the config.yaml file
my $config_file = File::Spec->catfile($script_dir, 'config.yaml');

# Check if the config file exists
if(!-e $config_file)
{
    die("Config file '$config_file' not found.\n");
}

# Create a new Config::IniFiles object
my $cfg = LoadFile($config_file);

# Check if the config file was successfully read
if(!$cfg) 
{
    die("Failed to read YAML file: $config_file\n");
}

# Retrieve the paths of the tools from the config file
my $bcftools = $cfg->{tools}->{bcftools};
my $rtg = $cfg->{tools}->{rtg};

# Temp folder #
my $tmp_folder = 'tmp';
unless(-d $tmp_folder)
{
	mkdir $tmp_folder;
}

my $infile = 'NA';
my $help = undef;
GetOptions(
        'i=s'   => \$infile,
        'help'      => \$help
) or die "Error in command line arguments\n";


# Print the help message #
if ($help)
{
    print_help1();
    exit;
}

if($infile eq 'NA')
{
        print "\nERROR : Please provide input vcf file !! \n\n";
        print_help1();
        exit;
}


my $sample = (split'\.',(split'\/',$infile)[-1],-1)[0];
my $resultfile1 = "$tmp_folder\/$sample\_single_allelic.vcf";
my $resultfile2 = "$tmp_folder\/$sample\_single_allelic_decomposed.vcf";
my $resultfile3 = "$sample\_preprocessed.vcf";

# Covert to single allelic #
my $status = 0;
$status = system("$bcftools norm --threads 64 -m- $infile -o $resultfile1");

if($status == 0)
{
	print "Successfully converted to singlealleic vcf file\n";
	$status = system("$rtg vcfdecompose --break-indels --break-mnps -Z -i $resultfile1 -o $resultfile2");
	
	if($status == 0)
	{
		print "Successfully decomposed\n";
		
		# Merge genotypes #
		&mergeGT($resultfile2,$resultfile3);
	}
	else
	{
		print "Error in rtg decompose\n";
		exit;
	}
}
else
{
	print "Error in converting to single allelic\n";
	exit;
}

# Cleaning temporary folder #
&cleanUp($tmp_folder);
exit;



#################### SUBROUTINES ###################################

sub mergeGT
{
	my $infile = shift;
	my $resultfile = shift;

	my @header = ();
	my @data = ();
	my @sampleindex = ();
	my %GenoTypeInfo = ();
	my $total = 0;
	my $merged = 0;

	open OUT,">$resultfile" or die "Can't open $resultfile for writing\n";
	open IN,"<$infile" or die "Can't open $infile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);

		if(/^\s*\#CHR.*/) # Extracting the samples from headers #
		{
			@header = split'\t',$_,-1;
		}

		if(/^\s*\#+.*/) # Skipping VCF headers #
		{
			print OUT "$_\n";
			next;
		}
		$total++;
		@data = split"\t",$_;
		my $chr = $data[0];
		my $startPos = $data[1];
		my $ref = $data[3];
		my $alt = $data[4];

		if($chr=~/CHM13/g) # formating the Chromosome number #
		{
			my ($Chmchr,$chrno) = split'\#',$chr,-1;
			$chr = $chrno;
		}
		
		# Printing previous variant if current line holds a new variant #
		if(! defined($GenoTypeInfo{"$chr\t$startPos\t$ref\t$alt"}))
		{
			foreach my $previosVar(keys %GenoTypeInfo)
			{
				$merged++;
				print OUT $GenoTypeInfo{$previosVar}{'INFO'}."\t";
				my @GTinfo = ();
				for(my $i=9;$i<=$#data;$i++)
				{
					if($i==9)
					{
						if($GenoTypeInfo{$previosVar}{$header[$i]}>1)
						{
							$GenoTypeInfo{$previosVar}{$header[$i]} = 1;
						}
					}
					my $gt = $GenoTypeInfo{$previosVar}{$header[$i]};
					push(@GTinfo,$gt);
				}
				my $gtinfo = join"\t",@GTinfo;
				print OUT "$gtinfo\n";
			}
			%GenoTypeInfo = (); # Re-initializing #
		}
		
		for(my $i=9;$i<=$#data;$i++)
		{
			if($i==9)
			{
				if($data[$i] eq '.')
				{
					$data[$i] = 0;
				}
				$GenoTypeInfo{"$chr\t$startPos\t$ref\t$alt"}{$header[$i]}+=$data[$i];
			}
			else
			{
				if(defined($GenoTypeInfo{"$chr\t$startPos\t$ref\t$alt"}{$header[$i]}))
				{
					my $preGT = $GenoTypeInfo{"$chr\t$startPos\t$ref\t$alt"}{$header[$i]};
					my ($phap1,$phap2) = split'\|',$preGT,-1;
					$phap1 = 0 if($phap1 eq '.');
					$phap2 = 0 if($phap2 eq '.');
					my($hap1,$hap2) = split'\|',$data[$i],-1;
					$hap1 = 0 if($hap1 eq '.');
					$hap2 = 0 if($hap2 eq '.');
					my $Hap1 = $hap1+$phap1;
					my $Hap2 = $hap2+$phap2;
					$Hap1=1 if($Hap1>1);
					$Hap2=1 if($Hap2>1);
					$data[$i] = "$Hap1|$Hap2";
				}
				$GenoTypeInfo{"$chr\t$startPos\t$ref\t$alt"}{$header[$i]}=$data[$i];
			}
			my $varinfo = join"\t",($chr,@data[1..8]);
			$GenoTypeInfo{"$chr\t$startPos\t$ref\t$alt"}{'INFO'} = $varinfo;
		}
	}
	foreach my $previosVar(keys %GenoTypeInfo)
	{
		$merged++;
		print OUT $GenoTypeInfo{$previosVar}{'INFO'}."\t";
		my @GTinfo = ();
		for(my $i=9;$i<=$#data;$i++)
		{
			if($i==9)
			{
				if($GenoTypeInfo{$previosVar}{$header[$i]}>1)
				{
					$GenoTypeInfo{$previosVar}{$header[$i]} = 1;
				}
			}
			my $gt = $GenoTypeInfo{$previosVar}{$header[$i]};
			push(@GTinfo,$gt);
		}
		my $gtinfo = join"\t",@GTinfo;
		print OUT "$gtinfo\n";
	}
	%GenoTypeInfo = (); # Re-initializing #

	print "Genotype merging completed successfully\n";
	print "Total Variants : $total\n";
	print "Total variants after Merging : $merged\n";

	close IN;
	close OUT;
}
