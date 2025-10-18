#!/usr/bin/perl -w
use strict;
use lib '.';
use Getopt::Long;
use perlModules::panscan qw(cleanUp print_uniq_var_help generateInput compareSNPDB compareINDELDB generateInputSV compareSV);

use Cwd 'getcwd', 'chdir';
my $orig_dir = getcwd();

#BEGIN#######################################################
# Version:   1.0
# Filename:  findUniqVariants.pl
# Author: Dr. Bipin Balan
# Date: 16/04/2024
# Purpose: Find uniq files in the given vcf file by compairing
#          with existing DBs.
#
# External modules used:
# Notes:
#    1.) To execute :
#        perl findUniqVariants.pl --i cohort.vcf --t threads --db SNP
#        execute perl findUniqVariants.pl to get complete arguments.
##############################################################
#       Modification Block
# Number #   Author    Change Description.      Date     Version
#
#########################################################END#

# Variable Declaration #
my $inputVcf = 'NA';
my $varType = 'NA';
my $db = 'ALL';
my $overlap = 80;
my $help = undef;
my $pbsv = 0;
my %types = ('SNP','SNP','INDEL','INDEL','SV','SV');
my $db_path = 'NA';

GetOptions(
	'i=s'   => \$inputVcf,
	't=s'   => \$varType,
	'pbsv'   => \$pbsv,
	'db=s'  => \$db,
	'db_path=s'  => \$db_path,
	'overlap=i'  => \$overlap,
	'help'      => \$help    
) or die "Error in command line arguments\n";


# Print the help message #
if ($help)
{
    print_uniq_var_help();
    exit;
}

$varType = uc($varType);
if($varType eq 'NA')
{
	print "\nERROR : Please provide variant type to compare(SNP/INDEL/SV) !! \n\n";
	print_uniq_var_help();
	exit;
}
elsif(!defined($types{$varType}))
{
	print "\nERROR : Unknown variant type <$varType>.Please provide variant type to compare(SNP/INDEL/SV) !! \n\n";
	print_uniq_var_help();
	exit;
}

if($inputVcf eq 'NA')
{
	print "ERROR: Please provide input VCF file !! \n";
	print_uniq_var_help();
	exit;
}

if($db_path eq 'NA')
{
	print "ERROR: Please provide the database folder path !! \n";
	print_uniq_var_help();
	exit;
}


# Result folders #
#my $tmpfolder = 'tmp_'.$varType;
#my $resultfolder = "$varType\_Comparison\_Result";

my $calldir = $ENV{'PANSCAN_CALLDIR'} // getcwd();

# Build absolute folders under the caller dir (instead of relative paths)
my $tmpfolder     = "$calldir/tmp_$varType";
my $resultfolder  = "$calldir/${varType}_Comparison_Result";


# generating resultfolders #
if(-d $tmpfolder)
{
	&cleanUp($tmpfolder);
	mkdir $tmpfolder
}
else
{
	mkdir $tmpfolder;
}

if(-d $resultfolder)
{
	my $prefix = (split'\_',$resultfolder,-1)[-1];
	my $fNo = 0;
	if($prefix eq 'RESULT')
	{
		$fNo=1;
	}
	elsif($prefix=~/^\d$/gi)
	{
		$fNo=$prefix+1;
	}
	$resultfolder = "$calldir/$varType\_Comparison\_Result\_$fNo";
	mkdir $resultfolder;
}
else
{
	mkdir $resultfolder;
}

# Database files #
my $dbSNP_SNP = "$db_path/DBSNP_SNP.txt";
my $dbSNP_INDEL = "$db_path/DBSNP_INDEL.txt";
my $dbSNP_SV_INS = "$db_path/DBSNP_SV_INS.txt";
my $dbSNP_SV_DEL = "$db_path/DBSNP_SV_DEL.txt";

my $gnomAD_SNP = "$db_path/GNOMAD_SNP.txt";
my $gnomAD_INDEL = "$db_path/GNOMAD_INDEL.txt";
my $gnomAD_SV_INS = "$db_path/GNOMAD_SV_INS.txt";
my $gnomAD_SV_DEL = "$db_path/GNOMAD_SV_DEL.txt";

my $Genomes_1000_SNP = "$db_path/1000Genome_SNP.txt";
my $Genomes_1000_INDEL = "$db_path/1000Genome_INDEL.txt";
my $Genomes_1000_SV_INS = "$db_path/1000Genome_SV_INS.txt";
my $Genomes_1000_SV_DEL = "$db_path/1000Genome_SV_DEL.txt";

my $HGSVC_SNP = "$db_path/HGSVC_SNP.txt";
my $HGSVC_INDEL = "$db_path/HGSVC_INDEL.txt";
my $HGSVC_SV_INS = "$db_path/HGSVC_SV_INS.txt";
my $HGSVC_SV_DEL = "$db_path/HGSVC_SV_DEL.txt";

my $GME_SNP = "$db_path/GME_SNP.txt";
my $GME_INDEL = "$db_path/GME_INDEL.txt";

my $DGV_SV_INS = "$db_path/DGV_SV_INS.txt";
my $DGV_SV_DEL = "$db_path/DGV_SV_DEL.txt";


if($varType eq 'SNP')
{
	# Calling subroutene to generate input files for SNP #
	my $inputFile = &generateInput($inputVcf,$tmpfolder,$varType,$resultfolder);

	my $sample = 'Input-cohort';
	my %DBs = ('dbSNP',$dbSNP_SNP,'gnomAD',$gnomAD_SNP,'1000Genomes',$Genomes_1000_SNP,'HGSVC',$HGSVC_SNP,'GME',$GME_SNP);
	my @db_order = ('dbSNP','gnomAD','1000Genomes','HGSVC','GME');
	my $totalDB = 0;
	if($db ne 'ALL')
	{
		@db_order = ();
		my %selectedDB = ();
		foreach my $eachInpDB(split',',$db)
		{
			if(defined($DBs{$eachInpDB}))
			{
				push(@db_order,$eachInpDB);
				$selectedDB{$eachInpDB} = $DBs{$eachInpDB};
			}
			else
			{
				print "Error : No database available for <$eachInpDB>\n";
				print "Please see the available $varType databases below : \n\n";
				print_uniq_var_help();
				exit;
			}
		}
		%DBs = %selectedDB;
		$totalDB = scalar(keys %DBs);
	}
	else
	{
		$totalDB = scalar(@db_order);
	}
	
	my $dbNo = 1;
	my $NovelSNPs = 0;
	foreach my $eachDB(@db_order)
	{
		&compareSNPDB(\$inputFile,$DBs{$eachDB},$sample,$eachDB,$resultfolder,$dbNo,$totalDB);
		$dbNo++;
	}
	print "The SNP comparison completed successfully and the results were generated in the folder : $resultfolder\n";
}
elsif($varType eq 'INDEL')
{
	# Calling subroutene to generate input files for INDEL #
	my $inputFile = &generateInput($inputVcf,$tmpfolder,$varType,$resultfolder);

	my $sample = 'Input-cohort';
        my %DBs = ('dbSNP',$dbSNP_INDEL,'gnomAD',$gnomAD_INDEL,'1000Genomes',$Genomes_1000_INDEL,'HGSVC',$HGSVC_INDEL,'GME',$GME_INDEL);
        my @db_order = ('dbSNP','gnomAD','1000Genomes','HGSVC','GME');
	my $totalDB = 0;
	if($db ne 'ALL')
	{
		@db_order = ();
		my %selectedDB = ();
		foreach my $eachInpDB(split',',$db)
		{
			if(defined($DBs{$eachInpDB}))
			{
				push(@db_order,$eachInpDB);
				$selectedDB{$eachInpDB} = $DBs{$eachInpDB};
			}
			else
			{
				print "Error : No database available for <$eachInpDB>\n";
				print "Please see the available $varType databases below : \n\n";
				print_uniq_var_help();
				&cleanUp($tmpfolder);
				&cleanUp($resultfolder);
				exit;
			}
		}
		%DBs = %selectedDB;
		$totalDB = scalar(keys %DBs);
	}
	else
	{
		$totalDB = scalar(@db_order);
	}
	
	my $dbNo = 1;
	foreach my $eachDB(@db_order)
	{
		&compareINDELDB(\$inputFile,$DBs{$eachDB},$sample,$eachDB,$resultfolder,$dbNo,$totalDB);
		$dbNo++;
	}
	print "The INDEL comparison completed successfully and the results were generated in the folder : $resultfolder\n";
}
elsif($varType eq 'SV')
{
	my %DB_Del = ('dbSNP',$dbSNP_SV_DEL,'gnomAD',$gnomAD_SV_DEL,'1000Genomes',$Genomes_1000_SV_DEL,'HGSVC',$HGSVC_SV_DEL,'DGV',$DGV_SV_DEL);
	my %DB_Ins = ('dbSNP',$dbSNP_SV_INS,'gnomAD',$gnomAD_SV_INS,'1000Genomes',$Genomes_1000_SV_INS,'HGSVC',$HGSVC_SV_INS,'DGV',$DGV_SV_INS);
	my @db_order = ('dbSNP','gnomAD','1000Genomes','HGSVC','DGV');
	my %selectedDB_INS = ();
	my %selectedDB_DEL = ();
	if($db ne 'ALL')
	{
		@db_order = ();
                foreach my $eachInpDB(split',',$db)
                {
			if(defined($DB_Del{$eachInpDB}))
			{
				push(@db_order,$eachInpDB);
				$selectedDB_DEL{$eachInpDB} = $DB_Del{$eachInpDB};
				$selectedDB_INS{$eachInpDB} = $DB_Ins{$eachInpDB};
			}
			else
			{
                                print "Error : No database available for <$eachInpDB>\n";
                                print "Please see the available $varType databases below : \n\n";
                                print_uniq_var_help();
				&cleanUp($tmpfolder);
				&cleanUp($resultfolder);
                                exit;
                        }
                }
		%DB_Del = %selectedDB_DEL;
		%DB_Ins = %selectedDB_INS;
	}

	# Calling subroutene to generate input files for SV #
	if($pbsv)
	{
		my ($svCount,$inputINS,$inputDEL,$inputDUP,$inputINV,$InsMap,$DelMap) = &generateInputSV($inputVcf,$tmpfolder,$resultfolder);
		if($svCount==0)
		{
			print "WARNING : There is no SVs present in the given vcf file. Pleas make sure the input is a pbsv or recommended SV vcf file\n";
		}
		else
		{
			if(-e $inputINS)
			{
				#&compareSV($inputINS,$tmpfolder,$resultfolder,$overlap,\@db_order,\%DB_Ins,$cfg,'Insertions',$InsMap);
				&compareSV($inputINS,$tmpfolder,$resultfolder,$overlap,\@db_order,\%DB_Ins,$db_path,'Insertions',$InsMap);
			}

			if(-e $inputDEL)
			{
				#&compareSV($inputDEL,$tmpfolder,$resultfolder,$overlap,\@db_order,\%DB_Del,$cfg,'Deletions',$DelMap);
				&compareSV($inputDEL,$tmpfolder,$resultfolder,$overlap,\@db_order,\%DB_Del,$db_path,'Deletions',$DelMap);
			}
		}
	}
	else
	{
		my ($svCount,$inputINS,$inputDEL,$InsMap,$DelMap) = &generateInput($inputVcf,$tmpfolder,$varType,$resultfolder);
		if($svCount==0)
		{
			print "WARNING : There is no SVs present in the given vcf file. Pleas make sure the input is a pbsv or recommended SV vcf file\n";
		}
		else
		{
			if(-e $inputINS)
			{
				#&compareSV($inputINS,$tmpfolder,$resultfolder,$overlap,\@db_order,\%DB_Ins,$cfg,'Insertions',$InsMap);
				&compareSV($inputINS,$tmpfolder,$resultfolder,$overlap,\@db_order,\%DB_Ins,$db_path,'Insertions',$InsMap);
			}

			if(-e $inputDEL)
			{
				#&compareSV($inputDEL,$tmpfolder,$resultfolder,$overlap,\@db_order,\%DB_Del,$cfg,'Deletions',$DelMap);
				&compareSV($inputDEL,$tmpfolder,$resultfolder,$overlap,\@db_order,\%DB_Del,$db_path,'Deletions',$DelMap);
			}
		}
	}
}

# Cleaning temporary folder #
&cleanUp($tmpfolder);

chdir($orig_dir) or warn "Could not return to original directory: $!";

exit;
