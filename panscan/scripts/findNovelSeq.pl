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
use strict;
use Cwd qw(getcwd abs_path);
use File::Copy;
use lib '.';
use YAML::XS 'LoadFile';
use Getopt::Long;
use perlModules::panscan qw(preprocessVCF mergeGT cleanUp print_help1);

# Specify the path to the config file
my $config_file = "config.yaml";

# Check if the config file exists
if(!-e $config_file) {
    die("Config file '$config_file' not found.\n");
}

# Load configuration from config.yaml
my $cfg = LoadFile($config_file);
if(!$cfg) {
    die("Failed to read YAML file: $config_file\n");
}

# Retrieve tool paths from config file
my $truvari    = $cfg->{tools}->{truvari};
my $bcftools   = $cfg->{tools}->{bcftools};
my $rtg        = $cfg->{tools}->{rtg};
my $bgzip      = $cfg->{tools}->{bgzip};
my $tabix      = $cfg->{tools}->{tabix};
my $cdhit_est  = $cfg->{tools}->{cdhit_est};
my $Rscript    = $cfg->{tools}->{Rscript};

# Variable Declaration
my $infile1 = undef;    # Input VCF file (pangenome)
my $infile2 = undef;    # Reference VCF file
my $pInfile1 = undef;   # Preprocessed input VCF file (optional)
my $pInfile2 = undef;   # Preprocessed reference VCF file (optional)
my $threads = 1;        # Number of threads
my $excludeSample = 'NA';
my $dpi = 600;          # DPI for ideogram figure
my $genome = 'HG38';
my $help = undef;
my $output = '';        # New output directory flag
# New flags can be added here if needed

GetOptions(
    'i=s'      => \$infile1,
    'r=s'      => \$infile2,
    't=s'      => \$threads,
    'dpi=i'    => \$dpi,
    'pInp=s'   => \$pInfile1,
    'pRef=s'   => \$pInfile2,
    'exclude=s'=> \$excludeSample,
    'genome=s' => \$genome,
    'op=s'     => \$output,  # New output directory flag
    'help'     => \$help
) or die "Error in command line arguments\n";

# Print help if requested
if ($help) {
    print_help1();
    exit;
}

my $karyo = "$cfg->{databases}->{dbpath}/HG38_karyotype.bed";
if($genome eq 'CHM13') {
    $karyo = "$cfg->{databases}->{dbpath}/CHM13_karyotype.bed";
}

# Determine result folder based on -op flag
my $resultfolder;
if($output) {
    $resultfolder = abs_path($output) or die "Cannot resolve output directory: $output\n";
    mkdir $resultfolder unless(-d $resultfolder);
} else {
    $resultfolder = 'NovelSeq_Results';
    mkdir $resultfolder unless(-d $resultfolder);
}

# Create temporary folder
my $tmpFolder = 'tmp';
if(-d $tmpFolder) {
    cleanUp($tmpFolder);
    mkdir $tmpFolder;
} else {
    mkdir $tmpFolder;
}

# Pre-processing of the VCF files
my $ppvcf1 = '';
my $ppvcf2 = '';
my %ExcludeSamples = ();
if($excludeSample ne 'NA') {
    foreach my $eachSamples(split /,/, $excludeSample, -1) {
        $ExcludeSamples{$eachSamples} = 0;
    }
}
if(defined($pInfile1)) {
    my $sample = (split /_/, (split /\./, (split /\//, $pInfile1)[-1], -1)[0], -1)[0];
    my $resultfile = "$tmpFolder/${sample}_preprocessed.vcf";
    $ppvcf1 = mergeGT($pInfile1, $resultfile, \%ExcludeSamples);
} else {
    $ppvcf1 = preprocessVCF($cfg, $infile1, $threads, $tmpFolder, \%ExcludeSamples);
}
if(defined($pInfile2)) {
    my $sample = (split /_/, (split /\./, (split /\//, $pInfile2)[-1], -1)[0], -1)[0];
    my $resultfile = "$tmpFolder/${sample}_preprocessed.vcf";
    $ppvcf2 = mergeGT($pInfile2, $resultfile, \%ExcludeSamples);
} else {
    $ppvcf2 = preprocessVCF($cfg, $infile2, $threads, $tmpFolder, \%ExcludeSamples);
}

# Extract SV Insertions
my $InsSV1 = extractSVIns($ppvcf1, $tmpFolder);
my $NovelSampleFile = getSVsampleInfo($InsSV1, $tmpFolder);
my $InsSV2 = extractSVIns($ppvcf2, $tmpFolder);

# Comparison of the VCFs using Truvari bench command
my $truvariResFolder = "$tmpFolder/Truvari";
my $status = system("$truvari bench -r 1000 -C 1000 -O 0.8 -p 0.8 -P 0.0 -s 50 -S 15 --sizemax 100000 -b $InsSV2 -c $InsSV1 -o $truvariResFolder");
if($status==0) {
    print "VCF comparison using Truvari completed \n";
    extractNovelInsertions($truvariResFolder, $tmpFolder, $resultfolder, $NovelSampleFile, $karyo, $Rscript, $dpi);
} else {
    print "Error in VCF comparison using truvari\n";
    exit;
}

# Removing temporary folders
cleanUp($tmpFolder);
exit;

#Subroutines
sub preprocessVCF
{
	my $infile = shift;
	my $resfolder = shift;
	my $thread = shift;
	my $sample = (split'\.',(split'\/',$infile)[-1],-1)[0];
	my $resultfile1 = "$resfolder\/$sample\_single_allelic.vcf";
	my $resultfile2 = "$resfolder\/$sample\_single_allelic_decomposed.vcf";
	my $resultfile3 = "$resfolder\/$sample\_single_allelic_decomposed_GT_Merged.vcf";

	print "Started pre-processing of vcf file $infile\n\n";

	# Covert to single allelic #
	my $status = 0;
	$status = system("bcftools norm --threads $thread -m- $infile -o $resultfile1");

	if($status == 0)
	{
		print "Successfully converted to singlealleic vcf file\n";
		$status = system("rtg vcfdecompose --break-indels --break-mnps -Z -i $resultfile1 -o $resultfile2");
		
		if($status == 0)
		{
			print "Successfully decomposed the vcf file\n";
			
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
			my $chr = (split'\.',$data[0],-1)[-1];
			
			my $startPos = $data[1];
			my $ref = $data[3];
			my $alt = $data[4];

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

	print "\nPreprocessing of input vcf file $infile completed successfully !!\n";
	print '######################################'."\n\n";
	return $resultfile3;
}

sub extractSVIns
{
	my $infile = shift;
	my $resfolder = shift;
	my $sample = (split'\_',(split'\.',(split'\/',$infile)[-1],-1)[0],-1)[0];
	my $resultfile1 = "$resfolder\/$sample\_SV_INS.vcf";
	my $resultfile2 = "$resfolder\/$sample\_SV_INS.vcf.gz";
	open OUT,">$resultfile1" or die "Can't open $resultfile1 for writing\n";
	open IN,"<$infile" or die "Can't open $infile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		if(/^\s*\#+.*/)
		{
			s/CHM13v2\.//gi;
			print OUT "$_\n";
			next;
		}
		my @data = split"\t",$_,-1;
		my $ref = $data[3];
		my $alt = $data[4];
		my $refLen = length($ref);
		my $altLen = length($alt);
		my $indelLength = abs($altLen-$refLen);
		if( ($indelLength>=50) && ($altLen>$refLen))
		{
			print OUT "$_\n";
		}
	}
	close IN;
	close OUT;
	print "Extracted SV Insertions from vcf file of sample $sample\n\n";
	my $status = system("bgzip -c $resultfile1 > $resultfile2");
	if($status==0)
	{
		print "Bgzip successful : $resultfile1\n";
		$status = system("tabix -p vcf $resultfile2");
		if($status==0)
		{
			print "Index generated for the file $resultfile2\n";
			return $resultfile2;
		}
		else
		{
			print "Error in indexing of the file : $resultfile2\n";
			exit;
		}
	}
	else
	{
		print "Error in bgzip of file : $resultfile1\n";
		exit;
	}
}

sub extractNovelInsertions
{
	my $infolder = shift;
	my $infile = "$infolder\/fp.vcf.gz";
	my $resultfile = "Novel_Insertions.fa";

	open OUT,">$resultfile" or die "Can't open $resultfile for writing\n";
	my $novelSeqNo = 0;
	if($infile=~/\.gz$/)
	{
		open(IN,"gunzip -c $infile |") or die "can't open $infile for reading";
	}
	else
	{
		open IN,"<$infile" or die "Can't open $infile for reading\n";
	}
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		next if(/^\s*\#+.*/);
		my @data = split"\t",$_,-1;
		my $chr = $data[0];
		my $seq = $data[4];
		$novelSeqNo++;
		my $seqHeader = "\>$novelSeqNo\_$chr\_$data[1]";
		print OUT "$seqHeader\n$seq\n";
	}
	close IN;
	close OUT;
	print "Total number of sequences before clustering : $novelSeqNo\n";
	
	&clusterNovelIns($resultfile,$infolder);
}

sub clusterNovelIns
{
	my $infile = shift;
	my $tmpFolder = shift;
	my $resultfile = 'NovelSequences.fa';
	my $resultfile2 = 'CDHIT_ClusterInfo.txt';
	my $uniqInsLen = 0;
	if(-f $resultfile)
	{
		unlink($resultfile);
	}
	if(-f $resultfile2)
	{
		unlink($resultfile2);
	}

	print "Clustering novel seq using CD-HIT\n";
	my %info = ();
	open IN,"<$infile" or die "Can't open $infile for reading\n";
	while(my $header = <IN>)
	{
		my $seq = <IN>;
		chomp($header,$seq);
		my @data = split'\_',$header,-1;
		my $group = join'_',@data[1..2];
		push(@{$info{$group}},"$header\n$seq");
	}
	close IN;

	open OUT,">$resultfile" or die "Can't open $resultfile for writing\n";
	open OUT2,">$resultfile2" or die "Can't open $resultfile2 for writing\n";
	foreach my $eachGroups(sort keys %info)
	{
		my $SeqNo = scalar(@{$info{$eachGroups}});
		if($SeqNo>1)
		{
			my $fname = $eachGroups;
			$fname=~tr/>//d;
			my $tmpfile = "$tmpFolder\/$fname\.fa";
			my $tmpresfile = "$tmpFolder\/$fname\_clustered\.fa";
			my $clusterfile = $tmpresfile.'.clstr';
			open TMP,">$tmpfile" or die "Can't open $tmpfile for writing\n";
			foreach my $eachSeq(@{$info{$eachGroups}})
			{
				print TMP "$eachSeq\n";
			}
			close TMP;
			my $retval = system("cd-hit-est \-i $tmpfile \-o $tmpresfile \-M 0 \-d 50 \-T 0 \-c 0.9");
			if($retval==0)
			{
				unlink($tmpfile);
				open IN1,"<$tmpresfile" or die "Can't open $tmpresfile for reading\n";
				while(<IN1>)
				{
					chomp;
					next if(/^\s*$/);
					if(/^\s*\>+.*/)
					{
						print OUT "$_\n";
					}
					else
					{
						my $l = length($_);
						$uniqInsLen+=$l;
						print OUT "$_\n";
					}
				}
				close IN1;
				unlink($tmpresfile);

				open IN2,"<$clusterfile" or die "Can't open $clusterfile for reading\n";
				while(<IN2>)
				{
					chomp;
					next if(/^\s*$/);
					print OUT2 "$_\n";
				}
				close IN2;
				unlink($clusterfile);

				print "Group $eachGroups : CD HIT run successfully !\n";
			}
			else
			{
				print "CD HIT run ERROR for the group <$eachGroups>\n";
			}
		}
		else
		{
			foreach my $eachSeq(@{$info{$eachGroups}})
			{
				foreach my $eachinfo(split'\n',$eachSeq,-1)
				{

					if($eachinfo=~/^\s*\>+.*/)
					{
						print OUT "$eachinfo\n";
					}
					else
					{
						my $l = length($eachinfo);
						$uniqInsLen+=$l;
						print OUT "$eachinfo\n";
					}
				}
			}
		}
	}
	close OUT;
	close OUT2;
	
	unlink $infile;
	unlink $resultfile2;

	print "############### Program Completed Successfully  ! ! ! #################\n\n";
	print "\n\nUnique Insertion Length After Clustering : $uniqInsLen bp\n";
}

sub deleteTmpFolder
{
	my $infolder = shift;
	foreach my $eachfiles(glob("$infolder\/*"))
	{
		unlink $eachfiles;
	}
	rmdir $infolder;
}
