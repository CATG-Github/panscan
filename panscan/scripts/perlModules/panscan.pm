package perlModules::panscan;

use strict;
use warnings;
use Cwd;
use Exporter 'import';
use Parallel::ForkManager;

our @EXPORT_OK = qw(preprocessVCF mergeGT cleanUp print_novel_seq_help print_uniq_var_help generateInput compareSNPDB compareINDELDB generateInputSV compareSVINS parseSVINS compareSVDEL parseSVDEL extractSVIns extractNovelInsertions clusterNovelIns getSVsampleInfo generateIdeogram runTruvari chrSplitVcf mergeTruvariVcfs);


sub cleanUp
{
        my $infolder = shift;
	foreach my $eachfiles(glob("$infolder\/*"))
	{
		if(-d $eachfiles)
		{
			&cleanUp($eachfiles);
		}
		else
		{
			unlink($eachfiles);
		}
	}
        rmdir $infolder;
}

sub print_novel_seq_help
{
print <<EOF;
        Usage: perl findNovelSeq.pl [options]

        Options:
        --i FILE       Specify the input pangenome vcf file. 
        --r FILE       Specify the reference pangenome vcf file to compare.
        --pInp FILE    Specify the pre-processed input pangenome vcf file. 
        --pRef FILE    Specify the pre-processed reference pangenome vcf file to compare.
        --exclude LIST Specify the sample/s to be excluded (',' separated)
        --t NUMBER     Specify the number of threads required.
			default : 1
        --dt NUMBER     Specify the number of threads for RTG decompose.
			default : 1 (Please note that decompose requires large memory and if
                                     increase the threads, please make sure that the system is 
				     having enough memory )
	--dpi          Specify the dpi of the Ideogram image.
			default : 600
	--genome       Specify reference genome ( HG38 or CHM13 )
			default : HG38
        --help         Display this help message.

EOF
}


sub print_uniq_var_help
{
print <<EOF;
        Usage: perl findUniqVariants.pl [options]

        Options:
        --i FILE       Specify the input vcf files.
	--pbsv         If input is pbsv vcf file for SV comparison.
		        default : off
        --t TEXT       Specify the variant type to compare (SNP/INDEL/SV).
        --db FILE      Specify the databases ( ',' sepearted for multiple databases ).
                       Available databases for
                                        SNP   : dbSNP,gnomAD,1000Genomes,GME
                                        InDel : gnomAD,1000Genomes,GME
                                        SV insertion    : DGV
                                        SV deletion     : DGV,1000Genomes
				default : ALL
        --overlap NUMBER     Specify the percentage of overlap for the SV comparison.
			default : 80
        --help         Display this help message.

EOF
}

sub preprocessVCF
{
	# Specify the path to the config file
	my $cfg = shift;
	my $infile = shift;
	my $threads = shift;
	my $tmp_folder = shift;
	my $dthreads = shift;
	my $refExcludeSamples = shift;

	my $sample = (split'\.',(split'\/',$infile)[-1],-1)[0];
	my $resultfile1 = "$tmp_folder\/$sample\_single_allelic.vcf";
	my $resultfile2 = "$tmp_folder\/$sample\_single_allelic_decomposed.vcf";
	my $resultfile3 = "$tmp_folder\/$sample\_preprocessed.vcf";

	# Retrieve the paths of the tools from the config file
	my $bcftools = $cfg->{tools}->{bcftools};
	my $rtg = $cfg->{tools}->{rtg};

	# Covert to single allelic #
	my $status = 0;
	$status = system("$bcftools norm --threads $threads -m- $infile -o $resultfile1");

	if($status == 0)
	{
		print "Successfully converted to singlealleic vcf file\n";

		&splitVcf($resultfile1,$tmp_folder,$dthreads); # Split input file w.r.t threads #
		&decomposeVcf($tmp_folder,$dthreads,$resultfile2,$rtg);
		my $retval = system("$bcftools sort \-T $tmp_folder $resultfile2 \-o $resultfile3");
		if($retval != 0 )
		{
			print "Error in sorting the devomposed vcf file <$resultfile2>\n";
			exit;
		}
	}
	else
	{
		print "Error in converting to single allelic\n";
		exit;
	}
	
	return $resultfile3;
}

sub mergeGT
{
	my $infile = shift;
	my $resultfile = shift;
	my $refExcludeSamples = shift;

	my @header = ();
	my @data = ();
	my @sampleindex = ();
	my %GenoTypeInfo = ();
	my $total = 0;
	my $merged = 0;
	my %refIndex = ();
	my %Chrs = map { $_ => $_ } (map { "chr$_" } (1..22, 'X', 'Y', 'M'));

	my %SkipSample = ();
	open OUT,">$resultfile" or die "Can't open $resultfile for writing\n";
	open IN,"<$infile" or die "Can't open $infile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);

		if(/^\s*\#CHR.*/) # Extracting the samples from headers #
		{
			@header = split'\t',$_,-1;
			for(my $i=0;$i<=$#header;$i++)
			{
				if(defined($$refExcludeSamples{$header[$i]}))
				{
					$SkipSample{$i} = 0;
				}
			}
		}

		if(/^\s*\#+.*/) # Skipping VCF headers #
		{
			print OUT "$_\n";
			next;
		}
		$total++;
		@data = split"\t",$_,-1;
		if(scalar(keys %SkipSample)>0)
		{
			my $skipGT = 0;
			my $sampleGT = 0;
			for(my $i=9;$i<=$#data;$i++)
			{
				foreach my $eachGT(split'\D',$data[$i],-1)
				{
					$eachGT = 0 if($eachGT=~/^\s*$/);
					if(defined($SkipSample{$i}))
					{
						$skipGT+=$eachGT;
					}
					else
					{
						$sampleGT+=$eachGT;
					}
				}
			}
			if( ($sampleGT==0) && ($skipGT>0))
			{
				next;
			}
		}
		my $chr = $data[0];
		unless(defined($Chrs{$chr}))
		{
			print "Warning : skipping <$chr>. Please make sure the chromosome names in the vcf files are in the format : chrNo (chr1,chr2,chrX etc.)\n";
			next;
		}

		my $startPos = $data[1];
		my $ref = $data[3];
		my $alt = $data[4];
		my $l1 = length($ref);
		my $l2 = length($alt);
		if( ($l1>$l2) || (abs($l2-$l1)<50) )
		{
			next;
		}

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
					if(defined($refIndex{$i}))
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
			if($data[$i]=~/\|+/g)
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
			else
			{
				$refIndex{$i} = 0;
				if($data[$i] eq '.')
				{
					$data[$i] = 0;
				}
				$GenoTypeInfo{"$chr\t$startPos\t$ref\t$alt"}{$header[$i]}+=$data[$i];
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
			if(defined($refIndex{$i}))
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
	return $resultfile
}

sub splitVcf
{
	my $infile = shift;
	my $tmpfolder = shift;
	my $threads = shift;

	my $line_count = `wc -l < $infile`;
	my $header_count = `grep -c "^#" < $infile`;
	chomp($line_count,$header_count);
	my $TotalRec = ($line_count-$header_count);
	my $RecPerfile = ($TotalRec/$threads);
	my $finalFileRec = 0;
	my @RecPerFile = ();
	if($RecPerfile=~/\./g)
	{
		$RecPerfile = int($RecPerfile);
		$finalFileRec = $TotalRec-($RecPerfile*($threads-1));
		@RecPerFile = (undef,($RecPerfile)x($threads-1),$finalFileRec);
	}
	else
	{
		@RecPerFile = (undef,($RecPerfile)x$threads);
	}

	my @fileHandlers = ();
	# Opening files #
	for(my $i=1;$i<=$threads;$i++)
	{
		my $tmpfile = "$tmpfolder\/in-vcf\_$i\.vcf";
		open my $OUT,">$tmpfile" or die "Can't open $tmpfile for writing\n";
		$fileHandlers[$i] = $OUT;
	}

	my $header = '';
	my $lineNo = 0;
	my $fileNo = 1;
	open IN,"<$infile" or die "Can't open $infile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		if(/^\s*\#+.*/)
		{
			if($header eq '')
			{
				$header	= $_;
			}
			else
			{
				$header.="\n$_";
			}
		}
		else
		{
			my $FH = $fileHandlers[$fileNo];
			my $TotRecToCheck = $RecPerFile[$fileNo];
			if($lineNo==0)
			{
				print $FH "$header\n";
			}

			my @data = split"\t",$_,-1;
			if(length($data[4])>=50)
			{
				print $FH "$_\n";
			}
			$lineNo++;

			if($lineNo==$RecPerFile[$fileNo])
			{
				close $FH;
				$fileNo++;
				$lineNo = 0;
			}
		}
	}
	close IN;
}

sub decomposeVcf
{
	my $inFolder = shift;
	my $threads = shift;
	my $resultfile = shift;
	my $rtg = shift;

	my $parellelFork = Parallel::ForkManager->new($threads);

	my $retval = 0;
	for(my $i=1;$i<=$threads;$i++)
	{
		$parellelFork->start and next;

		my $inVcf = "$inFolder\/in-vcf\_$i\.vcf";
		my $resVcf = "$inFolder\/decomposed-vcf\_$i\.vcf";

		if(-e $inVcf)
		{
			my $status = system("$rtg vcfdecompose --break-indels --break-mnps -Z -i $inVcf -o $resVcf");
			if($status!=0)
			{
				$retval = 1;
				print "Error in decomposition of the file <$inVcf>\n";
				exit;
			}
			else
			{
				unlink $inVcf;
			}
		}
		else
		{
			$retval = 1;
			print "Error ! the file <$inVcf> not found for decomposition.\n";
			exit;
		}
		$parellelFork->finish;
	}
	$parellelFork->wait_all_children;
	if($retval==0)
	{
		print "Decomposition of all files done !!\n";
		&mergeVcf($inFolder,$threads,$resultfile,'decomposed');
	}
	return $retval;
}

sub mergeVcf
{
	my $inFolder = shift;
	my $thread = shift;
	my $resultfile = shift;
	my $tag = shift;

	open OUT,">$resultfile" or die "Can't open $resultfile for writing\n";
	foreach my $i(1..$thread)
	{
		my $infile = "$inFolder\/$tag\-vcf\_$i\.vcf";
		open IN,"<$infile" or die "Can't open $infile for reading\n";
		while(<IN>)
		{
			chomp;
			next if(/^\s*$/);
			if(/^\s*#+.*/)
			{
				if($i==1)
				{
					print OUT "$_\n";
				}
				next;
			}
			my @data = split"\t",$_,-1;
			if(length($data[4])>=50)
			{
				print OUT "$_\n";
			}
		}
		close IN;
		unlink $infile;
	}
	close OUT;
}

#################### Uniq_variants SUBROUTINES ###################
sub generateInput
{
	# Accepting arguments #
	my $infile = shift;
	my $resfolder = shift;
	my $db = shift;
	my $rFolder = shift;

	# Resultfiles #
	my $resultfile = "$resfolder\/$db\.txt";
	my $resultfile1 = "$resfolder\/SV_INS.txt";
	my $resultfile2 = "$resfolder\/SV_DEL.txt";
	my $resultfile3 = "$rFolder\/SV_Stat.txt";


	# Opening the vcf file for reading #
	if($infile=~/\.gz$/)
	{
		open(IN,"gunzip -c $infile |") or die "can't open $infile for reading";
	}
	else
	{
		open IN,"<$infile" or die "Can't open $infile for reading\n";
	}
	
	open OUT,">$resultfile" or die "Can't open $resultfile for writing\n";
	open OUT1,">$resultfile1" or die "Can't open $resultfile1 for writing\n";
	open OUT2,">$resultfile2" or die "Can't open $resultfile2 for writing\n";
	open OUT3,">$resultfile3" or die "Can't open $resultfile3 for writing\n";
	my $totSVs = 0;
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/); # Skipping blank lines present, if any #
		next if(/^\s*\#+.*/); # Skipping headers #

		my @data = split"\t",$_,-1; # split & extract data #
		my $chr = $data[0];
		my $pos = $data[1];
		my $ref = $data[3];
		my $l1 = length($ref);
		my $alt = $data[4];
		
		foreach my $eachAlt(split'\,',$alt,-1) # Multi allelic variants #
		{
			my $l2 = length($eachAlt);
			my $diff = $l2-$l1;

			if( ($l1==1) && ($l2==$l1) ) # SNP #
			{
				if($db eq 'SNP')
				{
					print OUT "$chr\t$pos\t$ref\t$eachAlt\n";
				}
			}
			elsif(abs($diff)<50) # INDEL #
			{
				if($db eq 'INDEL')
				{
					print OUT "$chr\t$pos\t$ref\t$eachAlt\n";
				}
			}
			else # SV #
			{
				if($db eq 'SV')
				{
					my $chr = $data[0];
					$chr=~s/chr//gi;
					if($chr eq 'X')
					{
						$chr = 23;
					}
					if($chr eq 'Y')
					{
						$chr = 24;
					}
					if($diff>0)
					{
						$totSVs++;
						my $end = $pos+$diff;
						print OUT1 "$chr\t$pos\t$end\n";
					}
					else
					{
						$totSVs++;
						my $end = $pos+abs($diff);
						print OUT2 "$chr\t$pos\t$end\n";
					}
				}
			}
		}
	}
	close IN; # Closing the filehandler #
	close OUT;
	close OUT1;
	close OUT2;
	close OUT3;

	if($db ne 'SV')
	{
		unlink $resultfile1;
		unlink $resultfile2;
		unlink $resultfile3;
		return $resultfile;
	}
	else
	{
		unlink $resultfile3;
		return($totSVs,$resultfile1,$resultfile2,$resultfile3);
	}
}

sub compareSNPDB
{
	my $refInfile = shift;
	my $db = shift;
	my $sample = shift;
	my $dbName = shift;
	my $resfolder = shift;
	my $dbNo = shift;
	my $totalDB = shift;

	my $resultfile1 = '';
	if($dbNo==$totalDB)
	{
		$resultfile1 = "$resfolder\/After-comparison-with\-$dbName\_Final-Novel_SNPs.txt";
	}
	else
	{
		$resultfile1 = "$resfolder\/After-comparison-with\-$dbName\_Unique_SNPs.txt";
	}
	my $resultfile2 = "$resfolder\/$sample\_SNP\_Stat.txt";

	my %DB_Variants = ();
	my $total = 0;
	open IN1,"<$$refInfile" or die "Can't open $$refInfile for reading\n";
	while(<IN1>)
	{
		chomp;
		next if(/^\s*$/);
		$DB_Variants{$_} = 0;
		$total++;
	}
	close IN1;

	open IN2,"<$db" or die "Can't open $db for reading\n";
	while(<IN2>)
	{
		chomp;
		next if(/^\s*$/);
		if(defined($DB_Variants{$_}))
		{
			$DB_Variants{$_} = 1;
		}
	}
	close IN2;

	my $common = 0;
	open OUT1,">>$resultfile1" or die "Can't open $resultfile1 for writing\n";
	while(my($k,$v) = each %DB_Variants)
	{
		if($v==0)
		{
			print OUT1 "$k\n";
		}
		else
		{
			$common++;
		}
	}
	close OUT1;

	my $unique = $total-$common;
	open OUT2,">>$resultfile2" or die "Can't open $resultfile2 for writing\n";
	if($dbNo==1)
	{
		print OUT2 "Database\tTotal input SNPs\tCommon SNPs\tUnique SNPs\n";
	}
	print OUT2 "$dbName\t$total\t$common\t$unique\n";
	if($dbNo==$totalDB)
	{
		print OUT2 "\n\nTotal Novel SNPs : $unique\n";
	}
	close OUT2;
	
	$$refInfile = $resultfile1;
}

sub compareINDELDB
{
	my $refInfile = shift;
	my $db = shift;
	my $sample = shift;
	my $dbName = shift;
	my $resfolder = shift;
	my $dbNo = shift;
	my $totalDB = shift;

        my $resultfile1 = '';
        if($dbNo==$totalDB)
        {
                $resultfile1 = "$resfolder\/After-comparison-with\-$dbName\_Final-Novel_INDELs.txt";
        }
        else
        {
		$resultfile1 = "$resfolder\/After-comparison-with\-$dbName\_Unique_INDELs.txt";
        }

	my $resultfile2 = "$resfolder\/$sample\_INDEL\_Stat.txt";

	my %DB_Variants = ();
	my $total = 0;
	open IN1,"<$$refInfile" or die "Can't open $$refInfile for reading\n";
	while(<IN1>)
	{
		chomp;
		next if(/^\s*$/);
		$DB_Variants{$_} = 0;
		$total++;
	}
	close IN1;

	open IN2,"<$db" or die "Can't open $db for reading\n";
	while(<IN2>)
	{
		chomp;
		next if(/^\s*$/);
		if(defined($DB_Variants{$_}))
		{
			$DB_Variants{$_} = 1;
		}
	}
	close IN2;

	my $common = 0;
	open OUT1,">>$resultfile1" or die "Can't open $resultfile1 for writing\n";
	while(my($k,$v) = each %DB_Variants)
	{
		if($v==0)
		{
			print OUT1 "$k\n";
		}
		else
		{
			$common++;
		}
	}
	close OUT1;

	my $unique = $total-$common;
	open OUT2,">>$resultfile2" or die "Can't open $resultfile2 for writing\n";
	if($dbNo==1)
	{
		print OUT2 "Database\tTotal input INDELs\tCommon INDELs\tUnique INDELs\n";
	}
	print OUT2 "$dbName\t$total\t$common\t$unique\n";
        if($dbNo==$totalDB)
        {
                print OUT2 "\n\nTotal Novel INDELs : $unique\n";
        }
	close OUT2;
	
	$$refInfile = $resultfile1;
}

sub generateInputSV
{
	my $infile = shift;
	my $resFolder = shift;
	my $resFolder2 = shift;
	my $sample = 'Input-cohort';

	my $resultfile1 = "$resFolder\/$sample\_SV_INS.txt";
	my $resultfile2 = "$resFolder\/$sample\_SV_DEL.txt";
	my $resultfile3 = "$resFolder\/$sample\_SV_DUP.txt";
	my $resultfile4 = "$resFolder\/$sample\_SV_INV.txt";
	my $resultfile5 = "$resFolder2\/$sample\_SV_Stat.txt";

	my %VarInfo = ();
	my @header = ();
	open OUT1,">$resultfile1" or die "Can't open $resultfile1 for writing\n";
	open OUT2,">$resultfile2" or die "Can't open $resultfile2 for writing\n";
	open OUT3,">$resultfile3" or die "Can't open $resultfile3 for writing\n";
	open OUT4,">$resultfile4" or die "Can't open $resultfile4 for writing\n";
	open IN,"<$infile" or die "Can't open $infile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		if(/^\s*\#CHR.*/)
		{
			@header = split'\t',$_,-1;
		}
		next if(/^\s*\#+.*/);
		my @data = split"\t",$_,-1;
		if($data[6] eq 'PASS')
		{
			my ($SVType,$SVLen,$end) = ('NA','NA','NA');
			foreach my $eachInfo(split'\;',$data[7],-1)
			{
				my($k,$v) = split'\=',$eachInfo,-1;
				if($k eq 'SVTYPE')
				{
					$SVType = $v;
				}
				elsif($k eq 'END')
				{
					$end = $v;
				}
				elsif($k eq 'SVLEN')
				{
					$SVLen = $v;
				}
			}
			my $chr = $data[0];
			$chr=~s/chr//gi;
			if($chr eq 'X')
			{
				$chr = 23;
			}
			if($chr eq 'Y')
			{
				$chr = 24;
			}
			
			if($SVType eq 'INS')
			{
				$end = $data[1]+$SVLen;
				if(length($chr)==1)
				{
					print OUT1 "$chr\t$data[1]\t$end\n";
				}
			}
			elsif($SVType eq 'DEL')
			{
				if(length($chr)==1)
				{
					print OUT2 "$chr\t$data[1]\t$end\n";
				}
			}
			elsif($SVType eq 'DUP')
			{
				print OUT3 "$chr\t$data[1]\t$end\n";
			}
			elsif($SVType eq 'INV')
			{
				print OUT4 "$chr\t$data[1]\t$end\n";
			}
			
			foreach my $i(9..$#header)
			{
				my $gt = (split'\:',$data[$i],-1)[0];
				if( ($gt eq '1/0') ||($gt eq '1/1') ||($gt eq '0/1') ||($gt eq './1') || ($gt eq '1/.'))
				{
					$VarInfo{$header[$i]}{'TYPE'}{$SVType}++;
					$VarInfo{$header[$i]}{'TYPE'}{'TOTAL'}++;
				}
			}
		}
	}
	close IN;
	close OUT1;
	close OUT2;
	close OUT3;
	close OUT4;

	open OUT5,">$resultfile5" or die "Can't open $resultfile5 for writing\n";
	my @SVTypes = ('INS','DEL','DUP','INV','BND');
	print OUT5 "Sample\tTotal Variants\t",(join"\t",@SVTypes)."\n";
	my $totSVs = 0;
	foreach my $i(9..$#header)
	{
		my $sample = $header[$i];
		my @sample_result = ($sample);
		foreach my $eachType('TOTAL',@SVTypes)
		{
			if(defined($VarInfo{$header[$i]}{'TYPE'}{$eachType}))
			{
				$totSVs++;
				push(@sample_result,$VarInfo{$header[$i]}{'TYPE'}{$eachType});
			}
			else
			{
				push(@sample_result,0);
			}
		}
		my $eachResult = join"\t",@sample_result;
		print OUT5 "$eachResult\n";
	}
	close OUT5;
	if(-z $resultfile1)
	{
		unlink $resultfile1;
	}
	if(-z $resultfile2)
	{
		unlink $resultfile2;
	}
	if(-z $resultfile3)
	{
		unlink $resultfile3;
	}
	if(-z $resultfile4)
	{
		unlink $resultfile4;
	}
	return ($totSVs,$resultfile1,$resultfile2,$resultfile3,$resultfile4);
}

sub compareSVINS
{
	my $infile = shift;
	my $tmpfolder = shift;
	my $resultfolder = shift;
	my $svOverlap = shift;
	my $cfg = shift;

	my $DGV_SV_INS = "$cfg->{databases}->{dbpath}/DGV_SV_INS.txt";
	my $resultfile = "$tmpfolder\/After-comparison-with-DGV_Final-Unique-SVs.txt";
	my $retval = system("java $cfg->{classes}->{reciprocal_overlap} $infile $DGV_SV_INS $svOverlap > $resultfile");
	if($retval==0)
	{
		&parseSVINS($resultfile,$resultfolder);
		print "SV insertion comparison completed !!\n";
	}
	else
	{
		print "Problem in the reciprocal overlap\n";
		exit;
	}
}

sub parseSVINS
{
	my $infile = shift;
	my $resultfolder = shift;
	my $percentage = shift;

	my $resultfile1 = "$resultfolder\/After-comparison-with-DGV_Final-Unique-SV-Insertions.txt";
	my $resultfile2 = "$resultfolder\/SV-Insetion_Comparison_Stat.txt";

	my $total = 0;
	my $common = 0;
	my $unique = 0;

	open OUT1,">$resultfile1" or die "Can't open $resultfile1 for writng\n";
	open IN,"<$infile" or die "Can't open $infile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		$total++;
		my @data = split"\t",$_,-1;

		if($data[1]>0)
		{
			$common++;
		}
		elsif($data[1]==0)
		{
			my($chr,$cord) = split'\:',$data[0],-1;
			my($cord1,$cord2) = split'\-',$cord,-1;
			$chr=~tr/Chr//d;
			print OUT1 "$chr\t$cord1\t$cord2\n";
			$unique++;
		}
	}
	close IN;
	close OUT1;

	open OUT2,">$resultfile2" or die "Can't open $resultfile2 for writng\n";
	print OUT2 "Total insertions\tCommon insertions\tUnique insertions\n";
	print OUT2 "$total\t$common\t$unique\n";
	print OUT2 "\n\nTotal Novel SV Insertions : $unique\n";
	close OUT2;
}

sub compareSVDEL
{
	my $infile = shift;
	my $tmpfolder = shift;
	my $resultfolder = shift;
	my $svOverlap = shift;
	my $refDBorder = shift;
	my $refDBs = shift;
	my $cfg = shift;

	my $totalDBs = scalar(@$refDBorder);
	my $no = 1;
	my $retFile = '';
	foreach my $eachDb(@$refDBorder)
	{
		my $dbfile = $$refDBs{$eachDb};
		my $resultfile = '';
		$resultfile = "$tmpfolder\/After-comparison-with\-$eachDb\-SV-Deletions.txt";
		if($no==$totalDBs)
		{
			$resultfile = "$tmpfolder\/After-comparison-with\-$eachDb\_Final-Unique-SV-Deletions.txt";
		}
		
		my $retval = 'NA';

		if($no==1)
		{
			$retval = system("java $cfg->{classes}->{reciprocal_overlap} $infile $dbfile $svOverlap > $resultfile");
		}
		else
		{
			$retval = system("java $cfg->{classes}->{reciprocal_overlap} $retFile $dbfile $svOverlap > $resultfile");
		}
		if($retval==0)
		{
			$retFile = &parseSVDEL($resultfile,$resultfolder,$eachDb,$totalDBs,$no);
			if($no==$totalDBs)
			{
				print "SV Deletion comparison completed !!\n";
			}
		}
		else
		{
			print "Problem in the reciprocal overlap\n";
			exit;
		}
		$no++;
	}
}

sub parseSVDEL
{
	my $infile = shift;
	my $resultfolder = shift;
	my $db = shift;
	my $totalDB = shift;
	my $N = shift;

	my $resultfile1 = "$resultfolder\/After-comparison-with\-$db\_SV-Deletions.txt";
	if($N==$totalDB)
	{
		$resultfile1 = "$resultfolder\/After-comparison-with\-$db\_Final-Unique-SV-Deletions.txt";
	}
	my $resultfile2 = "$resultfolder\/SV-Deletion_Comparison_Stat.txt";

	my $total = 0;
	my $common = 0;
	my $unique = 0;

	open OUT1,">$resultfile1" or die "Can't open $resultfile1 for writing\n";
	open IN,"<$infile" or die "Can't open $infile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		$total++;
		my @data = split"\t",$_,-1;

		if($data[1]>0)
		{
			$common++;
		}
		elsif($data[1]==0)
		{
			my($chr,$cord) = split'\:',$data[0],-1;
			my($cord1,$cord2) = split'\-',$cord,-1;
			$chr=~tr/Chr//d;
			print OUT1 "$chr\t$cord1\t$cord2\n";
			$unique++;
		}
	}
	close IN;
	close OUT1;

	if($N==1)
	{
		open OUT2,">$resultfile2" or die "Can't open $resultfile2 for writng\n";
		print OUT2 "Database\tTotal deletions\tCommon deletions\tUnique deletions\n";
	}
	else
	{
		open OUT2,">>$resultfile2" or die "Can't open $resultfile2 for writng\n";
	}

	print OUT2 "$db\t$total\t$common\t$unique\n";
	if($N==$totalDB)
        {
                print OUT2 "\n\nTotal Novel SV Deletions : $unique\n";
        }
	close OUT2;

	return $resultfile1;
}

#################### novel_seq SUBROUTINES ###############################################

sub extractSVIns
{
	my $infile = shift;
	my $resfolder = shift;
	my $bgzip = shift;
	my $tabix = shift;

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
	my $status = system("$bgzip \-c $resultfile1 > $resultfile2");
	if($status==0)
	{
		print "Bgzip successful : $resultfile1\n";
		$status = system("$tabix \-p vcf $resultfile2");
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
	my $tmpfolder = shift;
	my $resfolder = shift;
	my $sampleInfofile = shift;
	my $karyo = shift;
	my $Rscript = shift;
	my $dpi = shift;
	my $threads = shift;
	my $cdhit_est = shift;
	
	my $infile = "$tmpfolder\/fp.vcf";
	my $resultfile = "$tmpfolder\/Novel_Insertions.fa";

	open OUT,">$resultfile" or die "Can't open $resultfile for writing\n";
	my $novelSeqNo = 0;
	open IN,"<$infile" or die "Can't open $infile for reading\n";
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
	
	&clusterNovelIns($resultfile,$tmpfolder,$resfolder,$sampleInfofile,$karyo,$Rscript,$dpi,$threads,$cdhit_est);
}

sub clusterNovelIns
{
	my $infile = shift;
	my $tmpFolder = shift;
	my $resFolder = shift;
	my $sampleInfoFile = shift;
	my $karyo = shift;
	my $Rscript = shift;
	my $dpi = shift;
	my $threads = shift;
	my $cdhit_est = shift;
	
	my %chrInfo = ("chr1",1,"chr2",2,"chr3",3,"chr4",4,"chr5",5,"chr6",6,"chr7",7,"chr8",8,"chr9",9,"chr10",10,"chr11",11,"chr12",12,"chr13",13,"chr14",14,"chr15",15,"chr16",16,"chr17",17,"chr18",18,"chr19",19,"chr20",20,"chr21",21,"chr22",22,"chrM",'M',"chrX",'X',"chrY",'Y');
	
	my $resultfile = "$resFolder\/NovelSequences.fa";
	my $summaryfile = "$resFolder\/NovelSeqInfo.txt";
	my $summaryfile1 = "$resFolder\/Sample_novel_Seq_count.txt";
	my $summaryfile2 = "$resFolder\/Sample_novel_Seq_summary.txt";
	my $resultfile1 = "$tmpFolder\/NovelSequences.fa";
	my $resultfile2 = "$tmpFolder\/CDHIT_ClusterInfo.txt";
	my $resultfile3 = "$tmpFolder\/Ideogram_input.txt";
	my $uniqInsLen = 0;
	if(-f $resultfile)
	{
		unlink($resultfile);
	}
	if(-f $resultfile2)
	{
		unlink($resultfile2);
	}

	# Read sample info file #
	my %Sinfo = ();
	open IN,"<$sampleInfoFile" or die "Can't open $sampleInfoFile for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		my @data = split"\t",$_,-1;
		$Sinfo{$data[2]} = "$data[0]\t$data[1]\t$data[3]\t$data[4]\t$data[5]";
	}
	close IN;

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

	open OUT1,">$resultfile1" or die "Can't open $resultfile1 for writing\n";
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
			my $retval = system("$cdhit_est \-i $tmpfile \-o $tmpresfile \-M 0 \-d 50 \-T $threads \-c 0.9");
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
						print OUT1 "$_\n";
					}
					else
					{
						my $l = length($_);
						$uniqInsLen+=$l;
						print OUT1 "$_\n";
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
						print OUT1 "$eachinfo\n";
					}
					else
					{
						my $l = length($eachinfo);
						$uniqInsLen+=$l;
						print OUT1 "$eachinfo\n";
					}
				}
			}
		}
	}
	close OUT1;
	close OUT2;

	open OUT3,">$resultfile" or die "Can't open $resultfile for writing\n";
	open OUT4,">$summaryfile" or die "Can't open $summaryfile for writing\n";
	print OUT4 "Novel Seq\tChr\tPos\tSV Insertion Length\tTotal Samples\tSample Names\n";
	open OUT41,">$resultfile3" or die "Can't open $resultfile3 for writing\n";
	print OUT41 "Chr\tStart\tEnd\tValue\n";
	open IN1,"<$resultfile1" or die "Can't open $resultfile1 for reading\n";
	my $sno = 0;
	my %SampleNoveSeqCnt = ();
	my %checkUniq = ();
	while(my $header = <IN1>)
	{
		my $seq = <IN1>;
		$sno++;
		chomp($header,$seq);
		my @h = split'\_',$header,-1;
		unless(defined($checkUniq{$h[1]}{$h[2]}))
		{
			print OUT41 $chrInfo{$h[1]}."\t$h[2]\t$h[2]\t1000\n";
		}
		$checkUniq{$h[1]}{$h[2]} = 0;
		$h[0] = '>'.$sno;
		$header = join'_',@h;

		print OUT3 "$header\n$seq\n";
		if(defined($Sinfo{$seq}))
		{
			print OUT4 "$header\t".$Sinfo{$seq}."\n";
			foreach my $eachSample(split'\,',(split"\t",$Sinfo{$seq},-1)[-1],-1)
			{
				$SampleNoveSeqCnt{$eachSample}++;
			}
		}
		
	}
	close IN1;
	close OUT3;
	close OUT4;
	close OUT41;

	open OUT5,">$summaryfile1" or die "Can't open $summaryfile1 for writing\n";
	print OUT5 "Sample\tTotal Number of Novel Sequences\n";
	foreach my $eachSampes(sort{$SampleNoveSeqCnt{$b}<=>$SampleNoveSeqCnt{$a};} keys %SampleNoveSeqCnt)
	{
		print OUT5 "$eachSampes\t".$SampleNoveSeqCnt{$eachSampes}."\n";
	}
	close OUT5;

	# Generate Ideogram #
	&generateIdeogram($resultfile3,$karyo,$resFolder,$tmpFolder,$Rscript,$dpi);
	
	unlink $infile;
	unlink $resultfile2;
	&cleanUp($tmpFolder);

	print "############### Program Completed Successfully  ! ! ! #################\n\n";

	open OUT6,">$summaryfile2" or die "Can't open $summaryfile2 for writing\n";
	print OUT6 "Total number of Novel Sequences : $sno\n";
	print OUT6 "Total novel sequence length : $uniqInsLen bp\n";
	close OUT6;
}

sub getSVsampleInfo
{
        my $infile = shift;
	my $tmpFolder = shift;
        my @samples = ();
	$infile=~s/\.gz$//;

	my $resultfile = "$tmpFolder\/SV_Ins_Sample_Info.txt";

	open OUT,">$resultfile" or die "Can't open $resultfile for writing\n";
        open IN,"<$infile" or die "Can't open $infile for reading\n";
        while(<IN>)
        {
                chomp;
                next if(/^\s*$/);
                next if(/^\s*\#+\#+.*/);
                my @data = split"\t",$_,-1;
                if(/^\s*\#+CHR.*/)
                {
                        @samples = @data;
                        next;
                }

                my @novelSamples = ();
                for(my $i=9;$i<=$#data;$i++)
                {
                        my $h1 = $data[$i];
                        my $h2 = $data[$i];
                        if($data[$i]=~/\|/g)
                        {
                                ($h1,$h2) = split'\|',$data[$i],-1;
                        }
                        if($data[$i]=~/\//g)
                        {
                                ($h1,$h2) = split'\/',$data[$i],-1;
                        }
                        $h1=0 if($h1 eq '.');
                        $h2=0 if($h2 eq '.');
                        if( ($h1>0) || ($h2>0))
                        {
                                push(@novelSamples,$samples[$i]);
                        }
                }
                my $totSamples = scalar(@novelSamples);
                my $novelSampleList = join',',@novelSamples;
		my $insLen = length($data[4]);
                print OUT "$data[0]\t$data[1]\t$data[4]\t$insLen\t$totSamples\t$novelSampleList\n";
        }
        close IN;
	return $resultfile;
}

sub generateIdeogram
{
        my $novel_seq = shift;
        my $karyotype = shift;
        my $path = shift;
        my $cdir = shift;
        my $Rscript = shift;
        my $dpi = shift;

	my $currdir = getcwd();
        my $resultfile = "$cdir\/Novel_Ideogram.sh";
        $novel_seq = "$currdir\/$novel_seq";
	my $resfile = "NovelSeq_Ideogram.tiff";

	my $range = 1000000;
	my $inpIdeogram = "$currdir\/$cdir\/Ideogram_Input_Shrinked.txt";
	my %Info = ();
	open IN,"<$novel_seq" or die "Can't open $novel_seq for reading\n";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		next if(/^\s*Chr.*/);
		my @data = split"\t",$_,-1;
		my($chr,$start,$end,$n1) = split"\t",$_,-1;
		my $n = int($start/$range);
		$Info{$chr}{$n}+=$n1;
	}
	close IN;

	open OUT1,">$inpIdeogram" or die "Can't open $inpIdeogram for writing\n";
	print OUT1 "Chr\tStart\tEnd\tValue\n";
	my $maxval = 0;
	foreach my $eachChr(sort keys %Info)
	{
		foreach my $cord(sort{$a<=>$b;} keys %{$Info{$eachChr}})
		{
			my $val = $Info{$eachChr}{$cord};
			my $start = ($cord*$range)+1;
			my $end = (($cord+1)*$range);
			my $N= int($range/$val);
			print OUT1 "$eachChr\t$start\t$end\t$val\n";
			if($maxval<$val)
			{
				$maxval = $val;
			}
		}
	}
	close OUT1;



        open OUT,">$resultfile" or die "Can't open $resultfile for writing\n";
print OUT <<"EOF";
setwd("$currdir\/$path")
require(RIdeogram)
human_karyotype <-read.table(\"$karyotype\", sep = \"\\t\", header = T, stringsAsFactors = F)
novel_info <-read.table(\"$inpIdeogram\", sep = \"\\t\", header = T, stringsAsFactors = F)
ideogram(karyotype=human_karyotype,synteny = NULL,label = NULL,label_type = NULL, overlaid = novel_info, colorset1 = c("#FFFFFF", "#D30808", "#D30808"), output = "NovelSeq_Ideogram.svg", width=200,Lx=10000)
svg2tiff("NovelSeq_Ideogram.svg",file = "$resfile", height=5.4,width=7,dpi=$dpi)
file.remove("NovelSeq_Ideogram.svg")
file.remove("Rplots.pdf")

EOF
        close OUT;
        my $retval = system("$Rscript $resultfile");
        if($retval==0)
        {
		unlink $inpIdeogram;
                print "Successfully generated Ideogram <$resfile>\n";
        }
        else
        {
                print "Error in running Rscript command\n";
        }
}

sub runTruvari
{
	my $infile = shift;
	my $reffile = shift;
	my $tmpFolder = shift;
	my $thread = shift;
	my $truvari = shift;
	my $bgzip = shift;
	my $tabix = shift;

	my %Chrs = map { $_ => $_ } (map { "chr$_" } (1..22, 'X', 'Y', 'M'));

	&chrSplitVcf($infile,$tmpFolder,$thread,'input',\%Chrs,$bgzip,$tabix);
	&chrSplitVcf($reffile,$tmpFolder,$thread,'ref',\%Chrs,$bgzip,$tabix);

	my $parellelFork = Parallel::ForkManager->new($thread);
	foreach my $eachChrs(sort{$a cmp $b;} keys %Chrs)
	{
		my $inVcf = "$tmpFolder\/input-vcf\_$eachChrs\.vcf.gz";
		my $refVcf = "$tmpFolder\/ref-vcf\_$eachChrs\.vcf.gz";
		my $resfolder = "$tmpFolder\/Truvari\_$eachChrs";
		$parellelFork->start and next;
		my $status = system("$truvari bench -r 1000 -C 1000 -O 0.8 -p 0.8 -P 0.0 -s 50 -S 15 --sizemax 100000 -b $refVcf -c $inVcf -o $resfolder");
		if($status==0)
		{
			unlink $inVcf;
			unlink "$inVcf\.tbi";
			unlink $refVcf;
			unlink "$refVcf\.tbi";
			print "Truvari ran successfully for $eachChrs\n";
		}
		else
		{
			print "Error in running Truvari on $eachChrs\n";
		}
		$parellelFork->finish;
	}
	$parellelFork->wait_all_children;
	print "Ran Truvari over all files !!\n";
	my $truvariFinalResultFile = "$tmpFolder\/fp.vcf";
	&mergeTruvariVcfs($truvariFinalResultFile,$tmpFolder);
}

sub chrSplitVcf
{
	my $infile = shift;
	my $tmpfolder = shift;
	my $threads = shift;
	my $tag = shift;
	my $refChrs = shift;
	my $bgzip = shift;
	my $tabix = shift;

	my %Chrs = %{$refChrs};

	my %fileHandlers = ();
	foreach my $eachChrs(sort{$a cmp $b;} keys %Chrs)
	{
		my $tmpfile = "$tmpfolder\/$tag\-vcf\_$eachChrs\.vcf";
		open my $OUT,">$tmpfile" or die "Can't open $tmpfile for writing\n";
		$fileHandlers{$eachChrs} = $OUT;
	}

	my $header = '';
	my $lineNo = 0;
	open(IN,"gunzip -c $infile |") or die "can't open $infile for reading";
	while(<IN>)
	{
		chomp;
		next if(/^\s*$/);
		if(/^\s*\#+.*/)
		{
			if($header eq '')
			{
				$header	= $_;
			}
			else
			{
				$header.="\n$_";
			}
			next;
		}
		else
		{
			if($lineNo==0)
			{
				foreach my $eachChrs(keys %Chrs)
				{
					my $FH = $fileHandlers{$eachChrs};
					print $FH "$header\n";
				}
			}

			my @data = split"\t",$_,-1;
			my $FH = $fileHandlers{$data[0]};
			print $FH "$_\n";
			$lineNo++;
		}
	}
	close IN;

	foreach my $eachChrs(keys %Chrs)
	{
		my $FH = $fileHandlers{$eachChrs};
		close $FH;
		my $tmpfile = "$tmpfolder\/$tag\-vcf\_$eachChrs\.vcf";
		my $tmpfile1 = "$tmpfolder\/$tag\-vcf\_$eachChrs\.vcf\.gz";
		my $status = system("$bgzip \-c $tmpfile > $tmpfile1");
		if($status!=0)
		{
			print "Error in bgzip <$tmpfile>\n";
			exit;
		}
		else
		{
			unlink $tmpfile;
			$status = system("$tabix \-p vcf $tmpfile1");
			if($status!=0)
			{
				print "Error in file indexing <$tmpfile1>\n";
				exit;
			}
		}
	}
}

sub mergeTruvariVcfs
{
	my $resultfile = shift;
	my $tmpFolder = shift;

	open OUT,">$resultfile" or die "Can't open $resultfile for writing\n";
	my $no = 0;
	foreach my $eachChrNo(1..22,'X','Y','M')
	{
		my $infile = "$tmpFolder\/Truvari\_chr$eachChrNo\/fp.vcf.gz";
		$no++;
		open(IN,"gunzip -c $infile |") or die "can't open $infile for reading";
		while(<IN>)
		{
			chomp;
			next if(/^\s*$/);
			if(/^\s*#+.*/)
			{
				if($no==1)
				{
					print OUT "$_\n";
				}
				next;
			}
			print OUT "$_\n";
		}
		close IN;
	}
	close OUT;
}

1;
