#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use YAML::XS 'LoadFile';
use perlModules::panscan qw(cleanUp print_help2);
use Cwd qw(getcwd abs_path);

# Use FindBin to determine the directory where this script resides.
my $script_dir = $FindBin::Bin;

# Set config file relative to the script directory.
my $config_file = "$script_dir/config.yaml";

# Check if the config file exists.
if ( ! -e $config_file ) {
    die("Config file '$config_file' not found.\n");
}

# Load the configuration from config.yaml.
my $cfg = LoadFile($config_file);
if ( !$cfg ) {
    die("Failed to read YAML file: $config_file\n");
}

# Variable declarations.
my $inputVcf         = 'NA';
my $varType          = 'NA';
my $db               = 'ALL';
my $overlap          = 80;
my $help             = undef;
my $pbsv             = 0;
my $db_path_override = '';  # New flag to override database path.
my $output           = '';  # New flag for output directory.
my %types            = (
    'SNP'   => 'SNP',
    'INDEL' => 'INDEL',
    'SV'    => 'SV'
);

# Get command-line options, including the new db-path and output flags.
GetOptions(
    'i=s'       => \$inputVcf,
    't=s'       => \$varType,
    'pbsv'      => \$pbsv,
    'db=s'      => \$db,
    'overlap=i' => \$overlap,
    'db-path=s' => \$db_path_override,
    'output=s'  => \$output,     # New output directory flag.
    'help'      => \$help    
) or die "Error in command line arguments\n";

# If help is requested, print help and exit.
if ($help) {
    print_help2();
    exit;
}

$varType = uc($varType);
if ($varType eq 'NA') {
    print "\nERROR: Please provide variant type to compare (SNP/INDEL/SV)!!\n\n";
    print_help2();
    exit;
} elsif ( !defined($types{$varType}) ) {
    print "\nERROR: Unknown variant type <$varType>. Please provide variant type to compare (SNP/INDEL/SV)!!\n\n";
    print_help2();
    exit;
}

if ($inputVcf eq 'NA') {
    print "ERROR: Please provide input VCF file!!\n";
    print_help2();
    exit;
}

# Determine the output directory.
my $resultfolder;
if ($output) {
    # Convert the provided output directory to an absolute path.
    $resultfolder = abs_path($output) or die "Cannot resolve output directory: $output\n";
    mkdir $resultfolder unless (-d $resultfolder);
} else {
    $resultfolder = "${varType}_Comparison_Result";
    mkdir $resultfolder unless (-d $resultfolder);
}

# Create temporary folder.
my $tmpfolder = 'tmp';
unless(-d $tmpfolder) {
    mkdir $tmpfolder;
}

# If result folder exists, clean it and re-create.
if(-d $resultfolder) {
    cleanUp($resultfolder);
    mkdir $resultfolder;
} else {
    mkdir $resultfolder;
}

# Setup database file paths.
# If the --db-path flag is provided, use that directory as the base;
# otherwise, use the relative paths from config.yaml.
my ($dbSNP, $gnomAD, $Genomes_1000, $GME,
    $GNOMAD_INDEL, $Genomes_1000_INDEL, $GME_INDEL,
    $DGV_SV_INS, $Genomes_1000_SV_DEL, $DGV_SV_DEL);

if ($db_path_override) {
    $dbSNP              = "$db_path_override/DBSNP_SNP.txt";
    $gnomAD             = "$db_path_override/GNOMAD_SNP.txt";
    $Genomes_1000       = "$db_path_override/1000Genome_SNP.txt";
    $GME                = "$db_path_override/GME_SNP.txt";
    $GNOMAD_INDEL       = "$db_path_override/GNOMAD_INDEL.txt";
    $Genomes_1000_INDEL = "$db_path_override/1000Genome_INDEL.txt";
    $GME_INDEL          = "$db_path_override/GME_INDEL.txt";
    $DGV_SV_INS         = "$db_path_override/DGV_SV_INS.txt";
    $Genomes_1000_SV_DEL = "$db_path_override/1000Genome_SV_DEL.txt";
    $DGV_SV_DEL         = "$db_path_override/DGV_SV_DEL.txt";
} else {
    $dbSNP              = $cfg->{databases}->{dbSNP};
    $gnomAD             = $cfg->{databases}->{gnomAD};
    $Genomes_1000       = $cfg->{databases}->{Genomes_1000};
    $GME                = $cfg->{databases}->{GME};
    $GNOMAD_INDEL       = $cfg->{databases}->{GNOMAD_INDEL};
    $Genomes_1000_INDEL = $cfg->{databases}->{Genomes_1000_INDEL};
    $GME_INDEL          = $cfg->{databases}->{GME_INDEL};
    $DGV_SV_INS         = $cfg->{databases}->{DGV_SV_INS};
    $Genomes_1000_SV_DEL = $cfg->{databases}->{Genomes_1000_SV_DEL};
    $DGV_SV_DEL         = $cfg->{databases}->{DGV_SV_DEL};
}

# Process by variant type.
if ($varType eq 'SNP') {
    my $inputFile = &generateInput($inputVcf, $tmpfolder, $varType, $resultfolder);
    my $sample    = 'Input';
    my %DBs       = (
        'dbSNP'       => $dbSNP,
        'gnomAD'      => $gnomAD,
        '1000Genomes' => $Genomes_1000,
        'GME'         => $GME,
    );
    my @db_order  = ('dbSNP','gnomAD','1000Genomes','GME');
    my $totalDB   = 0;
    if ($db ne 'ALL') {
        @db_order = ();
        my %selectedDB = ();
        foreach my $eachInpDB (split /,/, $db) {
            if ( defined($DBs{$eachInpDB}) ) {
                push(@db_order, $eachInpDB);
                $selectedDB{$eachInpDB} = $DBs{$eachInpDB};
            } else {
                print "Error: No database available for <$eachInpDB>\n";
                print "Please see the available $varType databases below:\n\n";
                print_help2();
                exit;
            }
        }
        %DBs = %selectedDB;
        $totalDB = scalar(keys %DBs);
    } else {
        $totalDB = scalar(@db_order);
    }
    
    my $dbNo = 1;
    foreach my $eachDB (@db_order) {
        &compareSNPDB(\$inputFile, $DBs{$eachDB}, $sample, $eachDB, $resultfolder, $dbNo, $totalDB);
        $dbNo++;
    }
    print "The SNP comparison completed successfully and the results were generated in the folder: $resultfolder\n";
}
elsif ($varType eq 'INDEL') {
    my $inputFile = &generateInput($inputVcf, $tmpfolder, $varType, $resultfolder);
    my $sample    = 'Input';
    my %DBs       = (
        'gnomAD'      => $GNOMAD_INDEL,
        '1000Genomes' => $Genomes_1000_INDEL,
        'GME'         => $GME_INDEL,
    );
    my @db_order  = ('gnomAD','1000Genomes','GME');
    my $totalDB   = 0;
    if ($db ne 'ALL') {
        @db_order = ();
        my %selectedDB = ();
        foreach my $eachInpDB (split /,/, $db) {
            if ( defined($DBs{$eachInpDB}) ) {
                push(@db_order, $eachInpDB);
                $selectedDB{$eachInpDB} = $DBs{$eachInpDB};
            } else {
                print "Error: No database available for <$eachInpDB>\n";
                print "Please see the available $varType databases below:\n\n";
                print_help2();
                cleanUp($tmpfolder);
                cleanUp($resultfolder);
                exit;
            }
        }
        %DBs = %selectedDB;
        $totalDB = scalar(keys %DBs);
    } else {
        $totalDB = scalar(@db_order);
    }
    
    my $dbNo = 1;
    foreach my $eachDB (@db_order) {
        &compareINDELDB(\$inputFile, $DBs{$eachDB}, $sample, $eachDB, $resultfolder, $dbNo, $totalDB);
        $dbNo++;
    }
    print "The INDEL comparison completed successfully and the results were generated in the folder: $resultfolder\n";
}
elsif ($varType eq 'SV') {
    my %DBs = (
        'DGV'         => $DGV_SV_DEL,
        '1000Genomes' => $Genomes_1000_SV_DEL,
    );
    my @db_order = ('DGV','1000Genomes');
    my %selectedDB = ();
    if ($db ne 'ALL') {
        @db_order = ();
        foreach my $eachInpDB (split /,/, $db) {
            if ( defined($DBs{$eachInpDB}) ) {
                push(@db_order, $eachInpDB);
                $selectedDB{$eachInpDB} = $DBs{$eachInpDB};
            } else {
                print "Error: No database available for <$eachInpDB>\n";
                print "Please see the available $varType databases below:\n\n";
                print_help2();
                cleanUp($tmpfolder);
                cleanUp($resultfolder);
                exit;
            }
        }
        %DBs = %selectedDB;
    }
    
    if ($pbsv) {
        my ($svCount, $inputINS, $inputDEL, $inputDUP, $inputINV) = &generateInputSV($inputVcf, $tmpfolder, $resultfolder);
        if ($svCount == 0) {
            print "WARNING: There are no SVs present in the given VCF file. Please make sure the input is a pbsv or recommended SV VCF file.\n";
        } else {
            if (-e $inputINS) {
                &compareSVINS($inputINS, $tmpfolder, $resultfolder, $overlap);
            }
            if (-e $inputDEL) {
                &compareSVDEL($inputDEL, $tmpfolder, $resultfolder, $overlap, \@db_order, \%DBs);
            }
        }
    } else {
        my ($svCount, $inputINS, $inputDEL) = &generateInput($inputVcf, $tmpfolder, $varType, $resultfolder);
        if ($svCount == 0) {
            print "WARNING: There are no SVs present in the given VCF file. Please make sure the input is a pbsv or recommended SV VCF file.\n";
        } else {
            if (-e $inputINS) {
                &compareSVINS($inputINS, $tmpfolder, $resultfolder, $overlap);
            }
            if (-e $inputDEL) {
                &compareSVDEL($inputDEL, $tmpfolder, $resultfolder, $overlap, \@db_order, \%DBs);
            }
        }
    }
}

# Clean temporary folder.
cleanUp($tmpfolder);
exit;

#################### SUBROUTINES ###################
sub generateInput {
    my ($infile, $resfolder, $db, $rFolder) = @_;
    my $resultfile = "$resfolder/$db.txt";
    my $resultfile1 = "$resfolder/SV_INS.txt";
    my $resultfile2 = "$resfolder/SV_DEL.txt";
    my $resultfile3 = "$rFolder/SV_Stat.txt";
    
    if ($infile =~ /\.gz$/) {
        open(IN, "gunzip -c $infile |") or die "can't open $infile for reading";
    } else {
        open(IN, "<$infile") or die "Can't open $infile for reading\n";
    }
    open(OUT, ">$resultfile") or die "Can't open $resultfile for writing\n";
    open(OUT1, ">$resultfile1") or die "Can't open $resultfile1 for writing\n";
    open(OUT2, ">$resultfile2") or die "Can't open $resultfile2 for writing\n";
    open(OUT3, ">$resultfile3") or die "Can't open $resultfile3 for writing\n";
    my $totSVs = 0;
    while (<IN>) {
        chomp;
        next if (/^\s*$/);
        next if (/^\s*\#+.*/);
        my @data = split "\t", $_, -1;
        my $chr = $data[0];
        my $pos = $data[1];
        my $ref = $data[3];
        my $l1 = length($ref);
        my $alt = $data[4];
        foreach my $eachAlt (split /,/, $alt, -1) {
            my $l2 = length($eachAlt);
            my $diff = $l2 - $l1;
            if (($l1 == 1) && ($l2 == $l1)) {
                if ($db eq 'SNP') {
                    print OUT "$chr\t$pos\t$ref\t$eachAlt\n";
                }
            } elsif (abs($diff) < 50) {
                if ($db eq 'INDEL') {
                    print OUT "$chr\t$pos\t$ref\t$eachAlt\n";
                }
            } else {
                if ($db eq 'SV') {
                    $chr =~ s/chr//gi;
                    if ($chr eq 'X') { $chr = 23; }
                    if ($chr eq 'Y') { $chr = 24; }
                    if ($diff > 0) {
                        $totSVs++;
                        my $end = $pos + $diff;
                        print OUT1 "$chr\t$pos\t$end\n";
                    } else {
                        $totSVs++;
                        my $end = $pos + abs($diff);
                        print OUT2 "$chr\t$pos\t$end\n";
                    }
                }
            }
        }
    }
    close(IN);
    close(OUT);
    close(OUT1);
    close(OUT2);
    close(OUT3);
    if ($db ne 'SV') {
        unlink $resultfile1;
        unlink $resultfile2;
        unlink $resultfile3;
        return $resultfile;
    } else {
        unlink $resultfile3;
        return ($totSVs, $resultfile1, $resultfile2, $resultfile3);
    }
}

sub compareSNPDB {
    my ($refInfile, $db, $sample, $dbName, $resfolder, $dbNo, $totalDB) = @_;
    my $resultfile1 = ($dbNo == $totalDB)
        ? "${resfolder}/After-comparison-with-${dbName}_Final-Novel_SNPs.txt"
        : "${resfolder}/After-comparison-with-${dbName}_Unique_SNPs.txt";
    my $resultfile2 = "${resfolder}/${sample}_SNP_Stat.txt";
    my %DB_Variants;
    my $total = 0;
    open(IN1, "<$$refInfile") or die "Can't open $$refInfile for reading\n";
    while (<IN1>) {
        chomp;
        next if (/^\s*$/);
        $DB_Variants{$_} = 0;
        $total++;
    }
    close(IN1);
    open(IN2, "<$db") or die "Can't open $db for reading\n";
    while (<IN2>) {
        chomp;
        next if (/^\s*$/);
        if (defined($DB_Variants{$_})) {
            $DB_Variants{$_} = 1;
        }
    }
    close(IN2);
    my $common = 0;
    open(OUT1, ">>$resultfile1") or die "Can't open $resultfile1 for writing\n";
    while (my ($k, $v) = each %DB_Variants) {
        if ($v == 0) {
            print OUT1 "$k\n";
        } else {
            $common++;
        }
    }
    close(OUT1);
    my $unique = $total - $common;
    open(OUT2, ">>$resultfile2") or die "Can't open $resultfile2 for writing\n";
    if ($dbNo == 1) {
        print OUT2 "Database\tTotal input SNPs\tCommon SNPs\tUnique SNPs\n";
    }
    print OUT2 "${dbName}\t$total\t$common\t$unique\n";
    if ($dbNo == $totalDB) {
        print OUT2 "\n\nTotal Novel SNPs : $unique\n";
    }
    close(OUT2);
    $$refInfile = $resultfile1;
}

sub compareINDELDB {
    my ($refInfile, $db, $sample, $dbName, $resfolder, $dbNo, $totalDB) = @_;
    my $resultfile1 = ($dbNo == $totalDB)
        ? "${resfolder}/After-comparison-with-${dbName}_Final-Novel_INDELs.txt"
        : "${resfolder}/After-comparison-with-${dbName}_Unique_INDELs.txt";
    my $resultfile2 = "${resfolder}/${sample}_INDEL_Stat.txt";
    my %DB_Variants;
    my $total = 0;
    open(IN1, "<$$refInfile") or die "Can't open $$refInfile for reading\n";
    while (<IN1>) {
        chomp;
        next if (/^\s*$/);
        $DB_Variants{$_} = 0;
        $total++;
    }
    close(IN1);
    open(IN2, "<$db") or die "Can't open $db for reading\n";
    while (<IN2>) {
        chomp;
        next if (/^\s*$/);
        if (defined($DB_Variants{$_})) {
            $DB_Variants{$_} = 1;
        }
    }
    close(IN2);
    my $common = 0;
    open(OUT1, ">>$resultfile1") or die "Can't open $resultfile1 for writing\n";
    while (my ($k, $v) = each %DB_Variants) {
        if ($v == 0) {
            print OUT1 "$k\n";
        } else {
            $common++;
        }
    }
    close(OUT1);
    my $unique = $total - $common;
    open(OUT2, ">>$resultfile2") or die "Can't open $resultfile2 for writing\n";
    if ($dbNo == 1) {
        print OUT2 "Database\tTotal input INDELs\tCommon INDELs\tUnique INDELs\n";
    }
    print OUT2 "${dbName}\t$total\t$common\t$unique\n";
    if ($dbNo == $totalDB) {
        print OUT2 "\n\nTotal Novel INDELs : $unique\n";
    }
    close(OUT2);
    $$refInfile = $resultfile1;
}

sub generateInputSV {
    my ($infile, $resFolder, $resFolder2) = @_;
    my $sample = 'Input';
    my $resultfile1 = "$resFolder/${sample}_SV_INS.txt";
    my $resultfile2 = "$resFolder/${sample}_SV_DEL.txt";
    my $resultfile3 = "$resFolder/${sample}_SV_DUP.txt";
    my $resultfile4 = "$resFolder/${sample}_SV_INV.txt";
    my $resultfile5 = "$resFolder2/${sample}_SV_Stat.txt";
    my %VarInfo;
    my @header;
    open(OUT1, ">$resultfile1") or die "Can't open $resultfile1 for writing\n";
    open(OUT2, ">$resultfile2") or die "Can't open $resultfile2 for writing\n";
    open(OUT3, ">$resultfile3") or die "Can't open $resultfile3 for writing\n";
    open(OUT4, ">$resultfile4") or die "Can't open $resultfile4 for writing\n";
    open(IN, "<$infile") or die "Can't open $infile for reading\n";
    while (<IN>) {
        chomp;
        next if (/^\s*$/);
        if (/^\s*\#CHR.*/) {
            @header = split '\t', $_, -1;
        }
        next if (/^\s*\#+.*/);
        my @data = split "\t", $_, -1;
        if ($data[6] eq 'PASS') {
            my ($SVType, $SVLen, $end) = ('NA', 'NA', 'NA');
            foreach my $eachInfo (split /;/, $data[7], -1) {
                my ($k, $v) = split /=/, $eachInfo, -1;
                if ($k eq 'SVTYPE') {
                    $SVType = $v;
                } elsif ($k eq 'END') {
                    $end = $v;
                } elsif ($k eq 'SVLEN') {
                    $SVLen = $v;
                }
            }
            my $chr = $data[0];
            $chr =~ s/chr//gi;
            if ($chr eq 'X') { $chr = 23; }
            if ($chr eq 'Y') { $chr = 24; }
            if ($SVType eq 'INS') {
                $end = $data[1] + $SVLen;
                if (length($chr) == 1) {
                    print OUT1 "$chr\t$data[1]\t$end\n";
                }
            } elsif ($SVType eq 'DEL') {
                if (length($chr) == 1) {
                    print OUT2 "$chr\t$data[1]\t$end\n";
                }
            } elsif ($SVType eq 'DUP') {
                print OUT3 "$chr\t$data[1]\t$end\n";
            } elsif ($SVType eq 'INV') {
                print OUT4 "$chr\t$data[1]\t$end\n";
            }
            foreach my $i (9 .. $#header) {
                my $gt = (split /:/, $data[$i], -1)[0];
                if (($gt eq '1/0') || ($gt eq '1/1') || ($gt eq '0/1') || ($gt eq './1') || ($gt eq '1/.')) {
                    $VarInfo{$header[$i]}{'TYPE'}{$SVType}++;
                    $VarInfo{$header[$i]}{'TYPE'}{'TOTAL'}++;
                }
            }
        }
    }
    close(IN);
    close(OUT1);
    close(OUT2);
    close(OUT3);
    close(OUT4);
    open(OUT5, ">$resultfile5") or die "Can't open $resultfile5 for writing\n";
    my @SVTypes = ('INS', 'DEL', 'DUP', 'INV', 'BND');
    print OUT5 "Sample\tTotal Variants\t" . (join "\t", @SVTypes) . "\n";
    my $totSVs = 0;
    foreach my $i (9 .. $#header) {
        my $sample = $header[$i];
        my @sample_result = ($sample);
        foreach my $eachType ('TOTAL', @SVTypes) {
            if (defined($VarInfo{$header[$i]}{'TYPE'}{$eachType})) {
                $totSVs++;
                push(@sample_result, $VarInfo{$header[$i]}{'TYPE'}{$eachType});
            } else {
                push(@sample_result, 0);
            }
        }
        my $eachResult = join "\t", @sample_result;
        print OUT5 "$eachResult\n";
    }
    close(OUT5);
    if (-z $resultfile1) { unlink $resultfile1; }
    if (-z $resultfile2) { unlink $resultfile2; }
    if (-z $resultfile3) { unlink $resultfile3; }
    if (-z $resultfile4) { unlink $resultfile4; }
    return ($totSVs, $resultfile1, $resultfile2, $resultfile3, $resultfile4);
}

sub compareSVINS {
    my ($infile, $tmpfolder, $resultfolder, $svOverlap) = @_;
    my $resultfile = "$tmpfolder/After-comparison-with-DGV_Final-Unique-SVs.txt";
    my $retval = system("java $cfg->{classes}->{reciprocal_overlap} $infile $DGV_SV_INS $svOverlap > $resultfile");
    if ($retval == 0) {
        &parseSVINS($resultfile, $resultfolder);
        print "SV insertion comparison completed!!\n";
    } else {
        print "Problem in the reciprocal overlap\n";
        exit;
    }
}

sub parseSVINS {
    my ($infile, $resultfolder) = @_;
    my $resultfile1 = "$resultfolder/After-comparison-with-DGV_Final-Unique-SV-Insertions.txt";
    my $resultfile2 = "$resultfolder/SV-Insetion_Comparison_Stat.txt";
    my $total = 0;
    my $common = 0;
    my $unique = 0;
    open(OUT1, ">$resultfile1") or die "Can't open $resultfile1 for writing\n";
    open(IN, "<$infile") or die "Can't open $infile for reading\n";
    while (<IN>) {
        chomp;
        next if (/^\s*$/);
        $total++;
        my @data = split "\t", $_, -1;
        if ($data[1] > 0) {
            $common++;
        } elsif ($data[1] == 0) {
            my ($chr, $cord) = split /:/, $data[0], -1;
            my ($cord1, $cord2) = split /-/, $cord, -1;
            $chr =~ tr/Chr//d;
            print OUT1 "$chr\t$cord1\t$cord2\n";
            $unique++;
        }
    }
    close(IN);
    close(OUT1);
    open(OUT2, ">$resultfile2") or die "Can't open $resultfile2 for writing\n";
    print OUT2 "Total insertions\tCommon insertions\tUnique insertions\n";
    print OUT2 "$total\t$common\t$unique\n";
    print OUT2 "\n\nTotal Novel SV Insertions : $unique\n";
    close(OUT2);
}

sub compareSVDEL {
    my ($infile, $tmpfolder, $resultfolder, $svOverlap, $refDBorder, $refDBs) = @_;
    my $totalDBs = scalar(@$refDBorder);
    my $no = 1;
    my $retFile = '';
    foreach my $eachDb (@$refDBorder) {
        my $dbfile = $$refDBs{$eachDb};
        my $resultfile = "$tmpfolder/After-comparison-with-${eachDb}-SV-Deletions.txt";
        if ($no == $totalDBs) {
            $resultfile = "$tmpfolder/After-comparison-with-${eachDb}_Final-Unique-SV-Deletions.txt";
        }
        my $retval = 'NA';
        if ($no == 1) {
            $retval = system("java $cfg->{classes}->{reciprocal_overlap} $infile $dbfile $svOverlap > $resultfile");
        } else {
            $retval = system("java $cfg->{classes}->{reciprocal_overlap} $retFile $dbfile $svOverlap > $resultfile");
        }
        if ($retval == 0) {
            $retFile = &parseSVDEL($resultfile, $resultfolder, $eachDb, $totalDBs, $no);
            if ($no == $totalDBs) {
                print "SV Deletion comparison completed!!\n";
            }
        } else {
            print "Problem in the reciprocal overlap\n";
            exit;
        }
        $no++;
    }
}

sub parseSVDEL {
    my ($infile, $resultfolder, $db, $totalDB, $N) = @_;
    my $resultfile1 = "${resultfolder}/After-comparison-with-${db}_SV-Deletions.txt";
    if ($N == $totalDB) {
        $resultfile1 = "${resultfolder}/After-comparison-with-${db}_Final-Unique-SV-Deletions.txt";
    }
    my $resultfile2 = "${resultfolder}/SV-Deletion_Comparison_Stat.txt";
    my $total = 0;
    my $common = 0;
    my $unique = 0;
    open(OUT1, ">$resultfile1") or die "Can't open $resultfile1 for writing\n";
    open(IN, "<$infile") or die "Can't open $infile for reading\n";
    while (<IN>) {
        chomp;
        next if (/^\s*$/);
        $total++;
        my @data = split "\t", $_, -1;
        if ($data[1] > 0) {
            $common++;
        } elsif ($data[1] == 0) {
            my ($chr, $cord) = split /:/, $data[0], -1;
            my ($cord1, $cord2) = split /-/, $cord, -1;
            $chr =~ tr/Chr//d;
            print OUT1 "$chr\t$cord1\t$cord2\n";
            $unique++;
        }
    }
    close(IN);
    close(OUT1);
    if ($N == 1) {
        open(OUT2, ">$resultfile2") or die "Can't open $resultfile2 for writing\n";
        print OUT2 "Database\tTotal deletions\tCommon deletions\tUnique deletions\n";
    } else {
        open(OUT2, ">>$resultfile2") or die "Can't open $resultfile2 for writing\n";
    }
    print OUT2 "${db}\t$total\t$common\t$unique\n";
    if ($N == $totalDB) {
        print OUT2 "\n\nTotal Novel SV Deletions : $unique\n";
    }
    close(OUT2);
    return $resultfile1;
}
