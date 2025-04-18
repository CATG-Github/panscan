#!/usr/bin/perl -w
use strict;
use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use YAML::XS 'LoadFile';
use Cwd qw(abs_path);
use perlModules::panscan qw(
    cleanUp 
    print_uniq_var_help 
    generateInput 
    compareSNPDB 
    compareINDELDB 
    generateInputSV 
    compareSVINS 
    parseSVINS 
    compareSVDEL 
    parseSVDEL
);

# Specify the path to the config file relative to this script’s directory.
my $config_file = "$FindBin::Bin/config.yaml";
if (! -e $config_file) {
    die("Config file '$config_file' not found.\n");
}
my $cfg = LoadFile($config_file);
if(!$cfg) {
    die("Failed to read YAML file: $config_file\n");
}

# Variable Declaration
my $inputVcf = 'NA';
my $varType  = 'NA';
my $db       = 'ALL';
my $overlap  = 80;
my $help     = undef;
my $pbsv     = 0;
my $output   = '';   # New output directory flag
my $db_path_override = '';  # New database base path override flag
my %types    = (
    'SNP'   => 'SNP',
    'INDEL' => 'INDEL',
    'SV'    => 'SV'
);

GetOptions(
    'i=s'       => \$inputVcf,
    't=s'       => \$varType,
    'pbsv'      => \$pbsv,
    'db=s'      => \$db,
    'overlap=i' => \$overlap,
    'op=s'      => \$output,         # New output option
    'db-path=s' => \$db_path_override, # New database path override option
    'help'      => \$help
) or die "Error in command line arguments\n";

if ($help) {
    print_uniq_var_help();
    exit;
}

$varType = uc($varType);
if ($varType eq 'NA') {
    print "\nERROR: Please provide variant type to compare (SNP/INDEL/SV)!!\n\n";
    print_uniq_var_help();
    exit;
} elsif (!defined($types{$varType})) {
    print "\nERROR: Unknown variant type <$varType>. Please provide variant type to compare (SNP/INDEL/SV)!!\n\n";
    print_uniq_var_help();
    exit;
}

if ($inputVcf eq 'NA') {
    print "ERROR: Please provide input VCF file!!\n";
    print_uniq_var_help();
    exit;
}

# Determine result folder.
my $resultfolder;
if ($output) {
    $resultfolder = abs_path($output) or die "Cannot resolve output directory: $output\n";
    mkdir $resultfolder unless (-d $resultfolder);
} else {
    $resultfolder = "${varType}_Comparison_Result";
    mkdir $resultfolder unless (-d $resultfolder);
}

# Temporary folder.
my $tmpfolder = 'tmp_' . $varType;
unless (-d $tmpfolder) {
    mkdir $tmpfolder;
} else {
    cleanUp($tmpfolder);
    mkdir $tmpfolder;
}

# Database files.
my ($dbSNP, $gnomAD, $Genomes_1000, $GME,
    $GNOMAD_INDEL, $Genomes_1000_INDEL, $GME_INDEL,
    $Genomes_1000_SV_DEL, $DGV_SV_DEL);

if ($db_path_override) {
    $dbSNP            = "$db_path_override/DBSNP_SNP.txt";
    $gnomAD           = "$db_path_override/GNOMAD_SNP.txt";
    $Genomes_1000     = "$db_path_override/1000Genome_SNP.txt";
    $GME              = "$db_path_override/GME_SNP.txt";
    $GNOMAD_INDEL     = "$db_path_override/GNOMAD_INDEL.txt";
    $Genomes_1000_INDEL = "$db_path_override/1000Genome_INDEL.txt";
    $GME_INDEL        = "$db_path_override/GME_INDEL.txt";
    $Genomes_1000_SV_DEL = "$db_path_override/1000Genome_SV_DEL.txt";
    $DGV_SV_DEL       = "$db_path_override/DGV_SV_DEL.txt";
} else {
    $dbSNP            = "$cfg->{databases}->{dbpath}/DBSNP_SNP.txt";
    $gnomAD           = "$cfg->{databases}->{dbpath}/GNOMAD_SNP.txt";
    $Genomes_1000     = "$cfg->{databases}->{dbpath}/1000Genome_SNP.txt";
    $GME              = "$cfg->{databases}->{dbpath}/GME_SNP.txt";
    $GNOMAD_INDEL     = "$cfg->{databases}->{dbpath}/GNOMAD_INDEL.txt";
    $Genomes_1000_INDEL = "$cfg->{databases}->{dbpath}/1000Genome_INDEL.txt";
    $GME_INDEL        = "$cfg->{databases}->{dbpath}/GME_INDEL.txt";
    $Genomes_1000_SV_DEL = "$cfg->{databases}->{dbpath}/1000Genome_SV_DEL.txt";
    $DGV_SV_DEL       = "$cfg->{databases}->{dbpath}/DGV_SV_DEL.txt";
}

if($varType eq 'SNP') {
    my $inputFile = generateInput($inputVcf, $tmpfolder, $varType, $resultfolder);
    my $sample = 'Input-cohort';
    my %DBs = (
        'dbSNP'       => $dbSNP,
        'gnomAD'      => $gnomAD,
        '1000Genomes' => $Genomes_1000,
        'GME'         => $GME
    );
    my @db_order = ('dbSNP','gnomAD','1000Genomes','GME');
    my $totalDB = ($db ne 'ALL')
      ? do {
            my %selectedDB;
            @db_order = ();
            foreach my $eachInpDB(split /,/, $db) {
                if(defined($DBs{$eachInpDB})) {
                    push(@db_order, $eachInpDB);
                    $selectedDB{$eachInpDB} = $DBs{$eachInpDB};
                } else {
                    print "Error: No database available for <$eachInpDB>\n";
                    print "Please see the available $varType databases below:\n\n";
                    print_uniq_var_help();
                    exit;
                }
            }
            %DBs = %selectedDB;
            scalar(keys %DBs);
        }
      : scalar(@db_order);
    my $dbNo = 1;
    foreach my $eachDB (@db_order) {
        compareSNPDB(\$inputFile, $DBs{$eachDB}, $sample, $eachDB, $resultfolder, $dbNo, $totalDB);
        $dbNo++;
    }
    print "The SNP comparison completed successfully and the results were generated in the folder: $resultfolder\n";
} elsif ($varType eq 'INDEL') {
    my $inputFile = generateInput($inputVcf, $tmpfolder, $varType, $resultfolder);
    my $sample = 'Input-cohort';
    my %DBs = (
        'gnomAD'      => $GNOMAD_INDEL,
        '1000Genomes' => $Genomes_1000_INDEL,
        'GME'         => $GME_INDEL
    );
    my @db_order = ('gnomAD','1000Genomes','GME');
    my $totalDB = ($db ne 'ALL')
      ? do {
            my %selectedDB;
            @db_order = ();
            foreach my $eachInpDB(split /,/, $db) {
                if(defined($DBs{$eachInpDB})) {
                    push(@db_order, $eachInpDB);
                    $selectedDB{$eachInpDB} = $DBs{$eachInpDB};
                } else {
                    print "Error: No database available for <$eachInpDB>\n";
                    print "Please see the available $varType databases below:\n\n";
                    print_uniq_var_help();
                    cleanUp($tmpfolder);
                    cleanUp($resultfolder);
                    exit;
                }
            }
            %DBs = %selectedDB;
            scalar(keys %DBs);
        }
      : scalar(@db_order);
    my $dbNo = 1;
    foreach my $eachDB (@db_order) {
        compareINDELDB(\$inputFile, $DBs{$eachDB}, $sample, $eachDB, $resultfolder, $dbNo, $totalDB);
        $dbNo++;
    }
    print "The INDEL comparison completed successfully and the results were generated in the folder: $resultfolder\n";
} elsif ($varType eq 'SV') {
    my %DBs = (
        'DGV'         => $DGV_SV_DEL,
        '1000Genomes' => $Genomes_1000_SV_DEL
    );
    my @db_order = ('DGV','1000Genomes');
    if ($db ne 'ALL') {
        my %selectedDB;
        @db_order = ();
        foreach my $eachInpDB(split /,/, $db) {
            if (defined($DBs{$eachInpDB})) {
                push(@db_order, $eachInpDB);
                $selectedDB{$eachInpDB} = $DBs{$eachInpDB};
            } else {
                print "Error: No database available for <$eachInpDB>\n";
                print "Please see the available $varType databases below:\n\n";
                print_uniq_var_help();
                cleanUp($tmpfolder);
                cleanUp($resultfolder);
                exit;
            }
        }
        %DBs = %selectedDB;
    }
    
    if ($pbsv) {
        my ($svCount, $inputINS, $inputDEL, $inputDUP, $inputINV) = generateInputSV($inputVcf, $tmpfolder, $resultfolder);
        if ($svCount == 0) {
            print "WARNING: There are no SVs present in the given VCF file. Please make sure the input is a pbsv or recommended SV VCF file.\n";
        } else {
            if (-e $inputINS) {
                compareSVINS($inputINS, $tmpfolder, $resultfolder, $overlap, $cfg);
            }
            if (-e $inputDEL) {
                compareSVDEL($inputDEL, $tmpfolder, $resultfolder, $overlap, \@db_order, \%DBs, $cfg);
            }
        }
    } else {
        my ($svCount, $inputINS, $inputDEL) = generateInput($inputVcf, $tmpfolder, $varType, $resultfolder);
        if ($svCount == 0) {
            print "WARNING: There are no SVs present in the given VCF file. Please make sure the input is a pbsv or recommended SV VCF file.\n";
        } else {
            if (-e $inputINS) {
                compareSVINS($inputINS, $tmpfolder, $resultfolder, $overlap);
            }
            if (-e $inputDEL) {
                compareSVDEL($inputDEL, $tmpfolder, $resultfolder, $overlap, \@db_order, \%DBs);
            }
        }
    }
}

cleanUp($tmpfolder);
exit;
