package perlModules::panscan;

use strict;
use warnings;
use Exporter 'import';

our @EXPORT_OK = qw(cleanUp print_help1 print_help2);


sub cleanUp
{
        my $infolder = shift;
	foreach my $eachfiles(glob("$infolder\/*"))
	{
		unlink($eachfiles);
	}
        rmdir $infolder;
}


sub print_help1
{
print <<EOF;
        Usage: perl preprocessVCF.pl [options]

        Options:
        --i FILE       Specify the input vcf files.
        --help         Display this help message.

EOF
}


sub print_help2
{
print <<EOF;
        Usage: perl findUniqVariants.pl [options]

        Options:
        --i FILE       Specify the input vcf files. ( For SV comparison, input should be pbsv vcf file )
        --t TEXT       Specify the variant type to compare (SNP/INDEL/SV).
        --db FILE      Specify the databases ( ',' sepearted for multiple databases ).
                       Available databases for
                                        SNP   : dbSNP,gnomAD,1000Genomes,GME
                                        InDel : gnomAD,1000Genomes,GME
                                        SV insertion    : DGV
                                        SV deletion     : DGV,1000Genomes
				default : ALL
        --op NUMBER     Specify the percentage of overlap for the SV comparison.
			default : 80
        --help         Display this help message.

EOF
}


1;
