#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	GTF2VCF           						                        #     
#									                                    #
#	Converts GTF files to VCF format                   	                        #
#									                                    #
#########################################################################

#############################
# Mike Boursnell            #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# mike.boursnell@aht.org.uk #
#############################


use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;

# VERSION OF SOFTWARE #
my $version						= "3";


####################
# Define variables #
####################

my $line_count					= 0;
my $array_size					= 0;
my $total_no_lines				= 0;

my $input_file					= "";
my $output_file					= "";
my $output_file_unsorted		= "";
my $single_line					= "";
my $answer						= "";
my $snp_type					= "";
my $snp_start					= "";
my $snp_end						= "";
my $snp_id						= "";
my $snp_alleles					= "";
my $chromosome					= "";
my $ref_allele					= "";
my $alt_allele					= "";
my $command						= "";
my $prefix						= "";

my @item						= ();


print color 'reset';
print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      GTF2VCF       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";

print "  - This program processes GTF files into VCF files.\n\n";

print color 'reset';

##########################
# Get name of input file #
##########################

until (-e $input_file)
{
	print "\nPlease input the name of your GTF file to be converted to VCF:      ";
	$input_file = <STDIN>;chomp $input_file;
	
	if ($input_file eq "ls"){print "\n";system ("ls *.gtf")}
	
	if ($input_file ne "ls")
	{
		if (! -e $input_file){print "\n\n>>>>>>>>  File $input_file not found.  Try again.  <<<<<<<<\n\n";}
	}
}


$prefix = &get_prefix ($input_file);
print "\nName of output file: (default= $prefix.vcf):                      ";

$output_file = <STDIN>;
chomp $output_file;

if ($output_file eq ""){$output_file = "$prefix".".vcf";}

$output_file_unsorted = "$prefix"."_unsorted.vcf";


until ($chromosome ne "")
{
	print "\nChromosome:  ";
	$chromosome = <STDIN>;
	chomp $chromosome;

	if ($chromosome eq "")
	{
		print "\n >>>>>>>   You must enter a chromosome number here!  (which must match the chromosome in the filename)  <<<<<<<<<\n"
	}
	if ((index($chromosome,"chr") == -1) && ($chromosome ne ""))
	{
		$chromosome = "chr"."$chromosome";
	}
}

print "Chromosome is $chromosome\n\n";

####################################################
# Open the list file to get the list of file names #
####################################################
open (IN, "$input_file") || die "Cannot open $input_file";
open (OUT, ">$output_file_unsorted")|| die "Cannot create output file: $output_file_unsorted";

print OUT "##fileformat=VCFv4.0\n";
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n";

print "Input file:            \t$input_file\n";
print "Unsorted output file:  \t$output_file_unsorted\n";
print "Sorted output file:    \t$output_file\n\n";

print color 'yellow';
print "Press 'return' to continue:  ";
$answer=<STDIN>;
print color 'reset';
$line_count=0;

while ($single_line = <IN> ) 
{
	$line_count = $line_count + 1;
	
	chomp $single_line;
	
	# Split line at TABS (should make or spaces)
	@item=split(/\s+/,$single_line);
	
	$array_size = scalar @item;
	
	$snp_type = $item[2];
	$snp_start = $item[3];
	$snp_end = $item[4];
	$snp_id = $item[9];
	$snp_alleles = $item[11];
	
	$snp_id = substr $snp_id, 1;
	$snp_id = substr $snp_id, 0, length($snp_id)-2;
	
	$ref_allele = substr $snp_alleles, 2, 1;
	$alt_allele = substr $snp_alleles, 1, 1;
	
	print "Line: $line_count\t$snp_alleles\t$ref_allele\t$alt_allele\n";
	
	#$answer = <STDIN>;
	############################
	# Now write to output file #
	############################
	
	#CHROM POS ID REF ALT QUAL FILTER INFO
	print OUT "$chromosome\t$snp_start\t$snp_id\t$ref_allele\t$alt_allele\t.\t.\t.\n";
	
} # end of while
	
$total_no_lines = $line_count;
	
	
####################
# Sort output file #
####################
print "\n\n\nSorting...\n\n\n";

$command = "(head -n 2 $output_file_unsorted; tail -n +3 $output_file_unsorted | sort -k2 -n) > $output_file";

print "$command\n";
system ("$command");

close OUT;
close IN;
close SORT;

print "\n\n";
print color 'bold yellow';
print "############\n";
print "# Finished #\n";
print "############\n\n\n";

print "Chromosome  \t$chromosome\n\n";

print "Unsorted file:\t$output_file_unsorted\n\n";
print "Sorted file:\t$output_file\n\n\n";
print color 'reset';
exit;


sub get_prefix
{
	my $filename = "";

	$filename = $_[0];
	if (index($filename,".") > 0)
	{
		$filename = substr($filename, 0, index($filename,"."));
	}
	if (index($filename,".") == -1)
	{
		$filename = $filename;
	}
}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	