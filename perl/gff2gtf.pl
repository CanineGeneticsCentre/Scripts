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
#use warnings::unused;
use Term::ANSIColor;

# VERSION OF SOFTWARE #
my $version						= "1";


####################
# Define variables #
####################

my $line_count					= 0;
my $array_count					= 0;
my $array_size					= 0;
my $total_no_lines				= 0;

my $mystring					= "";
my $input_file					= "";
my $output_file					= "";
my $output_file_sorted			= "";
my $output_file_unsorted		= "";
my $single_line					= "";
my $answer						= "";
my $snp_type					= "";
my $snp_start					= "";
my $snp_end						= "";
my $snp_id						= "";
my $snp_alleles					= "";
my $chromosome					= "";
my $source						= "";
my $ref_allele					= "";
my $alt_allele					= "";
my $command						= "";
my $prefix						= "";

my $feature						= "";
my $start						= "";
my $end							= "";
my $attribute					= "";
my $score						= "";
my $frame						= "";
my $strand						= "";
my $tag							= "";
my $value						= "";
my $tag_value_out				= "";
my $first_column				= "";

my @item						= ();
my @attribute_item				= ();
my @tag_value					= ();

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

print "  - This program processes GFF files into GTF files.\n\n";


##########################
# Get name of input file #
##########################

print "Name of input file:  ";
$input_file = <STDIN>;
chomp $input_file;

$prefix = &get_prefix ($input_file);

print "Name of output file:  ";
$output_file = <STDIN>;
chomp $output_file;

if ($output_file eq ""){$output_file = "$prefix".".vcf";}

$output_file_unsorted = "$prefix"."_unsorted.vcf";
$output_file_sorted = $output_file;

print "Chromosome: (i.e. what you want in the first column) ";
$chromosome = <STDIN>;
chomp $chromosome;

if (index($chromosome,"chr") == -1)
{
	$chromosome = "chr"."$chromosome";
}
print "Chromosome is $chromosome\n\n";

####################################################
# Open the list file to get the list of file names #
####################################################
open (IN, "$input_file") || die "Cannot open $input_file";
open (OUT, ">$output_file_unsorted")|| die "Cannot create output file: $output_file_unsorted";


print "Input file:\t$input_file\n";
print "Unsorted output file:\t$output_file_unsorted\n";
print "Sorted output file:\t$output_file\n";

print "\nPress 'return' to continue:  ";
$answer=<STDIN>;

$line_count=0;

while ($single_line = <IN> ) 
{
	$line_count = $line_count + 1;
	
	chomp $single_line;
	
	# Split line at TABS (should make or spaces)
	@item=split(/\s+/,$single_line);
	
	$array_size = scalar @item;
	
	$first_column = $item[0];
	$source = $item[1];
	$feature = $item[2];
	$start = $item[3];
	$end = $item[4];
	$score = $item[5];
	$strand = $item[6];
	$frame = $item[7];
	$attribute = $item[8];
		
	# Split attribute at semicolons 

	if (substr($first_column,0,1) ne "#")
	{	
		print OUT "$chromosome\t$source\t$feature\t$start\t$end\t$score\t$strand\t$frame\t";
		
		print "Line: $line_count\n";
		
		@attribute_item = split(";",$attribute);
		
		$array_size = (scalar @attribute_item) - 1;
	
		for ($array_count = 0; $array_count <= $array_size; $array_count++)
		{	
			@tag_value = split ("=",$attribute_item[$array_count]);
			
			$tag = $tag_value[0];
			$value = $tag_value[1];
			
			$tag_value_out = "$tag \"$value\"";
			
			print OUT $tag_value_out;
			print OUT ";";
		}
	
		# end of line
		print OUT "\n";
		
	}
	
	#$answer=<STDIN>;

	
} # end of while
	
$total_no_lines = $line_count;
	
	
####################
# Sort output file #
####################
print "\nSorting...\n\n\n";

$command = "(head -n 2 $output_file_unsorted; tail -n +3 $output_file_unsorted | sort -k2 -n) > $output_file";

print "$command\n";
system ("$command");

close OUT;
close IN;


print "\n\n";
print "############\n";
print "# Finished #\n";
print "############\n\n\n";

print "Chromosome  \t$chromosome\n\n";

print "Unsorted file:\t$output_file_unsorted\n";
print "Sorted file:\t$output_file_sorted\n\n";

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
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	