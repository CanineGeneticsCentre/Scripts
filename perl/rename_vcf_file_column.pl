#!/usr/bin/perl -w

##################################################################################
#									                                             #      
#	rename_vcf_file_column			                                             #     
#									                                             #
#	This PERL script changes the name in the column heading                 	 #
#									                                             #
##################################################################################

use strict;
use Term::ANSIColor;
use List::Util qw[min max];

#Strings
my $version						= "1";
my $vcf_file					= "";
my $vcf_output_file				= "";
my $vcf_line					= "";
my $prefix						= "";
my $answer						= "";
my $old_name					= "";
my $new_name					= "";
my $chromosome					= "";
my $position					= "";

#Numbers
my $no_of_samples				= 0;
my $no_of_samples_chrom_line	= 0;
my $chrom_line_array_size		= 0;

#Counters
my $array_count					= 0;
my $line_count					= 0;
my $column_count				= 0;

my $sample_count				= 0;

#Boolean
my $passed_header_lines			= ""; # true or false

#Arrays
my @chrom_line_array			= ();
my @vcf_line_array				= ();
my @sample_name_array			= ();
my @sample_name_array_new		= ();

#############################################################
print color 'reset';
print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "    rename VCF file column            \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";


print color 'yellow';
print "Version $version\n\n";
print "  - This PERL script changes the name of the column header in a VCF file\n\n";

###########################################################
# Get the name of the input VCF file                      #
###########################################################
$vcf_file = &get_file ("What is the name of the input VCF file (or .gVCF file)","vcf");

$prefix = &get_prefix($vcf_file);

if (index ($vcf_file,".vcf") > -1){$vcf_output_file = $prefix."_new.vcf";}
if (index ($vcf_file,".gVCF") > -1){$vcf_output_file = $prefix."_new.gVCF";}


###########################################################
# Open the input VCF file                                 #
# This is the MAIN LOOP of the script                     #
###########################################################
open (VCF, "$vcf_file") || die "Cannot open $vcf_file";

###########################################################
# Open the output tVCF file                               #
###########################################################
open (OUT, ">$vcf_output_file") || die "Cannot open $vcf_output_file";


while ($vcf_line = <VCF>) 
{
	chomp $vcf_line;
	$line_count = $line_count + 1;

	if (index($vcf_line,"#CHROM") > -1)
	{
		############################################################
		# Parsing of #CHROM data line to get a list of input files #
		############################################################
		##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Dennis5runs_merged

		$passed_header_lines = "true";

		@chrom_line_array = split(/\s+/,$vcf_line);
		$chrom_line_array_size = (scalar @chrom_line_array) - 1;
		$no_of_samples_chrom_line = $chrom_line_array_size - 8; # ignore first columns 0-9

		for ($array_count = 9; $array_count <= $chrom_line_array_size; $array_count++)
		{
			# Sample names go into sample_name_array_input_order
			$sample_name_array[$array_count - 8] = $chrom_line_array[$array_count];
			$column_count = $array_count - 8;
		}
		
		
        ############################################################
        # Total number of samples is defined by the number of      #
        # samples on the #CHROM line, $no_of_samples_chrom_line    #
        ############################################################
		$no_of_samples = $column_count;


		print "There are $no_of_samples columns in this VCF file\n\n";

		print "Enter the new names [press 'return' to keep the name unchanged]\n\n";

		for ($sample_count = 1; $sample_count <=$no_of_samples; $sample_count++)
		{
			print "\nSample $sample_count.\n";
			print "\tExisting sample name: $sample_name_array[$sample_count]  ";

			print "\tNew sample name:     >";
			$answer=<STDIN>;
			chomp $answer;

			if ($answer ne "")
			{
				$sample_name_array_new[$sample_count] = $answer;
			}
			else
			{
				$sample_name_array_new[$sample_count] = $sample_name_array[$sample_count];
			}
		}

		######################################
		# Print out first bit of #CHROM line #
		######################################
		print OUT "#CHROM      POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  ";

		######################################
		# Check changes                      #
		######################################

		&print_message("Check these changes carefully","message");

		print "Old name\tNew name\n";
		print "========\t========\n";
		for ($sample_count = 1; $sample_count <=$no_of_samples; $sample_count++)
		{
			# pad out disease status to 12 characters
			$old_name = $sample_name_array[$sample_count];
			$new_name = $sample_name_array_new[$sample_count];

			$old_name = sprintf("%-*s", 20, $old_name);
			$new_name = sprintf("%-*s", 20, $new_name);

			print "$old_name\t$new_name\n";
		}

		print "\nAre these changes OK (y/n)  ";
		$answer=<STDIN>;
		chomp $answer;

		if (lc $answer ne "y"){exit;}

		print "\n\n";

		######################################
		# Print out sample headers           #
		######################################
		for ($sample_count = 1; $sample_count <=$no_of_samples; $sample_count++)
		{
			print OUT "$sample_name_array_new[$sample_count]   ";
		}

		print OUT "\n";
	} # chrom line
	else
	{
		print OUT "$vcf_line\n";
	}

	if (($line_count % 1000000) == 0)
	{
		if ($passed_header_lines eq "true")
		{
			@vcf_line_array = split(/\s+/,$vcf_line);
			$chromosome = $vcf_line_array[0];
			$position = $vcf_line_array[1];

			print "Writing new VCF file.  Line: $line_count\tChr: $chromosome\tPos: $position\n";
		}
	}

} # while loop

close VCF;
close OUT;

&print_message ("Renaming completed","message");

print "Input file:          \t$vcf_file\n";
print "Output file renamed: \t$vcf_output_file\n\n";

exit;


#################################################################
# Subroutine to get a file name                                 #
# Argument 1: Text that user sees, asking for the file          #
# Argument 2: suffix of files to search for with ls e.g. ".bed" #
#################################################################
sub get_file
{
	my $input_message 	= $_[0];
	my $_file_type	 	= $_[1];
	my $_file_name		= "";
	my $_search_string	= "";

	if ($_file_type eq ""){$_file_type = "*"}

	until ((-e $_file_name) || (lc $_file_name eq "q"))
	{
		&print_message("$input_message","input");
		print "> ";

		$_file_name = <STDIN>;
		chomp $_file_name;

		# User types 'ls'
		if (($_file_name eq "ls") || ($_file_name eq "")) {print "\n";system ("ls *.vcf");system ("ls *.gVCF");}

		# Starts with 'ls' followed by search string
		if (($_file_name ne "ls") && (index ($_file_name,"ls") == 0))
		{
			$_search_string = substr($_file_name,3,99);
			print "\n";
			system ("ls *$_search_string*");
			$_file_name = "ls";
			print "\n";
		}

		if (($_file_name ne "ls")  && (lc $_file_name ne "q") && ($_file_name ne ""))
		{

			if (!-e $_file_name){print "\n  >>>>>>>>>>>>>  ERROR.  File $_file_name cannot be found <<<<<<<<<<<\n\n";}

			if (-e $_file_name)
			{
				$_file_name = $_file_name;
			}
		} # not ls

	} # until loop

	print "\n";
	if ($_file_type eq ".txt"){system("dos2unix $_file_name")}

	$_file_name = $_file_name;
} # get_file



######################################
# Subroutine to print screen message #
######################################
sub print_message
{
	my $_message_length 	= "";
	my $_pos_count		= 0;
	my $_char			= "";
	
	my $_message = $_[0];
	my $_style = $_[1];
	
	$_message_length = length($_message);
	
	if ($_style eq ""){$_char = "#"}
	if ($_style eq "input"){$_char = "~"}
	if ($_style eq "message"){$_char = "#"}
	if ($_style eq "warning"){$_char = "!"}
	
	print "\n\n";
	print color ' bold yellow';
	if ($_style eq "warning"){print color ' bold red'}
	if ($_style eq "input"){print color ' bold white'}
	
	for ($_pos_count = 1;$_pos_count <=($_message_length + 10);$_pos_count++){print $_char}
	
	print "\n$_char    $_message    $_char\n";
	
	for ($_pos_count = 1;$_pos_count <=($_message_length + 10);$_pos_count++){print $_char}
	
	print "\n\n";
	print color 'reset';

}#

###########################################################
# Subroutine to get filename before file type suffix      #
# (i.e. prefix) (e.g. to get "test" from "test.fasta")    #
# New version using rindex rather than index.  This means #
# it can deal with files like filename.something.vcf      #
###########################################################

sub get_prefix
{
	my $_filename 	= "";
	my $_dot_pos	= 0;

	$_filename = $_[0];
	$_dot_pos = rindex($_filename,".");

	if (rindex($_filename,".") > 0)
	{
		$_filename = substr($_filename, 0, rindex($_filename,"."));
	}
	if (rindex($_filename,".") == -1)
	{
		$_filename = $_filename;
	}

	$_filename = $_filename;

} # get_prefix