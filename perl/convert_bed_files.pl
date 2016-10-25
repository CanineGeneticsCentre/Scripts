#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	convert_bed_files            						                #     
#									                                    #
#	converts bed 12 to bed                      						#
#########################################################################

##############################
# Mike Boursnell June 2014   #
# Animal Health Trust        #
# Newmarket                  #
# UK                         #
# mike.boursnell@aht.org.uk  #
##############################

use strict;
use Getopt::Long;
use File::Basename;
use Term::ANSIColor;
use Cwd;

# VERSION OF SOFTWARE #
my $version						= "2";

# FILE NAMES
my $bed12_file					= ""; # File in BED12 format
my $bed_file					= ""; # File in BED 4-column format
my $command_log					= ""; # File with log of commands used

# VARIABLES
my $prefix						= "";
my $single_line					= "";
my $chromosome					= "";
my $region_start				= "";
my $region_end					= "";
my $label						= "";
my $strand						= "";					
my $block_count					= "";
my $block_sizes					= "";
my $block_starts				= "";
my $answer						= "";
my $exon_number					= "";
my $exon_name					= "";

# NUMBERS
my $no_of_columns				= 0;
my $no_of_blocks				= 0;
my $no_of_sizes					= 0;
my $no_of_starts				= 0;
my $exon_start					= 0;
my $exon_end					= 0;

# ARRAYS
my @item						= ();
my @sizes_array					= ();
my @starts_array				= ();


print color 'bold magenta';

print "\n\n";
print "################################\n";
print color 'bold white';
print "      convert_bed_files     \n";
print color 'bold magenta';
print "################################\n\n\n";

print color 'yellow';
print "Version $version\n\n";

print "  - This program converts a BED12 file to a 4-column BED file.\n\n";

print color 'reset';

##############################################################
# Get the name of the BED file - containing the exons        #
##############################################################
$bed12_file = &get_file("What is the name of your BED12 file?",".bed");

##############################################################
# Make up default name for bedtools output file              #
##############################################################
$prefix = &get_prefix($bed12_file);
$bed_file = $prefix."_converted.bed";

##############################################################
# Open command log file and write first details              #
##############################################################
&open_command_log;

open (BED12, "$bed12_file") || die "Cannot open $bed12_file";
open (BED, ">$bed_file")|| die "Cannot create output BED file: $bed_file";

while ($single_line = <BED12> ) 
{
	chomp $single_line;
	#Split line at white space into the array @item
	@item=split(/\s/,$single_line);

	$no_of_columns = (scalar @item) + 1;

	if ($no_of_columns >= 7)
	{
		$chromosome   = $item[0];
		$region_start = $item[1];
		$region_end   = $item[2];
		$label        = $item[3];
		$strand   	  = $item[5];

		$no_of_blocks = $item[9];
		$block_sizes  = $item[10];
		$block_starts = $item[11];
	}

	#print "$single_line\n\n";
	#print "$chromosome\t$region_start\t$region_end\t$label\t$strand\n";
	#print "$no_of_blocks\t$block_sizes\t$block_starts\n";

	#Now deal with exon starts and sizes
	@sizes_array = split(/,/,$block_sizes);
	@starts_array = split(/,/,$block_starts);

	$no_of_sizes = (scalar @sizes_array) + 1;
	$no_of_starts = (scalar @starts_array) + 1;


	for ( $block_count = 1; $block_count <= $no_of_blocks; $block_count++ )
	{
		#print "$block_count:\t$sizes_array[$block_count-1]\t$starts_array[$block_count-1]\n";

		if ( $strand eq "+" )
		{
			$exon_start = $region_start + $starts_array[$block_count-1]; # Removed + 1 from here
			$exon_end = $exon_start + $sizes_array[$block_count - 1]; # Removed - 1 from here

			if ($block_count < 10){$exon_number = "0".$block_count}else{$exon_number = $block_count}

			$exon_name = $label."_exon_".$exon_number;
		}

		if ( $strand eq "-" )
		{
			$exon_start = $region_start + $starts_array[$block_count-1]; # Removed + 1 from here
			$exon_end = $exon_start + $sizes_array[$block_count-1]; # Removed - 1 from here

			$exon_number = $no_of_blocks - $block_count + 1;

			if ($exon_number < 10){$exon_number = "0".$exon_number}

			$exon_name = $label."_exon_".$exon_number."_R";
		}

		print BED "$chromosome\t$exon_start\t$exon_end\t$exon_name\t0\t$strand\n";

		#print BED "$chromosome\t$exon_start\t$exon_end\t$exon_name\t0\t$strand\t$region_start\t$region_end\t$starts_array[$block_count-1]\t$sizes_array[$block_count-1]\n";
	}
	#$answer=<STDIN>;
}

&print_message("Finished converting BED12 file","message");

print "Input BED12 file:    \t$bed12_file\n\n";
print "Output BED file:     \t$bed_file\n\n\n\n";

exit;

###############################################################
# Subroutine to print screen message                          #
# Argument 1: Text of the message                             #
# Argument 2: Type of message "message", "input" or "warning" #
###############################################################

sub print_message
{
	my $message_length 	= "";
	my $pos_count		= 0;
	my $char			= "";
	
	my $message = $_[0];
	my $style = $_[1];
	
	if ($style eq ""){$style = ""}

	$message_length = length($message);
	
	if ($style eq ""){$char = "#"}
	if ($style eq "input"){$char = "~"}
	if ($style eq "message"){$char = "#"}
	if ($style eq "warning"){$char = "!"}
	
	print "\n\n";
	print color ' bold yellow';
	if ($style eq "warning"){print color ' bold red'}
	if ($style eq "input"){print color ' bold white'}
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n$char    $message    $char\n";
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n\n";
	print color 'reset';

}

#############################################
# Subroutine to execute unix command        #
# Argument 1: command to execute            #
#############################################
sub run_unix_command
{
	my $unix_command = "";
	$unix_command = $_[0];
		
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	print COMMAND_LOG "$unix_command\n";

	print("$unix_command\n\n");
	system("$unix_command");
}


####################################################################
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test" from "test.fasta")                           #
####################################################################
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


#################################################################
# Subroutine to get a file name                                 #
# Argument 1: Text that user sees, asking for the file          #
# Argument 2: suffix of files to search for with ls e.g. ".bed" #
#################################################################
sub get_file
{
	my $input_message 	= $_[0];
	my $file_type	 	= $_[1];
	my $_file_name		= "";

	if ($file_type eq ""){$file_type = "*"}

	until ((-e $_file_name) || (lc $_file_name eq "q"))
	{
		&print_message("$input_message","input");
		print "> ";

		$_file_name = <STDIN>;
		chomp $_file_name;

		# TEMPORARY DEFAULTS FOR TESTING
		if (($file_type eq ".bam") &&($_file_name eq "")){$_file_name = "AS_HC_AS_case_01_small.bam"}  # TEMP!!!
		if (($file_type eq ".bed") &&($_file_name eq "")){$_file_name = "AS_HC_small_2.bed"}  # TEMP!!!

		if ($_file_name eq "ls") {print "\n";system ("ls *"."$file_type")}
		if (($_file_name ne "ls")  && (lc $_file_name ne "q") && ($_file_name ne ""))
		{

			if (!-e $_file_name){print "/n  >>>>>>>>>>>>>  ERROR.  File $_file_name cannot be found <<<<<<<<<<<\n\n";}

			if (-e $_file_name)
			{
				$_file_name = $_file_name;
			}
		} # not ls
	} # until loop

	print "\n";
	if ($file_type eq ".txt"){system("dos2unix $_file_name")}

	$_file_name = $_file_name;
}

####################################################
# Subroutine to print to screen and to COMMAND_LOG #
####################################################
sub print_both
{
	my $message = $_[0];

	print "$message";
	print COMMAND_LOG "$message";
}

########################################################
# Subroutine to open COMMAND_LOG for bedtools_coverage #
########################################################

sub open_command_log
{
	$command_log = "convert_bedfiles_command_log.out";

	open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

	print COMMAND_LOG "COMMAND LOG for convert_bedfiles version $version\n\n";
	print COMMAND_LOG "BED12 file:                \t$bed12_file\n";
	print COMMAND_LOG "BED file:                  \t$bed_file\n";
}


##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

foreach (@_) {s/\n//g}  
foreach (@_) {s/\r//g}  

}
