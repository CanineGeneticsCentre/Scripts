#!/usr/bin/perl -w

#########################################################
#														#      
#	NGS CAPTURE SUCCESS ANALYSIS v 4.3					#     
#														#
#	THIS PERL SCRIPT WILL ANALYSE ILLUMINA NGS DATA		#
#														#
#########################################################

# The purpose of this script is to allow the quick analysis of whole genome versus targetted regions
# in order to calculate target capture efficiency.
# This version of the script will process multiple samples at the same time
# The output is a results folder containing seperate excel files for every bam file...
# ...as well as a single excel file into which all of these have been combined
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# There is the option of running the analyses
# with AND without duplicates.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#############################
# Louise Downs May 2012     #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# louise.downs@aht.org.uk   #
#############################

use strict;
use Getopt::Std ;
use File::Basename ;

####################
# Define variables #
####################
#Individual depth analyses

my $bait = "";
my $bam_file	= "";
my $bam_file_count ="";
my $command	= "";
my $count;
my $data = "";
my $data_type = "";
my $dup = "";
my $dup_ans = "";
my $filename_out_dup = "";
my $filename_out_no_dup = "";
my $i;
my $name	= "";
my $name_ans	= "";
my $no_bam_files	= "";
my $proceed	= "";
my $region_chr = "";
my $region_start = "";
my $region_end = "";
my $target_ans	= "";
my $target_interval_file = "";
my $to		= "";

my @bam_files = qw();		#qw = list of quoted words - don't use quotation marks
my @lines = ();

#Combining individual depth files into a single spreadsheet

my $bait_coverage = "";
my $bait_set = "";
my $bait_territory= "";
my $bases_align_wg = "";
my $bases_on_target = "";
my $filename_in = "";
my $h = "";
my $pct_bases_on_bait = "";
my $pct_bases_on_target = "";
my $pct_reads_align_wg = "";
my $reads_align_wg = "";
my $target_coverage = "";
my $target_territory = "";
my $ten_X = "";
my $two_X = "";
my $thirty_X = "";
my $total_reads = "";
my $twenty_X = "";

my @lines_in = qw();
my @target_line = qw();

########################
# Define non variables #
########################

my $title='Perl Mail demo';
my $from= 'louise.downs@aht.org.uk';
my $subject='NGS CAPTURE EFFICIENCY ANALYSIS';

#START TIMERS

BEGIN { our $start_run = time(); }

BEGIN { our $start_run2 = time(); }

BEGIN { our $start_run3 = time(); }

#TURN LOGGER ON

$| = 1;

open(STDOUT, "| tee log.rtf");

    use Term::ANSIColor;
    print color 'magenta';

print "\n";	
print "                 ##############################################\n";
print color 'bold white';
print "                 NGS data analysis: Capture efficiency pipeline\n";
print color 'reset';
print color 'magenta';
print "                 ##############################################\n";
print "\n";

print color 'bold white';

print "**************************************************************************************************\n";
print "* The purpose of this analysis is to compare the data statistics, including read depth and target*\n";
print "*               capture efficiency for a number of already-aligned bam files.                    *\n";
print "*                                                                                                *\n";
print "* ";
print color 'bold magenta';
print "Points to consider before starting:";
print color 'bold white';
print "                                                            *\n";
print "* 1. Input files required are Bait Interval file (.bed - see below for correct layout) and       *\n";
print "*    aligned binary read files (.bam).                                                           *\n";
print "* 2. Multiple bam files can be analysed in the same run, however:                                *\n";
print "*    a.  Bam files should all have been processed in the same way e.g. raw_align.bam files from *\n";
print "*        the NGS analysis pipeline.                                                              *\n";
print "*    b.  Ensure the bam files were originally created by aligning to the whole genome and not to *\n";
print "*        a region of interest. If in doubt use raw_align.bam files created during the NGS       *\n";
print "*        Analysis pipeline.                                                                      *\n";
print "*                                                                                                *\n";
print "* The bait interval file is the RNA bait file created during bait design, but must be edited to  *\n";
print "* look as follows:                                                                               *\n";
print "*                                                                                                *\n";
print '* @SQ SN:canfam2 LN:2531673953				          [copy this line verbatim]      *';
print "\n* chr#	bait_start	bait_end	+	bait_name                                        *\n";
print "* chr13	39648703	39648823	+	BI424281914_9     [example]                      *\n";
print "**************************************************************************************************\n\n";

print color 'reset';

&error_check;

#############################
# Name of Analysis and email   
#############################

print "\n~~~~~~~~~~~~~~~";
print "\nYour details...";
print "\n~~~~~~~~~~~~~~~\n";

print "\nPlease enter a name for this analysis (with no spaces):    ";
$name = <STDIN>;
chomp $name;

#################################################################
# Check if name of analysis has been used in the current folder 
#################################################################

if (-e "capture_results_$name/")
{ 
print color 'bold red';
print "\nA results folder with the name results_$name already exists!\n\n";
print color 'reset';
print "Would you like to specify a new name or use the existing folder?\n\n";
	
print "Enter 1 to specify a new name\n";
print "Enter 2 to use the existing folder\t\t"; 
	print color 'bold red';
	print "EXISTING FILES MAY BE OVERWRITTEN!\n";
	print color 'reset';
		
$name_ans = <STDIN>;
chomp $name_ans;

################################
# To specify a new name
#################################
	
if ($name_ans == 1)
	{
	print "\nPlease enter a new name for this analysis (with no spaces):    \n";
	$name = <STDIN>;
	chomp $name;

		if (-e "capture_results_$name/")
		{
		print color 'bold red';
		print "\nA results folder with the name results_$name also already exists!\n\n";
		print color 'reset';
		print "Would you like to start over or use the existing folder?\n\n";

		print "Enter 1 to start over\n";
		print "Enter 2 to use the existing folder\t\t"; 
		print color 'bold red';
		print "EXISTING FILES MAY BE OVERWRITTEN!\n";
		print color 'reset';
			
		$name_ans = <STDIN>;
		chomp $name_ans;
		
		###############################
		# To quit and start again
		##############################
		
		if ($name_ans == 1)
			{
			exit;
			}
		}
	}	
}

################################
# Create the results folder
##################################

unless (-e "capture_results_$name")
{
print "\nCreating new directory for results....\n";
$command = "mkdir capture_results_$name";
system("$command");
}

print("\nDirectory = 'capture_results_$name'\n");

###################################################
# Ask for email address
###################################################

print "\nPlease enter your email address:    ";
$to = <STDIN>;
chomp $to;

###############################################################
# Ask for the name of the Bait Interval File
################################################################

print "\nPlease enter the name of your Bait Interval file (with no spaces):    ";
$bait = <STDIN>;
chomp $bait;

#####################################
# Add bed if user hasn't added it #
#####################################

if (index($bait,".bed") == -1 )
{
	$bait = $bait.".bed";
}

################################################
# Ask how many bam files you want to run
###########################################

print "\nPlease input the number of Bam files you wish to analyse:   \n";
$no_bam_files = <STDIN>;
chomp $no_bam_files;

######################################
# Ask for the name(s) of bam input files
######################################

print "\nPlease input the name(s) of the bam files you wish to use, one at a time e.g. best_align_1234.bam:\n";

for $count ( 1 .. $no_bam_files )
{
print "\nBam file number $count :	";
$bam_file_count = <STDIN>;
chomp $bam_file_count;

#####################################
# Add bam if user hasn't added it #
#####################################
if (index($bam_file_count,".bam") == -1 )
	{
		$bam_file_count = $bam_file_count.".bam";
	}
####################################################
# Create an array of the bam files
###################################################

push(@bam_files, "$bam_file_count");
}

###########################################
# Check the bam files entered are correct
##########################################

print "\n\nThe Bam files entered are:\n";
foreach (@bam_files) {
	print "\t$_\n";
	}

print "\nAre these correct?";
print "\n\nPlease press enter to proceed (or 'N' to re-enter):	";
$proceed = <STDIN>;
chomp $proceed;

############################################
# If bam files are incorrect, re-enter
######################################

if (lc $proceed eq "n")
{
@bam_files = ();				#Empty array

print "\nPlease input the name of the bam files you wish to use, one at a time e.g. raw_align_1367.bam:\n";

for $count ( 1 .. $no_bam_files )
	{
	print "\nBam file number $count :	";
	$bam_file_count = <STDIN>;
	chomp $bam_file_count;

	#####################################
	# Add bam if user hasn't added it #
	#####################################
	if (index($bam_file_count,".bam") == -1 )
		{
			$bam_file_count = $bam_file_count.".bam";
		}
	####################################################
	# Create an array of the bam files
	###################################################

	push(@bam_files, "$bam_file_count");
	}

##################################################
# Check the bam files entered are correct
###############################################

print "\nThe Bam files entered are:\n";
foreach (@bam_files) {
	print "\t$_\n";
	}

print "\nAre these correct?";
print "\n\nPlease press enter to proceed (or 'Q' to quit):	";
$proceed = <STDIN>;
chomp $proceed;

############################################
# If bam files are incorrect again, quit and start again
########################################################

if (lc $proceed eq "q")
	{
	exit;
	}
}

#########################################################
# Ask if duplicates should be removed
########################################################

	print "\nWould you like to remove duplicates from the bam files?\n";
	print "\t(To compare metrices between files with and without duplicates select BOTH)\n\n";
	
	print "Enter 1 for Yes\n";
	print "Enter 2 for No\t\n"; 
	print "Enter 3 for Both\t\t";
		
	$dup_ans = <STDIN>;
	chomp $dup_ans;
		
if ($dup_ans == 1){$dup = "yes"}
if ($dup_ans == 2){$dup = "no"}
if ($dup_ans == 3){$dup = "both"}

if ($dup eq "yes" or "both")
{
####################################
#Ask if the data is SE or PE data?
####################################

print "\nDo you have a Paired-end or Single-end dataset?\n\n";

print "   Enter 1 for SE\n";
print "   Enter 2 for PE\t\t";

$data_type = <STDIN>;
chomp $data_type;

if ($data_type == 1){$data = "SE"}
if ($data_type == 2){$data = "PE"}

if ($data eq "SE")
	{
		print "\n\nYou have selected Single-end\n"
	}

if ($data eq "PE")
	{
		print "\n\nYou have selected Paired-end\n"
	}
}

###################################################
# Create a new target interval file
###################################################
###########################################################################
# How many regions would you like to target?                #
###########################################################################

print "\nPlease input the number of different regions you wish to target (e.g. 2):	\t";
$target_ans = <STDIN>;
chomp $target_ans;

###################################################
# Create target interval file
###################################################

open ( OUTFILE, ">target_interval_$name.txt" )
	or die "Couldn't create combined output file.";
$target_interval_file="target_interval_$name.txt";
print "\nThe Target Interval file is: \t$target_interval_file\n";
	print OUTFILE '@SQ SN:canfam2 LN:2531673953';
	print "\nHeading line added to $target_interval_file.\n";
close( OUTFILE );

#####################################################
# Define the chromosome, start and end for each region
#####################################################

for $count ( 1 .. $target_ans )
{
print "\nRegion $count chr:	";										# Define the chromosome
$region_chr = <STDIN>;
chomp $region_chr;

print "Region $count start location:	";							# Define the start location
$region_start = <STDIN>;
chomp $region_start;

print "Region $count end location:	";								# Define the end location
$region_end = <STDIN>;
chomp $region_end;

####################################################
# Add the region information to the target interval file
###################################################

open( OUTFILE, ">>$target_interval_file" );
	print OUTFILE "\nchr$region_chr\t$region_start\t$region_end\t+\tRegion_$count";
	print "\n\tRegion $count written to $target_interval_file.\n";
close( OUTFILE );

}

###################################################
# Create new combined output (excel) file and add headings
###################################################

if ($dup eq "no" or "both") {
open ( OUTFILE, ">capture_results_$name/Coverage_Summary_dup_$name.xls" )
	or die "Couldn't create combined output file.";
$filename_out_dup="capture_results_$name/Coverage_Summary_dup_$name.xls";
print "\nThe new depth summary file is: \t$filename_out_dup\n";
	print OUTFILE "Bam files\tTotal no of reads\tNo of reads aligned to whole genome\t";
	print OUTFILE "Percentage of reads aligned to whole genome\tBases mapped to whole genome\tBases mapped to target region\t";
	print OUTFILE "Mean coverage of baited region\tMean coverage of targetted region\tPercentage of bases mapped to baits\t";
	print OUTFILE "Pergentage of bases mapped to targets\tPercentage of bases with >2X coverage\t";
	print OUTFILE "Percentage of bases with >10X coverage\tPercentage of bases with >20X coverage\tPercentage of bases with >30X coverage\t\n";
	print "\n\tHeadings added $filename_out_dup.\n\n";
close( OUTFILE );
}

if ($dup eq "yes" or "both") {
open ( OUTFILE, ">capture_results_$name/Coverage_Summary_no_dup_$name.xls" )
	or die "Couldn't create combined output file.";
$filename_out_no_dup="capture_results_$name/Coverage_Summary_no_dup_$name.xls";
print "\nThe new depth summary file is: \t$filename_out_no_dup\n";
	print OUTFILE "Bam files\tTotal no of reads\tNo of reads aligned to whole genome\t";
	print OUTFILE "Percentage of reads aligned to whole genome\tBases mapped to whole genome\tBases mapped to target region\t";
	print OUTFILE "Mean coverage of baited region\tMean coverage of targetted region\tPercentage of bases mapped to baits\t";
	print OUTFILE "Pergentage of bases mapped to targets\tPercentage of bases with >2X coverage\t";
	print OUTFILE "Percentage of bases with >10X coverage\tPercentage of bases with >20X coverage\tPercentage of bases with >30X coverage\t\n";
	print "\n\tHeadings added $filename_out_no_dup.\n\n";
close( OUTFILE );
}

$command =  "clear";
system("$command");

##############################################################
# Print the details entered to the screen to check
###########################################################

print color "bold white";

print "\n================================================\n";
print "SUMMARY - PLEASE CHECK THESE DETAILS CAREFULLY!!\n";
print "================================================\n";

print color "bold green";

print "\nYOUR DETAILS :\n";

print color "bold white";

print "\nName of this analysis:\t$name\n\n"; 

print "Your email address:\t$to\n";

print color "bold green";

print "\nDATA FILES:\n\n";

print color "bold white";

print "Bam file(s):\n";

foreach (@bam_files) {
	print "\t\t\t$_\n";
	}

print "\nBait Interval File:\t";
print "$bait\n";
	
print "\nRegion(s) of interest:\n";

open ( INFILE, "<$target_interval_file" );
@lines = <INFILE>;
close ( INFILE );

for $i ( 1..$#lines ) {
	print ("$lines[ $i ]");
	}

print color "bold green";

print "\n\nSETTINGS :\n";

print color "bold white";

print "\nDuplicates:\t\t";

if ($dup eq "yes" )
{
print "Duplicates will be removed\n";

print "\nAnalysis:\t\t";

if  ($data eq "SE")
	{
	print "Single-end analysis\n";
	}

if  ($data eq "PE")
	{
	print "Paired-end analysis\n";
	}
}

if ($dup eq "both" )
{
print "Duplicates will be removed and metrices calculated for bam files with AND without duplicates\n";

print "\nAnalysis:\t\t";

if  ($data eq "SE")
	{
	print "Single-end analysis\n";
	}

if  ($data eq "PE")
	{
	print "Paired-end analysis\n";
	}
}

if ($dup eq "no")
{
print "Duplicates will no be removed\n";
}

print color "reset";

print "\n\nPlease press enter to proceed (or 'Q' to quit):      \n";

print color "bold red";
print "\nNOTE - RUN TIME MAY BE SEVERAL HOURS.\n"; 

$proceed = <STDIN>;
chomp $proceed; 
 
if (lc $proceed eq "q")
{
 exit;
}

print color "reset";

#############################################################
# Start analysis of bam files looping through array @bam_files
#################################################################

print "*****************************\n";
print "*Processing the Bam file(s).*\n";
print "*****************************\n";

foreach my $bam_files (@bam_files)
{
$bam_file = $bam_files;

print "\nThe bam file currently being processed is $bam_file\n\n";

$command = "cp $bam_file aligned_sorted.bam";
system("$command");

###############################################################
# Removal of PCR duplicates
###############################################################

if ($dup eq "yes" or "both")
	{
	open(STDERR, "| tee Duplicate_info.rtf");

	print "\nRemoving PCR duplicates.\n";

	if  ($data eq "PE")
		{
		$command =  "/opt/samtools/samtools rmdup aligned_sorted.bam duplicates_removed.bam";

		print("$command\n");
		system("$command");
		}
		
	if  ($data eq "SE")
		{
		$command =  "/opt/samtools/samtools rmdup -s aligned_sorted.bam duplicates_removed.bam";

		print("$command\n");
		system("$command");
		}

	print 	"\n Duplicates have been removed from $bam_file\n";

	print	"\nre-naming original file...\n\n";

	$command =  "mv aligned_sorted.bam aligned_plus_dup.bam";

	print("$command\n");
	system("$command");

	print	"\nre-naming new file...\n\n";

	$command =  "mv duplicates_removed.bam aligned_sorted.bam ";

	print("$command\n\n");
	system("$command");

	close(STDERR);
	}

######################################################################
# Run the metrics for the bam file
##########################################################

$command = "/usr/bin/java16 -jar -Xmx2g /opt/picard/CalculateHsMetrics.jar I=aligned_sorted.bam O=capture_eff BI=$bait  TI=$target_interval_file VALIDATION_STRINGENCY=LENIENT";
print("$command\n");
system("$command");

if ($dup eq "both")
	{
	$command = "/usr/bin/java16 -jar -Xmx2g /opt/picard/CalculateHsMetrics.jar I=aligned_plus_dup.bam O=capture_eff_dup BI=$bait  TI=$target_interval_file VALIDATION_STRINGENCY=LENIENT";
	print("$command\n");
	system("$command");
	}

######################################
# Move into results folder          #
######################################

print "\nThe file currently being processed is $bam_file...\n";
print "\nSorting results into new directory....\n";

	if ($dup eq "no") 
	{
	$command = "mv capture_eff capture_results_$name/Capture_efficiency_dup_$bam_file.txt";
	system("$command");
	}

	if ($dup eq "yes" or "both") 
	{
	$command = "mv capture_eff capture_results_$name/Capture_efficiency_no_dup_$bam_file.txt";
	system("$command");
	}

	if ($dup eq "both") 
	{
	$command = "mv capture_eff_dup capture_results_$name/Capture_efficiency_dup_$bam_file.txt";
	system("$command");
	}

##############################################################
# Add the relevant results to the Summary excel workbook
###########################################################

if ($dup eq "no" or "both")
	{
	print "\nAdding relevant data to $filename_out_dup....\n\n";

	###################################################
	# Read from individual capture efficiency files
	###############################################

	open( INFILE, "<capture_results_$name/Capture_efficiency_dup_$bam_file.txt" )
		or die( "Couldn't open file $filename_in: $!\n" );
	@lines_in = <INFILE>;
	close( INFILE );

	####################################################################
	# Locate and define the target value i.e. the total number of bases
	########################################################

		for $i ( 0 .. $#lines_in ) {
			chomp ( $lines_in[ $i ] );
			if ( $lines_in[ $i ] =~ m/BAIT_SET\t/ ) {					 	# Target the heading line, starting with BAIT_SET
			@target_line = split (/\t/, $lines_in[ $i + 1 ]);		# Split the next line into a new array
			for $h ( 0 .. $#target_line ) {
				chomp ( $target_line[ $h ] );
				$total_reads=$target_line[5];
				$reads_align_wg=$target_line[10];
				$pct_reads_align_wg=$target_line[11];
				$bases_align_wg=$target_line[12];
				$bases_on_target=$target_line[16];
				$bait_coverage=$target_line[20];
				$target_coverage=$target_line[21];
				$pct_bases_on_bait=$target_line[22];
				$pct_bases_on_target=$target_line[23];
				$two_X=$target_line[27];
				$ten_X=$target_line[28];
				$twenty_X=$target_line[29];
				$thirty_X=$target_line[30];
				}
			}
		}

	##############################################
	# Print details to combined excel depth file
	##########################################

	open( OUTFILE, ">>capture_results_$name/Coverage_Summary_dup_$name.xls");
		print OUTFILE "$bam_file\t$total_reads\t$reads_align_wg\t$pct_reads_align_wg\t";
		print OUTFILE "$bases_align_wg\t$bases_on_target\t$bait_coverage\t$target_coverage\t$pct_bases_on_bait\t";
		print OUTFILE "$pct_bases_on_target\t$two_X\t$ten_X\t$twenty_X\t$thirty_X\t\n";
		print "\tRelevant information added to $filename_out_dup.\n";
	close( OUTFILE );
	}

if ($dup eq "yes" or "both")
	{
	print "\nAdding relevant data to $filename_out_no_dup....\n\n";

	###################################################
	# Read from individual capture efficiency files
	###############################################

	open( INFILE, "<capture_results_$name/Capture_efficiency_no_dup_$bam_file.txt" )
		or die( "Couldn't open file $filename_in: $!\n" );
	@lines_in = <INFILE>;
	close( INFILE );

	####################################################################
	# Locate and define the target value i.e. the total number of bases
	########################################################

		for $i ( 0 .. $#lines_in ) {
			chomp ( $lines_in[ $i ] );
			if ( $lines_in[ $i ] =~ m/BAIT_SET\t/ ) {					 	# Target the heading line, starting with BAIT_SET
			@target_line = split (/\t/, $lines_in[ $i + 1 ]);		# Split the next line into a new array
			for $h ( 0 .. $#target_line ) {
				chomp ( $target_line[ $h ] );
				$total_reads=$target_line[5];
				$reads_align_wg=$target_line[10];
				$pct_reads_align_wg=$target_line[11];
				$bases_align_wg=$target_line[12];
				$bases_on_target=$target_line[16];
				$bait_coverage=$target_line[20];
				$target_coverage=$target_line[21];
				$pct_bases_on_bait=$target_line[22];
				$pct_bases_on_target=$target_line[23];
				$two_X=$target_line[27];
				$ten_X=$target_line[28];
				$twenty_X=$target_line[29];
				$thirty_X=$target_line[30];
				}
			}
		}

	##############################################
	# Print details to combined excel depth file
	##########################################

	open( OUTFILE, ">>capture_results_$name/Coverage_Summary_no_dup_$name.xls");
		print OUTFILE "$bam_file\t$total_reads\t$reads_align_wg\t$pct_reads_align_wg\t";
		print OUTFILE "$bases_align_wg\t$bases_on_target\t$bait_coverage\t$target_coverage\t$pct_bases_on_bait\t";
		print OUTFILE "$pct_bases_on_target\t$two_X\t$ten_X\t$twenty_X\t$thirty_X\t\n";
		print "\tRelevant information added to $filename_out_no_dup.\n";
	close( OUTFILE );
	}
}

##############################################
# Print generic details e.g. size of target and bait regions, to combined excel depth files
##########################################

if ($dup eq "no" or "both")
{
###################################################
# Read from individual capture efficiency files (with duplicates)
###############################################

open( INFILE, "<capture_results_$name/Capture_efficiency_dup_$bam_file.txt" )
	or die( "Couldn't open file $filename_in: $!\n" );
@lines_in = <INFILE>;
close( INFILE );

####################################################################
# Locate and define the target value i.e. the total number of bases
########################################################
	for $i ( 0 .. $#lines_in ) {
		chomp ( $lines_in[ $i ] );
		if ( $lines_in[ $i ] =~ m/BAIT_SET\t/ ) {					 	# Target the heading line, starting with BAIT_SET
		@target_line = split (/\t/, $lines_in[ $i + 1 ]);		# Split the next line into a new array
		for $h ( 0 .. $#target_line ) {
			chomp ( $target_line[ $h ] );
			$bait_set=$target_line[0];
			$bait_territory=$target_line[2];
			$target_territory=$target_line[3];
			}
		}
	}

open( OUTFILE, ">>capture_results_$name/Coverage_Summary_dup_$name.xls");
	print OUTFILE "\nBait file\t$bait_set\nSize of baited region\t$bait_territory\n";
	print OUTFILE "Size of targetted region\t$target_territory\n";
	print "\tRelevant information added to $filename_out_dup.\n";
close( OUTFILE );
}

if ($dup eq "yes" or "both")
{
###################################################
# Read from individual capture efficiency files
###############################################

open( INFILE, "<capture_results_$name/Capture_efficiency_no_dup_$bam_file.txt" )
	or die( "Couldn't open file $filename_in: $!\n" );
@lines_in = <INFILE>;
close( INFILE );

####################################################################
# Locate and define the target value i.e. the total number of bases
########################################################
	for $i ( 0 .. $#lines_in ) {
		chomp ( $lines_in[ $i ] );
		if ( $lines_in[ $i ] =~ m/BAIT_SET\t/ ) {					 	# Target the heading line, starting with BAIT_SET
		@target_line = split (/\t/, $lines_in[ $i + 1 ]);		# Split the next line into a new array
		for $h ( 0 .. $#target_line ) {
			chomp ( $target_line[ $h ] );
			$bait_set=$target_line[0];
			$bait_territory=$target_line[2];
			$target_territory=$target_line[3];
			}
		}
	}

open( OUTFILE, ">>capture_results_$name/Coverage_Summary_no_dup_$name.xls");
	print OUTFILE "\nBait file\t$bait_set\nSize of baited region\t$bait_territory\n";
	print OUTFILE "Size of targetted region\t$target_territory\n";
	print "\tRelevant information added to $filename_out_no_dup.\n";
close( OUTFILE );
}

#RUN TIMER END

my $end_run = time();
my $run_time = $end_run - our $start_run;

open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: $subject COMPLETE ($name)\n\n";
## Mail Body
print MAIL "Your next target capture efficiency analysis ($name) is complete\n\n";
print MAIL "For pipeline details see log.rtf\n\n";
print MAIL "Run time : $run_time seconds\n";
printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
close(MAIL);

print "\nANALYSIS COMPLETE! YOU HAVE BEEN NOTIFIED BY EMAIL\n\n";

print "Run time : $run_time seconds\n";
printf "%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];

print "XXXXXXXXXXXXXXXXXXXXXXX\n";
print "  ERRORS AND WARNINGS  \n";
print "XXXXXXXXXXXXXXXXXXXXXXX\n\n";

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print if /\bwarning\b/i;
}
close(LOG);

print "\n";

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print if /\berror\b/i;
}
close(LOG);

print "\nPlease check log.rtf for full details\n\n";

#############################
# Turn off logging          #
#############################

close(STDOUT);

#################################################
# Move log file to the results folder
################################################

$command = "mv log.rtf capture_results_$name/Log_$name.rtf";
system("$command");

$command = "rm $target_interval_file";
system("$command");

$command = "rm aligned_sorted.bam";
system("$command");

if ($dup eq "yes" or "both") {
$command = "rm aligned_plus_dup.bam";
system("$command");

$command = "rm Duplicate_info.rtf";
system("$command");
}

#########################################################
#
# Subroutine to stop at a point
#
######################################

sub error_check
{
my $proceed = "";

print "\n\nPlease press enter to proceed (or 'Q' to quit):	";
$proceed = <STDIN>;
chomp $proceed;

if (lc $proceed eq "q")
	{
	exit;
	}
}