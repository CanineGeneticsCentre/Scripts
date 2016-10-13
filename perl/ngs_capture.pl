#!/usr/bin/perl -w

#########################################################
#														#      
#	NGS CAPTURE SUCCESS ANALYSIS     					#     
#														#
#	THIS PERL SCRIPT WILL ANALYSE ILLUMINA NGS DATA		#
#														#
#########################################################

# The purpose of this script is to allow the quick analysis of whole genome versus targetted regions
# in order to calculate target capture efficiency.
# This version of the script will process multiple samples at the same time
# The output is a results folder containing separate output files for every bam file...
# ...as well as a single output file into which all of these have been combined
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# There is the option of running the analyses
# with OR without duplicates.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#######################################
# Louise Downs May 2012               #
# Animal Health Trust                 #
# Newmarket                           #
# UK                                  #
# louise.downs@aht.org.uk             #
# modified by Mike Boursnell May 2014 #
#######################################

use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;
use Cwd;

my $version					= "9";

####################
# Define variables #
####################

#File names
my $command_log				= "";
my $screen_log				= "";
my $list_file				= ""; # File containing the list of BAM files (usually .txt)
my $baits_file 				= ""; # Original baits file, not formatted correctly
my $baits_file_new			= ""; # Baits file formatted for picard
my $bam_file				= "";
my $bam_file_in_list_file	= "";
my $bam_file_without_path	= "";
my $aligned_sorted_bam		= "";
my $duplicates_removed_bam	= "";
my $capture_efficiency_file	= "";
my $coverage_summary_file	= "";
my $header_temp_file		= "";

#Directories
my $path_to_bam_file		= ""; # Can be set to /ahtlocal/local_storage/bam_files

# Numbers
my $list_count				= 0;
my $item_array_size			= 0;
my $chromosome_number		= 0;
my $target_region_count		= 0;
my $no_of_target_regions	= 0;
my $region_choice			= 0;
my $row_count				= 0;
my $no_of_rows				= 0; # no of rows of data to add to summary output file
my $no_of_exons				= 0; # if the baits file is for exome sequencing
my $exon_count				= 0; 
my $no_of_header_lines		= 0;
my $slash_pos				= 0;

#Boolean
my $show					= "false";
my $use_default_regions		= "yes";
my $some_files_not_found	= "false";
my $exome_sequencing		= "true"; # or false

#Strings
my $command					= "";
my $chromosome				= "";
my $position_string			= "";
my $prefix					= "";
my $answer					= "";
my $paired_ends_type 		= "";
my $paired_ends_type_type 	= "";
my $remove_duplicates 		= "";
my $remove_duplicates_ans 	= "";
my $filename_out_dup 		= "";
my $filename_out_no_dup 	= "";
my $region					= "";
my $read_file_method		= ""; # multiple or single: method of reading in the BAM files
my $run_title				= "";
my $run_title_ans			= "";
my $results_folder			= "";
my $proceed					= "";
my $region_chr 				= "";
my $region_start 			= "";
my $region_end 				= "";
my $target_ans				= "";
my $target_interval_file	= "";
my $email_address			= "";
my $bait_name				= "";
my $plus					= "";
my $explanation				= "";

my $no_of_files				= ""; # no of BAM files
my $count					= 0;
my $i						= 0;
my $capture_array_count 	= 0;

# Arrays
my @bam_file_array			= ();
my @bam_file_array_with_path= ();
my @bam_files 				= qw();		#qw = list of quoted words - don't use quotation marks
my @lines 					= ();
my @item					= ();  # single_line is split into array item
my @smallest_position		= ();
my @largest_position		= ();
my @target_region_array		= ();

# 2 dimensional array for storing all details of output for BAM files
my @file_number 					= ();
my @column_number 						= ();
my @summaryarray_2d 		= (\@file_number, \@column_number);
#Combining individual depth files into a single spreadsheet

my $baits_file_coverage 	= "";
my $baits_file_set 			= "";
my $baits_file_territory	= "";
my $bases_align_wg 			= "";
my $bases_on_target 		= "";
my $filename_in 			= "";
my $h 						= "";
my $pct_bases_on_bait 		= "";
my $pct_bases_on_target 	= "";
my $pct_reads_align_wg 		= "";
my $reads_align_wg 			= "";
my $target_coverage 		= "";
my $target_territory 		= "";
my $ten_X 					= "";
my $two_X 					= "";
my $thirty_X 				= "";
my $total_reads 			= "";
my $twenty_X 				= "";
my $single_line				= "";

my @target_file_array		= ();
my @header_file_array		= ();
my @baits_file_array		= ();
my @baits_file_new_array	= ();
my @capture_array 			= qw();
my @target_line 			= qw();

########################
# Define non variables #
########################


my $from = 'NGS_capture@aht.org.uk'; # Who e-mails come from

#START TIMERS

BEGIN { our $start_run = time(); }
BEGIN { our $start_run2 = time(); }
BEGIN { our $start_run3 = time(); }



print color 'magenta';

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "              ngs_capture   \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "    - This program compares read depth and target capture efficncy for BAM files\n\n";
print "    - Note: read these details below: \n\n";



print "********************************************************************************************************\n";
print "* The purpose of this analysis is to compare the data statistics, including read depth, target         *\n";
print "* capture efficiency and target region coverage for a number of already-aligned BAM files.             *\n";
print "*                                                                                                      *\n";
print "* The bait interval file is the RNA bait file created during bait design, but must be in this format:  *\n";
print "*                                                                                                      *\n";
print "*    chr20     015543005     015543125     BI426399621_0     1000     +                                *\n";
print "*                                                                                                      *\n";
print "*                                                                                                      *\n";
print "* If the baits file is missing, the target regions can be specified in the same format:                *\n";
print "*                                                                                                      *\n";
print "*    chr20     015543005     015543125     target_region_name     1000     +                           *\n";
print "*                                                                                                      *\n";
print "*    (Note: use a unique target region name, and a column with 1000 and a column with +)               *\n";
print "*                                                                                                      *\n";
print "*                                                                                                      *\n";
print "* Multiple BAM files can be analysed in the same run.  You must be sure that:                          *\n";
print "*                                                                                                      *\n";
print "*    a.  BAM files have all been processed in the same way e.g. final.bam files from                   *\n";
print "*        fastq2vcf (or intermediate f2b.bam files from fastq2vcf)                                      *\n";
print "*                                                                                                      *\n";
print "*    b.  Ensure the BAM files were originally created by aligning to the whole genome and not to       *\n";
print "*        a region of interest. If in doubt use f2b.bam files created during by fastq2vcf.              *\n";
print "*                                                                                                      *\n";
print "********************************************************************************************************\n\n";

print color 'reset';



#############################
# Name of Analysis and email   
#############################

&print_message("Please enter a name for this analysis (with no spaces)","input");

$run_title = <STDIN>;
chomp $run_title;

$results_folder = $run_title."_ngs_capture_results";
$command_log = "$run_title"."_ngs_capture_command_log.out";
$screen_log = "$run_title"."_ngs_capture_screen_log.out";
$header_temp_file = "$run_title"."_header_temp.txt";


##################
# TURN LOGGER ON #
##################

$| = 1;

open(STDOUT, "| tee $screen_log");

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

print COMMAND_LOG "COMMAND LOG for ngs_capture version $version\n\n";
print COMMAND_LOG "Analysis name:           \t$run_title\n";

if ($run_title ne "test")
{
	if (-e $results_folder){print "\n>>>> Folder $results_folder exists already. Choose a new name, or delete the folder.  <<<<\n\n";exit;}
}

&run_unix_command("mkdir -p $results_folder");  ## -p is temp!!!!!



##########################
# Set up some file names #
##########################
$coverage_summary_file = "$run_title"."_"."coverage_summary.out";
$target_interval_file=$run_title."_target_intervals.txt";


######################################
# E-mail address for notifications   #
######################################
&print_message("Please enter your email address if you want to be notified when the analysis finishes","input");
$email_address = <STDIN>;
chomp $email_address;

# Some short cuts
if ($email_address eq "m"){$email_address = 'mike.boursnell@aht.org.uk';}
if ($email_address eq "o"){$email_address = 'oliver.forman@aht.org.uk';}
if ($email_address eq "s"){$email_address = 'sally.ricketts@aht.org.uk';}
if ($email_address eq "b"){$email_address = 'rebekkah.hitti@aht.org.uk';}


#######################################################
# Ask where the BAM files are located                 #
#######################################################
&print_message("Where are the BAM files located?","input");

print "  <1> In the current directory (or path is specified in the input file)\n";
print "  <2> Somewhere else that you will specify a path to\n\n";

$answer = <STDIN>;chomp $answer;

if (substr($answer,0,1) eq "1"){$path_to_bam_file = ""}

if (substr($answer,0,1) eq "2")
{
	print "Type in the exact path to the BAM files, like this:\n\n";
	print "(e.g. /ahtlocal/local_storage/bam_files/bam_genome/ ) \n\n";

	$path_to_bam_file = <STDIN>;
	chomp $path_to_bam_file;
}



#######################################################
# Ask how you want to get the names of the BAM files  #
#######################################################
until ($read_file_method eq "multiple" || $read_file_method eq "single")
{
	&print_message("How do you want to read the input files?","input");

	print "  <1> Using a file of file names for your BAM files  [DEFAULT]\n";
	print "  <2> Entering the file names individually\n\n";

	$answer = <STDIN>;chomp $answer;
	if ($answer eq ""){$answer = "1"}

	if (substr($answer,0,1) eq "1"){$read_file_method = "multiple"}
	if (substr($answer,0,1) eq "2"){$read_file_method = "single"}
}


#######################################################
# Get the name of the file with the list of BAM files #
#######################################################
if ($read_file_method eq "multiple")
{
	&print_message("The input is a file with a list of the original BAM file names.","message");
	$list_file = &get_file("Name of the file with the list of BAM files",".txt");

	#############################################
	# Get the BAM file names from the list file #
	#############################################
	&get_bam_files_from_list_file;

	print "\nList of BAM files:\n\n";

	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		print "\tFile $list_count: $bam_file_array_with_path[$list_count]\n";
	}

} # read_file_method = "multiple"


#######################################################
# Get the name(s) of the BAM files individually       #
#######################################################
if ($read_file_method eq "single")
{
	$list_count = 0;

	until (lc $answer eq 'q')
	{
		$bam_file = &get_file("Name of the BAM file (type 'q' to finish entering BAM files)",".bam");

		if (lc $bam_file ne "q")
		{
			$list_count = $list_count + 1;
			$bam_file_array[$list_count] = $bam_file;
			$bam_file_array_with_path[$list_count] = $bam_file;
		}
		
		if (lc $bam_file eq "q")
		{
			print "Finished entering BAM files\n";
			last;
		}
	}

	$no_of_files = $list_count;

	print "\nList of BAM files (entered individually):\n\n";
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		print "\tFile $list_count: $bam_file_array[$list_count]\n";
	}

} # read_file_method = "multiple"


##################################
# Get the name of the Baits file #
##################################
$baits_file = &get_file("Please enter the name of your Bait Interval file",".bed");
$prefix = &get_prefix ($baits_file);
$baits_file_new = $run_title."_".$prefix."_reformatted.bed";


#############################################
# Ask if this is exome sequencing           #
#############################################
&print_message("What sort of capture is this?","input");

print " <1> Target capture using multiple oligos\n";
print " <2> Exome capture\n\n"; 

$answer = <STDIN>;
chomp $answer;
if (substr($answer,0,1) eq "1" ){$exome_sequencing = "false"}
if (substr($answer,0,1) eq "2" ){$exome_sequencing = "true"}



#############################################
# Ask if duplicates should be removed       #
#############################################
&print_message("Would you like to remove duplicates from the BAM files?","input");

print " <1> Yes\n";
print " <2> No     [DEFAULT]\n"; 

$answer = <STDIN>;
chomp $answer;

if ($answer eq ""){$answer = "2"}

if (substr($answer,0,1) eq "1" ){$remove_duplicates = "yes"}
if (substr($answer,0,1) eq "2" ){$remove_duplicates = "no"}


if ($remove_duplicates eq "yes")
{
	#####################################
	# Ask if the data is SE or PE data? #
	#####################################
	&print_message("Do you have a Single-end or Paired-end dataset?","input");
	print "  <1>  Single-end\n";
	print "  <2>  Paired-end    [DEFAULT]\n\n";

	$answer = <STDIN>;
	chomp $answer;
	if ($answer eq ""){$answer = "2"}

	if (substr($answer,0,1) eq "1" ){$paired_ends_type = "SE"}
	if (substr($answer,0,1) eq "2" ){$paired_ends_type = "PE"}
}

##############################################################
# Get header from BAM file to use for baits and target files #
##############################################################
&get_header_from_bam_file;


################################################################################
# Create a new target interval file and add the header from the BAM file first #
################################################################################
&create_new_baits_file;


################################################################################
# Create a new target interval file and add the header from the BAM file first #
################################################################################
&create_target_interval_file;


##################################################################
# Write re-formatted baits file                                  #
##################################################################
&reformat_baits_file;


######################################################
# Get guessed target regions from info in Baits file #
######################################################
if ($exome_sequencing eq "false") 
{&get_target_regions_from_baits_file;}
else
{&get_target_regions_from_baits_file_exome;}


####################################################
# Ask for target regions required                  #
####################################################
&get_target_regions_from_user;


#####################################################
# Add chosen target regions to Target Interval File #
#####################################################
&add_target_regions_to_file;

print "New Target Intervals file $target_interval_file created\n\n";


###########################################################
# Print the details entered to the screen to check        #
###########################################################

&print_message("PLEASE CHECK THESE DETAILS CAREFULLY!!","message");


&print_both("Baits File:                     \t$baits_file\n\n");
&print_both("Baits File (re-formatted):      \t$baits_file_new\n\n");	
&print_both("Target Interval File            \t$target_interval_file\n\n");

&print_both("Paired ends type:               \t$paired_ends_type\n\n");



if ($exome_sequencing eq "false")
{
	&print_both("Regions for target interval file:\n\n");

	for ($count = 1; $count<= $no_of_target_regions; $count++)
	{
		&print_both("\t$target_region_array[$count]\n");
	}
}
if ($exome_sequencing eq "true")
{
	&print_both("Some example regions for target interval file: (there are $no_of_target_regions regions)\n\n");

	for ($count = 1; $count<= $no_of_target_regions; $count++)
	{
		if ($count % 10000 == 0) {&print_both("\t$target_region_array[$count]\n");}
	}
}


if ($remove_duplicates eq "no"){&print_both("\n\nDuplicates will not be removed\n");}
if ($remove_duplicates eq "yes"){&print_both("\n\nDuplicates will  be removed\n");}

&print_message("RUN TIME MAY BE SEVERAL HOURS.  It might be best to use a unix detached screen (e.g. screen -S example)","message"); 

&pause;


#####################################################################
# Start analysis of bam files looping through array @bam_file_array #
#####################################################################

&print_message("Processing the BAM files","message");


for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	$bam_file = $bam_file_array_with_path[$list_count];
	$prefix = &get_prefix_without_path($bam_file);

	$aligned_sorted_bam =$bam_file_array_with_path[$list_count];


	$duplicates_removed_bam = "$run_title"."_"."$prefix"."_duplicates_removed.bam";
	$capture_efficiency_file = "$run_title"."_"."$prefix"."_capture_efficiency.out";

	&print_message("The BAM file currently being processed is $bam_file","message");


	##############################################################
	# Removal of PCR duplicates                                  #
	##############################################################
	if ($remove_duplicates eq "yes")
	{
		$duplicates_removed_bam = &remove_duplicates_with_rmdup ($aligned_sorted_bam, $duplicates_removed_bam);
	}

	##############################################################
	# Run CalculateHsMetrics for the bam file                    #
	##############################################################
	if ($remove_duplicates eq "no")
	{
		&run_unix_command("java -jar -Xmx8g /opt/picard/CalculateHsMetrics.jar I=$aligned_sorted_bam O=$capture_efficiency_file BI=$baits_file_new  TI=$target_interval_file VALIDATION_STRINGENCY=LENIENT");
	}
	if ($remove_duplicates eq "yes")
	{
		&run_unix_command("java -jar -Xmx8g /opt/picard/CalculateHsMetrics.jar I=$duplicates_removed_bam O=$capture_efficiency_file BI=$baits_file_new  TI=$target_interval_file VALIDATION_STRINGENCY=LENIENT");
	}

	######################################
	# Move into results folder           #
	######################################
	print "\nMoving results from $bam_file into new directory....\n\n";

	&run_unix_command("mv $capture_efficiency_file $results_folder/$capture_efficiency_file");


	###################################################
	# Read from individual capture efficiency files   #
	# so it can be added to the summary file          #
	###################################################
	open( INFILE, "<$results_folder/$capture_efficiency_file" ) or die( "Couldn't open file $capture_efficiency_file\n" );
	@capture_array = <INFILE>;
	close( INFILE );

	for $capture_array_count ( 0 .. scalar @capture_array ) 
	{
		chomp $capture_array[$capture_array_count];

		# Target the heading line, starting with BAIT_SET
		if ($capture_array[$capture_array_count] =~ m/BAIT_SET\t/ ) 
		{					 	
			@target_line = split (/\t/, $capture_array[$capture_array_count]);		# Split the next line into a new array
			for $h ( 0 .. $#target_line ) 
			{
				chomp $target_line[$h];
				$summaryarray_2d[0][$h] = $target_line[$h];
			}

			@target_line = split (/\t/, $capture_array[$capture_array_count + 1]);		# Split the next line into a new array
			$no_of_rows = scalar @target_line;
			for $h ( 0 .. $#target_line ) 
			{
				chomp $target_line[$h];
				$summaryarray_2d[$list_count][$h] = $target_line[$h];
			}
			last;
		}
	} # capture_array loop

} # bam_file_loop



#############################################################
# Print results to summary results file (for all BAM files) #
#############################################################

open(OUTFILE, ">>$results_folder/$coverage_summary_file");

print OUTFILE "Results from picard CalculateHsMetrics\n";

print OUTFILE "Bait file\t$baits_file\n";

print OUTFILE "Remove duplicates\t$remove_duplicates\n";
for ($count = 1; $count<= $no_of_target_regions; $count++)
{
	print OUTFILE "Region $count\t$target_region_array[$count]\n";
}

print OUTFILE "METRIC\t";

#BAM file names
for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print OUTFILE "$bam_file_array[$list_count]\t";
}

print OUTFILE "\n";
for ($row_count=1;$row_count<$no_of_rows;$row_count++)
{
	for ($list_count=0;$list_count<=$no_of_files;$list_count++)
	{
		print OUTFILE "$summaryarray_2d[$list_count][$row_count]\t";
	} # file loop

	$explanation = &get_explanation($summaryarray_2d[0][$row_count]);

	print OUTFILE "$explanation\n";
} # row loop
close OUTFILE;

#RUN TIMER END

my $end_run = time();
my $run_time = $end_run - our $start_run;

open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $email_address\n";
print MAIL "From: $from\n";
print MAIL "Subject: NGS capture analysis completed ($run_title)\n\n";
## Mail Body
print MAIL "Your target capture efficiency analysis ($run_title) is complete\n\n";
print MAIL "All the results are in the folder $results_folder\n\n";

if ($no_of_files == 1)
{print MAIL "Capture efficiency file for BAM file:        \t$capture_efficiency_file\n";}
else
{print MAIL "Capture efficiency file for last BAM file:        \t$capture_efficiency_file\n";}

print MAIL "Coverage summary file (all BAM files combined):   \t$coverage_summary_file\n\n";

print MAIL "For pipeline details see $screen_log and $command_log\n\n";
print MAIL "Run time : $run_time seconds\n";
printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
close(MAIL);

&print_message("NGS CAPTURE ANALYSIS COMPLETE!","message");

print "Run time : $run_time seconds\n";
printf "%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];

print "\nAll the results are in the folder $results_folder\n\n";

if ($no_of_files == 1)
{
	print  "Capture efficiency file for BAM file:        \t$capture_efficiency_file\n\n";
	print  "Coverage summary file for single BAM file:   \t$coverage_summary_file\n\n";
}
else
{
	print  "Capture efficiency file for last BAM file:        \t$capture_efficiency_file\n\n";
	print  "Coverage summary file (all BAM files combined):   \t$coverage_summary_file\n\n";
}

print "Please check $command_log and $screen_log for full details of the run\n\n";


#################################################
# Move log file to the results folder
################################################

$command = "mv $screen_log $results_folder/$screen_log";
system("$command");
$command = "mv $command_log $results_folder/$command_log";
system("$command");


###########################################
# Delete files which are no longer needed #
###########################################
&delete_file("$header_temp_file");
&delete_file("$baits_file_new");
&delete_file("$target_interval_file");

if ($remove_duplicates eq "yes")
{
	for ($list_count=0;$list_count<=$no_of_files;$list_count++)
	{
		$bam_file = $bam_file_array[$list_count];
		$prefix = &get_prefix_without_path($bam_file);
		$duplicates_removed_bam = "$run_title"."_"."$prefix"."_duplicates_removed.bam";

		&delete_file("$duplicates_removed_bam");
	} 
}


#############################
# Turn off logging          #
#############################
close(STDOUT);


######################################
# Subroutine to stop at a point      #
######################################

sub pause
{
	my $proceed = "";
	print "\nPlease press enter to proceed (or 'Q' to quit):	";
	$proceed = <STDIN>;
}




######################################
# Subroutine to print screen message #
######################################

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


########################################################
# Open the list file to get the list of BAM file names #
########################################################
sub get_bam_files_from_list_file
{
	open (LIST, "$list_file") || die "Cannot open $list_file";
	$list_count=0;

	print "\n";

	while ($bam_file_in_list_file = <LIST> ) 
	{
		chomp $bam_file_in_list_file;
		
		$list_count=$list_count + 1;

		$slash_pos = index($bam_file_in_list_file,"/");

		#############################################
		# If path is specified in the file name     #
		#############################################
		if ($slash_pos > -1)
		{
			$bam_file_array_with_path[$list_count]=$bam_file_in_list_file;
			$bam_file = &get_prefix_without_path($bam_file_in_list_file);
			$bam_file = $bam_file.".bam";
			$bam_file_array[$list_count] = $bam_file;

			print "$bam_file_in_list_file\tPath specified in the file name\n";
		}

		#############################################
		# If path is not specified in the file name #
		#############################################
		if (index($bam_file_in_list_file,"/") == -1)
		{
			####################################
			# If path is specified by user     #
			####################################
			if ($path_to_bam_file ne "")
			{
				$bam_file_array_with_path[$list_count]=$path_to_bam_file.$bam_file_in_list_file;
				$bam_file_array[$list_count] = $bam_file_in_list_file;

				print "$bam_file_in_list_file\tPath specified separately by user\n";
			}

			####################################
			# If path is NOT specified by user #
			####################################
			if ($path_to_bam_file eq "")
			{
				$bam_file_array_with_path[$list_count]=$bam_file_in_list_file;
				$bam_file_array[$list_count] = $bam_file_in_list_file;

				print "$bam_file_in_list_file\tNo path is specified\n";
			}

		}

		#############################################
		# If path is specified in the file name     #
		# and by user then warn user                #
		#############################################
		if ((index($bam_file_in_list_file,"/") > -1) && ($path_to_bam_file ne ""))
		{
			&print_message("You mustn't specify the path in the input file and also as a user-entered path","warning");

			print "User entered path: $path_to_bam_file\n\n";
			print "File name in input file: $bam_file\n\n";

			exit;
		}

		$prefix = &get_prefix("$bam_file");

		if (! -e $bam_file_array_with_path[$list_count])
		{
			print "\t >>>>>>>>>>   $bam_file_array_with_path[$list_count]\tNOT FOUND\n";
			$some_files_not_found = "true";
		}
			
	}

	close LIST;

	$no_of_files = $list_count;

	print "\nNo. of files: $no_of_files\n\n";

	if ($some_files_not_found eq "true")
	{
		&print_message("Not all the BAM files were found","warning");
		print "Please check your list and try again\n\n";
		exit;
	}
} # get_bam_files_from_list_file


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

###########################################################
# Subroutine to get filename before file type suffix      #
# (i.e. prefix) (e.g. to get "test" from "test.fasta")    #
# New version using rindex rather than index.  This means #
# it can deal with files like filename.something.vcf      #
###########################################################

sub get_prefix_without_path
{
	my $_filename 	= "";
	my $_dot_pos	= 0;
	my $_slash_pos  = 0;


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

	if (rindex($_filename,"/") > 0)
	{
		$_filename = substr($_filename, rindex($_filename,"/") + 1);
	}

	$_filename = $_filename;

} # get_prefix_without_path

#
# Remove duplicates from BAM file with samtools rmdup
sub remove_duplicates_with_rmdup
{
	my $filein = "";
	my $fileout = "";
	$filein = $_[0];
	$fileout = $_[1];

	if  ($paired_ends_type eq "PE")
	{
		&run_unix_command("/opt/samtools/samtools rmdup $filein $fileout");
	}
		
	if  ($paired_ends_type eq "SE")
	{
		&run_unix_command("/opt/samtools/samtools rmdup -s $filein $fileout");
	}

	$fileout = $fileout;

}# remove_duplicates_with_rmdup


sub get_file
{
	my $input_message 	= $_[0];
	my $file_type	 	= $_[1];
	my $file_sub		= "";

	if ($file_type eq ""){$file_type = "*"}

	until ((-e $file_sub) || (lc $file_sub eq "q"))
	{
		&print_message("$input_message","input");
		print "> ";

		$file_sub = <STDIN>;
		chomp $file_sub;

		if ($file_sub eq ""){$file_sub = "exome_baits.bed"}


		if ($file_sub eq "ls") 
		{
			if (($path_to_bam_file ne "") && ($file_type eq ".bam"))
			{
				print "\n";system ("ls $path_to_bam_file*"."$file_type");
			}
			else
			{
				print "\n";system ("ls *"."$file_type");
			}
		}

		if (($file_sub ne "ls")  && (lc $file_sub ne "q"))
		{

			if (!-e $file_sub){print "/n  >>>>>>>>>>>>>  ERROR.  File $file_sub cannot be found <<<<<<<<<<<\n\n";}

			if (-e $file_sub)
			{
				$file_sub = $file_sub;
			}
		} # not ls
	} # until loop

	print "\n";
	$command = "dos2unix $file_sub";
	if ($file_type eq ".txt"){system("$command")}

	$file_sub = $file_sub;
}


###############################################
# Add target regions to target intervals file #
###############################################
sub add_target_regions_to_file 
{
	for ($count = 1; $count<= $no_of_target_regions; $count++)
	{
		$region = $target_region_array[$count];

		@item = split (/:/, $region);
		$chromosome = $item[0];
		$position_string = $item[1];

		@item = split (/-/, $position_string);
		$region_start = $item[0];
		$region_end = $item[1];

		if ($chromosome eq "chr39"){$chromosome = "chrX"}

		print TARGETS "$chromosome\t$region_start\t$region_end\t+\tRegion $count\n";
	}
	close TARGETS;

	print "Target regions added to file\n\n";
} #add_target_regions_to_file

sub get_header_from_bam_file 
{
	&run_unix_command("samtools view -H $bam_file_array_with_path[1] > $header_temp_file");

	#Read this into header_file_array
	open (HEADER, "$header_temp_file" ) or die "Couldn't open $header_temp_file";
		@header_file_array = <HEADER>;
	close HEADER;

} #get_header_from_bam_file

sub create_new_baits_file 
{
	#Open new Baits file and write headers to it
	open ( BAITS_NEW, ">$baits_file_new" ) or die "Couldn't create new BAITS file $baits_file_new\n";
	foreach $single_line (@header_file_array)
	{
		print BAITS_NEW $single_line;
	}
} #create_new_baits_file

sub create_target_interval_file 
{
	$no_of_header_lines = 0;

	#Open new target Intervals file and write headers to it
	open ( TARGETS, ">$target_interval_file" ) or die "Couldn't create new BAITS file $target_interval_file\n";
	foreach $single_line (@header_file_array)
	{
		print TARGETS $single_line;

		$no_of_header_lines = $no_of_header_lines + 1;
	}
} # create_target_interval_file



##############################################################
# Open old Baits file and make a new re-formatted Baits file #
##############################################################
sub reformat_baits_file 
{
	open (BAITS, "<$baits_file" ) or die "Couldn't open baits file $baits_file";
	@baits_file_array = <BAITS>;
	close BAITS;

	&print_message("Creating new re-formatted BAITS file...","message");
	$chromosome = "0";


	####################################################
	# Set up arrays for smallest and largest positions #
	####################################################
	for ($count = 1; $count<40; $count++)
	{
		$smallest_position[$count] = 1000000000;
		$largest_position[$count] = 1;
	}

	foreach $single_line (@baits_file_array)
	{
		chomp $single_line;
		@item=split(/\t/,$single_line);

		$item_array_size = scalar @item;

		####################################################
		# Only read the line if it doesn't start with an @ #
		####################################################
		if (index($single_line,"@") == -1)
		{
			print "Line doesn't contain @  $single_line\n";
			print "Item array size = $item_array_size\n";

			if ($item_array_size > 5)
			{
				$chromosome = $item[0];
				$region_start = $item[1];
				$region_end = $item[2];
				$bait_name = $item[3];
				$plus = $item[5];

				if ($chromosome eq "chr39"){$chromosome = "chrX"}

				print BAITS_NEW "$chromosome\t$region_start\t$region_end\t+\t$bait_name\n";

				print  "Writing to BAITS_NEW: $chromosome\t$region_start\t$region_end\t+\t$bait_name\n";

				$chromosome_number = substr($chromosome,3,99);
				if ($chromosome_number eq "X"){$chromosome_number = 39} # dog only

				############################################################
				# Check. Baits file should be this format:                 #
				# chr20	015543005	015543125	BI426399621_0	1000	+  #
				############################################################
				if ($plus ne "+")
				{
					&print_message("The sixth column of the unformatted BAITS file should be a plus sign '+'","warning");

					print "The format of each line should be:\n\n";
					print "   chr20      015543005     015543125     BI426399621_0     1000     +\n\n";

					print "Your line looks like this:\n\n";

					print "$single_line\n\n";

					print "If the unformatted BAITS file is not in this format the program won't work\n\n";
					exit;
				}

				############################################################
				# Store smallest and largest positions for each chromosome #
				############################################################
				if ($chromosome_number > 0) 
				{
					if ($region_start < $smallest_position[$chromosome_number])
					{
						$smallest_position[$chromosome_number] = $region_start;
					}
					if ($region_end > $largest_position[$chromosome_number])
					{
						$largest_position[$chromosome_number] = $region_end;
					}
				} # if chr > 0

			}# If array size > 5

		}# if line doesn't contain '@'

	} # for each baits file array

	close BAITS_NEW;

	print "\tRe-formatted Baits file: $baits_file_new created\n\n";

} # reformat_baits_file


##################################################################
# Look at Baits file and guess which target regions are required #
##################################################################
sub get_target_regions_from_baits_file 
{

	&print_message("Choose your target regions","message");

	$target_region_count =0;

	for ($count = 1; $count<40; $count++)
	{
		if ($largest_position[$count] > $smallest_position[$count])
		{
			# Save in array
			$target_region_count = $target_region_count + 1;
			$target_region_array[$target_region_count] = "chr$count:$smallest_position[$count]-$largest_position[$count]";
		}
	}

	$no_of_target_regions = $target_region_count;

	print "These are the $no_of_target_regions regions found from the Baits file:\n\n";

	for ($count = 1; $count<= $no_of_target_regions; $count++)
	{
		print "\t$target_region_array[$count]\n";
	}
} # get_target_regions_from_baits_file


##################################################################
# EXOME                                                          #
# Look at Baits file and guess which target regions are required #
##################################################################
sub get_target_regions_from_baits_file_exome 
{

	&print_message("Getting exons from re-formatted baits file $baits_file_new","message");

	###############################################################
	# Open reformatted baits file and store all lines in an array #
	###############################################################
	open (BAITS_NEW, "<$baits_file_new" ) or die "Couldn't open reformatted baits file $baits_file_new";
	@baits_file_new_array = <BAITS_NEW>;
	close BAITS_NEW;

	$no_of_exons = scalar @baits_file_new_array;

	$target_region_count =0;

	for ($exon_count = 1; $exon_count<=$no_of_exons; $exon_count++)
	{
		$single_line = $baits_file_new_array[$exon_count];

		if (index($single_line,"chr") == 0)
		{

			@item=split(/\t/,$single_line);

			$item_array_size = scalar @item;

			if ($item_array_size == 5)
			{
				$chromosome = $item[0];
				$region_start = $item[1];
				$region_end = $item[2];

				# Save in array
				if (index($chromosome,"chr") == 0)
				{
					$target_region_count = $target_region_count + 1;
					$target_region_array[$target_region_count] = "$chromosome:$region_start-$region_end";
				}
			}

		} # line has 'chr' in it.

	} # exon_count loop

	$no_of_target_regions = $target_region_count;

	print "These are some target regions (actually exons) found from the reformatted Baits file:\n\n";

	for ($exon_count = 1; $exon_count< $no_of_exons; $exon_count++)
	{
		if ($exon_count % 10000 == 0){print "\t$exon_count: $target_region_array[$exon_count]\n"}
	}
} # get_target_regions_from_baits_file_exome


################################
# Get target regions from user #
################################
sub get_target_regions_from_user 
{
	
	&print_message("Do you want to use these regions as they are?","input");

	print "  <1> Yes, use these regions     [DEFAULT]\n";
	print "  <2> No, edit these regions\n";
	print "  <3> Enter new regions\n\n";

	$answer = <STDIN>;
	chomp $answer;
	if ($answer eq ""){$answer = "1"}

	if (substr($answer,0,1) eq "1" ){$use_default_regions = "yes"; $region_choice = 1}
	if (substr($answer,0,1) eq "2" ){$use_default_regions = "no"; $region_choice = 2}
	if (substr($answer,0,1) eq "3" ){$use_default_regions = "no"; $region_choice = 3}

	if ($region_choice == 2)
	{
		########################################################
		# Define the chromosome, start and end for each region #
		########################################################
		$target_region_count = 0;
		print "\nEnter regions in the form chr15:23100455-245860120 ('q' to quit entereing regions)\n\n";

		for $count ( 1 .. 100 )
		{
			print "Region $count:  (default = $target_region_array[$count])  ";
			$region = <STDIN>;
			chomp $region;
			if ($region eq ""){$region = $target_region_array[$count]} # use default
			if (lc $region eq "q"){last;}
			# Save in array
			$target_region_count = $target_region_count + 1;
			$target_region_array[$target_region_count] = $region;
		}
		$no_of_target_regions = $target_region_count;

	} # region_choice ==3 . Enter new regions


	if ($region_choice == 3)
	{
		########################################################
		# Define the chromosome, start and end for each region #
		########################################################
		$target_region_count = 0;
		print "Enter regions in the form chr15:23100455-245860120 ('q' to quit entering regions)\n\n";

		for $count ( 1 .. 100 )
		{
			print "Region $count :  ";
			$region = <STDIN>;
			chomp $region;

			if (($region eq "") || ($region eq "q")){last;}
			# Save in array
			$target_region_count = $target_region_count + 1;
			$target_region_array[$target_region_count] = $region;
		}
		$no_of_target_regions = $target_region_count;

	} # region_choice ==3 . Enter new regions

} # get_target_regions_from_user

#############################################
# Subroutine to delete files                #
#############################################

sub delete_file
{
	my $file_to_be_deleted = "";	
	$file_to_be_deleted = $_[0];
	
	if (! -e "$file_to_be_deleted")
	{
		print "\n$file_to_be_deleted could not be found to be deleted\n";
		print COMMAND_LOG "\n$file_to_be_deleted could not be found to be deleted\n";
	}
	
	if (-e "$file_to_be_deleted")
	{
		$command = "rm  $file_to_be_deleted";
		system("$command");
		print COMMAND_LOG "\n$file_to_be_deleted was deleted\n";
	}
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

sub get_explanation
{
	my $metric = $_[0];

	if ($metric eq "BAIT_SET") {$explanation = " The name of the bait set used in the hybrid selection."}
	if ($metric eq "GENOME_SIZE") {$explanation = " The number of bases in the reference genome used for alignment."}
	if ($metric eq "BAIT_TERRITORY") {$explanation = " The number of bases which have one or more baits on top of them."}
	if ($metric eq "TARGET_TERRITORY") {$explanation = " The unique number of target bases in the experiment where target is usually exons etc."}
	if ($metric eq "BAIT_DESIGN_EFFICIENCY") {$explanation = " Target terrirtoy / bait territory. 1 == perfectly efficient, 0.5 = half of baited bases are not target."}
	if ($metric eq "TOTAL_READS") {$explanation = " The total number of reads in the SAM or BAM file examine."}
	if ($metric eq "PF_READS") {$explanation = " The number of reads that pass the vendor's filter."}
	if ($metric eq "PF_UNIQUE_READS") {$explanation = " The number of PF reads that are not marked as duplicates."}
	if ($metric eq "PCT_PF_READS") {$explanation = " PF reads / total reads. The percent of reads passing filter."}
	if ($metric eq "PCT_PF_UQ_READS") {$explanation = ' PF Unique Reads / Total Reads.'}
	if ($metric eq "PF_UQ_READS_ALIGNED") {$explanation = " The number of PF unique reads that are aligned with mapping score > 0 to the reference genome."}
	if ($metric eq "PCT_PF_UQ_READS_ALIGNED") {$explanation = " PF Reads Aligned / PF Reads."}
	if ($metric eq "PF_UQ_BASES_ALIGNED") {$explanation = " The number of bases in the PF aligned reads that are mapped to a reference base. Accounts for clipping and gaps."}
	if ($metric eq "ON_BAIT_BASES") {$explanation = " The number of PF aligned bases that mapped to a baited region of the genome."}
	if ($metric eq "NEAR_BAIT_BASES") {$explanation = " The number of PF aligned bases that mapped to within a fixed interval of a baited region, but not on a baited region."}
	if ($metric eq "OFF_BAIT_BASES") {$explanation = " The number of PF aligned bases that mapped to neither on or near a bait."}
	if ($metric eq "ON_TARGET_BASES") {$explanation = " The number of PF aligned bases that mapped to a targeted region of the genome."}
	if ($metric eq "PCT_SELECTED_BASES") {$explanation = " On+Near Bait Bases / PF Bases Aligned."}
	if ($metric eq "PCT_OFF_BAIT") {$explanation = " The percentage of aligned PF bases that mapped neither on or near a bait."}
	if ($metric eq "ON_BAIT_VS_SELECTED") {$explanation = " The percentage of on+near bait bases that are on as opposed to near."}
	if ($metric eq "MEAN_BAIT_COVERAGE") {$explanation = " The mean coverage of all baits in the experiment."}
	if ($metric eq "MEAN_TARGET_COVERAGE") {$explanation = " The mean coverage of targets that received at least coverage depth = 2 at one base."}
	if ($metric eq "PCT_USABLE_BASES_ON_BAIT") {$explanation = " The number of aligned, de-duped, on-bait bases out of the PF bases available."}
	if ($metric eq "PCT_USABLE_BASES_ON_TARGET") {$explanation = " The number of aligned, de-duped, on-target bases out of the PF bases available."}
	if ($metric eq "FOLD_ENRICHMENT") {$explanation = " The fold by which the baited region has been amplified above genomic background."}
	if ($metric eq "ZERO_CVG_TARGETS_PCT") {$explanation = " The number of targets that did not reach coverage=2 (Read Depth=2) over any base."}
	if ($metric eq "FOLD_80_BASE_PENALTY") {$explanation = " The fold over-coverage necessary to raise 80% of bases in 'non-zero-cvg' targets to the mean coverage level in those targets."}
	if ($metric eq "PCT_TARGET_BASES_2X") {$explanation = " The percentage of ALL target bases achieving 2X or greater coverage (Read Depth)."}
	if ($metric eq "PCT_TARGET_BASES_10X") {$explanation = " The percentage of ALL target bases achieving 10X or greater coverage (Read Depth)."}
	if ($metric eq "PCT_TARGET_BASES_20X") {$explanation = " The percentage of ALL target bases achieving 20X or greater coverage (Read Depth)."}
	if ($metric eq "PCT_TARGET_BASES_30X") {$explanation = " The percentage of ALL target bases achieving 30X or greater coverage (Read Depth)."}
	if ($metric eq "PCT_TARGET_BASES_40X") {$explanation = " The percentage of ALL target bases achieving 40X or greater coverage (Read Depth)."}
	if ($metric eq "PCT_TARGET_BASES_50X") {$explanation = " The percentage of ALL target bases achieving 50X or greater coverage (Read Depth)."}
	if ($metric eq "PCT_TARGET_BASES_100X") {$explanation = " The percentage of ALL target bases achieving 100X or greater coverage (Read Depth)."}
	if ($metric eq "HS_LIBRARY_SIZE") {$explanation = " The estimated number of unique molecules in the selected part of the library."}
	if ($metric eq "HS_PENALTY_10X") {$explanation = " The 'hybrid selection penalty' incurred to get 80% of target bases to 10X. This metric should be interpreted as"}
	if ($metric eq "HS_PENALTY_20X") {$explanation = " The 'hybrid selection penalty' incurred to get 80% of target bases to 20X. This metric should be interpreted as"}
	if ($metric eq "HS_PENALTY_30X") {$explanation = " The 'hybrid selection penalty' incurred to get 80% of target bases to 30X. This metric should be interpreted as"}
	if ($metric eq "HS_PENALTY_40X") {$explanation = " The 'hybrid selection penalty' incurred to get 80% of target bases to 40X. This metric should be interpreted as"}
	if ($metric eq "HS_PENALTY_50X") {$explanation = " The 'hybrid selection penalty' incurred to get 80% of target bases to 50X. This metric should be interpreted as"}
	if ($metric eq "HS_PENALTY_100X") {$explanation = " The 'hybrid selection penalty' incurred to get 80% of target bases to 100X. This metric should be interpreted as"}
	if ($metric eq "AT_DROPOUT") {$explanation = " A measure of how undercovered <= 50% GC regions are relative to the mean. For each GC bin [0..50] we calculate a = % of target territory, and b = % of aligned reads aligned to these targets. AT DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total reads that should have mapped to GC<=50% regions mapped elsewhere."}
	if ($metric eq "GC_DROPOUT") {$explanation = " A measure of how undercovered >= 50% GC regions are relative to the mean. For each GC bin [50..100] we calculate a = % of target territory, and b = % of aligned reads aligned to these targets. GC DROPOUT is then abs(sum(a-b when a-b < 0)). E.g. if the value is 5% this implies that 5% of total reads that should have mapped to GC>=50% regions mapped elsewhere"}

	$explanation = $explanation;

}#get_explanation