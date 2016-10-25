 #!/usr/bin/perl -w

############################################################################
#									                                             
#	coverage_plot  		                                                            
#									                                       
#	THIS PERL SCRIPT WILL compare coverage of cases and controls           
#									                                       
#   It works in three stages:
#
# Stage 1: Runs GATK/DepthOfCoverage.jar to get the basic coverage data
#          from your BAM files.
#
# Stage 2: The perl script then makes two files from this:
#          File 1: a coverage file with only the most significant 
#                  differences between cases and controls.
#          File 2: a SINE data output file for the next stage.
#
# Stage 3: The sine data file from stage 2 is analysed to look for
#          possible inserted SINEs
#
# Output files are 1) A list of possible SINEs
#                  2) A file to plot SINEs
#
############################################################################

#############################
# Mike Boursnell Dec 2012   #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use Getopt::Std ;
use File::Basename ;
use Statistics::Distributions;
use Statistics::Basic;

#Version
my $version						= "13";

# Constants
my $gatk_directory				= "gatk";
my $use_samtools				= "no";
my $no_output_columns			= 4; # This is true for GATK DepthOfCoverage files
my $show 						= "false"; # For debugging 
my $show_during					= "false"; # For debugging 
my $show_sine					= "false"; # For debugging
my $include_individual_columns	= "true";
my $verbose_out					= "false";
my $debugging_mode				= "true";
my $normalise_coverage			= "true";

#Thresholds for coverage_all.out file #
# There is a small amount of filtering for this file if you set filter_large = "true" only
my $coverage_threshold_for_all	= 5; # Used for filtering when creating coverage_all.out
my $filter_large				= "false"; # Do you filter when making coverage_all.out file

# Thresholds for coverage_assoc.out         #
# Only filtered positions go into this file #
my $t_test_threshold			= 3.5;
my $covar_threshold				= 0.5; # threshold for coverage_assoc.out
my $coverage_threshold_for_assoc= 4;   # if mean coverage is less than this then it is not written to the ASSOC file
my $filter_assoc				= "true"; # Do you want to filter when making the "assoc" file
my $chi_squared_threshold		= 0; # Not used anymore

# Look for SINEs - thresholds for filtering #
my $sine_stringency				= "";
my $flat_before_threshold		= 4; # the number of positions before the SINE where coverage is flat
my $flat_after_threshold		= 4; # the number of positions after the SINE where coverage is flat
my $sine_length_upper_threshold	= 17;
my $sine_length_lower_threshold	= 10;
my $sine_score_threshold		= 9.7; # threshold for filtering SINEs on SINE SCORE
my $covar_whole_window_threshold= 0.26; # coefficient of variance threshold for possible_sines.out. Less than this
my $crc_flat_threshold			= 0.08; # crc (coverage ratio change) must be less than this value in regions before and after the cliffs
my $crc_flat_during_threshold 	= 0.07; # change must be less than this value between the cliffs
my $crc_cliff_threshold			= 0.15; # change must be more than this for the "cliff" at the start and end of the SINE
my $cliff_difference_threshold	= 0.14; # cliff up and cliff down must nearly cancel out. Less than this
my $mean_sine_coverage_threshold= 6; # mean threshold coverage in a SINE

my $saved_by_sine_score			= "false";
my $saved_by_inclusion			= "false";
my $saved_by_strict				= "false";
my $saved_by_sine_score_count	= 0;
my $saved_by_inclusion_count	= 0;
my $saved_by_strict_count		= 0;


#SINES
my $flat_before_count			= 0;
my $flat_during_count			= 0;
my $flat_after_count			= 0;
my $window_size					= 0;
my $check_count					= 0;
my $plot_count					= 0;
my $cliff_1_count				= 0; # The "cliff 1" is the start edge found when DP changes suddenly
my $cliff_2_count				= 0; # The "cliff 2" is the end edge found when DP changes back
my $sine_repeat_length			= 0; # length of repeat at end of SINEs which shows up as an increase in coverage
my $end_of_sine					= 0; # position at end of raised coverage (or lowered)
my $cliff_1_size				= 0;
my $cliff_2_size				= 0;
my $cliff_1_size_abs			= 0;
my $cliff_2_size_abs			= 0;
my $cliff_difference			= 0;
my $coverage_ratio_change_signed = 0; # still has + or - sign (not ABS as $coverage_ratio_change)
my $sine_score					= 0;
my $total_change_before			= 0;
my $total_change_during			= 0;
my $total_change_after  		= 0;
my $mean_change_before			= 0;
my $mean_change_during			= 0;
my $mean_change_after  			= 0;
my $inclusion_score				= 0;
my $inclusion_threshold 		= 0; # Calculated from various constants above

# SINE filtering counts #
my $sine_coverage_filtered_count 		= 0;
my $sine_cliff_filtered_count 			= 0;
my $sine_cliffdiff_filtered_count 		= 0;
my $sine_repeat_length_filtered_count 	= 0;
my $sine_covar_filtered_count 			= 0;
my $sine_flat_filtered_count 			= 0;
my $sine_flat_during_filtered_count 	= 0;

		
# Other variables
my $list_count					= 0;
my $checked_count				= 0;
my $no_of_files					= 0;
my $no_of_lines_coverage_all_out= 0;
my $no_of_lines_large_out		= 0;
my $no_of_lines_sine_data		= 0;
my $array_size					= 0;
my $file_count					= 0;
my $no_of_blocks				= 0;
my $block_count					= 0;
my $affected_total_coverage		= 0;
my $normal_total_coverage		= 0;
my $affected_total_coverage_normalised = 0;
my $normal_total_coverage_normalised 	= 0;
my $no_of_affecteds				= 0;
my $no_of_normals				= 0;
my $mean_coverage_affected		= 0;
my $mean_coverage_normal		= 0;
my $mean_coverage_both			= 0;
my $line_count					= 0;
my $bad_line_count				= 0;
my $expected_affected			= 0;
my $expected_normal				= 0;
my $chi_squared				= 0;
my $sum_squares_affected	= 0;
my $sum_squares_normal		= 0;
my $variance_affected		= 0;
my $variance_normal			= 0;
my $sd_affected				= 0;
my $sd_normal				= 0;
my $covar_affected			= 0; # Coefficient of variation (dimensionless so better than SD)
my $covar_normal			= 0;
my $covar_all_samples		= 0;
my $covar_whole_window			= 0; # Sum of Coefficient of variation for the whole window
my $mean_covar_whole_window		= 0; # Mean Coefficient of variation for the whole window
my $covar_during				= 0; # Sum of Coefficient of variation between cliffs
my $mean_covar_during			= 0; # Mean Coefficient of variation between cliffs
my $coverage_ratio				= 0;
my $last_coverage_ratio			= 0;
my $coverage_ratio_change		= 0;
my $chi_filtered_count			= 0;
my $covar_filtered_count		= 0; # count of positions filtered from assoc file because of coeff of variation
my $t_test_filtered_count		= 0; # count of positions filtered from assoc file because of T-test
my $coverage_filtered_count		= 0; # count of positions filtered from assoc file because of mean coverage
my $keep_for_assoc_count		= 0;
my $crc_lose_count				= 0;
my $t_test						= 0;
my $percentage_done				= 0;
my $total_covar					= 0; # mean coefficient of varation across region
my $percentage_inclusion		= 0;
my $percentage_sine_score		= 0;
my $percentage_strict			= 0;
my $percentage_t_test_lost		= 0;
my $percentage_covar_lost		= 0;
my $percentage_coverage_lost	= 0;
my $percentage_kept_all			= 0;
my $percentage_kept_assoc		= 0;
my $keep_for_coverage_all_count	= 0;
my $lose_for_coverage_all_count	= 0;
my $sum_of_means				= 0;
my $mean_of_means				= 0;


# File names
my $config_file				= "";
my $bam_file				= "";
my $list_file				= "";
my $prefix					= "";
my $command					= "";
my $coverage_output_file	= ""; # This is the output file from GATK DepthOfCoverage
my $coverage_mpileup_file	= ""; # This is the output file from samtools mpileup (probably not used)
my $coverage_all_large_out	= ""; # The merged coverage file made with the UNIX paste command (called coverage_all_large.out)
my $coverage_all_out		= ""; # This file is a simplified version of the coverage_all_large.out file
my $output_assoc_file		= ""; # This is a filtered version of coverage_all_out with only significant lines retained
my $output_sine_data_file	= ""; # This is the file with the relevant info for searching for SINE insertions
my $possible_sines_out		= ""; # This is the output file wit the results from the SINE search (possible sines out)
my $possible_sines_plot_data_out = ""; # File holding data for plotting SINE histograms

my $position				= "";
my $coverage				= "";
my $input_string			= "";
my $run_title				= "";
my $sample_name				= "";
my $answer					= "";
my $ref						= "";
my $ref_seq_name			= "";
my $use_defined_region		= "yes";
my $region					= "";
my $mem						= "-Xmx4g";
my $command_log				= "";
my $single_line				= "";
my $read_length				= "";
my $chromosome				= "";
my $GATK_region_string		= "";
my $GATK_validation_stringency = "LENIENT";
my $samtools_region_string	= "";
my $temp_dir_string			= "";
my $tempdir					= "";
my $prefix_name				= "";
my $second_column_found		= "";
my $status					= "";
my $species					= "";
my $position_string			= "";
my $files_to_use			= ""; # Start with BAM files or coverage.out files
my $coverage_string			= "";
my $keep					= ""; # flag to decide what goes into the assoc file
my $keep_chi				= "";
my $keep_t_test				= "true";
my $keep_covar				= "true";
my $keep_coverage			= "true";
my $keep_crc				= "true";
my $sine_found				= "false";
my $cliff_position			= "";


my @bam_file_array			= ();
my @sample_name_array		= ();
my @status_array			= ();
my @item					= ();
my @item_position			= ();
my @coverage_array			= ();
my @total_coverage_array	= ();
my @mean_coverage_array		= ();
my @file_array				= ();
my @weighting_array			= ();

use Term::ANSIColor;
print color 'bold magenta';

print "\n\n";
print "                 ##################################################\n";
print color 'bold white';
print "                                      coverage_plot                    \n";
print " \n";
print "                     Compares coverage between cases and controls        \n";
print "                     Also specifically looks for SINE insertions         \n";
print color 'reset';
print color 'bold magenta';
print "                 ##################################################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  This program uses GATK DepthOfCoverage to get read depths.\n";
print "  Then it compares cases and controls to find significant differences\n";
print "  across the region specified.\n\n\n";

print "  It also looks for SINE insertions.\n\n";

print "  Significant differences are output to a file ending:  \t _coverage_assoc.out.\n";
print "  Possible SINE insertions are output to a file ending: \t _possible_sines.out.\n";
print "  Plotting data for SINEs are output to a file ending:  \t _possible_sines_plot_data.out.\n\n\n";

print color 'reset';


print "~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Which do you want to do?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   <1>  Use BAM files to start analysis from the beginning\n";
print "   <2>  Use 'coverage_all_large.out' file to analyse the depth of coverage figures.\n";
print "   <3>  Use 'coverage_all.out' (filtered) file to analyse the depth of coverage figures.\n";
print "   <4>  Use 'sine_data.out' file to look for SINE elements.\n\n";

$answer=<STDIN>;
chomp $answer;
if (substr($answer,0,1) eq "1" ){$files_to_use = "bam"}
if (substr($answer,0,1) eq "2" ){$files_to_use = "coverage_all_large"}
if (substr($answer,0,1) eq "3" ){$files_to_use = "coverage_all"} # should work but doesn't
if (substr($answer,0,1) eq "4" ){$files_to_use = "sine"}



if ($files_to_use eq "bam")
{
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "Please enter a name for this analysis (with no spaces) \n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
	$run_title = <STDIN>;
	chomp $run_title;
}

if ($files_to_use eq "coverage_all_large")
{
	until ((-e "$run_title"."_coverage_all_large.out") && ($run_title ne ""))
	{
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		print "Please enter the name for this analysis \n";
		print "so that the coverage_all_large file can be found \n";
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
		$run_title = <STDIN>;
		chomp $run_title;
		$coverage_all_large_out = "$run_title"."_coverage_all_large.out";
		$output_sine_data_file = "$run_title"."_sine_data.out";
		
		if ($run_title eq "ls"){print "\n";system ("ls *_coverage_all_large.out");print "\n";}
		if ($run_title ne "ls")
		{
			if (! -e "$coverage_all_large_out"){print "\n\n>>>>>>>>  File $coverage_all_large_out not found.  Try again.  <<<<<<<<\n\n";}
			if (-e "$coverage_all_large_out"){print "\n\nFile $coverage_all_large_out found.\n\n";}
		}
	}
}


if ($files_to_use eq "coverage_all")
{
	until ((-e "$run_title"."_coverage_all.out") && ($run_title ne ""))
	{
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		print "Please enter the name for this analysis   \n";
		print "so that the coverage_all file can be found\n";
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
		$run_title = <STDIN>;
		chomp $run_title;
		$coverage_all_out = "$run_title"."_coverage_all.out";
		$output_sine_data_file = "$run_title"."_sine_data.out";
		
		if ($run_title eq "ls"){print "\n";system ("ls *_coverage_all.out");print "\n";}
		if ($run_title ne "ls")
		{
			if (! -e "$coverage_all_out"){print "\n\n>>>>>>>>  File $coverage_all_out not found.  Try again.  <<<<<<<<\n\n";}
		}
	}
}

if ($files_to_use eq "sine")
{
	until ((-e "$run_title"."_sine_data.out") && ($run_title ne ""))
	{
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
		print "Please enter the name for this analysis\n";
		print "so that the sine_data file can be found\n";
		print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
		$run_title = <STDIN>;
		chomp $run_title;
		$coverage_all_out = "$run_title"."_coverage_all.out";
		$output_sine_data_file = "$run_title"."_sine_data.out";
		$possible_sines_out = "$run_title"."_possible_sines.out";
		if ($run_title eq "ls"){print "\n";system ("ls *_coverage_all.out")}
		if ($run_title ne "ls")
		{
			if (! -e "$output_sine_data_file"){print "\n\n>>>>>>>>  File $output_sine_data_file not found.  Try again.  <<<<<<<<\n\n";}
		}
	}
}


#######################
# Set some file names #
#######################
$coverage_all_out = "$run_title"."_coverage_all.out";
$output_assoc_file = "$run_title"."_coverage_assoc.out";
$output_sine_data_file = "$run_title"."_sine_data.out";
$possible_sines_out = "$run_title"."_possible_sines.out";
$possible_sines_plot_data_out = "$run_title"."_possible_sines_plot_data.out";
$coverage_all_large_out = "$run_title"."_coverage_all_large.out";

if ($files_to_use eq "coverage_all_large")
{
	print "\n\nPre-existing coverage_all_large file found:  \t$coverage_all_large_out\n\n";
}

if ($files_to_use eq "coverage_all")
{
	print "\n\nPre-existing coverage_all file found:  \t$coverage_all_out\n\n";
}

if ($files_to_use eq "sine")
{
	print "\n\nPre-existing coverage file found:   \t$coverage_all_out\n\n";
	print "Pre-existing SINE data file found:  \t$output_sine_data_file\n\n";
}

until ($sine_stringency ne "")
{
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "What level of stringency do you want for searching for SINE elements?\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
	print "   <1>  Normal\n";
	print "   <2>  Low (find more SINEs)\n";
	print "   <3>  High (find fewer SINEs)\n";
	print "   <4>  Help on criteria\n";

	$answer=<STDIN>;
	chomp $answer;
	if (substr($answer,0,1) eq "1" ){$sine_stringency = "normal"}
	if (substr($answer,0,1) eq "2" ){$sine_stringency = "low"}
	if (substr($answer,0,1) eq "3" ){$sine_stringency | "high"}

	if (substr($answer,0,1) eq "4" )
	{
		print "\n\n";
		print "###################################################################\n";
		print "# The following thresholds are used to search for possible SINEs. #\n";
		print "###################################################################\n\n";
		
		print "A SINE is assumed to have flat regions before and after, and a 'cliff' at the start and end.\n\n";

		print "     |||||||||||||| \n";
		print "     |||||||||||||| \n";
		print "     |||||||||||||| \n";
		print "|||||||||||||||||||||||\n";
		print "|||||||||||||||||||||||\n";
		print "|||||||||||||||||||||||\n";
		print "|||||||||||||||||||||||\n";
		print "<FLAT>  <REPEAT>   <FLAT>\n";
		print "     ^            ^    \n";
		print "  <CLIFF>      <CLIFF> \n";
		
		print "\n\n~~~~~~~~~~~~~~~~~~~\n";
		print "FILTERING CRITERIA:\n";
		print "~~~~~~~~~~~~~~~~~~~\n\n";
		
		print "\nFor the following criteria the score must be greater or equal to the threshold:\n\n";

		print "\tFlat before threshold:                       \t$flat_before_threshold\n";
		print "\tFlat after threshold:                        \t$flat_after_threshold\n";
		print "\tMinimum SINE length:                         \t$sine_length_lower_threshold\n";
		print "\tCliff size threshold:                        \t$crc_cliff_threshold\n";
		print "\tMean coverage threshold:                     \t$mean_sine_coverage_threshold\n";
		print "\tSINE score threshold:                        \t$sine_score_threshold\n\n\n";


		print "For the following criteria the score must be less than or equal to the threshold:\n\n";

		print "\tVariance over whole window threshold:        \t$covar_whole_window_threshold\n";
		print "\tMaximum SINE length:                         \t$sine_length_upper_threshold\n";
		print "\tCoverage ratio change in flat sections:      \t$crc_flat_threshold\n";
		print "\tCoverage ratio change in central section:    \t$crc_flat_during_threshold\n";
		print "\tCliff difference threshold:                  \t$cliff_difference_threshold\n\n\n";

		print "\n\nPress return to continue\n\n";
		$answer=<STDIN>;
		
	} # help section
	
} # until $sine_stringency ne ""

##################################
# Set SINE stringency thresholds #
##################################
if ($sine_stringency eq "normal")
{
	# Look for SINEs - thresholds for filtering #
	$flat_before_threshold			= 4; # the number of positions before the SINE where coverage is flat
	$flat_after_threshold			= 4; # the number of positions after the SINE where coverage is flat
	$sine_length_upper_threshold	= 17;
	$sine_length_lower_threshold	= 10;
	$sine_score_threshold			= 9.7; # threshold for filtering SINEs on SINE SCORE
	$covar_whole_window_threshold	= 0.30; # coefficient of variance threshold for possible_sines.out. Less than this
	$crc_flat_threshold				= 0.08; # crc (coverage ratio change) must be less than this value in regions before and after the cliffs
	$crc_flat_during_threshold 		= 0.07; # change must be less than this value between the cliffs
	$crc_cliff_threshold			= 0.15; # change must be more than this for the "cliff" at the start and end of the SINE
	$cliff_difference_threshold		= 0.14; # cliff up and cliff down must nearly cancel out. Less than this
	$mean_sine_coverage_threshold	= 6; # mean threshold coverage in a SINE
}	

if ($sine_stringency eq "low")
{
	# Look for SINEs - thresholds for filtering #
	$flat_before_threshold			= 3; # the number of positions before the SINE where coverage is flat
	$flat_after_threshold			= 3; # the number of positions after the SINE where coverage is flat
	$sine_length_upper_threshold	= 21;
	$sine_length_lower_threshold	= 8;
	$sine_score_threshold			= 9.5; # threshold for filtering SINEs on SINE SCORE
	$covar_whole_window_threshold	= 0.33; # coefficient of variance threshold for possible_sines.out. Less than this
	$crc_flat_threshold				= 0.09; # crc (coverage ratio change) must be less than this value in regions before and after the cliffs
	$crc_flat_during_threshold 		= 0.09; # change must be less than this value between the cliffs
	$crc_cliff_threshold			= 0.14; # change must be more than this for the "cliff" at the start and end of the SINE
	$cliff_difference_threshold		= 0.15; # cliff up and cliff down must nearly cancel out. Less than this
	$mean_sine_coverage_threshold	= 5; # mean threshold coverage in a SINE
}	

if ($sine_stringency eq "high") # Not changed yet. May use as extra low
{
	# Look for SINEs - thresholds for filtering #
	$flat_before_threshold			= 4; # the number of positions before the SINE where coverage is flat
	$flat_after_threshold			= 4; # the number of positions after the SINE where coverage is flat
	$sine_length_upper_threshold	= 17;
	$sine_length_lower_threshold	= 10;
	$sine_score_threshold			= 9.7; # threshold for filtering SINEs on SINE SCORE
	$covar_whole_window_threshold	= 0.26; # coefficient of variance threshold for possible_sines.out. Less than this
	$crc_flat_threshold				= 0.08; # crc (coverage ratio change) must be less than this value in regions before and after the cliffs
	$crc_flat_during_threshold 		= 0.07; # change must be less than this value between the cliffs
	$crc_cliff_threshold			= 0.15; # change must be more than this for the "cliff" at the start and end of the SINE
	$cliff_difference_threshold		= 0.14; # cliff up and cliff down must nearly cancel out. Less than this
	$mean_sine_coverage_threshold	= 6; # mean threshold coverage in a SINE
}	


until (-e "$list_file")
{
	print "\n\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "Please input the name of your disease status file with a list of file names of the BAM files. \n";
	print "NOTE: this file has a second column (tab delimited) to specify Affected, Carrier or Normal.\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
	
	$list_file = <STDIN>;
	chomp $list_file;
	
	if ($list_file eq "")
	{
		$list_file = "test.txt";
	}
	
	if (($list_file ne "ls" ) && (index($list_file,".txt") == -1 )){$list_file = $list_file.".txt"}

	
	if ($list_file eq "ls")
	{
		print "\n";
		system ("ls *.txt");
		print "\n";
	}
	
	
	if ($list_file ne "ls")
	{
		if (! -e "$list_file"){print "\nFile doesn't exist. Try again...    \n";}
	}
}

#############################################
# Make sure the list file is in Unix format #
#############################################

$command = "dos2unix $list_file";
system("$command");


####################################################
# Open the list file to get the list of file names #
####################################################
open (LIST, "$list_file") || die "Cannot open $list_file";
$list_count=1;

while ($single_line = <LIST> ) 
{
	chomp $single_line;
	
	@item=split(/\t/,$single_line);
		
	$array_size = scalar @item;
	
	if ($array_size == 1)
	{
		print "\n\n";
		print "######################\n";
		print "#  INPUT FILE ERROR  #\n";
		print "######################\n\n";
		
		print "There should be a second column containing the statuses - Affected, Carrier or Normal!\n\n";

		exit;
	}
	
	if ($array_size == 2)
	{
		$bam_file = $item[0];
		$status=lc ($item[1]);
		$second_column_found = "true";
		
		if ($status eq "case") {$status="affected"}
		
		if ($status eq "normal") {$status="control"}
		if ($status eq "clear") {$status="control"}
		
		if ($status eq "affected"){$no_of_affecteds = $no_of_affecteds + 1}
		if ($status eq "control") {$no_of_normals = $no_of_normals + 1}
	}
	
	#####################################################
	# Get prefix of bam file (ie omitting .bam suffix ) #
	#####################################################
	$prefix_name = &get_prefix($bam_file);
	
	
	######################################################
	# Add file type suffix .bam if user hasn't added it  #
	######################################################

	if (index($bam_file,".bam") == -1 ){$bam_file = $bam_file.".bam"}


	$bam_file_array[$list_count]=$bam_file;
	$status_array[$list_count]=$status;
	
	$list_count = $list_count + 1;

} # end of while loop for reading list_file


close LIST;

$no_of_files=$list_count - 1;


###################
# List file names #
###################
print "\n\nThere are $no_of_files BAM files in this file of file names.\n\n";

if ($second_column_found eq "true")
{
	print "The second column in the input file has the disease statuses.\n\n";
}

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print "File $list_count	\t$bam_file_array[$list_count]\t\t$status_array[$list_count]\n";
}

print "\n\nNumber of affecteds: \t$no_of_affecteds\n";
print "Number of normals:   \t$no_of_normals\n";

print "\n\nThe order should be Affected - Carrier - Control\n\n";

print "\nIf these file names and statuses look OK, press enter to proceed (or 'Q' to quit):      ";
$answer = <STDIN>;
chomp $answer; 
 
if (lc $answer eq "q"){exit;} 
	

	
################################
# Make the temp java directory #
################################
$tempdir = "$ENV{HOME}/javatempdir";
$temp_dir_string = " -Djava.io.tmpdir=$tempdir";
if (! -e $tempdir)
{
	unless(mkdir $tempdir){die "Unable to create temporary Java directory $tempdir";}
	$temp_dir_string = " -Djava.io.tmpdir=$tempdir";	
}

#############################################
# Use BAM files to start from the beginning #
#############################################
if ($files_to_use eq "bam")
{

	##################################
	# Define data files              #
	##################################

	print "\n\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print " Which reference sequence do you want to use?\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
	print "   Enter 1 for CanFam3\n";
	print "   Enter 2 for CanFam2\n";


	$answer = <STDIN>;
	chomp $answer;

	if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3"; $species = "canis_familiaris"}
	if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/canfam2/canfam2.fasta"; $ref_seq_name = "canfam2"; $species = "canis_familiaris"}


	print "\n\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "You must define the genome region you are using\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

	$use_defined_region = "yes";

	##################################
	# If a defined region is chosen  #
	##################################
	if ($use_defined_region eq "yes")
	{
		print "\nPlease define your region of interest (eg 'chr5:21000000-23000000'):      ";
		$region = <STDIN>;
		chomp $region;
		
		###################################################
		# If only chromosome  is given (e.g. 15 or chr15) #
		###################################################
		if (index($region,":") == -1)
		{
			if (index($region,"chr") == -1){$region = "chr"."$region";}
		}
		
		
		#################################################
		# If full region given chr15:34000000-390000000 #
		#################################################
		if (index($region,":") > -1)
		{
			if (index($region,"chr") > -1)
			{
				$chromosome = substr($region,3,index($region,":")-3);
			}
			if (index($region,"chr") == -1)
			{
				$region = "chr"."$chromosome";
			}
		}
		if (index($region,":") == -1)
		{
			if (index($region,"chr") == 0)
			{
				$chromosome = substr($region,3,99);
			}
		}
		
		if (index($chromosome,"chr") == -1){$chromosome = "chr"."$chromosome";}
		
		$GATK_region_string = "-L $region";
		$samtools_region_string = " -r $region";

	} # if ($use_defined_region eq "yes")

} # end of if $files_to_use eq "bam"

	
$command_log = "coverage_command_log.out";
	
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "COMMAND LOG for coverage.pl $version\n\n";

print COMMAND_LOG "Running with files_to_use = $files_to_use\n\n";

print COMMAND_LOG "Reference sequence:\t$ref\n";
print COMMAND_LOG "Region:            \t$region\n";
print COMMAND_LOG "Input file:        \t$list_file\n\n";
print COMMAND_LOG "List of BAM files:\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print COMMAND_LOG "File $list_count	\t$bam_file_array[$list_count]\t\t$status_array[$list_count]\n";
}
print COMMAND_LOG "\n\n";

print COMMAND_LOG "Thresholds for writing _coverage_all.out:\n\n";

if ($filter_large eq "true")
{
	print COMMAND_LOG "\tMean coverage threshold:                   \t$coverage_threshold_for_all (for either cases or controls)\n\n\n";
}

print COMMAND_LOG "Thresholds for writing _coverage_assoc.out:\n\n";

print COMMAND_LOG "\tVariation within cases or controls:           \t$covar_threshold\n";
print COMMAND_LOG "\tT-test statistic threshold:                   \t$t_test_threshold\n\n\n";


print COMMAND_LOG "Thresholds for looking for possible SINES\n\n";

print COMMAND_LOG "\tThresholds for lengths:\n\n";

print COMMAND_LOG "\t\tLength of flat region before:              \t$flat_before_threshold\n";
print COMMAND_LOG "\t\tLength of flat region after:               \t$flat_after_threshold\n";
print COMMAND_LOG "\t\tSINE length upper threshold:               \t$sine_length_upper_threshold\n\n\n";
print COMMAND_LOG "\t\tSINE length lower threshold:               \t$sine_length_lower_threshold\n\n";


print COMMAND_LOG "\tThresholds for change in Coverage Ratio:\n\n";

print COMMAND_LOG "\t\tMax coverage ratio change in before/after: \t$crc_flat_threshold\n";
print COMMAND_LOG "\t\tMax coverage ratio change central repeat:  \t$crc_flat_during_threshold\n";
print COMMAND_LOG "\t\tMin coverage ratio change for cliff:       \t$crc_cliff_threshold\n";
print COMMAND_LOG "\t\tMax cliff difference (start/end) threshold:\t$cliff_difference_threshold\n\n";

print COMMAND_LOG "\t\tNOTE: the 'cliff' means the sharp change in coverage at the start and end of the SINE repeat\n\n";


if ($files_to_use eq "bam")
{
	######################
	# DOC ANALYSIS LOOP  #
	######################
		
	print "\n\n";
	print "#-----------------------------------#\n";
	print "# Calculating depth of coverage     #\n";
	print "#-----------------------------------#\n\n";

	##############
	# Affecteds  #
	##############
	for ($file_count=1;$file_count <=$no_of_files;$file_count++)
	{
		$bam_file = $bam_file_array[$file_count];
		$prefix=&get_prefix($bam_file);
		$coverage_output_file = "$run_title"."_"."$prefix"."_coverage.out";
		$coverage_mpileup_file = "$run_title"."_"."$prefix"."_mpileup.out";
		
		if ($use_samtools eq "no")
		{
			$input_string = $input_string." $coverage_output_file";
			&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T DepthOfCoverage $GATK_region_string -R $ref -I $bam_file -o $coverage_output_file -S $GATK_validation_stringency");
		}
		if ($use_samtools eq "yes")
		{
			$input_string = $input_string." $coverage_mpileup_file";
			&run_unix_command ("/opt/samtools/samtools mpileup $bam_file -q 0 -D $samtools_region_string > $coverage_mpileup_file");
		}

		
		################################
		# Delete various metrics files #
		################################
		if ($use_samtools eq "no")
		{
			&delete_file ("$coverage_output_file".".sample_cumulative_coverage_counts");
			&delete_file ("$coverage_output_file".".sample_cumulative_coverage_proportions");
			&delete_file ("$coverage_output_file".".sample_interval_statistics");
			&delete_file ("$coverage_output_file".".sample_interval_summary");
			&delete_file ("$coverage_output_file".".sample_statistics");
			&delete_file ("$coverage_output_file".".sample_summary");
		}
		
	} # end of file_count loop

	$input_string = $input_string." ";

	#################################
	# Combine output files into one #
	# multi-column file             #
	#################################
	print "\n\n";
	print "#----------------------------------------------------------------------------------#\n";
	print "# Using UNIX 'paste' command to combine output files into 'coverage_all_large.out' #\n";
	print "#----------------------------------------------------------------------------------#\n\n";

	&run_unix_command("paste $input_string > $coverage_all_large_out");


}# if $files_to_use eq "bam"


##################################################################################
# If the individual coverage.out files have been made using GATK DepthOfCoverage #
# and the combined coverage_all_large.out file has also been made                #
##################################################################################
if (($files_to_use eq "bam") || ($files_to_use eq "coverage_all_large"))
{

	####################################################################################
	# Simplify coverage_all_large.out file to produce the main "coverage_all.out" file #
	# This just has the positions and the coverage for each BAM file                   #
	# whereas the coverage_all_large.out file has more columns which can be removed    #
	####################################################################################

	print "\n\n";
	print "#----------------------------------------------------------------------------------------#\n";
	print "# Simplifying 'coverage_all_large.out' file to create combined 'coverage_all.out' file   #\n";
	if ($filter_large eq "true")
	{
	print "# (this is also filtered to remove regions of low coverage)                              #\n";
	}
	print "#----------------------------------------------------------------------------------------#\n\n";

	print COMMAND_LOG "\n\nSimplifying $coverage_all_large_out to create $coverage_all_out\n\n";
	
	open (LARGE, "$coverage_all_large_out") || die "Cannot open $coverage_all_large_out";
	open (OUT_ALL, ">$coverage_all_out") || die "Cannot open output_file";

	###########
	# Headers #
	###########
	print OUT_ALL "CHR\tPOSITION";

	
	#########################
	# Print column headings #
	#########################
	for ($file_count=1;$file_count <=$no_of_files;$file_count++)
	{
		print OUT_ALL "\t$bam_file_array[$file_count]";
	}
	print OUT_ALL "\n";

	
	##############################
	# Set coverage array to zero #
	##############################
	for ($file_count=1;$file_count <=$no_of_files;$file_count++)
	{
		$total_coverage_array[$file_count] = 0;
		$mean_coverage_array[$file_count] = 0;
	}
	
	
	####################################
	# Read coverage_all_large.out file #
	####################################
	while ($single_line = <LARGE> ) 
	{
		$line_count = $line_count + 1;
		chomp $single_line;
		
		if ($line_count > 1)  
		{
			@item=split(/\t/,$single_line);
				
			$array_size = scalar @item;
			$position_string = $item[0];
			
			@item_position = split(/:/,$position_string);
			$chromosome=$item_position[0];
			$position=$item_position[1];
				
			if (($line_count % 200000) == 0){print "Simplifying coverage_all_large_out. \tLine: $line_count\tPosition: $position\n";}

			$no_of_blocks = $array_size / $no_output_columns;
			
			################################
			# These two should be the same #
			# This is just a reality check #
			################################
			if ($no_of_blocks == $no_of_files)
			{
				$affected_total_coverage = 0;
				$normal_total_coverage = 0;
				
				$coverage_string = "";
				
				
				for ($block_count=1;$block_count <=$no_of_blocks;$block_count++)
				{
					$coverage_array[$block_count] = $item[($block_count * $no_output_columns)-2];
					$coverage_array[$block_count] = sprintf('%.0f', $coverage_array[$block_count]); 
					
					if ($status_array[$block_count] eq "affected")
					{
						$affected_total_coverage = $affected_total_coverage + $coverage_array[$block_count];
					}
					if ($status_array[$block_count] eq "control")
					{
						$normal_total_coverage = $normal_total_coverage + $coverage_array[$block_count];
					}
					
					#Store coverage columns until you see if it is worth writing this line
					$coverage_string = "$coverage_string"."\t"."$coverage_array[$block_count]";
					
					######################################
					# Store total coverage for each file #
					# to get the mean coverage later     #
					######################################
					$total_coverage_array[$block_count] = $total_coverage_array[$block_count] + $coverage_array[$block_count];
					
				}# $block_count loop
				
				###################################
				# Write to coverage_all.out file  #
				###################################
				print OUT_ALL "$chromosome\t$position"; 
				print OUT_ALL "\t$coverage_string\n";
						
			} # if ($no_of_blocks == $no_of_files)
			
			
			###########################################################
			# If no of blocks isn't correct then record as a bad line #
			###########################################################
			if ($no_of_blocks != $no_of_files)
			{
				$bad_line_count = $bad_line_count + 1;
			} # if ($no_of_blocks != $no_of_files)
			
		} # end of if $line_count > 1

			#if ($line_count > 1000000){goto HERE1};
			
	} # end of while $single_line
	HERE1:
	close LARGE;
	close OUT_ALL;

	print COMMAND_LOG "\n\nFile $coverage_all_out created from $coverage_all_large_out\n\n";
	
	$no_of_lines_large_out = $line_count;
	
	
	####################################################
	# Calculate mean coverage for each file            #
	# (This uses data from reading coverage_all_large) #
	####################################################
	print "\nMean coverages in $coverage_all_large_out: \n\n";
	print COMMAND_LOG "\nMean coverages in $coverage_all_large_out: \n\n";
	
	$sum_of_means = 0;
	
	for ($file_count=1;$file_count <=$no_of_files;$file_count++)
	{
		if ($no_of_lines_large_out > 0 )
		{
			$mean_coverage_array[$file_count] = $total_coverage_array[$file_count] / $no_of_lines_large_out;
			$sum_of_means = $sum_of_means + $mean_coverage_array[$file_count];
			
		}
	}
	$mean_of_means = $sum_of_means / $no_of_files;
	
	
	##########################################################
	# Calculate weighting factors                            #
	# (To compensate for the fact that some samples may have #
	# overall better coverage than others for some reason    #
	##########################################################
	print             "\nMean coverages and weighting factors from $coverage_all_large_out: \n\n";
	print COMMAND_LOG "\nMean coverages and weighting factors from $coverage_all_large_out: \n\n";
	
	print             "File No \tName\t\tMean coverage\tWeighting\n\n";
	print COMMAND_LOG "File No \tName\t\tMean coverage\tWeighting\n\n";
	
	for ($file_count=1;$file_count <=$no_of_files;$file_count++)
	{
		$weighting_array[$file_count] = $mean_coverage_array[$file_count] / $mean_of_means;
		
		# Format numbers #
		$mean_coverage_array[$file_count] = sprintf('%.0f', $mean_coverage_array[$file_count]); 
		$weighting_array[$file_count] = sprintf('%.2f', $weighting_array[$file_count]);
		
		print "File $file_count: \t$bam_file_array[$file_count]\t$mean_coverage_array[$file_count]\t$weighting_array[$file_count]\n";
		print COMMAND_LOG "File $file_count: \t$bam_file_array[$file_count]\t$mean_coverage_array[$file_count]\t$weighting_array[$file_count]\n";
	}
			
} # if $files_to_use eq "bam" OR "coverage_all_large"



########################################################################################
# This part of the analysis can be used on previously created "coverage_all.out" files #
# $files_to_use eq "coverage_all" (this is the filtered version of coverage_all_large) #                                                         #
########################################################################################

if (($files_to_use eq "bam") || ($files_to_use eq "coverage_all_large") || ($files_to_use eq "coverage_all"))
{

	print "\n\n";
	print "#-----------------------------------#\n";
	print "# Reading coverage_all.out file...  #\n";
	print "#                                   #\n";
	print "# Calculating coverage statistics   #\n";
	print "# and creating coverage_assoc.out   #\n";
	print "#-----------------------------------#\n\n";

	#########################################
	# Read all of coverage file into memory #
	#########################################
	print "Reading $coverage_all_out ... \n\n";
	print COMMAND_LOG "\nReading $coverage_all_out ... \n\n";

	open (IN_ALL, "$coverage_all_out") || die "Cannot open $coverage_all_out\n\n";
		@file_array = <IN_ALL>;
		$no_of_lines_coverage_all_out = scalar @file_array;
	close IN_ALL;

	#print "File $coverage_all_out has been read into the array file_array\n\n";

	#######################################################
	# Open file for filtered coverage data and statistics #
	#######################################################
	open (OUT_ASSOC, ">$output_assoc_file") || die "Cannot open $output_assoc_file\n\n";

	###################################################
	# Open file to save data for later SINE searching #
	###################################################
	open (SINE_DATA, ">$output_sine_data_file") || die "Cannot open $output_sine_data_file\n\n";

	###########
	# Headers #
	###########
	print OUT_ASSOC "CHR\tPOSITION\tTYPE";
	print SINE_DATA "CHR\tPOSITION\tAV_A\tAV_N\tCOVAR\tRATIO\tCHANGE\n";

	########################################################
	# Print column headings, affected first, then controls #
	# if $include_individual_columns is "false" then these #
	# aren't written (makes a smaller output file)         #
	########################################################

	if ($include_individual_columns eq "true")
	{
		for ($file_count=1;$file_count <=$no_of_files;$file_count++)
		{
			print OUT_ASSOC "\t$bam_file_array[$file_count]";
		}
	} # if you want indivial DP data columns


	print OUT_ASSOC "\tAV_A\tAV_N\tCOV_A\tCOV_N\tRATIO\tCHANGE\tT_TEST\n";

	$line_count =0;

	print "No of lines: $no_of_lines_coverage_all_out in $coverage_all_out\n\n";
	print COMMAND_LOG "No of lines: $no_of_lines_coverage_all_out in $coverage_all_out\n\n";
	
	
	###################################################
	# Looping through file_array for coverage_all.out #
	# to create ASSOC out file and SINE out file      #
	###################################################
	for ($line_count = 0;$line_count<=$no_of_lines_coverage_all_out - 1; $line_count++)
	{
		if ($line_count > 0)
		{
			$single_line = $file_array[$line_count];
			@item=split(/\t/,$single_line);
			$array_size = scalar @item;

			$chromosome=$item[0];
			$position=$item[1];
			
			if (($line_count % 200000) == 0){print "Creating filtered assoc.out file.  Line: $line_count  \tPosition: $position\n";}

			$no_of_blocks = $array_size - 3;
			
			if ($no_of_blocks == $no_of_files)
			{
				$affected_total_coverage = 0;
				$normal_total_coverage = 0;
				$affected_total_coverage_normalised = 0;
				$normal_total_coverage_normalised = 0;
				
				$coverage_string = "";
				
				#################################################
				# Put together string of coverages to show user #
				#################################################
				for ($block_count=1;$block_count <=$no_of_blocks;$block_count++)
				{
					$coverage_array[$block_count] = $item[$block_count + 2];
					$coverage_array[$block_count] = sprintf('%.0f', $coverage_array[$block_count]); 

					$coverage_string = "$coverage_string"."\t"."$coverage_array[$block_count]";
					
					if ($status_array[$block_count] eq "affected")
					{
						#print "Block count: $block_count\n";
						#print "Status: $status_array[$block_count]\n";
						#print "coverage_array[$block_count]: \t$coverage_array[$block_count]\n";
						#print "weighting_array[$block_count]: \t$weighting_array[$block_count]\n";
						
						$affected_total_coverage = $affected_total_coverage + $coverage_array[$block_count];
						$affected_total_coverage_normalised = $affected_total_coverage_normalised + ($coverage_array[$block_count] * $weighting_array[$block_count]);
					}
					
					if ($status_array[$block_count] eq "control")
					{
						#print "Block count: $block_count\n";
						#print "Status: $status_array[$block_count]\n";
						#print "coverage_array[$block_count]: \t$coverage_array[$block_count]\n";
						#print "weighting_array[$block_count]: \t$weighting_array[$block_count]\n";
						
						$normal_total_coverage = $normal_total_coverage + $coverage_array[$block_count];
						$normal_total_coverage_normalised = $normal_total_coverage_normalised + ($coverage_array[$block_count] * $weighting_array[$block_count]);
					}
				}

				###################################################
				# Mean coverages on each line of coverage_all.out #
				###################################################
				if ($normalise_coverage eq "false")
				{
					$mean_coverage_affected = $affected_total_coverage / $no_of_affecteds;
					$mean_coverage_normal = $normal_total_coverage / $no_of_normals;
					$mean_coverage_both = ($affected_total_coverage + $normal_total_coverage) / ($no_of_affecteds + $no_of_normals);
				
				}
				
				if ($normalise_coverage eq "true")
				{
					$mean_coverage_affected = $affected_total_coverage_normalised / $no_of_affecteds;
					$mean_coverage_normal = $normal_total_coverage_normalised / $no_of_normals;
					$mean_coverage_both = ($affected_total_coverage_normalised + $normal_total_coverage_normalised) / ($no_of_affecteds + $no_of_normals);
				
				}

				
				################################################################
				# Normalised coverages (using mean for each sample to correct) #
				################################################################
					
				############
				# Variance #
				############
				$sum_squares_affected = 0;
				$sum_squares_normal = 0;
				$variance_affected = 0;
				$variance_normal = 0;
				$covar_affected = 0;
				$covar_normal = 0;
				$coverage_ratio = 0;

				for ($block_count=1;$block_count <=$no_of_blocks;$block_count++)
				{
					if ($status_array[$block_count] eq "affected")
					{
						$sum_squares_affected = $sum_squares_affected + (($coverage_array[$block_count] - $mean_coverage_affected) ** 2);
					}
					if ($status_array[$block_count] eq "control")
					{
						$sum_squares_normal = $sum_squares_normal + (($coverage_array[$block_count] - $mean_coverage_normal) ** 2);
					}
				}
				
				$variance_affected = $sum_squares_affected / $no_of_affecteds;
				$variance_normal = $sum_squares_normal / $no_of_normals;

				
				######################
				# Standard Deviation #
				######################
				$sd_affected = sqrt ($variance_affected);
				$sd_normal = sqrt ($variance_normal);
			
			
				##########
				# T-test #
				##########
				$t_test = 0;
				if (($variance_affected > 0) && ($variance_normal > 0))
				{
					$t_test = ($mean_coverage_affected - $mean_coverage_normal) / (sqrt((($variance_affected / $no_of_affecteds) + ($variance_normal / $no_of_normals))));
					$t_test = abs($t_test); # t_statistic can be positve or negative
				}
				
				############################
				# Coefficient of variation #
				############################
				if ($mean_coverage_affected > 0 ) {$covar_affected = $sd_affected / $mean_coverage_affected;}
				if ($mean_coverage_normal > 0) {$covar_normal = $sd_normal / $mean_coverage_normal;}

				$covar_all_samples = ($covar_affected + $covar_normal ) / 2;
				
				##################
				# Coverage ratio #
				##################
				if ($mean_coverage_normal > 0) {$coverage_ratio = $mean_coverage_affected / $mean_coverage_normal;}


				################################
				# Format to two decimal places #
				################################ 
				$covar_affected = sprintf('%.2f', $covar_affected); 
				$covar_normal = sprintf('%.2f', $covar_normal); 
				$coverage_ratio = sprintf('%.2f', $coverage_ratio );
				$t_test = sprintf('%.2f', $t_test );
				$covar_all_samples = sprintf('%.2f', $covar_all_samples );
				
			} # if no_blocks = no_files
			
			if ($no_of_blocks != $no_of_files)
			{
				$bad_line_count = $bad_line_count + 1;
			}
			
			###################################
			# Check change in coverage ratios #
			###################################

			if ($coverage_ratio > 0)
			{
				$coverage_ratio_change_signed = (1 - ($last_coverage_ratio / $coverage_ratio));
				$coverage_ratio_change = abs($coverage_ratio_change_signed);	
			}

			
			################################
			# Format to two decimal places #
			################################
			$mean_coverage_affected = sprintf('%.2f', $mean_coverage_affected);
			$mean_coverage_normal = sprintf('%.2f', $mean_coverage_normal);
			$coverage_ratio_change = sprintf('%.2f', $coverage_ratio_change);
			$coverage_ratio_change_signed = sprintf('%.2f', $coverage_ratio_change_signed);
			

			###################################################
			# Make assoc.out file with only filtered lines    #
			###################################################
			$keep = "false";
			$keep_t_test = "true";
			$keep_covar = "true";
			$keep_coverage = "true";
			$keep_crc = "true";		

			# If coefficient of variance is too great then the line is filtered out #
			if (($covar_affected > $covar_threshold) || ($covar_normal > $covar_threshold)){$covar_filtered_count = $covar_filtered_count + 1}

			# If t_test difference is too small then the line is filtered out from the assoc file #
			if ($t_test <= $t_test_threshold){$t_test_filtered_count = $t_test_filtered_count + 1}

			# If mean_coverage difference is too small then the line is filtered out #
			if ($mean_coverage_both <= $coverage_threshold_for_assoc){$coverage_filtered_count = $coverage_filtered_count + 1}
			
			# if filter eq "true" means if the global variable $filter_assoc is "false" then output everything #
			if ($filter_assoc eq "true")
			{
				if ($t_test <= $t_test_threshold)
				{
					$keep_t_test = "false";
				}

				if (($covar_affected > $covar_threshold) || ($covar_normal > $covar_threshold)) 
				{
					$keep_covar = "false";
				}
				
				if (($mean_coverage_both <= $coverage_threshold_for_assoc) || ($mean_coverage_normal <= $coverage_threshold_for_assoc)) 
				{
					$keep_coverage = "false";
				}

			} # filtering on

		
			#######################################################
			# print to SINE_DATA file for SINE data               #
			# (no filtering here as we want consecutive positions #
			#######################################################
			print SINE_DATA "$chromosome\t$position\t$mean_coverage_affected\t$mean_coverage_normal\t$covar_all_samples\t$coverage_ratio\t$coverage_ratio_change_signed\n";

			
			#############################################################################
			# Only write to the assoc.out file if the t-test difference is great enough  #
			# and if the variation is not too great for either cases or controls        #
			#############################################################################
			if (($keep_t_test eq "true") && ($keep_covar eq "true") && ($keep_coverage eq "true"))
			{
				print OUT_ASSOC "$chromosome\t$position\tCOVERAGE";
				if ($include_individual_columns eq "true"){print OUT_ASSOC "$coverage_string";}
				print OUT_ASSOC "\t$mean_coverage_affected\t$mean_coverage_normal";
				print OUT_ASSOC "\t$covar_affected\t$covar_normal";
				print OUT_ASSOC "\t$coverage_ratio\t$coverage_ratio_change"; # CHANGE column is used in coverage
				print OUT_ASSOC "\t$t_test";
				print OUT_ASSOC "\n";
				
				$keep_for_assoc_count = $keep_for_assoc_count + 1;
			}

	#if ($keep_for_assoc_count > 100000) {close OUT_ASSOC; close SINE_DATA; print "\n\n\nTEMP EXIT ON 1000000\n\n\n"; goto HERE;}


			#if ($filter_assoc eq "false")
			#{
			#	print OUT_ASSOC "$chromosome\t$position";
			#	print OUT_ASSOC "$coverage_string";
			#	print OUT_ASSOC "\t$mean_coverage_affected\t$mean_coverage_normal\t$coverage_ratio\t$coverage_ratio_change";
			#	print OUT_ASSOC "\t$covar_affected\t$covar_normal\t$t_test";
			#	print OUT_ASSOC "\n";
			#} # filtering off

			###############################################
			# Store the last coverage ratio that was used #
			###############################################
			$last_coverage_ratio = $coverage_ratio;
				
		} # end of if $line_count > 0
		
	} # $line_count LOOP reading coverage_all_large.out file

	close OUT_ASSOC;
	close SINE_DATA;

	HERE:

} # if files_to_use eq coverage bam or LARGE (NOT sine)

print "\nBad line count in coverage.out file: $bad_line_count\n\n";
print COMMAND_LOG "\nBad line count in coverage.out file: $bad_line_count\n\n";

##################
# LOOK FOR SINES #
##################


##########################
# Read in SINE data file #
##########################

print "Reading in sine data file $output_sine_data_file...\n\n";
print COMMAND_LOG "Reading in sine data file $output_sine_data_file...\n\n";

open (IN, "$output_sine_data_file") || die "Cannot open $output_sine_data_file";
	@file_array = <IN>;
	$no_of_lines_sine_data = scalar @file_array;
close IN;

print "No of lines in $output_sine_data_file: $no_of_lines_sine_data\n\n\n";
print COMMAND_LOG "No of lines in $output_sine_data_file: $no_of_lines_sine_data\n\n\n";



###################################################
# Look for SINES in sine_data file  (SINE LOOP)   #
###################################################
$window_size = $flat_before_threshold + 1 + $sine_length_upper_threshold + $flat_after_threshold;
$coverage_ratio_change = 0;
$cliff_position = 1;


print "################################\n";
print "#  Checking for SINE elements  #\n";
print "################################\n\n";

print COMMAND_LOG "Checking for SINE elements...\n\n";

open (SINES, ">$possible_sines_out") || die "Cannot open $possible_sines_out";
open (SINES_PLOT, ">$possible_sines_plot_data_out") || die "Cannot open $possible_sines_plot_data_out";

print SINES "===================================================\n";
print SINES "Possible SINES output file from Coverage version $version\n";
print SINES "===================================================\n\n";

print SINES "Thresholds for statistics:\n\n";

print SINES "\tVariation within cases or controls threshold:\t$covar_threshold\n";
print SINES "\tT-test statistic threshold:                  \t $t_test_threshold\n\n\n";


print SINES "Thresholds for lengths:\n\n";

print SINES "\tLength of flat region before:              \t$flat_before_threshold\n";
print SINES "\tLength of flat region after:               \t$flat_after_threshold\n";
print SINES "\tSINE length lower threshold:               \t$sine_length_lower_threshold\n";
print SINES "\tSINE length upper threshold:               \t$sine_length_upper_threshold\n\n\n";

print SINES "Thresholds for change in Coverage Ratio:\n\n";

print SINES "\tCoverage ratio change in before/after:   \t$crc_flat_threshold\n";
print SINES "\tCoverage ratio change in central repeat: \t$crc_flat_during_threshold\n";
print SINES "\tCoverage ratio change for cliff:         \t$crc_cliff_threshold\n";
print SINES "\tCliff difference (start/end) threshold:  \t$cliff_difference_threshold\n\n";

if ($debugging_mode eq "false"){print SINES "CHR\tPOSITION\tSINE_SCORE\tSELECTION\n";}

if ($debugging_mode eq "true")
{
	print SINES "CHR\tPOSITION\tSINE_SCORE\tSELECTION";

	print SINES "\tMCB\tMCD\tMCA\tCD\tCS1\tMCOV\tMCOV_d\tSRL\tMAFF\tMNO\tFBC\tFAC\tC1\tC2\n";
}	

#############################################################
# This score is used to select regions that score well for  #
# the things that make a SINE: flat coverage before etc etc #
# This is pretty strict but I have subtracted 1 at the      #
# end to allow a few more things through                    #
#############################################################
$inclusion_threshold = $flat_before_threshold + 1 + 1 + $flat_after_threshold;
		
$checked_count = 0;
		
for ($line_count = 0;$line_count<=$no_of_lines_sine_data +1; $line_count++)
{
	if (($line_count > 1) && ($line_count < $no_of_lines_sine_data - $window_size))
	{
		if (($line_count % 200000) == 0)
		{
			$percentage_done = ($line_count / $no_of_lines_sine_data) * 100;
			$percentage_done = sprintf('%.1f', $percentage_done);
			
			print "Checking for SINEs.  Percent done: $percentage_done %\tPosition: $cliff_position\n";
		}
		
		##########################################
		# Start looking forward from line_count  #
		# (moving window)                        #
		##########################################
		
		############################################################
		# Make string to show when you find SINE                   #
		#                                                          #
		# Also calculate Coefficient of Variation for whole window #
		# --> gets mean_covar_whole_window	                       #
		############################################################
		$coverage_string = "";
		$covar_whole_window = 0;
		
		for ($check_count = 1;$check_count <= $window_size; $check_count++)
		{
			$single_line = $file_array[$line_count + $check_count - 1];
			chomp $single_line;
			@item=split(/\t/,$single_line);
			
			$covar_whole_window = $covar_whole_window + $item[4];
			$coverage_ratio_change_signed = $item[6];
			$coverage_string = "$coverage_string "."$coverage_ratio_change_signed";
		}
		
		$mean_covar_whole_window = $covar_whole_window / $window_size;
		$mean_covar_whole_window = sprintf('%.2f', $mean_covar_whole_window);
		
		
		##############################################################
		# Get data at the position where you are looking for a cliff #
		##############################################################
		$single_line = $file_array[$line_count + $flat_before_threshold];
		chomp $single_line;
		@item=split(/\t/,$single_line);
		$chromosome = $item[0];
		$cliff_position = $item[1];
	
	
		#########################
		# Set variables to zero #
		#########################
		$sine_score = 0;
		
		$cliff_1_count = 0;
		$cliff_2_count = 0;
		$cliff_1_size = 0;
		$cliff_2_size = 0;
		$cliff_1_size_abs = 0;
		$cliff_2_size_abs = 0;
		$cliff_difference = 0;
		
		$flat_before_count = 0;
		$flat_after_count =  0;
		$flat_during_count = 0;
		$sine_repeat_length =0;
		
		$total_change_before = 0;
		$total_change_during = 0;
		$total_change_after = 0;
		
		$mean_change_before = 0;
		$mean_change_during = 0;
		$mean_change_after = 0;
		
		$covar_all_samples = 0;
		$total_covar = 0;
		
		
		#################################################################
		# Check positions BEFORE SINE (actually before the first cliff) #
		# --> gets flat_before_count                                    #
		#################################################################
		for ($check_count = 1;$check_count <= $flat_before_threshold; $check_count++)
		{
			$single_line = $file_array[$line_count + $check_count - 1];
			chomp $single_line;
			
			@item=split(/\t/,$single_line);

			$coverage_ratio = $item[5];
			$coverage_ratio_change_signed = $item[6];
			$coverage_ratio_change = abs ($coverage_ratio_change_signed);
			
			$total_change_before = $total_change_before + abs($coverage_ratio_change);
			
			if ($coverage_ratio_change <=  $crc_flat_threshold){$flat_before_count = $flat_before_count + 1}
			
			if ($show_during eq "true")
			{#TEMPORARY!!!!
				print "Before Cliff. Line: $line_count\tCheck: $check_count\n";
				print "$coverage_string\n";
				print "CRC: $coverage_ratio_change_signed\tTotal change before: $total_change_before\tFlat before count: $flat_before_count\n\n";
				$answer=<STDIN>;
			}
			
		} # before cliff
		
		
		#############################################################
		# Check actual "cliff" i.e. where coverage suddenly changes #
		# --> gets cliff_1_size, cliff_1_count, mean_coverage_affected,
		#     mean_coverage_normal
		#############################################################
		$single_line = $file_array[$line_count + $flat_before_threshold ];
		chomp $single_line;
		@item=split(/\t/,$single_line);
		
		$mean_coverage_affected = $item[2];
		$mean_coverage_normal = $item[3];
		$coverage_ratio = $item[5];
		$coverage_ratio_change_signed = $item[6];
		$coverage_ratio_change = abs ($coverage_ratio_change_signed);
		
		$covar_all_samples = $item[4];
		
		if ($coverage_ratio_change >  $crc_cliff_threshold)
		{
			$cliff_1_count = $cliff_1_count + 1;
			$cliff_1_size = $coverage_ratio_change_signed;
		}
		
		
		#####################################################################
		# Check positions DURING SINE (between the first and second cliffs) #
		# --> gets sine_repeat_length
		#####################################################################
		$covar_during = 0;
		for ($check_count = $flat_before_threshold + 2;$check_count < $flat_before_threshold  + 2 + $sine_length_upper_threshold; $check_count++)
		{
			$sine_repeat_length = $sine_repeat_length + 1;
			$single_line = $file_array[$line_count + $check_count - 1];
			chomp $single_line;
			@item=split(/\t/,$single_line);
			
			$covar_all_samples = $item[4]; # <-- covar across all samples at this position
			
			$covar_during = $covar_during + $covar_all_samples;
			
			$coverage_ratio = $item[5];
			$coverage_ratio_change_signed = $item[6];
			$coverage_ratio_change = abs ($coverage_ratio_change_signed);

			############################################################
			# Keep total change score (but only up until second cliff) #
			# --> gets flat_during_count
			############################################################
			if ($coverage_ratio_change <= $crc_cliff_threshold)
			{
				$total_change_during = $total_change_during + abs($coverage_ratio_change);
			}
			
			if ($coverage_ratio_change <=  $crc_flat_during_threshold){$flat_during_count = $flat_during_count + 1}
		
			if ($show_during eq "true")
			{#TEMPORARY!!!!
				print "During SINE. Line: $line_count\tCheck: $check_count\n";
				print "$coverage_string\n";
				print "SINE repeat length: $sine_repeat_length\tCRC: $coverage_ratio_change_signed\tTotal change during: $total_change_during\n\n";
				$answer=<STDIN>;
			}
			
			
			###############################################################
			# End SINE if there is another cliff
			# (at the moment this is either positive or negative. Change? #
			# Could look at cliff_difference?
			# --> gets cliff_2_count
			###############################################################
			$end_of_sine = $check_count;
			if ($coverage_ratio_change > $crc_cliff_threshold)
			{
				$end_of_sine = $check_count;
				$cliff_2_count = $cliff_2_count + 1;
				$cliff_2_size = $coverage_ratio_change_signed;
				goto END_OF_SINE;
			}
		} # during SINE
		
		END_OF_SINE:
		
		################################
		# Calculate mean_covar_during  #
		################################
		$mean_covar_during = $covar_during/$sine_repeat_length;
		$mean_covar_during = sprintf('%.2f', $mean_covar_during);
		
		
		###################################################
		# Check positions AFTER SINE (after second cliff) #
		# --> gets flat_after_count
		###################################################
		for ($check_count = $end_of_sine + 1;$check_count <= ($end_of_sine + $flat_after_threshold); $check_count++)
		{
			$single_line = $file_array[$line_count + $check_count - 1];
			chomp $single_line;
			@item=split(/\t/,$single_line);

			$coverage_ratio = $item[5];
			$coverage_ratio_change_signed = $item[6];
			$coverage_ratio_change = abs ($coverage_ratio_change_signed);

			$total_change_after = $total_change_after + abs($coverage_ratio_change);
				
			if ($coverage_ratio_change <=  $crc_flat_threshold){$flat_after_count = $flat_after_count + 1}
			
			if ($show_during eq "true")
			{#TEMPORARY!!!!
				print "After Cliff. Line: $line_count\tCheck: $check_count\n";
				print "$coverage_string\n";
				print "CRC: $coverage_ratio_change_signed\tTotal change after: $total_change_after\tFlat after count: $flat_after_count\n\n";
				$answer=<STDIN>;
			}
			
		} # after cliff
		
		
		##########################################################
		# If there are too many flat afters, reduce to threshold #
		##########################################################
		if ($flat_after_count > $flat_after_threshold){$flat_after_count = $flat_after_threshold}
		
		
		#############################
		# --> gets cliff_difference #
		#############################
		$cliff_difference = abs($cliff_1_size + $cliff_2_size);
		$cliff_1_size_abs = abs($cliff_1_size);
		$cliff_2_size_abs = abs($cliff_2_size);
		
		
		###################################################
		# Calculate mean changes before, during and after #
		###################################################
		if ($flat_before_threshold > 0)
		{
			$mean_change_before = ($total_change_before/$flat_before_threshold);
		}
		if ($sine_repeat_length > 0)
		{
			$mean_change_during = ($total_change_during/$sine_repeat_length);
		}
		if ($flat_after_threshold > 0)
		{
			$mean_change_after = ($total_change_after/$flat_after_threshold);
		}
		
		$mean_change_before = sprintf('%.2f', $mean_change_before);
		$mean_change_during = sprintf('%.2f', $mean_change_during);
		$mean_change_after = sprintf('%.2f', $mean_change_after);
		
		
		#######################################################################
		# Calculate SINE SCORE - This is another way of selecting interesting #
		# regions and uses a quality score rather then a strict threshold     #
		#######################################################################
		
		$sine_score = 10 - $mean_change_before - $mean_change_during - 
		$mean_change_after - $cliff_difference + ($cliff_1_size_abs/2) - ($mean_covar_during * 2);
		
		
		#######################################################
		# If there is no cliff we don't want to select anyway #
		# so set the sine score to 1                          #
		#######################################################
		if ($cliff_1_size_abs == 0){$sine_score = 1;}
		
		$sine_score = sprintf('%.2f', $sine_score);
		
		
		########################################################################
		# Calculate INCLUSION SCORE - This is used to select for regions       #
		# that have a good score for each of the various requirements          #
		# such as flat before etc.  Saved if greater than inclusion_threshold. #   
		# Note: flat_during_count is not used as this can vary 	               #
		########################################################################
		$inclusion_score = $flat_before_count + $cliff_1_count + $cliff_2_count + $flat_after_count;
		
		
		####################################################
		# Save if SINE score is above SINE_score_threshold #
		####################################################
		$saved_by_sine_score = "false";
		
		if ($sine_score >= $sine_score_threshold)
		{
			$saved_by_sine_score = "true";
			$saved_by_sine_score_count = $saved_by_sine_score_count + 1;
			
			if ($verbose_out eq "true")
			{
				print "  >>>>>>  Selected on SINE score of $sine_score\n";
				
				print SINES "Selected on SINE score\n\n";
				
				print SINES "CHR: $chromosome\tStart of cliff:  $cliff_position\tSINE score: $sine_score\n\n";
									
				print SINES "Line in file: $line_count\n";
				print SINES "Affected mean: $mean_coverage_affected\t\tNormal mean: $mean_coverage_normal\n\n";
				print SINES "Changes in coverage ratio in the region:\n\n";
				print SINES "$coverage_string\n\n";
				print SINES "Change thresholds.   Before and after cliffs: $crc_flat_threshold    Between cliffs: $crc_flat_during_threshold   Size of cliff: $crc_cliff_threshold\n";
				print SINES "\tFlat before count:     \t$flat_before_count (out of $flat_before_threshold)\n";
				print SINES "\tCliff 1 count:         \t$cliff_1_count\n";
				print SINES "\tFlat during count:     \t$flat_during_count (Range is $sine_length_lower_threshold - $sine_length_upper_threshold)\n";
				print SINES "\tCliff 2 count:         \t$cliff_1_count\n";
				print SINES "\tFlat after count:      \t$flat_after_count (out of $flat_after_threshold)\n";
				print SINES "\tCliff difference:      \t$cliff_difference (threshold is $cliff_difference_threshold)\n";
				print SINES "\tSine repeat length:    \t$sine_repeat_length\n\n";
				print SINES "\tCoeff of Var:          \t$covar_whole_window\n\n";
				print SINES "Inclusion criteria:\n\n";
				print SINES "Counts:    \t abs($flat_before_count + $cliff_1_count + $flat_during_count + $cliff_2_count + $flat_after_count ))\n";
				print SINES "Threshold: \t abs($flat_before_threshold  + 1 + $sine_length_lower_threshold + 1 + $flat_after_threshold)\n\n";

				print SINES "===================================================================\n\n";
			}# if verbose
			
			# SINE SCORE
			if ($debugging_mode eq "false")
			{print SINES "$chromosome\t$cliff_position\t$sine_score\tSINE_SCORE\n";}	
			 
			 if ($debugging_mode eq "true")
			{
				print SINES "$chromosome\t$cliff_position\t$sine_score\tSINE_SCORE";
				
				print SINES "\t$mean_change_before\t$mean_change_during\t$mean_change_after\t$cliff_difference";
				
				print SINES "\t$cliff_1_size_abs\t$mean_covar_whole_window\t$mean_covar_during\t$sine_repeat_length";
					
				print SINES "\t$mean_coverage_affected\t$mean_coverage_normal\t$flat_before_count\t$flat_after_count\t$cliff_1_count\t$cliff_2_count\n";
					
			}
			
			
		}# sine score threshold
		
		###############################################################
		# Save if INCLUSION_SCORE is greater then inclusion_threshold #
		# OMIT THIS METHOD OF SELECTION FOR THE MOMENT                #
		###############################################################
		$saved_by_inclusion = "false";
		
		#if (($mean_affected > $mean_sine_coverage_threshold) && ($mean_normal > $mean_sine_coverage_threshold))
		#{
		#	if ($inclusion_score >= ($inclusion_threshold ))
		#	{
		#		$saved_by_inclusion = "true";
		#		$saved_by_inclusion_count = $saved_by_inclusion_count + 1;
		#		
		#		# INCLUSION #
		#		if ($debugging_mode eq "false")
		#		{print SINES "$chromosome\t$cliff_position\t$sine_score\tINCLUSION\n";}	
		#		 
		#		 if ($debugging_mode eq "true")
		#		{
		#			print SINES "$chromosome\t$cliff_position\t$sine_score\tINCLUSION";
		#			
		#			print SINES "\t$total_change_before\t$total_change_during\t$total_change_after\t$cliff_difference";
		#			
		#			print SINES "\t$cliff_1_size_abs\t$mean_covar_whole_window\t$mean_covar_during\t$sine_repeat_length";
		#			
		#			print SINES "\t$mean_coverage_affected\t$mean_coverage_normal\t$flat_before_count\t$flat_after_count\t$cliff_1_count\t$cliff_2_count\n";
		#			
		#		}
		#		
		#	} # inclusion threshold
		#} # above coverage thresholds
		
		########################################################################
		# STRICT selection
		#
		# This section only saves the region if it strictly passes the 
		# criterion for each characteristic.
		# Some of these are SINE related, e.g. having two cliffs, and some of
		# them are general quality measures such as having a reasonable
		# depth of coverage (greater than mean_sine_threshold)
		#
		# NOTE: currently doesn't use any flat during measure
		########################################################################
		$saved_by_strict = "false";
		
		
		######################################################
		# Count how many are filtered when looking for SINEs #
		######################################################
		
		if (($mean_coverage_affected <= $mean_sine_coverage_threshold) || ($mean_coverage_normal <= $mean_sine_coverage_threshold))
		{$sine_coverage_filtered_count = $sine_coverage_filtered_count + 1}
		
		if (($cliff_1_count == 0) || ($cliff_2_count == 0))
		{$sine_cliff_filtered_count = $sine_cliff_filtered_count + 1}
		
		if ($cliff_difference > $cliff_difference_threshold)
		{$sine_cliffdiff_filtered_count = $sine_cliffdiff_filtered_count + 1}
		
		if (($sine_repeat_length < $sine_length_lower_threshold) || ($sine_repeat_length > $sine_length_upper_threshold))
		{$sine_repeat_length_filtered_count = $sine_repeat_length_filtered_count + 1}
		
		if ($mean_covar_whole_window > $covar_whole_window_threshold)
		{$sine_covar_filtered_count = $sine_covar_filtered_count + 1}
		
		if (($flat_before_count < $flat_before_threshold) || ($flat_after_count < $flat_after_threshold)) 
		{$sine_flat_filtered_count = $sine_flat_filtered_count + 1}
		
		if ($flat_during_count < ($sine_repeat_length - 1))
		{$sine_flat_during_filtered_count = $sine_flat_during_filtered_count + 1}
								
		if (($mean_coverage_affected > $mean_sine_coverage_threshold) && ($mean_coverage_normal > $mean_sine_coverage_threshold))
		{
			if (($cliff_1_count > 0) && ($cliff_2_count > 0)) # MUST have both cliffs
			{
				if ($cliff_difference <= $cliff_difference_threshold)  # cliffs must be of roughly the same size
				{
					if (($sine_repeat_length >= $sine_length_lower_threshold) && ($sine_repeat_length <= $sine_length_upper_threshold))
					{
						if ($mean_covar_whole_window <= $covar_whole_window_threshold) # not too much variation across samples
						{
							if (($flat_before_count >= $flat_before_threshold) && ($flat_after_count >= $flat_after_threshold)) # flat before and after
							{
								if ($flat_during_count >= ($sine_repeat_length - 1))
								{
								
									$saved_by_strict = "true";
									$saved_by_strict_count = $saved_by_strict_count + 1;
									
									if ($verbose_out eq "true")
									{
										print SINES "CHR: $chromosome\tStart of cliff:  $cliff_position\tSINE score: $sine_score\n\n";
										
										print SINES "Line in file: $line_count\n";
										print SINES "Affected mean: $mean_coverage_affected\t\tNormal mean: $mean_coverage_normal\n\n";
										print SINES "Changes in coverage ratio in the region:\n\n";
										print SINES "$coverage_string\n\n";
										print SINES "Change thresholds.   Before and after cliffs: $crc_flat_threshold    Between cliffs: $crc_flat_during_threshold   Size of cliff: $crc_cliff_threshold\n";
										print SINES "\tFlat before count:     \t$flat_before_count (out of $flat_before_threshold)\n";
										print SINES "\tCliff 1 count:         \t$cliff_1_count\n";
										print SINES "\tFlat during count:     \t$flat_during_count (Range is $sine_length_lower_threshold - $sine_length_upper_threshold)\n";
										print SINES "\tCliff 2 count:         \t$cliff_1_count\n";
										print SINES "\tFlat after count:      \t$flat_after_count (out of $flat_after_threshold)\n";
										print SINES "\tCliff difference:      \t$cliff_difference (threshold is $cliff_difference_threshold)\n";
										print SINES "\tSine repeat length:    \t$sine_repeat_length\n\n";
										
										print SINES "\tEnd of sine:           \t$end_of_sine\n";
										print SINES "\tWindow size:           \t$window_size\n";
										print SINES "\tCoeff of Var:          \t$covar_whole_window\n\n";
										print SINES "Inclusion criteria:\n\n";
										print SINES "Counts:    \t abs($flat_before_count + $cliff_1_count + $flat_during_count + $cliff_2_count + $flat_after_count ))\n";
										print SINES "Threshold: \t abs($flat_before_threshold  + 1 + $sine_length_lower_threshold + 1 + $flat_after_threshold)\n\n";

										print SINES "-------------------------------------------------------------------\n\n";
									} # verbose out
									
									# STRICT #
									if ($debugging_mode eq "false")
									{print SINES "$chromosome\t$cliff_position\t$sine_score\tSTRICT\n";}	
									 
									 if ($debugging_mode eq "true")
									{
										print SINES "$chromosome\t$cliff_position\t$sine_score\tSTRICT";
										
										print SINES "\t$total_change_before\t$total_change_during\t$total_change_after\t$cliff_difference";
										
										print SINES "\t$cliff_1_size_abs\t$mean_covar_whole_window\t$mean_covar_during\t$sine_repeat_length";
						
										print SINES "\t$mean_coverage_affected\t$mean_coverage_normal\t$flat_before_count\t$flat_after_count\t$cliff_1_count\t$cliff_2_count\n";
						
									}
				
									##############################################
									# Write to sines_plot_out for later plotting #
									##############################################
									print SINES_PLOT "POSSIBLE SINE\n";
									for($plot_count=0;$plot_count<=$window_size;$plot_count++)
									{
										$single_line = $file_array[$line_count + $plot_count];
										chomp $single_line;
										@item=split(/\t/,$single_line);
										$position = $item[1];
										$coverage_ratio = $item[5];
										
										print SINES_PLOT "$position\t$coverage_ratio\n";
									}
								
								}# flat during count
		
							} # flat counts
						
						} # Covar threshold
						
					} # sine repeat length 
					
				} # cliff difference
			
			} # cliff counts
			
		}# If over mean threshold
		
		
		if ($cliff_position > 99964974111) #################  111
		{

			print "Line: $line_count\t$cliff_position\n";
			print "$coverage_string\n\n";
			print "CRC flat threshold: $crc_flat_threshold\tCRC flat during threshold: $crc_flat_during_threshold\n\n";
			print "\tFlat before count:     \t$flat_before_count (out of $flat_before_threshold)\n";
			print "\tCliff 1 count:         \t$cliff_1_count\n";
			print "\tFlat during count:     \t$flat_during_count ($sine_length_lower_threshold - $sine_length_upper_threshold)\n";
			print "\tCliff 2 count:         \t$cliff_1_count\n";
			print "\tFlat after count:      \t$flat_after_count (out of $flat_after_threshold)\n";
			print "\tCliff 1:               \t$cliff_1_size (must be more than $crc_cliff_threshold)\n";
			print "\tCliff 2:               \t$cliff_2_size (must be more than $crc_cliff_threshold)\n";
			print "\tCliff difference:      \t$cliff_difference (must be less then $cliff_difference_threshold)\n";
			print "\tSine repeat length:    \t$sine_repeat_length\n";
			print "\tCovar whole window:    \t$covar_whole_window (mean is $mean_covar_whole_window)\n";
			print "\tCovar during:			\t$covar_during (mean is $mean_covar_during)\n";
			print "\tCovar all samples:     \t$covar_all_samples\n\n";
	
			print "SINE SCORE criteria:\n\n";
			print "\tMean change before:    \t$mean_change_before\n";
			print "\tMean change during:    \t$mean_change_during\n";
			print "\tMean change after:     \t$mean_change_after\n";
			print "\tCliff difference:      \t$cliff_difference\n";
			print "\tCliff size 1 abs:      \t$cliff_1_size_abs (divided by 2)\n";
			print "\tCliff size 2 abs:      \t$cliff_2_size_abs (divided by 2)\n";
			print "\tMean covar during:     \t$mean_covar_during\n";
			print "\tSine score:            \t$sine_score (must be >= $sine_score_threshold\n\n";
			if ($saved_by_sine_score eq "true"){print "\t\t>>>>> Saved by SINE_SCORE. $sine_score >= $sine_score_threshold\n\n";}
			
			print "INCLUSION criteria:\n\n";
			
			print "\tFlat before count:     \t$flat_before_count\n";
			print "\tCliff 1 count:         \t$cliff_1_count\n";
			print "\tFlat during count:     \t$flat_during_count\n";
			print "\tCliff 2 count:         \t$cliff_1_count\n";
			print "\tFlat after count:      \t$flat_after_count\n\n";
			
			print "\tInclusion score:       \t$inclusion_score\n";
			print "\tInclusion threshold:   \t$inclusion_threshold\n\n";
			if ($saved_by_inclusion eq "true"){print "\t\t>>>>> Saved by INCLUSION: $inclusion_score >= $inclusion_threshold\n\n";}
	
			print "STRICT criteria:\n\n";
			print "\tMean affected:         \t$mean_coverage_affected (threshold is $mean_sine_coverage_threshold)\n";
			print "\tMean normal:           \t$mean_coverage_normal (threshold is $mean_sine_coverage_threshold)\n";
			print "\tCliff 1 count:         \t$cliff_1_count (must be 1)\n";
			print "\tCliff 2 count:         \t$cliff_2_count (must be 1)\n";
			print "\tCliff difference:      \t$cliff_difference (threshold is $cliff_difference_threshold)\n";
			print "\tSine repeat length:    \t$sine_repeat_length ($sine_length_lower_threshold - $sine_length_upper_threshold)\n";
			print "\tMean covar whole wind: \t$mean_covar_whole_window (must be less than or equal to $covar_whole_window_threshold)\n";
			print "\tMean covar SINE:       \t$mean_covar_during (not used in STRICT currently)\n";

			print "\tSTRICT also uses the inclusion score threshold:\n\n";
			print "\tInclusion score:       \t$inclusion_score\n";
			print "\tInclusion threshold:   \t$inclusion_threshold\n\n";
			
			if ($saved_by_strict eq "true"){print "\t\t>>>>> Saved by STRICT\n\n";}
			
			$answer=<STDIN>;
			
		} # if cliff position > 6497114 etc
		
		
	} # if $line_count > 1
}

close IN;
close SINES;
close SINES_PLOT;


#################################
# Calculate various percentages #
#################################

if ($no_of_lines_large_out > 0 ){$percentage_kept_all = sprintf('%.1f', (100* ($keep_for_coverage_all_count / $no_of_lines_large_out)));}

if ($no_of_lines_coverage_all_out > 0 ){
	$percentage_t_test_lost = sprintf('%.6f', (100 * ($t_test_filtered_count / $no_of_lines_coverage_all_out)));
	$percentage_covar_lost = sprintf('%.6f', (100 * ($covar_filtered_count / $no_of_lines_coverage_all_out)));
	$percentage_coverage_lost = sprintf('%.6f', (100 * ($coverage_filtered_count / $no_of_lines_coverage_all_out)));
	$percentage_kept_assoc = sprintf('%.6f', (100* ($keep_for_assoc_count / $no_of_lines_coverage_all_out)));
}


$percentage_strict = 100 * ($saved_by_strict_count / $no_of_lines_sine_data);
$percentage_strict = sprintf('%.5f', $percentage_strict);
$percentage_inclusion = 100 * ($saved_by_inclusion_count / $no_of_lines_sine_data);
$percentage_inclusion = sprintf('%.3f', $percentage_inclusion);
$percentage_sine_score = 100 * ($saved_by_sine_score_count / $no_of_lines_sine_data);
$percentage_sine_score = sprintf('%.3f', $percentage_sine_score);

print "\n\n";
print "########################\n";
print "#  Coverage FINISHED   #\n";
print "########################\n\n";

print "Input list file:                      \t$list_file\n";
print "Output coverage all file:             \t$coverage_all_out\n\n\n";
print "Output coverage association file:     \t$output_assoc_file\n\n\n";

if (($files_to_use eq "bam") || ($files_to_use eq "coverage_all_large"))
{
	if ($filter_large eq "true")
	{
		print "Filtering to create $coverage_all_out:\n\n";
		print "\tNo. retained because of coverage threshold:    \t$keep_for_coverage_all_count ($percentage_kept_all%)\n";
		
		print "\tNo. kept for coverage_all.out:                 \t$keep_for_coverage_all_count\n";
		print "\tNo. lost for coverage_all.out:                 \t$lose_for_coverage_all_count\n";
		print "\tNo. of lines in coverage_all_large.out:                     \t$no_of_lines_large_out\n\n";
	}
}

if ($files_to_use ne "sine")
{
	print "\nFiltering to select what goes into assoc file:\n\n";
	print "\tOutput assoc file:     \t$output_assoc_file\n\n";

	print "\tTotal number of lines:                            \t$no_of_lines_coverage_all_out\n";
	print "\tNo. filtered out because of T-test threshold:     \t$t_test_filtered_count ($percentage_t_test_lost%)\n";
	print "\tNo. filtered out because of variation threshold:  \t$covar_filtered_count ($percentage_covar_lost%)\n";
	print "\tNo. filtered out because of Coverage threshold:   \t$coverage_filtered_count ($percentage_coverage_lost%)\n";
	print "\tNo. retained for coverage_assoc.out file:         \t$keep_for_assoc_count ($percentage_kept_assoc%)\n\n\n";
}

print "Filtering to select best possible SINES for ASSOC file:\n\n";
print "\tPossible SINES file:  \t$possible_sines_out\n\n";

print "\tRetained because of sine score:                  \t$saved_by_sine_score_count ($percentage_sine_score %)\n";
print "\tRetained because of strict criteria:             \t$saved_by_strict_count ($percentage_strict %)\n\n";
print "\tSINE stringency:                                 \t$sine_stringency\n\n";

print "\n\nDetails of which filtering methods filtered the SINE results:\n\n";
print "\tTotal number of lines of SINE data: \t$no_of_lines_sine_data\n";
print "\tsine_coverage_filtered_count:     \t$sine_coverage_filtered_count\n";
print "\tsine_cliff_filtered_count:        \t$sine_cliff_filtered_count\n";
print "\tsine_cliffdiff_filtered_count:    \t$sine_cliffdiff_filtered_count\n";
print "\tsine_repeat_length_filtered_count:\t$sine_repeat_length_filtered_count\n";
print "\tsine_covar_filtered_count:        \t$sine_covar_filtered_count\n";
print "\tsine_flat_filtered_count:         \t$sine_flat_filtered_count\n";
print "\tsine_flat_during_filtered_count:  \t$sine_flat_during_filtered_count\n\n";



###############################
# COMMAND_LOG final additions #
###############################

print COMMAND_LOG "Input list file:                  \t$list_file\n";
print COMMAND_LOG "Output coverage all file:         \t$coverage_all_out\n\n\n";

if (($files_to_use eq "bam") || ($files_to_use eq "coverage_all_large"))
{
	if ($filter_large eq "true")
	{
		print COMMAND_LOG "Filtering to create $coverage_all_out:\n\n";
		print COMMAND_LOG "\tNo. retained because of coverage threshold:    \t$keep_for_coverage_all_count ($percentage_kept_all%)\n";
		
		print COMMAND_LOG "\tNo. kept for coverage_all.out:                 \t$keep_for_coverage_all_count\n";
		print COMMAND_LOG "\tNo. lost for coverage_all.out:                 \t$lose_for_coverage_all_count\n";
		print COMMAND_LOG "\tNo. of lines in coverage_all_large.out:                     \t$no_of_lines_large_out\n\n";
	}
	if ($filter_large eq "false")
	{
		print COMMAND_LOG "Filtering to create $coverage_all_out:   NONE\n\n";
	}
}

print COMMAND_LOG "\nFiltering to select what goes into assoc file:\n\n";
print COMMAND_LOG "\tOutput assoc file:     \t$output_assoc_file\n\n";

print COMMAND_LOG "\tTotal number of lines:                            \t$no_of_lines_coverage_all_out\n";
print COMMAND_LOG "\tNo. filtered out because of T-test threshold:     \t$t_test_filtered_count ($percentage_t_test_lost%)\n";
print COMMAND_LOG "\tNo. filtered out because of variation threshold:  \t$covar_filtered_count ($percentage_covar_lost%)\n";
print COMMAND_LOG "\tNo. filtered out because of Coverage threshold:   \t$coverage_filtered_count ($percentage_coverage_lost%)\n";
print COMMAND_LOG "\tNo. retained for coverage_assoc.out file:         \t$keep_for_assoc_count ($percentage_kept_assoc%)\n\n\n";


print COMMAND_LOG "Filtering to select best possible SINES for SINE file:\n\n";
print COMMAND_LOG "\tPossible SINES file:  \t$possible_sines_out\n\n";

print COMMAND_LOG "\tRetained because of sine score:                  \t$saved_by_sine_score_count ($percentage_sine_score %)\n";
print COMMAND_LOG "\tRetained because of strict criteria:             \t$saved_by_strict_count ($percentage_strict %)\n";
print COMMAND_LOG "\tSINE stringency:                                 \t$sine_stringency\n\n";

print COMMAND_LOG "\n\nDetails of which filtering methods filtered the SINE results:\n\n";
print COMMAND_LOG "\tTotal number of lines of SINE data: \t$no_of_lines_sine_data\n";
print COMMAND_LOG "\tsine_coverage_filtered_count:     \t$sine_coverage_filtered_count\n";
print COMMAND_LOG "\tsine_cliff_filtered_count:        \t$sine_cliff_filtered_count\n";
print COMMAND_LOG "\tsine_cliffdiff_filtered_count:    \t$sine_cliffdiff_filtered_count\n";
print COMMAND_LOG "\tsine_repeat_length_filtered_count:\t$sine_repeat_length_filtered_count\n";
print COMMAND_LOG "\tsine_covar_filtered_count:        \t$sine_covar_filtered_count\n";
print COMMAND_LOG "\tsine_flat_filtered_count:         \t$sine_flat_filtered_count\n";
print COMMAND_LOG "\tsine_flat_during_filtered_count:  \t$sine_flat_during_filtered_count\n\n";


close COMMAND_LOG;

exit;


#############################################
# Subroutine to execute unix command        #
#############################################

sub run_unix_command
{
	my $unix_command = "";
	$unix_command = $_[0];
	print "\n";
	print("$unix_command\n");
	system("$unix_command");
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	
	print COMMAND_LOG "$unix_command\n";

}

####################################################################
#                                                                  #
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test" from "test.fasta")                           #
#                                                                  #
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


#############################################
# Subroutine to delete files                #
#############################################

sub delete_file
{
	my $file_to_be_deleted = "";	

	$file_to_be_deleted = $_[0];
	
	if (-e "$file_to_be_deleted")
	{
		$command = "rm  $file_to_be_deleted";
		system("$command");
	}

}

