# !/usr/bin/perl -w

#########################################################################
#									                                    #      
#	deletion_plot      						                        	#     
#									                                    #
#	Looks at deletions in BAm files (from the CIGAR string)             #
#   across your cases and controls, and provides data for               #
#	an association plot    										        #
#									                                    #
#########################################################################

#############################
# Mike Boursnell Feb 2013   #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use warnings;

use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;
use Cwd;
my $version							= "3";

#Various Fixed Parameters 
my $deletion_length_threshold		= 20;
my $default_window_size				= 500;
my $include_individual_columns		= "true";
my $t_test_threshold				= 2; # Only T-test values over this go into the ASSOC file
my $show							= "false"; #for debugging
my $ref								= "/home/genetics/canfam2/canfam2.fasta";


#new for deletion_plot
my $cigar							= "";
my $deletion_length					= 0;
my $D_position						= 0; # position of "D" in CIGAR string
my $position_start					= 0;
my $position_end					= 0;
my $deletion_size					= 0;


my $window_size						= 0;
my $files_to_use					= "";
my $prefix							= "";
my $last_prefix						= "";
my $last_chromosome					= "";
my $command_log						= "";
my $single_line						= "";
my $array_size						= "";
my $list_count						= 0;
my $no_of_files						= 0;
my $status							= "";
my $no_of_affecteds					= 0;
my $no_of_normals					= 0;
my $second_column_found				= "false";
my $line_count						= 0;
my $position						= 0;

my $no_of_lines_deletions_out		= 0;
my $no_of_lines_windows_file		= 0;
my $no_of_lines_sam_file			= 0;
my $no_of_lines_all_file			= 0;
my $input_string					= "";
my $first_position					= 0;
my $last_position					= 0;
my $window_start					= 0;
my $window_end						= 0;
my $window_count					= 0; # Maybe replace with deletion_count ???
my $last_line_checked				= 0;
my $last_line_checked_even			= 0;
my $last_line_checked_odd			= 0;
my $position_count					= 0;
my $deletion_count					= 0; # number of deletion in a window
my $something_in_window				= "false";
my $no_of_blocks					= 0;
my $block_count						= 0;
my $affected_total_deletion			= 0;
my $normal_total_deletion			= 0;
my $normalise_deletion				= "false";
my $bad_line_count					= 0;
my $even_check						= 0;
my $even							= "false";
my $check_count						= 0;
my $next_character					= "";

#Statistics
my $mean_deletion_affected		= 0;
my $mean_deletion_normal			= 0;
my $mean_deletion_both			= 0;
my $sum_squares_affected			= 0;
my $sum_squares_normal 				= 0;
my $variance_affected				= 0;
my $variance_normal					= 0;
my $covar_affected					= 0;
my $covar_normal					= 0;
my $deletion_ratio				= 0;
my $sd_affected						= 0;
my $sd_normal						= 0;
my $t_test							= 0;
my $covar_all_samples				= 0;


# FIle names #
my $sam_file						= "";
my $deletions_out					= ""; # deletion output file
my $deletion_bam					= "";
my $deletion_all_out				= "";
my $deletion_window_out				= "";
my $list_file						= "";
my $command							= "";
my $bam_file						= "";
my $bai_file						= "";
my $bam_file_in_results				= "";
my $bai_file_in_results				= "";
my $output_assoc_file				= "";
my $bam_file_in_folder				= "";
my $bam_files_directory				= "";
my $bam_file_in_bam_files			= "";

my $vcf_merge_string				= "";
my $run_title						= "";
my $chromosome						= "";
my $UG_region_string				= ""; #This is the string for the UnifiedGenotyper command line -L chr12
my $sample_name						= "";
my $answer							= "";
my $calls							= "";
my $species							= "";
my $ref_seq_name					= "";
my $other_ref_sequence 				= "false";
my $current_directory 				= "";
my $folder 							= "";
my $deletion_string				= "";


my @bam_file_array					= ();
my @sample_name_array				= ();
my @status_array					= ();
my @item							= ();
my @file_array						= ();
my @line_array						= (); # for making merged file
my @deletion_array				= ();

#######################################
# Create a name for the BAM directory #
#######################################
$bam_files_directory = "$ENV{HOME}/bam_files";

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "#######################################\n";
print color 'bold white';
print "            deletion_plot                \n";
print color 'bold magenta';
print "#######################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This program looks for deletions in a number of BAM files\n";
print "    and then plots an association plot for differences between\n";
print "    cases and controls.\n";

print color 'reset';
print "\n\n";
print "The input is a file with a list of the original BAM file names, with disease status in the second column (tab-delimited).\n\n";



until (-e $list_file)
{
	print "Name of the file of BAM file names:      ";
	$list_file = <STDIN>;
	chomp $list_file;
	
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

print "\n\n~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Which do you want to do?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   <1>  Use BAM files to start analysis from the beginning\n";
print "   <2>  Use SAM files (if you have already created them)\n";
#print "   <3>  Use 'deletions.out files derived from the SAM files.\n";
#print "   <4>  Use 'windows.out' files which have been grouped into windows.\n";
#print "   <5>  Use 'deletion_all.out' files which has multicolumns of deletion lists.\n\n";

$answer=<STDIN>;
chomp $answer;
if (substr($answer,0,1) eq "1" ){$files_to_use = "bam"}
if (substr($answer,0,1) eq "2" ){$files_to_use = "sam"}
if (substr($answer,0,1) eq "3" ){$files_to_use = "deletions"}
if (substr($answer,0,1) eq "4" ){$files_to_use = "windows"}
if (substr($answer,0,1) eq "5" ){$files_to_use = "deletion_all_out"}


print "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Please enter a name for this analysis (with no spaces) \n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
$run_title = <STDIN>;
chomp $run_title;


print "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Please specify a minimum size for deletions to look for (default=1) \n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
$deletion_length_threshold = <STDIN>;
chomp $deletion_length_threshold;

if ($deletion_length_threshold eq "") { $deletion_length_threshold = 1}


print "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Please specify a window size for grouping deletions (default=1000) \n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
$window_size = <STDIN>;
chomp $window_size;

if ($window_size eq "") { $window_size = 1000}



##############
# File names #
##############
$deletion_all_out = "$run_title"."_deletion_all.out";

#############################################
# Make sure the list file is in Unix format #
#############################################

$command = "dos2unix $list_file";
#print("\n$command\n");
system("$command");
print "\n";


####################################################
# Open the list file to get the list of file names #
####################################################
open (LIST, "$list_file") || die "Cannot open $list_file";

while ($single_line = <LIST> ) 
{
	$list_count = $list_count + 1;
	
	chomp $single_line;
	
	@item=split(/\t/,$single_line);
		
	$array_size = scalar @item;
	
	if ($array_size == 1)
	{
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
	$prefix = &get_prefix($bam_file);
	
	
	######################################################
	# Add file type suffix .bam if user hasn't added it  #
	######################################################

	if (index($bam_file,".bam") == -1 ){$bam_file = $bam_file.".bam"}

	#################################################################
	# Check file exists in current directory or bam_files directory #
	#################################################################
	$bam_file_in_bam_files = "$bam_files_directory/$bam_file";
	$bam_file_array[$list_count]=$bam_file;
	$status_array[$list_count]=$status;
	
	if (-e $bam_file)
	{
		print "\tFile $bam_file exists in the current directory\n";
		#$bam_file_array[$list_count]=$bam_file;
	}
	if ((!-e $bam_file) && (-e $bam_file_in_bam_files))
	{
		print "\tFile $bam_file exists in the 'bam_files' directory\n";
		#$bam_file_array[$list_count]=$bam_file_in_bam_files;
	}
	if ((!-e $bam_file) && (!-e $bam_file_in_bam_files))
	{
		print "\nERROR: File $bam_file not in the current directory or the 'bam_files' directory\n\n";
		exit;
	}

	
} # end of while loop for reading list_file (file of file names of BAM files)

close LIST;

$no_of_files = $list_count;

print "\n\n";
print "#################\n";
print "# Check details #\n";
print "#################\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print "$bam_file_array[$list_count]\t\t$status_array[$list_count]\n";
}

print "\n\nNo of BAM files: $no_of_files\n\n";

print "Deletion length threshold: \t\t$deletion_length_threshold\n\n";

print "Window size for grouping:  \t\t$window_size\n\n";

print "  > Press return to continue.\n\n";
$answer=<STDIN>;


$command_log = $run_title."deletion_command_log.out";
	
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

print COMMAND_LOG "COMMAND LOG for deletion_plot.pl $version\n\n";

print COMMAND_LOG "Run title:                \t$run_title\n\n";

print COMMAND_LOG "Files_to_use:             \t$files_to_use\n\n";

print COMMAND_LOG "Deletion length threshold:\t$deletion_length_threshold\n\n";

print COMMAND_LOG "Window size:              \t$window_size\n\n";

print COMMAND_LOG "Input file:               \t$list_file\n\n";

print COMMAND_LOG "List of BAM files:\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print COMMAND_LOG "File $list_count	\t$bam_file_array[$list_count]\t\t$status_array[$list_count]\n";
}
print COMMAND_LOG "\n\n";


if ($files_to_use eq "bam")
{
	print "\n\n";
	print "#------------------------------------------------------------#\n";
	print "# Using samtools to make SAM files from the BAM files...     #\n";
	print "#------------------------------------------------------------#\n";

	#####################################################
	# Stage 1                                           #
	# deletion_plot  makes SAM files from each BAM file #                   
	#####################################################
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		$bam_file = $bam_file_array[$list_count];
		$prefix = &get_prefix($bam_file_array[$list_count]);
		$sam_file = $prefix.".sam";
		$deletion_bam = $prefix."_deletion.bam";
		#$deletion_bam_txt = $prefix."_deletion_bam.txt";
		$bam_file_in_bam_files = "$bam_files_directory/$bam_file";
		
		$deletions_out = "$run_title"."_"."$prefix"."_deletion.out";
		
		print "Input BAM file:     \t$bam_file\n";
		print "Output deletion file:  \t $deletions_out\n\n";
		
		if (-e $bam_file)
		{
			&run_unix_command("/opt/samtools/samtools view $bam_file -o $sam_file");
		}
	
	
		if ((!-e $bam_file) && (-e $bam_file_in_bam_files))
		{
			&run_unix_command("/opt/samtools/samtools view $bam_file_in_bam_files -o $sam_file");
		}
		print "\n";
	}
}


if (($files_to_use eq "bam") || ($files_to_use eq "sam"))
{
	
	
	print "\n\n";
	print "#------------------------------------------------------------------#\n";
	print "# Using SAM files to look for deletions in CIGAR strings           #\n";
	print "#------------------------------------------------------------------#\n\n";

	#################################################
	# Stage 2                                       #
	# Use SAM files to look at CIGAR strings        #
	# and make deletions list files (deletions.out) #
	#################################################
	$first_position = 100000000;
	$last_position = 0;
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		$bam_file = $bam_file_array[$list_count];
		$prefix = &get_prefix($bam_file_array[$list_count]);
		$sam_file = $prefix.".sam";

		$deletions_out = "$run_title"."_"."$prefix"."_deletions.out";
		
		if (! -e $sam_file)
		{
			print "File $sam_file could not be found.\n\n\n";
			exit;
		}
	
		print "\n---------------------------------------------------------\n";
		print "Reading file $list_count:\t$sam_file into memory...\n";

		open (IN_SAM, "$sam_file") || die "Cannot open $sam_file\n\n";
			@file_array = <IN_SAM>;
			$no_of_lines_sam_file = scalar @file_array;
		close IN_SAM;
		
		open (OUT_DELETIONS, ">$deletions_out") || die "Cannot open $deletions_out\n\n";
		
		print "\tNo. of lines in SAM file $sam_file: $no_of_lines_sam_file\n\n";
		
		###########################################
		# Looping through file_array for SAM file #
		###########################################
		print "Writing file $list_count:\t$deletions_out\n";
		
		for ($line_count = 0;$line_count<=$no_of_lines_sam_file - 1; $line_count=$line_count + 2)
		{
			$single_line = $file_array[$line_count];
			chomp $single_line;
			
			@item=split(/\t/,$single_line);
			$array_size = scalar @item;

			$chromosome      = $item[2];
			$position        = $item[3];
			$cigar   		 = $item[5];

			if ($position < $first_position){$first_position = $position};
			if ($position > $last_position){$last_position = $position};
			
			if (($line_count % 100000) == 0){print "\tSearching SAM file $list_count/$no_of_files\tLine: $line_count/$no_of_lines_sam_file\n";}
			
			$D_position = index($cigar,"D");
			$next_character = substr($cigar,$D_position + 2,1);
			
			if ($D_position > -1)
			{
				#$deletion_length = substr($cigar,$D_position + 1,1);
				#print "CIGAR: $cigar\tDeletion_length: $deletion_length\n";
				

				if ($next_character =~ /^\d+?$/)
				{
					$deletion_length = substr($cigar,$D_position + 1,2);
				}
				else
				{
					$deletion_length = substr($cigar,$D_position + 1,1);
				}
				
				#print "CIGAR: $cigar\tDEL: $deletion_length\t$next_character\n";
				
				if ($deletion_length > $deletion_length_threshold)
				{
					print OUT_DELETIONS "$chromosome\t$position\t$cigar\n";
				}
			}
		}
		
		print "\n\nFor file $sam_file:\n\n";
		print "\tLowest position:  \t$first_position\n";
		print "\tHighest position: \t$last_position\n\n";

		close OUT_DELETIONS;
	}
}



############################################
# These are the overall lowest and highest #
# positions from all the deletions files   #
############################################

print "\tLowest position in all SAM files:  \t$first_position\n";
print "\tHighest position in all SAM files: \t$last_position\n";


if (($files_to_use eq "bam") || ($files_to_use eq "sam") || ($files_to_use eq "deletions"))
{
	###################################################
	# Stage 3                                         #
	# Reducing the deletions text file          	  #
	# to a file which contains a list of windows      #
	# and the number of deletions within each window  #
	###################################################
	
	print "\n\n";
	print "#---------------------------------------------------------#\n";
	print "# Grouping data from deletions.out files into windows...  #\n";
	print "#---------------------------------------------------------#\n\n";
	
	$window_size=$default_window_size;
	
	##################################
	# Loop through the list of files #
	##################################
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		$prefix = &get_prefix($bam_file_array[$list_count]);
		$deletions_out = "$run_title"."_"."$prefix"."_deletions.out";
		$deletion_window_out = "$run_title"."_"."$prefix"."_window.out";
		
		open (IN_DELETIONS, "$deletions_out") || die "Cannot open $deletions_out\n\n";
			@file_array = <IN_DELETIONS>;
			$no_of_lines_deletions_out = scalar @file_array;
		close IN_DELETIONS;
		
		open (OUT_WINDOW,">$deletion_window_out")|| die "Cannot create output file: $deletion_window_out";
		
		print OUT_WINDOW "FILE_NAME\tCHROMOSOME\tPOS_START\tWINDOW_START\tNO_IN_WINDOW\n";
		
		print "\n\n-------------------------------------------------------------------\n";
		print "Grouping file $list_count: $deletions_out into windows of $window_size...\n\n";
		
		#########################
		# Set up initial states #
		#########################
		$window_start = $first_position;
		$window_end = $window_start + $window_size;
		$last_line_checked = 1;
		$last_line_checked_even = 1;
		$last_line_checked_odd = 1;
		
		$even_check= 0; # first $even = "false";
		
		##########################################################
		# Move forward a window at a time                        #
		# (overlapping each time by half a window)               #
		# (using file stored in @file_array)                     #
		#                                                        #
		# At each window position look down the whole of         #
		# the deletions file and see which lines fit within it.  #
		##########################################################
		$check_count = 0;
		
		# This looop is a loop of windows
		for ($position_count=$first_position;$position_count<=$last_position;$position_count = $position_count + ($window_size/2))
		{
			##########################################
			# Set start and end points of the window #
			##########################################
			$window_start = $position_count;
			$window_end = $position_count + $window_size;
			
			$window_count = 0;
			
			$check_count = $check_count + 1;
			
			if (($check_count % 1000) == 0){print "\tGrouping file $list_count/$no_of_files into windows. Position: $position_count/$last_position\n";}
			
			# Work out if this an even or an odd window #
			$even_check = $even_check + 1;
			
			if (($even_check % 2) == 0) {$even="true"} else {$even = "false"};
			
			if ($even eq "true"){$last_line_checked = $last_line_checked_even;}
			if ($even eq "false"){$last_line_checked = $last_line_checked_odd;}
				
			#################################################################################
			# Now loop through the deletions file to find any positions that are in         #
			# this window  (There may be some windows for which there are no positions)     #
			# Once past window_end then break out of the loop                               #
			#################################################################################
			
			### SHOULD WE USE LAST_LINE_CHECKED HERE?????  ###
			
			$something_in_window="false";
			
			for ($line_count = $last_line_checked;$line_count<$no_of_lines_deletions_out;$line_count++)
			{
				$single_line = $file_array[$line_count];
				chomp $single_line;
			
				@item=split(/\t/,$single_line);
				$chromosome = $item[0];
				$position = $item[1];
			
				if ($even eq "true"){$last_line_checked_even = $line_count;}
				if ($even eq "false"){$last_line_checked_odd = $line_count;}
				
				#print "\tLine:$line_count\tPos: $position\tWind: $window_start-$window_end\n";
				
				if ($position < $window_end)
				{
					##########################################
					# Check if position is within the window #
					##########################################
					if (($position >= $window_start) && ($position < $window_end))
					{
						$window_count = $window_count + 1;
						
						#print "\t\tFOUND Line:$line_count\tPos: $position\tWind: $window_start-$window_end\t$window_count\n";
						
						$something_in_window="true";
						#$position_found = $position_start;
						
						if ($even eq "true"){$last_line_checked_even = $line_count;}
						if ($even eq "false"){$last_line_checked_odd = $line_count;}
						
						#print "File: $list_count\tLine: $line_count\tPosition: $position_start $window_start > $position_start > $window_end\t$window_count\n";
						
					}

				} # if position_start is < window_end
				
				#####################################################
				# If position is beyond the end of the window then  #
				# write to the output window file.  Then jump out   #
				# of the loop (no point going any further)          #
				#####################################################
				if ($position >= $window_end)
				{
					#print "END_LOOP.  $position > $window_end\n\n";
					goto END_LOOP;
				}
				
				
			} # $line_count loop
			
			END_LOOP: # jump out of $line_count loop to here. Past end of this window#
			
			if ($something_in_window eq "true")
			{
				print OUT_WINDOW "$prefix\t$chromosome\t$window_start\t$window_count\n";
			}
			
			if ($something_in_window eq "false")
			{
				print OUT_WINDOW "$last_prefix\t$last_chromosome\t$window_start\t0\n";
			}
			
			$something_in_window = "false";
			
			$window_count = 0;
			$last_line_checked = $line_count;
			$last_prefix = $prefix;
			$last_chromosome = $chromosome;
			
			
		} # Position_count loop
		
		close OUT_WINDOW;
		

	} # $list_count loop	
}


if (($files_to_use eq "bam") || ($files_to_use eq "sam") || ($files_to_use eq "deletions") || ($files_to_use eq "windows"))
{
	print "\n\n";
	print "#---------------------------------------------------#\n";
	print "# Merging windows files into deletions_all.out...   #\n";
	print "#---------------------------------------------------#\n\n";
	
	######################################################
	# Stage 4                                            #
	# Load in files one at a time and add to @file_array #
	# to make the combined deletion_all.out file         #
	######################################################
	
	open (OUT_ALL,">$deletion_all_out")|| die "Cannot create output file: $deletion_all_out";
	
	# Write Headers to first line #
	print OUT_ALL "CHR\tWINDOW_START";
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		#$bam_file = $bam_file_array[$list_count];
		$prefix = &get_prefix($bam_file_array[$list_count]);
		$deletion_window_out = "$run_title"."_"."$prefix"."_window.out";
		
		#Column header for this file #
		print OUT_ALL "\t$prefix";
		
		# Read file into file_array #
		print "\n\n---------------------------------------------------------\n";
		print "\tReading file $list_count: $deletion_window_out into memory...\n";
		
		open (IN_WINDOW, "$deletion_window_out") || die "Cannot open $deletion_window_out\n\n";
			@file_array = <IN_WINDOW>;
			$no_of_lines_windows_file = scalar @file_array;
		close IN_WINDOW;
		
		# Work down file storing data into @line_array #
		for ($line_count = 1;$line_count <$no_of_lines_windows_file;$line_count++)
		{
			$single_line = $file_array[$line_count];
			chomp $single_line;
		
			@item=split(/\t/,$single_line);
			$prefix = $item[0];
			$chromosome = $item[1];
			$position = $item[2];
			$deletion_count = $item[3];
		
			if ($list_count == 1)
			{
				$line_array[$line_count] = "$chromosome\t$position\t$deletion_count";
			}
			if ($list_count > 1)
			{
				$line_array[$line_count] = "$line_array[$line_count]\t$deletion_count";
			}
			
		}
		
	} # $list_count loop
	
	# Work down file writing data from @line_array to OUT_ALL #
	print "\n\t\t==> Writing combined file $deletion_all_out...\n";
	print OUT_ALL "\n";
	
	for ($line_count = 1;$line_count<$no_of_lines_windows_file;$line_count++)
	{
		print OUT_ALL "$line_array[$line_count]\n";
	}
	
} # end if if $files_to_use eq etc (Stag e 4)



if (($files_to_use eq "bam") || ($files_to_use eq "sam") || ($files_to_use eq "deletion") || ($files_to_use eq "windows") || ($files_to_use eq "deletion_all_out"))
{
	################################################
	# Stage 5                                      #
	# Doing statistics and creating assoc.out file #
	################################################
	print "\n\n";
	

	print "#------------------------------------#\n";
	print "# Reading deletion_all.out file...   #\n";
	print "#                                    #\n";
	print "# Calculating deletion statistics    #\n";
	print "# and creating deletion_assoc.out    #\n";
	print "#------------------------------------#\n\n";
	# Read file into file_array #
	print "\tReading file $deletion_all_out\n";
	open (IN_ALL, "$deletion_all_out") || die "Cannot open $deletion_all_out\n\n";
		@file_array = <IN_ALL>;
		$no_of_lines_all_file = scalar @file_array;
	close IN_ALL;
	
	print "\n\tNo. of lines in $deletion_all_out:\t$no_of_lines_all_file\n\n";
	
	#######################################################
	# Open file for filtered deletion data and statistics #
	#######################################################
	$output_assoc_file="$run_title"."_deletion_w"."$window_size"."_assoc.out";
	
	open (OUT_ASSOC, ">$output_assoc_file") || die "Cannot open $output_assoc_file\n\n";
	
	###########
	# Headers #
	###########
	print OUT_ASSOC "CHR\tPOSITION";
	
	if ($include_individual_columns eq "true")
	{
		for ($list_count=1;$list_count <=$no_of_files;$list_count++)
		{
			$prefix = &get_prefix($bam_file_array[$list_count]);
			print OUT_ASSOC "\t$prefix";
		}
	} # if you want indivial DP data columns


	print OUT_ASSOC "\tAV_A\tAV_N\tCOV_A\tCOV_N\tRATIO\tT_TEST\n";

				
	#####################################################
	# Looping through file_array for deletion_all.out #
	#####################################################
	for ($line_count = 0;$line_count<$no_of_lines_all_file; $line_count++)
	{
		if ($line_count > 0)
		{
			$single_line = $file_array[$line_count];
			@item=split(/\t/,$single_line);
			$array_size = scalar @item;

			$chromosome=$item[0];
			$position=$item[1];
			
			if (($line_count % 1000) == 0){print "\tCreating assoc.out file.  Line: $line_count  \tPosition: $position\n";}

			$no_of_blocks = $array_size - 2;
			
			if ($no_of_blocks == $no_of_files)
			{
				$affected_total_deletion = 0;
				$normal_total_deletion = 0;
				#$affected_total_deletion_normalised = 0;
				#$normal_total_deletion_normalised = 0;
				
				$deletion_string = "";
				
				##################################################
				# Put together string of deletion to show user #
				##################################################
				for ($block_count=1;$block_count <=$no_of_blocks;$block_count++)
				{
					$deletion_array[$block_count] = $item[$block_count + 1];
					$deletion_array[$block_count] = sprintf('%.0f', $deletion_array[$block_count]); 

					$deletion_string = "$deletion_string"."\t"."$deletion_array[$block_count]";
					
					if ($status_array[$block_count] eq "affected")
					{
						$affected_total_deletion = $affected_total_deletion + $deletion_array[$block_count];
						#$affected_total_deletion_normalised = $affected_total_deletion_normalised + ($deletion_array[$block_count] * $weighting_array[$block_count]);
					}
					
					if ($status_array[$block_count] eq "control")
					{
						$normal_total_deletion = $normal_total_deletion + $deletion_array[$block_count];
						#$normal_total_deletion_normalised = $normal_total_deletion_normalised + ($deletion_array[$block_count] * $weighting_array[$block_count]);
					}
				} # $block_count loop
				
				############################################################
				# Mean deletion counts on each line of deletion_all.out #
				############################################################
				if ($normalise_deletion eq "false")
				{
					$mean_deletion_affected = $affected_total_deletion / $no_of_affecteds;
					$mean_deletion_normal = $normal_total_deletion / $no_of_normals;
					$mean_deletion_both = ($affected_total_deletion + $normal_total_deletion) / ($no_of_affecteds + $no_of_normals);
				}
				
				############
				# Variance #
				############
				$sum_squares_affected = 0;
				$sum_squares_normal = 0;
				$variance_affected = 0;
				$variance_normal = 0;
				$covar_affected = 0;
				$covar_normal = 0;
				$deletion_ratio = 0;

				for ($block_count=1;$block_count <=$no_of_blocks;$block_count++)
				{
					if ($status_array[$block_count] eq "affected")
					{
						$sum_squares_affected = $sum_squares_affected + (($deletion_array[$block_count] - $mean_deletion_affected) ** 2);
					}
					if ($status_array[$block_count] eq "control")
					{
						$sum_squares_normal = $sum_squares_normal + (($deletion_array[$block_count] - $mean_deletion_normal) ** 2);
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
					$t_test = ($mean_deletion_affected - $mean_deletion_normal) / (sqrt((($variance_affected / $no_of_affecteds) + ($variance_normal / $no_of_normals))));
					$t_test = abs($t_test); # t_statistic can be positve or negative
				}
				
				############################
				# Coefficient of variation #
				############################
				if ($mean_deletion_affected > 0 ) {$covar_affected = $sd_affected / $mean_deletion_affected;}
				if ($mean_deletion_normal > 0) {$covar_normal = $sd_normal / $mean_deletion_normal;}

				$covar_all_samples = ($covar_affected + $covar_normal ) / 2;
				
				
				####################
				# deletion ratio #
				####################
				if ($mean_deletion_normal > 0) {$deletion_ratio = $mean_deletion_affected / $mean_deletion_normal;}


				################################
				# Format to two decimal places #
				################################ 
				$covar_affected = sprintf('%.2f', $covar_affected); 
				$covar_normal = sprintf('%.2f', $covar_normal); 
				$deletion_ratio = sprintf('%.2f', $deletion_ratio );
				$t_test = sprintf('%.2f', $t_test );
				$covar_all_samples = sprintf('%.2f', $covar_all_samples );
				$mean_deletion_affected = sprintf('%.2f', $mean_deletion_affected);
				$mean_deletion_normal = sprintf('%.2f', $mean_deletion_normal);
			
			
				######################################
				# Write to output ASSOC file         #
				# (if greater than T-test threshold) #
				######################################
				if ($t_test >= $t_test_threshold)
				{
					print OUT_ASSOC "$chromosome\t$position";
					if ($include_individual_columns eq "true"){print OUT_ASSOC "$deletion_string";}
					print OUT_ASSOC "\t$mean_deletion_affected\t$mean_deletion_normal";
					print OUT_ASSOC "\t$covar_affected\t$covar_normal";
					print OUT_ASSOC "\t$deletion_ratio\t$t_test";
					print OUT_ASSOC "\n";
				}
			} # if $no_of_blocks eq $no_of_files
			
			if ($no_of_blocks != $no_of_files)
			{
				$bad_line_count = $bad_line_count + 1;
			}
		}
		
	} # $line_count loop
	
	close OUT_ASSOC;
} # end if if $files_to_use eq etc

close COMMAND_LOG;
close OUT_ALL;

print "\n\n";
print "###############\n";
print "#  FINISHED   #\n";
print "###############\n\n";

print "Output files:\n\n";

print "File with all the deletion data:\n\n";  
print "\t$deletion_all_out\n\n\n";

print "deletion association file, for plotting association graph: \n\n";   	   
print "\t$output_assoc_file\n\n";

print "Window size used:\t$window_size\n\n\n";

exit;

#############################################
#                                           #
# Subroutine to execute unix command        #
#                                           #
#############################################

sub run_unix_command
{
	my $unix_command = "";
	my $step = "";	
	$unix_command = $_[0];
		
	print("   $unix_command\n");
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


sub find_file
{
	my $file_to_find 			= "";
	my $current_directory		= "";
	my $folder					= "";
	
	$file_to_find = $_[0];
	
	$current_directory=getcwd;
	#chomp $current_directory;
	opendir (DIR, $current_directory);

	while (my $folder = readdir(DIR)) 
	{
			if ((-d "$current_directory/$folder") && !($folder =~ m/^\./))
			{
				if (-e "$folder/$file_to_find")
				{
					$file_to_find = "$folder/$file_to_find";
					return ($file_to_find);
				}
			}

	}

	closedir (DIR);

		
}

#############################################
# Subroutine to delete files                #
#############################################

sub delete_file
{
	my $file_to_be_deleted = "";	

	$file_to_be_deleted = $_[0];
	$command = "rm  $file_to_be_deleted";
	print("Deleting file   $file_to_be_deleted\n");
	system("$command");

}



