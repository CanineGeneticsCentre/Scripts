#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	delly_plot      						                        	#     
#									                                    #
#	Extracts data from DELLY files and plots it					        #
#									                                    #
#########################################################################

#############################
# Mike Boursnell Feb 2013   #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;
use Cwd;
my $version							= "3";


#Various Fixed Parameters 
my $s_value							= 8;
my $default_window_size				= 2000;
my $include_individual_columns		= "true";
my $t_test_threshold				= 0.1; # Only T-test values over this go into the ASSOC file
my $show							= "false"; #for debugging
my $ref								= "/home/genetics/canfam2/canfam2.fasta";
my $deletion_size_threshold        	= 100; # for awk filtering
my $supporting_reads_threshold		= 3;   # for awk filtering
my $mapping_quality_threshold       = 5;  # for awk filtering

#new for delly
my $pass							= 0; # pass 1 = delly, pass 2 = invy, pass3 = duppy
my $delly_program					= ""; # This is delly, invy or duppy
my $position_start					= 0;
my $position_end					= 0;
my $deletion_size					= 0;
my $supporting						= 0; # No supporting reads
my $mapping_quality					= 0;
my $deletion_id						= "";
my $window_size						= 0;
my $files_to_use					= "";
my $delly_bam_txt 					= "";
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
my $no_of_lines_delly_bam_txt		= 0;
my $line_count						= 0;
my $position						= 0;
my $no_of_lines_filtered_out		= 0;
my $no_of_lines_simplified_out		= 0;
my $no_of_lines_windows_file		= 0;
my $no_of_lines_all_file			= 0;
my $input_string					= "";
my $first_position					= 0;
my $last_position					= 0;
my $window_start					= 0;
my $window_end						= 0;
my $window_count					= 0; # Maybe replace with delly_count ???
my $last_line_checked				= 0;
my $last_line_checked_even			= 0;
my $last_line_checked_odd			= 0;
my $position_count					= 0;
my $delly_count						= 0; # number of delly in a window
my $something_in_window				= "false";
my $no_of_blocks					= 0;
my $block_count						= 0;
my $affected_total_delly			= 0;
my $normal_total_delly				= 0;
my $normalise_delly					= "false";
my $bad_line_count					= 0;
my $even_check						= 0;
my $even							= "false";

#Statistics
my $mean_delly_affected				= 0;
my $mean_delly_normal				= 0;
my $mean_delly_both					= 0;
my $sum_squares_affected			= 0;
my $sum_squares_normal 				= 0;
my $variance_affected				= 0;
my $variance_normal					= 0;
my $covar_affected					= 0;
my $covar_normal					= 0;
my $delly_ratio						= 0;
my $sd_affected						= 0;
my $sd_normal						= 0;
my $t_test							= 0;
my $covar_all_samples				= 0;


# FIle names #
my $delly_out						= ""; # DELLY output file
my $delly_simplified_out			= ""; # simplified version of $delly_out
my $delly_filtered_out				= ""; # filtered version of delly_out

my $delly_bam						= "";
my $delly_all_out					= "";
my $delly_window_out				= "";
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
my $delly_string					= "";


my @bam_file_array					= ();
my @sample_name_array				= ();
my @status_array					= ();
my @item							= ();
my @file_array						= ();
my @line_array						= (); # for making merged file
my @delly_array				= ();

#############################################
# Create a name for the BAM files directory #
# (This is if you have all your BAM files   #
# in one place                              #
#############################################
$bam_files_directory = "$ENV{HOME}/bam_files";

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "#######################################\n";
print color 'bold white';
print "            delly_plot                \n";
print color 'bold magenta';
print "#######################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This program uses the output files from the program DELLY.\n";
print "    DELLY looks for deletions in NGS sequence by looking for\n";
print "    read pairs that are too far apart.\n";

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
print "   <1>  Use BAM files to start analysis from the beginning.\n";
print "   <2>  Use '_delly.out' files already made with DELLY.\n";
print "   <3>  Use '_filtered.out' files which have been filtered from the DELLY out files.\n";
#print "   <4>  Use '_simplified.out' files which have been simplified from the DELLY filtered files.\n";
#print "   <5>  Use '_window.out' files which have been created from the DELLY simplified files.\n\n";

$answer=<STDIN>;
chomp $answer;
if (substr($answer,0,1) eq "1" ){$files_to_use = "bam"}
if (substr($answer,0,1) eq "2" ){$files_to_use = "delly_out"}
if (substr($answer,0,1) eq "3" ){$files_to_use = "filtered_out"}
if (substr($answer,0,1) eq "4" ){$files_to_use = "simplified_out"}
if (substr($answer,0,1) eq "5" ){$files_to_use = "window_out"}


print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Which DELLY program do you want to use?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   <1>  DELLY for deletions\n";
print "   <2>  INVY for inversions\n";
print "   <3>  DUPPY for tandem duplications\n";
print "   <4>  All three (not yet implemented)\n\n";

$answer=<STDIN>;
chomp $answer;
if (substr($answer,0,1) eq "1" ){$delly_program = "delly"}
if (substr($answer,0,1) eq "2" ){$delly_program = "invy"}
if (substr($answer,0,1) eq "3" ){$delly_program = "duppy"}
if (substr($answer,0,1) eq "4" ){$delly_program = "all"}




print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Please enter a name for this analysis (with no spaces) \n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
$run_title = <STDIN>;
chomp $run_title;




#############################################
# Make sure the list file is in Unix format #
#############################################

$command = "dos2unix $list_file";
print("\n$command\n");
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
		print "File $bam_file exists in the current directory\n";
		#$bam_file_array[$list_count]=$bam_file;
	}
	if ((!-e $bam_file) && (-e $bam_file_in_bam_files))
	{
		print "File $bam_file exists in the 'bam_files' directory\n";
		#$bam_file_array[$list_count]=$bam_file_in_bam_files;
	}
	if ((!-e $bam_file) && (!-e $bam_file_in_bam_files))
	{
		print "\nERROR: File $bam_file not in the current directory or the 'bam_files' directory\n\n";
		exit;
	}

	
} # end of while loop for reading list_file (file of file names of BAM files)

$no_of_files = $list_count;

print "\nNo of BAM files: $no_of_files\n";
close LIST;


$command_log = "delly_command_log.out";
	
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "COMMAND LOG for delly.pl $version\n\n";

print COMMAND_LOG "Running with files_to_use = $files_to_use\n\n";

print COMMAND_LOG "Input file:        \t$list_file\n\n";
print COMMAND_LOG "List of BAM files:\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print COMMAND_LOG "File $list_count	\t$bam_file_array[$list_count]\t\t$status_array[$list_count]\n";
}
print COMMAND_LOG "\n\n";

#############################
# Three passes              #
# 1=delly, 2=invy, 3=duppy  #
#############################
for ($pass = 1;$pass <=3;$pass++)
{

	if ($pass == 1){$delly_program="delly"}
	if ($pass == 2){$delly_program="invy"}
	if ($pass == 3){$delly_program="duppy"}

	#########################
	# Output all file names #
	#########################
	if ($delly_program eq "delly"){$delly_all_out = "$run_title"."_delly_all.out"}
	if ($delly_program eq "invy"){$delly_all_out = "$run_title"."_invy_all.out"}
	if ($delly_program eq "duppy"){$delly_all_out = "$run_title"."_duppy_all.out"}


	if ($files_to_use eq "bam")
	{
		print "\n\n";
		if ($delly_program eq "delly")
		{
			print "#-------------------------------------------------------#\n";
			print "# Using DELLY to look for deletions in the BAM files...   #\n";
			print "#-------------------------------------------------------#\n";
		}
		if ($delly_program eq "invy")
		{
			print "#-------------------------------------------------------#\n";
			print "# Using INVY to look for inversions in the BAM files... #\n";
			print "#-------------------------------------------------------#\n";
		}
		if ($delly_program eq "duppy")
		{
			print "#-----------------------------------------------------------------#\n";
			print "# Using DUPPY to look for tandem duplications in the BAM files... #\n";
			print "#-----------------------------------------------------------------#\n";
		}
		if ($delly_program eq "all")
		{
			print "#-------------------------------------------------------#\n";
			print "# Using DELLY, INVY and DUPPY to look for deletions,    #\n";
			print "# inversions and duplications in the BAM files...       #\n";
			print "# (NOT YET IMPLEMENTED)                         .       #\n
			print "#-------------------------------------------------------#\n";
		}
		
		###############################################
		# Stage 1                                     #
		# DELLY reads the BAM files and creates       #
		# several output files
		###############################################
		for ($list_count=1;$list_count<=$no_of_files;$list_count++)
		{
			$bam_file = $bam_file_array[$list_count];
			$prefix = &get_prefix($bam_file_array[$list_count]);
			$bam_file_in_bam_files = "$bam_files_directory/$bam_file";
			
			
			$delly_bam = "$run_title"."_"."$prefix"."_"."$delly_program".".bam";
			$delly_bam_txt = "$run_title"."_"."$prefix"."_"."$delly_program"."_bam.txt";
			$delly_out = "$run_title"."_"."$prefix"."_"."$delly_program".".out";
			
			print "Input BAM file:     \t$bam_file\n";
			print "Output DELLY file:  \t $delly_out\n\n";
			
			if (-e $bam_file)
			{
				if ($delly_program eq "delly")
				{&run_unix_command("/opt/delly/delly -i $run_title -s $s_value -o $delly_out  $bam_file")}
				
				if ($delly_program eq "invy")
				{&run_unix_command("/opt/delly/invy -i $run_title -o $delly_out  $bam_file")}
				
				if ($delly_program eq "duppy")
				{&run_unix_command("/opt/delly/duppy -i $run_title -o $delly_out $bam_file")}
			}
		
		
			if ((!-e $bam_file) && (-e $bam_file_in_bam_files))
			{

				if ($delly_program eq "delly")
				{&run_unix_command("/opt/delly/delly -i $run_title -s $s_value -o $delly_out  $bam_file_in_bam_files")}
				
				if ($delly_program eq "invy")
				{&run_unix_command("/opt/delly/invy -i $run_title -o $delly_out  $bam_file_in_bam_files")}
				
				if ($delly_program eq "duppy")
				{&run_unix_command("/opt/delly/duppy -i $run_title -o $delly_out $bam_file_in_bam_files")}
				
			}
			
			print "\n";
		}
	}

	if (($files_to_use eq "bam") || ($files_to_use eq "delly_out"))
	{
		
		
		print "\n\n";
		print "#------------------------------------------------------------------#\n";
		print "# Filtering delly_out files using cat, grep and awk (unix stuff!)  #\n";
		print "#------------------------------------------------------------------#\n";

		###############################################
		# Stage 2                                     #
		# Use linux awk to filter the output files    #
		###############################################
		for ($list_count=1;$list_count<=$no_of_files;$list_count++)
		{
			$bam_file = $bam_file_array[$list_count];
			$prefix = &get_prefix($bam_file_array[$list_count]);
			
			$delly_out = "$run_title"."_"."$prefix"."_"."$delly_program".".out";
			$delly_filtered_out = "$run_title"."_"."$prefix"."_"."$delly_program"."_filtered.out";
			
			if (! -e $delly_out)
			{
				print "File $delly_out could not be found.\n\n\n";
				exit;
			}
		
			print "File $list_count:\t$delly_out\n";
			print "\t";

			
			&run_unix_command("cat $delly_out | grep '>' | awk '\$4>=$deletion_size_threshold && \$5>$supporting_reads_threshold && \$6>$mapping_quality_threshold' > $delly_filtered_out");

		}
	}

	if (($files_to_use eq "bam") || ($files_to_use eq "delly_out")|| ($files_to_use eq "filtered_out") )
	{
		###############################################
		# Stage 3                                     #
		# Simplify these delly filtered files so that #
		# they only contain FILENAME, CHR and POS     #
		###############################################
		
		print "\n\n";
		print "#------------------------------------------#\n";
		print "# Simplifying filtered.out files...        #\n";
		print "#------------------------------------------#\n\n";

		$first_position = 100000000;
		$last_position = 0;
		
		for ($list_count=1;$list_count<=$no_of_files;$list_count++)
		{
			$bam_file = $bam_file_array[$list_count];
			$prefix = &get_prefix($bam_file_array[$list_count]);

			$delly_filtered_out = "$run_title"."_"."$prefix"."_"."$delly_program"."_filtered.out";
			$delly_simplified_out = "$run_title"."_"."$prefix"."_"."$delly_program"."_simplified.out";
			
			
			print "\tSimplifying file $list_count: $delly_filtered_out...\n";
			
			open (IN_ALL, "$delly_filtered_out") || die "Cannot open $delly_bam_txt\n\n";
				@file_array = <IN_ALL>;
				$no_of_lines_filtered_out = scalar @file_array;
			close IN_ALL;
			
			open (OUT, ">$delly_simplified_out")|| die "Cannot create output file: $delly_simplified_out";
			
			print OUT "FILE\tCHR\tSTART\tEND\tSIZE\tQUAL\tID\n";
			
			#######################################################
			# Looping through file_array for delly_simplified_out #
			#######################################################
			for ($line_count = 0;$line_count<=$no_of_lines_filtered_out - 1; $line_count++)
			{
				$single_line = $file_array[$line_count];
				chomp $single_line;
				
				if (substr ($single_line, 0, 3) eq "chr")
				{
					@item=split(/\t/,$single_line);
					$array_size = scalar @item;

					$chromosome      = $item[0];
					$position_start  = $item[1];
					$position_end    = $item[2];
					$deletion_size   = $item[3];
					$supporting      = $item[4]; # not used in output file?
					$mapping_quality = $item[5];
					$deletion_id     = $item[6];
					
				} # if chr line
				
				if ($position_start < $first_position){$first_position = $position_start};
				if ($position_end > $last_position){$last_position = $position_end};
				
				print OUT "$prefix\t$chromosome\t$position_start\t$position_end\t$deletion_size\t$mapping_quality\t$deletion_id\n";	
			} # $line_count loop

		} # $list_count loop to simplfy delly text files
	}
	close OUT;

	print "\n\n";

	############################################
	# These are the overall lowest and highest #
	# positions from all the input files       #
	############################################
	print "These are calculated during the filtered--> simplified stage:\n\n";
	print "\tLowest position:  \t$first_position\n";
	print "\tHighest position: \t$last_position\n";


	if (($files_to_use eq "bam") || ($files_to_use eq "delly_out") || ($files_to_use eq "filtered_out") || ($files_to_use eq "simplified_out"))
	{
		###################################################
		# Stage 4                                         #
		# Reducing the simplified delly text file     	  #
		# to a file which contains a list of windows      #
		# and the number of deletions within each window  #
		###################################################
		
		print "\n\n";
		print "#---------------------------------------------------------#\n";
		print "# Grouping data from simplified_out files into windows... #\n";
		print "#---------------------------------------------------------#\n\n";
		
		$window_size=$default_window_size;
		
		##################################
		# Loop through the list of files #
		##################################
		for ($list_count=1;$list_count<=$no_of_files;$list_count++)
		{
			$prefix = &get_prefix($bam_file_array[$list_count]);
			
			$delly_simplified_out = "$run_title"."_"."$prefix"."_"."$delly_program"."_simplified.out";
			$delly_window_out = "$run_title"."_"."$prefix"."_"."$delly_program"."_window.out";
			
			open (IN_SIMPLIFIED, "$delly_simplified_out") || die "Cannot open $delly_simplified_out\n\n";
				@file_array = <IN_SIMPLIFIED>;
				$no_of_lines_simplified_out = scalar @file_array;
			close IN_SIMPLIFIED;
			
			open (OUT_WINDOW,">$delly_window_out")|| die "Cannot create output file: $delly_window_out";
			
			print OUT_WINDOW "FILE_NAME\tCHROMOSOME\tWINDOW_START\tNO_IN_WINDOW\n";
			
			print "\tGrouping file $list_count: $delly_simplified_out into windows of $window_size...\n";
			
			# Set up initial states #
			$window_start = $first_position;
			$window_end = $window_start + $window_size;
			$last_line_checked = 1;
			$last_line_checked_even = 1;
			$last_line_checked_odd = 1;
			
			#print "\tNo of lines: $no_of_lines_filtered_out\n\n";
			$even_check= 0; # first $even = "false";
			
			##########################################################
			# Move forward a window at a time (not overlapping)      #
			# (using file stored in @file_array)                     #
			#                                                        #
			# At each window position look down the whole of         #
			# the simplified file and see which lines fit within it. #
			##########################################################
			
			for ($position_count=$first_position;$position_count<=$last_position;$position_count = $position_count + ($window_size/2))
			{
				##########################################
				# Set start and end points of the window #
				##########################################
				$window_start = $position_count;
				$window_end = $position_count + $window_size;
				
				#print OUT_WINDOW "CHECKING WINDOW $window_start - $window_end\n";
				
				$window_count = 0;
				
				if (($position_count % 1000) == 0){print "Position: $position_count   Creating _windows.out file\n";}
				
				# Work out if this an even or an odd window #
				$even_check = $even_check + 1;
				if (($even_check % 2) == 0) {$even="true"} else {$even = "false"};
				
				if ($even eq "true")
				{
					$last_line_checked = $last_line_checked_even;
				}
				if ($even eq "false")
				{
					$last_line_checked = $last_line_checked_odd;
				}
						
				##################################################################################
				# Now loop through simplified file to find any positions that are in this window #
				# (There may be some windows for which there are no positions)                   #
				##################################################################################
				$something_in_window="false";
				
				for ($line_count = 1;$line_count<$no_of_lines_simplified_out;$line_count++)
				{
					$single_line = $file_array[$line_count];
					chomp $single_line;
				
					@item=split(/\t/,$single_line);
					$prefix = $item[0];
					$chromosome = $item[1];
					$position_start = $item[2];
				
					if ($even eq "true"){$last_line_checked_even = $line_count;}
					if ($even eq "false"){$last_line_checked_odd = $line_count;}
					
					if ($position_start < $window_end)
					{
						##########################################
						# Check if position is within the window #
						##########################################
						if (($position_start >= $window_start) && ($position_start < $window_end))
						{
							$window_count = $window_count + 1;
							
							$something_in_window="true";
							#$position_found = $position_start;
							
							if ($even eq "true")
							{
								$last_line_checked_even = $line_count;
							}
							if ($even eq "false")
							{
								$last_line_checked_odd = $line_count;
							}
							
							#print "File: $list_count\tLine: $line_count\tPosition: $position_start $window_start > $position_start > $window_end\t$window_count\n";
							
						}

					} # if position_start is < window_end
					
					#####################################################
					# If position is beyond the end of the window then  #
					# write to the output window file.  Then jump out   #
					# of the loop (no point going any further)          #
					#####################################################
					if ($position >= $window_end) ### SHOULD THESE BE POSITION_START??????????
					{
						goto END_LOOP;
					}
					
					
				} # $line_count loop
				
				END_LOOP: # jump out of loop to here. Past end of this window#
				
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


	if (($files_to_use eq "bam") || ($files_to_use eq "delly_out") || ($files_to_use eq "filtered_out") || ($files_to_use eq "simplified_out") || ($files_to_use eq "window_out"))
	{
		print "\n\n";
		print "#-----------------------------------------#\n";
		print "# Merging delly_out window files...   #\n";
		print "#-----------------------------------------#\n\n";
		
		######################################################
		# Stage 4                                            #
		# Load in files one at a time and add to @file_array #
		# to make the combined delly_all.out file            #
		######################################################
		
		print "Opening file $delly_all_out\n\n";
		
		open (OUT_ALL,">$delly_all_out")|| die "Cannot create output file: $delly_all_out";
		
		# Write Headers to first line #
		print OUT_ALL "CHR\tWINDOW_START";
		
		for ($list_count=1;$list_count<=$no_of_files;$list_count++)
		{
			#$bam_file = $bam_file_array[$list_count];
			$prefix = &get_prefix($bam_file_array[$list_count]);
			$delly_window_out = "$run_title"."_"."$prefix"."_"."$delly_program"."_window.out";
			
			#Column header for this file #
			print OUT_ALL "\t$prefix";
			
			# Read file into file_array #
			print "\tReading file $list_count: $delly_window_out\n";
			open (IN, "$delly_window_out") || die "Cannot open $delly_window_out\n\n";
				@file_array = <IN>;
				$no_of_lines_windows_file = scalar @file_array;
			close IN;
			
			# Work down file storing data into @line_array #
			for ($line_count = 1;$line_count <$no_of_lines_windows_file;$line_count++)
			{
				$single_line = $file_array[$line_count];
				chomp $single_line;
			
				@item=split(/\t/,$single_line);
				$prefix = $item[0];
				$chromosome = $item[1];
				$position = $item[2];
				$delly_count = $item[3];
			
				if ($list_count == 1)
				{
					$line_array[$line_count] = "$chromosome\t$position\t$delly_count";
				}
				if ($list_count > 1)
				{
					$line_array[$line_count] = "$line_array[$line_count]\t$delly_count";
				}
				
			}
			
		} # $list_count loop
		
		# Work down file writing data from @line_array to OUT_ALL #
		print "\n\t\t==> Writing combined file $delly_all_out...\n";
		print OUT_ALL "\n";
		
		for ($line_count = 1;$line_count<$no_of_lines_windows_file;$line_count++)
		{
			print OUT_ALL "$line_array[$line_count]\n";
		}
		
		close OUT_ALL;
		
	} # end if if $files_to_use eq etc (Stage 4)


	if (($files_to_use eq "bam") || ($files_to_use eq "delly_out") || ($files_to_use eq "filtered_out") || ($files_to_use eq "simplified_out") || ($files_to_use eq "window_out") || ($files_to_use eq "delly_all_out"))
	{
		################################################
		# Stage 5                                      #
		# Doing statistics and creating assoc.out file #
		################################################
		print "\n\n";
		

		print "#------------------------------------#\n";
		print "# Reading delly_all.out file... #\n";
		print "#                                    #\n";
		print "# Calculating delly statistics  #\n";
		print "# and creating delly_assoc.out  #\n";
		print "#------------------------------------#\n\n";
		# Read file into file_array #
		print "\tReading file $delly_all_out\n";
		open (IN_ALL, "$delly_all_out") || die "Cannot open $delly_all_out\n\n";
			@file_array = <IN_ALL>;
			$no_of_lines_all_file = scalar @file_array;
		close IN_ALL;
		
		print "\n\tNo. of lines in $delly_all_out:\t$no_of_lines_all_file\n\n";
		
		#######################################################
		# Open file for filtered delly data and statistics #
		#######################################################
		$output_assoc_file="$run_title"."_"."$delly_program"."_w"."$window_size"."_assoc.out";
		
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
		# Looping through file_array for delly_all.out #
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
					$affected_total_delly = 0;
					$normal_total_delly = 0;
					#$affected_total_delly_normalised = 0;
					#$normal_total_delly_normalised = 0;
					
					$delly_string = "";
					
					##################################################
					# Put together string of delly to show user #
					##################################################
					for ($block_count=1;$block_count <=$no_of_blocks;$block_count++)
					{
						$delly_array[$block_count] = $item[$block_count + 1];
						$delly_array[$block_count] = sprintf('%.0f', $delly_array[$block_count]); 

						$delly_string = "$delly_string"."\t"."$delly_array[$block_count]";
						
						if ($status_array[$block_count] eq "affected")
						{
							$affected_total_delly = $affected_total_delly + $delly_array[$block_count];
							#$affected_total_delly_normalised = $affected_total_delly_normalised + ($delly_array[$block_count] * $weighting_array[$block_count]);
						}
						
						if ($status_array[$block_count] eq "control")
						{
							$normal_total_delly = $normal_total_delly + $delly_array[$block_count];
							#$normal_total_delly_normalised = $normal_total_delly_normalised + ($delly_array[$block_count] * $weighting_array[$block_count]);
						}
					} # $block_count loop
					
					############################################################
					# Mean delly counts on each line of delly_all.out #
					############################################################
					if ($normalise_delly eq "false")
					{
						$mean_delly_affected = $affected_total_delly / $no_of_affecteds;
						$mean_delly_normal = $normal_total_delly / $no_of_normals;
						$mean_delly_both = ($affected_total_delly + $normal_total_delly) / ($no_of_affecteds + $no_of_normals);
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
					$delly_ratio = 0;

					for ($block_count=1;$block_count <=$no_of_blocks;$block_count++)
					{
						if ($status_array[$block_count] eq "affected")
						{
							$sum_squares_affected = $sum_squares_affected + (($delly_array[$block_count] - $mean_delly_affected) ** 2);
						}
						if ($status_array[$block_count] eq "control")
						{
							$sum_squares_normal = $sum_squares_normal + (($delly_array[$block_count] - $mean_delly_normal) ** 2);
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
						$t_test = ($mean_delly_affected - $mean_delly_normal) / (sqrt((($variance_affected / $no_of_affecteds) + ($variance_normal / $no_of_normals))));
						$t_test = abs($t_test); # t_statistic can be positve or negative
					}
					
					############################
					# Coefficient of variation #
					############################
					if ($mean_delly_affected > 0 ) {$covar_affected = $sd_affected / $mean_delly_affected;}
					if ($mean_delly_normal > 0) {$covar_normal = $sd_normal / $mean_delly_normal;}

					$covar_all_samples = ($covar_affected + $covar_normal ) / 2;
					
					
					####################
					# delly ratio #
					####################
					if ($mean_delly_normal > 0) {$delly_ratio = $mean_delly_affected / $mean_delly_normal;}


					################################
					# Format to two decimal places #
					################################ 
					$covar_affected = sprintf('%.2f', $covar_affected); 
					$covar_normal = sprintf('%.2f', $covar_normal); 
					$delly_ratio = sprintf('%.2f', $delly_ratio );
					$t_test = sprintf('%.2f', $t_test );
					$covar_all_samples = sprintf('%.2f', $covar_all_samples );
					$mean_delly_affected = sprintf('%.2f', $mean_delly_affected);
					$mean_delly_normal = sprintf('%.2f', $mean_delly_normal);
				
				
					######################################
					# Write to output ASSOC file         #
					# (if greater than T-test threshold) #
					######################################
					if ($t_test >= $t_test_threshold)
					{
						print OUT_ASSOC "$chromosome\t$position";
						if ($include_individual_columns eq "true"){print OUT_ASSOC "$delly_string";}
						print OUT_ASSOC "\t$mean_delly_affected\t$mean_delly_normal";
						print OUT_ASSOC "\t$covar_affected\t$covar_normal";
						print OUT_ASSOC "\t$delly_ratio\t$t_test";
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

} # passes loop
close COMMAND_LOG;
close OUT_ALL;

print "\n\n";
print "###############\n";
print "#  FINISHED   #\n";
print "###############\n\n";

print "Output files:\n\n";

print "File with all the delly data:\n\n";  
print "\t$delly_all_out\n\n\n";

print "delly association file, for plotting association graph: \n\n";   	   
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



