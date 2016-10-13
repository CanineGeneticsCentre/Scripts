#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	singletons_plot   						                            #     
#									                                    #
#	Looks at a batch of BAM files, some cases and some controls.        #
#									                                    #
#   It looks for singletons in the BAM files and records whether there  #
#   are significantly more in the cases than in the controls.           #
#									                                    #
#   It then produces an output which can be plotted like an association #
#   plot									                            #
#########################################################################

#############################
# Mike Boursnell July 2012  #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;
use Cwd;
my $version							= "6";


#Various Fixed Parameters 
my $default_window_size				= 1000;
my $include_individual_columns		= "true";
my $t_test_threshold				= 1; # Only T-test values over this go into the ASSOC file
my $show							= "false"; #for debugging

#new for singleton
my $window_size						= 0;
my $chromosome_to_use				= ""; # currently can't work on multiple chromosomes so you have to choose
my $files_to_use					= "";
my $singleton_bam_txt 				= "";
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
my $no_of_lines_singletons_bam_txt	= 0;
my $line_count						= 0;
my $position						= 0;
my $no_of_lines_simplified_singleton_out		= 0;
my $no_of_lines_windows_file		= 0;
my $no_of_lines_all_file			= 0;
my $input_string					= "";
my $first_position					= 0;
my $last_position					= 0;
my $window_start					= 0;
my $window_end						= 0;
my $window_count					= 0; # Maybe replace with singletons_count ???
my $last_line_checked				= 0;
my $last_line_checked_even			= 0;
my $last_line_checked_odd			= 0;
my $position_count					= 0;
my $singletons_count				= 0; # number of singletons in a window
my $something_in_window				= "false";
my $no_of_blocks					= 0;
my $block_count						= 0;
my $affected_total_singletons		= 0;
my $normal_total_singletons			= 0;
my $normalise_singletons			= "false";
my $bad_line_count					= 0;
my $even_check						= 0;
my $even							= "false";
my $start_row						= 0;
my $end_row							= 0;

#Statistics
my $mean_singletons_affected		= 0;
my $mean_singletons_normal			= 0;
my $mean_singletons_both			= 0;
my $sum_squares_affected			= 0;
my $sum_squares_normal 				= 0;
my $variance_affected				= 0;
my $variance_normal					= 0;
my $covar_affected					= 0;
my $covar_normal					= 0;
my $singletons_ratio				= 0;
my $sd_affected						= 0;
my $sd_normal						= 0;
my $t_test							= 0;
my $covar_all_samples				= 0;


# FIle names #
my $singleton_bam					= "";
my $simplified_singleton_out					= "";
my $singletons_all_out				= "";
my $singleton_window_out			= "";
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
my $singletons_string				= "";


my @bam_file_array					= ();
my @sample_name_array				= ();
my @status_array					= ();
my @item							= ();
my @file_array						= ();
my @simplified_file_array			= ();
my @line_array						= (); # for making merged file
my @singletons_array				= ();

######################################
# Create a new for the BAM directory #
######################################
$bam_files_directory = "$ENV{HOME}/bam_files";

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "#######################################\n";
print color 'bold white';
print "            singletons_plot             \n";
print color 'bold magenta';
print "#######################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This program analyses the unmapped pairs (singletons) in BAM files\n";
print "    (singletons may indicate an insertion relative to the reference sequence)\n\n";

print "    The program compares cases and controls and produces an output file \n";
print "    which can be plotted as an association plot.\n\n";

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

&print_message("Which do you want to do?","input");

print "   <1>  Use BAM files to start analysis from the beginning\n";
print "   <2>  Use '_singleton_bam.txt' files already made with samtools.\n";
#print "   <3>  Use 'simplified_singleton_out' simplified files already made by this program.\n";


$answer=<STDIN>;
chomp $answer;
if (substr($answer,0,1) eq "1" ){$files_to_use = "bam"}
if (substr($answer,0,1) eq "2" ){$files_to_use = "singleton_bam_txt"}
if (substr($answer,0,1) eq "3" ){$files_to_use = "simplified_singleton_out"}
if (substr($answer,0,1) eq "4" ){$files_to_use = "singleton_windows_out"}
if (substr($answer,0,1) eq "5" ){$files_to_use = "singletons_all_out"}

print "\nFiles to use: $files_to_use\n\n";

&print_message("Please enter a name for this analysis (with no spaces)","input");

$run_title = <STDIN>;
chomp $run_title;


&print_message("Which chromosome are you looking at?","input");

print "Enter chromosome number:  ";
$chromosome_to_use = <STDIN>;
chomp $chromosome_to_use;

if (index($chromosome_to_use,"chr") == -1){$chromosome_to_use = "chr"."$chromosome_to_use"}

print "Chromosome chosen: $chromosome_to_use\n\n";


##############
# File names #
##############
$singletons_all_out = "$run_title"."_singletons_all.out";

#############################################
# Make sure the list file is in Unix format #
#############################################

$command = "dos2unix $list_file";
print("\n$command\n");
system("$command");
print "\n";

#############################################
# Get window size for the plot              #
#############################################
if ($window_size == 0)
{
	print "Window size for the plot (default = $default_window_size):   ";
	$answer=<STDIN>;
	chomp $answer;
	if ($answer eq ""){$window_size=$default_window_size}else{$window_size = $answer}
	print "\n";
}
	
	
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
	$singleton_bam_txt = $prefix."_singleton_bam.txt";


	#################################################################
	# Check file exists in current directory or bam_files directory #
	#################################################################
	$bam_file_in_bam_files = "$bam_files_directory/$bam_file";
	$bam_file_array[$list_count]=$bam_file;
	$status_array[$list_count]=$status;
	

	#####################################
	# Check for BAM files (if required) #
	#####################################
	if ($files_to_use eq "bam")
	{
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

	} # if using BAM files

	###################################################
	# Check for singleton_bam_txt files (if required) #
	###################################################
	if ($files_to_use eq "singleton_bam_txt")
	{
		if (-e $singleton_bam_txt)
		{
			print "File $singleton_bam_txt exists in the current directory\n";
			#$bam_file_array[$list_count]=$bam_file;
		}

		if (!-e $singleton_bam_txt)
		{
			print "\nERROR: File $singleton_bam_txt not in the current directory\n\n";
			exit;
		}

	} # if using singleton_bam_txt files

	
} # end of while loop for reading list_file (file of file names of BAM files)

$no_of_files = $list_count;

if ($files_to_use eq "bam"){print "\nNo of BAM files: $no_of_files\n"}
if ($files_to_use eq "singleton_bam_txt"){print "\nNo of singleton_bam_txt files: $no_of_files\n"}

close LIST;


$command_log = "singletons_command_log.out";
	
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "COMMAND LOG for singletons_plot.pl $version\n\n";

print COMMAND_LOG "Running with files_to_use = $files_to_use\n\n";

print COMMAND_LOG "Input file:        \t$list_file\n\n";
print COMMAND_LOG "List of BAM files:\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print COMMAND_LOG "File $list_count	\t$bam_file_array[$list_count]\t\t$status_array[$list_count]\n";
}
print COMMAND_LOG "\n\n";


if ($files_to_use eq "bam")
{
	print "\n";

	&print_message("Using samtools to extract singleton reads from the BAM files...","message");


	###############################################
	# Stage 1                                     #
	# Create output files with singleton reads    #
	# First a new BAM file and then a text file   #
	# (looks for BAM files in directory bam_files #
	# if it can't be found in current directory)  #
	###############################################
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		$bam_file = $bam_file_array[$list_count];
		$prefix = &get_prefix($bam_file_array[$list_count]);
		$singleton_bam = $prefix."_singleton.bam";
		$singleton_bam_txt = $prefix."_singleton_bam.txt";
		$bam_file_in_bam_files = "$bam_files_directory/$bam_file";
		
		print "File: $bam_file\n";
		print "\nCreating $singleton_bam and $singleton_bam_txt\n";
		
		if (-e $bam_file)
		{
			&run_unix_command("samtools view -u -f 8 -F 260 $bam_file > $singleton_bam");
			&run_unix_command("samtools view $singleton_bam > $singleton_bam_txt");
		}
		
		if ((!-e $bam_file) && (-e $bam_file_in_bam_files))
		{
			&run_unix_command("samtools view -u -f 8 -F 260 $bam_file_in_bam_files > $singleton_bam");
			&run_unix_command("samtools view $singleton_bam > $singleton_bam_txt");
		}
		
		# Delete intermediate BAM file #
		&delete_file("$singleton_bam");
		
		print "\n";
	}
}

if (($files_to_use eq "bam") || ($files_to_use eq "singleton_bam_txt") )
{
	###############################################
	# Stage 2                                     #
	# Simplify these singleton text files so that #
	# they only contain FILENAME, CHR and POS     #
	###############################################
	
	print "\n\n";
	&print_message("Simplifying singleton_bam.txt files..","message");

	print COMMAND_LOG "\nSimplifying singleton_bam.txt files..\n";

	$first_position = 1000000000;
	$last_position = 0;
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		$bam_file = $bam_file_array[$list_count];
		$prefix = &get_prefix($bam_file_array[$list_count]);
		$singleton_bam_txt = $prefix."_singleton_bam.txt";
		$simplified_singleton_out = $prefix."_simplified_singleton.out";
		
		&print_both ("\tSimplifying file $list_count: $singleton_bam_txt...\n");
		
		open (IN_ALL, "$singleton_bam_txt") || die "Cannot open $singleton_bam_txt\n\n";
			@file_array = <IN_ALL>;
			$no_of_lines_simplified_singleton_out = scalar @file_array;
		close IN_ALL;
		
		open (OUT, ">$simplified_singleton_out")|| die "Cannot create output file: $simplified_singleton_out";
		
		print OUT "FILE\tCHR\tPOS\n";
		
		########################################################
		# Looping through file_array of singleton_bam_txt file #
		########################################################
		for ($line_count = 1;$line_count<=$no_of_lines_simplified_singleton_out - 1; $line_count++)
		{
			$single_line = $file_array[$line_count];
			chomp $single_line;
			
			@item=split(/\t/,$single_line);
			$array_size = scalar @item;

			$chromosome=$item[2];
			$position=$item[3];
			
			if ($position < $first_position){$first_position = $position};
			if ($position > $last_position){$last_position = $position};
			
			print OUT "$prefix\t$chromosome\t$position\n";	
		}
		#print "File $list_count\tFirst position:\t$first_position\tLast position: \t$last_position\n";

	} # $list_count loop to simplfy singleton text files
}
close OUT;

&print_message("Do you want to delete the singleton_bam_txt files as the simplified files have all the essential data?","input");

#$answer = <STDIN>;

print "\n\n";

############################################
# These are the overall lowest and highest #
# positions from all the input files       #
############################################
print "\tLowest position:  \t$first_position\n";
print "\tHighest position: \t$last_position\n";

if ($files_to_use eq "simplified_singleton_out")
{
	print "When you have chose option simplified_singleton_out then I'm not sure these work...\n\n";
	print "\tLowest position:  \t$first_position\n";
	print "\tHighest position: \t$last_position\n";

	$answer=<STDIN>;
}



if (($files_to_use eq "bam") || ($files_to_use eq "singleton_bam_txt") || ($files_to_use eq "simplified_singleton_out"))
{
	###################################################
	# Stage 3                                         #
	# Reducing the simplified singleton text file     #
	# to a file which contains a list of windows      #
	# and the number of singletons within each window #
	###################################################
	
	print "\n\n";
	&print_message("Grouping data from simplified_singleton_out files into windows...","message");

	
	if ($window_size == 0)
	{
		print "Window size (default = $default_window_size):   ";
		$answer=<STDIN>;
		chomp $answer;
		if ($answer eq ""){$window_size=$default_window_size}else{$window_size = $answer}
		print "\n";
	}
	
	##########################################
	# Work down the list of simplified files #
	##########################################
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		$prefix = &get_prefix($bam_file_array[$list_count]);
		$simplified_singleton_out = $prefix."_simplified_singleton.out";
		$singleton_window_out = $prefix."_window_singleton.out";
		
		#Read each file into an array
		open (IN_ALL, "$simplified_singleton_out") || die "Cannot open $simplified_singleton_out for input\n\n";
			@simplified_file_array = <IN_ALL>;
			$no_of_lines_simplified_singleton_out = scalar @simplified_file_array;
		close IN_ALL;
		
		open (OUT_WINDOW,">$singleton_window_out")|| die "Cannot create output window file: $singleton_window_out";
		
		print OUT_WINDOW "FILE\tCHR\tWINDOW_START\tNO_IN_WINDOW\n";
		
		print "\tGrouping file $list_count: $simplified_singleton_out into windows of $window_size...\n";
		
		$window_start = $first_position;
		$window_end = $window_start + $window_size;
		$last_line_checked = 1;
		$last_line_checked_even = 1;
		$last_line_checked_odd = 1;
		
		#print "\tNo of lines: $no_of_lines_simplified_singleton_out\n\n";
		
		####################################################################
		# Find first and last line of chosen chromosome in simplified file #
		####################################################################
		$start_row = 0;
		for ($line_count = 1; $line_count <= $no_of_lines_simplified_singleton_out; $line_count++)
		{
			$single_line = $simplified_file_array[$line_count];
			chomp $single_line;
			@item=split(/\t/,$single_line);
			$chromosome = $item[1];



			if ($chromosome eq $chromosome_to_use)
			{
				$start_row = $line_count;
				print "Chr $chromosome_to_use forund on line $start_row\n";
				$answer=<STDIN>;
			}
		}


		#####################################################
		# Move forward a window at a time (overlapping 50%) #
		# (using file in @simplified_file_array)            #
		#####################################################
		$even_check= 0; # first $even = "false";
		
		for ($position_count = $first_position; $position_count<=$last_position; $position_count = $position_count + ($window_size/2))
		{
			$even_check = $even_check + 1;
			if (($even_check % 2) == 0) {$even="true"} else {$even = "false"};
			
			
			$window_start = $position_count;
			$window_end = $position_count + $window_size;
			
			$something_in_window="false";
			
			if ($even eq "true")
			{
				$last_line_checked = $last_line_checked_even;
			}
			if ($even eq "false")
			{
				$last_line_checked = $last_line_checked_odd;
			}
			
			print "Position count: $position_count\n";
			

			


			##################################################################################
			# Now loop through simplified file to find any positions that are in this window #
			# (There may be some windows for which there are no positions)                   #
			# NOTE: chr has to match as well                                                 #
			##################################################################################
			for ($line_count = $last_line_checked; $line_count <= $no_of_lines_simplified_singleton_out; $line_count++)
			{
				$single_line = $simplified_file_array[$line_count];
				chomp $single_line;
			
				@item=split(/\t/,$single_line);
				$prefix = $item[0];
				$chromosome = $item[1];
				$position = $item[2];
				
				
				#############################################
				# Check chromosome is the chosen chromosome #
				#############################################
				if ($chromosome eq $chromosome_to_use)
				{
					##########################################
					# Check if position is within the window #
					##########################################
					if (($position >= $window_start) && ($position < $window_end))
					{
						$window_count = $window_count + 1;
						
						$something_in_window="true";
						
						$last_line_checked = $line_count;
						
						if ($even eq "true")
						{
							$last_line_checked_even = $line_count;
						}
						if ($even eq "false")
						{
							$last_line_checked_odd = $line_count;
						}
					}
				
					# DEBUGGING
					if ($line_count > 4247)
					{
						print "Line:           \t$line_count\n\n";
						print "prefix:         \t$prefix\n";
						print "chromosome:     \t$chromosome\n";
						print "position:       \t$position\n";

						print "window_start:   \t$window_start\n";
						print "window_end:     \t$window_end\n";
						
						$answer=<STDIN>;
					}

				} # if chr is chromosome_to_use


				#####################################################
				# If position is beyond the end of the window then  #
				# write to the output window file.  Then jump out   #
				# of the loop (no point going any further)          #
				#####################################################
				if ($position >= $window_end)
				{
					goto END_LOOP;
				}
				
			} # $line_count loop
			
			END_LOOP: # jump out of loop to here #
			
			if ($something_in_window eq "true")
			{
				print OUT_WINDOW "$prefix\t$chromosome\t$window_start\t$window_count\n";

				print "$prefix\t$chromosome\t$window_start\t$window_count\n";
				$answer=<STDIn>;
			}
			
			if ($something_in_window eq "false")
			{
				print OUT_WINDOW "$last_prefix\t$last_chromosome\t$window_start\t0\n";
			}
			$window_count = 0;
			$last_line_checked = $line_count;
			$last_prefix = $prefix;
			$last_chromosome = $chromosome;

		} # position_count loop for moving along windows
		
		close OUT_WINDOW;
		
	} # $list_count loop	
}


if (($files_to_use eq "bam") || ($files_to_use eq "singleton_bam_txt") || ($files_to_use eq "simplified_singleton_out") || ($files_to_use eq "singleton_windows_out"))
{
	print "\n\n";
	print "#-----------------------------------------#\n";
	print "# Merging singleton_out window files...   #\n";
	print "#-----------------------------------------#\n\n";
	
	######################################################
	# Stage 4                                            #
	# Load in files one at a time and add to @file_array #
	# to make the combined singleton_all.out file        #
	######################################################
	
	open (OUT_ALL,">$singletons_all_out")|| die "Cannot create output file: $singletons_all_out";
	
	# Write Headers to first line #
	print OUT_ALL "CHR\tWINDOW_START";
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		#$bam_file = $bam_file_array[$list_count];
		$prefix = &get_prefix($bam_file_array[$list_count]);
		$singleton_window_out = $prefix."_window_singleton.out";
		
		#Column header for this file #
		print OUT_ALL "\t$prefix";
		
		# Read file into file_array #
		print "\tReading file $list_count: $singleton_window_out\n";
		open (WINDOW_IN, "$singleton_window_out") || die "Cannot open $singleton_window_out\n\n";
			@file_array = <WINDOW_IN>;
			$no_of_lines_windows_file = scalar @file_array;
		close WINDOW_IN;
		
		# Work down file storing data into @line_array #
		for ($line_count = 1;$line_count <$no_of_lines_windows_file;$line_count++)
		{
			$single_line = $file_array[$line_count];
			chomp $single_line;
		
			@item=split(/\t/,$single_line);
			$prefix = $item[0];
			$chromosome = $item[1];
			$position = $item[2];
			$singletons_count = $item[3];
		
			if ($list_count == 1)
			{
				$line_array[$line_count] = "$chromosome\t$position\t$singletons_count";
			}
			if ($list_count > 1)
			{
				$line_array[$line_count] = "$line_array[$line_count]\t$singletons_count";
			}
			
		}
		
	} # $list_count loop
	
	# Work down file writing data from @line_array to OUT_ALL #
	print "\n\t\t==> Writing combined file $singletons_all_out...\n";
	print OUT_ALL "\n";
	
	for ($line_count = 1;$line_count<$no_of_lines_windows_file;$line_count++)
	{
		print OUT_ALL "$line_array[$line_count]\n";
	}
	
} # end if if $files_to_use eq etc (Stage 4)



if (($files_to_use eq "bam") || ($files_to_use eq "singleton_bam_txt") || ($files_to_use eq "simplified_singleton_out") || ($files_to_use eq "singleton_windows_out") || ($files_to_use eq "singletons_all_out"))
{
	################################################
	# Stage 5                                      #
	# Doing statistics and creating assoc.out file #
	################################################
	print "\n\n";
	

	print "#------------------------------------#\n";
	print "# Reading singletons_all.out file... #\n";
	print "#                                    #\n";
	print "# Calculating singletons statistics  #\n";
	print "# and creating singletons_assoc.out  #\n";
	print "#------------------------------------#\n\n";
	
	#########################################################
	# Read whole of singletons_all_out file into file_array #
	#########################################################
	print "\tReading file $singletons_all_out\n";
	open (IN_ALL, "$singletons_all_out") || die "Cannot open $singletons_all_out\n\n";
		@file_array = <IN_ALL>;
		$no_of_lines_all_file = scalar @file_array;
	close IN_ALL;
	
	print "\n\tNo. of lines in $singletons_all_out:\t$no_of_lines_all_file\n\n";
	
	#########################################################
	# Open file for filtered singletons data and statistics #
	#########################################################
	$output_assoc_file="$run_title"."_singletons_w"."$window_size"."_assoc.out";
	
	open (OUT_ASSOC, ">$output_assoc_file") || die "Cannot open $output_assoc_file\n\n";
	
	###########
	# Headers #
	###########
	print OUT_ASSOC "CHR\tPOS\tTYPE";
	
	if ($include_individual_columns eq "true")
	{
		for ($list_count=1;$list_count <=$no_of_files;$list_count++)
		{
			$prefix = &get_prefix($bam_file_array[$list_count]);
			print OUT_ASSOC "\t$prefix";
		}
	} # if you want indivial DP data columns ($include_individual_columns eq "true")


	print OUT_ASSOC "\tAV_A\tAV_N\tCOV_A\tCOV_N\tRATIO\tCHANGE\tT_TEST\n";

	#print OUT_ASSOC "$chromosome\t$position";
	#if ($include_individual_columns eq "true"){print OUT_ASSOC "$singletons_string";}
	#print OUT_ASSOC "\t$mean_singletons_affected\t$mean_singletons_normal\t$singletons_ratio";
	#print OUT_ASSOC "\t$covar_affected\t$covar_normal\t$t_test";
	#print OUT_ASSOC "\n";
				
	#####################################################
	# Looping through file_array for singletons_all.out #
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
			
			if (($line_count % 10000) == 0){print "\tCreating assoc.out file.  Line: $line_count  \tPosition: $position\n";}

			$no_of_blocks = $array_size - 2;
			if ($no_of_blocks == $no_of_files)
			{
				$affected_total_singletons = 0;
				$normal_total_singletons = 0;
				#$affected_total_singletons_normalised = 0;
				#$normal_total_singletons_normalised = 0;
				
				$singletons_string = "";
				
				##################################################
				# Put together string of singletons to show user #
				##################################################
				for ($block_count=1;$block_count <=$no_of_blocks;$block_count++)
				{
					$singletons_array[$block_count] = $item[$block_count + 1];
					$singletons_array[$block_count] = sprintf('%.0f', $singletons_array[$block_count]); 

					$singletons_string = "$singletons_string"."\t"."$singletons_array[$block_count]";
					
					if ($status_array[$block_count] eq "affected")
					{
						$affected_total_singletons = $affected_total_singletons + $singletons_array[$block_count];
						#$affected_total_singletons_normalised = $affected_total_singletons_normalised + ($singletons_array[$block_count] * $weighting_array[$block_count]);
					}
					
					if ($status_array[$block_count] eq "control")
					{
						$normal_total_singletons = $normal_total_singletons + $singletons_array[$block_count];
						#$normal_total_singletons_normalised = $normal_total_singletons_normalised + ($singletons_array[$block_count] * $weighting_array[$block_count]);
					}
				} # $block_count loop
				
				############################################################
				# Mean singleton counts on each line of singletons_all.out #
				############################################################
				if ($normalise_singletons eq "false")
				{
					$mean_singletons_affected = $affected_total_singletons / $no_of_affecteds;
					$mean_singletons_normal = $normal_total_singletons / $no_of_normals;
					$mean_singletons_both = ($affected_total_singletons + $normal_total_singletons) / ($no_of_affecteds + $no_of_normals);
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
				$singletons_ratio = 0;

				
				############################
				# Calculate sun of squares #
				############################
				for ($block_count=1;$block_count <=$no_of_blocks;$block_count++)
				{
					if ($status_array[$block_count] eq "affected")
					{
						$sum_squares_affected = $sum_squares_affected + (($singletons_array[$block_count] - $mean_singletons_affected) ** 2);
					}
					if ($status_array[$block_count] eq "control")
					{
						$sum_squares_normal = $sum_squares_normal + (($singletons_array[$block_count] - $mean_singletons_normal) ** 2);
					}
				}
				
				############################
				# Calculate Variance       #
				############################
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
					$t_test = ($mean_singletons_affected - $mean_singletons_normal) / (sqrt((($variance_affected / $no_of_affecteds) + ($variance_normal / $no_of_normals))));
					$t_test = abs($t_test); # t_statistic can be positve or negative
				}
				
				
				############################
				# Coefficient of variation #
				############################
				if ($mean_singletons_affected > 0 ) {$covar_affected = $sd_affected / $mean_singletons_affected;}
				if ($mean_singletons_normal > 0) {$covar_normal = $sd_normal / $mean_singletons_normal;}

				$covar_all_samples = ($covar_affected + $covar_normal ) / 2;
				
				
				####################
				# Singletons ratio #
				####################
				if ($mean_singletons_normal > 0) {$singletons_ratio = $mean_singletons_affected / $mean_singletons_normal;}


				################################
				# Format to two decimal places #
				################################ 
				$covar_affected = sprintf('%.2f', $covar_affected); 
				$covar_normal = sprintf('%.2f', $covar_normal); 
				$singletons_ratio = sprintf('%.2f', $singletons_ratio );
				$t_test = sprintf('%.2f', $t_test );
				$covar_all_samples = sprintf('%.2f', $covar_all_samples );
				$mean_singletons_affected = sprintf('%.2f', $mean_singletons_affected);
				$mean_singletons_normal = sprintf('%.2f', $mean_singletons_normal);
			
			
				######################################
				# Write to output ASSOC file         #
				# (if greater than T-test threshold) #
				######################################
				if ($t_test >= $t_test_threshold)
				{
					print OUT_ASSOC "$chromosome\t$position\tsingletons";
					
					if ($include_individual_columns eq "true"){print OUT_ASSOC "$singletons_string";}
					
					print OUT_ASSOC "\t$mean_singletons_affected\t$mean_singletons_normal";
					print OUT_ASSOC "\t$covar_affected\t$covar_normal";
					print OUT_ASSOC "\t$singletons_ratio";
					print OUT_ASSOC "\tNA"; # Nothing to go in CHANGE column
					print OUT_ASSOC "\t$t_test";
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

print "File with all the singletons data:\n\n";  
print "\t$singletons_all_out\n\n\n";

print "Singleton association file, for plotting association graph: \n\n";   	   
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

####################################################
# Subroutine to print to screen and to COMMAND_LOG #
####################################################

sub print_both
{
	my $message = $_[0];

	print "$message";
	print COMMAND_LOG "$message";
}
