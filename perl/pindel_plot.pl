#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	pindel_plot      						                        	#     
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
use warnings;

use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;
use Cwd;
my $version							= "5";


#Various Fixed Parameters 
my $s_value							= 8;
my $default_window_size				= 500;
my $include_individual_columns		= "true";
my $t_test_threshold				= 2; # Only T-test values over this go into the ASSOC file
my $show							= "false"; #for debugging
my $proportion_supporting_threshold	= 0.25;


my $ref								= "/home/genetics/canfam2/canfam2.fasta";
my $deletion_size_threshold        	= 100; # for awk filtering
my $supporting_reads_threshold		= 3;   # for awk filtering
my $mapping_quality_threshold       = 5;  # for awk filtering

#
my $region							= "";
my $region_start					= 0;
my $region_end						= 0;
my $region_OK						= "false";
		
#new for pindel_plot
my $proportion_supporting			= 0;
my $total_supporting				= 0;
my $total_reads						= 0;
my $no_reads						= 0;
my $short_line						= "";
my $read_line_count					= 0;
my $read_1							= "";
my $read_2							= "";
my $start_position					= 0;
my $ref_seq							= "";
my $hash_line						= 0;
my $data_line						= 0;
my $ref_line						= 0;
my $first_read_line					= 0;
my $sv_type							= "";
my $sv_length						= 0;
my $sv_start						= 0;
my $expected_data_array_size		= 0; # (This is no_of_files * 7 ) + 31
my $no_supporting					= 0;
my $pass							= 0;
my $file_type						= ""; # LI, INV etc
my $sv_end							= "";
my $check_count						= 0;

my $cigar							= "";
my $position_start					= 0;
my $position_end					= 0;
my $deletion_size					= 0;
my $supporting						= 0; # No supporting reads
my $mapping_quality					= 0;
my $deletion_id						= "";
my $window_size						= 0;
my $files_to_use					= "";
my $deletion_bam_txt 				= "";
my $prefix							= "";
my $last_prefix						= "";
my $last_chromosome					= "";
my $command_log						= "";
my $single_line						= "";
my $array_size						= "";
my $list_count						= 0;
my $file_count						= 0;
my $no_of_files						= 0;
my $status							= "";
my $no_of_affecteds					= 0;
my $no_of_normals					= 0;
my $second_column_found				= "false";
my $no_of_lines_deletion_bam_txt	= 0;
my $line_count						= 0;
my $position						= 0;
my $no_of_lines_filtered_out		= 0;
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
my $deletion_count					= 0; # number of deletions in a window
my $something_in_window				= "false";
my $no_of_blocks					= 0;
my $block_count						= 0;
my $affected_total_supports			= 0;
my $normal_total_supports			= 0;
my $normalise_supports				= "false";
my $bad_line_count					= 0;
my $even_check						= 0;
my $even							= "false";

#Statistics
my $mean_supports_affected			= 0;
my $mean_supports_normal			= 0;
my $mean_supports_both				= 0;
my $sum_squares_affected			= 0;
my $sum_squares_normal 				= 0;
my $variance_affected				= 0;
my $variance_normal					= 0;
my $covar_affected					= 0;
my $covar_normal					= 0;
my $supports_ratio					= 0;
my $sd_affected						= 0;
my $sd_normal						= 0;
my $t_test							= 0;
my $covar_all_samples				= 0;
my $expected_affected				= 0;
my $expected_normal					= 0;
my $chi_squared						= 0;

# FIle names #
my $pindel_D						= "";
my $config_file						= ""; # pindel config file
my $pindel_out						= ""; # pindel output
my $pindel_all_out					= ""; # Out put from this program after pindel has been run
my $sam_file						= "";
my $deletions_out					= ""; # deletion output file
my $deletion_simplified_out			= ""; # simplified version of $deletions_out
my $deletion_filtered_out			= ""; # filtered version of deletion_out
my $file_to_process					= ""; # pindel out file to process (e.g. _D, _INV, _LI etc
my $pindel_concat_out				= ""; # Final output file, concatenated from all the _all_out files

my $deletion_bam					= "";
my $deletion_all_out				= "";
my $deletion_window_out				= "";
my $list_file						= "";
my $command							= "";
my $bam_file						= "";
my $bai_file						= "";
my $bam_bai_file					= "";
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
my $deletion_string					= "";


my @bam_file_array					= ();
my @bam_file_name_array				= ();
my @status_array					= ();
my @item							= ();
my @file_array						= ();
my @line_array						= (); # for making merged file
my @supports_array					= ();
my @pindel_column_name_array		= ();
my @unique_upstream_array			= ();
my @unique_downstream_array			= ();
my @total_upstream_array			= ();
my @total_downstream_array			= ();

				

#######################################
# Create a name for the BAM directory #
#######################################
$bam_files_directory = "$ENV{HOME}/bam_files";

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "#######################################\n";
print color 'bold white';
print "            pindel_plot                \n";
print color 'bold magenta';
print "#######################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This program uses the output files from the program pindel\n";
print "    and compares the results across cases and controls,\n";
print "    providing data for a Structural Variation (SV) association plot.\n\n";

print "    This version allows the user to choose a threshold for \n";
print "    the required proportion of reads supporting any SV.\n\n";

print "    Any SV with less than this threshold does not go into the output files\n\n";

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
print "   <1>  Use BAM files to start pindel analysis from the beginning\n";
print "   <2>  Use pindel output files such as _D (deletion file) files\n";


$answer=<STDIN>;
chomp $answer;
if (substr($answer,0,1) eq "1" ){$files_to_use = "bam"}
if (substr($answer,0,1) eq "2" ){$files_to_use = "output"}


print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Please enter a name for this analysis (with no spaces) \n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
$run_title = <STDIN>;
chomp $run_title;

##############
# File names #
##############
$deletion_all_out = "$run_title"."_supports_all.out";
$config_file = $run_title."_config.txt";
$pindel_out = $run_title."_pindel";

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
	$bam_file_name_array[$list_count] = $prefix;
	
	$bai_file = $prefix.".bai";
	$bam_bai_file = $prefix.".bam.bai";
	
	
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
		#print "File $bam_file exists in the current directory\n";
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

	####################################################
	# Check that bai file exists in for name.bam.bai   #
	# If they are in form name.bai make a renamed copy #
	####################################################
	if (! -e $bam_bai_file)
	{
		print "The index file $bam_bai_file cannot be found\n\n";
		
		if (-e $bai_file)
		{
			print "\tCopying $bai_file to $bam_bai_file";
			system("cp $bai_file $bam_bai_file");
		}
		if (! -e $bai_file)
		{
			print "\n\n";
			print "#########################\n";
			print "#  INDEX FILE MISSING   #\n";
			print "#########################\n\n";
			
			print "Can't find an index file for the BAM file $bam_file\n\n";
			print "This should be $bai_file or $bam_bai_file\n\n";
			
			print "Create an Index file using Samtools, and try again.\n\n\n";
			
			exit;
		}
	}
	
} # end of while loop for reading list_file (file of file names of BAM files)

$no_of_files = $list_count;

print "\nNo of BAM files: $no_of_files\n\n\n";
close LIST;


######################################################## TEMP!!!!!!!!!!!

if ($files_to_use eq "deletions")
{
 goto stage3;
}


print "\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " Which reference sequence do you want to use?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   Enter 1 for CanFam3\n";
print "   Enter 2 for CanFam2\n";
print "   Enter 3 for EquCab2\n";
print "   Enter 4 for Human\n\n";

print "   Enter 5 for Strep. equi\n";
print "   Enter 6 for Strep. zoo\n\n";

print "   Enter 9 for other\n\n";


$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3"; $species = "canis_familiaris"}
if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/canfam2/canfam2.fasta"; $ref_seq_name = "canfam2"; $species = "canis_familiaris"}
if (substr($answer,0,1) eq "3" ){$ref = "/home/genetics/equcab2/equcab2.fasta"; $ref_seq_name = "equcab2"; $species = "equus_caballus"}
if (substr($answer,0,1) eq "4" ){$ref = "/home/genetics/human/human.fasta"; $ref_seq_name = "human"; $species = "homo_sapiens"}

if (substr($answer,0,1) eq "5" ){$ref = "/home/genetics/strep_equi/strep_equi.fasta"; $ref_seq_name = "s_equi"; $species = "streptococcus_equi"}
if (substr($answer,0,1) eq "6" ){$ref = "/home/genetics/strep_zoo/strep_zoo.fasta"; $ref_seq_name = "s_zoo"; $species = "streptococcus_zoo"}

print "\n\n";




print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Please define your region of interest (eg 'chr5:21000000-23000000' or just 'chr5')\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

until ($region_OK eq "true")
{
	print "> ";
	$region = <STDIN>;
	chomp $region;

	######################################
	# Check for presence of chr, : and - #
	######################################
	if ((index($region,":") == -1) || (index($region,"chr") == -1) || (index($region,"-") == -1))
	{
		print "\nYou must give the region in this format: chr15:53000000-55000000\n\n";
		$region_OK = "false";
	}

	#################################################
	# If full region given chr15:34000000-390000000 #
	#################################################
	if ((index($region,"chr") > -1) && (index($region,":") > -1) && (index($region,"-") > -1))
	{
		$chromosome = substr($region,3,index($region,":")-3);
		$region_start = substr($region,index($region,":")+1,index($region,"-") - index($region,":") - 1 );
		$region_end = substr($region,index($region,"-")+1,99);
		$region_OK = "true";
	}

}
print "Region:        \t$region\n";
print "Chromosome:    \t$chromosome\n";
print "Region start:  \t$region_start\n";
print "Region end:    \t$region_end\n\n";

print "\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " What proportion of reads must support the SV? (default = $proportion_supporting_threshold)\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

$answer = <STDIN>;
chomp $answer;

if ($answer ne "")
{
	$proportion_supporting_threshold = $answer;
}

if ($files_to_use eq "bam")
{


	#############################
	# Stage 1                   #
	# Create pindel config file #
	#############################
	open (CONFIG, ">$config_file") || die "Cannot create output file: $config_file";

	print "\nCreating pindel CONFIG file...\n\n";

	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		$bam_file = $bam_file_array[$list_count];
		$prefix = &get_prefix($bam_file);
		
		$bam_file_in_bam_files = "$bam_files_directory/$bam_file";
		if (-e $bam_file)
		{
			print CONFIG "$bam_file\t200\t$prefix\n";
			print "\t$bam_file\t200\t$prefix\n";
		}

		if ((!-e $bam_file) && (-e $bam_file_in_bam_files))
		{
			print CONFIG "$bam_file_in_bam_files\t200\t$prefix\n";
			print "\t$bam_file\t200\t$prefix\n";
		}
	}
	close CONFIG;

} # if files_to_use eq "bam"

$command_log = "$run_title"."_pindel_plot_command_log.out";
	
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "COMMAND LOG for pindel_plot.pl $version\n\n";

print COMMAND_LOG "Running with files_to_use = $files_to_use\n\n";

print COMMAND_LOG "Input file:        \t$list_file\n\n";
print COMMAND_LOG "List of BAM files:\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print COMMAND_LOG "File $list_count	\t$bam_file_array[$list_count]\t\t$status_array[$list_count]\n";
}
print COMMAND_LOG "\n\n";


###########################################
# Stage 2                                 #
# Running pindel on the list of BAM files #
###########################################
if ($files_to_use eq "bam")
{
	print "\n\n";
	print "#---------------------------------------------#\n";
	print "# Running pindel using $config_file	...   \n";
	print "#---------------------------------------------#\n\n";

	#####################################################
	# Stage 1                                           #
	# pindel makes several output files #                   
	#####################################################

	&run_unix_command("/opt/pindel/pindel -f $ref -i $config_file -c $region -o $pindel_out");
}



###########################################################
# Stage 3: Parse pindel output file and for each feature  #
# check which which samples contribute to that feature    #
###########################################################
stage3:



if (($files_to_use eq "bam") || ($files_to_use eq "output"))
{
	print "\n\n";
	print "#---------------------------------------------#\n";
	print "# Processing pindel output files              #\n";
	print "#---------------------------------------------#\n\n";
	
	#########################################################
	# Open grand concatenated output file with all SV types #
	#########################################################
	$pindel_concat_out = "$run_title"."_pindel_concat_all.out";
	
	open (OUT_CONCAT,">$pindel_concat_out")|| die "Cannot create output file: $pindel_concat_out";
	
	#####################################
	# headers for pindel_concat_all_out #
	#####################################
	print OUT_CONCAT "POS\tTYPE";
	for ($file_count = 1; $file_count <=$no_of_files;$file_count++)
	{
		$prefix = &get_prefix($bam_file_array[$file_count]);
		print OUT_CONCAT "\t$prefix";
	}
				
	print OUT_CONCAT "\tT_test\n";
			
	
	###########################################################
	# Each pass processes a different pindel output file type #
	###########################################################
	print "\n";
	for ($pass = 1;$pass<=5;$pass++)
	{
		if ($pass ==1){$file_type="D"}
		if ($pass ==2){$file_type="SI"}
		if ($pass ==3){$file_type="INV"}
		if ($pass ==4){$file_type="TD"}
		if ($pass ==5){$file_type="LI"}

		if ($file_type ne "LI")
		{
			$expected_data_array_size = ($no_of_files * 7) + 31;
		}
		if ($file_type eq "LI")
		{
			$expected_data_array_size = ($no_of_files * 5) + 10;
		}
		
		###################################
		# Make names for in and out files #
		###################################
		$file_to_process = "$run_title"."_pindel_"."$file_type";
		$pindel_all_out = "$run_title"."_pindel_"."$file_type"."_all.out";
		
		print "Processing pindel output file: $file_to_process \t--> \t$pindel_all_out\n";
		
		if (!-e "$file_to_process"){print "\nFile $file_to_process not found\n\n";}
	
		if (-e $file_to_process)
		{
			
			#print "\nReading pindel file $file_to_process\n\n";
			
			open (PINDEL_OUT, "$file_to_process") || die "Cannot open $file_to_process\n\n";
				@file_array = <PINDEL_OUT>;
				$no_of_lines_deletions_out = scalar @file_array;
			close PINDEL_OUT;
			
			open (OUT_ALL,">$pindel_all_out")|| die "Cannot create output file: $pindel_all_out";
			
			
			#############################################################
			# headers for pindel_all_out                                #
			# These are in the order in the original input file         #
			# (i.e. probably affecteds followed by controls).           #
			# This is not necessarily the same as the alphabetic order  #
			# that pindel uses for output, so our output has to correct #
			# for that.                                                 #
			#############################################################
			print OUT_ALL "POS\tTYPE";
			for ($file_count = 1; $file_count <=$no_of_files;$file_count++)
			{
				$prefix = &get_prefix($bam_file_array[$file_count]);
				print OUT_ALL "\t$prefix";
			}
			print OUT_ALL "\tT_test\n"	;
			
			
			$first_read_line = 1000000;
			
			for ($line_count = 0;$line_count < $no_of_lines_deletions_out;$line_count++)
			{
				$single_line = $file_array[$line_count];
				chomp $single_line;
			
				$short_line = substr($single_line,0,110);
				
				######################################
				# This reads the line of hashes and  #
				# then sets up which lines are which #
				######################################
				if (substr($single_line,1,5) eq "#####")
				{
					$hash_line = $line_count;
					$data_line = $line_count + 1;
					$ref_line = $line_count + 2;
					$first_read_line = $line_count + 3;
				}
				
				############################
				# This reads the data line #
				############################
				if ($line_count == $data_line)
				{
					@item=split(/\s+/,$single_line);
					$array_size = scalar @item;
					
					if ($file_type ne "LI")
					{
						$sv_type = $item[1];
						$sv_length = $item[2];
						$sv_start = $item[9];
						$sv_end = $item[10];
					}
					
					# LI file has a different format #
					if ($file_type eq "LI")
					{
						$sv_type = $item[1];
						$sv_start = $item[4];
						$sv_end = $item[7];
						$sv_length = $sv_end - $sv_start;
					}
							
					if ($array_size != $expected_data_array_size)
					{
						print "\t\tLine: $line_count  Array size is $array_size NOT $expected_data_array_size\n\n";
					}
					
					
					
					$affected_total_supports = 0;
					$normal_total_supports = 0;
				
					#######################################################################
					# Loop through the files in the original input file of BAM file names #
					# to count the total number of supporting reads (up and downstream)   #
					#######################################################################
					for ($file_count = 1; $file_count <=$no_of_files;$file_count++)
					{
						if ($file_type ne "LI")
						{
							$pindel_column_name_array[$file_count] = $item[24 + ($file_count * 7)    ];
							$unique_upstream_array[$file_count]    = $item[24 + ($file_count * 7) + 4];
							$unique_downstream_array[$file_count]  = $item[24 + ($file_count * 7) + 6];
							
							$total_upstream_array[$file_count]    = $item[24 + ($file_count * 7) + 1];
							$total_downstream_array[$file_count]  = $item[24 + ($file_count * 7) + 2];
							
						}
						
						if ($file_type eq "LI")
						{
							$pindel_column_name_array[$file_count] = $item[5 + ($file_count * 5)    ];
							$unique_upstream_array[$file_count]    = $item[5 + ($file_count * 5) + 2];
							$unique_downstream_array[$file_count]  = $item[5 + ($file_count * 5) + 4];
							
							#########################################################
							# No info for total reads in LI files so use unique     #
							# This will give 100% proportion and allow all through  #
							#########################################################
							$total_upstream_array[$file_count]    = $unique_upstream_array[$file_count];
							$total_downstream_array[$file_count]  = $unique_downstream_array[$file_count];
						}

						
						#print "Sample data. File $file_count\t$sv_start\tType: $file_type\tExp Array: $expected_data_array_size\n";
						#print "name:       \t$pindel_column_name_array[$file_count]\n";
						#print "upstream:   \t$unique_upstream_array[$file_count]\n";
						#print "downstream: \t$unique_downstream_array[$file_count]\n";
						
						#$answer=<STDIN>;
						$no_supporting = $unique_upstream_array[$file_count] + $unique_downstream_array[$file_count];
						$no_reads = $total_upstream_array[$file_count] + $total_downstream_array[$file_count];
						
						#############################
						# Totals across all samples #
						#############################					
						$total_supporting = $total_supporting + $no_supporting;
						$total_reads = $total_reads + $no_reads;
						
						
						# Store in array for later sum of squares calculation #
						$supports_array[$file_count] = $no_supporting;
						
						
						if ($status_array[$file_count] eq "affected")
						{
							$affected_total_supports = $affected_total_supports + $no_supporting;
						}
						if ($status_array[$file_count] eq "control")
						{
							$normal_total_supports = $normal_total_supports + $no_supporting;
						}
					
					} # $file_count loop 2
					
					##################################################################
					# Calculate proportion of supporting to reads across all samples #
					##################################################################
					if (($total_reads + $total_supporting) > 0 )
					{
						$proportion_supporting = $total_supporting / ($total_reads + $total_supporting);
					}
					if (($total_reads + $total_supporting) == 0 )
					{
						$proportion_supporting = 0;
					}
					
					
					######################
					# Set totals to zero #
					######################
					$total_supporting = 0;
					$total_reads = 0;
					
					##########################################################
					# Again only write to output file if it is within region #
					# and if the proportion of supporting reads is higher    #
					# than the threshold                                     #
					##########################################################
					
					#print "$sv_start\tProportion proportion_supporting: $proportion_supporting\t(threshold is $proportion_supporting_threshold)\n";
						#$answer=<STDIN>;
						
					if (($sv_start >= $region_start ) && ( $sv_start <= $region_end))
					{
						
						
						if ($proportion_supporting > $proportion_supporting_threshold)
						{
							
							#####################################
							# First write the first two columns #
							#####################################
							
							print OUT_ALL "$sv_start\t$sv_type";
							print OUT_CONCAT "$sv_start\t$sv_type";
							
						
							###################################################
							# Then write in the SAME ORDER as the input files #
							# (default order is alphanumeric!!)               #
							###################################################
							for ($file_count = 1; $file_count <=$no_of_files;$file_count++)
							{
							
								#########################################################
								# Check through all four columns to find the one        #
								# that matches the BAM file name in bam_file_name_array #
								# then write the data FOR THAT COLUMN to output files   #
								#########################################################
								for ($check_count = 1; $check_count <=$no_of_files;$check_count++)
								{
									
									#$answer = <STDIN>;
									
									if ($bam_file_name_array[$file_count] eq $pindel_column_name_array[$check_count])
									{
										#print "Position $sv_start\n\n";
										print OUT_ALL "\t$supports_array[$check_count]";
										print OUT_CONCAT "\t$supports_array[$check_count]";
										
										
										#print "File count: $file_count/$no_of_files \tCheck_count: $check_count/$no_of_files\n";
										#print "Sample name array $file_count: $bam_file_name_array[$file_count]";
										#print "\tFile name array $check_count: $pindel_column_name_array[$check_count]\n\n";
										#print "Supports_array(check_count): $supports_array[$check_count]\n";
										#print "Supports_array(file_count): $supports_array[$file_count]\n";
									
										#print "PAUSE\n";
										#$answer=<STDIN>;
										
										last; # exit loop
										
									}
								} # check_count loop
		
							} # file_count loop
						
						} # if proportion is high enough
						
					} # if in region
					
					#################################################################
					# Do statistics (but only if there are enough supporting reads) #
					#################################################################
					if ($proportion_supporting > $proportion_supporting_threshold)
					{
						############################################################
						# Mean support counts on each line of pindel_all.out       #
						############################################################
						if ($normalise_supports eq "false")
						{
							$mean_supports_affected = $affected_total_supports / $no_of_affecteds;
							$mean_supports_normal = $normal_total_supports / $no_of_normals;
							$mean_supports_both = ($affected_total_supports + $normal_total_supports) / ($no_of_affecteds + $no_of_normals);
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
						$supports_ratio = 0;
							
						##################################################
						# Calculate sum of squares across all samples    #
						# (use $block_count as $file_count is being used #
						##################################################
						for ($block_count=1;$block_count <=$no_of_files;$block_count++)
						{
							if ($status_array[$block_count] eq "affected")
							{
								$sum_squares_affected = $sum_squares_affected + (($supports_array[$block_count] - $mean_supports_affected) ** 2);
							}
							if ($status_array[$block_count] eq "control")
							{
								$sum_squares_normal = $sum_squares_normal + (($supports_array[$block_count] - $mean_supports_normal) ** 2);
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
						if (($variance_affected > 0) || ($variance_normal > 0))
						{
							$t_test = ($mean_supports_affected - $mean_supports_normal) / (sqrt((($variance_affected / $no_of_affecteds) + ($variance_normal / $no_of_normals))));
							$t_test = abs($t_test); # t_statistic can be positve or negative
						}
						
						############################
						# Coefficient of variation #
						############################
						if ($mean_supports_affected > 0 ) {$covar_affected = $sd_affected / $mean_supports_affected;}
						if ($mean_supports_normal > 0) {$covar_normal = $sd_normal / $mean_supports_normal;}

						$covar_all_samples = ($covar_affected + $covar_normal ) / 2;
						
						
						##########################
						# Supporting reads ratio #
						##########################
						if ($mean_supports_normal > 0) {$supports_ratio = $mean_supports_affected / $mean_supports_normal;}

						
						################################
						# Format to two decimal places #
						################################ 
						$covar_affected = sprintf('%.2f', $covar_affected); 
						$covar_normal = sprintf('%.2f', $covar_normal); 
						$supports_ratio = sprintf('%.2f', $supports_ratio );
						$t_test = sprintf('%.3f', $t_test );
						$covar_all_samples = sprintf('%.2f', $covar_all_samples );
						$mean_supports_affected = sprintf('%.2f', $mean_supports_affected);
						$mean_supports_normal = sprintf('%.2f', $mean_supports_normal);
						$sd_affected = sprintf('%.2f', $sd_affected);
						$sd_normal = sprintf('%.2f', $sd_normal);
				
				
						####################################################
						# Only write to output file if it is within region #
						####################################################
						if (($sv_start >= $region_start ) && ( $sv_start <= $region_end))
						{
							print OUT_ALL "\t$t_test\n";
							print OUT_CONCAT "\t$t_test\n";
						}
					
					} # if if ($proportion_supporting > $proportion_supporting_threshold)
					
				} # if data line
				if ($line_count == $ref_line)
				{
					$ref_seq= $single_line;
				}
				
				#######################################
				# Now you are reading the reads lines #
				#######################################
				if ($line_count >= $first_read_line)
				{
					@item=split(/\s+/,$single_line);
					$array_size = scalar @item;
					

					$read_1=$item[1];
					$read_2=$item[2];
					$start_position = $item[4];
					$sample_name = $item[6];
					
					$read_line_count = $line_count - $first_read_line + 1;
					
					#print "Read line $read_line_count (array size: $array_size)\n";
					#print "read_1: $read_1\n";
					#print "read_2: $read_2\n";
					#print "start_position: $start_position\n";
					#print "Sample: $sample_name\n\n";
					
					
					#print "\t$read_1\t$read_2\t$start_position\t$sample_name\n";
				}
				
				
			} # $line_count loop
			
			close OUT_ALL;
					
		} # if file exists
	
		
	} # pass loop
	
	close OUT_CONCAT;
}

print "\n\nOutput file:  \t$pindel_concat_out\n\n";

exit;

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
		
		print "Grouping file $list_count: $deletion_simplified_out into windows of $window_size...\n";
		
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
			
			#print "\nCHECKING WINDOW $window_start - $window_end\n\n";
			
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
			# Now loop through deletions file to find any positions that are in this window  #
			# (There may be some windows for which there are no positions)                   #
			##################################################################################
			$something_in_window="false";
			
			for ($line_count = 1;$line_count<$no_of_lines_deletions_out;$line_count++)
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
				if ($position >= $window_end)
				{
					#print "END_LOOP.  $position > $window_end\n\n";
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
		print "\tReading file $list_count: $deletion_window_out\n";
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
	$output_assoc_file="$run_title"."_supports_w"."$window_size"."_assoc.out";
	
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
				$affected_total_supports = 0;
				$normal_total_supports = 0;
				#$affected_total_supports_normalised = 0;
				#$normal_total_supports_normalised = 0;
				
				$deletion_string = "";
				
				##################################################
				# Put together string of deletion to show user #
				##################################################
				for ($block_count=1;$block_count <=$no_of_blocks;$block_count++)
				{
					$supports_array[$block_count] = $item[$block_count + 1];
					$supports_array[$block_count] = sprintf('%.0f', $supports_array[$block_count]); 

					$deletion_string = "$deletion_string"."\t"."$supports_array[$block_count]";
					
					if ($status_array[$block_count] eq "affected")
					{
						$affected_total_supports = $affected_total_supports + $supports_array[$block_count];
						#$affected_total_supports_normalised = $affected_total_supports_normalised + ($supports_array[$block_count] * $weighting_array[$block_count]);
					}
					
					if ($status_array[$block_count] eq "control")
					{
						$normal_total_supports = $normal_total_supports + $supports_array[$block_count];
						#$normal_total_supports_normalised = $normal_total_supports_normalised + ($supports_array[$block_count] * $weighting_array[$block_count]);
					}
				} # $block_count loop
				
				############################################################
				# Mean deletion counts on each line of deletion_all.out #
				############################################################
				if ($normalise_supports eq "false")
				{
					$mean_supports_affected = $affected_total_supports / $no_of_affecteds;
					$mean_supports_normal = $normal_total_supports / $no_of_normals;
					$mean_supports_both = ($affected_total_supports + $normal_total_supports) / ($no_of_affecteds + $no_of_normals);
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
				$supports_ratio = 0;

				for ($block_count=1;$block_count <=$no_of_blocks;$block_count++)
				{
					if ($status_array[$block_count] eq "affected")
					{
						$sum_squares_affected = $sum_squares_affected + (($supports_array[$block_count] - $mean_supports_affected) ** 2);
					}
					if ($status_array[$block_count] eq "control")
					{
						$sum_squares_normal = $sum_squares_normal + (($supports_array[$block_count] - $mean_supports_normal) ** 2);
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
					$t_test = ($mean_supports_affected - $mean_supports_normal) / (sqrt((($variance_affected / $no_of_affecteds) + ($variance_normal / $no_of_normals))));
					$t_test = abs($t_test); # t_statistic can be positve or negative
				}
				
				############################
				# Coefficient of variation #
				############################
				if ($mean_supports_affected > 0 ) {$covar_affected = $sd_affected / $mean_supports_affected;}
				if ($mean_supports_normal > 0) {$covar_normal = $sd_normal / $mean_supports_normal;}

				$covar_all_samples = ($covar_affected + $covar_normal ) / 2;
				
				
				####################
				# deletion ratio #
				####################
				if ($mean_supports_normal > 0) {$supports_ratio = $mean_supports_affected / $mean_supports_normal;}


				################################
				# Format to two decimal places #
				################################ 
				$covar_affected = sprintf('%.2f', $covar_affected); 
				$covar_normal = sprintf('%.2f', $covar_normal); 
				$supports_ratio = sprintf('%.2f', $supports_ratio );
				$t_test = sprintf('%.2f', $t_test );
				$covar_all_samples = sprintf('%.2f', $covar_all_samples );
				$mean_supports_affected = sprintf('%.2f', $mean_supports_affected);
				$mean_supports_normal = sprintf('%.2f', $mean_supports_normal);
			
			
				######################################
				# Write to output ASSOC file         #
				# (if greater than T-test threshold) #
				######################################
				if ($t_test >= $t_test_threshold)
				{
					print OUT_ASSOC "$chromosome\t$position";
					if ($include_individual_columns eq "true"){print OUT_ASSOC "$deletion_string";}
					print OUT_ASSOC "\t$mean_supports_affected\t$mean_supports_normal";
					print OUT_ASSOC "\t$covar_affected\t$covar_normal";
					print OUT_ASSOC "\t$supports_ratio\t$t_test";
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



