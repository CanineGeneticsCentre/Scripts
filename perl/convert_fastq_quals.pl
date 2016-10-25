#!/usr/bin/perl -w
use strict;
use Term::ANSIColor;

############################################################
# Converts fastq files from old Illumina quality scores    #
# to Sanger quality scores                                 #
############################################################
my $version					= "5";

#################################################
# Sanger scores subtract 33 from the ASCII code #
# Old Illumina subtract 64 from the ASCII code  #
# so the conversion_amount is 31                #
#################################################
my $conversion_amount		= 31; # Subtract this from Illumina ASCII scores to get lower Sanger scores
my $no_of_levels			= 40; # levels in histogram of scores.
my $some_less_than_64 			= "false";

#File names
my $list_file				= "";
my $fastq_file         		= "";
my $fastq_file_new			= "";

#String variables
my $single_line				= "";
my $answer					= "";
my $id_line					= "";
my $seq_line				= "";
my $plus_line				= "";
my $qual_line				= "";
my $qual_line_new			= "";
my $char					= "";
my $char_new				= "";
my $read_file_method		= "";
my $reads					= "";
my $reads2					= "";
my $command					= "";
my $second_column_found		= "";
my $third_column_found		= "";
my $new_sample_name			= "";
my $prefix					= "";
my $testing_mode			= "no";
my $conversion_mode			= "";
my $process_file 			= ""; # yes or no
my $looks_like_illumina 	= "";
my $looks_like_sanger 		= "";


#Number variables
my $no_paired_files			= 0; # This is 2 if multiple paired files are used, and 1 otherwsie
my $line_count				= 0;
my $list_count				= 0;
my $loop_count				= 0;
my $pair_count				= 0;
my $no_of_files				= 0;
my $array_count				= 0;
my $mod						= 0;
my $level					= 0;
my $block_count				= 0;
my $pos						= 0;
my $sanger_more_than_40_count	= 0;
my $sanger_less_than_64_count	= 0;
my $sanger_41_count			= 0;
my $level_size				= 0;
my $threshold				= 0;
my $qual_score				= 0;
my $ascii_score				= 0;
my $ascii_score_new			= 0;
my $old_score				= 0;
my $new_score				= 0;
my $max_qual_score			= 0;
my $min_qual_score			= 1000;
my $max_ascii_score			= 0;
my $min_ascii_score			= 1000;
my $array_size				= 0;
my $max_ascii_score_array	= 0;
my $min_new_score			= 1000;
my $max_new_score			= 0;
my $min_old_score			= 1000;
my $max_old_score			= 0;


my @ascii_score_array		= ();
my @new_score_array			= ();
my @old_score_array			= ();

my @item					= ();
#my @fastq_file_array		= ();
my @fastq_file_1_array		= ();
my @fastq_file_2_array		= ();
my @fastq_file_1_new_array	= ();
my @fastq_file_2_new_array	= ();
my @sample_name_array		= ();

print "\n\n";
print color 'bold magenta';

print "################################\n";
print color 'bold white';
print "      convert_fastq_quals       \n";
print color 'bold magenta';
print "################################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This PERL script converts FASTQ file quality scores from Old to New\n\n";
print "    Some old Illumina FASTQ files have a different way of coding the quality scores.\n";
print "    The GATK programs will not work if some files have the Old scores.\n\n";
print "    (All recent data should be in the New system).\n\n";
print color 'reset';


&print_message("Which direction do you want to convert the quality scores?","input");

print "   <1>\tOld --> New\n";
print "   <2>\tNew --> Old\n\n";
print "   ";
$answer = <STDIN>;
chomp $answer;
if (substr($answer,0,1) eq "1" ){$conversion_mode = "old_to_new"}
if (substr($answer,0,1) eq "2" ){$conversion_mode = "new_to_old"}


&print_message("How do you want to read the input files?","input");

print "   <1>\tMULTIPLE FILES (using a file of file names)\n";
print "   <2>\tSINGLE FILE\n\n";
print "   ";
$answer = <STDIN>;
chomp $answer;
if ($answer eq ""){$answer = "1"} # default

$answer = substr($answer,0,1);

if (substr($answer,0,1) eq "1" ){$read_file_method = "multiple"}
if (substr($answer,0,1) eq "2" ){$read_file_method = "single"}



###################################
# Input file names                #
###################################

if ($read_file_method eq "single")
{
	until (-e $fastq_file)
	{
		print "\nInput FASTQ file:               \t";
		
		$fastq_file = <STDIN>;
		chomp $fastq_file;
		
		if ($fastq_file eq "ls"){print "\n";system ("ls *.fastq")}
		
		if ($fastq_file ne "ls")
		{
			if (! -e $fastq_file){print "\n\n>>>>>>>>  File $fastq_file not found.  Try again.  <<<<<<<<\n\n";}
		}
	}
	

		

	##############################################################
	# Assign to array 1 of file names but just use first element #
	# (array 2 is not used because no_of_pairs = 1)              #
	##############################################################
	$fastq_file_1_array[1]=$fastq_file;
	$no_of_files=1;
	$no_paired_files=1;
	
	print   "Output file (new FASTQ file):   \t";
	$fastq_file_new = <STDIN>;
	chomp $fastq_file_new;
	if (index($fastq_file_new,".fastq") == -1 ){$fastq_file_new = $fastq_file_new.".fastq"}
	
	$fastq_file_1_new_array[1]=$fastq_file_new;

} # End of read = single


if ($read_file_method eq "multiple")
{
	until (-e "$list_file")
	{
		print "\nPlease input the name of your file with a list of file names of the FASTQ files:    (type 'ls' to get a list of .txt files)  ";
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
			print "##########################\n";
			print "#   FILE FORMAT ERROR    #\n";
			print "##########################\n\n";
			print "This file only appears to have a single column\n\n";
			print "For paired-end (PE) files you need two columns separated by a tab.\n\n";
			close LIST;
			exit;
			$no_paired_files=1;
		}
		
		if ($array_size == 2)
		{
			$reads = $item[0];
			$reads2 = $item[1];
			
			#####################################################################
			# Add file type suffixes .fastq if user hasn't added it in the file #
			#####################################################################
			if (index($reads,".fastq") == -1 ){$reads = $reads.".fastq"}
			if (index($reads2,".fastq") == -1 ){$reads2 = $reads2.".fastq"}
			
			$fastq_file_1_array[$list_count]=$reads;
			$fastq_file_2_array[$list_count]=$reads2;
			
			#Output file names
			$prefix = &get_prefix ($reads);
			$fastq_file_1_new_array[$list_count] = "$prefix"."_newquals.fastq";
			
			$prefix = &get_prefix ($reads2);
			$fastq_file_2_new_array[$list_count] = "$prefix"."_newquals.fastq";
			
			$list_count = $list_count + 1;
			
			$no_paired_files=2;
			
		}
		
		# If there is a third column to use as the sample name
		if ($array_size == 3)
		{
			$reads = $item[0];
			$reads2 = $item[1];
			$new_sample_name = $item[2];
			
			#####################################################################
			# Add file type suffixes .fastq if user hasn't added it in the file #
			#####################################################################
			if (index($reads,".fastq") == -1 ){$reads = $reads.".fastq"}
			if (index($reads2,".fastq") == -1 ){$reads2 = $reads2.".fastq"}
			
			$fastq_file_1_array[$list_count] = $reads;
			$fastq_file_2_array[$list_count] = $reads2;
			$sample_name_array[$list_count] = $new_sample_name;
			
			#Output file names
			$fastq_file_1_new_array[$list_count] = "$new_sample_name"."_newquals_1.fastq";
			$fastq_file_2_new_array[$list_count] = "$new_sample_name"."_newquals_2.fastq";
			
			$third_column_found = "true";
			
			$list_count = $list_count + 1;
			
			$no_paired_files=2;
		}
			
	} # while

	close LIST;

	$no_of_files=$list_count - 1;
	
	
	###################
	# List file names #
	###################
	print "\nThere are $no_of_files pairs of 'reads' files in this file of file names.\n\n";
	
	if ($third_column_found eq "true")
	{
		print "There is also a third column in the input file which will be used for rthe output file name.\n\n";
	}
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		if ($third_column_found eq "true")
		{
			print "Pair $list_count	\t$fastq_file_1_array[$list_count]  \t$fastq_file_2_array[$list_count]\t\t$sample_name_array[$list_count]\n";
			print "   Output file 1: $fastq_file_1_new_array[$list_count]\n";
			print "   Output file 2: $fastq_file_2_new_array[$list_count]\n\n";
		}
		if ($third_column_found eq "false")
		{
			print "Pair $list_count	\t$fastq_file_1_array[$list_count]  \t$fastq_file_2_array[$list_count]\n";
		}
	}

} # End of if ($read_file_method eq "multiple")


###############################
# Tell user what is happening #
###############################
if ($conversion_mode eq "old_to_new"){&print_message("Converting quality scores from OLD to NEW...","message");}
if ($conversion_mode eq "new_to_old"){&print_message("Converting quality scores from NEW to OLD...","message");}


if ($testing_mode eq "yes") {$answer = <STDIN>}

for ($loop_count=1;$loop_count <=$no_of_files;$loop_count++)
{

	for ($pair_count = 1;$pair_count <=$no_paired_files;$pair_count++)
	{
		$process_file = "no";
		
		
		if ($read_file_method eq "single")
		{
			$fastq_file = $fastq_file_1_array[$loop_count];
			$fastq_file_new = $fastq_file_1_new_array[$loop_count];
			
			if ($fastq_file ne ""){$process_file = "yes";}
		}
		
		if ($read_file_method eq "multiple")
		{
		
			if ($pair_count == 1) 
			{
				$fastq_file = $fastq_file_1_array[$loop_count];
				$fastq_file_new = $fastq_file_1_new_array[$loop_count];
			
				print "\n\n";
				print "Loop $loop_count  File: $pair_count\n\n";
				print "FASTQ file:     $fastq_file\n";
				print "New FASTQ file: $fastq_file_new\n\n";
				
				if ($fastq_file ne ""){$process_file = "yes";}
			}
		
			if ($pair_count == 2) 
			{
				$fastq_file = $fastq_file_2_array[$loop_count];
				$fastq_file_new = $fastq_file_2_new_array[$loop_count];
				
				print "\n\n";
				print "Loop $loop_count  File: $pair_count\n\n";
				print "FASTQ file:     $fastq_file\n";
				print "New FASTQ file: $fastq_file_new\n\n";
				
				if ($fastq_file ne ""){$process_file = "yes";}
			}
		
		} # read file method is multiple
		
		if ($process_file eq "yes")
		{
			if ($testing_mode eq "yes") {$answer = <STDIN>}
			
			if ($testing_mode eq "no")
			{
			
				open (IN, "$fastq_file") || die "Cannot open $fastq_file";
				open (OUT, ">$fastq_file_new") || die "Cannot create $fastq_file_new";

				######################
				# Set up score array #
				######################
				for ($array_count = 0;$array_count <=126;$array_count++)
				{
					$ascii_score_array[$array_count] = 0;
					$old_score_array[$array_count] = 0;
					$new_score_array[$array_count] = 0;
				}


				############################################
				# Read through FASTQ file a line at a time #
				############################################
				while ($single_line = <IN>)
				{
					$line_count = $line_count + 1;
					$mod = $line_count % 4;
					
					chomp $single_line;
					
					if ($mod == 1)
					{
						$id_line = $single_line;
						print OUT "$id_line\n";
					}
					if ($mod == 2)
					{
						$seq_line = $single_line;
						print OUT "$seq_line\n";
					}
					if ($mod == 3)
					{
						$plus_line = $single_line;
						print OUT "$plus_line\n";
					}

						
					###############################################################
					# Last line of block of four lines is the quality scores line #
					###############################################################
					if ($mod == 0)
					{
						$qual_line = $single_line;
						
						$block_count = $block_count	 + 1;
						
						if (($block_count % 10000) == 0) 
						{
							if ($read_file_method eq "multiple"){print "File: $loop_count\tPair: $pair_count\t$fastq_file\t$block_count\t$conversion_mode\n"}
							if ($read_file_method eq "single"){print "File: $fastq_file\t$block_count\t$conversion_mode\n"}
						}

						$qual_line_new = "";
						
						
						##########################################
						# Run along the string of quality scores #
						##########################################
						for ($pos = 0; $pos < length $qual_line; $pos++)
						{
							$char = substr ($qual_line, $pos, 1);
							
							if (length ($char) > 0) 
							{
								#####################################################
								#  Work out what old or new quality scores would be #
								#####################################################
								$ascii_score = ord($char);
								
								$old_score = ord($char)- 64; # Before Illumina 1.8
								$new_score = ord($char)- 33; # Sanger


								###############################################
								# Any ASCII less than 64 (means it is Sanger) #
								###############################################
								if ($ascii_score < 64){$some_less_than_64 = "true";}


								##################################################
								# Record New scores over 40 (the GATK threshold) #
								##################################################
								if ($new_score > 40)
								{
									$sanger_more_than_40_count = $sanger_more_than_40_count + 1;
								}
								if ($new_score == 41)
								{
									$sanger_41_count = $sanger_41_count + 1;
								}

								########################################
								# Store scores in arrays for Histogram #
								########################################
								if ($ascii_score >= 0 ) {$ascii_score_array[$ascii_score] = $ascii_score_array[$ascii_score] + 1}
								if ($new_score >= 0 ) {$new_score_array[$new_score] = $new_score_array[$new_score] + 1}
								if ($old_score >= 0 ) {$old_score_array[$old_score] = $old_score_array[$old_score] + 1}


								################################################
								#  Create new ASCII value for new quality line #
								################################################

								if ($conversion_mode eq "old_to_new")
								{
									$ascii_score_new = $ascii_score - $conversion_amount;
								}
								if ($conversion_mode eq "new_to_old")
								{
									$ascii_score_new = $ascii_score + $conversion_amount;
								}
								
								$char_new = chr($ascii_score_new);
								

								######################################
								# Continue building new quality line #
								######################################
								$qual_line_new = $qual_line_new.$char_new;
								

							} # if char <> ""
							

							# Get max and min quals #
							if ($ascii_score > $max_ascii_score)
							{
								$max_ascii_score = $ascii_score;
								
								$max_old_score = $max_ascii_score- 64;
								$max_new_score = $max_ascii_score- 33;
								
							}
							if ($ascii_score < $min_ascii_score)
							{
								$min_ascii_score = $ascii_score;
								
								$min_old_score = $min_ascii_score- 64;
								$min_new_score = $min_ascii_score- 33;
							}
								
							
						} # End of for pos loop
						
						print OUT "$qual_line_new\n";
						
					} # if $mod == 0 (ie a quality scores line)
					
				} # while $single_line <IN>

				#####################################
				# List distribution of ASCII scores #
				#####################################

				print "\n\n\n";
				print "Distribution of ASCII scores.\n\n";
				print "=============================\n\n";
			
				print "ASCII\tNumber\n";
				
				$max_ascii_score_array = 0;
				
				for ($array_count = 30; $array_count <=126; $array_count++)
				{
					if ($ascii_score_array[$array_count] > 0)
					{
						print "$array_count\t$ascii_score_array[$array_count]\n";
						if ($ascii_score_array[$array_count] > $max_ascii_score_array)
						{
							$max_ascii_score_array = $ascii_score_array[$array_count];
						}
					}
					
				}
					


				#############################
				# Histogram of ASCII values #
				#############################

				$level_size = $max_ascii_score_array / $no_of_levels;
				print "\n\n\n\n\n\n\n";

				for ($level=0;$level<$no_of_levels;$level++)
				{
					$threshold= ($no_of_levels - $level) * $level_size;

					for ($array_count = 30; $array_count <=126;$array_count++)
					{		
						if ($ascii_score_array[$array_count] >= $threshold){print "|";}else{print " ";}	
					}
					print "\n";
				}



				print "\n30........40........50........60........70........80........90........100.......110.......120.......130\n";

				print "\n\n";
				
				$min_qual_score = $min_ascii_score - 33;
				$max_qual_score = $max_ascii_score - 33;
				
		
				print "Histogram of PRE-CONVERSION ASCII scores. \tASCII scores $min_ascii_score-$max_ascii_score\n\n";

				print " >> For OLD quality scores subtract 64, for NEW quality scores subtract 33 (from the ASCII value)\n\n";


				print "    If quality scores in this file are NEW they would go from:       \t$min_new_score-$max_new_score\n";
				print "    If quality scores in this file are OLD they would go from:       \t$min_old_score-$max_old_score\n\n";
				
				#Decide which system
				print "Checks on pre-conversion file:\n\n";

				print "\nFactors making pre-conversion file it look like New scores\n\n";
				
				if ($max_new_score <= 40)
				{
					print "\tThere were no New scores over 40\n";
					$looks_like_sanger = "true";
				}
				
				if ($some_less_than_64 eq "true")
				{
					print "\tThere were some ASCII scores less than 64.\n";
					$looks_like_sanger = "true";
				}
				
				if ($looks_like_sanger eq "")
				{
					print "\tNONE\n";
				}
				
				print "\n\nFactors making pre-conversion file look like Old scores\n\n";
				
				
				if (($max_new_score > 40) && ($max_old_score <=40))
				{
					print "\tThe maximum New score was over 40 and the maximum Old score was not over 40\n";
					$looks_like_illumina = "true";
				}
				
				
				if ($sanger_more_than_40_count > 0)
				{
					print "\tThere are $sanger_more_than_40_count scores more than 40. ($sanger_41_count scores of 41).  This will cause a problem in GATK\n";
					$looks_like_illumina = "true";
				}
				
				if ($some_less_than_64 eq "false")
				{
					print "\tThere are no ASCII scores less than 64.\n";
					$looks_like_illumina = "true";
				}

				if ($looks_like_illumina eq "")
				{
					print "\tNONE\n";
				}
				
				&print_message("CONVERSION of $fastq_file FINISHED","message");
				print "\n";

				if ($conversion_mode eq "old_to_new")
				{
					print "Original FASTQ file (Old quality scores):      \t$fastq_file\n";
					print "Converted FASTQ file (New quality scores):     \t$fastq_file_new\n\n";

					print "(QUAL scores have been reduced by $conversion_amount in the output file.)\n\n";
				}
				if ($conversion_mode eq "new_to_old")
				{
					print "Original FASTQ file (New quality scores):      \t$fastq_file\n";
					print "Converted FASTQ file (Old quality scores):     \t$fastq_file_new\n\n";

					print "(QUAL scores have been increased by $conversion_amount in the output file.)\n\n";
				}


				
			
			} # end of 'if testing_mode eq "yes" 
		
		
		} # if $process_file eq "yes"
	
	} # end of pair loop
	
	
	if ($read_file_method eq "multiple")
	{
		print_message("End of loop $loop_count","message");
	}
	
} # end of main loop

exit;

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
