#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	fastq_merge  	                                                    #
#   merges fastq files with unix cat command                            #
#									                                    #
#   Details:                                                            #
#									                                    #
#   When we get FASTQ files from Whole Genome Sequencing they often     #
#   come as two separate files.  These have to be merged, and we use    #
#   unix 'cat' command.                                                 #
#									                                    #
#   If this was done manually, and if the wrong files were used,        #
#   there would be no way of tracing the error.  So it is important     #
#   to use this script, which flags up any erros, and stores a detailed #
#   log file, which records all the file sizes as well.                 #
#########################################################################

###############################
# Mike Boursnell Nov 2014     #
# Animal Health Trust         #
# Newmarket                   #
# UK                          #
# mike.boursnelln@aht.org.uk  #
###############################

use strict;
use File::Basename ;
use Term::ANSIColor;

# Version
my $version						= "13";

# Constants
my $from 						= 'fastq_merge@gen-x1404-ws01.aht.org.uk'; # Who e-mails come from

#FILE NAMES
my $list_file					= ""; # file with list of BAM file names
my $fastq_file					= "";
my $fastq_file_1_gzipped		= "";
my $fastq_file_2_gzipped		= "";
my $fastq_files_to_merge		= ""; #string of FASTQ files to merge
my $log_file					= "";

#FILE NAMES 2-dimensional method
my $fastq_file_merged			= "";
my $fastq_file_merged_gzipped	= "";

#STRINGS
my $command						= ""; # unix command
my $merge_string		    	= ""; # string of input files for picard
my $answer						= "";
my $tempdir						= ""; # java temporary directory
my $temp_dir_string				= "";
my $prefix						= ""; # Part of file name before the dot.
my $single_line					= "";
my $start_time					= "";
my $end_time					= "";
my $run_time					= "";
my $email_address				= "";

#NUMBERS
my $row_count					= 0; # rows in input file
my $no_of_files					= 0;
my $no_of_fastq_files_to_merge	= 0;
my $array_size					= 0;
my $filesize_1					= 0;
my $filesize_2					= 0;
my $filesize_merged				= 0;
my $filesize_merged_estimate	= 0;
my $column_count				= 0;
my $file_count					= 0; # column_count + 1
my $line_count					= 0; # There should only be 2 lines in the file
my $merged_name_column			= 0; # Number of last column which should contain the name for the merged file

#BOOLEAN
my $some_index_files_missing	="false";
my $make_index_files			= "no";
my $merge_sequence_dictionaries = "no";
my $not_all_files_found			= ""; # true or false
my $overwrite_previous_merged_file = ""; # true or false
my $overwrite_previous_gzipped_file = ""; # true or false

#ARRAYS
# 2 dimensional arrays for storing file names and sizes
my @file_number 				= ();
my @row_number 					= ();
my @fastq_file_array_2d 		= (\@file_number, \@row_number);
my @file_number_2				= ();
my @row_number_2 					= ();
my @fastq_filesize_array_2d		= (\@file_number_2, \@row_number_2); # is this OK?

#ARRAYS
my @item						= ();
my @fastq_files_string_array	= ();
my @merged_file_array			= ();
my @fastq_files_size_array		= ();
my @merged_filesize_array		= ();
my @merged_filesize_array_calculated		= ();

my @fastq_file_column_array			= (); # If there are more then 2 files to merge
my @fastq_file_column_array_gzipped	= (); # If there are more then 2 files to merge

#File sizes
my @fastq_filesize_1_array		= ();
my @fastq_filesize_2_array		= ();




print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "              fastq_merge   \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "    - This program runs the unix 'cat' command to merge several FASTQ files into one\n\n";
print "    - All the details are recorded in a log file, which should be kept.\n\n";

print "The input file contains a list of FASTQ files to be merged\n\n";

print "The format is tab or space-separated columns with the input FASTQ files and the name of the output merged FASTQ file on each of two lines\n\n";

print "For example:  ABC_01_R1.fastq   ABC_02_R1.fastq   ABC_03_R1.fastq   ABC_1.fastq \n";
print "              ABC_01_R2.fastq   ABC_02_R2.fastq   ABC_03_R2.fastq   ABC_2.fastq \n\n";

print color 'reset';

until (-e "$list_file")
{
	&print_message("Please input the name of your file with a list of file names of the FASTQ files to merge","input");

	print ">  ";

	$list_file = <STDIN>;
	chomp $list_file;

	if (($list_file eq "ls") || ($list_file eq ""))
	{
		print "\nList of .txt files...\n\n";
		system ("ls *.txt");
		print "\n";
	}
	
	if (($list_file ne "ls") && ($list_file ne ""))
	{
		if (! -e "$list_file"){print "\nFile doesn't exist. Try again...    \n";}
	}
}

$prefix = &get_prefix($list_file);

$log_file = $prefix."_fastq_merge_log.out";
open (LOG, ">$log_file")|| die "Cannot create output file: $log_file";

print LOG "Log file of fastq_merge version $version\n\n";

print LOG "File containing list of FASTQ files to merge: $list_file\n\n";

print LOG "Contents of this input file:\n\n";


#############################################
# Make sure the list file is in Unix format #
#############################################
print "\n";
$command = "dos2unix $list_file";
system("$command");



####################################################
# Open the list file to get the list of file names #
####################################################
open (LIST, "$list_file") || die "Cannot open $list_file";
$row_count = 0;
$line_count = 0;

while ($single_line = <LIST> ) 
{
	$line_count = $line_count + 1;

	chomp $single_line;

	print LOG "\t$single_line\n";

	@item=split(/\s+/,$single_line);
				
	$array_size = scalar @item;
	
	if ($array_size == 1)
	{
		print "\n\n";
		print "##########################\n";
		print "#   FILE FORMAT ERROR    #\n";
		print "##########################\n\n";
		print "This file only appears to have a single column\n\n";
		print "For paired-end (PE) files you need at least two columns separated by a tab, plus a third for the merged name\n\n";
		close LIST;
		exit;
	}# array_size = 1
	
	if ($array_size == 2)
	{
		print "\n\n";
		print "##########################\n";
		print "#   FILE FORMAT ERROR    #\n";
		print "##########################\n\n";
		print "This file only appears to have a two columns\n\n";
		print "For paired-end (PE) files you need at least two columns separated by a tab, plus a third for the merged name\n\n";
		close LIST;
		exit;
	}# array_size = 2


	#############################################################
	# This now deals with any number of FASTQ files to merge    #
	# There are two rows in the files and a number of columns   #
	# The number of columns is the number of files to be merged #
	# plus the name of the merged file.                         #
	#############################################################
	if ($array_size > 2)
	{
		$row_count = $row_count + 1; # The files has two rows, 1 each of the paired end file list

		$no_of_fastq_files_to_merge = $array_size - 1;
		$merged_name_column = $array_size;



		###################################
		# Loop for FASTQ files to merge...#
		###################################
		for ($column_count = 0; $column_count <  $no_of_fastq_files_to_merge; $column_count++)
		{
			#############################################################
			# file_count is just the columns ( In unix the first column #
			# is zero, so we make it 1 which seems more intuitive)      #
			#############################################################
			$file_count = $column_count + 1;

			$fastq_file_array_2d[$file_count][$row_count] = $item[$column_count];

			#####################################################
			# If the input filenames refer to the gzipped files #
			#####################################################
			if (index($item[$column_count],".gz") > -1)
			{
				print "Column $file_count: $fastq_file_column_array[$file_count]\tPLEASE GUNZIP THIS FILE\n";

				&print_message("All files have to be gunzipped before merging","warning");

				exit;
			}


			#####################################################
			# Make string for all FASTQ files                   #
			#####################################################
			if ($file_count == 1)
			{
				$fastq_files_string_array[$row_count] = "$fastq_file_array_2d[$file_count][$row_count]";
			}
			if ($file_count > 1)
			{
				$fastq_files_string_array[$row_count] = "$fastq_files_string_array[$row_count]\t$fastq_file_array_2d[$file_count][$row_count]";
			}

		} # column_count loop 


		####################################
		# Merged file is last in the line  #
		####################################
		$merged_file_array[$row_count] = $item[$array_size - 1];

		if (index($merged_file_array[$row_count],".fastq") == -1)
		{
			$merged_file_array[$row_count] = $merged_file_array[$row_count].".fastq";
		}
		
	}# array_size > 2

} # While $single_line (going down the lines in input file - probably 2 lines)

close LIST;

########################################################
# This is the number of columns (not counting the last #
# column which contains the file to be merged into     #
########################################################
$no_of_files=$file_count;



######################################
# E-mail address for notifications   #
######################################

&print_message("Please enter your email address","input");
$email_address = <STDIN>;
chomp $email_address;

# Some short cuts
if ($email_address eq "m"){$email_address = 'mike.boursnell@aht.org.uk';}
if ($email_address eq "o"){$email_address = 'oliver.forman@aht.org.uk';}
if ($email_address eq "s"){$email_address = 'sally.ricketts@aht.org.uk';}
if ($email_address eq "k"){$email_address = 'karen.steward@aht.org.uk';}
if ($email_address eq "a"){$email_address = 'amy.charbonneau@aht.org.uk';}
if ($email_address eq "b"){$email_address = 'rebekkah.hitti@aht.org.uk';}
if ($email_address eq "g"){$email_address = 'graham.newland@aht.org.uk';}
if ($email_address eq "l"){$email_address = 'louise.pettitt@aht.org.uk';}

print LOG "\n\nE-mail address: \t$email_address\n\n";


########################################################
# Now check everything and make string for running cat #
########################################################

for ($row_count=1;$row_count<=2;$row_count++)
{
	for ($file_count=1;$file_count<=$no_of_files;$file_count++)
	{

		# Check file exists and record size
		if (-e $fastq_file_array_2d[$file_count][$row_count])
		{
			#Record size
			$fastq_filesize_array_2d[$file_count][$row_count] = -s "$fastq_file_array_2d[$file_count][$row_count]";

			#Make string
			if ($file_count == 1)
			{
				$fastq_files_string_array[$row_count] = "$fastq_file_array_2d[$file_count][$row_count]";
			}
			else
			{
				$fastq_files_string_array[$row_count] = "$fastq_files_string_array[$row_count]\t$fastq_file_array_2d[$file_count][$row_count]";
			}

		}
		else
		{
			&print_message("FASTQ file $fastq_file_array_2d[$file_count][$row_count] cannot be found!","warning");
			print LOG "FASTQ file $fastq_file_array_2d[$file_count][$row_count] cannot be found!\n\n";
			$not_all_files_found = "true";
		}
	}

	#
	# Checked if merged file already exists #

}

####################################
# Warn if not all files were found #
####################################
if ($not_all_files_found eq "true")
{
	&print_message("Not all the FASTQ files were found!","warning");
	
	print "Please check the input file $list_file carefully\n\n";

	print LOG "Not all the FASTQ files were found!\n";
	print LOG "Script terminated\n";

	exit;
}


######################################
# List file names for user to check  #
######################################
&print_message("PLEASE CHECK THE FOLLOWING *VERY* CAREFULLY","message");

for ($row_count=1;$row_count<=2;$row_count++)
{
	print "Files on row $row_count\n\n";

	for ($file_count=1;$file_count<=$no_of_files;$file_count++)
	{
		print "   File $file_count: $fastq_file_array_2d[$file_count][$row_count]\n";
	}
	print "\n";
	print "      ==> To be merged to: $merged_file_array[$row_count]\n\n";
}

print "\nStrings for unix 'cat' command\n\n";
for ($row_count=1;$row_count<=2;$row_count++)
{
	print "  $fastq_files_string_array[$row_count]";
	print " > $merged_file_array[$row_count]\n\n";
}


print "\n\nAre these the correct files to be merged, and is the output file name correct?\n\n";

print "If all is ";
print color 'bold red';
print "DEFINITELY";
print color 'reset';
print " correct, press 'return' otherwise press 'Q' to quit\n";



$answer=<STDIN>;
chomp $answer;
if (lc $answer eq "q"){exit;}


######################################
# Warn if merged file already exists #
######################################
for ($row_count=1;$row_count<=2;$row_count++)
{
	$fastq_file_merged = $merged_file_array[$row_count];
	$fastq_file_merged_gzipped = $fastq_file_merged.".gz";

	if (-e $fastq_file_merged)
	{
		&print_message("Merged file $fastq_file_merged already exists","warning");

		print "Do you want to overwrite it? (y/n)";

		$answer=<STDIN>;
		chomp $answer;

		if (lc $answer eq "y")
		{
			$overwrite_previous_merged_file = "true";
		}
		else
		{
			$overwrite_previous_merged_file = "false";
			print "\nOK.  Exiting program\n\n";
			exit;
		}
	} # if merged file exists

	if (-e $fastq_file_merged_gzipped)
	{
		&print_message("Merged and gzipped file $fastq_file_merged_gzipped already exists","warning");

		print "Do you want to overwrite it? (y/n)";

		$answer=<STDIN>;
		chomp $answer;

		if (lc $answer eq "y")
		{
			$overwrite_previous_gzipped_file = "true";
		}
		else
		{
			$overwrite_previous_gzipped_file = "false";
			print "\nOK.  Exiting program\n\n";
			exit;
		}
	} # if merged gzipped file exists

} # row_count loop


&print_message("Starting to merge FASTQ files","message");

print "\nCalculating the expected size of the merged file...\n\n";
for ($row_count=1;$row_count<=$no_of_files;$row_count++)
{	
	$merged_filesize_array_calculated[$row_count] = 0;

	for ($file_count=1;$file_count<=$no_of_files;$file_count++)
	{
		$merged_filesize_array_calculated[$row_count] = $merged_filesize_array_calculated[$row_count] + $fastq_filesize_array_2d[$file_count][$row_count];
	}
}


################################################
# Use unix 'cat' command to merge              #
################################################
for ($row_count=1;$row_count<=$no_of_files;$row_count++)
{
	
	&print_message("File $row_count/$no_of_files    Merging...","message");

	######################################
	# Get file names from the arrays     #
	######################################
	$fastq_files_to_merge = $fastq_files_string_array[$row_count];
	$fastq_file_merged = $merged_file_array[$row_count];
	
	$start_time = time();
	print "\n";

	########################################################
	# Run the actual unix 'cat' command to merge the files #
	########################################################
	&run_unix_command("cat $fastq_files_to_merge > $fastq_file_merged");

	
	###################################
	# Get size of merged FASTQ file   #
	###################################
	if (-e $fastq_file_merged)
	{
		$merged_filesize_array[$row_count] = -s "$fastq_file_merged";
	}

	$end_time = time();
	$run_time = $end_time - $start_time;
	&log_time($run_time, "\nTime taken to run cat on files $fastq_files_to_merge");

	&print_message("File $row_count/$no_of_files.  Finished merging files","message");

	print LOG "\n\n###########################################################################n\n";

} # list count loop for merging with cat

print LOG "Merging of FASTQ files with fastq_merge has finished\n";


print "\nThese are the sizes of all the files after merging:\n\n";
for ($row_count=1;$row_count<=$no_of_files;$row_count++)
{	
	for ($file_count=1;$file_count<=$no_of_files;$file_count++)
	{
		print "   $fastq_file_array_2d[$file_count][$row_count]  \tSize: $fastq_filesize_array_2d[$file_count][$row_count]\n";
	}
	print "\n";
	print "  $merged_file_array[$row_count]\tSize: $merged_filesize_array[$row_count]\n\n";
	print "  Expected size for merged file: $merged_filesize_array_calculated[$row_count]\n\n";
}


################################################
# gzip the final FASTQ files                   #
################################################

&print_message("Gzipping final merged FASTQ file","message");

for ($row_count=1;$row_count<=$no_of_files;$row_count++)
{
	$fastq_file_merged = $merged_file_array[$row_count];

	if ($overwrite_previous_gzipped_file eq "true")
	{
		&run_unix_command("gzip -f $fastq_file_merged");
	}
	else
	{
		&run_unix_command("gzip $fastq_file_merged");
	}

} # list count loop for gzipping

print LOG "Gzipping of FASTQ files with fastq_merge has finished\n";

close LOG;


#####################################################
# Copy log to /home/genetics/command_logs directory #
#####################################################
&run_unix_command("cp $log_file /home/genetics/command_logs/$log_file");


&print_message("FINISHED MERGING FASTQ FILES","message");

print "Your merging of FASTQ files with fastq_merge has finished\n\n";
print "All the details are recorded in the log file:  $log_file\n\n";



#####################################################
# Send a mail to let user know it has finished      #
#####################################################
open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $email_address\n";
print MAIL "From: $from\n";
print MAIL "Subject: fastq_merge has finished\n\n";

## Mail Body
print MAIL "Script:     fastq_merge version $version\n\n";

print MAIL "Your merging of FASTQ files with fastq_merge has finished\n\n";

for ($row_count=1;$row_count<=$no_of_files;$row_count++)
{
	$fastq_file_merged = $merged_file_array[$row_count];

	print MAIL "  $fastq_file_merged\n";
} # list count loop for gzipping


print MAIL "\n\nThese are the sizes of all the files after merging:\n\n";
for ($row_count=1;$row_count<=$no_of_files;$row_count++)
{	
	for ($file_count=1;$file_count<=$no_of_files;$file_count++)
	{
		print MAIL "  $fastq_file_array_2d[$file_count][$row_count]  \tSize: $fastq_filesize_array_2d[$file_count][$row_count]\n";
	}
	print MAIL "\n";
	print MAIL "   $merged_file_array[$row_count] \tSize: $merged_filesize_array[$row_count]\n\n";
	print MAIL "   Expected size for merged file: \t$merged_filesize_array_calculated[$row_count]\n\n\n";
}

close(MAIL);

exit;

########################################################################################################
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

}#
 
###########################################################
# Subroutine to get filename before file type suffix      #
# (i.e. prefix) (e.g. to get "test" from "test.fasta")    #
###########################################################

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

##########################################################
# Subroutine to pause until user hits 'return'           #
##########################################################
sub pause
{
	my $_answer = "";
	print "\n Press RETURN to continue\n";
	$_answer=<STDIN>;
}

##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{
	foreach (@_) {s/\n//g}  
	foreach (@_) {s/\r//g}  
}



#############################################
# Subroutine to execute unix command        #
# Argument 1: command to execute            #
#############################################
sub run_unix_command
{
	my $_unix_command = "";
	$_unix_command = $_[0];
		
	print(" >> $_unix_command\n\n");
	system("$_unix_command");
}
	

##############################################
# Subroutine to record size of INPUT file    #
# (without sending an e-mail if it is zero)  #
##############################################
sub record_input_file_size
{
	my $_inputfile = "";
	my $_filesize = "";	

	$_inputfile = $_[0];
	
	if (-e $_inputfile)
	{
		$_filesize = -s "$_inputfile";
		print "Input file: $_inputfile\t\tSize: $_filesize\n";
		print LOG "Input file: $_inputfile\t\tSize: $_filesize\n";
	}
	else
	{
		$_filesize=0;
		print "Input file: $_inputfile\t\tSize: Not found\n";
		print LOG "Input file: $_inputfile\t\tSize: Not found\n";
	}
} # record_input_file_size



##############################################
# Subroutine to record size of output file   #
# (without sending an e-mail if it is zero)  #
##############################################

sub record_output_file_size
{
	my $_outputfile = "";
	my $_filesize = "";	

	$_outputfile = $_[0];
	
	if (-e $_outputfile)
	{
		$_filesize = -s "$_outputfile";
		print "Output file: $_outputfile\t\tSize: $_filesize\n";
		print LOG "Output file: $_outputfile\t\tSize: $_filesize\n";
	}
	else
	{
		$_filesize=0;
		print "Output file: $_outputfile\t\tSize: Not found\n";
		print LOG "Output file: $_outputfile\t\tSize: Not found\n";
	}
} # record_output_file_size



sub log_time
{
	my $time_to_log = "";
	my $description = "";
	$time_to_log = $_[0];
	$description = $_[1];

	printf "$description: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $time_to_log)[7,2,1,0];
}
	

	####################################################
# Subroutine to print to screen and to COMMAND_LOG #
####################################################
sub print_both
{
	my $_message = $_[0];

	print "$_message";
	print LOG "$_message";
}