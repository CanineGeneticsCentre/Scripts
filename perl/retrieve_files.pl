#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	retrieve_files  	                                                    #
#   bRetrieves files from the Workstation working area                    #
#   to the Workstation local_storage                                    #
#									                                    #
#   Details:                                                            #
#									                                    #
#   Because the Workstation has limited disc space, we need to back up  #
#   files regularly.  It is very importnat that this whole process      #
#   is recorded somewhere.  This program does the bcak up, and then     #
#   stores a log file in /home/genetics/comand_logs                     #
#########################################################################

###############################
# Mike Boursnell Nov 2014     #
# Animal Health Trust         #
# Newmarket                   #
# UK                          #
# mike.boursnelln@aht.org.uk  #
###############################

use strict;
use Getopt::Long;
use File::Basename;
use Term::ANSIColor;
use Cwd;

# VERSION OF SOFTWARE #
my $version						= "9";

#Directories
my $local_storage_path			= "/ahtlocal/local_storage/";

#File names
my $list_file					= "";
my $command_log					= "";
my $file_in_local				= "";
my $full_file_path				= "";

#Strings
my $command						= "";
my $file_to_be_retrieved		= "";
my $email_address				= "";
my $answer						= "";
my $from 						= 'retrieve_files@unix.aht.org.uk'; # Who e-mails come from
my $source_directory			= "/ahtlocal/local_storage"; # This is a constant but could be changed
my $destination_directory		= ""; 
my $read_file_method			= ""; # single or multiple
my $start_time					= "";
my $end_time					= "";
my $run_time					= "";

#Numbers
my $no_of_lines					= 0;
my $file_count					= 0;
my $file_exist_count			= 0;
my $total_no_files				= 0;
my $file_size_source 			= 0;
my $file_size_destination 		= 0;
my $dir_count					= 0;
my $total_dir_count				= 0;
my $number_of_dirs_1			= 0;
my $number_of_dirs_2			= 0;
my $number_of_dirs_3			= 0;

#Arrays
my @file_array					= ();
my @dir_array					= ();

#Boolean
my $gunzip_files					= "false"; # true or false

#Date
my $date						= 'date';

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900; 



print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      retrieve_files       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';

if ( index($version,"w") == -1 ){ print "Version $version (for samba64)\n\n"; } else { print "Version $version (for Workstation)\n\n"; }

print "  - This program retrieves files from /ahtlocal/local_storage\n";
print "    A list of files to be retrieved are provided in a text file\n";
print "    You are e-mailed when it's finished\n\n";

print "  - The files are COPIED not MOVED so no files are deleted.\n\n";


$command_log = "retrieve_log"."_"."$hour"."."."$min"."_"."$mday"."_"."$mon"."_"."$year".".out";
$destination_directory=getcwd; # i.e. current directory

print "Command log file: $command_log\n\n";

######################################
# E-mail address for notifications   #
######################################

&print_message("Please enter your email address to be notified when files are retrieved","input");
$email_address = <STDIN>;
chomp $email_address;


# Some short cuts for e-mails
if ($email_address eq "a"){$email_address = 'amy.charbonneau@aht.org.uk';}
if ($email_address eq "b"){$email_address = 'rebekkah.hitti@aht.org.uk';}
if ($email_address eq "c"){$email_address = 'chris.jenkinsd@aht.org.uk';}
if ($email_address eq "g"){$email_address = 'graham.newland@aht.org.uk';}
if ($email_address eq "k"){$email_address = 'karen.steward@aht.org.uk';}
if ($email_address eq "m"){$email_address = 'mike.boursnell@aht.org.uk';}
if ($email_address eq "o"){$email_address = 'oliver.forman@aht.org.uk';}
if ($email_address eq "s"){$email_address = 'sally.ricketts@aht.org.uk';}



#####################################################################
# ASK IF YOU WANT TO READ THE FILENAMES FROM A "FILE OF FILE NAMES" #
#####################################################################
$read_file_method = "single";


######################################
# Get file with list of files        #
######################################
if ($read_file_method eq "multiple")
{
	######################################
	# Get file with list of files        #
	######################################
	until (-e "$list_file")
	{
		&print_message("Please input the name of your file with a list of files to be retrieved to $destination_directory","input");
		
		$list_file = <STDIN>;chomp $list_file;
		
		if ($list_file eq "ls"){print "\n";system ("ls *.txt"); print "\n";}
	}
	#############################################
	# Make sure the list file is in Unix format #
	#############################################

	$command = "dos2unix $list_file";
	system("$command");
	print "\n\n";

	####################################################
	# Open the list file to get the list of files      #
	####################################################
	open (LIST, "$list_file") || die "Cannot open $list_file";
	$no_of_lines=0;
	while ($file_to_be_retrieved = <LIST> ) 
	{
		chomp $file_to_be_retrieved;

		if (-e "$source_directory/$file_to_be_retrieved")
		{
			$file_count = $file_count + 1;
			$file_array[$file_count] = $file_to_be_retrieved;
		}
		if (! -e "$source_directory/$file_to_be_retrieved")
		{
			print "File $file_to_be_retrieved cannot be found\n";
		}
	}
} # if read_file_method eq multiple


######################################
# Get files one by one               #
######################################
if ($read_file_method eq "single")
{
	$file_count = 1;

	####################################################
	# Ask user which directory to use                  #
	####################################################
	&get_list_of_all_directories;

	&print_message("Which directory do you want to retrieve the files from?","input");

	for ($dir_count = 1; $dir_count <= $total_dir_count; $dir_count++)
	{
			print "<$dir_count>\t$dir_array[$dir_count]\n";
	}

	print "\n>  ";
	$answer=<STDIN>;
	chomp $answer;

	$source_directory = $dir_array[$answer];
	if (-e $source_directory)
	{
		print "\nList of files in $source_directory:\n\n";
	}
	else
	{
		print "\n\nERROR: $source_directory can't be found\n\n";
		exit;
	}


	######################
	# Show list of files #
	######################
	system ("ls $source_directory");

	&print_message("Please input the names of your files to be retrieved (q to end list)","input");
	print "Source directory:      \t$source_directory\n";
	print "Destination directory: \t$destination_directory\n\n";

	until (lc $file_to_be_retrieved eq "q")
	{
		
		print "File $file_count:    ";
		$file_to_be_retrieved = <STDIN>;
		chomp $file_to_be_retrieved;


		print "\n";

		if ($file_to_be_retrieved eq "ls"){print "\n";system ("ls $source_directory/"); print "\n";}

		if (($file_to_be_retrieved ne "ls") && (lc $file_to_be_retrieved ne "q"))
		{
			if (! -e "$source_directory/$file_to_be_retrieved")
			{
				print "File $file_to_be_retrieved cannot be found\n\n";
			}
			if (-e "$source_directory/$file_to_be_retrieved")
			{
				$file_array[$file_count] = $file_to_be_retrieved;
				$file_count = $file_count + 1;
				$file_exist_count = $file_exist_count + 1;
			}
		}
	}
} # if read_file_method eq single

$total_no_files = $file_exist_count;

print "Total number of files: $total_no_files\n\n";


####################################################
# Ask user if the files should be gzipped          #
####################################################
&print_message("Do you want the files gun-zipped when they get to the destination folder?","input");

print "   <1> YES - gunzip the files\n";
print "   <2> NO\n\n";
$answer=<STDIN>;
if (substr($answer,0,1) eq "1" ){$gunzip_files = "true"} else {$gunzip_files = "false"}


####################################################
# Open the list file to SHOW the list of files     #
####################################################
&print_message("List of $total_no_files files to be retrieved from $source_directory","message");


for ($file_count = 1; $file_count <=$total_no_files; $file_count++)
{
	print "File $file_count: $file_array[$file_count]\n";
}


print "\nPress 'Return' if these files are OK, or press 'Q' to quit\n\n";
$answer=<STDIN>;
chomp $answer;

if (lc $answer eq "q"){exit;}

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

print COMMAND_LOG "Command log of perl script retrieve_files, version $version\n\n";

print COMMAND_LOG "List of files retrieved\n\n";



################
# START TIMERS #
################
$start_time = time();



####################################################
# Open the list file to get the list of files      #
####################################################
$no_of_lines=0;
for ($file_count = 1; $file_count <=$total_no_files; $file_count++)
{
	$file_to_be_retrieved = $file_array[$file_count];

	if (-e "$source_directory/$file_to_be_retrieved")
	{
		$command = "cp $source_directory/$file_to_be_retrieved $destination_directory/$file_to_be_retrieved";
		print "$command\n";
		
		print COMMAND_LOG "\n=====================================================================\n";
		print COMMAND_LOG "FILE:    $file_to_be_retrieved\n";
		print COMMAND_LOG "COMMAND: $command\n\n";
		system ("$command");
		

		#######################################
		# Write file sizes to the command log #
		#######################################

		$file_size_source = -s "$source_directory/$file_to_be_retrieved";
		$file_size_destination = -s "$destination_directory/$file_to_be_retrieved";

		print COMMAND_LOG "Check on file sizes:\n\n";

		print COMMAND_LOG "File: $file_to_be_retrieved\n";
		print COMMAND_LOG "Source file size: $file_size_source\n";
		print COMMAND_LOG "Source file size: $file_size_destination\n\n";


		if ($gunzip_files eq "true")
		{
			$command = "gunzip $destination_directory/$file_to_be_retrieved";
			print "$command\n";
			print COMMAND_LOG "GUNZIP:    $command\n\n";
			system ("$command");
		}
	}
}


####################################################################
# Calculate the elapsed total run time at this point (end of loop) #
####################################################################
$end_time = time();
$run_time = $end_time - $start_time;


##########################################
# Send an e-mail when finished           #
##########################################
open(MAIL, "|/usr/sbin/sendmail -t");
 
## Mail Header
print MAIL "To: $email_address\n";
print MAIL "From: command_list\n";
print MAIL "Subject: Retrieving of files in $list_file to $destination_directory\n\n";
## Mail Body
print MAIL "Program: retrieve_files\tVersion $version\n\n";

print MAIL "This program running the retrieving of your files from $source_directory has finished\n\n";

if ($read_file_method eq "multiple"){print MAIL "File with list of files to be retrieved:     \t$list_file\n\n";}

print MAIL "Files to be retrieved:\n\n";

for ($file_count = 1; $file_count <=$total_no_files; $file_count++)
{
	print MAIL "\t$file_array[$file_count]\n";
}

print MAIL "\nTime taken for retrieving of files : $run_time seconds\n";
printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];

close MAIL;





######################################################
# Copy command log to command log folder in genetics #
######################################################

&run_unix_command("cp $command_log /home/genetics/command_logs/$command_log","Copy retrieve_files command log to genetics folder");

print COMMAND_LOG "\nTime taken for retrieving of files : $run_time seconds\n";
printf COMMAND_LOG"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];

close COMMAND_LOG;

&print_message("Retrieving program completed","message");

print "All files in $list_file have been copied to $destination_directory\n\n";




exit;

######################################
# Subroutine to print screen message #
######################################
sub print_message
{
	my $_message_length 	= "";
	my $_pos_count		= 0;
	my $_char			= "";
	
	my $_message = $_[0];
	my $_style = $_[1];
	
	$_message_length = length($_message);
	
	if ($_style eq ""){$_char = "#"}
	if ($_style eq "input"){$_char = "~"}
	if ($_style eq "message"){$_char = "#"}
	if ($_style eq "warning"){$_char = "!"}
	
	print "\n\n";
	print color ' bold yellow';
	if ($_style eq "warning"){print color ' bold red'}
	if ($_style eq "input"){print color ' bold white'}
	
	for ($_pos_count = 1;$_pos_count <=($_message_length + 10);$_pos_count++){print $_char}
	
	print "\n$_char    $_message    $_char\n";
	
	for ($_pos_count = 1;$_pos_count <=($_message_length + 10);$_pos_count++){print $_char}
	
	print "\n\n";
	print color 'reset';

}#

#############################################
# Subroutine to execute unix command        #
#############################################

sub run_unix_command($;$)
{
	my $unix_command = "";
	my $step = "";	
	$unix_command = $_[0];
	$step = $_[1];
	print "\n";
	print("$unix_command\n");
	
	print COMMAND_LOG "$unix_command\n";

	my $returnCode = system("$unix_command");
	if ($returnCode != 0) {
		print COMMAND_LOG "run_unix_command detected an ERROR ($returnCode) when running command: $unix_command";
		&send_email_err_alert($returnCode,$unix_command);
	}
}

sub get_list_of_all_directories
{
	opendir (DIR, $local_storage_path);
	$dir_count = 0;
	$total_dir_count = 0;
	while ($file_in_local = readdir(DIR)) 
	{
		$full_file_path = "$local_storage_path"."$file_in_local";

		if (-e $full_file_path)
		{
			if ((-d $full_file_path)  && (index($full_file_path,".") == -1))
			{
				$total_dir_count = $total_dir_count + 1;
				$dir_array[$total_dir_count] = $full_file_path;
			}
		}# if file exists
	}

	$number_of_dirs_1 = $total_dir_count;
	close DIR;

	#######################################
	# Now look at sub-directories level 2 #
	#######################################
	for ($dir_count = 1; $dir_count <= $number_of_dirs_1; $dir_count++)
	{
		opendir (DIR2, $dir_array[$dir_count]);
		while ($file_in_local = readdir(DIR2)) 
		{
			$full_file_path = $dir_array[$dir_count]."/".$file_in_local;

			if (-e $full_file_path)
			{
					if ((-d $full_file_path) && (index ($file_in_local,".") != 0))
					{
						$total_dir_count = $total_dir_count + 1;
						$dir_array[$total_dir_count] = $full_file_path;
					}
			}
		}
	}
	close DIR2;
	$number_of_dirs_2 = $total_dir_count;

	#######################################
	# Now look at sub-directories level 3 #
	#######################################
	for ($dir_count = $number_of_dirs_1 + 1; $dir_count <= $number_of_dirs_2; $dir_count++)
	{
		opendir (DIR3, $dir_array[$dir_count]);
		while ($file_in_local = readdir(DIR3)) 
		{
			$full_file_path = $dir_array[$dir_count]."/".$file_in_local;
			if (-e $full_file_path)
			{
					if ((-d $full_file_path) && (index ($file_in_local,".") != 0))
					{
						$total_dir_count = $total_dir_count + 1;
						$dir_array[$total_dir_count] = $full_file_path;
					}
			}
		}
	}

	close DIR3;
} # get_list_of_all_directories