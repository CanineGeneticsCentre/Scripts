#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	delete_files  	                                                    #
#   deletes files from the Workstation working area                     #
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
my $version						= "w5";

#File names
my $list_file					= "";
my $command_log					= "";

#Strings
my $command						= "";
my $file_to_be_deleted		= "";
my $email_address				= "";
my $answer						= "";
my $from 						= 'delete_files@unix.aht.org.uk'; # Who e-mails come from
my $source_directory			= "";
my $local_storage_directory		= "/ahtlocal/local_storage"; # This is a constant but could be changed
my $read_file_method			= ""; # single or multiple
my $padding						= "";

#Numbers
my $no_of_lines					= 0;
my $file_count					= 0;
my $file_exist_count			= 0;
my $total_no_files				= 0;
my $file_size_source 			= 0;
my $file_size_destination 		= 0;

#Arrays
my @file_array					= ();
my @file_found_in_storage_array = (); # if file has been backed up

#Boolean
my $gzip_files					= "false"; # true or false
my $file_found					= "false";

#Date
my $date						= 'date';

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900; 



print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      delete_files       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';

if ( index($version,"w") == -1 ){ print "Version $version (for samba64)\n\n"; } else { print "Version $version (for Workstation)\n\n"; }

print "  - This program deletes a list of files provided in a text file\n";
print "  - and then e-mails you when it's finished\n\n";



$command_log = "delete_log"."_"."$hour"."."."$min"."_"."$mday"."_"."$mon"."_"."$year".".out";
$source_directory=getcwd;

print "Command log file: $command_log\n\n";

######################################
# E-mail address for notifications   #
######################################

&print_message("Please enter your email address to be notified when files are deleted","input");
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


#####################################################################
# ASK IF YOU WANT TO READ THE FILENAMES FROM A "FILE OF FILE NAMES" #
#####################################################################

$answer = "";

until ($read_file_method eq "multiple" || $read_file_method eq "single")
{
		&print_message("How do you want to read the input files?","input");

	print "   <1> Using a file of file names for your files to be deleted\n";
	print "   <2> Add the names of the files one by one\n\n";

	$answer = <STDIN>;chomp $answer;
	if (substr($answer,0,1) eq "1"){$read_file_method = "multiple"}
	if (substr($answer,0,1) eq "2"){$read_file_method = "single"}
}

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
		&print_message("Please input the name of your file with a list of files to be deleted","input");
		
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
	$file_count =0;
	$file_exist_count = 0;
	while ($file_to_be_deleted = <LIST> ) 
	{
		chomp $file_to_be_deleted;

		if (-e "$file_to_be_deleted")
		{
			$file_count = $file_count + 1;
			$file_exist_count = $file_exist_count + 1;
			$file_array[$file_count] = $file_to_be_deleted;
			
			print "File $file_to_be_deleted found OK\n";
		}
		if (! -e "$file_to_be_deleted")
		{
			print "File $file_to_be_deleted cannot be found\n";
			
		}
	}
	
	$total_no_files = $file_exist_count;
	
} # if read_file_method eq multiple


######################################
# Get files one by one               #
######################################
if ($read_file_method eq "single")
{
	$file_count =  1;
	$file_exist_count =  0;

	&print_message("Please input the name of your file(s) to be backed up to $local_storage_directory (q to end list)","input");

	until (lc $file_to_be_deleted eq "q")
	{
		
		print "File $file_count:    ";
		$file_to_be_deleted = <STDIN>;
		chomp $file_to_be_deleted;
		print "\n";

		if ($file_to_be_deleted eq "ls"){print "\n";system ("ls"); print "\n";}

		if (($file_to_be_deleted ne "ls") && (lc $file_to_be_deleted ne "q"))
		{
			if (! -e "$file_to_be_deleted")
			{
				print "File $file_to_be_deleted cannot be found\n\n";
			}
			if (-e "$file_to_be_deleted")
			{
				$file_array[$file_count] = $file_to_be_deleted;
				$file_count = $file_count + 1;
				$file_exist_count = $file_exist_count + 1;
			}
		}
	}
	
	$total_no_files = $file_exist_count;
	
} # if read_file_method eq single



print "\nThere are $file_exist_count files\n\n";



####################################################
# Open the list file to SHOW the list of files     #
####################################################
&print_message("List of files to be deleted","message");


for ($file_count = 1; $file_count <=$total_no_files; $file_count++)
{
	print "File $file_count: $file_array[$file_count]";
	$file_found = "false";

	# Look in local_storage
	if (-e "$local_storage_directory/$file_array[$file_count]"){$file_found = "true"}

	if (-e "$local_storage_directory/bam_files/$file_array[$file_count]"){$file_found = "true"}

	if (-e "$local_storage_directory/fastq_files/$file_array[$file_count]"){$file_found = "true"}

	if (-e "$local_storage_directory/fastq_files/original_fastq_files/$file_array[$file_count]"){$file_found = "true"}

	if (-e "$local_storage_directory/fastq_files/merged_fastq_files/$file_array[$file_count]"){$file_found = "true"}

	$padding = &add_space($file_array[$file_count]);

	if ($file_found eq "true")
	{
		print "$padding file found in local_storage\n";
		$file_found_in_storage_array[$file_count] = "Found";
	}
	else
	{
		print "$padding file not found in local_storage so will NOT be deleted\n";
		$file_found_in_storage_array[$file_count] = "Not found";
	}
} # show user list of files and whether they are on local_storage


print "\nPress 'Y' if these files are OK, or press 'Q' to quit\n\n";
$answer=<STDIN>;
chomp $answer;

if (lc $answer ne "y"){exit;}


###############################
# Add this to the COMMAND_LOG #
###############################
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "Command log of perl script delete_files, version $version\n\n";
print COMMAND_LOG "List of files to be deleted\n\n";


for ($file_count = 1; $file_count <=$total_no_files; $file_count++)
{
	$padding = &add_space($file_array[$file_count]);

	if ($file_found_in_storage_array[$file_count] eq "Found")
	{
		print COMMAND_LOG "File $file_count: $file_array[$file_count] $padding will be deleted\n";
	}
	else
	{
		print COMMAND_LOG "File $file_count: $file_array[$file_count] $padding will NOT be deleted\n";
	}

}
print COMMAND_LOG "\n\n";


####################################################
# Now delete the relevant files                    #
####################################################
$no_of_lines=0;
&print_message("Deleting files...","message");

for ($file_count = 1; $file_count <=$total_no_files; $file_count++)
{
	$file_to_be_deleted = $file_array[$file_count];

	if (-e "$file_to_be_deleted")
	{
		if ($file_found_in_storage_array[$file_count] eq "Found")
		{

			$command = "rm -I $source_directory/$file_to_be_deleted";
			print COMMAND_LOG ("$command\n");

			system ("$command");
			&print_both ("   File $file_count: $file_array[$file_count] $padding was deleted\n\n");
		}
		else
		{
			&print_both ("   File $file_count: $file_array[$file_count] $padding was NOT  deleted\n\n");
		}	

	} # if file exists
	else
	{
		&print_both ("File $file_count: $file_array[$file_count] $padding was not found\n\n");
	}
}


##########################################
# Send an e-mail when finished           #
##########################################
open(MAIL, "|/usr/sbin/sendmail -t");
 
## Mail Header
print MAIL "To: $email_address\n";
print MAIL "From: delete_files\n";
print MAIL "Subject: Deleting of files in $list_file\n\n";

## Mail Body
print MAIL "Program: delete_files\tVersion $version\n\n";

print MAIL "This program running the deleting of your files has finished\n\n";
print MAIL "(Only files that wer found in local_stoarge were deleted)\n\n";

if ($read_file_method eq "multiple"){print MAIL "File with list of files to be backed up:     \t$list_file\n\n";}

print MAIL "Files to be backed up:\n\n";

for ($file_count = 1; $file_count <=$total_no_files; $file_count++)
{
	print MAIL "\t$file_array[$file_count]";
	$padding = &add_space($file_array[$file_count]);

	if ($file_found_in_storage_array[$file_count] eq "Found")
	{
		print MAIL "$padding File was deleted\n";
	}
	else
	{
		print MAIL "$padding File was NOT deleted\n";
	}
}

close MAIL;

######################################################
# Copy command log to command log folder in genetics #
######################################################

&run_unix_command("cp $command_log /home/genetics/command_logs/$command_log","Copy delete_files command log to genetics folder");

close COMMAND_LOG;

&print_message("Delete_files program completed","message");

print "All files have been deleted (but only if they were backed up to local_storage)\n\n";

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


sub add_space
{
	my $string = $_[0];
	my $space = "";
	my $count = 0;

	for($count=1; $count < 40-length($string) ; $count++)
	{
		$space = $space." ";
	}

	$space = $space;
}


####################################################
# Subroutine to print to screen and to COMMAND_LOG #
####################################################
sub print_both
{
	my $_message = $_[0];

	print "$_message";
	print COMMAND_LOG "$_message";
}