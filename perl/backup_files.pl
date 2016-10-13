#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	backup_files  	                                                    #
#   backs up files from the Workstation working area                    #
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
my $version						= "w9";

#File names
my $list_file					= "";
my $command_log					= "";
my $bai_file					= "";
my $bam_bai_file				= "";
my $file_to_find				= "";
my $file_in_local				= "";
my $full_file_path				= "";

#Directories
my $local_storage_path			= "/ahtlocal/local_storage/";


#Strings
my $command						= "";
my $file_to_be_backed_up		= "";
my $email_address				= "";
my $answer						= "";
my $from 						= 'backup_files@unix.aht.org.uk'; # Who e-mails come from
my $source_directory			= "";
my $destination_directory		= "/ahtlocal/local_storage"; # This is a constant but could be changed
my $read_file_method			= ""; # single or multiple
my $prefix						= "";

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
my $gzip_files					= "false"; # true or false
my $some_bam_files				= "false"; # if there are any BAM files
my $backup_bai					= "true";

#Date
my $date						= 'date';

my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900; 



print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      backup_files       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';

if ( index($version,"w") == -1 ){ print "Version $version (for samba64)\n\n"; } else { print "Version $version (for Workstation)\n\n"; }

print "  - This program backs up a list of files provided in a text file\n";
print "  - and then e-mails you when it's finished\n\n";

print "  - The files are COPIED not MOVED so no files are deleted.\n\n";

print color 'reset';

$command_log = "backup_log"."_"."$hour"."."."$min"."_"."$mday"."_"."$mon"."_"."$year".".out";
$source_directory=getcwd;


####################################
# Get directories in local_storage #
####################################
&get_list_of_all_directories;


######################################
# E-mail address for notifications   #
######################################

&print_message("Please enter your email address to be notified when files are backed up","input");
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

	print "   <1> Add the names of the files one by one           [DEFAULT]\n";
	print "   <2> Using a file of file names for your files to be backed up\n\n";
	

	$answer = <STDIN>;chomp $answer;
	if ($answer eq ""){$answer = "1"} # DEFAULT
	if (substr($answer,0,1) eq "1"){$read_file_method = "single"}
	if (substr($answer,0,1) eq "2"){$read_file_method = "multiple"}
	
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
		&print_message("Please input the name of your file with a list of files to be backed up to $destination_directory","input");
		
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
	while ($file_to_be_backed_up = <LIST> ) 
	{
		chomp $file_to_be_backed_up;

		if (-e "$file_to_be_backed_up")
		{
			$file_count = $file_count + 1;
			$file_exist_count = $file_exist_count + 1;
			$file_array[$file_count] = $file_to_be_backed_up;
			
			print "File $file_to_be_backed_up     found\n";

			#############################
			# Check if it is a BAM file #
			#############################
			if (index ($file_to_be_backed_up,"bam") > -1){$some_bam_files = "true"}
		}
		if (! -e "$file_to_be_backed_up")
		{
			print "File $file_to_be_backed_up CANNOT be found\n";
			
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

	&print_message("Please input the name of your file(s) to be backed up to $destination_directory (q to end list)","input");

	until (lc $file_to_be_backed_up eq "q")
	{
		
		print "File $file_count:    ";
		$file_to_be_backed_up = <STDIN>;
		chomp $file_to_be_backed_up;
		print "\n";

		if ($file_to_be_backed_up eq "ls"){print "\n";system ("ls"); print "\n";}

		if (($file_to_be_backed_up ne "ls") && (lc $file_to_be_backed_up ne "q"))
		{
			if (! -e "$file_to_be_backed_up")
			{
				print "File $file_to_be_backed_up cannot be found\n\n";
			}
			if (-e "$file_to_be_backed_up")
			{
				$file_array[$file_count] = $file_to_be_backed_up;
				$file_count = $file_count + 1;
				$file_exist_count = $file_exist_count + 1;
			}
		}
	}
	
	$total_no_files = $file_exist_count;
	
} # if read_file_method eq single



print "\nThere are $file_exist_count files\n\n";


####################################################
# Ask user if the files should be gzipped          #
####################################################
&print_message("Do you want the files gzipped when they get to the destination folder?","input");

print "   <1> YES - gzip the files\n";
print "   <2> NO                               [DEFAULT]\n\n";
$answer=<STDIN>;
if (substr($answer,0,1) eq "1" ){$gzip_files = "true"} else {$gzip_files = "false"}


####################################################
# Ask user which directory to use                  #
####################################################
&print_message("Which directory do you want to use for backup?","input");

for ($dir_count = 1; $dir_count <= $total_dir_count; $dir_count++)
{
		print "<$dir_count>\t$dir_array[$dir_count]\n";
}

print "\n>  ";
$answer=<STDIN>;
chomp $answer;

$destination_directory = $dir_array[$answer];
if (-e $destination_directory)
{
	#print "$destination_directory exists\n";
}
else
{
	print "\n\nERROR: $destination_directory can't be found\n\n";
	exit;
}

if ($some_bam_files eq "true")
{
	####################################################
	# Ask user if the files should be gzipped          #
	####################################################
	&print_message("Do you want any BAI files backed up along with the BAM?","input");

	print "   <1> YES - back up BAI as well         [DEFAULT]\n";
	print "   <2> NO\n\n";
	$answer=<STDIN>;

	if (substr($answer,0,1) eq "2" ){$backup_bai = "false"} else {$backup_bai = "true"}
}

####################################################
# Open the list file to SHOW the list of files     #
####################################################
&print_message("List of files to be backed-up to $destination_directory","message");


for ($file_count = 1; $file_count <=$total_no_files; $file_count++)
{
	print "File $file_count: $file_array[$file_count]\n";
}


print "\nPress 'Return' if these files are OK, or press 'Q' to quit\n\n";
$answer=<STDIN>;
chomp $answer;

if (lc $answer eq "q"){exit;}

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

print COMMAND_LOG "Command log of perl script backup_files, version $version\n\n";

print COMMAND_LOG "List of files backed up\n\n";




####################################################
# Open the list file to get the list of files      #
####################################################
$no_of_lines=0;
for ($file_count = 1; $file_count <=$total_no_files; $file_count++)
{
	$file_to_be_backed_up = $file_array[$file_count];

	if (-e "$file_to_be_backed_up")
	{
		$command = "cp $source_directory/$file_to_be_backed_up $destination_directory/$file_to_be_backed_up";
		print "$command\n";
		
		print COMMAND_LOG "FILE:    $file_to_be_backed_up\n\n";
		print COMMAND_LOG "COMMAND: $command\n\n";
		system ("$command");
		
		#
		# Look for BAI #
		################
		if ($backup_bai eq "truexx")
		{
			$prefix = &get_prefix($file_to_be_backed_up);
			$bai_file = $prefix.".bai";
			$bam_bai_file = $prefix.".bam.bai";

			if (-e $bai_file)
			{
				$command = "cp $source_directory/$bai_file $destination_directory/$bai_file";
				print "$command\n";
				
				print COMMAND_LOG "FILE:    $bai_file\n\n";
				print COMMAND_LOG "COMMAND: $command\n\n";
				system ("$command");
			}
			if (-e $bam_bai_file)
			{
				$command = "cp $source_directory/$bam_bai_file $destination_directory/$bam_bai_file";
				print "$command\n";
				
				print COMMAND_LOG "FILE:    $bam_bai_file\n\n";
				print COMMAND_LOG "COMMAND: $command\n\n";
				system ("$command");
			}
		}
		#######################################
		# Write file sizes to the command log #
		#######################################

		$file_size_source = -s "$source_directory/$file_to_be_backed_up";
		$file_size_destination = -s "$destination_directory/$file_to_be_backed_up";

		print COMMAND_LOG "Check on file sizes:\n\n";

		print COMMAND_LOG "  File: $file_to_be_backed_up\n";
		print COMMAND_LOG "  Source file size:      \t$file_size_source\n";
		print COMMAND_LOG "  Destination file size: \t$file_size_destination\n\n";


		if ($gzip_files eq "true")
		{
			$command = "gzip $destination_directory/$file_to_be_backed_up";
			print "$command\n\ns";
			print COMMAND_LOG "GZIP:    $command\n\n";
			system ("$command");
		}

		print COMMAND_LOG "========================================================================\n\n";
	}
}


##########################################
# Send an e-mail when finished           #
##########################################
open(MAIL, "|/usr/sbin/sendmail -t");
 
## Mail Header
print MAIL "To: $email_address\n";
print MAIL "From: backup_files\n";
print MAIL "Subject: Back up of files in $list_file to $destination_directory\n\n";
## Mail Body
print MAIL "Program: backup_files\tVersion $version\n\n";

print MAIL "This program running the backup of your files has finished\n\n";
print MAIL "Check all the files have been copied correctly before deleting the originals\n\n";

if ($read_file_method eq "multiple"){print MAIL "File with list of files to be backed up:     \t$list_file\n\n";}

print MAIL "Destination directory: $destination_directory\n\n";

print MAIL "Files backed up:\n\n";

for ($file_count = 1; $file_count <=$total_no_files; $file_count++)
{
	print MAIL "\t$file_array[$file_count]\n";
}

close MAIL;

######################################################
# Copy command log to command log folder in genetics #
######################################################

&run_unix_command("cp $command_log /home/genetics/command_logs/$command_log","Copy backup_files command log to genetics folder");

close COMMAND_LOG;

&print_message("Backup program completed","message");

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

###########################################################
# Subroutine to get filename before file type suffix      #
# (i.e. prefix) (e.g. to get "test" from "test.fasta")    #
###########################################################

sub get_prefix
{
	my $_filename = "";

	$_filename = $_[0];
	
	if (index($_filename,".") > 0)
	{
		$_filename = substr($_filename, 0, index($_filename,"."));
	}
	if (index($_filename,".") == -1)
	{
		$_filename = $_filename;
	}
} # get_prefix


###########################################################
# Subroutine to get file suffix (e.g. '.txt', '.bam')     #
###########################################################

sub get_suffix
{
	my $_filename = "";

	$_filename = $_[0];
	
	if (index($_filename,".") > 0)
	{
		$_filename = substr($_filename,index($_filename,".") + 1, 99);
	}
	if (index($_filename,".") == -1)
	{
		$_filename = $_filename;
	}

	print "Suffix: $_filename\n";

} # get_suffix

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