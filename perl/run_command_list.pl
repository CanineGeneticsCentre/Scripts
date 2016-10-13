#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename;
use Term::ANSIColor;
use Cwd;

# VERSION OF SOFTWARE #
my $version						= "3";

#File names
my $list_file					= "";
my $command_log					= "run_command_list_command_log.out";

#Strings
my $command						= "";
my $single_line					= "";
my $email_address				= "";
my $answer						= "";
my $from 						= 'command_list@samba64.aht.org.uk'; # Who e-mails come from

#Numbers
my $list_count					= 0;
my $no_of_lines					= 0;

#Timers
my $start_time 					= "";
my $end_time 					= "";
my $run_time					= "";
my $start_stage_time 			= "";
my $end_stage_time 				= "";
my $stage_time 					= "";

print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      run_command_list       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';

if ( index($version,"w") == -1 ){ print "Version $version (for samba64)\n\n"; } else { print "Version $version (for Workstation)\n\n"; }

print "  - This program runs a list of commands provided in a text file\n";
print "  - and then e-mails you when it's finished\n\n";


######################################
# Gte file with list of commands     #
######################################
until (-e "$list_file")
{
	&print_message("Please input the name of your file with a list of unix commands","input");
	
	$list_file = <STDIN>;chomp $list_file;
	
	if ($list_file eq "ls"){print "\n";system ("ls *.txt"); print "\n";}
}


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


#############################################
# Make sure the list file is in Unix format #
#############################################

$command = "dos2unix $list_file";
system("$command");
print "\n\n";


####################################################
# Open the list file to SHOW the list of commands  #
####################################################
&print_message("Here are the commands in this file","message");

open (LIST, "$list_file") || die "Cannot open $list_file";
while ($single_line = <LIST> ) 
{
	chomp $single_line;
	print "$single_line\n";
}

close LIST;

print "\nPress 'Return' if these commands are OK, or press 'Q' to quit\n\n";
$answer=<STDIN>;
chomp $answer;

if (lc $answer eq "q"){exit;}

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

$start_time = time();

####################################################
# Open the list file to get the list of commands   #
####################################################
open (LIST, "$list_file") || die "Cannot open $list_file";
$list_count=1;
$no_of_lines=0;
while ($single_line = <LIST> ) 
{
	chomp $single_line;

	print "\nCommand > $single_line\n\n";

	$start_stage_time = time();

	system("$single_line");

	$end_stage_time = time();
	$stage_time = $end_stage_time - $start_stage_time;
	printf "  Time for stage:\t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $stage_time)[7,2,1,0];

	print COMMAND_LOG "\nCommand > $single_line\n\n";
	printf COMMAND_LOG "  Time for stage:\t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $stage_time)[7,2,1,0];
	print COMMAND_LOG "============================================================\n\n";
}

$end_time = time();
$run_time = $end_time - $start_time;

printf "  Time for all commands:\t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
printf COMMAND_LOG "  Time for all commands:\t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
print COMMAND_LOG "============================================================\n";
print COMMAND_LOG "============================================================\n\n";
close LIST;
close COMMAND_LOG;


##########################################
# Send an e-mail at the end of each loop #
##########################################
open(MAIL, "|/usr/sbin/sendmail -t");

	## Mail Header
	print MAIL "To: $email_address\n";
	print MAIL "From: command_list\n";
	print MAIL "Subject: Your list of unix commands has been completed (File: $list_file)\n\n";
	## Mail Body
	print MAIL "Command list file:     \t$list_file\n\n";
	print MAIL "Commands:\n\n";

	open (LIST, "$list_file") || die "Cannot open $list_file";
	while ($single_line = <LIST> ) 
	{
		chomp $single_line;
		print MAIL "  $single_line\n";
	}

	close LIST;


close(MAIL);

&print_message("List of commands completed","message");

print "List of commands in $list_file completed\n\n";

print "An e-mail has been sent to $email_address\n\n";

print "There is a command_log file with times for each stage: $command_log\n\n";
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