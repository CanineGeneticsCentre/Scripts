 #!/usr/bin/perl -w

############################################################################
#									                                       #      
#	Compress gVCF files               		                               #     
#									                                       #
#	Compresses gVCF files in the correct way for GATK gVCFs			   	   #
#									                                       #
############################################################################

#############################
# Mike Boursnell Nov 2015   #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use Getopt::Long ;
use File::Basename ;
use Term::ANSIColor;

#Constants
my $version						= "1";

# File names
my $gVCF_file					= "";
my $gVCF_file_compressed		= ""; # 
my $gVCF_file_index				= "";
my $command_log					= "";

#Strings
my $prefix						= "";
my $command						= "";
my $answer						= "";
my $start_time					= "";
my $end_time					= "";
my $run_time					= "";

###############################
# Process flags               #
###############################

GetOptions("file:s"=>\$gVCF_file);

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "                 compress_gVCF_files        \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "    - This program compresses GATK gVCF files,\n";
print "      using bgzip and tabix in the correct way.\n\n";

print color 'reset';

if ($gVCF_file eq "")
{
	$gVCF_file = &get_file("What is the name of the gVCF file?","gVCF");
}

#############################
# Make up output file names #
#############################
$prefix = &get_prefix("$gVCF_file");
$gVCF_file_compressed = $gVCF_file.".gz";
$gVCF_file_index = $gVCF_file_compressed.".tbi";

$command_log = $prefix."_compress_gVCF_files_command_log.out";

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "COMMAND LOG for comparess_gVCF_files version $version\n\n";


print COMMAND_LOG "Input gVCFfile:             \t$gVCF_file\n\n";
print COMMAND_LOG "Output compressed file:     \t$gVCF_file_compressed\n\n";
print COMMAND_LOG "Output index file:          \t$gVCF_file_index\n\n";

$start_time = time();

##############################
# First compress with bgzip #
#############################
&print_message("Compressing $gVCF_file with bgzip","message");

$command = "bgzip $gVCF_file";

&run_unix_command_single("bgzip $gVCF_file");

##############################
# First compress with bgzip #
#############################
&print_message("Creating index of $gVCF_file_compressed with tabix","message");

&run_unix_command_single("tabix -p vcf $gVCF_file_compressed");


$end_time = time();
$run_time = $end_time - $start_time;

&print_message("Compression of $gVCF_file finished","message");

print  "InputgVCF file:            \t$gVCF_file\n";
print  "Output compressed file:    \t$gVCF_file_compressed\n";
print  "Output index file:         \t$gVCF_file_index\n\n";

print COMMAND_LOG "InputgVCF file:            \t$gVCF_file\n";
print COMMAND_LOG "Output compressed file:    \t$gVCF_file_compressed\n";
print COMMAND_LOG "Output index file:         \t$gVCF_file_index\n\n";

printf COMMAND_LOG "\nRun time: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];

close COMMAND_LOG;

exit;

sub get_prefix
{
	my $_filename 	= "";
	my $_dot_pos	= 0;

	$_filename = $_[0];
	$_dot_pos = rindex($_filename,".");

	if (rindex($_filename,".") > 0)
	{
		$_filename = substr($_filename, 0, rindex($_filename,"."));
	}
	if (rindex($_filename,".") == -1)
	{
		$_filename = $_filename;
	}

	$_filename = $_filename;

} # get_prefix



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


#################################################################
# Subroutine to get a file name                                 #
# Argument 1: Text that user sees, asking for the file          #
# Argument 2: suffix of files to search for with ls e.g. ".bed" #
#################################################################
sub get_file
{
	my $input_message 	= $_[0];
	my $_file_type	 	= $_[1];
	my $_file_name		= "";
	my $_search_string	= "";

	if ($_file_type eq ""){$_file_type = "*"}

	until ((-e $_file_name) || (lc $_file_name eq "q"))
	{
		&print_message("$input_message","input");
		print "> ";

		$_file_name = <STDIN>;
		chomp $_file_name;

		# User types 'ls'
		if ($_file_name eq "ls") {print "\n";system ("ls *"."$_file_type")}

		# Starts with 'ls' followed by search string
		if (($_file_name ne "ls") && (index ($_file_name,"ls") == 0))
		{
			$_search_string = substr($_file_name,3,99);
			print "\n";
			system ("ls *$_search_string*");
			$_file_name = "ls";
			print "\n";
		}

		if (($_file_name ne "ls")  && (lc $_file_name ne "q") && ($_file_name ne ""))
		{

			if (!-e $_file_name){print "\n  >>>>>>>>>>>>>  ERROR.  File $_file_name cannot be found <<<<<<<<<<<\n\n";}

			if (-e $_file_name)
			{
				$_file_name = $_file_name;
			}
		} # not ls

	} # until loop

	print "\n";
	if ($_file_type eq ".txt"){system("dos2unix $_file_name")}

	$_file_name = $_file_name;
} # get_file

#############################################
# Subroutine to execute unix command        #
# Argument 1: command to execute            #
#############################################
sub run_unix_command_single
{
	my $unix_command = "";
	$unix_command = $_[0];
		
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	print COMMAND_LOG "$unix_command\n";

	print("$unix_command\n\n");
	system("$unix_command");
}
