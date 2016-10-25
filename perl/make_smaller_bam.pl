	#!/usr/bin/perl

#############################
# Mike Boursnell            #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# mike.boursnell@aht.org.uk #
#############################
use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;


my $version						= "1";
my $mem							= "-Xmx4g";
my $bam_file					= "";
my $bai_file					= "";
my $prefix						= "";
my $region						= "";
my $bam_file_region				= "";
my $bai_file_region				= "";
my $answer						= "";

print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "       make_smaller_bam       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This program makes a smaller BAM file focussed on a chosen region\n\n";


print color 'bold cyan';

print "    (The input files are a BAM and its BAI file)\n\n";

print color 'reset';

until (-e $bam_file)
{
	&print_message("What is the prefix for your input BAM file:  ","input");

	$prefix = <STDIN>;
	chomp $prefix;

	$bam_file = $prefix.".bam";
	$bai_file = $prefix.".bai";

	if ($prefix eq "ls"){print "\n";system ("ls *.bam")}
			
	if ($prefix ne "ls")
	{
		if (! -e $bam_file){print "\n\n>>>>>>>>  File $bam_file not found.  Try again.  <<<<<<<<\n";}
		if (! -e $bai_file){print "\n\n>>>>>>>>  File $bai_file not found.  Try again.  <<<<<<<<\n\n";}
	}

	
}

$bam_file_region = $prefix."_region.bam";
$bai_file_region = $prefix."_region.bai";

##################################
# Define your region of analysis #
##################################
&print_message("Please define your region of interest (eg 'chr5:21000000-23000000')","input");
$region = <STDIN>;
chomp $region;


&run_unix_command("/opt/samtools/samtools view $bam_file $region -b -o $bam_file_region","Make smaller BAM file to specified region");



&record_input_file_size ("$bam_file");
&record_output_file_size ("$bam_file_region");	

&print_message("Finished making Bam file for region $region","message");

############################################################################
# FASTQ2VCF Step 15:  Now make an index file for this new smaller BAM file #
############################################################################
&print_message("New Index for region-only BAM file being created","message");

&run_unix_command("java $mem -jar /opt/picard/BuildBamIndex.jar I=$bam_file_region O=$bai_file_region VALIDATION_STRINGENCY=LENIENT","Make BAM index");

&record_input_file_size ("$bam_file_region");
&record_output_file_size ("$bai_file_region");

&print_message("New Index for region-only BAM file has been created","message");
		
&record_output_file_size ("$bai_file_region");	


&print_message("New sammler BAM file has been created","message");


print "Input BAM file:  \t$bam_file\n\n";
print "Output BAM file: \t$bam_file_region\n\n";
exit;


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

##############################################
# Subroutine to record size of output file   #
# (without sending an e-mail if it is zero)  #
##############################################

sub record_output_file_size
{
	my $outputfile = "";
	my $filesize = "";	

	$outputfile = $_[0];
	
	if (-e $outputfile)
	{
		$filesize = -s "$outputfile";
		print COMMAND_LOG "\n  Output file: $outputfile\t\tSize: $filesize\n";
		print "\n  Output file: $outputfile\t\tSize: $filesize\n\n"
	}
	else
	{
		$filesize=0;
		print COMMAND_LOG "\n  Output file: $outputfile\t\tSize: Not found\n";
		print "\n  Output file: $outputfile\t\tSize: Not found\n\n"
	}


}

##############################################
# Subroutine to record size of INPUT file    #
# (without sending an e-mail if it is zero)  #
##############################################
sub record_input_file_size
{
	my $inputfile = "";
	my $filesize = "";	

	$inputfile = $_[0];
	
	if (-e $inputfile)
	{
		$filesize = -s "$inputfile";
		print COMMAND_LOG "\n  Input file: $inputfile\t\tSize: $filesize\n";
		print "\n  Input file: $inputfile\t\tSize: $filesize\n\n"
	}
	else
	{
		$filesize=0;
		print COMMAND_LOG "\n  Input file: $inputfile\t\tSize: Not found\n";
		print "\n  Input file: $inputfile\t\tSize: Not found\n\n"
	}
}
	
#############################################
# Subroutine to execute unix command        #
#############################################

sub run_unix_command
{
	my $unix_command = "";

	$unix_command = $_[0];

	print "\n";
	print("$unix_command\n");
	system("$unix_command");
}