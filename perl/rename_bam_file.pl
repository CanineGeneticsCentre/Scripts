#!/usr/bin/perl -w

##################################################################################
#									                                             #      
#	rename_bam		                                                             #     
#									                                             #
#	This PERL script changes the name of a BAM file (and records it in a log)    #
#    It can also changes the name imbedded in the headers of the BAM File	     #
#									                                             #
##################################################################################

use strict;
use Term::ANSIColor;
use List::Util qw[min max];

#Constants
my $version						= "1";
my $picard_validation_stringency = "STRICT";
my $memory_string				= "-Xmx35g";

#Strings
my $vcf_line					= "";
my $prefix						= "";
my $prefix_new					= "";
my $answer						= "";
my $old_name					= "";
my $new_name					= "";
my $chromosome					= "";
my $position					= "";
my $tempdir						= "";
my $temp_dir_string				= "";
my $run_title					= "";
my $lib							= "";
my $sample_name					= "";
my $header_line					= "";
my $PU_string					= "";
my $SM_string					= "";
my $SM_string_new				= "";
my $LB_string					= "";
my $ID_string					= "";
my $datestring					= "";

#File names
my $temp_header_file			= "";
my $bam_file					= "";
my $bam_output_file				= "";
my $command_log					= "";

#Numbers
my $no_of_samples				= 0;
my $no_of_samples_chrom_line	= 0;
my $chrom_line_array_size		= 0;
my $pos_PU						= 0;
my $pos_SM						= 0;
my $pos_ID						= 0;
my $pos_LB						= 0;
my $space_pos					= 0;
my $input_file_size				= 0;
my $output_file_size			= 0;

#Counters
my $array_count					= 0;
my $line_count					= 0;
my $column_count				= 0;

my $sample_count				= 0;

#Boolean
my $passed_header_lines			= ""; # true or false

#Arrays
my @chrom_line_array			= ();
my @vcf_line_array				= ();
my @sample_name_array			= ();
my @sample_name_array_new		= ();

#############################################################
print color 'reset';
print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      rename BAM            \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";


print color 'yellow';
print "Version $version\n\n";
print "  - This PERL script changes the name of a BAM file, and [importantly]\n";
print "    records this change in a log file.\n\n";
print "  - It will also change the name imbedded in the BAM header\n\n";

###########################################################
# Get the name of the input VCF file                      #
###########################################################
$bam_file = &get_file ("What is the name of the input BAM file","bam");

$prefix = &get_prefix($bam_file);
$temp_header_file = $prefix."_temp_header.txt";

&print_message("What name would you like for the renamed file","input");

$bam_output_file = <STDIN>;
chomp $bam_output_file;

if (index($bam_output_file,".bam") == -1)
{
	$bam_output_file = $bam_output_file.".bam"
}

&print_message("Please check this name change","message");

print "Existing BAM file name:  $bam_file\n\n";
print "New BAM file name:       $bam_output_file\n\n";

print "Continue with renaming? (y/n)   ";
$answer=<STDIN>;
chomp $answer;
if (lc $answer ne "y"){exit;}

$prefix_new = &get_prefix($bam_output_file);

&run_unix_command_single("samtools view -H $bam_file > $temp_header_file");

&run_unix_command_single("cat $temp_header_file");

###########################################################
# Open the temporary header file                          #
# To get the Read Groups information                      #
# @RG     ID:V1   PL:ILLUMINA   PU:V_23005  LB:canfam3  SM:V_23005
###########################################################
open (TEMP_HEADER, "$temp_header_file") || die "Cannot open $temp_header_file";
while ($header_line = <TEMP_HEADER>) 
{
	chomp $header_line;
	if (index($header_line,"RG") == 1)
	{
		$pos_PU = index($header_line,"PU:");
		$pos_SM = index($header_line,"SM:");
		$pos_ID = index($header_line,"ID:");
		$pos_LB = index($header_line,"LB:");

		print "\n\n$header_line\n";

		if ($pos_ID > 0)
		{
			# Get first tab (or space) after ID
			$space_pos = index($header_line,"\t",$pos_ID);
			if ($space_pos == -1){$space_pos = index($header_line," ",$pos_ID);}

			if ($space_pos > -1){$ID_string = substr($header_line,$pos_ID + 3,($space_pos - $pos_ID - 3));}else{$ID_string = substr($header_line,$pos_ID + 3)}
		}

		if ($pos_LB > 0)
		{
			# Get first tab (or space) after LB
			$space_pos = index($header_line,"\t",$pos_LB);
			if ($space_pos == -1){$space_pos = index($header_line," ",$pos_LB);}

			if ($space_pos > -1){$LB_string = substr($header_line,$pos_LB + 3,($space_pos - $pos_LB - 3));}else{$LB_string = substr($header_line,$pos_LB + 3)}
		}

# @RG     ID:V1   PL:ILLUMINA   PU:V_23005  LB:canfam3  SM:V_23005

		if ($pos_PU > 0)
		{
			# Get first tab (or space) after PU
			$space_pos = index($header_line,"\t",$pos_PU);
			if ($space_pos == -1){$space_pos = index($header_line," ",$pos_PU);}

			if ($space_pos > -1){$PU_string = substr($header_line,$pos_PU + 3,($space_pos - $pos_PU - 3));}else{$PU_string = substr($header_line,$pos_PU + 3)}
		}

		if ($pos_SM > 0)
		{
			# Get first tab (or space) after SM
			$space_pos = index($header_line,"\t",$pos_SM);
			if ($space_pos == -1){$space_pos = index($header_line," ",$pos_SM);}

			if ($space_pos > -1){$SM_string = substr($header_line,$pos_SM + 3, ($space_pos - $pos_SM - 3));}else{$SM_string = substr($header_line,$pos_SM+3)}
		}

		&print_message("Naming in the BAM header lines","message");

		print "ID string: $ID_string\n";
		print "LB string: $LB_string\n\n";

		print "PU string: $PU_string\n";
		print "SM string: $SM_string\n\n";

		print "We usually give the two fields (PU and SM) the same name (ie the sample name).\n\n";
		print "What new name would you like? (Press 'return' to keep the same name)\n\n";

		$answer=<STDIN>;
		chomp $answer;
		if ($answer ne ""){$SM_string_new = $answer} else {$SM_string_new = $SM_string}

		&print_message("Please check this header name change","message");

		print "Existing header sample name:  $SM_string\n";
		print "New header sample name:       $SM_string_new\n\n";

		print "Continue with renaming? (y/n)   ";
		$answer=<STDIN>;
		chomp $answer;
		if (lc $answer ne "y"){exit;}
	}
}

close TEMP_HEADER;


####################
# Command Log name #
####################
$command_log = "BAM_name_change_".$prefix."_".$prefix_new."_command_log.out";



################################
# Make the temp java directory #
################################
$tempdir = "$ENV{HOME}/javatempdir"; 
$temp_dir_string = " -Djava.io.tmpdir=$tempdir";
if (! -e $tempdir)
{
	unless(mkdir $tempdir){die "Unable to create temporary Java directory $tempdir";}
	$temp_dir_string = " -Djava.io.tmpdir=$tempdir";	
}

$datestring = localtime();

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

print COMMAND_LOG "Log of BAM file name change\n\n";
print COMMAND_LOG "Date and time $datestring\n\n";
print COMMAND_LOG "Existing BAM file name:       \t$bam_file\n";
print COMMAND_LOG "New BAM file name:            \t$bam_output_file\n\n";
print COMMAND_LOG "Existing header sample name:  \t$SM_string\n";
print COMMAND_LOG "New header sample name:       \t$SM_string_new\n\n";


print COMMAND_LOG "Unix command used:\n\n";;

&run_unix_command_single("java $memory_string $temp_dir_string -jar /opt/picard/AddOrReplaceReadGroups.jar I=$bam_file O=$bam_output_file rgID=$ID_string LB=$LB_string PL='ILLUMINA' PU=$SM_string_new SM=$SM_string_new SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=$picard_validation_stringency");	
	
# removed rgID=$run_title LB=$lib
#&run_unix_command_single("java -Xmx45g  -Djava.io.tmpdir=/home/LANPARK/mboursnell/javatempdir -jar /opt/picard/AddOrReplaceReadGroups.jar I=V_21485_V_21485.fastq_aligned.sam O=V_21485_V_21485.fastq_aligned_sorted_rg.bam rgID=V_21485 LB=canfam3 PL='ILLUMINA' PU=V_21485.fastq SM=V_21485.fastq SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT");
&print_message ("Renaming completed","message");

$input_file_size = -s "$bam_file";
$output_file_size = -s "$bam_output_file";

print COMMAND_LOG "\n\n";
print COMMAND_LOG "Input file size:              \t$input_file_size\n";
print COMMAND_LOG "Output file size:             \t$output_file_size\n";
close COMMAND_LOG;

&run_unix_command_single("cp $command_log /home/genetics/command_logs/$command_log");

print "Input file:          \t$bam_file\n";
print "Output file renamed: \t$bam_output_file\n\n";

print "Input file size:      \t$input_file_size\n";
print "Output file size:     \t$output_file_size\n\n\n";

exit;


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
		if (($_file_name eq "ls") || ($_file_name eq "")) {print "\n";system ("ls *$_file_type");;}

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

###########################################################
# Subroutine to get filename before file type suffix      #
# (i.e. prefix) (e.g. to get "test" from "test.fasta")    #
# New version using rindex rather than index.  This means #
# it can deal with files like filename.something.vcf      #
###########################################################

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

#############################################
# Subroutine to execute unix command        #
# Argument 1: command to execute            #
#############################################
sub run_unix_command_single
{
	my $unix_command = "";
	$unix_command = $_[0];
		
	print COMMAND_LOG "$unix_command\n";

	print("$unix_command\n\n");
	system("$unix_command");
}
