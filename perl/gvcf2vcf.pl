#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	gvcf2vcf 						                                    #     
#									                                    #
#	Converts a group of gVCF files into a single multi-column VCF file  #
#									                                    #
#########################################################################

#############################
# Mike Boursnell Jan 2015   #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;
use Cwd;

my $version							= "w8";
my $testing_mode					= "off";

##################################################
#  GATK memory suggestions for variant calling   #
#------------------------------------------------#
#Tool						 UG                  #
#Available modes 	     NT,NCT,SG               #
#CPU threads (-nct)      3 / 6 / 24              #
#Data threads (-nt)      8 / 4 / 1               #
#Memory (Gb) 		    32 / 16 / 4              #
##################################################

############################################################
# Memory and threading settings (can be altered later on)  #
#----------------------------------------------------------#
#                              Workstation      samba64    #
#                                                          #
# Java memory                       60             4       #
#                                                          #
# $no_of_threads_gatk_ug_nct        1              1       #
# $no_of_threads_gatk_ug_nt         16             2       #
#                                                          #
############################################################


my $memory_in_Gb				= "4"; # memory setting for Java in gigabytes
my $no_of_threads_bwa			= "2";

my $no_of_threads_gatk_ug_nt	= "2";  # UnifiedGenotyper -nt (or HaplotypeCaller)

my $workstation 				= "unknown";
my $e_mail_from 				= ""; # Who e-mails come from

# Workstation version - (if 'w' in version name)
if (index($version,"w") > -1)
{
	$memory_in_Gb				= "45"; # memory setting for Java in gigabytes
	$no_of_threads_gatk_ug_nt	= "1";  # UnifiedGenotyper -nt (but NOT HaplotypeCaller)
	$workstation 				= "true";
	$e_mail_from 				= 'NGS_analysis@gen-x1404-ws01.aht.org.uk'; # Who e-mails come from
}
else
# Samba64 version (if no 'w' in version name)
{
	$memory_in_Gb				= "4"; # memory setting for Java in gigabytes
	$no_of_threads_gatk_ug_nt	= "1";  # UnifiedGenotyper -nt (but NOT HaplotypeCaller)
	$workstation 				= "false";
	$e_mail_from 				= 'NGS_analysis@samba64.org.uk'; # Who e-mails come from
}



#Various Parameters (e.g. for the Unified Genotyper)
my $gatk_directory					= "gatk";
my $GATK_validation_stringency		= "STRICT";
my $stand_emit_conf					= 30;
my $stand_call_conf					= 30;
my $min_reads_platypus				= 2; # Minimum number of reads required for a call in platypus (2 is the default)
my $max_alt_alleles					= 6; # This is how many alleles the INDEL caller in UnifiedGenotyper can allow

# Integers
my $list_count						= 0;
my $no_of_files						= 0;
my $dot_pos							= 0;
my $size							= 0;
my $chromosome_count				= 0;

# Boolean
my $use_default_stand_values		= "true"; # default values of stand_emit_conf and stand_call_conf
my $use_defined_region				= "";
my $run_chromosomes_separately		= "false";
my $has_threading_option			= "true"; # true or false. Can it use multiple CPU threads
my $has_memory_option				= "true"; # true or false. Can it use memory setting like -Xmx4g
my $choice_ok						= "false";
my $fix_misencoded_qual_scores		= ""; # 'yes' or 'no'

# File names
my $vcf_file						= ""; # output VCF file
my $list_file						= ""; # file of file names
my $gVCF_file						= "";
my $bai_file						= "";
my $gVCF_file_in_results				= "";
my $bai_file_in_results				= "";
my $output_vcf						= ""; # Output VCF file
my $output_gvcf						= ""; # Output GVCF from HaplotypeCaller GVCF method
my $output_SNPs_vcf					= ""; # Output VCF in separate SNPs
my $output_Indels_vcf				= ""; # Output VCF if separate Indels
my $command_log_file				= "";
my $screen_log_file					= ""; # log file for screen output

#Strings
my $quality_scores					= ""; # Can be "New" or "Old". Assigned by user during program
my $GATK_fix_quals_string			= "";
my $ref								= ""; # The reference sequence (e.g. canfam3)
my $temp_dir_string					= "";
my $tempdir							= "";
my $staticdata						= "";
my $input							= "";
my $memory_string					= "";

my $default_value					= "";
my $region							= "";
my $command							= "";
my $vcf_merge_string				= "";
my $prefix							= ""; # Part of sfile name before the dot
my $GATK_input_string				= ""; # string of -I bam_file_1 -I bam_file_2 for GATK
my $GVCF_input_string				= ""; # string of -V GVC_file_ -V GVCF_file_2 for GATK
my $run_title						= "";
my $chromosome						= "";
my $GATK_chromosome_only_string		= ""; #This is the chr only string for the UnifiedGenotyper command line -L chr12
my $GATK_region_string				= ""; # -L chr12:50000000-60000000
my $freeBayes_region_string			= ""; # --region chr12:50000000-60000000
my $samtools_region_string			= ""; # This is the region string in samtools format: chr12:50000000-60000000
my $platypus_region_string			= ""; # This is the region string in platypus format: --regions chr12:50000000-60000000
my $sample_name						= "";
my $cleaned_sorted_bam				= "";
my $answer							= "";
my $calls							= "";
my $gVCF_file_for_unified_genotyper 	= "";
my $species							= "";
my $ref_seq_name					= "";
my $other_ref_sequence 				= "false";
my $current_directory 				= "";
my $folder 							= "";
my $gVCF_file_in_folder				= "";
my $variant_caller					= "HaplotypeCallerGVCF"; # fixed for gVCF 2 VCF
my $search_string					= "";

#Timers
my $start_time						= "";
my $end_time						= "";
my $current_time					= "";

my @gVCF_file_array					= ();
my @chromosome_end_array			= (); # for running all chromosomes separately

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "              gvcf2vcf   \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';

if ( index ($version,"w") == -1 ){ print "Version $version\n\n"; } else { print "Version $version (for Workstation)\n\n"; }

print "    - This program runs GATK programs to convert multiple gVCF files to a single multi-column VCF file\n\n";

print color 'reset';


################################
# Make the temp java directory #
# Check if staticdata exists   #
################################
$tempdir = "$ENV{HOME}/javatempdir"; # Individual's HOME space (if staticdata isn't online)
$temp_dir_string = " -Djava.io.tmpdir=$tempdir";

&print_message("The input is a file with a list of the gVCF file names.","input");


##########################################################
# Get the file that contains a list of the gVCF filenames #
##########################################################
until (-e $list_file)
{
	print "Name of the file of file names:      ";
	$list_file = <STDIN>;
	chomp $list_file;
	
	if ($list_file eq ""){$list_file = "AS_HC_small_bam_input.txt"} # TEMP!!!!
	if ($list_file eq "ls")
	{
		print "\n";
		system ("ls *.txt");
		print "\n";
	}
	# Starts with 'ls'
	if (($list_file ne "ls") && (index ($list_file,"ls") == 0))
	{
		$search_string = substr($list_file,3,99);
		print "\n";
		system ("ls *$search_string*");
		$list_file = "ls";
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
print "\n";


########################################################
# Open the list file to get the list of gVCF file names #
########################################################
open (LIST, "$list_file") || die "Cannot open $list_file";
$list_count=1;

while ($gVCF_file = <LIST> ) 
{
	chomp $gVCF_file;
	
	$prefix = &get_prefix("$gVCF_file");
	
	print "$list_count\t$gVCF_file\t";
	
	if (-e $gVCF_file){print "File found\n"}

	if (! -e $gVCF_file)
	{
		print "File not found\n";
		#print "\n\n";
		#print "###########################################################################################\n";
		#print "$gVCF_file cannot be found\n\n";
		#print "Correct the input file of file names, or the locations of the input gVCF files and try again\n";
		#print "###########################################################################################\n\n";
		#close LIST;
		#exit;
	}
	
	$gVCF_file_array[$list_count]=$gVCF_file;
	$list_count=$list_count + 1;
}

$no_of_files = $list_count - 1;

close LIST;


##############################################
# Get name of output VCF file                #
##############################################
$vcf_file = "ls";
until ($vcf_file ne "ls")
{
	&print_message("Type a name for the output VCF file:","input");

	$vcf_file = <STDIN>; chomp $vcf_file;

	if ($vcf_file eq "ls")
	{
		print "\n";
		system ("ls *.vcf");
		print "\n";
	}
}

if (index($vcf_file,".vcf") == -1)
{
	print "\nAdded vcf suffix to $vcf_file ";

	$vcf_file = $vcf_file.".vcf";

	print "to give $vcf_file\n\n";
}

##############################################
# Make run_title out of output VCF file name #
##############################################
$run_title = &get_prefix($vcf_file);


############################
# Open the screen log file #
############################
$screen_log_file = "$run_title"."_gvcf2vcf_screen_log.out";
$| = 1;
open(STDERR, "| tee $screen_log_file");


##################################
# Define data files              #
##################################
&print_message("Which reference sequence do you want to use?","input");

print "   <1>  Dog - CanFam3\n";
print "   <2>  Horse - EquCab2\n";
print "   <3>  Human  \n";
print "   <4>  Cat - FelCat5\n\n";

print "   <5>  Strep. equi\n";
print "   <6>  Strep. zoo\n\n";

print "   <9>  Other\n\n";

$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3";}

if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/equcab2/equcab2.fasta"; $ref_seq_name = "equcab2";}
if (substr($answer,0,1) eq "3" ){$ref = "/home/genetics/human/human.fasta"; $ref_seq_name = "human";}
if (substr($answer,0,1) eq "4" ){$ref = "/home/genetics/felcat5/felCat5.fasta"; $ref_seq_name = "felcat5";}

if (substr($answer,0,1) eq "5" ){$ref = "/home/genetics/strep_equi/strep_equi.fasta"; $ref_seq_name = "s_equi";}
if (substr($answer,0,1) eq "6" ){$ref = "/home/genetics/strep_zoo/strep_zoo.fasta"; $ref_seq_name = "s_zoo";}

if ($ref eq ""){print "\n\nYou have to choose a reference sequence\n\n";exit;}

if (! -e $ref)
{
	&print_message("Reference sequence $ref_seq_name not found on this server","warning");

		print "$ref not found\n\n";
		exit;
}


######################################
# Choose your own reference sequence #
######################################
if (substr($answer,0,1) eq "9" )
{
	until (-e "$ref")
	{
		print "What is the full path of your reference sequence?\n";
		print "e.g. /home/genetics/canfam3/canfam3.fasta\n\n";
		print " >   ";
		$ref = <STDIN>;
		chomp $ref;
		if (! -e $ref){print "\n\n>>>>>>    $ref not found!.  Try again   <<<<<<\n\n";}
	}
	
	print "\n\nChoose a short name for this reference (e.g. 'dog')\n\n";
	print " >  ";
	$ref_seq_name = <STDIN>;
	chomp $ref_seq_name;
	$other_ref_sequence = "true";
}



###################################################################################
# Filtering preferences for GATK                                                  #
# If you want to change the default values of stand_call_conf and stand_emit_conf #
###################################################################################
if (($variant_caller eq "UnifiedGenotyper") || ($variant_caller eq "HaplotypeCaller") || ($variant_caller eq "HaplotypeCallerGVCF"))
{
	print "\n\n";

	&print_message("Would you like to use the default values of stand_emit_conf and stand_call_conf?","input");
	print "  (The GATK variant callers use these to control quality thresholds)\n\n";


	print "   <1>  YES   [DEFAULT]\n";
	print "   <2>  NO\n\n";

	$answer = <STDIN>;
	chomp $answer;

	if ($answer eq ""){$answer = "1"} # default

	if (substr($answer,0,1) eq "1")
	{
		$use_default_stand_values = "true";
	}

	if (substr($answer,0,1) eq "2")
	{
		$use_default_stand_values = "false";

		print "New value for stand_call_conf (default is 30 - lower produces more variants):    ";
		$stand_call_conf = <STDIN>;
		chomp $stand_call_conf;
		
		print "New value for stand_emit_conf (default is 30 - lower produces more variants):    ";
		$stand_emit_conf = <STDIN>;
		chomp $stand_emit_conf;
	} # New values

} # if variant caller is GATK



####################################################
# Get user input on all memory and threads options #
####################################################
$memory_string = "-Xmx".$memory_in_Gb."g";

 while ($choice_ok eq "false")
{
		&print_message("Current memory and CPU thread settings","message");

		print "Memory setting for java steps (Gigabytes):                     \t$memory_in_Gb\n";
		print "Number of data threads (-nt) for GATK variant calling:         \t$no_of_threads_gatk_ug_nt\n";

		print "\nWould you like to change these? (y/n)   [default = 'n'] ";
		$answer=<STDIN>;chomp $answer;$answer = lc $answer;

		if ($answer eq "y")
		{
			&print_message("Enter new values (press 'return' to keep existing value)","input");

			$choice_ok = "false";

			print "  Memory setting for java steps (in Gigabytes)                  [current value = $memory_in_Gb]:      ";
			$input = <STDIN>;chomp $input;
			if ($input ne ""){$memory_in_Gb = $input}
			if ($memory_in_Gb > 60){$memory_in_Gb = "60"}


			print "  Number of data threads (-nt) for GATK variant calling?        [current value = $no_of_threads_gatk_ug_nt]:      ";
			$input = <STDIN>;chomp $input;
			if ($input ne ""){$no_of_threads_gatk_ug_nt = $input}
			if ($no_of_threads_gatk_ug_nt > 16){$no_of_threads_gatk_ug_nt = "16"}
			

			if ($no_of_threads_gatk_ug_nt > 1)
			{
				&print_message("Number of threads warning","warning");

				print "\tIf you increase the number of threads then, although the process should run faster,\n";
				print "\tthe memory required will also increase (by the number of threads used)\n\n";

				print "\tSo if you are combining a lot of large gVCF files, you should keep -nt low\n\n";

				print "Press 'return' to continue\n\n";
				$answer=<STDIN>;
			}
		}# Yes to change memory options
		else
		{
			$choice_ok = "true";
		}

} # while choice_ok eq "false"


# Make up java memory_string
$memory_string = "-Xmx"."$memory_in_Gb"."g";

#########################
# open Command Log file #
#########################
$command_log_file = "$run_title"."_gvcf2vcf_command_log.out";
open (COMMAND_LOG, ">$command_log_file")|| die "Cannot create output file: $command_log_file";
print COMMAND_LOG "COMMAND LOG for gvcf2vcf version $version\n";


################################
# Show all the summary details #
################################
&print_message_both("SUMMARY details","message");

&print_both("Run title:                  \t$run_title\n");
&print_both("Reference sequence:         \t$ref\n\n");

&print_both("GATK region string:         \t$GATK_region_string\n\n");

&print_both("GATK constants:\n\n");

&print_both("    --stand_emit_conf:      \t$stand_emit_conf\n");
&print_both("    --stand_call_conf:      \t$stand_call_conf\n");
&print_both("    --max_alt_alleles:      \t$max_alt_alleles\n\n");

&print_both("    --validation_stringency \t$GATK_validation_stringency\n\n");

&print_both("Memory setting:             \t$memory_string\n");

&print_both("No. of data threads:        \t$no_of_threads_gatk_ug_nt\n\n");


if ($use_defined_region eq "yes")
{	
	&print_both("Chromosome:                 \t$chromosome\n");

	&print_both("GATK region string:         \t$GATK_region_string\n");
}

&print_both("\nThere are $no_of_files gVCF files\n\n");

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	&print_both("File $list_count	\t$gVCF_file_array[$list_count]\n");
}

&print_both("\nOutput VCF file:   \t$vcf_file\n\n");

print "\n\n  >> Press return to continue\n\n";
$answer=<STDIN>;
chomp $answer;
if ($answer eq "q" || $answer eq "Q"){exit;}



######################################
# First get a list of the gVCF files #
######################################
for ($list_count=1;$list_count <=$no_of_files;$list_count++)
{
	$gVCF_file = $gVCF_file_array[$list_count];
	$sample_name = &get_prefix ($gVCF_file);

	print "File: $list_count \t$gVCF_file\n";
		
	$GATK_input_string = $GATK_input_string." -I $gVCF_file";

} # list_count loop

$GATK_input_string = $GATK_input_string." ";


##########################################################
#  Make a VCF file using all the gVCF files in parallel  #
##########################################################

&print_message_both("Making VCF from gVCF files, using GATK GenotypeGVCFs","message");
$start_time = time();
printf COMMAND_LOG "  Start time:    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $start_time)[7,2,1,0];


#########################################################################################
# Loop through the gVCF files - creating an input string for GATK GenotypeGVCFs         #
#########################################################################################

for ($list_count=1;$list_count <=$no_of_files;$list_count++)
{
	$gVCF_file = $gVCF_file_array[$list_count];

	$GVCF_input_string = $GVCF_input_string." -V $gVCF_file";
} # list_count loop

$GVCF_input_string = $GVCF_input_string." ";

print "gVCF input string: $GVCF_input_string\n\n";


########################################################
# GenotypeGVCFs - run on all the files at once         #
########################################################

&print_message("Running GenotypeGVCFs on all $no_of_files files... ","message");

&print_both("Input string: $GVCF_input_string\n\n");

$output_vcf = "$run_title"."_variants_gVCF.vcf";

if ($use_default_stand_values eq "true")
{
	&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T GenotypeGVCFs -nt $no_of_threads_gatk_ug_nt  $GATK_region_string  $GVCF_input_string -o $output_vcf -S $GATK_validation_stringency","Running GenotypeGVCFs on all VCF files - default filters");	
}
if ($use_default_stand_values eq "false")
{
	&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T GenotypeGVCFs -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -nt $no_of_threads_gatk_ug_nt  $GATK_region_string  $GVCF_input_string -o $output_vcf -S $GATK_validation_stringency","Running GenotypeGVCFs on all VCF files - user-specified filters");	
}

&record_output_file_size("$output_vcf");

&print_message("The VCF files from Haplotype_Caller gVCF mode have been merged into a single VCF file $output_vcf","message");


$current_time = time();
$end_time = $current_time - $start_time;
printf COMMAND_LOG "\n\nStart time:    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $start_time)[7,2,1,0];
printf COMMAND_LOG "\nEnd time:      \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $current_time)[7,2,1,0];
printf COMMAND_LOG "\nTotal time:    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $end_time)[7,2,1,0];

close (COMMAND_LOG);
close (STDERR);

print "\n\n";

######################################################
# Copy command log to command log folder in genetics #
######################################################
&run_unix_command("cp $command_log_file /home/genetics/command_logs/$command_log_file","Copy command log to genetics folder");


####################
# Finished program #
####################
&print_message("FINISHED RUNNING gvcf2vcf","message");

&print_both("\nThere were $no_of_files gVCF files\n\n");

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
&print_both("File $list_count	\t$gVCF_file_array[$list_count]\n");
}

&print_both("\nOutput VCF file:   \t$vcf_file\n\n");

print "For information on this run look at the file with a list of all commands:  \t$command_log_file\n\n";
print "Also look at the screen log file which captures some of the screen output: \t$screen_log_file\n\n";
	
exit;


#############################################
#                                           #
# Subroutine to execute unix command        #
#                                           #
#############################################
sub run_unix_command
{
	my $unix_command = "";
	$unix_command = $_[0];
		
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	print COMMAND_LOG "$unix_command\n";

	print("$unix_command\n\n");
	system("$unix_command");
}


####################################################################
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test" from "test.fasta")                           #
####################################################################
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

##################################################################
# Subroutine to print screen message and add to command log file #
##################################################################

sub print_message_both
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
	print COMMAND_LOG "\n\n";

	print color ' bold yellow';
	if ($style eq "warning"){print color ' bold red'}
	if ($style eq "input"){print color ' bold white'}
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char; print COMMAND_LOG $char}
	
	print "\n$char    $message    $char\n";
	print COMMAND_LOG "\n$char    $message    $char\n";

	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char; print COMMAND_LOG $char}
	
	print "\n\n";
	print COMMAND_LOG "\n\n";

	print color 'reset';

}

####################################################
# Subroutine to print to screen and to COMMAND_LOG #
####################################################

sub print_both
{
	my $message = $_[0];

	print "$message";
	print COMMAND_LOG "$message";
}



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
	
	if ($testing_mode eq "on")
	{
		print "\nPRESS 'RETURN' TO CONTINUE\n";
		$answer = <STDIN>;
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

##########################################################
# Subroutine to record size of INPUT and OUTPUT files    #
##########################################################
sub record_input_and_output_file_size
{
	my $inputfile = "";
	my $outputfile = "";
	my $filesize_in = "";
	my $filesize_out = "";	

	$inputfile = $_[0];
	$outputfile = $_[1];
	
	if (-e $inputfile)
	{
		$filesize_in = -s "$inputfile";
		print COMMAND_LOG "\n  Input file: $inputfile\t\tSize: $filesize_in\n";
		print "\n  Input file: $inputfile\t\tSize: $filesize_in\n"
	}
	else
	{
		$filesize_in = 0;
		print COMMAND_LOG "\n  Input file: $inputfile\t\tSize: Not found\n";
		print "\n  Input file: $inputfile\t\tSize: Not found\n"
	}
	if (-e $outputfile)
	{
		$filesize_out = -s "$outputfile";
		print COMMAND_LOG "  Output file: $outputfile\t\tSize: $filesize_out\n\n";
		print "  Output file: $outputfile\t\tSize: $filesize_out\n\n"
	}
	else
	{
		$filesize_out = 0;
		print COMMAND_LOG "  Output file: $outputfile\t\tSize: Not found\n\n";
		print "  Output file: $outputfile\t\tSize: Not found\n\n"
	}
}

sub get_chr_start_and_end
{
	$chromosome_end_array[1]=122678785;
	$chromosome_end_array[10]=69331447;
	$chromosome_end_array[11]=74389097;
	$chromosome_end_array[12]=72498081;
	$chromosome_end_array[13]=63241923;
	$chromosome_end_array[14]=60966679;
	$chromosome_end_array[15]=64190966;
	$chromosome_end_array[16]=59632846;
	$chromosome_end_array[17]=64289059;
	$chromosome_end_array[18]=55844845;
	$chromosome_end_array[19]=53741614;
	$chromosome_end_array[2]=85426708;
	$chromosome_end_array[20]=58134056;
	$chromosome_end_array[21]=50858623;
	$chromosome_end_array[22]=61439934;
	$chromosome_end_array[23]=52294480;
	$chromosome_end_array[24]=47698779;
	$chromosome_end_array[25]=51628933;
	$chromosome_end_array[26]=38964690;
	$chromosome_end_array[27]=45876710;
	$chromosome_end_array[28]=41182112;
	$chromosome_end_array[29]=41845238;
	$chromosome_end_array[3]=91889043;
	$chromosome_end_array[30]=40214260;
	$chromosome_end_array[31]=39895921;
	$chromosome_end_array[32]=38810281;
	$chromosome_end_array[33]=31377067;
	$chromosome_end_array[34]=42124431;
	$chromosome_end_array[35]=26524999;
	$chromosome_end_array[36]=30810995;
	$chromosome_end_array[37]=30902991;
	$chromosome_end_array[38]=23914537;
	$chromosome_end_array[4]=88276631;
	$chromosome_end_array[5]=88915250;
	$chromosome_end_array[6]=77573801;
	$chromosome_end_array[7]=80974532;
	$chromosome_end_array[8]=74330416;
	$chromosome_end_array[9]=61074082;
	$chromosome_end_array[39]=123869142;
}
