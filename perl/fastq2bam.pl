#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	FASTQ2BAM           						                        #     
#									                                    #
#	PROCESS FASTQ FILES TO BAM FILES           	                        #
#									                                    #
#########################################################################

#############################
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# oliver.forman@aht.org.uk  #
# mike.boursnell@aht.org.uk #
#############################

#####################################################
# From an original script by Oliver Forman Sep 2010 #
# May 2012 - Modifed my Mike Boursnell              #
#####################################################

use strict;
use Getopt::Std ;
use File::Basename ;
#use warnings::unused;

# VERSION OF SOFTWARE #
my $version						= "20";


####################
# Define variables #
####################
#paths for software
my $bwa_path					= "/opt/bwa/bwa"; # this runs the current default version of bwa (NEW VERSION)

#Constants
my $testing_mode				= "off";
my $picard_validation_stringency = "LENIENT";
my $delete_intermediate_files 	= "yes";
my $max_open_temp_files			= 7900;

my $list_count					= 0; #Counter for the list of input files
my $no_of_files					= 0; #No of files in the list of input files (or no of lines if Paired Ends)
my $loop_count					= 0; #Loops round once for each of the files in the list of files
my $run_time					= 0;
my $next_folder					= 0;
my $start_folder				= 1;
my $array_size					= 0; # Size of item array to check there are two columns in the file
my $keep_all_for_first_file		= "yes"; # This keeps all intermediate files for first file, for error-checking

my $results_fastq2bam_folder			= "";  # Folder for all results (no sub folders any more)
my $use_other_ref_sequence		= "no";
my $dummy_dbsnp_file			= "";
my $species						= "";
my $command						= "";
my $ref							= "";
my $ref_prefix					= "";
my $reads						= "";
my $reads2						= "";
my $mem							= "";
my $proceed						= "";
my $data						= "";
my $email_address				= "";
my $run_title					= "";
my $read_file_method			= "";
my $answer						= "";
my $single_line					= "";
my $ref_seq_name				= "";  # Name of referecne sequence for pindel analysis
my $log_file					= "";
my $sample_name					= "";  # Name of each sample. Used to give the .vcf file an individual sample name
my $new_sample_name				= ""; #sample name given in input file (optional)
my $start_time					= "";
my $end_time					= "";
my $current_time				= "";
my $last_current_time			= "";
my $stage_time					= ""; # Time for each stage of the pipeline
my $third_column_found			= "false";
my $second_column_found			= "false";
my $lib							= ""; #This is used by the ReadGroups section and i have made it let us record the reference sequence
my $fix_quals_string			= ""; # Used if it is necessary to fix old-style Illumina (pre 1.8) quality scores
my $bwa_alignment_method		= ""; # this can be aln or mem (see bwa documentation)

#### File names ####

my $list_file					= ""; # input file of file names of fastq files
my $aln_sa_sai					= "";
my $aln_sa1_sai					= "";
my $aln_sa2_sai					= "";
my $aligned_sam					= "";
my $aligned_sorted_bam  		= "";
my $aligned_sorted_bai			= "";
my $aligned_sorted_rg_bam		= "";
my $aligned_sorted_rg_bai		= "";
my $validate_out				= "";
my $aligned_sorted_viewed_bam	= "";
my $aligned_sorted_viewed_bai	= "";
my $final_bam 					= "";
my $final_bai					= "";
my $command_log					= "";
my $title						= "";
my $flagstat_out				= "";
my $readme_file					= "";
my $bam_input_file				= ""; # File of file names for bam2vcf

my @reads1_file_array			= ();
my @reads2_file_array			= ();
my @sample_name_array			= ();
my @item						= ();


########################
# Define non variables #
########################

$title='Perl Mail demo';
my $from= 'NGS_analysis@samba64.aht.org.uk';
			
			
################
# START TIMERS #
################
BEGIN { our $start_run = time(); }
BEGIN { our $start_run2 = time(); }
BEGIN { our $start_run3 = time(); }

$start_time = time();

$command = "date";
system("$command");

#################
#TURN LOGGER ON #
#################



use Term::ANSIColor;
print color 'bold cyan';


if ($testing_mode eq "on"){print"\n\nTESTING MODE ON\n\n";}

print color 'reset';

use Term::ANSIColor;
print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      FASTQ2BAM       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";

print "  - This program processes FASTQ files into 'raw' BAM files.\n\n";
print "    BAM files can then be processed with bam2vcf.\n\n";
print color 'bold magenta';
print "NOTE: this version uses the new version of bwa (bwa-07)\n\n";
print color 'reset';


#############################
# Name the analysis run     #
#############################

&print_message("Your details...","input");

print "\nPlease enter a name for this analysis.\n";
print "(Keep if fairly short as it becomes part of the file name - and with no spaces):\n\n";
$run_title = <STDIN>;
chomp $run_title;

# Name some files etc
$log_file = "$run_title"."_fastq2bam_log.out";
$results_fastq2bam_folder = "results_fastq2bam_$run_title";


$| = 1;  open(STDERR, "| tee $log_file");


#########################################################
# Check if results folder with this name exists already #
#########################################################

if (-e "$results_fastq2bam_folder")
{ 
	print "\nA results folder with the name $results_fastq2bam_folder already exists!\n\n";
	print "Please start again and choose another name, or delete the existing folder\n\n";
	exit;
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
if ($email_address eq "k"){$email_address = 'karen.steward@aht.org.uk';}
if ($email_address eq "s"){$email_address = 'sally.ricketts@aht.org.uk';}
if ($email_address eq "a"){$email_address = 'amy.charbonneau@aht.org.uk';}


#####################################
# Ask if the data is SE or PE data? #
#####################################

&print_message("Do you have a Single-end or Paired-end dataset?","input");
print "   <1> Single-end (SE)\n";
print "   <2> Paired-end (PE)\n\n";

$answer = <STDIN>; chomp $answer;

if ($answer eq ""){$answer = "2"};

if (substr($answer,0,1) eq "1" ){$data = "SE"}
if (substr($answer,0,1) eq "2" ){$data = "PE"}


##################################
# Define data files              #
##################################

&print_message("Which reference sequence do you want to use?","input");

print "   <1> CanFam3\n";
print "   <2> CanFam3nu (unknow chromosomes removed)\n";
print "   <3> CanFam2\n";
print "   <4> EquCab2\n";
print "   <5> Human\n\n";

print "   <6> Strep. equi\n";
print "   <7> Strep. zoo\n\n";

print "   <9> other\n\n";

$answer = <STDIN>; chomp $answer;

if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3/canfam3_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/canfam3nu/canfam3nu.fasta"; $ref_seq_name = "canfam3nu"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3/canfam3_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "3" ){$ref = "/home/genetics/canfam2/canfam2.fasta"; $ref_seq_name = "canfam2"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam2/canfam2_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "4" ){$ref = "/home/genetics/equcab2/equcab2.fasta"; $ref_seq_name = "equcab2"; $species = "equus_caballus";$dummy_dbsnp_file = "/home/genetics/equcab2/equcab2_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "5" ){$ref = "/home/genetics/human/human.fasta"; $ref_seq_name = "human"; $species = "homo_sapiens";$dummy_dbsnp_file = "/home/genetics/human/human_dummy_DBSNP.vcf"}

if (substr($answer,0,1) eq "6" ){$ref = "/home/genetics/strep_equi/strep_equi.fasta"; $ref_seq_name = "s_equi"; $species = "streptococcus_equi";$dummy_dbsnp_file = "/home/genetics/strep_equi/strep_equi_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "7" ){$ref = "/home/genetics/strep_zoo/strep_zoo.fasta"; $ref_seq_name = "s_zoo"; $species = "streptococcus_zoo";$dummy_dbsnp_file = "/home/genetics/strep_zoo/strep_zoo_dummy_DBSNP.vcf"}




######################################
# Choose your own reference sequence #
######################################
if (substr($answer,0,1) eq "9" )
{
	$use_other_ref_sequence = "yes";
	 
	until (-e "$ref")
	{
		&print_message("Enter the full path of your reference sequence","input");
		print "e.g. /home/genetics/canfam3/canfam3.fasta\n\n";
		print " >   ";
		$ref = <STDIN>;
		chomp $ref;
		if (! -e $ref){print "\n\n>>>>>>    $ref not found!.  Try again   <<<<<<\n\n";}
	}
	
	$ref_prefix = &get_prefix($ref);
	
	$dummy_dbsnp_file = "$ref_prefix"."_dummy_DBSNP.vcf";
	
	
	&print_message("Choose a short name for this reference (e.g. dog)","input");
	print " >  ";
	$ref_seq_name = <STDIN>;
	chomp $ref_seq_name;
}


#########################################################
# Check if REF file exists (if a REF file is specified) #
########################################################

if (! -e "$ref")
{ 
	&print_message("File $ref does not exist","warning");
	print "  You might need to make a new_bwa version of the file.\n\n";
	exit;
}

$lib = $ref_seq_name;

###################################################################
# Get bit before name of ref sequence to delete "test.dict" later #
###################################################################
$ref_prefix = &get_prefix ($ref);



###################################################################
# Check if REF.DICT files exist (if a pre-indexed file is chosen) #
###################################################################

if ((! -e "$ref") || (! -e "$ref.bwt") || (! -e "$ref.pac") || (! -e "$ref.sa") || (! -e "$ref.amb") || (! -e "$ref.dict") || (! -e "$ref.ann") || (! -e "$ref.fai"))
{ 
	print "##############################\n";
	print "#  REFERENCE SEQUENCE ERROR  #\n";
	print "##############################\n\n";
	print "\nNot all the correct files for an indexed refererence sequence exist.\n\n";
	print "You need the following files:  .dict, .bwt, .pac, .sa, .amb, .ann, .fai\n\n";
	print "This suggests that there is not already a pre-indexed reference file\n\n";
	
	if(! -e "$ref"){print "$ref does not exist\n";}
	if(! -e "$ref.bwt"){print "$ref.bwt does not exist\n";}
	if(! -e "$ref.pac"){print "$ref pac does not exist\n";}
	if(! -e "$ref.sa"){print "$ref.sa does not exist\n";}
	if(! -e "$ref.amb"){print "$ref.amb does not exist\n";}
	if(! -e "$ref.dict"){print "$ref.dict does not exist\n";}
	if(! -e "$ref.ann"){print "$ref.ann does not exist\n";}
	if(! -e "$ref.fai"){print "$ref.fai does not exist\n";}
	
	exit;
} 



##########################################################
# Ask which bwa method you want to use aln or mem        #
##########################################################

&print_message("Which bwa alignment method do you want to use?","input");

print "   <1>  mem - this is a new method avaialble with the latest version of bwa (for reads over 70bp)\n";
print "   <2>  aln - this is the method we have historically been using\n\n";

$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1" ){$bwa_alignment_method = "mem";}
if (substr($answer,0,1) eq "2" ){$bwa_alignment_method = "aln";}


##########################################################
# Ask if the fastq file has old Illumina quality scores  #
##########################################################

&print_message("Are your Quality Scores of the correct Sanger type (i.e. new style)?","input");

print "  Illumina scores prior to Illumina 1.8 were 'old-style'              \n";
print "  All current FASTQ files should be OK.                               \n\n";

print "   <1> New style\n";
print "   <2> Old style\n\n";

$answer = <STDIN>;
chomp $answer;

if ($answer eq ""){$answer = "1"} # default is New style

if (substr($answer,0,1) eq "1" ){$fix_quals_string = ""}
if (substr($answer,0,1) eq "2" ){$fix_quals_string = "-I ";}

if (($bwa_alignment_method eq "mem") && ($answer eq "2"))
{
	&print_message("The mem alignment method doesn't have an option for old-style scores","message");
	print "You'll have to fix them using convert_fastq_quals or when you are running bam2vcf\n\n";
	$fix_quals_string = "";
	exit;
}


#####################################################################
# ASK IF YOU WANT TO READ THE FILENAMES FROM A "FILE OF FILE NAMES" #
#####################################################################
$answer = "";

until ($read_file_method eq "multiple" || $read_file_method eq "single")
{
	&print_message("How do you want to read the input files?","input");
	print "   <1> MULTIPLE FILES (using a file of file names)\n";
	print "   <2> SINGLE FILE\n\n";

	$answer = <STDIN>; chomp $answer;	
	if (substr($answer,0,1) eq "1" ){$read_file_method = "multiple"}
	if (substr($answer,0,1) eq "2" ){$read_file_method = "single"}

}


#####################################
# Input file names for fastq option #
#####################################

if ($read_file_method eq "single")
{
	if ($data eq "SE")
	{
		print "\nPlease input the name of your fastq sequence file:      ";
		$reads = <STDIN>;
		chomp $reads;

		
		# Add .fastq suffix if necessary
		if (index($reads,".fastq") == -1){$reads = $reads.".fastq"}

		# Assign to array of file names but just use first element #
		$reads1_file_array[1]=$reads;
		$no_of_files = 1;
	}
	
	if ($data eq "PE")
	{
		print "\nPlease input the name of your 1st fastq sequence file:      ";
		$reads = <STDIN>;
		chomp $reads;
		print "\nPlease input the name of your 2nd fastq sequence file:      ";
		$reads2 = <STDIN>;
		chomp $reads2;
		
		# Assign to array of file names but just use first element #
		$reads1_file_array[1]=$reads;
		$reads2_file_array[1]=$reads2;
		$no_of_files = 1;
	}
} #  if ($read_file_method eq "single")


if ($read_file_method eq "multiple")
{
	until (-e "$list_file")
	{
		if ($data eq "SE")
		{
			print "\nPlease input the name of your file with a list of file names of the FASTQ files:      ";
		}
		if ($data eq "PE")
		{
			print "\nPlease input the name of your file with a list of file names of the FASTQ files\n\n";
			print "(The two paired-end convert_fastq_quals file names must be on the same line separated by a TAB):      ";
		}
		
		$list_file = <STDIN>;
		chomp $list_file;

		if ($list_file ne "ls")
		{
			if (! -e $list_file){print "\n\n>>>>>>>>  File $list_file not found.  Try again.  <<<<<<<<\n\n";}
		}
		
		if ($list_file eq "ls"){print "\n";system ("ls *.txt")}
	}
	
	#############################################
	# Make sure the list file is in Unix format #
	#############################################

	$command = "dos2unix $list_file";
	#print("\n$command\n");
	system("$command");
	

	####################################################
	# Open the list file to get the list of file names #
	####################################################
	open (LIST, "$list_file") || die "Cannot open $list_file";
	$list_count=1;
	$no_of_files=0;

	while ($single_line = <LIST> ) 
	{
		chomp $single_line;

		if (length $single_line > 1)  # so if user has left blank lines at the end of the file of file names it won't matter
		{
			if ($data eq "SE")
			{
				@item=split(/\t/,$single_line);
				$array_size = scalar @item;
				
				if ($array_size == 1)
				{
					$reads = $single_line;
					if (index($reads,".fastq") == -1 ){$reads = $reads.".fastq"}
					$reads1_file_array[$list_count]=$reads;
				} # array size 1
				
				if ($array_size == 2)
				{
					$reads = $item[0];
					$new_sample_name = $item[1];
					
					if (index($reads,".fastq") == -1 ){$reads = $reads.".fastq"}
					$reads1_file_array[$list_count]=$reads;
					
					$sample_name_array[$list_count] = $new_sample_name;
					$second_column_found = "true";
					
				} # array size 2
			}
			
			if ($data eq "PE") # Read two columns separated by a TAB for the two files reads1 and reads2
			{
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
					
					$reads1_file_array[$list_count]=$reads;
					$reads2_file_array[$list_count]=$reads2;
					$sample_name_array[$list_count] = "";
				
				}
				
				#
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
					
					$reads1_file_array[$list_count] = $reads;
					$reads2_file_array[$list_count] = $reads2;
					$sample_name_array[$list_count] = $new_sample_name;
					$third_column_found = "true";
				}

			} # PE

			$list_count=$list_count + 1;

		} # if $single_line ne ""
	}

	close LIST;

	$no_of_files=$list_count - 1;
	
	
	########################
	# Check if files exist #
	########################
	
	&print_message("Checking input files...","message");

	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		if ($data eq "SE")
		{
			if (! -e "$reads1_file_array[$list_count]")
			{ 
				print "\nFile $reads1_file_array[$list_count] does not exist.\n\n";
				exit;
			} 
		}
		if ($data eq "PE")
		{
			if (! -e "$reads1_file_array[$list_count]")
			{ 
				print "\nFile $reads1_file_array[$list_count] does not exist..\n\n";
				exit;
			} 
			
			if (! -e "$reads2_file_array[$list_count]")
			{ 
				print "\nFile $reads2_file_array[$list_count] does not exist...\n\n";
				exit;
			} 
		}
	}
	

	
	###################
	# List file names #
	###################
	if ($data eq "SE"){print "\nThere are $no_of_files FASTQ files in this file of file names.\n\n";}
	if ($data eq "PE"){print "\nThere are $no_of_files pairs of FASTQ files in this file of file names.\n\n";}
	
	if ($third_column_found eq "true")
	{
		print "There is also a third column in the input file which will be used for renaming the files.\n\n";
	}
	if ($second_column_found eq "true")
	{
		print "There is also a second column in the input file which will be used for renaming the files.\n\n";
	}
	
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		if ($data eq "SE")
		{
			if ($second_column_found eq "true")
			{
				print "File $list_count	\t$reads1_file_array[$list_count]\t\t$sample_name_array[$list_count]\n";
			}
			if ($second_column_found eq "false")
			{
				print "File $list_count	\t$reads1_file_array[$list_count]\n";
			}
	
		}
		if ($data eq "PE")
		{
			if ($third_column_found eq "true")
			{
				print "Pair $list_count	\t$reads1_file_array[$list_count]  \t$reads2_file_array[$list_count]\t\t$sample_name_array[$list_count]\n";
			}
			if ($third_column_found eq "false")
			{
				print "Pair $list_count	\t$reads1_file_array[$list_count]  \t$reads2_file_array[$list_count]\n";
			}
		}
	}

	print "\nIf these file names look OK, press enter to proceed (or 'Q' to quit):      ";
	$proceed = <STDIN>;
	chomp $proceed; 
	 
	if (lc $proceed eq "q"){exit;} 

} # End of if ($read_file_method eq "multiple")

	

&print_message("Please enter the memory setting for this analysis","input");

print "Memory setting (default is -Xmx4g):      ";
$mem = <STDIN>;
chomp $mem;
 
if ($mem eq ""){$mem = "-Xmx4g"}
	


#################################
# Show details to user to check #
#################################
&print_message("SUMMARY - PLEASE CHECK THESE DETAILS CAREFULLY!!","message");

&print_message("YOUR DETAILS","message");

print "  Name of this analysis:         \t$run_title\n\n"; 
print "  Your email address:            \t$email_address\n\n";

&print_message("SETTINGS","message");

if ($data eq "SE")        {print "Single-end analysis            \tYES\n\n";}
if ($data eq "PE")        {print "Paired-end analysis            \tYES\n\n";}

print "  Memory setting:                \t$mem\n\n";
print "  BWA alignment option:          \t$bwa_alignment_method\n\n";


&print_message("DATA FILES","message");

print "  Reference sequence:            \t$ref\n\n";

print "  Number of FASTQ files to analyse:\t$no_of_files\n\n";
	
for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	if ($data eq "SE"){print "   Reads file $list_count	\t$reads1_file_array[$list_count]\n";}
	if ($data eq "PE"){print "   Pair of reads files $list_count	\t$reads1_file_array[$list_count]\t$reads2_file_array[$list_count]\n";}
}


&print_message("Please press 'ENTER' to proceed with the analysis run (or 'Q' to quit)","input");

$proceed = <STDIN>;
chomp $proceed; 
if (lc $proceed eq "q"){exit;}
	

###########################
# Make up some file names #
###########################	 
$bam_input_file = "$run_title"."_bam_input.txt";
$command_log = "$run_title"."_fastq2bam_command_log.out";


#######################################
# Open file of file names for bam2vcf #
#######################################
open (BAM_FOF, ">$bam_input_file")|| die "Cannot create output file: $bam_input_file";



#########################
# open Command Log file #
#########################

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "COMMAND LOG for fastq2bam version$version\n\n";
print COMMAND_LOG "Analysis name:        \t$run_title\n\n";
print COMMAND_LOG "Reference sequence:   \t$ref_seq_name\n\n";

print COMMAND_LOG "bwa path used:        \t$bwa_path\n";
print COMMAND_LOG "bwa alignment method: \t$bwa_alignment_method\n\n";

print COMMAND_LOG "List of input files\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	if ($data eq "SE")
	{
		print COMMAND_LOG "FASTQ FILE $list_count	\t$reads1_file_array[$list_count]\n";
	}
	if ($data eq "PE")
	{
		print COMMAND_LOG "FASTQ FILES $list_count	\t$reads1_file_array[$list_count]\t$reads2_file_array[$list_count]\n";
	}
}

print COMMAND_LOG "\n";

##############################################################################
##############################################################################
# Start MAIN LOOP.  This runs through each of the files in the list of files #
##############################################################################
##############################################################################

$last_current_time = time();

for ($loop_count=1;$loop_count <=$no_of_files;$loop_count++)
{
	######################################################
	# Assign some file names which change for each loop  #
	######################################################
	
	$reads = $reads1_file_array[$loop_count];
	$reads2 = $reads2_file_array[$loop_count];
	
	
	######################################################
	# Create sample name (to use for name BAM files etc) #
	######################################################
	
	$sample_name = &get_prefix ($reads);
	

	#############################################################################
	# If a sample name has been specified in the third column of the input file #	
	#############################################################################
	if ($sample_name_array[$loop_count] ne "")
	{
		$sample_name = $sample_name_array[$loop_count];
	}
	if (substr($sample_name,length($sample_name)-2,2) eq "_1")
	{
		$sample_name = substr($sample_name,0,length($sample_name)-2);
	}
	
	&print_message("Loop count: $loop_count   Sample name: $sample_name","message");
	
	
	#########################################################
	# Set up various File names which vary for each loop    #
	#########################################################

	$aln_sa_sai = "$run_title"."_"."$sample_name"."_aln_sa.sai";
	$aln_sa1_sai = "$run_title"."_"."$sample_name"."_aln_sa1.sai";
	$aln_sa2_sai = "$run_title"."_"."$sample_name"."_aln_sa2.sai";
	$aligned_sam = "$run_title"."_"."$sample_name"."_aligned.sam";
	$aligned_sorted_bam = "$run_title"."_"."$sample_name"."_aligned_sorted.bam";
	$aligned_sorted_bai = "$run_title"."_"."$sample_name"."_aligned_sorted.bai";
	$aligned_sorted_rg_bam = "$run_title"."_"."$sample_name"."_aligned_sorted_rg.bam";
	$aligned_sorted_rg_bai = "$run_title"."_"."$sample_name"."_aligned_sorted_rg.bai";
	$aligned_sorted_viewed_bam = "$run_title"."_"."$sample_name"."_aligned_sorted_viewed.bam";
	$aligned_sorted_viewed_bai = "$run_title"."_"."$sample_name"."_aligned_sorted_viewed.bai";
	$final_bam = "$run_title"."_"."$sample_name"."_f2b.bam";
	$final_bai = "$run_title"."_"."$sample_name"."_f2b.bai";
	$validate_out = "$run_title"."_"."$sample_name"."_validate_summary.out";
	$flagstat_out = "$run_title"."_"."$sample_name"."_flagstat.out";
	$readme_file = "$run_title"."_"."$sample_name"."_readme.out";
	
	
	####################################################################
	# 1 Run bwa                                                        #
	####################################################################

	#####################
	# SINGLE END DATA ! #
	#####################

	if  ($data eq "SE")
	{		

		if ($bwa_alignment_method eq "aln")
		{
			print"\n====================================================================================";
			print "\nFile $loop_count/$no_of_files   Step 1/5 Step 1 Running BWA to produce SAI file - SE";
			print "\n===================================================================================\n\n\n";
			
			#Step 1 Running BWA to produce aligned BAM file - SE

			&print_message("File $loop_count/$no_of_files   Step 1/5 Running BWA to produce aligned BAM file - SE","message");

			&run_unix_command("$bwa_path $bwa_alignment_method $fix_quals_string $ref $reads > $aln_sa_sai","1");

			&record_input_file_size("$reads");
			&record_output_file_size("$aln_sa_sai");


			print"\n====================================================================================";
			print "\nFile $loop_count/$no_of_files   Step 1/5 Step 1 COMPLETED Running BWA to produce SAI file - SE";
			print "\n====================================================================================\n\n\n";

			&test_mode_subroutine;
		} # bwa aln

		if ($bwa_alignment_method eq "mem")
		{
			print"\n====================================================================================";
			print "\nFile $loop_count/$no_of_files   Step 1/5 Looking for suffix array coordinates for good sequence hits";
			print "\n===================================================================================\n\n\n";

			&print_message("File $loop_count/$no_of_files  Step 1/5 Running BWA mem to produce aligned SAM file - SE","message");

			&run_unix_command("$bwa_path mem -M $ref $reads > $aligned_sam","Converting FASTQ to SAM using bwa mem");

			&record_input_file_size("$reads");
			&record_output_file_size("$aligned_sam");


			print"\n====================================================================================";
			print "\nFile $loop_count/$no_of_files   Step 1/5 COMPLETED Running BWA mem to produce aligned SAM file";
			print "\n====================================================================================\n\n\n";

			&test_mode_subroutine;
		} # bwa mem

	} # SE

	#####################
	# PAIRED END DATA ! #
	#####################

	if  ($data eq "PE")
	{
	
		if ($bwa_alignment_method eq "aln")
		{

			##########
			# Read 1 #
			##########
			
			&print_message("File $loop_count/$no_of_files   Step 1a/5 Running BWA ALN to produce SAI file - PE1","message");

			&run_unix_command("$bwa_path $bwa_alignment_method $fix_quals_string $ref $reads > $aln_sa1_sai","1a");

			&record_input_file_size("$reads");
			&record_output_file_size("$aln_sa1_sai");
			
			&print_message("File $loop_count/$no_of_files   Step 1a/5 COMPLETED Running BWA ALN to produce SAI file - PE1","message");

			&test_mode_subroutine;
			
			
			##########
			# Read 2 #
			##########
			
			&print_message("File $loop_count/$no_of_files   Step 1b/5 Running BWA ALN to produce SAI file - PE2","message");
			
			&run_unix_command("$bwa_path aln $fix_quals_string $ref $reads2 > $aln_sa2_sai","1b");
			
			&record_input_file_size("$reads2");
			&record_output_file_size("$aln_sa2_sai");

			&print_message("File $loop_count/$no_of_files   Step 1b/5 COMPLETED Running BWA ALN to produce SAI file - PE2","message");

			&test_mode_subroutine;

		} #bwa aln

		if ($bwa_alignment_method eq "mem")
		{

			################################################
			# Read 1 and Read 2 are done together with mem #
			################################################
			
			&print_message("File $loop_count/$no_of_files   Step 1/5 Running BWA MEM to produce SAM file - PE","message");
			
			&run_unix_command("$bwa_path mem -M $ref $reads $reads2 > $aligned_sam","Converting FASTQ to SAM using bwa mem");

			&record_input_file_size("$reads");
			&record_input_file_size("$reads2");
			&record_output_file_size("$aligned_sam");
			
			&print_message("File $loop_count/$no_of_files   Step 1/5 COMPLETED Running BWA MEM to produce SAM file - PE","message");

			&test_mode_subroutine;
			
		} #bwa mem

	} # PE

	&problem_check;


	##########################################################################
	# 2 Convert SAI files into SAM files  (aln only, not required for mem)   #
	##########################################################################
	if ($bwa_alignment_method eq "aln")
	{
		if  ($data eq "SE")
		{
			print "\n\n===================================================================================";
			print "\nFile $loop_count/$no_of_files   Step 2/5 Converting suffix array coordinates into chromosomal coordinates";
			print "\n=====================================================================================\n";

			&run_unix_command("$bwa_path samse $ref $aln_sa_sai $reads > $aligned_sam","Converting SAI to SAM");

			&record_input_file_size("$aln_sa_sai");
			&record_output_file_size("$aligned_sam");

			&test_mode_subroutine;
		
			&problem_check;
		
			print "\n==============================================================================================";
			print "\nFile $loop_count/$no_of_files   Step 2/5 Suffix array coordinates have been converted into chromosomal coordinates";
			print "\n==============================================================================================\n\n\n";
		}

		if  ($data eq "PE")
		{
			print "\n\n==============================================================================================";
			print "\nFile $loop_count/$no_of_files   Step 2/5 Converting suffix array coordinates into chromosomal coordinates";
			print "\n==============================================================================================\n\n";

			&run_unix_command("$bwa_path sampe $ref $aln_sa1_sai $aln_sa2_sai $reads $reads2 > $aligned_sam","Converting SAI to SAM");

			&record_input_file_size("$aln_sa1_sai");
			&record_input_file_size("$aln_sa2_sai");
			&record_output_file_size("$aligned_sam");
			
			&test_mode_subroutine;
		
			&problem_check;
		
			print "\n==============================================================================================";
			print "\nFile $loop_count/$no_of_files   Step 2/5 Suffix array coordinates have been converted into chromosomal coordinates";
			print "\n==============================================================================================\n\n\n";
		}

	} # aln only

	print " * * * \n";



	##########################################################
	# Run all the picard stages in one                       #
	#  - Sort the SAM file generated by BWA                  #
	#  - Convert the SAM file (from BWA) to a BAM file       #
	#  - Add Read Group information to the BAM file          #	
	##########################################################


	&print_message("File $loop_count/$no_of_files   Step 3/5 Picard stages","message");
	print "  - ReadGroup info to be added to the BAM file\n";
	print "  - The aligned reads to be sorted\n";
	print "  - The file to be converted to BAM format\n\n";


	&run_unix_command("java $mem -jar /opt/picard/AddOrReplaceReadGroups.jar I=$aligned_sam O=$aligned_sorted_rg_bam rgID=$run_title LB=$lib PL='ILLUMINA' PU=$sample_name SM=$sample_name SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=$picard_validation_stringency","3");	
	
	# Previous individual commands
	#&run_unix_command("java $mem -jar /opt/picard/AddOrReplaceReadGroups.jar I=$aligned_sorted_viewed_bam O=$final_bam rgID=$run_title LB=$lib PL='ILLUMINA' PU=$sample_name SM=$sample_name","6");	
	
	#&run_unix_command("java $mem -jar /opt/picard/SortSam.jar I=$aligned_sam O=$aligned_sorted_sam SO=coordinate VALIDATION_STRINGENCY=$picard_validation_stringency","3");

	#&run_unix_command("java $mem -jar /opt/picard/SamFormatConverter.jar I=$aligned_sorted_sam O=$aligned_sorted_bam VALIDATION_STRINGENCY=$picard_validation_stringency","4");

	&record_input_file_size("$aligned_sam");
	&record_output_file_size("$aligned_sorted_rg_bam");

	&problem_check;
	
	&print_message("File $loop_count/$no_of_files   Step 3/5 Picard stages COMPLETED","message");

	

	################################################
	# 4 Use picard ValidateSamFile Summary         #
	################################################
	
	&print_message("File $loop_count/$no_of_files   Step 4/5 Picard ValidateSamFile","message");
	
	&run_unix_command("java $mem -jar /opt/picard/ValidateSamFile.jar I=$aligned_sorted_rg_bam O=$validate_out MODE=SUMMARY MAX_OPEN_TEMP_FILES=$max_open_temp_files VALIDATION_STRINGENCY=SILENT","4");

	&record_input_file_size("$aligned_sorted_rg_bam");
	&record_output_file_size("$validate_out");

	&print_message("File $loop_count/$no_of_files   Step 4/5 Picard ValidateSamFile COMPLETED","message");
	
	

	################################################
	# 5 Use samtools flagstat                     #
	################################################
	
	&print_message("File $loop_count/$no_of_files   Step 5/5 samtools flagstat","message");
	
	&run_unix_command("/opt/samtools/samtools flagstat $aligned_sorted_rg_bam > $flagstat_out","5");

	&record_input_file_size("$aligned_sorted_rg_bam");
	&record_output_file_size("$flagstat_out");

	&print_message("File $loop_count/$no_of_files   Step 5/5 samtools flagstat COMPLETED","message");
	
	
	
	######################################
	# Create the results folder          #
	######################################

	&run_unix_command("mkdir $results_fastq2bam_folder","Make results folder");


	##########################################
	# Move various files to a Results Folder #
	##########################################
	
	&move_to_results_folder ("$aligned_sorted_rg_bam");
	&move_to_results_folder ("$aligned_sorted_rg_bai");
	
	&move_to_results_folder ("$validate_out");
	&move_to_results_folder ("$flagstat_out");
	&move_to_results_folder ("$readme_file");


	###########################################
	# Rename as final.bam (in results folder) #
	###########################################

	&run_unix_command("mv $results_fastq2bam_folder/$aligned_sorted_rg_bam $results_fastq2bam_folder/$final_bam","Rename BAM file");
	&run_unix_command("mv $results_fastq2bam_folder/$aligned_sorted_rg_bai $results_fastq2bam_folder/$final_bai","Rename BAI file");
	
	&record_output_file_size("$results_fastq2bam_folder/$final_bam");
	&record_output_file_size("$results_fastq2bam_folder/$final_bai");
	
	print "\n\n";
	

	############################################################################
	# Add name of BAM file for bam file to file names for next stage (bam2vcf) #
	############################################################################
	print BAM_FOF ("$final_bam\t$sample_name\n");


	##########
	# README #
	##########

	open (READMEFILE, ">$readme_file"); 

	print READMEFILE "##################################\n";
	print READMEFILE "Summary of fastq2bam results files\n";
	print READMEFILE "##################################\n\n";

	print READMEFILE "Program:\t\t\tfastq2bam\tVersion: $version\n\n";
		
	print READMEFILE "Run title:\t\t\t$run_title\n\n";
	
	print READMEFILE "Reference sequence:\t$ref_seq_name\n\n";
	
	print READMEFILE "\nBAM ALIGNMENT FILE (Raw alignments to the reference)\n\n";

	print READMEFILE "\t$final_bam\n\n\n";
	
	print READMEFILE "To check this BAM file look at these two files:\n\n";
	
	print READMEFILE "\t$validate_out\n\n";
	print READMEFILE "\t$flagstat_out\n\n\n";
	
	print READMEFILE "This BAM file now needs to go through the rest of the NGS pipeline using bam2vcf\n\n";

	close (READMEFILE);

	&move_to_results_folder ("$readme_file");


	##############################################
	# Delete these files at the end of each loop #
	##############################################
	if (-e "$aln_sa_sai"){&delete_file ("$aln_sa_sai")}
	if (-e "$aln_sa1_sai"){&delete_file ("$aln_sa1_sai")}
	if (-e "$aln_sa2_sai"){&delete_file ("$aln_sa2_sai")}

	if ($delete_intermediate_files eq "yes") 
	{
		&delete_file ("$aligned_sam");
	}

	########################################
	# Calculate the run time for this loop #
	########################################
	$current_time = time();
	$run_time = $current_time - $start_time;


	##########################################
	# Send an e-mail at the end of each loop #
	##########################################
	open(MAIL, "|/usr/sbin/sendmail -t");

	## Mail Header
	print MAIL "To: $email_address\n";
	print MAIL "From: $from\n";
	print MAIL "Subject: FASTQ2BAM: Sample $sample_name in run $run_title has finished\n\n";
	## Mail Body
	print MAIL "fastq2bam PERL script version $version\n\n";
	print MAIL "The processing of FASTQ files into BAM files, of $sample_name in run $run_title is complete\n\n";
	print MAIL "Run time so far : $run_time seconds\n";
	printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
	close(MAIL);

	print COMMAND_LOG "\n\n";
	print COMMAND_LOG "========================================================================\n";
	print COMMAND_LOG "                       End of loop $loop_count                          \n";
	print COMMAND_LOG "========================================================================\n\n\n";
} # End of MAIN LOOP

close BAM_FOF; # bam file of file names



#####################################################
# Delete these files after all loops have completed #
#                                                   #
# (But not if you are using pre-indexed reference   #
# files such as CanFam2)                            #
#####################################################

$current_time = time();
$run_time = $current_time - $start_time;

open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $email_address\n";
print MAIL "From: $from\n";
print MAIL "Subject: FASTQ2BAM: Run $run_title has finished\n\n";
## Mail Body
print MAIL "fastq2bam script version $version\n\n";
print MAIL "Your next generation sequence analysis ($run_title) is complete\n\n";
print MAIL "For pipeline details see $log_file\n\n";
print MAIL "For individual file information see readme files\n\n";
print MAIL "For a list of commands see $command_log and $log_file\n\n";
print MAIL "Also check _validate_out and _flagstat_out for checks on the BAM files\n\n";

printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
close(MAIL);

print "\n\n";
print "################################################################\n";
print "# FASTQ2BAM ANALYSIS COMPLETE! YOU HAVE BEEN NOTIFIED BY EMAIL #\n";
print "################################################################\n\n";

print "Run title: $run_title\n\n";
print "Run time: ";
printf "%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];

print "\n\nFinal results folder:\t$results_fastq2bam_folder";

print "\n\nFor details of the analysis run:\tcheck $command_log and $log_file for full details\n\n";

print "For details of the BAM file produced:\tcheck $validate_out and $flagstat_out for checks on the BAM file\n\n";


#############################
# Turn off logging          #
#############################

close(STDERR);
close (COMMAND_LOG);

&move_to_results_folder ("$log_file");
&move_to_results_folder ("$command_log");
&move_to_results_folder("$bam_input_file"); # FOR for bam2vcf

exit;
#############################################
# Subroutine to move file to results folder #
#############################################
sub move_to_results_folder
{
	my $file_to_be_moved = "";	
	$file_to_be_moved = $_[0];
	$command = "mv  $file_to_be_moved $results_fastq2bam_folder/$file_to_be_moved";
	print("$command\n");
	print COMMAND_LOG ("$command\n");
	system("$command");
}


#############################################
# Subroutine to delete files                #
#############################################

sub delete_file
{
	my $file_to_be_deleted = "";	

	$file_to_be_deleted = $_[0];
	
	if (! -e "$file_to_be_deleted")
	{
		print "\n$file_to_be_deleted could not be found to be deleted\n";
		print COMMAND_LOG "\n$file_to_be_deleted could not be found to be deleted\n";
	}
	
	if (-e "$file_to_be_deleted")
	{
		$command = "rm  $file_to_be_deleted";
		system("$command");
		print COMMAND_LOG "\n$file_to_be_deleted was deleted\n";
	}
}

#############################################
# Subroutine to execute unix command        #
#############################################
sub run_unix_command
{
	my $unix_command = "";
	my $step = "";	
	$unix_command = $_[0];
	$step = $_[1];
	print "\n";
	print("$unix_command\n");
	system("$unix_command");
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	print COMMAND_LOG "File: $loop_count/$no_of_files \tStep: $step/5\n";
	print COMMAND_LOG "$unix_command\n";

}

##############################################
# Subroutine to record size of output file   #
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
		print "\n  Output file: $outputfile\t\tSize: $filesize\n\n";
	}
	else
	{
		$filesize=0;
		print COMMAND_LOG "\n  Output file: $outputfile\t\tSize: Not found\n";
		print "\n  Output file: $outputfile\t\tSize: Not found\n\n";
	}


	$current_time = time();
	$end_time = $current_time - $start_time;
	$stage_time	= $current_time - $last_current_time;
	printf COMMAND_LOG "\n  Total time:    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $end_time)[7,2,1,0];
	printf COMMAND_LOG "  Time for stage:\t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $stage_time)[7,2,1,0];
	
	$last_current_time = $current_time;
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
# Subroutine to move log     			    #
#############################################
sub move_log
{

	my $file_to_be_moved = "";	

	$file_to_be_moved = $_[0];
	$command = "mv  $file_to_be_moved $results_fastq2bam_folder/$file_to_be_moved";
	system("$command");

}

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



###########################
# Carries out error check #
# and e-mails user        #
###########################
sub problem_check
{
	#print "\nERR checking switched OFF\n\n";
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

sub test_mode_subroutine
{
		if ($testing_mode eq "on")
		{
			print "\n\n>>>>>>>>>>>>>  Testing mode is set to ON              <<<<<<<<<<<<<<<<<\n";
			print ">>>>>>>>>>>>>  Press return to see a list of files    <<<<<<<<<<<<<<<<<\n\n";
			$answer = <STDIN>;
			$command = "ls -lh $sample_name*";
			system("$command");
			print "\n\n>>>>>>>>>>>>>  Press return to continue to the next stage    <<<<<<<<<<<<<<<<<\n\n";
			$answer = <STDIN>;
		}

}		

##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{
	foreach (@_) {s/\n//g}  
	foreach (@_) {s/\r//g}  
}