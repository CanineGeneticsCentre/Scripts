#!/usr/bin/perl -w

##################################################################################
#									                                             #      
#	NGS_pipeline  						                                         #     
#									                                             #
#	This PERL script runs the Broad Pipeline using PERL rather than Queue.jar	 #
#									                                             #
##################################################################################

#############################
# Mike Boursnell July 2012  #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# mike.boursnell@aht.org.uk #
#                           #
# based on an original PERL #
# script by Oliver Forman   #
#############################



use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;

# VERSION OF SOFTWARE #
my $version						= "11";


####################
# Define variables #
####################

#Constants
my $gatk_directory				= "gatk";
my $testing_mode				= "off";
my $ValidationStringency		= "LENIENT";
my $dummy_dbsnp_file			= "";
my $fix_mate_pairs				= "no";
my $max_open_temp_files			= "7900";
my $temp_dir_string				= " -Djava.io.tmpdir=javatempdir";
my $tempdir						= "javatempdir";
my $delete_intermediates		= "yes";
my $validate_sam_files			= "no";
my $make_smaller_bam_file		= "yes";
my $use_s_option_in_bwa_sampe	= "yes"; #This runs bwa sampe with the '-s' option which disables Smith-Waterman alignment for the unmapped mate
my $remove_duplicates			= "no";
my $use_samtools_rmdup			= "yes"; # This means use samtools rmdup rather than picard/MarkDuplicates
my $platform					= "Illumina";

my $list_count					= 0; #Counter for the list of input files
my $no_of_files					= 0; #No of files in the list of input files (or no of lines if Paired Ends)
my $loop_count					= 0; #Loops round once for each of the files in the list of files
my $read_length					= 0;
my $run_time					= 0;
my $results_folder_count		= 0; # Counts existing results folder to see if run has already partially run.
my $next_folder					= 0;
my $start_folder				= 1;

#Various Parameters (e.g. for the Unified Genotyper)
my $stand_emit_conf				= 30;
my $stand_call_conf				= 30;
my $indels_caller				= "UnifiedGenotyper"; # alternatives are UnifiedGenotyper OR SomaticIndelDetector

my $word						= "";
my $command						= "";
my $ref							= "";
my $ref_prefix					= "";
my $reads						= "";
my $mem							= "";
my $proceed						= "";
my $data_type					= "";
my $data						= "PE"; # PE of SE
my $reads2						= "";
my $email_address				= "";
my $run_title					= "";
my $platform_type 				= "";
my $reads3 						= "";
my $suffix						= "";
my $annotate_variants			= "";
my $maq							= "";
my $ref_size					= "";
my $ref_genome					= "";
my $structural_variant_analysis	= "";
my $chromosome					= "";
my $chromosome_string			= "";
my $insert						= "";
my $region						= "";
my $use_defined_region			= "";
my $input_file					= "";
my $bam_file					= "";
my $bai_file					= "";
my $bam_bai_file				= "";
my $bam_file_mates_fixed		= "";
my $filename					= "";
my $list_file					= "";
my $read_file_method			= "";
my $answer						= "";
my $single_line					= "";
my $results_folder				= "";
my $file_to_be_moved			= "";
my $input_string				= "";
my $species						= "";
my $ref_seq_name				= "";  # Name of reference sequence for pindel analysis
my $log_file					= "log.rtf";
my $exit_on_problem				= "no";
my $sample_name					= "";  # Name of each sample. Used to give the .vcf file an individual sample name
my $use_preindexed_reference	= "";  # can be 'yes' or 'no'  Flags whether you want to use pred-indexed ref seq such as CanFam2 or EquCab1
my $make_parallel_VCF			= ""; # This makes use of multiple BAM files in GATK to make a combined VCF output file

#### File names ####

my $aln_sa_sai					= "";
my $aln_sa1_sai					= "";
my $aln_sa2_sai					= "";
my $cleaned_bam					= "";
my $aligned_sam					= "";
my $aligned_sorted_sam  		= "";
my $aligned_sorted_bam  		= "";
my $aligned_sorted_RG_bam		= "";
my $aligned_sorted_bai			= "";
my $aligned_sorted_bam_bai		= "";
my $duplicates_removed_bam		= "";
my $cleaned_sorted_bam			= "";
my $cleaned_sorted_bai			= "";
my $aligned_plus_dup_bam		= "";
my $region_only_bam				= "";
my $region_only_bai				= "";
my $queryname_sorted_bam		= "";
my $recal_bam					= "";
my $recal_bai					= "";
my $reads_csv					= "";
my $snps_raw_vcf				= "";
my $snps_raw_vcf_idx			= "";
my $snps_parallel_vcf			= "";
my $snps_final_vcf				= "";
my $detailed_output_bed			= "";
my $indels_raw_vcf				= "";
my $indels_raw_vcf_alt			= "";
my $indels_raw_vcf_idx			= "";
my $indels_parallel_vcf			= "";
my $indels_final_vcf			= "";
my $indels_raw_bed				= "";
my $all_reads_bam				= ""; # Only used if you focus in on a particular region
my $all_reads_bai				= ""; # Only used if you focus in on a particular region
my $annotated_indels_vcf		= "";
my $annotated_snps_vcf			= "";
my $ref_dict					= ""; # This is the file like canfam3.dict
my $insert_size_pdf				= "";
my $gc_bias_pdf					= "";
my $command_log					= "";


# Broad file names
my $sai_coordinates_1			= "";
my $sai_coordinates_2			= "";
my $reverted_bam				= "";
my $reverted_bai				= "";
my $realigned_sam				= "";
my $realigned_bam				= "";
my $realigned_bai				= "";
my $realigned_rg_bam			= "";
my $realigned_rg_bai			= "";
my $runtitle_clean_dedup_bam 	= "";
my $runtitle_clean_dedup_bai 	= "";
my $runtitle_clean_dedup_recal_bam 	= "";
my $runtitle_clean_dedup_recal_bai 	= "";
my $bam_file_for_unified_genotyper = "";
my $runtitle_sample_pre_validation	= "";
my $runtitle_sample_post_validation = "";
my $runtitle_bam_intervals		= "";
my $runtitle_clean_bam			= "";
my $runtitle_clean_bai			= "";

my $runtitle_sample_metrics		= "";
my $runtitle_sample_pre_recal_csv = "";
my $runtitle_sample_post_recal_csv = "";
my $runtitle_sample_pre_recal_grp = "";
my $runtitle_sample_post_recal_grp = "";
my $pdf_output_dir_pre			= "";
my $pdf_output_dir_post			= "";

# Fixed preferences #
my $check_depth			= "no";

my @reads1_file_array	= ();
my @reads2_file_array	= ();
my @bam_file_array		= ();
my @item				= ();


########################
# Define non variables #
########################

my $title='Perl Mail demo';
my $from= 'NGS_analysis@samba64.aht.org.uk';
my $subject='NGS ANALYSIS';
my $reads4 = 'out.fastq';

			
			
################
# START TIMERS #
################
BEGIN { our $start_run = time(); }
BEGIN { our $start_run2 = time(); }
BEGIN { our $start_run3 = time(); }


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



#################
#TURN LOGGER ON #
#################

$| = 1;

open(STDOUT, "| tee log.rtf");


print color 'bold cyan';

$command =  "clear";
system("$command");

if ($testing_mode eq "on"){print"\n\nTESTING MODE ON\n\n";}

print color 'reset';

use Term::ANSIColor;
print color 'bold magenta';

print "\n\n";
print "                 ##################################################\n";
print color 'bold white';
print "                                   NGS_pipeline                   \n";
print " \n";
print "                      NGS data analysis pipeline for samba64\n";
print color 'reset';
print color 'bold magenta';
print "                 ##################################################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  This version is based on the Broad Pipeline, as outlined in the document:\n\n";
print "  'Best Practice Variant Detection with the GATK' \n\n";

print "  This version doesn't used the 'retired' CounCovariates\n";
print "  but uses the new 'BaseRecalibrator' instead\n\n";


print color 'reset';


#############################
# Name the analysis run     #
#############################

print "\n~~~~~~~~~~~~~~~";
print "\nYour details...";
print "\n~~~~~~~~~~~~~~~\n";

print "\nPlease enter a name for this analysis (with no spaces):    ";
$run_title = <STDIN>;
chomp $run_title;



######################################
# E-mail address for notifications   #
######################################

print "\nPlease enter your email address:    ";
$email_address = <STDIN>;
chomp $email_address;

if ($email_address eq "m"){$email_address = 'mike.boursnell@aht.org.uk';}


#####################################
# Ask how long the reads are        #
#####################################


#print "\n\n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
#print " How long are the reads from the sequencing?    \n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

#print "   Length: ";
#$read_length = <STDIN>;
#chomp $read_length;


##################################
# Define data files              #
##################################

print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " Which reference sequence do you want to use?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   Enter 1 for CanFam3\n";
print "   Enter 2 for CanFam2\n";
print "   Enter 3 for EquCab2\n";
print "   Enter 4 for Test\n\n";

$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3/canfam3_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/canfam2/canfam2.fasta"; $ref_seq_name = "canfam2"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam2/canfam2_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "3" ){$ref = "/home/genetics/equcab2/equcab2.fasta"; $ref_seq_name = "equcab2"; $species = "equus_caballus";$dummy_dbsnp_file = "/home/genetics/equcab2/equcab2_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "4" ){$ref = "/home/genetics/test/test.fasta"; $ref_seq_name = "test"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/test/test_dummy_DBSNP.vcf"}

$ref_dict="/home/genetics/$ref_seq_name/$ref_seq_name.dict";

###################################################################
# Check if REF.DICT files exist (if a pre-indexed file is chosen) #
###################################################################

if ((! -e "$ref") || (! -e "$ref.bwt") || (! -e "$ref.rbwt") || (! -e "$ref.rpac") || (! -e "$ref.pac") || (! -e "$ref.rsa") || (! -e "$ref.sa") || (! -e "$ref.amb") || (! -e "$ref_dict") || (! -e "$ref.ann") || (! -e "$ref.fai"))
{ 
	print "##############################\n";
	print "#  REFERENCE SEQUENCE ERROR  #\n";
	print "##############################\n\n";
	print "\nNot all the correct files for an indexed refererence sequence exist.\n\n";
	print "You need the following files:  .dict, .rbwt, .bwt, .rpac, .pac, .rsa, .sa, .amb, .ann, .fai, .dict\n\n";
	print "This suggests that there is not already a pre-indexed reference file\n\n";
	
	if(! -e "$ref"){print "$ref does not exist\n";}
	if(! -e "$ref.bwt"){print "$ref.bwt does not exist\n";}
	if(! -e "$ref.rbwt"){print "$ref.rbwt does not exist\n";}
	if(! -e "$ref.rpac"){print "$ref.rpac does not exist\n";}
	if(! -e "$ref.pac"){print "$ref pac does not exist\n";}
	if(! -e "$ref.rsa"){print "$ref.rsa does not exist\n";}
	if(! -e "$ref.sa"){print "$ref.sa does not exist\n";}
	if(! -e "$ref.amb"){print "$ref.amb does not exist\n";}
	if(! -e "$ref_dict"){print "$ref_dict does not exist\n";}
	if(! -e "$ref.ann"){print "$ref.ann does not exist\n";}
	if(! -e "$ref.fai"){print "$ref.fai does not exist\n";}
	
	exit;
} 



###########################################################################
# Region preferences                                                      #
# Would you like to focus the analysis on a specific region of the genome #
###########################################################################


print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Would you like to focus on alignments to specific region of the genome?\n";
print "(This is strongly advised as it speeds up the Variant Calling process)\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

print "   Enter 1 for YES\n";
print "   Enter 2 for NO\n\n";

$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1"){$use_defined_region = "yes"}
if (substr($answer,0,1) eq "2"){$use_defined_region = "no"}

if ($use_defined_region eq "yes")
{
	print "\nPlease define your region of interest (eg chr5:21000000-23000000):      ";
	$region = <STDIN>;
	chomp $region;
	
	# Use this for chromosome in GATK -L command line option
	if (index($region,":") == -1)
	{
		$chromosome = $region;
		if (index($chromosome,"chr") == -1){$chromosome = "chr"."$chromosome";}
	}
	
	if (index($region,":") > -1)
	{
		if (index($region,"chr") > -1)
		{
			$chromosome = substr($region,3,index($region,":")-3);
		}
	}
	
	$chromosome_string = " -L $region";
	
	
	
}

#########################################################
# Ask if you want to make a multiple parallel VCF file  #
#########################################################
print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Would you like to make an extra single VCF file by analysing all the BAM files in parallel?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

print "   <1>  Make a single VCF by analysing all BAM files in parallel.\n";
print "   <2>  Make multiple VCFs by processing the BAM files consecutively (as before)\n";
print "   <3>  Make parallel and consecutive VCF files to compare them\n\n";

$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1"){$make_parallel_VCF = "parallel"}
if (substr($answer,0,1) eq "2"){$make_parallel_VCF = "consecutive"}
if (substr($answer,0,1) eq "3"){$make_parallel_VCF = "parallel and consecutive"}


#####################################################################
# ASK IF YOU WANT TO READ THE FILENAMES FROM A "FILE OF FILE NAMES" #
#####################################################################
print "\n\n";

$answer = "";

until ($answer eq "1" || $answer eq "2")
{
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print " How do you want to read the input files?\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
	print "   Enter 1 for MULTIPLE FILES (using a file of file names)\n";
	print "   Enter 2 for SINGLE FILE\n\n";

	$answer = <STDIN>;
	chomp $answer;

	$answer = substr($answer,0,1);
	
	if ($answer eq "1"){$read_file_method = "multiple"}
	if ($answer eq "2"){$read_file_method = "single"}

}


###################################
# Input file names                #
###################################

if ($read_file_method eq "single")
{
	until (-e $bam_file)
	{
		print "\nPlease input the name of your aligned reads file (aligned.bam):      ";
		$bam_file = <STDIN>;
		chomp $bam_file;
	}
	
	######################################################
	# Add file type suffix .bam if user hasn't added it  #
	######################################################

	if (index($bam_file,".bam") == -1 ){$bam_file = $bam_file.".bam"}

	# Assign to array of file names but just use first element #
	$bam_file_array[1]=$bam_file;
	$no_of_files=1;
	
} # End of read = single

if ($read_file_method eq "multiple")
{
	until (-e "$list_file")
	{
		print "\nPlease input the name of your file with a list of file names of the .bam files:      ";
		$list_file = <STDIN>;
		chomp $list_file;

		if ($list_file eq "ls")
		{
			print "\n";
			system ("ls *.txt");
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
	
	####################################################
	# Open the list file to get the list of file names #
	####################################################
	open (LIST, "$list_file") || die "Cannot open $list_file";
	$list_count=1;
	
	while ($bam_file = <LIST> ) 
	{
		chomp $bam_file;
		
		######################################################
		# Add file type suffix .bam if user hasn't added it  #
		######################################################

		if (index($bam_file,".bam") == -1 ){$bam_file = $bam_file.".bam"}
	
		$bam_file_array[$list_count]=$bam_file;
		$list_count=$list_count + 1;
	}

	close LIST;

	$no_of_files=$list_count - 1;
	
	
	###################
	# List file names #
	###################
	print "\n\nThere are $no_of_files BAM files in this file of file names.\n\n";
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		print "File $list_count	\t$bam_file_array[$list_count]\n";
	}

	print "\nIf these file names look OK, press enter to proceed (or 'Q' to quit):      ";
	$proceed = <STDIN>;
	chomp $proceed; 
	 
	if (lc $proceed eq "q"){exit;} 
	
	
	
	#########################################################
	# Check if results folder with this name exists already #
	#########################################################

	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		$results_folder = "results_"."$run_title"."_"."$list_count";
		
		if (-e $results_folder)
		{
			$results_folder_count = $results_folder_count + 1;
		}
		print "File $list_count	\t$bam_file_array[$list_count]\n";
	}
	
	# Record penultimate folder (most likely to be complete)
	$next_folder = $results_folder_count + 1;
	
	if ($results_folder_count > 1)
	{ 
		print "\n\nResults folder(s) for this analysis name ($run_title) already exist:\n\n";
		
		for ($list_count=1;$list_count<=$no_of_files;$list_count++)
		{
			$results_folder = "results_"."$run_title"."_"."$list_count";
			
			if (-e $results_folder)
			{
				print "$list_count:  Folder $results_folder already exists\n";
			}
		}
		if ($results_folder_count == $no_of_files)
		{
			print "\nChoose another analysis name (or else all these previous folders will be OVERWRITTEN!):  ";
			$run_title=<STDIN>;
			chomp $run_title;
		}
		if ($results_folder_count < $no_of_files)
		{
			print "\nDo you want to start this analysis at a file other than number 1?\n";
			print "(For example if the run has been interrupted)\n\n";
			
			print "If so, check the contents of the results folders (the last one may not have all the required files)\n";
			print "and then choose to start at $results_folder_count or $next_folder";
			
			print "\n\nNumber to start at: ";
			$start_folder = <STDIN>;
			chomp $start_folder;
			
			if ($start_folder eq ""){$start_folder = $next_folder;}
			
			print "\n\nStart folder has been set at $start_folder\n\n";
		}
	}



} # End of if ($read_file_method eq "multiple")




print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Please enter the memory setting for this analysis\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

print "Memory setting (default is -Xmx4g):      ";
$mem = <STDIN>;
chomp $mem;
 
if ($mem eq ""){$mem = "-Xmx4g"}



#################################
# Annotate variants preferences #
#################################

print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " Would you like to annotate your SNP and INDEL calls using Ensembl?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";


print "   Enter 1 for YES\n";
print "   Enter 2 for NO\n\n";

$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1"){$annotate_variants = "yes"}
if (substr($answer,0,1) eq "2"){$annotate_variants = "no"}


###############################
# PINDEL analysis preferences #
###############################

if ($data eq "PE")
{

	print "\n\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print " Would you like to perform Structural Variant analysis?\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

	print "   Enter 1 for YES\n";
	print "   Enter 2 for NO\n\n";


	$answer = <STDIN>;
	chomp $answer;

	if (substr($answer,0,1) eq "1"){$structural_variant_analysis = "yes"}
	if (substr($answer,0,1) eq "2"){$structural_variant_analysis = "no"}


	if ($structural_variant_analysis eq "yes")
	{
		print "\nSTRUCTURAL VARIANT ANALYSIS WILL BE INCLUDED\n\n";
		print "\nPlease input the expected average insert size for your library (eg 200):      ";
		$insert = <STDIN>;
		chomp $insert;
		print "\n\n";
	}

	if ($structural_variant_analysis eq "no")
	{
		print "\nSTRUCTURAL VARIANT ANALYSIS WILL NOT BE INCLUDED\n";
	}



	if ($structural_variant_analysis eq "yes"){
		if ($ref_genome eq "partial"){
			print "You have indicated that the reference is a partial genome sequence\n";

			print "\nPlease enter the chromosome number you are interested in:    ";
			$chromosome = <STDIN>;
			chomp $chromosome;
		} # if ($structural_variant_analysis eq "yes")
	}

} # End of if ($data eq "PE")

	
$command =  "clear";
system("$command");


print color "bold white";

print "==================================================\n";
print "SUMMARY 1 - PLEASE CHECK THESE DETAILS CAREFULLY!!\n";
print "==================================================\n\n";

print color "bold green";
print "YOUR DETAILS :\n\n";

print color "bold white";
print "Name of this analysis:\t$run_title\n\n"; 

print "Your email address:\t$email_address\n\n";
print color "reset";

print "Please press enter to proceed (or 'Q' to quit):      ";
$proceed = <STDIN>;
chomp $proceed; 
 
if (lc $proceed eq "q"){exit;} 


print "\n\n";
print color "bold white";

print "==================================================\n";
print "SUMMARY 2 - PLEASE CHECK THESE DETAILS CAREFULLY!!\n";
print "==================================================\n\n";

print color "bold green";
print "SETTINGS :\n\n";
print color "bold white";


if ($data eq "SE")        {print "Single-end analysis    \t\t\tNO\n\n";}
if ($data eq "PE")        {print "Paired-end analysis    \t\t\tYES\n\n";}

if ($check_depth eq "yes"){print "Check depth of coverage\t\t\tYES\n\n";}
if ($check_depth eq "no") {print "Check depth of coverage\t\t\tNO\n\n";}


if ($maq eq "yes")        {print "Include MAQ analysis   \t\t\tYES\n\n"}
if ($maq eq "no")         {print "Include MAQ analysis   \t\t\tNO\n\n"}

if ($annotate_variants eq "yes")        {print "Annotate variants     \t\t\tYES\n\n";}
if ($annotate_variants eq "no")         {print "Annotate variants     \t\t\tNO\n\n";}

if ($structural_variant_analysis eq "no")          {print "Structural variant analysis\t\tNO\n\n";}
	
if ($structural_variant_analysis eq "yes")
{
	{print "Structural variant analysis\t\tYES\n\n";}

	if ($ref_genome eq "partial"){print "(Chromosome $chromosome)\n\n"}
	if ($ref_genome eq "whole"){print "\n\n"}
}
print "bwa settings:  use -s in bwa sampe:\t$use_s_option_in_bwa_sampe\n\n";

print "Memory setting: \t\t\t$mem\n\n";

print "UnifiedGenotyper constants:\n\n";

print "    -stand_emit_conf:\t\t\t$stand_emit_conf\n";
print "    -stand_call_conf:\t\t\t$stand_call_conf\n\n";

print color "reset";

print "Please press enter to proceed (or 'Q' to quit):      ";
$proceed = <STDIN>;
chomp $proceed; 
 
if (lc $proceed eq "q"){exit;} 

print "\n\n";

print color "bold white";

print "==================================================\n";
print "SUMMARY 3 - PLEASE CHECK THESE DETAILS CAREFULLY!!\n";
print "==================================================\n\n";

print color "bold green";
print "DATA FILES :\n\n";
print color "bold white";

if ($ref_genome eq "whole"){print "Whole genome reference\n\n";}

if ($ref_genome eq "partial"){print "Partial genome reference.  Chromosome $chromosome\n\n";}

print "Use pre-indexed Reference file at $ref\n\n";

print "Number of samples to analyse:\t\t$no_of_files\n\n";
	
for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print "Bam file $list_count	\t$bam_file_array[$list_count]\n";
}

print color "reset";

print "\nPlease press enter to proceed (or 'Q' to quit):      ";

print color "bold red";
print "\nNOTE - RUN TIME MAY BE SEVERAL HOURS.\n\n"; 

$proceed = <STDIN>;
chomp $proceed; 
 
if (lc $proceed eq "q")
{
	exit;
} 
		  


print color "reset";

#########################
# open Command Log file #
#########################
$command_log = "$run_title"."_command_log.txt";

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "COMMAND LOG for NGS pipeline version $version\n\n";
print COMMAND_LOG "Analysis name: $run_title\n\n";

print COMMAND_LOG "PREFERENCES: \n\n";

print COMMAND_LOG "Reference sequence:      \t$ref\n";
print COMMAND_LOG "Name of this analysis:   \t$run_title\n"; 
print COMMAND_LOG "E-mail address:          \t$email_address\n";
print COMMAND_LOG "Data type:               \t$data\n";
print COMMAND_LOG "Check depth of coverage: \t$check_depth\n";
print COMMAND_LOG "MAQ analysis:            \t$maq\n";
print COMMAND_LOG "Annotate variants:       \t$annotate_variants\n";
print COMMAND_LOG "Structural variants:     \t$structural_variant_analysis\n";
print COMMAND_LOG "Ref genome partial?      \t$ref_genome\n";
print COMMAND_LOG "Region:                  \t$region\n";
print COMMAND_LOG "Chromosome:              \t$chromosome\n";
print COMMAND_LOG "Memory setting:          \t$mem\n\n";
print COMMAND_LOG "UnifiedGenotyper constants:\n\n";
print COMMAND_LOG "    -stand_emit_conf:\t\t\t$stand_emit_conf\n";
print COMMAND_LOG "    -stand_call_conf:\t\t\t$stand_call_conf\n\n";
print COMMAND_LOG "Make parallel VCF:		\t$make_parallel_VCF\n";

print COMMAND_LOG "\nDATA FILES :\n\n";

print COMMAND_LOG "Number of samples to analyse:\t\t$no_of_files\n\n";
	
for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print COMMAND_LOG "Bam file $list_count	\t$bam_file_array[$list_count]\n";
}

	
##############################################################################
##############################################################################
# Start MAIN LOOP.  This runs through each of the files in the list of files #
##############################################################################
##############################################################################


for ($loop_count=$start_folder;$loop_count <=$no_of_files;$loop_count++)
{

	######################################
	# Create the results folder          #
	######################################

	$results_folder = "results_"."$run_title"."_"."$loop_count";
	$command = "mkdir $results_folder";

	print("$command\n");
	system("$command");

	
	
	######################################################
	# Assign correct name to $bam_file                   #
	######################################################
	$bam_file = $bam_file_array[$loop_count];
	
	
	######################################################
	# Create sample name (to use for name VCF files etc) #
	######################################################
	$sample_name = &get_prefix ($bam_file);
	
	

	#########################################################
	# Set up various File names (not all necessarily used)  #
	#########################################################
	$bai_file= "$sample_name".".bai";
	$bam_bai_file = "$sample_name".".bam.bai";
	$cleaned_sorted_bam = "$sample_name"."_cleaned_sorted.bam";
	$cleaned_sorted_bai = "$sample_name"."_cleaned_sorted.bai";
	$aligned_sam = "$sample_name"."_aligned.sam";
	$aligned_sorted_sam = "$sample_name"."_aligned_sorted.sam";
	$aligned_sorted_bam = "$sample_name"."_aligned_sorted.bam";
	$region_only_bam = "$sample_name"."_region.bam";
	$region_only_bai = "$sample_name"."_region.bai";
	$aligned_sorted_bai = "$sample_name"."_aligned_sorted.bai";
	$aligned_sorted_bam_bai = "$sample_name"."_aligned_sorted.bam.bai";
	$aligned_sorted_RG_bam = "$sample_name"."_aligned_sorted_RG.bam";
	$duplicates_removed_bam = "$sample_name"."_duplicates_removed.bam";
	$aligned_plus_dup_bam = "$sample_name"."_aligned_plus_dup.bam";
	$all_reads_bam = "$sample_name"."_all_reads.bam";
	$all_reads_bai = "$sample_name"."_all_reads.bai";
	$reads_csv = "$sample_name"."_reads.csv";
	$recal_bam = "$sample_name"."_recal.bam";
	$recal_bai = "$sample_name"."_recal.bai";
	$cleaned_bam = "$sample_name"."_cleaned.bam";
	$snps_raw_vcf = "$sample_name"."_snps_raw.vcf";
	$snps_raw_vcf_idx = "$sample_name"."_snps_raw.vcf.idx";
	$snps_parallel_vcf = "SNPS_"."$run_title".".vcf";
	$snps_final_vcf = "SNPS_"."$sample_name".".vcf";
	$indels_raw_vcf = "$sample_name"."_indels_raw.vcf";
	$indels_raw_vcf_alt = "$sample_name"."_indels_raw_alt.vcf";
	$indels_raw_vcf_idx = "$sample_name"."_indels_raw.vcf.idx";
	$indels_parallel_vcf = "INDELS_"."$run_title".".vcf";
	$indels_final_vcf = "INDELS_"."$sample_name".".vcf";
	$indels_raw_bed = "$sample_name"."_indels_raw.bed";
	$queryname_sorted_bam = "$sample_name"."_queryname_sorted.bam";
	$detailed_output_bed = "$sample_name"."_detailed_output.bed";
	$bam_file_mates_fixed = "$sample_name"."_mates_fixed.bam";
	$annotated_indels_vcf	= "annotated_indels_"."$sample_name.vcf";
	$annotated_snps_vcf	= "annotated_snps_"."$sample_name.vcf";
	$insert_size_pdf = "$sample_name"."_insert_size.pdf";
	$gc_bias_pdf = "$sample_name"."_gc_bias.pdf";
	
	#Broad file names specified
	$sai_coordinates_1 = "$sample_name.1.1.sai";
	$sai_coordinates_2 = "$sample_name.1.2.sai";
	$reverted_bam = "$sample_name.reverted.bam";
	$realigned_sam = "$sample_name.1.realigned.sam";
	$realigned_bam = "$sample_name.1.realigned.bam";
	$realigned_bai = "$sample_name.1.realigned.bai";
	$realigned_rg_bam  = "$sample_name.1.realigned.rg.bam";
	$realigned_rg_bai  = "$sample_name.1.realigned.rg.bai";
	$runtitle_sample_pre_validation = "$run_title"."_"."$sample_name"."_pre_validation.out";
	$runtitle_sample_post_validation = "$run_title"."_"."$sample_name".".post_validation.out";
	$runtitle_bam_intervals = "$run_title"."_"."$sample_name".".bam.intervals";
	$runtitle_clean_bam = "$run_title"."_"."$sample_name".".clean.bam";
	$runtitle_clean_bai = "$run_title"."_"."$sample_name".".clean.bai";
	$runtitle_clean_dedup_bam = "$run_title"."_"."$sample_name".".clean.dedup.bam";
	$runtitle_clean_dedup_bai = "$run_title"."_"."$sample_name".".clean.dedup.bai";
	$runtitle_clean_dedup_recal_bam =  "$run_title"."_"."$sample_name".".clean.dedup.recal.bam";
	$runtitle_clean_dedup_recal_bai =  "$run_title"."_"."$sample_name".".clean.dedup.recal.bai";
	$runtitle_sample_metrics = "$run_title"."_"."$sample_name".".metrics";
	$runtitle_sample_pre_recal_csv  = "$run_title"."_"."$sample_name".".pre_recal.csv";
	$runtitle_sample_post_recal_csv  = "$run_title"."_"."$sample_name".".post_recal.csv";
	$runtitle_sample_pre_recal_grp  = "$run_title"."_"."$sample_name".".pre_recal.grp";
	$runtitle_sample_post_recal_grp  = "$run_title"."_"."$sample_name".".post_recal.grp";
	$pdf_output_dir_pre = "$results_folder/$run_title"."_"."$sample_name"."_PDF_files_pre";
	$pdf_output_dir_post = "$results_folder/$run_title"."_"."$sample_name"."_PDF_files_post";
	
	

	##########################################################
	# If there isn't a BAM index file (.bai) then create one #
	##########################################################
	
	if (! -e "$bai_file"){print "$bai_file doesn't exist\n";}
	if (-e "$bai_file"){print "$bai_file exists\n";}
	if (! -e "$bam_bai_file"){print "$bam_bai_file doesn't exist\n";}
	if (-e "$bam_bai_file"){print "$bam_bai_file exists\n";}
	
	
	if ((! -e "$bai_file") && (! -e "$bam_bai_file"))
	{
		print "\n===================================================";
		print "\nFile $loop_count/$no_of_files   New Index for BAM file being created";
		print "\n===================================================\n\n\n";
		
		&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/BuildBamIndex.jar I=$bam_file VALIDATION_STRINGENCY=LENIENT","Make index on original BAM file");
		
		print "\n===============================================";
		print "\nFile $loop_count/$no_of_files   New Index for BAM file created";
		print "\n===============================================\n\n\n";
			
	}
	
	
	################################################################
	# If you are using a defined region then create a new smaller  #
	# BAM file with only this region in it                         #
	################################################################
	if ($use_defined_region eq "yes")
	{
		if ($make_smaller_bam_file eq "yes")
		{
			&run_unix_command("/opt/samtools/samtools view $bam_file $region -b -o $region_only_bam","0");

			&record_output_file_size ("$region_only_bam");	
			
			print "\n========================================================";
			print "\nFile $loop_count/$no_of_files   Step 0:  Bam file for region $region created";
			print "\n========================================================\n\n\n";
			
			
			########################################################
			# Now make an index file for this new smaller BAM file #
			########################################################
			
			print "\n============================================================";
			print "\nFile $loop_count/$no_of_files   New Index for region-only BAM file being created";
			print "\n============================================================\n\n\n";
		
		
			&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/BuildBamIndex.jar I=$region_only_bam VALIDATION_STRINGENCY=LENIENT","Make index");

			&record_output_file_size ("$region_only_bai");
			
			print "\n============================================================";
			print "\nFile $loop_count/$no_of_files   Step 0:  Bam Index file for region $region created";
			print "\n============================================================\n\n\n";
		
		} # end of if ($make_smaller_bam_file eq "yes")
		
		if ($make_smaller_bam_file eq "no")
		{
			$region_only_bam = $bam_file;
		}
	} # if ($use_defined_region eq "yes"}
	
	if ($use_defined_region eq "no")
	{
		$region_only_bam = $bam_file;
	}
	
	###############################################################
	# Broad stage 1:  Revert Bam files using picard/RevertBam.jar #
	###############################################################

	
	print	"\n====================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   Broad stage 1/17 Using picard/RevertSam to revert the BAM file to its original settings";
	print 	"\n====================================================================================================\n\n\n";
	
	
	&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/RevertSam.jar INPUT=$region_only_bam OUTPUT=$reverted_bam VALIDATION_STRINGENCY=SILENT  SO=queryname  CREATE_INDEX=true","1");
	
	
	print	"\n======================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   Broad stage 1/17 COMPLETED picard/RevertSam reverted the BAM file to its original settings";
	print 	"\n======================================================================================================\n\n";
			
	&record_output_file_size ("$reverted_bam");	
	&test_mode_subroutine;
	
	
	if ($data eq "PE")
	{
	##################################################################
	# Broad stage 2:  Generating SA coordinates of the input reads   #
	##################################################################
	
		print	"\n============================================================================";
		print 	"\nFile $loop_count/$no_of_files   Broad stage 2/17 generating SA coordinates of the input reads_1";
		print 	"\n============================================================================\n\n\n";
		
		&run_unix_command ("bwa aln -t 1 -q 5 $ref -b1 $reverted_bam > $sai_coordinates_1","2");

		print	"\n======================================================================================";
		print 	"\nFile $loop_count/$no_of_files   Broad stage 2/17 COMPLETED generating SA coordinates of the input reads_1";
		print 	"\n======================================================================================\n\n\n";
	
	&record_output_file_size ("$sai_coordinates_1");	
	&test_mode_subroutine;
	 
	##################################################################
	# Broad stage 3:  Generating SA coordinates of the input reads   #
	##################################################################
	
		print	"\n===========================================================================";
		print 	"\nFile $loop_count/$no_of_files   Broad stage 3/17 generating SA coordinates of the input reads_2";
		print 	"\n===========================================================================\n\n\n";
	
	
		&run_unix_command ("bwa aln -t 1 -q 5 $ref -b2 $reverted_bam > $sai_coordinates_2","3");

		print	"\n======================================================================================";
		print 	"\nFile $loop_count/$no_of_files   Broad stage 3/17 COMPLETED generating SA coordinates of the input reads_2";
		print 	"\n======================================================================================\n\n\n";
	
	&record_output_file_size ("$sai_coordinates_2");	
	&test_mode_subroutine;


	#####################################################################################
	# Broad stage 4: Run bwa sampe to generate alignments in SAM format (paired ends)   #
	#####################################################################################

		print	"\n====================================================================================";
		print 	"\nFile $loop_count/$no_of_files   Broad stage 4/17 Running bwa sampe to generate alignments in SAM format";
		print   "\n (running with 'use s option in bwa sampe' set to $use_s_option_in_bwa_sampe)";
		print 	"\n====================================================================================\n\n\n";
		
		if ($use_s_option_in_bwa_sampe eq "yes")
		{
			&run_unix_command ("bwa sampe -s $ref $sai_coordinates_1 $sai_coordinates_2 $reverted_bam $reverted_bam  > $realigned_sam","4");
		}
		if ($use_s_option_in_bwa_sampe eq "no")
		{
			&run_unix_command ("bwa sampe $ref $sai_coordinates_1 $sai_coordinates_2 $reverted_bam $reverted_bam  > $realigned_sam","4");
		}
		
		print	"\n====================================================================================";
		print 	"\nFile $loop_count/$no_of_files   Broad stage 4/17 COMPLETED bwa sampe generated alignments in SAM format";
		print   "\n (ran with 'use s option in bwa sampe' set to $use_s_option_in_bwa_sampe)";
		print 	"\n====================================================================================\n\n\n";
	
	&record_output_file_size ("$realigned_sam");
	&test_mode_subroutine;
	 
	} # End of if ($data eq "PE")
	
		
#####################################################################
# Broad stage 5: Use picard/SortSam to sort SAM file by coordinate  #
#####################################################################

print	"\n=================================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 5/17 Using picard/SortSam to sort SAM file by coordinate";
print 	"\n=================================================================================\n\n\n";

&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/SortSam.jar I=$realigned_sam O=$realigned_bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true","5");

print	"\n====================================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 5/17 COMPLETED picard/SortSam sorted SAM file by coordinate";
print 	"\n====================================================================================\n\n\n";
	
&record_output_file_size ("$realigned_bam");
&test_mode_subroutine;

 
##########################################################################################
# Broad stage 6: Use picard/AddOrReplaceReadGroups to update ReadGroup info in BAM file  #
##########################################################################################

print	"\n======================================================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 6/17 Using picard/AddOrReplaceReadGroups to update ReadGroup info in BAM file";
print 	"\n======================================================================================================\n\n\n";

&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/AddOrReplaceReadGroups.jar I=$realigned_bam O=$realigned_rg_bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true RGID=$sample_name  RGLB=$sample_name  RGPL=illumina  RGPU=$sample_name  RGSM=$sample_name","6");

print	"\n=======================================================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 6/17 COMPLETED picard/AddOrReplaceReadGroups updated ReadGroup info in BAM file";
print 	"\n=======================================================================================================\n\n\n";

&record_output_file_size ("$realigned_rg_bam");
&test_mode_subroutine;

 
##########################################################################################
# Broad stage 7: Use picard/ValidateSamFile to validate the BAM file                     #
##########################################################################################

if ($validate_sam_files eq "yes")
{
	print	"\n==================================================================================";
	print 	"\nFile $loop_count/$no_of_files   Broad stage 7/17 Using picard/ValidateSamFile to validate the BAM file";
	print 	"\n==================================================================================\n\n\n";
	
	&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/ValidateSamFile.jar INPUT=$realigned_rg_bam OUTPUT=$runtitle_sample_pre_validation VALIDATION_STRINGENCY=SILENT  MODE=VERBOSE  MAX_OUTPUT=100  MAX_OPEN_TEMP_FILES=$max_open_temp_files","7");

	print	"\n==================================================================================";
	print 	"\nFile $loop_count/$no_of_files   Broad stage 7/17 COMPLETED picard/ValidateSamFile validated the BAM file";
	print 	"\n==================================================================================\n\n\n";

	&record_output_file_size ("$runtitle_sample_pre_validation");
	&test_mode_subroutine;
}

##########################################################################################################
# Broad stage 8: Use GenomeAnalysisTK/RealignerTargetCreator to identify intervals in need of realignment #
##########################################################################################################

print	"\n=======================================================================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 8/17 Use GenomeAnalysisTK/RealignerTargetCreator to identify intervals in need of realignment";
print 	"\n=======================================================================================================================\n\n\n";


#&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref -I $realigned_rg_bam -o $runtitle_bam_intervals -mismatch 0.0","8");
&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T RealignerTargetCreator -R $ref $chromosome_string -I $realigned_rg_bam -o $runtitle_bam_intervals -mismatch 0.0","8");


print	"\n===============================================================================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 8/17 COMPLETED GenomeAnalysisTK/RealignerTargetCreator has identified intervals in need of realignment";
print 	"\n===============================================================================================================================\n\n\n";

&record_output_file_size ("$runtitle_bam_intervals");
&test_mode_subroutine;

	
##########################################################################################################
# Broad stage 9: Use GenomeAnalysisTK/IndelRealigner to Realign Indels                                   #
##########################################################################################################

print	"\n====================================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 9/17 Using GenomeAnalysisTK/IndelRealigner to Realign Indels";
print 	"\n====================================================================================\n\n\n";


&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T IndelRealigner -R $ref -I $realigned_rg_bam -targetIntervals $runtitle_bam_intervals -o $runtitle_clean_bam -compress 0 -model USE_READS","9");

print	"\n==============================================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 9/17 COMPLETED using GenomeAnalysisTK/IndelRealigner to Realign Indels";
print 	"\n==============================================================================================\n\n\n";

&record_output_file_size ("$runtitle_clean_bam");
&test_mode_subroutine;


##########################################################################################################
# Broad stage 10: Use picard/MarkDuplicates to mark duplicates                                            #
#                Cleaned BAM --> Dedup BAM
##########################################################################################################

if ($remove_duplicates eq "yes")
{

	if ($use_samtools_rmdup eq "no")
	{
		print	"\n==========================================================================";
		print 	"\nFile $loop_count/$no_of_files   Broad stage 10/17 Use picard/MarkDuplicates to mark duplicates";
		print 	"\n==========================================================================\n\n\n";

		&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/MarkDuplicates.jar I=$runtitle_clean_bam O=$runtitle_clean_dedup_bam  VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true M=$runtitle_sample_metrics","10");

		print	"\n======================================================================================";
		print 	"\nFile $loop_count/$no_of_files   Broad stage 10/17 COMPLETED using picard/MarkDuplicates to mark duplicates";
		print 	"\n======================================================================================\n\n\n";
		
	} # don't use samtools rmdup - use MArkDuplicates
	
		if ($use_samtools_rmdup eq "yes")
	{
		print	"\n==========================================================================";
		print 	"\nFile $loop_count/$no_of_files   (non) Broad stage 10/17 Use samtools rmdup to remove duplicates";
		print 	"\n==========================================================================\n\n\n";

		&run_unix_command (" /opt/samtools/samtools rmdup $runtitle_clean_bam runtitle_clean_dedup_bam","10");

		print	"\n======================================================================================";
		print 	"\nFile $loop_count/$no_of_files   (non) Broad stage 10/17 COMPLETED using samtools rmdup to remove duplicates";
		print 	"\n======================================================================================\n\n\n";
		
	} # do use samtools rather than picard to remove duplicates

}
if ($remove_duplicates eq "no")
{
	$runtitle_clean_dedup_bam = $runtitle_clean_bam;
}

&record_output_file_size ("$runtitle_clean_dedup_bam");
&test_mode_subroutine;

	
##########################################################################################################
# Broad stage 11: Use GenomeAnalysisTK/BaseRecalibrator to generate a recalibration table   (pre)        #
##########################################################################################################

&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T BaseRecalibrator -R $ref -I $runtitle_clean_dedup_bam -knownSites $dummy_dbsnp_file -o  $runtitle_sample_pre_recal_grp","11");

print	"\n================================================================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 11/17 Used GenomeAnalysisTK/CountCovariates to count covariates ";
print 	"\n================================================================================================================\n\n\n";

&record_output_file_size ("$runtitle_sample_pre_recal_csv");
&test_mode_subroutine;



##########################################################################################################
# Broad stage 12: Use GenomeAnalysisTK/TableRecalibration to update the base quality scores              #
##########################################################################################################

&run_unix_command ("java  $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T TableRecalibration -I $runtitle_clean_dedup_bam  -R $ref -baq CALCULATE_AS_NECESSARY  -recalFile $runtitle_sample_pre_recal_csv  -o $runtitle_clean_dedup_recal_bam","12");

print	"\n===============================================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 12/17 Used GenomeAnalysisTK/TableRecalibration to update the base quality scores ";
print 	"\n================================================================================================\n\n\n";

&record_output_file_size ("$runtitle_clean_dedup_recal_bam");
&test_mode_subroutine;

	
##########################################################################################################
# Broad stage 13: Use AnalyzeCovariates to plot some PDF files (PRE)                                     #
##########################################################################################################

&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/AnalyzeCovariates.jar -recalFile $runtitle_sample_pre_recal_csv -outputDir $pdf_output_dir_pre","13");

print	"\n========================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 13/17 Used AnalyzeCovariates to plot some PDF files (pre)";
print 	"\n========================================================================\n\n\n";

&record_output_file_size ("$pdf_output_dir_pre");
&test_mode_subroutine;

	
##########################################################################################################
# Broad stage 14: Use picard/ValidateSamFile to validate the cleaned BAM file                            #
##########################################################################################################

if ($validate_sam_files eq "yes")
{
	&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/ValidateSamFile.jar INPUT=$runtitle_clean_dedup_recal_bam OUTPUT=$runtitle_sample_post_validation VALIDATION_STRINGENCY=SILENT MODE=VERBOSE MAX_OUTPUT=100 MAX_OPEN_TEMP_FILES=$max_open_temp_files","14");

	print	"\n===============================================================================";
	print 	"\nFile $loop_count/$no_of_files   Broad stage 14/17 Used picard/ValidateSamFile to validate the cleaned BAM file";
	print 	"\n===============================================================================\n\n\n";

	&record_output_file_size ("$runtitle_sample_post_validation");
	&test_mode_subroutine;
}
	
##########################################################################################################
# Broad stage 15: Use GenomeAnalysisTK/CountCovariates to count covariates (post)                        #
##########################################################################################################

&run_unix_command ("java $mem $temp_dir_string -jar /opt/gatk_v1/GenomeAnalysisTK.jar -T CountCovariates -I $runtitle_clean_dedup_recal_bam -R $ref -knownSites $dummy_dbsnp_file -recalFile $runtitle_sample_post_recal_csv -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate","15");
 
print	"\n========================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 15/17 Used GenomeAnalysisTK/CountCovariates to count covariates (post)";
print 	"\n========================================================================\n\n\n";

&record_output_file_size ("$runtitle_sample_post_recal_csv");
&test_mode_subroutine;


##########################################################################################################
# Broad stage 16: Do something with the list of files                                                    #
##########################################################################################################

#.qlog/QTEST2.cohort.list.bamList: List(QTEST2.test1.clean.dedup.recal.bam) > List(QTEST2.cohort.list)

print	"\n========================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 16/17 Should have done something with the list of files (not yet fixed) ";
print	"\n========================================================================";

&test_mode_subroutine;


##########################################################################################################
# Broad stage 17: Use AnalyzeCovariates to plot some PDF files again (POST)                              #
##########################################################################################################

&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/AnalyzeCovariates.jar -recalFile $runtitle_sample_post_recal_csv -outputDir $pdf_output_dir_post","17");

print	"\n=========================================================================";
print 	"\nFile $loop_count/$no_of_files   Broad stage 17/17 Used AnalyzeCovariates to plot some PDF files (post)";
print 	"\n=========================================================================\n\n\n";

&record_output_file_size ("$pdf_output_dir_post");
&test_mode_subroutine;



########################################################################################
# 18 Making INDEL calls  
########################################################################################

if ($make_parallel_VCF eq "consecutive" || $make_parallel_VCF eq "parallel and consecutive")
{

	if ($indels_caller eq "SomaticIndelDetector")
	{
		
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T SomaticIndelDetector -R $ref -I $runtitle_clean_dedup_recal_bam -o $indels_raw_vcf  -ws 300 --unpaired","INDEL CALLS");

		&record_output_file_size ("$indels_raw_vcf");
		&test_mode_subroutine;

	}
	
	if ($indels_caller eq "UnifiedGenotyper")
	{

		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T UnifiedGenotyper $chromosome_string -glm INDEL -R $ref -I $runtitle_clean_dedup_recal_bam -o $indels_raw_vcf  -S $ValidationStringency","INDEL CALLS");

		&record_output_file_size ("$indels_raw_vcf");
		&test_mode_subroutine;

	}
	

	#~~~~~~~~~~~~~~~#
	#ERROR CHECK v2!#
	#~~~~~~~~~~~~~~~#
	&problem_check;


	print "\n================================================";
	print "\nFile $loop_count/$no_of_files   Step 18 Indel calls have been made";
	print "\n================================================\n\n\n";

	&test_mode_subroutine;
	
} # end of if ($make_parallel_VCF eq "consecutive" || $make_parallel_VCF eq "parallel and consecutive")


########################################################################################
# 19 Making SNP calls   
########################################################################################


if ($make_parallel_VCF eq "consecutive" || $make_parallel_VCF eq "parallel and consecutive")
{

	&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $chromosome_string -glm SNP -I $runtitle_clean_dedup_recal_bam -stand_call_conf $stand_call_conf -stand_emit_conf $stand_emit_conf -o $snps_raw_vcf  -S $ValidationStringency");


	print "\n==============================================";
	print "\nFile $loop_count/$no_of_files   Step 19 SNP Calls have been made";
	print "\n==============================================\n";

	&record_output_file_size ("$snps_raw_vcf");
	&test_mode_subroutine;

} # End of if ($make_parallel_VCF eq "consecutive" || $make_parallel_VCF eq "parallel and consecutive")


	########################################
	# Structural variant analysis - Pindel #
	########################################


	if ($structural_variant_analysis eq "yes")
	{

		if ($ref_genome eq "partial")
		{

			print "\n\n";
			print "===============================================\n";
			print "Starting script to find structural variants...\n";
			print "===============================================\n\n";

			print "Sorting Bam file by query name...";

			&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/SortSam.jar I=$cleaned_sorted_bam O=$queryname_sorted_bam SO=queryname VALIDATION_STRINGENCY=LENIENT","Sorting BAM file");

			print "\n==============================================\n";
			print "Bam sorted, now making input files for pindel \n";
			print "==============================================\n\n";

			&record_output_file_size ("$queryname_sorted_bam");
			&test_mode_subroutine;
	
	
			$command =  "perl /opt/pindel/bam2pindel.pl -i $queryname_sorted_bam -o sv -s $sample_name -pi $insert";

			print("$command\n");
			system("$command");
			print COMMAND_LOG "$command\n";

			print "\n\nRunning Pindel\n\n";

			
			$command =  "/opt/pindel/pindel022 -f $ref -p sv_chr$chromosome.txt -c chr$chromosome -o pindel";

			print("$command\n");
			system("$command");
			print COMMAND_LOG "$command\n";

		}


		if ($ref_genome eq "whole")
		{

			print "\n\nStarting script to find structural variants...";
			print "\n\nSorting Bam file by query name...";

			$command =  "java $mem $temp_dir_string -jar /opt/picard/SortSam.jar I=$cleaned_sorted_bam O=$queryname_sorted_bam SO=queryname VALIDATION_STRINGENCY=LENIENT";
			print("$command\n");
			system("$command");
			print COMMAND_LOG "$command\n";

			print "\n\nBam sorted, now making input files for pindel\n\n";

			$command =  "perl /opt/pindel/bam2pindel.pl -i $queryname_sorted_bam -o sv -s $sample_name -pi $insert";
			print("$command\n");
			system("$command");
			print COMMAND_LOG "$command\n";

			print "\n\nMerging all Pindel files\n\n";

			$command =  "cat sv_chr1.txt sv_chr2.txt sv_chr3.txt sv_chr4.txt sv_chr5.txt sv_chr6.txt sv_chr7.txt sv_chr8.txt sv_chr9.txt sv_chr10.txt sv_chr11.txt sv_chr12.txt sv_chr13.txt sv_chr14.txt sv_chr15.txt sv_chr16.txt sv_chr17.txt sv_chr18.txt sv_chr19.txt sv_chr20.txt sv_chr21.txt sv_chr22.txt sv_chr23.txt sv_chr24.txt sv_chr25.txt sv_chr26.txt sv_chr27.txt sv_chr28.txt sv_chr29.txt sv_chr30.txt sv_chr31.txt sv_chr32.txt sv_chr33.txt sv_chr34.txt sv_chr35.txt sv_chr36.txt sv_chr37.txt sv_chr38.txt sv_chrX.txt sv_chrM.txt sv_chrUn.txt>all.txt";
			print("$command\n");
			system("$command");
			print COMMAND_LOG "$command\n";

			print "\n\nRunning Pindel\n\n";

			$command =  "/opt/pindel/pindel022 -f $ref -p all.txt -c ALL -o pindel";
			print("$command\n");
			system("$command");
			print COMMAND_LOG "$command\n";

		}

	} # End of if ($structural_variant_analysis eq "yes")


	###############################
	# RUN PINDEL TO VCF CONVERTER #
	# (CANFAM ONLY AT PRESENT)    #
	###############################


	if ($structural_variant_analysis eq "yes")
	{

		$command =  "/opt/pindel/pindel2vcf -p pindel_BP  -r $ref -R canfam2 -d 2006";
		print("$command\n");
		system("$command");
		print COMMAND_LOG "$command\n";

		$command =  "/opt/pindel/pindel2vcf -p pindel_D  -r $ref -R canfam2 -d 2006";
		print("$command\n");
		system("$command");
		print COMMAND_LOG "$command\n";

		$command =  "/opt/pindel/pindel2vcf -p pindel_INV  -r $ref -R canfam2 -d 2006";
		print("$command\n");
		system("$command");
		print COMMAND_LOG "$command\n";

		$command =  "/opt/pindel/pindel2vcf -p pindel_LI  -r $ref -R canfam2 -d 2006";
		print("$command\n");
		system("$command");
		print COMMAND_LOG "$command\n";

		$command =  "/opt/pindel/pindel2vcf -p pindel_SI  -r $ref -R canfam2 -d 2006";
		print("$command\n");
		system("$command");
		print COMMAND_LOG "$command\n";

		$command =  "/opt/pindel/pindel2vcf -p pindel_TD  -r $ref -R canfam2 -d 2006";
		print("$command\n");
		system("$command");
		print COMMAND_LOG "$command\n";

		print "\n\nVCF FILES MADE FOR PINDEL RESULTS\n\n";

	} # End of if ($structural_variant_analysis eq "yes")



	####################################################################
	# Checking depth of coverage, average insert size and GC bias
	####################################################################

	if ($check_depth eq "yes")
	{
		print "\n Checking depth of coverage....\n\n";

		$command = "java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T DepthOfCoverage -R $ref -I $runtitle_clean_dedup_recal_bam -o depth  -omitBaseOutput -S $ValidationStringency";

		print("$command\n");
		system("$command");
		print COMMAND_LOG "$command\n";
	}

	print "\n========================================";
	print "\n Calculating GC Bias....";
	print "\n========================================\n\n";

	&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/CollectGcBiasMetrics.jar REFERENCE_SEQUENCE=$ref INPUT=$runtitle_clean_dedup_recal_bam O=out1.junk CHART_OUTPUT=$gc_bias_pdf VALIDATION_STRINGENCY=LENIENT","GC_BIAS");

	&record_output_file_size("$gc_bias_pdf");
	
	print "\n========================================";
	print "\nFile $loop_count/$no_of_files   GC bias has been calculated";
	print "\n========================================\n";


	
	##############################
	# Plot Insert Size Histogram #
	##############################
	print "\n======================================================";
	print "\n Plotting insert size histogram (for PE data only)...";
	print "\n======================================================\n\n";


	&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/CollectInsertSizeMetrics.jar INPUT=$runtitle_clean_dedup_recal_bam O=out2.junk HISTOGRAM_FILE=$insert_size_pdf VALIDATION_STRINGENCY=LENIENT","SIZE_METRICS");

	&record_output_file_size("$insert_size_pdf");



	print "\n=========================================";
	print "\nFile $loop_count/$no_of_files   Insert size histogram plotted";
	print "\n=========================================\n";




	if ($structural_variant_analysis eq "yes")
	{
		$command = "mkdir $results_folder/Raw_SV_data";

		print("$command\n");
		system("$command");
		print COMMAND_LOG "$command\n";
	}


	print "\nSorting results into new directory....\n\n";



	########################################
	# Some files don't need to be kept     #
	# (keep files for first loop to check) #
	########################################
	
	if ($loop_count == 1){$delete_intermediates = "no"}else{$delete_intermediates = "yes"}
	
	if ($delete_intermediates eq "yes")
	{
		&delete_file("$realigned_sam");
		&delete_file("$reverted_bam");
		&delete_file("$reverted_bai");
		&delete_file("$realigned_bam");
		&delete_file("$realigned_bai");
		&delete_file("$runtitle_clean_dedup_bam");
		&delete_file("$runtitle_clean_dedup_bai");
		&delete_file("$runtitle_clean_bam");
		&delete_file("$runtitle_clean_bai");
		&delete_file("$realigned_rg_bam");
		&delete_file("$realigned_rg_bai");
		
		&delete_file("$sai_coordinates_1");
		&delete_file("$sai_coordinates_2");
		
		&delete_file("$runtitle_bam_intervals");

	}
	
	##########################################
	# Move various files to a Results Folder #
	##########################################
	
	&move_to_results_folder ("$runtitle_clean_dedup_recal_bam");
	&move_to_results_folder ("$runtitle_clean_dedup_recal_bai");
	&move_to_results_folder ("$runtitle_sample_pre_validation");
	&move_to_results_folder ("$runtitle_sample_post_validation");
	&move_to_results_folder ("$runtitle_sample_metrics");
	&move_to_results_folder ("$runtitle_sample_pre_recal_csv");
	&move_to_results_folder ("$runtitle_sample_post_recal_csv");

	&move_to_results_folder ("$snps_raw_vcf");
	&move_to_results_folder ("$aligned_sorted_bai");
	&move_to_results_folder ("$aligned_sorted_bam");
	&move_to_results_folder ("$indels_raw_vcf");
	&move_to_results_folder ("depth");
	&move_to_results_folder ("depth.sample_summary");
	&move_to_results_folder ("$gc_bias_pdf");
	&move_to_results_folder ("$insert_size_pdf");
	&move_to_results_folder ("Duplicate_info.rtf");
	&move_to_results_folder ("Alignment_Summary.xls");



	#####################################
	# Change the names of various files #
	#####################################

	$command = "mv $results_folder/$indels_raw_vcf $results_folder/$indels_final_vcf";
	system("$command");

	$command = "mv $results_folder/$snps_raw_vcf $results_folder/$snps_final_vcf";
	system("$command");

	$command = "mv $results_folder/depth $results_folder/depth_bp.xls";
	system("$command");

	$command = "mv $results_folder/depth.sample_summary $results_folder/depth_summary.xls";
	system("$command");


	##############################################################################
	# If Structural Variants has been done, move the files to the Results Folder #
	##############################################################################
	if ($structural_variant_analysis eq "yes")
	{

		&delete_file ("queryname_sorted.bam");

		&move_to_results_folder ("pindel_BP");
		&move_to_results_folder ("pindel_D");
		&move_to_results_folder ("pindel_INV");
		&move_to_results_folder ("pindel_LI");
		&move_to_results_folder ("pindel_SI");
		&move_to_results_folder ("pindel_TD");


		$command = "mv $results_folder/pindel_BP $results_folder/Raw_SV_data/SV_Break_points.txt";
		system("$command");

		$command = "mv $results_folder/pindel_D $results_folder/Raw_SV_data/SV_Deletions.txt";
		system("$command");

		$command = "mv $results_folder/pindel_INV $results_folder/Raw_SV_data/SV_Inversions.txt";
		system("$command");

		$command = "mv $results_folder/pindel_LI $results_folder/Raw_SV_data/SV_Long_insertions.txt";
		system("$command");

		$command = "mv $results_folder/pindel_SI $results_folder/Raw_SV_data/SV_Non_template_seq_in_deleletion.txt";
		system("$command");

		$command = "mv $results_folder/pindel_TD $results_folder/Raw_SV_data/SV_Tandom_Dup.txt";
		system("$command");


		&move_to_results_folder ("pindel_BP.vcf");
		&move_to_results_folder ("pindel_D.vcf");
		&move_to_results_folder ("pindel_INV.vcf");
		&move_to_results_folder ("pindel_LI.vcf");
		&move_to_results_folder ("pindel_SI.vcf");
		&move_to_results_folder ("pindel_TD.vcf");


		$command = "mv $results_folder/pindel_BP.vcf $results_folder/SV_Break_points.vcf";
		system("$command");

		$command = "mv $results_folder/pindel_D.vcf $results_folder/SV_Deletions.vcf";
		system("$command");

		$command = "mv $results_folder/pindel_INV.vcf $results_folder/SV_Inversions.vcf";
		system("$command");

		$command = "mv $results_folder/pindel_LI.vcf $results_folder/SV_Long_insertions.vcf";
		system("$command");

		$command = "mv $results_folder/pindel_SI.vcf $results_folder/SV_Non_template_seq_in_deleletion.vcf";
		system("$command");

		$command = "mv $results_folder/pindel_TD.vcf $results_folder/SV_Tandom_Dup.vcf";
		system("$command");


		if ($ref_genome eq "partial")
		{
			&delete_file ("sv_chr$chromosome.txt");
		}

		if ($ref_genome eq "whole")
		{
			&delete_file ("all.txt");
			&delete_file ("sv_chr1.txt");
			&delete_file ("sv_chr2.txt");
			&delete_file ("sv_chr3.txt");
			&delete_file ("sv_chr4.txt");
			&delete_file ("sv_chr5.txt");
			&delete_file ("sv_chr6.txt");
			&delete_file ("sv_chr7.txt");
			&delete_file ("sv_chr8.txt");
			&delete_file ("sv_chr9.txt");
			&delete_file ("sv_chr10.txt");
			&delete_file ("sv_chr11.txt");
			&delete_file ("sv_chr12.txt");
			&delete_file ("sv_chr13.txt");
			&delete_file ("sv_chr14.txt");
			&delete_file ("sv_chr15.txt");
			&delete_file ("sv_chr16.txt");
			&delete_file ("sv_chr17.txt");
			&delete_file ("sv_chr18.txt");
			&delete_file ("sv_chr19.txt");
			&delete_file ("sv_chr20.txt");
			&delete_file ("sv_chr21.txt");
			&delete_file ("sv_chr22.txt");
			&delete_file ("sv_chr23.txt");
			&delete_file ("sv_chr24.txt");
			&delete_file ("sv_chr25.txt");
			&delete_file ("sv_chr26.txt");
			&delete_file ("sv_chr27.txt");
			&delete_file ("sv_chr28.txt");
			&delete_file ("sv_chr29.txt");
			&delete_file ("sv_chr30.txt");
			&delete_file ("sv_chr31.txt");
			&delete_file ("sv_chr32.txt");
			&delete_file ("sv_chr33.txt");
			&delete_file ("sv_chr34.txt");
			&delete_file ("sv_chr35.txt");
			&delete_file ("sv_chr36.txt");
			&delete_file ("sv_chr37.txt");
			&delete_file ("sv_chr38.txt");
			&delete_file ("sv_chrX.txt");
			&delete_file ("sv_chrM.txt");
			&delete_file ("sv_chrUn.txt");

		}

	} # End of if ($structural_variant_analysis eq "yes")

	
	##########
	# README #
	##########

	open (MYFILE, '>>README.rtf'); 

	print	MYFILE	"========================\n";
	print	MYFILE	"Summary of results files\n";
	print	MYFILE	"========================\n\n\n";

	print	MYFILE	"PDF FILES\n\n";

	print	MYFILE	"$gc_bias_pdf          \tHistogram of GC content of aligned reads\n";
	print	MYFILE	"$insert_size_pdf      \tInsert size histogram\n\n";

	print	MYFILE	"ALIGNMENT FILES\n\n";

	print	MYFILE	"$runtitle_clean_dedup_recal_bam	Best alignments to the reference after processing by GATK\n\n";

	print 	MYFILE	"SNP AND INDEL CALLS\n\n";

	print	MYFILE	"$indels_final_vcf      \tList of InDels\n";
	print	MYFILE	"$snps_final_vcf        \tList of SNPs\n";
	print	MYFILE	"annotated_indels.xls   \tList of InDels annotated using the ensembl database\n";
	print	MYFILE	"annotated_snp.xls      \tList of SNPs annotated using the ensembl database\n\n";


	print 	MYFILE	"STRUCTURAL VARIANT FILES\n\n";

	print 	MYFILE	"SV_Break_points.txt\n";
	print 	MYFILE	"SV_Deletions.txt\n";
	print 	MYFILE	"SV_Inversions.txt\n";
	print 	MYFILE	"SV_Long_insertions.txt\n";
	print 	MYFILE	"SV_Non_template_seq_in_deleletion.txt\n";
	print 	MYFILE	"SV_Tandom_Dup.txt\n\n";

	print	MYFILE	"INFORMATION FILES\n\n";	

	print	MYFILE	"Depth_summary	\tSummary of reads depth across the target region\n";
	print	MYFILE	"Log.rtf	\t\tRun log for the NGS pipeline\n";
	print	MYFILE	"Duplicate_info.rtf	fraction of PCR duplicates reads in the dataset\n\n";

	print	MYFILE	"NOTES\n\n"; 
	print	MYFILE	".bam files must have and associated .bai index file for loading into IGV\n";
	print	MYFILE	"Alignment_Summary.xls can be used to calculate success of target enrichment\n\n";

	close (MYFILE);

	&move_to_results_folder ("README.rtf");

	#########################################################################################


	##############################################
	# Delete these files at the end of each loop #
	##############################################
	&delete_file ("$aligned_sam");
	&delete_file ("$aln_sa_sai");
	&delete_file ("$aligned_sorted_sam");
	&delete_file ("$reads_csv");
	&delete_file ("$recal_bam");
	&delete_file ("$recal_bai");
	&delete_file ("forRealigner.intervals");
	&delete_file ("$cleaned_bam");
	&delete_file ("$detailed_output_bed");
	&delete_file ("$aln_sa2_sai");
	&delete_file ("$aln_sa1_sai");
	&delete_file ("depth.sample_cumulative_coverage_counts");
	&delete_file ("depth.sample_cumulative_coverage_proportions");
	&delete_file ("depth.sample_statistics");
	##&delete_file ("*.dict");  # <<<<<<<<<<<<<<<<<<<<<< check!
	&delete_file ("*.junk");
	&delete_file ("$snps_raw_vcf_idx");
	&delete_file ("$indels_raw_vcf_idx");



	###########################
	#                         # 
	#   ANNOTATING VARIANTS   #
	#                         #
	###########################


	if ($annotate_variants eq "yes")
	{
		print "\n\n";
		print "############################\n";
		print "#  ANNOTATING VARIANTS...  #\n";
		print "############################\n\n";

		$command = ". ~/.bashrc";

		print("$command\n\n");
		system("$command");

		#################
		# Annotate SNPs #
		#################
		&run_unix_command("perl /opt/ensembl/variant_effect_predictor.pl --species $species --input $results_folder/$snps_final_vcf --format VCF --output_file $annotated_snps_vcf -v","ANNOTATE SNPS");

		&record_output_file_size ("$annotated_snps_vcf");
			
		###################
		# Annotate Indels #
		###################
		&run_unix_command("perl /opt/ensembl/variant_effect_predictor.pl --species $species --input $results_folder/$indels_final_vcf --format pileup --output_file annotated_indels_vcf -v","ANNOTATE INDELS");

		&record_output_file_size ("$annotated_indels_vcf");
		
		print "\n\nVariant annotation complete\n\n";

		&move_to_results_folder ("$annotated_snps_vcf");
		&move_to_results_folder ("$annotated_indels_vcf");

	}

	#RUN TIMER END

	
	##########################################
	# Send an e-mail at the end of each loop #
	##########################################
	open(MAIL, "|/usr/sbin/sendmail -t");

	## Mail Header
	print MAIL "To: $email_address\n";
	print MAIL "From: $from\n";
	print MAIL "Subject: NGS ANALYSIS: Sample $sample_name in run $run_title has finished\n\n";
	## Mail Body
	print MAIL "NGS PERL script version $version\n\n";
	print MAIL "Your next generation sequence analysis of $sample_name in run $run_title is complete\n\n";
	print MAIL "Run time so far : $run_time seconds\n";
	printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
	close(MAIL);


} # End of MAIN LOOP for ($loop_count=1;$loop_count <=$no_of_files;$loop_count++)


##########################################################
#  Make a VCF file using all the BAM files in parallel   #
##########################################################

####################################################
# Open the list file to get the list of file names #
# (seem to lose the file names so try this)        #
####################################################

print "\n\nInput file: $input_file\n\n";
print "List file: $list_file\n\n";

if ($input_file eq "bam")
{
	print "Opening file of file names...\n\n";
	
	open (LIST, "$list_file") || die "Cannot open $list_file";
	$list_count=1;

	while ($bam_file = <LIST> ) 
	{
		chomp $bam_file;
		
		######################################################
		# Add file type suffix .bam if user hasn't added it  #
		######################################################
		$bam_file_array[$list_count]=$bam_file;
		$list_count=$list_count + 1;
	}

close LIST;
}


###################
# List file names #
###################
print "\n\nThere are $no_of_files files in this file of file names.\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print "File $list_count	\t$bam_file_array[$list_count]\n";
}

print "\nIf these file names look OK, press enter to proceed (or 'Q' to quit):      ";
#$proceed = <STDIN>;
$proceed = "";
chomp $proceed; 
 
if (lc $proceed eq "q"){exit;} 

		
print "\n\n#--------------------------------------------#\n";
print "#  Making VCF with all BAM files in parallel #\n";
print "#--------------------------------------------#\n\n";
if ($make_parallel_VCF eq "parallel" || $make_parallel_VCF eq "parallel and consecutive")
{

	# First get a list of the BAM files #
	for ($list_count=1;$list_count <=$no_of_files;$list_count++)
	{
	
		print "LIST COUNT: $list_count\n";
		
		$results_folder = "results_"."$run_title"."_"."$list_count";
		print "RESULTS FOLDER: $results_folder\n";
		
		$bam_file = $bam_file_array[$list_count];
		print "BAM FILE: $bam_file\n";
		
		$sample_name = &get_prefix ($bam_file);
		print "SAMPLE NAME: $sample_name\n";
		

		################################################################
		# Make up name of cleaned dedupped recalibrated final bam file #
		################################################################
		$runtitle_clean_dedup_recal_bam =  "$run_title"."_"."$sample_name".".clean.dedup.recal.bam";
			
		$bam_file_for_unified_genotyper = "$results_folder/$runtitle_clean_dedup_recal_bam";
		
		$input_string = $input_string." -I $bam_file_for_unified_genotyper";
	}
	$input_string = $input_string." ";

	&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $chromosome_string -glm SNP $input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o $snps_parallel_vcf  -S $ValidationStringency","PARALLEL SNPS VCF");	

	&record_output_file_size ("$snps_parallel_vcf");	
	
	print "\n===================================================================";
	print "\nA SNP VCF file using all the BAM files in parallel has been created";
	print "\n===================================================================\n\n\n";

	&test_mode_subroutine;
	
	&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $chromosome_string -glm INDEL $input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o  $indels_parallel_vcf -S $ValidationStringency","PARALLEL INDELS VCF");
	
	&record_output_file_size ("$indels_parallel_vcf");

	
	print "\n=====================================================================";
	print "\nAn Indel VCF file using all the BAM files in parallel has been created";
	print "\n=====================================================================\n\n\n";

	&test_mode_subroutine;
	
	#####################################################
	# Move parallel VCF file to the last results folder #
	#####################################################
	&move_to_results_folder ("$snps_parallel_vcf");
	&move_to_results_folder ("$indels_parallel_vcf");
	
	
} # End of if ($make_parallel_VCF eq "parallel" || $make_parallel_VCF eq "parallel and consecutive")



my $end_run = time();
$run_time = $end_run - our $start_run;

open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $email_address\n";
print MAIL "From: $from\n";
print MAIL "Subject: $subject COMPLETE ($run_title)\n\n";
## Mail Body
print MAIL "NGS PERL script version $version\n\n";
print MAIL "Your next generation sequence analysis ($run_title) is complete\n\n";
print MAIL "For pipeline details see log.rtf\n\n";
print MAIL "For file information see README.rtf\n\n";
print MAIL "For a list of commands see $command_log\n\n";
print MAIL "Run time : $run_time seconds\n";
printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
close(MAIL);

print "\n\n";
print "######################################################\n";
print "# ANALYSIS COMPLETE! YOU HAVE BEEN NOTIFIED BY EMAIL #\n";
print "######################################################\n\n";

print "Run time : $run_time seconds\n";
printf "%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];


print "\n\nPlease check log.rtf for full details\n\n";


#############################
# Turn off logging          #
#############################

close(STDOUT);
close (COMMAND_LOG);

&move_log ("log.rtf");
&move_log ("$command_log");



#############################################
#                                           #
# Subroutine to move file to results folder #
#                                           #
#############################################

sub move_to_results_folder
{

	my $file_to_be_moved = "";	

	$file_to_be_moved = $_[0];
	$command = "mv  $file_to_be_moved $results_folder/$file_to_be_moved";
	print("$command\n");
	system("$command");

}

#############################################
#                                           #
# Subroutine to delete files              r #
#                                           #
#############################################

sub delete_file
{

	my $file_to_be_deleted = "";	

	$file_to_be_deleted = $_[0];
	$command = "rm  $file_to_be_deleted";
	print("Deleting file   $file_to_be_deleted\n");
	system("$command");

}


#############################################
#                                           #
# Subroutine to execute unix command        #
#                                           #
#############################################

sub run_unix_command
{
	my $unix_command = "";
	my $step = "";	
	$unix_command = $_[0];
	$step = $_[1];
		
	print("$unix_command\n\n");

	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	print COMMAND_LOG "File: $loop_count/$no_of_files \tStep: $step/17        Time:  time();\n";
	print COMMAND_LOG "$unix_command\n";
	
	system("$unix_command");

}

##############################################
#                                            #
# Subroutine to record size of output file   #
#                                            #
##############################################

sub record_output_file_size
{
	my $outputfile = "";
	my $filesize = "";	

	$outputfile = $_[0];
	
	if (-e $outputfile)
	{
		$filesize = -s "$outputfile";
	}
	else
	{
		$filesize="Not found";
	}
	print COMMAND_LOG "\n  Output file: $outputfile\t\tSize: $filesize\n";
	print "\n  Output file: $outputfile\t\tSize: $filesize\n\n";
	
	####################################
	# Send e-mail if file size is zero #
	####################################
	if ($filesize < 1)
	{
		open(MAIL, "|/usr/sbin/sendmail -t");
		## Mail Header
		print MAIL "To: $email_address\n";
		print MAIL "From: $from\n";
		print MAIL "Subject: NGS ANALYSIS ZERO FILE SIZE: Sample $sample_name in run $run_title. File $outputfile\n\n";
		## Mail Body
		print MAIL "NGS_pipeline script version $version\n\n";
		print MAIL "The output file $outputfile has zero file size.\n\n";
		print MAIL "Something may be wrong.  Please check\n\n";
		close(MAIL);
	}

}

		
		
#############################################
#                                           #
# Subroutine to move log     			    #
#                                           #
#############################################

sub move_log
{

	my $file_to_be_moved = "";	

	$file_to_be_moved = $_[0];
	$command = "mv  $file_to_be_moved $results_folder/$file_to_be_moved";
	system("$command");

}

####################################################################
#                                                                  #
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test" from "test.fasta")                           #
#                                                                  #
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
	print "\nERR checking switched OFF\n\n";
	#$answer = <STDIN>;
	
	#open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
	
	#while($single_line = <LOG>)
	#{
	#	chomp $single_line;
		
	#	if (index(lc $single_line,"a user error") != -1)
	#	{
	#		print "\nA PROBLEM HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n";
	#		print "\n>> $single_line\n\n";
	#		print"\n\nContinue?\n\n";

	#		open (MAIL, "|/usr/sbin/sendmail -t");
	#		print MAIL "To: $email_address\n";
	#		print MAIL "From: $from\n";
	#		print MAIL "Subject: Problem in NGS pipeline version $version ($run_title)\n\n";
	#		print MAIL "A PROBLEM HAS OCCURED IN NGS ANALYSIS! ($run_title)\n\n$single_line\n\nPlease check the log.rtf for details.\n";

	#		close(MAIL);
			
	#		if ($exit_on_problem eq "yes"){print "\nExit on problem\n\n";}
	#	}

	#}
	
	#close(LOG);


}
	

sub test_mode_subroutine
{
		if ($testing_mode eq "on")
		{
			print "\n\n";
			print ">>>>>>>>>>>>>       Testing mode is set to ON         <<<<<<<<<<<<<<<<<\n";
			print ">>>>>>>>>>>>>  Press return to see a list of files    <<<<<<<<<<<<<<<<<\n\n";
			$answer = <STDIN>;
			$command = "ls -lh $sample_name*";
			system("$command");
			print "\n\n>>>>>>>>>>>>>  Press return to continue to the next stage    <<<<<<<<<<<<<<<<<\n\n";
			$answer = <STDIN>;
		}

}		