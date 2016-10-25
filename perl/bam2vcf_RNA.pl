#!/usr/bin/perl -w

##################################################################################
#									                                             #      
#	bam2vcf_RNA			                                                         #     
#									                                             #
#	This PERL script converts BAM files to VCF files, using the Broad pipeline	 #
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
use File::Basename ;
use Term::ANSIColor;
#use warnings::unused;

# VERSION OF SOFTWARE #
my $version						= "28";

####################
# Define variables #
####################


#Constants
my $gatk_directory				= "gatk"; # "gatk_v1" is the previous GATK version and "gatk" is the current GATK version 2.0
my $testing_mode				= "off";
my $GATK_validation_stringency	= "LENIENT";
my $picard_validation_stringency = "SILENT";
my $max_open_temp_files			= "1000"; # SET THIS BASED ON ULIMIT
my $max_file_handles			= "1000"; # SET THIS BASED ON ULIMIT
my $max_records_in_ram			= "600000";
my $delete_intermediate_files	= "yes";
my $validate_bam_files			= "yes";
my $check_depth					= "no";
my $ensembl_directory			= "/opt/ensembl.new";
my $keep_all_for_first_file		= "yes"; # This keeps all intermediate files for first file, for error-checking
my $debug_mode					= "true";
my $simple_pipeline				= "false";

#Various Parameters (e.g. for the Unified Genotyper)
my $stand_emit_conf				= 30;
my $stand_call_conf				= 30;
my $max_alt_alleles				= 6; # This is how many alleles the INDEL caller in UnifiedGenotyper can allow
my $indels_caller				= "UnifiedGenotyper"; # alternatives are UnifiedGenotyper OR SomaticIndelDetector


my $temp_dir_string				= " -Djava.io.tmpdir=javatempdir";
my $tempdir						= "javatempdir";

my $use_s_option_in_bwa_sampe	= "no"; #This runs bwa sampe with the '-s' option which disables Smith-Waterman alignment for the unmapped mate
my $remove_duplicates			= "yes";
my $use_samtools_rmdup			= "no"; # This means use samtools rmdup rather than picard/MarkDuplicates

my $fix_misencoded_qual_scores	= ""; # Yes or No This fixes old Illumina quality scores which are too high for GATK
my $misencoded_quals_fixed		= ""; # True or False.  Set to True when they are fixed by -fixMisencodedQuals
my $fix_quals_string			= ""; # This is the string for the command line option in GATK

my $ulimit						= 0; #Taken from the value given by "ulimit -n". Used to set max_file_handles and  max_open_temp_files
my $list_count					= 0; #Counter for the list of input files
my $no_of_files					= 0; #No of files in the list of input files (or no of lines if Paired Ends)
my $loop_count					= 0; #Loops round once for each of the files in the list of files
my $run_time					= 0;
my $results_folder_count		= 0; # Counts existing results folder to see if run has already partially run.
my $next_folder					= 0;
my $start_folder				= 1;
my $zero_error_count			= 0;
my $bam_missing_for_UG_count	= 0;
my $array_size					= 0; # Size of item array to check there are two columns in the file

my $ref_prefix					= "";
my $zero_size_reported			= "";
my $command						= "";
my $ref							= "";
my $mem							= "";
my $answer						= "";
my $data						= "PE"; # PE of SE
my $email_address				= "";
my $run_title					= "";
my $annotate_variants			= "";
my $maq							= "";
my $ref_genome					= "";
my $structural_variant_analysis	= "";
my $chromosome					= "";
my $insert						= "";
my $region						= "";
my $use_defined_region			= "";
my $stop_on_missing_bai			= "yes"; # Leave this as 'yes' as a starting value
my $start_time					= "";
my $end_time					= "";
my $current_time				= "";
my $last_current_time			= "";
my $stage_time					= "";
my $second_column_found			= "false"; # If second column is found in input file of fil e names (used as sample name)
my $single_line					= "";
my $bam_files_from_tophat		= "true"; # Bam files from TopHat need an extra ReorderSam stage

my $bam_file_mates_fixed		= "";
my $list_file					= "";
my $read_file_method			= "";
my $results_folder				= "";
my $input_string				= "";
my $running_input_string		= "";
my $species						= "";
my $ref_seq_name				= "";  # Name of reference sequence for pindel analysis
my $log_file					= "";
my $sample_name					= "";  # Name of each sample. Used to give the .vcf file an individual sample name
my $prefix						= "";  # Name of bam_file ommitting the suffix .bam
my $make_parallel_VCF			= ""; # This makes use of multiple BAM files in GATK to make a combined VCF output file
my $lib							= ""; # to be used as RGLB by AddOrReplaceReadGroups
my $GATK_region_string			= ""; # String to be used with the -L option of GATK tools such as UnifiedGenotyper
my $samtools_region_string		= ""; # String to be used with samtools
my $UG_region_string 			= ""; # String to be used with Unified Genotyper (just narrows down to chromosome e.g. -L chr15)
my $GATK_known_sites_string		= ""; # String for known SNPs for BaseRecalibrator
my $use_other_ref_sequence 		= "no";

#### File names ####
my $bam_file						= ""; # Original BAM file from fastq2bam
my $bai_file						= "";
my $bam_bai_file					= "";

my $rg_bam							= ""; # After AddReadGroups
my $rg_bai							= "";

my $runtitle_dedup_bam 				= ""; # After MarkDuplicates (dedup)
my $runtitle_dedup_bai 				= "";
my $runtitle_sample_metrics			= ""; # Output from MarkDuplicates

my $runtitle_dedup_bam_intervals 	= ""; # After RealignerTargetCreator (clean)

my $runtitle_clean_dedup_bam 		= ""; # After IndelRealigner (clean)
my $runtitle_clean_dedup_bai 		= "";

my $runtitle_sample_recal_csv 		= ""; # After CountCovariates

my $runtitle_clean_dedup_recal_bam 	= ""; # After TableRecalibration (recal)
my $runtitle_clean_dedup_recal_bai 	= "";

my $input_for_validate_bam			= ""; # Files ready to be validated (depends on whether simple or full)
my $input_for_validate_bai			= ""; # Files ready to be validated (depends on whether simple or full)

my $region_only_bam					= ""; # After reducing BAM to a region
my $region_only_bai					= "";

my $final_bam						= ""; # Final BAM file after pipeline has run
my $final_bai						= "";

my $readme_file						= "README_bam2vcf_RNA.out";

### File names new for RNA
my $rs_bam							= ""; # After ReorderSam
my $rs_bai							= "";


#VCF files for SNPs and Indels
my $snps_vcf					= "";
my $snps_vcf_idx				= "";
my $snps_parallel_vcf			= "";
my $snps_final_vcf				= "";
my $indels_vcf					= "";
my $indels_vcf_alt				= "";
my $indels_vcf_idx				= "";
my $indels_parallel_vcf			= "";
my $indels_final_vcf			= "";
my $indels_bed					= "";

#LOG files
my $validate_out				= "";
my $flagstat_out				= "";

#Other files
my $annotated_indels_parallel_vcf		= ""; # Version for multi-sample parallel VCf file
my $annotated_snps_parallel_vcf			= ""; # Version for multi-sample parallel VCf file
my $ref_dict					= ""; # This is the file like can fam 3.dict
my $insert_size_pdf				= ""; # Used by CollectInsertSizeMetrics
my $gc_bias_pdf					= ""; # Used by CollectGcBiasMetrics
my $command_log					= "";
my $bam_file_for_unified_genotyper = ""; # Used in parallel VCF section
my $dummy_dbsnp_file			= ""; # The dummy dbSNP file used by....
my $actual_dbsnp_file			= ""; # Correct dbSNP file for this chromosome (if it exists)
my $dbsnp_file					= ""; # dbsnp file used by GATK

my $runtitle_sample_recal_grp = ""; # For BaseRealigner

my @bam_file_array				= ();
my @sample_name_array			= ();
my @item						= ();

########################
# Define non variables #
########################

my $from= 'NGS_analysis@samba64.aht.org.uk';

			
			
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

##############
# Get ulimit #
##############

$ulimit = qx{echo `ulimit -n`};
print "ulimit 1: $ulimit\n\n";

#$ulimit = qx{ulimit -n};
#print "ulimit 2: $ulimit\n\n";

#$ulimit = system("ulimit -n");
#print "ulimit 3: $ulimit\n\n";


#################
#TURN LOGGER ON #
#################




print color 'bold cyan';

print "\n\n\n\n";

#$command =  "clear";
#system("$command");

if ($testing_mode eq "on"){print"\n\nTESTING MODE ON\n\n";}

print color 'reset';

use Term::ANSIColor;
print color 'bold magenta';

print "\n\n";
print "                 ##################################################\n";
print color 'bold white';
print "                                      bam2vcf_RNA                    \n";
print " \n";
print "                    Converts BAM files from TopHat into VCF files        \n";
print color 'reset';
print color 'bold magenta';
print "                 ##################################################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  This program is based on the Broad pipeline for processing BAM files.\n\n";
print "  This version is slightly modified for BAM files created using TopHat\n\n";
print "  (i.e. from RNA-sequencing data)\n\n";

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

$log_file = "$run_title"."_bam2vcf_log.out";

######################################
# E-mail address for notifications   #
######################################

print "\nPlease enter your email address:    ";
$email_address = <STDIN>;
chomp $email_address;

if ($email_address eq "m"){$email_address = 'mike.boursnell@aht.org.uk';}
if ($email_address eq "o"){$email_address = 'oliver.forman@aht.org.uk';}
if ($email_address eq "k"){$email_address = 'karen.steward@aht.org.uk';}


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
print "   Enter 4 for Human\n\n";

print "   Enter 5 for Strep. equi\n";
print "   Enter 6 for Strep. zoo\n\n";

print "   Enter 9 for other\n\n";


$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3/canfam3_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/canfam2/canfam2.fasta"; $ref_seq_name = "canfam2"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam2/canfam2_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "3" ){$ref = "/home/genetics/equcab2/equcab2.fasta"; $ref_seq_name = "equcab2"; $species = "equus_caballus";$dummy_dbsnp_file = "/home/genetics/equcab2/equcab2_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "4" ){$ref = "/home/genetics/human/human.fasta"; $ref_seq_name = "human"; $species = "homo_sapiens";$dummy_dbsnp_file = "/home/genetics/human/human_dummy_DBSNP.vcf"}

if (substr($answer,0,1) eq "5" ){$ref = "/home/genetics/strep_equi/strep_equi.fasta"; $ref_seq_name = "s_equi"; $species = "streptococcus_equi";$dummy_dbsnp_file = "/home/genetics/strep_equi/strep_equi_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "6" ){$ref = "/home/genetics/strep_zoo/strep_zoo.fasta"; $ref_seq_name = "s_zoo"; $species = "streptococcus_zoo";$dummy_dbsnp_file = "/home/genetics/strep_zoo/strep_zoo_dummy_DBSNP.vcf"}


if (substr($answer,0,1) eq "u" ){$ref = "/home/genetics/canfam3nu/canfam3nu.fasta"; $ref_seq_name = "canfam3nu"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3nu/canfam3nu_dummy_DBSNP.vcf"}


######################################
# Choose your own reference sequence #
######################################
if (substr($answer,0,1) eq "9" )
{
	$use_other_ref_sequence = "yes";
	
	until (-e "$ref")
	{
		print "What is the full path of your reference sequence?\n";
		print "e.g. /home/genetics/canfam3/canfam3.fasta\n\n";
		print " >   ";
		$ref = <STDIN>;
		chomp $ref;
		if (! -e $ref){print "\n\n>>>>>>    $ref not found!.  Try again   <<<<<<\n\n";}
	}
	
	$ref_prefix = &get_prefix($ref);
	
	$dummy_dbsnp_file = "$ref_prefix"."_dummy_DBSNP.vcf";
	
	print "Dummy DBSNP file: $dummy_dbsnp_file\n\n";
	
	print "\n\nChoose a short name for this reference (e.g. dog)\n\n";
	print " >  ";
	$ref_seq_name = <STDIN>;
	chomp $ref_seq_name;

}


$lib = $ref_seq_name;

$ref_dict="/home/genetics/$ref_seq_name/$ref_seq_name.dict";

###################################################################
# Check if REF.DICT files exist (if a pre-indexed file is chosen) #
###################################################################

if ((! -e "$ref") || (! -e "$ref.bwt") || (! -e "$ref.rbwt") || (! -e "$ref.rpac") || (! -e "$ref.pac") || (! -e "$ref.rsa") || (! -e "$ref.sa") || (! -e "$ref.amb") || (! -e "$ref.dict") || (! -e "$ref.ann") || (! -e "$ref.fai"))
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
	if(! -e "$ref.dict"){print "$ref.dict does not exist\n";}
	if(! -e "$ref.ann"){print "$ref.ann does not exist\n";}
	if(! -e "$ref.fai"){print "$ref.fai does not exist\n";}
	
	exit;
} 

###########################################################################
# Do you need to fix old-style Quality Scores (ie prior to Illumina 1.8)  #
# These will be over 41 and too high for the GATK programs                #
# the GATK option -fixMisendcodedQuals will correct these scores          #
###########################################################################


print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Are your Quality Scores of the correct Sanger type (i.e. new style)?   \n\n";
print "  Illumina scores prior to Illumina 1.8 were 'old-style'               \n";
print "  All current BAM files should be OK.                                  \n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";


print "   Enter 1 for YES - new style\n";
print "   Enter 2 for NO  - old style\n\n";

$answer = <STDIN>;
chomp $answer;

if ($answer eq ""){$answer = "1"} # default

if (substr($answer,0,1) eq "1"){$fix_misencoded_qual_scores = "no";$misencoded_quals_fixed = "true"}
if (substr($answer,0,1) eq "2"){$fix_misencoded_qual_scores = "yes";$misencoded_quals_fixed = "false"}


if ($fix_misencoded_qual_scores eq "yes"){$fix_quals_string = "-fixMisencodedQuals "} else {$fix_quals_string = ""}


###########################################################################
# Simple version of the pipeline                                                     #
# Would you like to focus the analysis on a specific region of the genome #
###########################################################################


print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Would you like to use the simple version of the pipeline?        \n";
print "(This omits the Indel Realigning and the Base Recalibrating)     \n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

print "   Enter 1 for YES\n";
print "   Enter 2 for NO\n\n";

$answer = <STDIN>;
chomp $answer;

if ($answer eq ""){$answer = "2"} # default

if (substr($answer,0,1) eq "1"){$simple_pipeline = "true"}
if (substr($answer,0,1) eq "2"){$simple_pipeline = "false"}


###########################################################################
# Region preferences                                                      #
# Would you like to focus the analysis on a specific region of the genome #
###########################################################################


#print "\n\n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
#print "Would you like to focus on alignments to specific region of the genome?\n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";


#print "   Enter 1 for YES\n";
#print "   Enter 2 for NO\n\n";

#$answer = <STDIN>;
#chomp $answer;

#if ($answer eq ""){$answer = "1"} # default

#if (substr($answer,0,1) eq "1"){$use_defined_region = "yes"}
#if (substr($answer,0,1) eq "2"){$use_defined_region = "no"}
$use_defined_region = "no";



##########################################################################################
# If no defined region is chosen we need to specify a SNP file for GATK BaseRecalibrator #
##########################################################################################
if ($use_defined_region eq "no")
{
	############################################
	# If the reference sequence is 'other' you # 
	# will have to specify a dummy SNP file    #
	############################################
	if ($use_other_ref_sequence eq "yes")
	{
		until (-e "$dummy_dbsnp_file")
		{
			print "If you don't choose a defined region, you still need a dummy SNP file.\n\n";
			
			print "Your reference sequence is $ref\n\n";

			print "What is the full path of your dummy dbsnp file?\n";
			print "e.g. /home/genetics/canfam3/canfam3_dummy_DBSNP.vcf\n\n";
			print " >   ";
			$dummy_dbsnp_file = <STDIN>;
			chomp $dummy_dbsnp_file;
			if (! -e $dummy_dbsnp_file)
			{
				print "\n\n";
				print "########################################################################\n";
				print ">>>>>>    $dummy_dbsnp_file not found!.  Try again   <<<<<<\n\n";
				print "The format of the file should be REFSEQ_dummy_DBSNP.vcf\n";
				print "where REFSEQ is the root name of your reference sequence.\n\n";
				
				print "For example if your reference sequence is canfam3.fasta,\n";
				print "the dummy SNP file is canfam3_dummy_DBSNP.vcf\n\n";
				print "########################################################################\n\n";
			}
		}
		$dbsnp_file = $dummy_dbsnp_file;
		$GATK_known_sites_string = " -knownSites $dbsnp_file";
	
	} #use_other_ref_sequence eq "yes"
	
	if ($use_other_ref_sequence eq "no")
	{
		$dbsnp_file = $dummy_dbsnp_file;
		$GATK_known_sites_string = " -knownSites $dbsnp_file";
	}	#use_other_ref_sequence eq "no"
	
} # if ($use_defined_region eq "no")


######################################################################################################
# If a defined region is chosen we look for a SNP file for this chromosome for GATK BaseRecalibrator #
######################################################################################################
if ($use_defined_region eq "yes")
{
	print "\nPlease define your region of interest (eg 'chr5:21000000-23000000' or just 'chr5'):      ";
	$region = <STDIN>;
	chomp $region;
	
	###################################################
	# If only chromosome  is given (e.g. 15 or chr15) #
	###################################################
	if (index($region,":") == -1)
	{
		if (index($region,"chr") == -1){$region = "chr"."$region";}
	}
	
	
	#################################################
	# If full region given chr15:34000000-390000000 #
	#################################################
	if (index($region,":") > -1)
	{
		if (index($region,"chr") > -1)
		{
			$chromosome = substr($region,3,index($region,":")-3);
		}
		if (index($region,"chr") == -1)
		{
			$region = "chr"."$chromosome";
		}
	}
	if (index($region,":") == -1)
	{
		if (index($region,"chr") == 0)
		{
			$chromosome = substr($region,3,99);
		}
	}
	
	if (index($chromosome,"chr") == -1){$chromosome = "chr"."$chromosome";}
	
	$GATK_region_string = " -L $region";
	$samtools_region_string = "$region";
	$UG_region_string = "-L $chromosome";
	$GATK_known_sites_string = "";
	
	
	####################################################
	# Check if there is a SNP file for this chromosome #
	####################################################
	if ($chromosome ne "")
	{
		$actual_dbsnp_file = "/home/genetics/$ref_seq_name/"."$chromosome"."snps.vcf";
		
		# File exists
		if (-e $actual_dbsnp_file)
		{
			print "\n\n";
			print "##############################################\n";
			print "A SNP file for this chromosome has been found.\n";
			print "SNP file $actual_dbsnp_file exists\n";
			print "##############################################\n\n";
			
			$dbsnp_file = $actual_dbsnp_file;
			$GATK_known_sites_string = " -knownSites $dbsnp_file";
		}
		
		# File not found
		if (!-e $actual_dbsnp_file)
		{
			print "\n\n";
			print "#######################################################\n";
			print "A SNP file for this chromosome has not been found.     \n";
			print "SNP file $actual_dbsnp_file doesn't exist              \n";
			print "Dummy SNP file $dummy_dbsnp_file will be used instead  \n";
			print "#######################################################\n\n";
			
			print "Do you want to stop and create a SNP file for this chromosome now?   (y/n)        ";
			
			$answer = <STDIN>;
			chomp $answer;
			if (lc $answer eq "y"){exit;} 
			
			##################################
			# Ask for name of dummy SNP file #
			##################################
			until (-e "$dummy_dbsnp_file")
			{
				print "Your reference sequence is $ref\n\n";
				
				print "What is the full path of your dummy dbsnp file?\n";
				print "e.g. /home/genetics/canfam3/canfam3_dummy_DBSNP.vcf\n\n";
				print " >   ";
				$dummy_dbsnp_file = <STDIN>;
				chomp $dummy_dbsnp_file;
				if (! -e $dummy_dbsnp_file){
					print "\n\n";
					print "########################################################################\n";
					print ">>>>>>    $dummy_dbsnp_file not found!.  Try again   <<<<<<\n\n";
					print "The format of the file should be REFSEQ_dummy_DBSNP.vcf\n";
					print "where REFSEQ is the root name of your reference sequence.\n\n";
					
					print "For example if your reference sequence is canfam3.fasta,\n";
					print "the dummy SNP file is canfam3_dummy_DBSNP.vcf\n\n";
					print "########################################################################\n\n";
				}
			}
			$dbsnp_file = $dummy_dbsnp_file;
			$GATK_known_sites_string = " -knownSites $dbsnp_file";
			
		} # actual DBSNP not found
	}
	
} # if ($use_defined_region eq "yes")

#########################################################
# Ask if you want to make a multiple parallel VCF file  #
#########################################################

###########################################################
# The parallel VCF option is now the only one recommended #
###########################################################

$make_parallel_VCF = "parallel";


######################################################################
# ASK IF YOU WANT TO READ THE FILENAMES FROM A "FILE OF FIL E NAMES" #
######################################################################
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
	if ($answer eq ""){$answer = "1"} # default
	
	$answer = substr($answer,0,1);
	
	if ($answer eq "1"){$read_file_method = "multiple"}
	if ($answer eq "2"){$read_file_method = "single"}

}


###################################
# Input fil e names               #
###################################

if ($read_file_method eq "single")
{
	until (-e $bam_file)
	{
		print "\nPlease input the name of your BAM file:      ";
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
	$sample_name = &get_prefix ($bam_file);
	$prefix = $sample_name;
	$second_column_found = "false";
	$sample_name_array[1] = $sample_name;
	
	#########################################################
	# Check sample_name doesn't have the string "bam" in it #
	#########################################################
	if (index($sample_name,"bam") > -1)
	{
		print "\n\n";
		print "############\n";
		print "# WARNING! #\n";
		print "############\n\n";
		
		print "File name prefix: $sample_name\n\n";
		
		print "File names cannot have the string 'bam' in the prefix\n\n";
		print "This will prevent samtools from working.\n\n";
		
		print "Please rename your files without the string 'bam' in the prefix\n\n";
		
		exit;
	}
		
} # End of read = single

if ($read_file_method eq "multiple")
{
	until (-e "$list_file")
	{
		print "\nPlease input the name of your file with a list of file names of the .bam files:    (type 'ls' to get a list of .txt files)  ";
		$list_file = <STDIN>;
		chomp $list_file;
		
		if (($list_file ne "ls" ) && (index($list_file,".txt") == -1 )){$list_file = $list_file.".txt"}

		
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
	
	while ($single_line = <LIST> ) 
	{
		chomp $single_line;
		
		@item=split(/\t/,$single_line);
			
		$array_size = scalar @item;
		
		if ($array_size == 1)
		{
			$bam_file=$single_line;
			$sample_name = &get_prefix ($bam_file);
		}
		
		if ($array_size == 2)
		{
			$bam_file = $item[0];
			$sample_name=$item[1];
			$second_column_found = "true";
		}
		
		#####################################################
		# Get prefix of bam file (ie omitting .bam suffix ) #
		#####################################################
		$prefix = &get_prefix($bam_file);
		
		
		######################################################
		# Add file type suffix .bam if user hasn't added it  #
		######################################################

		if (index($bam_file,".bam") == -1 ){$bam_file = $bam_file.".bam"}
	
		
		#########################################################
		# Check sample_name doesn't have the string "bam" in it #
		#########################################################
		if (index($sample_name,"bam") > -1)
		{
			print "\n\n";
			print "############\n";
			print "# WARNING! #\n";
			print "############\n\n";
			
			print "File name prefix: $sample_name\n\n";
			
			print "File names cannot have the string 'bam' in the prefix\n\n";
			print "This will prevent samtools from working.\n\n";
			
			print "Please rename your files without the string 'bam' in the prefix\n\n";
			
			exit;
		}
		
		$bam_file_array[$list_count]=$bam_file;
		$sample_name_array[$list_count]=$sample_name;
		
		$list_count = $list_count + 1;
	}

	close LIST;

	$no_of_files=$list_count - 1;
	
	
	###################
	# List file names #
	###################
	print "\n\nThere are $no_of_files BAM files in this file of file names.\n\n";
	
	if ($second_column_found eq "true")
	{
		print "There is also a second column in the input file which will be used for renaming the files.\n\n";
	}
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		if ($second_column_found eq "false") {print "File $list_count	\t$bam_file_array[$list_count]\n";}
		if ($second_column_found eq "true") {print "File $list_count	\t$bam_file_array[$list_count]\t$sample_name_array[$list_count]\n";}
	}

	print "\nIf these file names look OK, press enter to proceed (or 'Q' to quit):      ";
	$answer = <STDIN>;
	chomp $answer; 
	 
	if (lc $answer eq "q"){exit;} 
	
	
	##########################################################
	# If there isn't a BAM index file (.bai) then create one #
	##########################################################
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		$bam_file = $bam_file_array[$list_count];
		$prefix = &get_prefix ($bam_file);
		$bai_file= "$prefix".".bai";
		$bam_bai_file = "$prefix".".bam.bai";
	
		if (! -e "$bai_file"){print "$bai_file doesn't exist\n";}
		if (-e "$bai_file"){print "$bai_file exists\n";}
		if (! -e "$bam_bai_file"){print "$bam_bai_file doesn't exist\n";}
		if (-e "$bam_bai_file"){print "$bam_bai_file exists\n";}
		
		
		if ((! -e "$bai_file") && (! -e "$bam_bai_file"))
		{
			if ($stop_on_missing_bai eq "yes")
			{
				print "######################################\n";
				print "#  No index file for BAM file found! #\n";
				print "######################################\n\n";
				
				print COMMAND_LOG "\nNo index file for $bam_file found\n\n";
				
				print "Each BAM file should have a matching BAI file.\n\n";
				print "If you continue the BAI files will be created, but using the original ones is better, if you have them.\n\n";
				
				print "Do you want to stop now, and put the BAI files in this same directory? (y/n)  ";
				
				$answer = <STDIN>;
				if (lc $answer eq "y"){exit;} 
				if (lc $answer eq "n"){$stop_on_missing_bai = "no"} 
			}
				
		}
	
	} # check if BAI files are present
	
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
		#print "File $list_count	\t$bam_file_array[$list_count]\n";
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

print color 'bold magenta';
print "\n  NOTE: This is currently done from NGS SNP Handler\n\n";
print color 'reset';

print "       Press return to continue\n\n";

#print "   Enter 1 for YES\n";
#print "   Enter 2 for NO\n\n";

$answer = <STDIN>;
chomp $answer;


#if (substr($answer,0,1) eq "1"){$annotate_variants = "yes"}
#if (substr($answer,0,1) eq "2"){$annotate_variants = "no"}
#if ($answer eq ""){$annotate_variants = "yes"}

$annotate_variants = "no";



###############################
# PINDEL analysis preferences #
###############################

if ($data eq "PExxxx")
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
$structural_variant_analysis="no";
	
$command =  "clear";
system("$command");

$| = 1;

open(STDOUT, "| tee $log_file");


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
$answer = <STDIN>;
chomp $answer; 
 
if (lc $answer eq "q"){exit;} 


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

if ($remove_duplicates eq "yes"){print "Remove duplicates\t\t\tYES\n\n";}
if ($remove_duplicates eq "no") {print "Do not remove duplicates\t\t\tNO\n\n";}

if ($annotate_variants eq "yes")        {print "Annotate variants     \t\t\tYES\n\n";}
if ($annotate_variants eq "no")         {print "Annotate variants     \t\t\tNO\n\n";}

if ($structural_variant_analysis eq "no")          {print "Structural variant analysis\t\tNO\n\n";}
	
if ($structural_variant_analysis eq "yes")
{
	{print "Structural variant analysis\t\tYES\n\n";}

	if ($ref_genome eq "partial"){print "(Chromosome $chromosome)\n\n"}
	if ($ref_genome eq "whole"){print "\n\n"}
}

if ($use_defined_region eq "yes")
{
	print "Defined region:       \t\t\t$region\n\n";
	print "  GATK region string:   \t\t\t$GATK_region_string\n";
	print "  UG region string:     \t\t\t$UG_region_string\n";
	print "  Samtools region string:\t\t\t$samtools_region_string\n\n";
}
if ($use_defined_region eq "no")
{
	print "Defined region:\t\t\t\tNO\n\n";
}

print "SNP file for BaseRecalibrator:\t\t$dbsnp_file\n\n";

print "GATK knownsites string:       \t\t$GATK_known_sites_string\n\n";

#print "bwa settings:  use -s in bwa sampe:\t$use_s_option_in_bwa_sampe\n\n";

print "Memory setting: \t\t\t$mem\n\n";

print "UnifiedGenotyper constants:\n\n";

print "    --stand_emit_conf:\t\t\t$stand_emit_conf\n";
print "    --stand_call_conf:\t\t\t$stand_call_conf\n";
print "    --max_alt_alleles:\t\t\t$max_alt_alleles (for calling Indels)\n\n";


print color "reset";

print "Please press enter to proceed (or 'Q' to quit):      ";
$answer = <STDIN>;
chomp $answer; 
 
if (lc $answer eq "q"){exit;} 

print "\n\n";

print color "bold white";

print "==================================================\n";
print "SUMMARY 3 - PLEASE CHECK THESE DETAILS CAREFULLY!!\n";
print "==================================================\n\n";

print color "bold green";
print "DATA FILES :\n\n";
print color "bold white";

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

$answer = <STDIN>;
chomp $answer; 
 
if (lc $answer eq "q")
{
	exit;
} 
		  

$start_time = time();

print color "reset";

#########################
# open Command Log file #
#########################
$command_log = "$run_title"."_bam2vcf_command_log.out";

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "COMMAND LOG for BAM2VCF_RNA version $version\n\n";
print COMMAND_LOG "Analysis name (run title): $run_title\n\n";

print COMMAND_LOG "PREFERENCES: \n\n";

print COMMAND_LOG "Reference sequence:      \t$ref\n";
print COMMAND_LOG "Name of this analysis:   \t$run_title\n"; 
print COMMAND_LOG "E-mail address:          \t$email_address\n";
print COMMAND_LOG "Data type (PE or SE)     \t$data\n";
print COMMAND_LOG "Check depth of coverage: \t$check_depth\n";
print COMMAND_LOG "Remove duplicates:       \t$remove_duplicates\n";
print COMMAND_LOG "Annotate variants:       \t$annotate_variants\n";
print COMMAND_LOG "Use defined region:      \t$use_defined_region\n";
print COMMAND_LOG "Region:                  \t$region\n";
print COMMAND_LOG "Chromosome:              \t$chromosome\n";
print COMMAND_LOG "Fix quality scores:      \t$fix_misencoded_qual_scores\n";
print COMMAND_LOG "Simple pipeline:         \t$simple_pipeline\n";
print COMMAND_LOG "Memory setting:          \t$mem\n\n";
print COMMAND_LOG "UnifiedGenotyper constants:\n\n";
print COMMAND_LOG "    -stand_emit_conf:\t\t\t$stand_emit_conf\n";
print COMMAND_LOG "    -stand_call_conf:\t\t\t$stand_call_conf\n\n";
print COMMAND_LOG "Make parallel VCF:		\t$make_parallel_VCF\n";
print COMMAND_LOG "Method of file name entry is $read_file_method\n\n";

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

$last_current_time = time();

for ($loop_count=$start_folder;$loop_count <=$no_of_files;$loop_count++)
{

	##########################################################################
	# Fixing misencoded quality scores must be set at the start of each loop #
	##########################################################################
	if ($fix_misencoded_qual_scores eq "yes"){$fix_quals_string = "-fixMisencodedQuals "} else {$fix_quals_string = ""}
	

	################################################################
	# Set flag to "no" for if any zero file size has been reported #
	# (we only want it for the first time in each loop)            #
	################################################################
	$zero_size_reported = "no";
	
	
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
	$prefix = &get_prefix ($bam_file);
	
	######################################################
	# Create sample name (to use for name VCF files etc) #
	######################################################
	$sample_name = $sample_name_array[$loop_count];
	
	
	#########################################################
	# Set up various File names (not all necessarily used)  #
	#########################################################
	
	$bai_file= "$prefix".".bai";
	$bam_bai_file = "$prefix".".bam.bai";
	
	# After AddReadGroups
	$rg_bam  = "$run_title"."_"."$sample_name"."_rg.bam";
	$rg_bai  = "$run_title"."_"."$sample_name"."_rg.bai";
	
	# After ReorderSam (new for RNA)
	$rs_bam  = "$run_title"."_"."$sample_name"."_rs.bam";
	$rs_bai  = "$run_title"."_"."$sample_name"."_rs.bai";
	
	# After MarkDuplicates (deduped)
	$runtitle_dedup_bam = "$run_title"."_"."$sample_name".".dedup.bam";
	$runtitle_dedup_bai = "$run_title"."_"."$sample_name".".dedup.bai";
	
	# After RealignerTargetCreator
	$runtitle_dedup_bam_intervals = "$run_title"."_"."$sample_name".".dedup.bam.intervals";
	
	# After IndelRealigner (cleaned)
	$runtitle_clean_dedup_bam = "$run_title"."_"."$sample_name".".clean.dedup.bam";
	$runtitle_clean_dedup_bai = "$run_title"."_"."$sample_name".".clean.dedup.bai";
	
	# After CountCovariates
	$runtitle_sample_recal_csv  = "$run_title"."_"."$sample_name".".recal.csv";
	$runtitle_sample_recal_grp  = "$run_title"."_"."$sample_name".".recal.grp"; # Use for BaseRecalibrator
	
	# After TableRecalibration (recal)
	$runtitle_clean_dedup_recal_bam =  "$run_title"."_"."$sample_name".".clean.dedup.recal.bam";
	$runtitle_clean_dedup_recal_bai =  "$run_title"."_"."$sample_name".".clean.dedup.recal.bai";
	
	# After reducing BAM to a region
	$region_only_bam = "$run_title"."_"."$sample_name"."_region.bam";
	$region_only_bai = "$run_title"."_"."$sample_name"."_region.bai";
	
	# After all stages have finsished, the BAM file is renamed to final_bam
	$final_bam = "$run_title"."_"."$sample_name"."_final.bam";
	$final_bai = "$run_title"."_"."$sample_name"."_final.bai";
	
	
	#SNP VCF files
	$snps_vcf = "$run_title"."_"."$sample_name"."_snps.vcf";
	$snps_vcf_idx = "$run_title"."_"."$sample_name"."_snps.vcf.idx";
	$snps_final_vcf = "SNPS_"."$run_title"."_"."$sample_name".".vcf";
	$snps_parallel_vcf = "SNPS_"."$run_title".".vcf";
	$annotated_snps_parallel_vcf = "annotated_SNPS_"."$run_title".".vcf";
	
	#Indel VCF files
	$indels_vcf = "$run_title"."_"."$sample_name"."_indels.vcf";
	$indels_vcf_alt = "$run_title"."_"."$sample_name"."_indels_alt.vcf";
	$indels_vcf_idx = "$run_title"."_"."$sample_name"."_indels.vcf.idx";
	$indels_bed = "$run_title"."_"."$sample_name"."_indels.bed";
	$indels_final_vcf = "INDELS_"."$run_title"."_"."$sample_name".".vcf";
	$indels_parallel_vcf = "INDELS_"."$run_title".".vcf";
	$annotated_indels_parallel_vcf = "annotated_INDELS_"."$run_title".".vcf";
	
	$bam_file_mates_fixed = "$sample_name"."_mates_fixed.bam";
	$insert_size_pdf = "$sample_name"."_insert_size.pdf";
	$gc_bias_pdf = "$sample_name"."_gc_bias.pdf";
	$validate_out = "$run_title"."_"."$sample_name"."_validate.out";
	$flagstat_out = "$run_title"."_"."$sample_name"."_flagstat.out";
	
	#file names
	
	$runtitle_sample_metrics = "$run_title"."_"."$sample_name".".metrics";
	
	
	

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
		
		&record_output_file_size ("$bai_file");	
		
		print "\n===============================================";
		print "\nFile $loop_count/$no_of_files   New Index for BAM file created";
		print "\n===============================================\n\n\n";
			
	}
	
	

##########################################################################################
# bam2vcf_RNA stage 1: Use picard/AddOrReplaceReadGroups to update ReadGroup info in BAM file  #
##########################################################################################

print	"\n============================================================================================================";
print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 1/13 Using picard/AddOrReplaceReadGroups to update ReadGroup info in BAM file";
print 	"\n============================================================================================================\n\n\n";

&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/AddOrReplaceReadGroups.jar I=$bam_file O=$rg_bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true RGID=$sample_name  RGLB=$lib  RGPL=illumina  RGPU=$sample_name  RGSM=$sample_name","1");

&record_output_file_size ("$rg_bam");

print	"\n=============================================================================================================";
print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 1/13 COMPLETED picard/AddOrReplaceReadGroups updated ReadGroup info in BAM file";
print 	"\n=============================================================================================================\n\n\n";

print "\n* * * \n";


##########################################################################################
# bam2vcf_RNA stage 2: Use picard/ReorderSam to make sure contigs are in the right order    #
# (only required in the BAM files have come from TopHat)                                 #
##########################################################################################

if ($bam_files_from_tophat eq "true")
{
	print	"\n=============================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 2/13 Using picard/ReorderSam to update contig order in BAM file";
	print 	"\n=============================================================================================\n\n\n";

	&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/ReorderSam.jar I=$rg_bam O=$rs_bam R=$ref VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true","2");

	&record_output_file_size ("$rs_bam");

	print	"\n=================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 2/13 COMPLETED picard/ReorderSam to update contig order in BAM file";
	print 	"\n=================================================================================================\n\n\n";

	print "\n* * * \n";

} # bam file is from TopHat


#########################################
# If it is not from TopHat then ...
#########################################
if ($bam_files_from_tophat eq "false")
{
	$rs_bam = $rg_bam;
	&record_output_file_size ("$rs_bam");
}


##########################################################################################################
# bam2vcf_RNA stage 3: Use picard/MarkDuplicates to mark duplicates                                           #
##########################################################################################################

if ($remove_duplicates eq "yes")
{

	if ($use_samtools_rmdup eq "no")
	{
		print	"\n\n==============================================================================";
		print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 3/13 Use picard/MarkDuplicates to mark duplicates";
		print 	"\n==============================================================================\n\n\n";

		&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/MarkDuplicates.jar I=$rs_bam O=$runtitle_dedup_bam  VALIDATION_STRINGENCY=$picard_validation_stringency CREATE_INDEX=true M=$runtitle_sample_metrics MAX_FILE_HANDLES=$max_file_handles MAX_RECORDS_IN_RAM=$max_records_in_ram","3");

		&record_output_file_size ("$runtitle_dedup_bam");
		
		print	"\n==========================================================================================";
		print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 3/13 COMPLETED using picard/MarkDuplicates to mark duplicates";
		print 	"\n==========================================================================================\n\n";
		
		print "\n* * * \n";
		
	} # don't use samtools rmdup - use MarkDuplicates
	
	if ($use_samtools_rmdup eq "yes")
	{
		print	"\n\n==========================================================================";
		print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 3/13 Use samtools rmdup to remove duplicates";
		print 	"\n==========================================================================\n\n\n";

		&run_unix_command (" /opt/samtools/samtools rmdup $rs_bam $runtitle_dedup_bam","3 samtools");

		&record_output_file_size ("$runtitle_dedup_bam");
		
		print	"\n======================================================================================";
		print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 3/13 COMPLETED using samtools rmdup to remove duplicates";
		print 	"\n======================================================================================\n\n";
		
		print "\n* * * \n";
		
	} # do use samtools rather than picard to remove duplicates

}


##############################################################
# If you are not removing duplicates you have to rename the  #
# BAM file so it can be processed correctly by the next step #
##############################################################
if ($remove_duplicates eq "no")
{
	$runtitle_dedup_bam = $rs_bam;
	&record_output_file_size ("$runtitle_dedup_bam");
}




#################################################################################
# Only carry out these next 4 stages if you have not chosen the simple pipeline #
#################################################################################
if ($simple_pipeline eq "false")
{
	#################################################################################################################
	# bam2vcf_RNA stage 4: Use GenomeAnalysisTK/RealignerTargetCreator to identify intervals in need of realignment #
	#################################################################################################################

	print	"\n\n==========================================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 4/13 Use GenomeAnalysisTK/RealignerTargetCreator to identify intervals in need of realignment";
	print 	"\n==========================================================================================================================\n\n\n";

	&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T RealignerTargetCreator -I $runtitle_dedup_bam -o $runtitle_dedup_bam_intervals $GATK_region_string -R $ref -mismatch 0.0 -S $GATK_validation_stringency $fix_quals_string","4");

	&record_output_file_size ("$runtitle_dedup_bam_intervals");

	print	"\n===================================================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 4/13 COMPLETED GenomeAnalysisTK/RealignerTargetCreator has identified intervals in need of realignment";
	print 	"\n===================================================================================================================================\n\n";

	print "\n* * * \n";
	
	
	##########################################################################################################
	# bam2vcf_RNA stage 5: Use GenomeAnalysisTK/IndelRealigner to Realign Indels                             #
	##########################################################################################################

	print	"\n\n=========================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 5/13 Using GenomeAnalysisTK/IndelRealigner to Realign Indels";
	print 	"\n=========================================================================================\n\n\n";

	&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T IndelRealigner -I $runtitle_dedup_bam  -o $runtitle_clean_dedup_bam $GATK_region_string -R $ref -targetIntervals $runtitle_dedup_bam_intervals -compress 0 -model USE_READS -S $GATK_validation_stringency $fix_quals_string","5");

	&record_output_file_size ("$runtitle_clean_dedup_bam");

	print	"\n===================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 5/13 COMPLETED using GenomeAnalysisTK/IndelRealigner to Realign Indels";
	print 	"\n===================================================================================================\n\n";

	print "\n* * * \n";

	
	#####################################################################################
	# Record that the misencoded base qualities have been fixed (if they required it)   #
	# (The GATK step inputs a BAM file and outputs a BAM file with the fixed qualities) #
	#####################################################################################
	$misencoded_quals_fixed = "true";
	$fix_quals_string = "";
	
	
	##########################################################################################################
	# bam2vcf_RNA stage 6: Use GenomeAnalysisTK/BaseRecalibrator to generate a recalibration table           #
	##########################################################################################################

	print	"\n\n===========================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 6/13 Using GenomeAnalysisTK/BaseRecalibrator to generate a Recalibration Table ";
	print 	"\n===========================================================================================================\n\n\n";


	&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T BaseRecalibrator -I $runtitle_clean_dedup_bam -o $runtitle_sample_recal_grp $GATK_region_string -R $ref $GATK_known_sites_string $fix_quals_string","6");

	&record_output_file_size ("$runtitle_sample_recal_grp");

	print	"\n===========================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 6/13 Used GenomeAnalysisTK/BaseRecalibrator to generate a Recalibration Table ";
	print 	"\n===========================================================================================================\n\n";


	print "\n* * * \n";
	
	
	##########################################################################################################
	# bam2vcf_RNA stage 7: Use GenomeAnalysisTK/PrintReads to update the base quality scores                 #
	##########################################################################################################

	print	"\n\n=====================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 7/13 Using GenomeAnalysisTK/PrintReads to update the base quality scores ";
	print 	"\n=====================================================================================================\n\n\n";


	&run_unix_command ("java  $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T PrintReads  -I $runtitle_clean_dedup_bam -o $runtitle_clean_dedup_recal_bam $GATK_region_string -R $ref  -BQSR $runtitle_sample_recal_grp -S $GATK_validation_stringency $fix_quals_string","7");

	&record_output_file_size ("$runtitle_clean_dedup_recal_bam");

	print	"\n====================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 7/13 Used GenomeAnalysisTK/PrintReads to update the base quality scores ";
	print 	"\n====================================================================================================\n\n";

	print "\n* * * \n";
	
} # if $simple_pipeline eq "false"


##############################################################################
# If misencoded quality scores have been fixed, change the string for GATK   #
# Note: they are only fixed in the full pipeline because the simple pipeline #
# has not yet had any GATK steps which write a new BAM file                  #
##############################################################################
if ($fix_misencoded_qual_scores eq "yes")
{
	if ($misencoded_quals_fixed eq "false"){$fix_quals_string = "-fixMisencodedQuals "} else {$fix_quals_string = ""}
}



####################################################################
# Whether you have chosen the simple or the full pipeline          #
# will affect which file you take into the next stage (validating) #
####################################################################
if ($simple_pipeline eq "true")
{
	$input_for_validate_bam = $runtitle_dedup_bam;
	$input_for_validate_bai = $runtitle_dedup_bai;
}
if ($simple_pipeline eq "false")
{
	$input_for_validate_bam = $runtitle_clean_dedup_recal_bam;
	$input_for_validate_bai = $runtitle_clean_dedup_recal_bai;
}

	
	
##########################################################################################################
# bam2vcf_RNA stage 8: Use picard/ValidateSamFile to validate the cleaned BAM file (Summary only)             #
##########################################################################################################

if ($validate_bam_files eq "yes")
{

	print	"\n===============================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 8/13 Using picard/ValidateSamFile to validate the cleaned BAM file";
	print 	"\n===============================================================================================\n\n\n";

	
	&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/ValidateSamFile.jar INPUT=$input_for_validate_bam OUTPUT=$validate_out VALIDATION_STRINGENCY=$picard_validation_stringency MODE=SUMMARY MAX_OUTPUT=100 MAX_OPEN_TEMP_FILES=$max_open_temp_files","8");

	#&record_output_file_size ("$validate_out"); # not required for flow of ipeline. Seems to get zero size when it isn't
	
	&move_to_results_folder ("$validate_out");
	
	print	"\n=========================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 8/13 COMPLETED using picard/ValidateSamFile to validate the cleaned BAM file";
	print 	"\n=========================================================================================================\n\n\n";

	print "\n* * * \n";
	

#########################################################################################################
# bam2vcf_RNA stage 9: Use samtools flagstat to get simple stats on bam file                                 #
#########################################################################################################


	print "\n=====================================================================";
	print "\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 9/13 samtools flagstat to be carried out";
	print "\n=====================================================================\n\n";


	&run_unix_command("/opt/samtools/samtools flagstat $input_for_validate_bam > $flagstat_out","9");

	&record_output_file_size("$flagstat_out");

	&move_to_results_folder ("$flagstat_out");
		
	print "\n=========================================================================";
	print "\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 9/13 COMPLETED samtools flagstat carried out";
	print "\n=========================================================================\n\n";
	
	print "\n* * * \n";
	
} # if validate_bam_files



######################################################################################################################
# bam2vcf_RNA stage 10: If you are using a defined region then create a new smaller BAM file with only this region in it   #
######################################################################################################################
if ($use_defined_region eq "yes")
{

	print "\n\n=========================================================";
	print "\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 10a/13:  Making Bam file for region";
	print "\n    Region: $region";
	print "\n===========================================================\n\n";
	
	&run_unix_command("/opt/samtools/samtools view $runtitle_clean_dedup_recal_bam $region -b -o $region_only_bam","9a");

	&record_output_file_size ("$region_only_bam");	
	
	print "\n============================================================";
	print "\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 10a/13:  Bam file for region created";
	print "\n    Region: $region";
	print "\n============================================================\n\n\n";
	
	print "\n* * * \n";
	
	
	########################################################
	# Now make an index file for this new smaller BAM file #
	########################################################
	
	print "\n\n============================================================";
	print "\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 10b/13:  New Index for region-only BAM file being created";
	print "\n============================================================\n\n";


	&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/BuildBamIndex.jar I=$region_only_bam O=$region_only_bai VALIDATION_STRINGENCY=LENIENT","9b");

	&record_output_file_size ("$region_only_bai");
	
	print "\n============================================================";
	print "\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 10b/13:  Bam Index file for region $region created";
	print "\n============================================================\n\n\n";
	
	print "\n* * * \n";
	
} # if ($use_defined_region eq "yes"}



#########################################################
# Whether you have chosen a defined region              #
# will affect which file you name as the final bam file #
#########################################################
if ($use_defined_region eq "yes")
{
	&run_unix_command("mv $region_only_bam $final_bam","Change name to final_bam");
	&run_unix_command("mv $region_only_bai $final_bai","Change name to final_bai");
}
if ($use_defined_region eq "no")
{
	&run_unix_command("mv $input_for_validate_bam $final_bam","Change name to final_bam");
	&run_unix_command("mv $input_for_validate_bai $final_bai","Change name to final_bai");
}

	


###########################################################################
# Create input string for parallel use of the UnifiedGenotyper at the end #
###########################################################################
$running_input_string = $running_input_string." -I $results_folder/$final_bam";


########################################################################################
# bam2vcf_RNA stage 11: Making SNP calls   
########################################################################################


if ($make_parallel_VCF eq "consecutive" || $make_parallel_VCF eq "parallel and consecutive")
{

	print "\n======================================================";
	print "\nFile $loop_count/$no_of_files   Step 11/13 Making SNP calls";
	print "\n======================================================\n\n";

	&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $UG_region_string -glm SNP -I $final_bam -stand_call_conf $stand_call_conf -stand_emit_conf $stand_emit_conf -o $snps_vcf -S $GATK_validation_stringency $fix_quals_string","11");

	&record_output_file_size ("$snps_vcf");
	
	
	print "\n==========================================================";
	print "\nFile $loop_count/$no_of_files   Step 11/13 COMPLETED: SNP calls have been made";
	print "\n==========================================================\n";

	print "\n* * * \n";
	
} # End of if ($make_parallel_VCF eq "consecutive" || $make_parallel_VCF eq "parallel and consecutive")



########################################################################################
# bam2vcf_RNA stage 12:  Making INDEL calls   
########################################################################################

if ($make_parallel_VCF eq "consecutive" || $make_parallel_VCF eq "parallel and consecutive")
{

	print "\n==================================================";
	print "\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 12/13 Making INDEL calls";
	print "\n==================================================\n\n";
	
	if ($indels_caller eq "SomaticIndelDetector")
	{
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T SomaticIndelDetector -R $ref -I $final_bam -o $indels_vcf  -ws 300 --unpaired $fix_quals_string","12");
	}
	
	if ($indels_caller eq "UnifiedGenotyper")
	{
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T UnifiedGenotyper $UG_region_string -glm INDEL --max_alternate_alleles $max_alt_alleles -R $ref -I $final_bam -o $indels_vcf  -S $GATK_validation_stringency $fix_quals_string","12");
	}
	
	&record_output_file_size ("$indels_vcf");	
		
		
	print "\n===========================================================";
	print "\nFile $loop_count/$no_of_files   bam2vcf_RNA stage 10/13 COMPLETED making INDEL calls";
	print "\n===========================================================\n\n\n";


	print "\n* * * \n";
	
	
	##############################################
	# Now run variant effect predictor if chosen #
	##############################################


} # end of if ($make_parallel_VCF eq "consecutive" || $make_parallel_VCF eq "parallel and consecutive")





####################################################################
# Checking depth of coverage, average insert size and GC bias
####################################################################

if ($check_depth eq "yes")
{
	print "\n Checking depth of coverage....\n\n";

	&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T DepthOfCoverage -R $ref -I $final_bam -o depth  -omitBaseOutput -S $GATK_validation_stringency $fix_quals_string","depth");

}

print "\n====================================";
print "\nFile $loop_count/$no_of_files   Calculating GC Bias....";
print "\n====================================\n\n";

&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/CollectGcBiasMetrics.jar REFERENCE_SEQUENCE=$ref INPUT=$final_bam O=out1.junk CHART_OUTPUT=$gc_bias_pdf VALIDATION_STRINGENCY=LENIENT","GC_BIAS");

&record_output_file_size("$gc_bias_pdf");

print "\n========================================";
print "\nFile $loop_count/$no_of_files   GC bias has been calculated";
print "\n========================================\n";

print "\n* * * \n";


##############################
# Plot Insert Size Histogram #
##############################
print "\n======================================================";
print "\n Plotting insert size histogram (for PE data only)...";
print "\n======================================================\n\n";


&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/CollectInsertSizeMetrics.jar INPUT=$final_bam O=out2.junk HISTOGRAM_FILE=$insert_size_pdf VALIDATION_STRINGENCY=LENIENT","SIZE_METRICS");

&record_output_file_size_no_mail("$insert_size_pdf");



print "\n=========================================";
print "\nFile $loop_count/$no_of_files   Insert size histogram plotted";
print "\n=========================================\n\n\n";

print "\n* * * \n";

	##########################################
	# Move various files to a Results Folder #
	##########################################
	
	&move_to_results_folder ("$final_bam");
	&move_to_results_folder ("$final_bai");
	&move_to_results_folder ("$runtitle_sample_metrics"); # MarkDuplicates info
	
	if ($make_parallel_VCF eq "consecutive" || $make_parallel_VCF eq "parallel and consecutive")
	{
		&move_to_results_folder ("$snps_vcf");
		&move_to_results_folder ("$indels_vcf");
	}
	
	if ($check_depth eq "yes")
	{
		&move_to_results_folder ("depth");
		&move_to_results_folder ("depth.sample_summary");
	}
	
	&move_to_results_folder ("$gc_bias_pdf");
	&move_to_results_folder ("$insert_size_pdf");
	
	##############################################################
	# For first loop copy over all files (for checking purposes) #
	# (only if $keep_all_for_first_file = "yes")                 #
	##############################################################
	if (($loop_count == 1) && ($keep_all_for_first_file eq "yes"))
	{
		&move_to_results_folder ("$rg_bam");
		&move_to_results_folder ("$rg_bai");
		&move_to_results_folder ("$rs_bam");
		&move_to_results_folder ("$rs_bai");
		&move_to_results_folder ("$runtitle_dedup_bam");
		&move_to_results_folder ("$runtitle_dedup_bai");
		&move_to_results_folder ("$runtitle_clean_dedup_bam");
		&move_to_results_folder ("$runtitle_clean_dedup_bai");
		&move_to_results_folder ("$runtitle_clean_dedup_recal_bam");
		&move_to_results_folder ("$runtitle_clean_dedup_recal_bai");
		&move_to_results_folder ("$region_only_bam");
		&move_to_results_folder ("$region_only_bai");
		&move_to_results_folder ("$runtitle_sample_recal_csv");
	}
	
	########################################
	# Some files don't need to be kept     #
	# (keep files for first loop to check) #
	########################################
	
	if ($delete_intermediate_files eq "yes")
	{
		# Do not delete files on first loop (for debugging)
		if ($loop_count > 1)
		{
			#Simple pipeline only
			&delete_file("$rs_bam");
			&delete_file("$rs_bai");
			&delete_file("$runtitle_dedup_bam");
			&delete_file("$runtitle_dedup_bai");
			
			if ($use_defined_region eq "yes")
			{
				&delete_file("$region_only_bam");
				&delete_file("$region_only_bai");
			}
			
			#Full pipeline (some extras to delete)
			if ($simple_pipeline eq "false")
			{
				&delete_file("$runtitle_dedup_bam_intervals");
				&delete_file("$runtitle_clean_dedup_bam");
				&delete_file("$runtitle_clean_dedup_bai");
				&delete_file("$runtitle_sample_recal_grp");
				&delete_file("$runtitle_clean_dedup_recal_bam");
				&delete_file("$runtitle_clean_dedup_recal_bai");		
			}
			
		}
	}
	

	#####################################
	# Change the names of various files #
	#####################################

	$command = "mv $results_folder/$indels_vcf $results_folder/$indels_final_vcf";
	system("$command");

	$command = "mv $results_folder/$snps_vcf $results_folder/$snps_final_vcf";
	system("$command");

	$command = "mv $results_folder/depth $results_folder/depth_bp.xls";
	system("$command");

	$command = "mv $results_folder/depth.sample_summary $results_folder/depth_summary.xls";
	system("$command");



	
	##########
	# README #
	##########

	open (READMEFILE, ">$readme_file"); 

	print	READMEFILE	"============================================\n";
	print	READMEFILE	"Summary of bam2vcf_RNA results files\n";
	print	READMEFILE	"============================================\n\n\n";


	print	READMEFILE	"ALIGNMENT FILES\n\n";

	print	READMEFILE	"\t$final_bam	   \tBest alignments to the reference after processing by GATK\n\n\n";

	
	if ($annotate_variants eq "yes")
	{
		print 	READMEFILE	"ANNOTATION FILES\n\n";
		
		print	READMEFILE	"\tannotated_indels.xls   \tList of Indels annotated using the ensembl database\n";
		print	READMEFILE	"\tannotated_snp.xls      \tList of SNPs annotated using the ensembl database\n\n\n";
	}
	
	print	READMEFILE	"PDF FILES\n\n";

	print	READMEFILE	"\t$gc_bias_pdf          \tHistogram of GC content of aligned reads\n";
	print	READMEFILE	"\t$insert_size_pdf      \tInsert size histogram\n\n\n";
	
	if ($structural_variant_analysis eq "yes")
	{
		print 	READMEFILE	"STRUCTURAL VARIANT FILES\n\n";
		print 	READMEFILE	"\tSV_Break_points.txt\n";
		print 	READMEFILE	"\tSV_Deletions.txt\n";
		print 	READMEFILE	"\tSV_Inversions.txt\n";
		print 	READMEFILE	"\tSV_Long_insertions.txt\n";
		print 	READMEFILE	"\tSV_Non_template_seq_in_deleletion.txt\n";
		print 	READMEFILE	"\tSV_Tandom_Dup.txt\n\n\n";
	}
	print	READMEFILE	"INFORMATION FILES\n\n";	

	if ($check_depth eq "yes") {print	READMEFILE	"Depth_summary	\tSummary of reads depth across the target region\n\n\n";}
	
	print	READMEFILE	"LOG FILES\n\n";
	
	print	READMEFILE	"\t$command_log  \tCommand log for bam2vcf_RNA (very useful for tracking errors)\n";
	print	READMEFILE	"\t$log_file     \t\tRun log for bam2vcf_RNA\n\n\n";


	print	READMEFILE	"NOTES\n\n"; 
	print	READMEFILE	"\t.bam files must have an associated .bai index file for loading into IGV\n";

	close (READMEFILE);

	&move_to_results_folder ("$readme_file");

	#########################################################################################


	##############################################
	# Delete these files at the end of each loop #
	##############################################
	&delete_file ("forRealigner.intervals");
	&delete_file ("depth.sample_cumulative_coverage_counts");
	&delete_file ("depth.sample_cumulative_coverage_proportions");
	&delete_file ("depth.sample_statistics");
	&delete_file ("*.junk");
	&delete_file ("$snps_vcf_idx");
	&delete_file ("$indels_vcf_idx");



	#RUN TIMER END

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
	print MAIL "Subject: BAM2VCF_RNA Run: $run_title  Sample: $sample_name has finished\n\n";
	## Mail Body
	print MAIL "bam2vcf_RNA PERL script version $version\n\n";
	print MAIL "Your next generation sequence analysis of $sample_name in run $run_title is complete\n\n";
	print MAIL "Run time so far : $run_time seconds\n";
	printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
	close(MAIL);

	

} # End of MAIN LOOP for ($loop_count=1;$loop_count <=$no_of_files;$loop_count++)


##############################################################
# Delete files at end of all loops instead of for each loop? #
##############################################################

for ($loop_count=$start_folder;$loop_count <=$no_of_files;$loop_count++)
{
	$bam_file = $bam_file_array[$loop_count];
	$sample_name = $sample_name_array[$loop_count];
	
	# File names
	$rg_bam  = "$run_title"."_"."$sample_name"."_rg.bam";
	$rg_bai  = "$run_title"."_"."$sample_name"."_rg.bai";
	
	if ($delete_intermediate_files eq "yes")
	{
		print "\n\nThis only occurs if delete_intermediate_files is yes.  It is $delete_intermediate_files\n\n";
		&delete_file("$rg_bam"); 
		&delete_file("$rg_bai");
		#&delete_file("$runtitle_dedup_bam");
		#&delete_file("$runtitle_dedup_bai");
		#&delete_file("$runtitle_dedup_bam_intervals");
		#&delete_file("$runtitle_clean_dedup_bam");
		#&delete_file("$runtitle_clean_dedup_bai");
		#&delete_file("$runtitle_clean_dedup_recal_bam");
		#&delete_file("$runtitle_clean_dedup_recal_bai");
		#&delete_file("$region_only_bam");
		#&delete_file("$region_only_bai");
		
		# Recal files
		&delete_file ("$runtitle_sample_recal_grp");
	}
	
} # end of deleting loop

##########################################################
#  Make a VCF file using all the BAM files in parallel   #
##########################################################

####################################################
# Open the list file to get the list of file names #
# (seem to lose the file names so try this)        #
####################################################


print "List file of file names: $list_file\n\n";

if ($read_file_method eq "multiple")
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



	###################
	# List file names #
	###################
	print "\n\nThere are $no_of_files files in this file of file names.\n\n";

	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		print "File $list_count	\t$bam_file_array[$list_count]\n";
	}


			
	print "\n\n";
	print "=========================================================\n";
	print "Stage 13: Making VCF files with all BAM files in parallel    \n";
	print "=========================================================\n\n";
	if ($make_parallel_VCF eq "parallel" || $make_parallel_VCF eq "parallel and consecutive")
	{

		# First get a list of the BAM files #
		for ($list_count=1;$list_count <=$no_of_files;$list_count++)
		{
		
			print "File count: $list_count\n";
			
			$results_folder = "results_"."$run_title"."_"."$list_count";
			print "\tResults folder: \t$results_folder\n";
			
			$bam_file = $bam_file_array[$list_count];
			print "\tBAM file:       \t$bam_file\n";
			
			$sample_name = $sample_name_array[$list_count];
			print "\tSample name:    \t$sample_name\n\n";
			

			################################################################
			# Make up name of cleaned dedupped recalibrated final bam file #
			################################################################
			$final_bam =  "$run_title"."_"."$sample_name"."_final.bam";
				
			$bam_file_for_unified_genotyper = "$results_folder/$final_bam";
			
			if (! -e "$bam_file_for_unified_genotyper")
			{
				$bam_missing_for_UG_count = $bam_missing_for_UG_count + 1;
			}
			
			$input_string = $input_string." -I $bam_file_for_unified_genotyper";
		}
		$input_string = $input_string." ";

		print "\n\n";
		print "========================================================\n";
		print "Stage 13a: Making SNP VCF with all BAM files in parallel \n";
		print "========================================================\n\n";
	
	
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $UG_region_string -glm SNP $input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o $snps_parallel_vcf  -S $GATK_validation_stringency $fix_quals_string","PARALLEL SNPS VCF");	

		&record_output_file_size ("$snps_parallel_vcf");	
		
		if ($bam_missing_for_UG_count == 0)
		{
			print "\n===================================================================";
			print "\nA SNP VCF file using all the BAM files in parallel has been created";
			print "\n===================================================================\n\n\n";
			
			print "\n* * * \n";
		}
		
		print "\n\n";
		print "\n\n";
		print "==========================================================\n";
		print "Stage 13b: Making Indel VCF with all BAM files in parallel \n";
		print "==========================================================\n\n";
		
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $UG_region_string -glm INDEL $input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf --max_alternate_alleles $max_alt_alleles -o $indels_parallel_vcf -S $GATK_validation_stringency $fix_quals_string","PARALLEL INDELS VCF");
		
		&record_output_file_size ("$indels_parallel_vcf");

		if ($bam_missing_for_UG_count == 0)
		{
			print "\n======================================================================";
			print "\nAn Indel VCF file using all the BAM files in parallel has been created";
			print "\n======================================================================\n\n\n";
			
			print "\n* * * \n";
		}
		
		if ($bam_missing_for_UG_count > 0)
		{
			print "\n=====================================================================";
			print "\nSNP and Indel files were not created as not all files were found";
			print "\nPlease check their location and use VCF_PARALLEL later";
			print "\n=====================================================================\n\n\n";
			
			print "\n* * * \n";
		}
		
		
		##############################################################
		# Annotate SNPs and Indels with the variant effect predictor #
		##############################################################
		
		if ($annotate_variants eq "yes")
		{
			###################
			# Annotate SNPs   #
			###################
		
			print "\n=====================";
			print "\n Annotating SNP calls";
			print "\n=====================\n\n";

			if ($use_defined_region eq "yes")
			{
				&run_unix_command("perl $ensembl_directory/variant_effect_predictor.pl --species $species --chr $chromosome --input $snps_parallel_vcf --format vcf --output_file $annotated_snps_parallel_vcf -verbose","VEP_SNPs");
			}
			if ($use_defined_region eq "no")
			{
				&run_unix_command("perl $ensembl_directory/variant_effect_predictor.pl --species $species --input $snps_parallel_vcf --format vcf --output_file $annotated_snps_parallel_vcf -verbose","VEP_SNPs");
			}
		
			&record_output_file_size ("$annotated_snps_parallel_vcf");
			
			print "\n==========================================";
			print "\n COMPLETED: SNP annotations have been made";
			print "\n==========================================\n";

			print "\n* * * \n";
			
			&move_to_results_folder ("$annotated_snps_parallel_vcf");


			###################
			# Annotate Indels #
			###################
		
			print "\n=======================";
			print "\n Annotating Indel calls";
			print "\n=======================\n\n";


			if ($use_defined_region eq "yes")
			{
				&run_unix_command("perl $ensembl_directory/variant_effect_predictor.pl --species $species --chr $chromosome --input $indels_parallel_vcf --format vcf --output_file $annotated_indels_parallel_vcf -verbose","VEP_Indels");
			}
			if ($use_defined_region eq "no")
			{
				&run_unix_command("perl $ensembl_directory/variant_effect_predictor.pl --species $species --input $indels_parallel_vcf --format vcf --output_file $annotated_indels_parallel_vcf -verbose","VEP_Indels");
			}
		
			&record_output_file_size ("$annotated_indels_parallel_vcf");
			
			print "\n============================================";
			print "\n COMPLETED: Indel annotations have been made";
			print "\n============================================\n";

			print "\n* * * \n";
			
			&move_to_results_folder ("$annotated_indels_parallel_vcf");
			

		} # if annotate_variants eq "yes"
		
		
		#####################################################
		# Move parallel VCF file to the last results folder #
		#####################################################
		&move_to_results_folder ("$snps_parallel_vcf");
		&move_to_results_folder ("$indels_parallel_vcf");
		
		
	} # End of if ($make_parallel_VCF eq "parallel" || $make_parallel_VCF eq "parallel and consecutive")

} # End of if ($read_file_method eq "multiple")

my $end_run = time();
$run_time = $end_run - our $start_run;

open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $email_address\n";
print MAIL "From: $from\n";
print MAIL "Subject: BAM2VCF_RNA: Run $run_title has finished\n\n";
## Mail Body
print MAIL "bam2vcf_RNA PERL script version $version\n\n";
print MAIL "Your next generation sequence analysis ($run_title) is complete\n\n";
print MAIL "For a command log see $command_log (very useful for debugging)\n\n";
print MAIL "For pipeline details see $log_file\n\n";
print MAIL "For file information see $readme_file\n\n";
print MAIL "Run time: ";
printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
close(MAIL);

print "\n\n";
print "######################################################\n";
print "# ANALYSIS COMPLETE! YOU HAVE BEEN NOTIFIED BY EMAIL #\n";
print "######################################################\n\n";

print "Run title: $run_title\n\n";
print "Run time: ";
printf "%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];


print "\n\nPlease check $command_log and $log_file for full details\n\n";

if ($bam_missing_for_UG_count > 0)
{
	print "\n#########################################################################";
	print "\nWhen the Unified Genotyper ran to make SNP and Indel files from multiple";
	print "\nBAM files, the VCF files were NOT created as not all files were found";
	print "\nPlease check their location and use VCF_PARALLEL later";
	print "\n##########################################################################\n\n\n";
}



#############################
# Turn off logging          #
#############################

close(STDOUT);
close (COMMAND_LOG);

&move_log ("$log_file");
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
	
	if (-e $file_to_be_moved)
	{
		$command = "mv $file_to_be_moved $results_folder/$file_to_be_moved";
		print COMMAND_LOG "$command\n";
		system("$command");
	}
	if (! -e $file_to_be_moved)
	{
		$command = "mv $file_to_be_moved $results_folder/$file_to_be_moved";
		print COMMAND_LOG "COULD NOT FIND $file_to_be_moved TO MOVE TO RESULTS FOLDER\n";
		print "COULD NOT FIND $file_to_be_moved TO MOVE TO RESULTS FOLDER\n";
	}

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
	print COMMAND_LOG "File: $loop_count/$no_of_files \tStep: $step/13  \n";
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
	printf COMMAND_LOG "  Total time:    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $end_time)[7,2,1,0];
	printf COMMAND_LOG "  Time for stage:\t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $stage_time)[7,2,1,0];
	
	$last_current_time = $current_time;
	
	####################################
	# Send e-mail if file size is zero #
	####################################
	if ($zero_size_reported eq "no")
	{
		if ($filesize < 1)
		{
			$zero_error_count = $zero_error_count + 1;
			
			$zero_size_reported = "yes";
			
			open(MAIL, "|/usr/sbin/sendmail -t");
			## Mail Header
			print MAIL "To: $email_address\n";
			print MAIL "From: $from\n";
			print MAIL "Subject: NGS ANALYSIS ZERO FILE SIZE: Program: bam2vcf_RNA  Run $run_title.  Sample $sample_name. File $outputfile.\n\n";
			## Mail Body
			print MAIL "bam2vcf_RNA script version $version\n\n";
			print MAIL "Run:    \t$run_title\n";
			print MAIL "Sample: \t$sample_name\n";
			print MAIL "File:   \t$outputfile\n\n";
			print MAIL "The output file $outputfile has zero file size or can't be found.\n\n";
			print MAIL "Something may be wrong.  Please check\n\n";
			close(MAIL);
		}
		
	} # if ($zero_size_reported eq "no";)
	
	if ($testing_mode eq "on")
	{
		print "\nPRESS 'RETURN' TO CONTINUE\n";
		$answer = <STDIN>;
	}

}


##############################################
#                                            #
# Subroutine to record size of output file   #
# (without sending an e-mail if it is zero)  #
##############################################

sub record_output_file_size_no_mail
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
