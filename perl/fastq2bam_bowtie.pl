#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	FASTQ2bam_bowtie        						                    #     
#									                                    #
#	PROCESS FASTQ FILES TO BAM FILES   using bowtie2 instead of bwa     #
#									                                    #
#########################################################################

#############################
# Mike Boursnell July 2013  #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# mike.boursnell@aht.org.uk #
#############################


use strict;
use Getopt::Std ;
use File::Basename ;
#]use warnings::unused;

# VERSION OF SOFTWARE #
my $version						= "1";


####################
# Define constants #
####################
#paths for software
my $bowtie_path					= "/opt/bowtie2/bowtie2"; # Note "bowtie" means "bowtie2"
my $tophat_path					= "/opt/tophat/tophat";
my $bowtie_index				= "";

my $gatk_directory				= "gatk"; # "gatk_v1" is the previous GATK version and "gatk" is the current GATK version 2.0
my $testing_mode				= "off";
my $use_base_recalibrator		= "yes";  # If this is set to "no" then the old CountCovariates is used instead
my $use_samtools_rmdup			= "no"; # This means use samtools rmdup rather than picard/MarkDuplicates

my $ValidationStringency		= "LENIENT";
my $picard_validation_stringency = "SILENT";
my $GATK_validation_stringency	= "LENIENT";
my $ensembl_directory			= "/opt/ensembl.new";
my $validate_bam_files			= "yes";
my $check_depth					= "no";
my $delete_intermediate_files	= "yes";
my $max_open_temp_files			= 7900;
my $platform					= "Illumina";
my $remove_duplicates			= "yes";
my $keep_all_for_first_file		= "yes"; # This keeps all intermediate files for first file, for error-checking

#Various Parameters (e.g. for the Unified Genotyper)
my $stand_emit_conf				= 30;
my $stand_call_conf				= 30;
my $max_alt_alleles				= 6; # This is how many alleles the INDEL caller in UnifiedGenotyper can allow
my $indels_caller				= "UnifiedGenotyper"; # alternatives are UnifiedGenotyper OR SomaticIndelDetector

####################
# Define variables #
####################
my $list_count					= 0; #Counter for the list of input files
my $no_of_files					= 0; #No of files in the list of input files (or no of lines if Paired Ends)
my $loop_count					= 0; #Loops round once for each of the files in the list of files
my $run_time					= 0;
my $next_folder					= 0;
my $start_folder				= 1;
my $array_size					= 0; # Size of item array to check there are two columns in the file
my $bam_missing_for_UG_count	= 0;
my $no_of_steps					= 0; # Number of steps in a stage
my $zero_error_count			= 0;

# New for RNA #
my $insert_size					= 0;
my $tophat_output_dir			= "";
my $tophat_string				= "";
my $tophat_string_1				= "";
my $tophat_string_2				= "";
my $read_length					= 0;
my $mate_inner_dist				= 0;
my $mate_std_dev				= 0;
my $calculate_inner_mate_distance = "";
my $mate_inner_string			= "";
my $zero_size_reported			= "";
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
my $input_file					= "";
my $list_file					= "";
my $read_file_method			= "";
my $answer						= "";
my $single_line					= "";
my $results_folder_bam2vcf		= "";  # Folder for bam2vcf output files
my $results_all_folder			= "";  # Folder for fastq2bam files (and also contains bam2vcf folders)
my $ref_seq_name				= "";  # Name of referecne sequence for pindel analysis
my $log_file					= "";
my $sample_name					= "";  # Name of each sample. Used to give the .vcf file an individual sample name
my $new_sample_name				= "";  # Sample name given in input file (optional)
my $use_preindexed_reference	= "";  # Can be 'yes' or 'no'  Flags whether you want to use pred-indexed ref seq such as Can Fam2 or EquCab1
my $start_time					= "";
my $end_time					= "";
my $current_time				= "";
my $last_current_time			= "";
my $stage_time					= ""; # Time for each stage of the pipeline
my $third_column_found			= "false";
my $second_column_found			= "false";
my $lib							= ""; #This is used by the ReadGroups section and i have made it let us record the reference sequence
my $input_string				= "";
my $running_input_string		= ""; # string of BAM file na mes for UG
my $species						= "";

my $prefix_name					= "";  # Name of bam_file ommitting the suffix .bam
my $bam_file_mates_fixed		= "";
my $temp_dir_string				= " -Djava.io.tmpdir=javatempdir";
my $tempdir						= "javatempdir";
my $use_defined_region			= "";
my $region						= "";
my $GATK_region_string			= ""; # String to be used with the -L option of GATK tools such as UnifiedGenotyper
my $samtools_region_string		= ""; # String to be used with samtools
my $UG_region_string 			= ""; # String to be used with Unified Genotyper (just narrows down to chromosome e.g. -L chr15)
my $GATK_known_sites_string		= ""; # String for known SNPs for BaseRecalibrator/CountCovariates
my $make_parallel_VCF			= ""; # This makes use of multiple BAM files in GATK to make a combined VCF output file
my $structural_variant_analysis	= "";
my $ref_genome					= "";
my $chromosome					= "";
my $annotate_variants			= "";
my $use_other_ref_sequence 		= "no";

#### File names FASTQ2BAM ####

my $accepted_hits_bam			= "";
my $accepted_hits_bai			= "";
my $accepted_hits_rg_bam		= "";
my $accepted_hits_rg_bai		= "";
my $deletions_bed				= "";
my $insertions_bed				= "";
my $junctions_bed				= "";
my $unmapped_bam				= "";

my $aln_sa_sai					= "";
my $aln_sa1_sai					= "";
my $aln_sa2_sai					= "";
my $aligned_sam					= "";
my $aligned_sorted_bam  		= "";
my $aligned_sorted_bai			= "";
my $aligned_sorted_rg_bam		= "";
my $aligned_sorted_rg_bai		= "";
my $aligned_sorted_viewed_bam	= "";
my $aligned_sorted_viewed_bai	= "";
my $command_log					= "";
my $title						= "";


#### File names BAM2VCF ####
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

my $region_only_bam					= ""; # After reducing BAM to a region
my $region_only_bai					= "";

my $final_bam						= ""; # Final BAM file after BAM2VCF pipeline has run
my $final_bai						= "";

my $final_tophat_bam				= ""; # Final BAM file after FASTQ2BAM pipeline has run
my $final_tophat_bai				= "";

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
my $validate_tophat_out			= "";
my $flagstat_tophat_out			= "";
my $validate_out				= "";
my $flagstat_out				= "";

#Other files
my $annotated_indels_vcf		= "";
my $annotated_snps_vcf			= "";
my $annotated_indels_parallel_vcf		= ""; # Version for multi-sample parallel VCf file
my $annotated_snps_parallel_vcf			= ""; # Version for multi-sample parallel VCf file
my $ref_dict					= ""; # This is the file like can fam 3.dict
my $insert_size_pdf				= ""; # Used by CollectInsertSizeMetrics
my $gc_bias_pdf					= ""; # Used by CollectGcBiasMetrics
my $bam_file_for_unified_genotyper = ""; # Used in parallel VCF section
my $dummy_dbsnp_file			= ""; # The dummy dbSNP file used by....
my $actual_dbsnp_file			= ""; # Correct dbSNP file for this chromosome (if it exists)
my $dbsnp_file					= ""; # dbsnp file used by GATK

my $runtitle_sample_recal_grp 	= ""; # For BaseRealigner
my $recal_file_count_covariates = ""; # For CountCovariates

my @bam_file_array				= ();
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


use Term::ANSIColor;
print color 'bold cyan';


if ($testing_mode eq "on"){print"\n\nTESTING MODE ON\n\n";}

print color 'reset';

use Term::ANSIColor;
print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      fastq2bam_bowtie      \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";

print "  - This program processes FASTQ files into BAM files using bowtie2 instaed of bwa.\n\n";


print "  - NOT YET EDITED FROM RNA VERSION\n\n";

exit;



print color 'reset';


#############################
# Name the analysis run     #
#############################

print "\n~~~~~~~~~~~~~~~";
print "\nYour details...";
print "\n~~~~~~~~~~~~~~~\n";

print "\nPlease enter a name for this analysis.\n";
print "(Keep if fairly short as it becomes part of the file name - and with no spaces):\n\n";
$run_title = <STDIN>;
chomp $run_title;


#########################################################
# Check if results folder with this name exists already #
#########################################################

$results_all_folder = "results_all_"."$run_title";
	
if (-e $results_all_folder)
{
	print "\n\n";
	print "###########\n";
	print "# WARNING #\n";
	print "###########\n\n";
	print "A results folder called $results_all_folder already exists\n\n";
	print "Choose a new run title or delete the existing folder.\n\n";
	
	exit;
}

$log_file = "$run_title"."_fastq2bam_RNA_log.out";

$| = 1;

open(STDERR, "| tee $log_file");


######################################
# E-mail address for notifications   #
######################################

print "\nPlease enter your email address:    ";
$email_address = <STDIN>;
chomp $email_address;

# Some short cuts
if ($email_address eq "m"){$email_address = 'mike.boursnell@aht.org.uk';}
if ($email_address eq "o"){$email_address = 'oliver.forman@aht.org.uk';}
if ($email_address eq "k"){$email_address = 'karen.steward@aht.org.uk';}


$platform = "Illumina";


##########################################
# ASK IF STARTING FILES ARE BAM OR FASTQ #
##########################################
$input_file = "fastq";



#####################################
# Ask if the data is SE or PE data? #
# ASSUME PE from now on             #
#####################################

#print "\n\n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
#print " Do you have a Paired-end or Single-end dataset?\n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
#print "   Enter 1 for SE\n";
#print "   Enter 2 for PE\n\n";
#$answer = <STDIN>;
#chomp $answer;
#if (substr($answer,0,1) eq "1" ){$data = "SE"}
#if (substr($answer,0,1) eq "2" ){$data = "PE"}

#ASSUME PE
$data = "PE";

############################################
# Data to calculate inner mate distance    #
############################################
if ($data eq "PE")
{
	print "\n\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "For TopHat, the 'Mate Inner Distance' can be specified\n";
	print "(This is the distance between the inner ends of the two reads)\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

	print "Do you want to enter data so this can be calculated?\n\n";
	print "(Note: the defaults work pretty well so it is not essential)\n\n";
	
	print "   Enter 1 for YES\n";
	print "   Enter 2 for NO\n\n";
	
	$answer = <STDIN>;
	chomp $answer;
	if (substr($answer,0,1) eq "1" ){$calculate_inner_mate_distance = "yes"}
	if (substr($answer,0,1) eq "2" ){$calculate_inner_mate_distance = "no"}

	if ($calculate_inner_mate_distance eq "yes")
	{
		##################################
		# Distance between read pairs    #
		##################################
		print "\n";
		print "Enter sequencing read length:     \t\t  ";

		$read_length = <STDIN>;
		chomp $read_length;

		##################################
		# Distance between read pairs    #
		##################################
		print "\n";
		print "Enter mean insert size:           \t\t";

		$insert_size = <STDIN>;
		chomp $insert_size;


		##################################
		# Distance between read pairs    #
		##################################
		print "\n";
		print "Enter standard deviation for insert size: (default=20):  ";;

		$mate_std_dev = <STDIN>;
		chomp $mate_std_dev;

		$mate_inner_dist = int ($insert_size - ($read_length * 2));
		
		print "\nThese values will be used:\n\n";
		print "\t--mate-inner-dist: \t$mate_inner_dist  (calculated as $insert_size - (2 x $read_length))\n";
		print "\t--mate-std-dev:    \t$mate_std_dev\n\n";
		
		$mate_inner_string = "--mate-inner-dist $mate_inner_dist --mate-std-dev $mate_std_dev";
		
	}
	
	if ($calculate_inner_mate_distance eq "no")
	{
		print "\nDefault values will be used (works pretty well)\n\n";
		print "\t--mate-inner-dist: \t50\n";
		print "\t--mate-std-dev:    \t20\n\n";
		
		$mate_inner_string = "";
	}
	
} # if PE


##################################
# Define data files              #
##################################

$use_preindexed_reference = "yes";

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

if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3/canfam3_dummy_DBSNP.vcf";$bowtie_index="/home/genetics/canfam3/bowtie_files/canfam3"}
if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/canfam2/canfam2.fasta"; $ref_seq_name = "canfam2"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam2/canfam2_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "3" ){$ref = "/home/genetics/equcab2/equcab2.fasta"; $ref_seq_name = "equcab2"; $species = "equus_caballus";$dummy_dbsnp_file = "/home/genetics/equcab2/equcab2_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "4" ){$ref = "/home/genetics/human/human.fasta"; $ref_seq_name = "human"; $species = "homo_sapiens";$dummy_dbsnp_file = "/home/genetics/human/human_dummy_DBSNP.vcf"}

if (substr($answer,0,1) eq "5" ){$ref = "/home/genetics/strep_equi/strep_equi.fasta"; $ref_seq_name = "s_equi"; $species = "streptococcus_equi";$dummy_dbsnp_file = "/home/genetics/strep_equi/strep_equi_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "6" ){$ref = "/home/genetics/strep_zoo/strep_zoo.fasta"; $ref_seq_name = "s_zoo"; $species = "streptococcus_zoo";$dummy_dbsnp_file = "/home/genetics/strep_zoo/strep_zoo_dummy_DBSNP.vcf"}


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

###################################################################
# Get bit before name of ref sequence to delete "test.dict" later #
###################################################################
$ref_prefix = &get_prefix ($ref);



###################################################################
# Check if REF.DICT files exist (if a pre-indexed file is chosen) #
###################################################################

if ((! -e "$ref") || (! -e "$ref.bwt") || (! -e "$ref.rbwt") || (! -e "$ref.rpac") || (! -e "$ref.pac") || (! -e "$ref.rsa") || (! -e "$ref.sa") || (! -e "$ref.amb") || (! -e "$ref.dict") || (! -e "$ref.ann") || (! -e "$ref.fai"))
{ 
	print "##############################\n";
	print "#  REFERENCE SEQUENCE ERROR  #\n";
	print "##############################\n\n";
	print "\nNot all the correct files for an indexed refererence sequence exist.\n\n";
	print "You need the following files:  .dict, .rbwt, .bwt, .rpac, .pac, .rsa, .sa, .amb, .ann, .fai\n\n";
	print "This suggests that there is not already a pre-indexed reference file\n\n";
	
	exit;
} 


############################################################################
# Remove duplicates                                                        #
# Would you like to remove duplicates - noramlly YES but might not want to #
############################################################################


print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Would you like to remove duplicate reads?\n";
print "(This is normally strongly advised, but for special purposes you might not want to)\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";


print "   Enter 1 for YES\n";
print "   Enter 2 for NO\n\n";

$answer = <STDIN>;
chomp $answer;

if ($answer eq ""){$answer = "1"} # default

if (substr($answer,0,1) eq "1"){$remove_duplicates = "yes"}
if (substr($answer,0,1) eq "2"){$remove_duplicates = "no"}


###########################################################################
# Region preferences                                                      #
# Would you like to focus the analysis on a specific region of the genome #
###########################################################################


#print "\n\n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
#print "Would you like to focus on alignments to specific region of the genome?\n";
#print "(e.g. a sequence-captured region of a larger genome)\n";
#print "(This is strongly advised as it speeds up the whole process)\n\n";
#print "  NB: For bacterial genomes use 'chr1'\n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";#


#print "   Enter 1 for YES\n";
#print "   Enter 2 for NO\n\n";

#$answer = <STDIN>;
#chomp $answer;

#if ($answer eq ""){$answer = "1"} # default

#if (substr($answer,0,1) eq "1"){$use_defined_region = "yes"}
#if (substr($answer,0,1) eq "2"){$use_defined_region = "no"}





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
	
	$GATK_region_string = "-L $region";
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


#####################################
# Input file names for fastq option #
#####################################


if ($read_file_method eq "single")
{
	if ($data eq "SE")
	{
		print "\nPlease input the name of your fastq sequence file (eg reads):      ";
		$reads = <STDIN>;
		chomp $reads;

		# Assign to array of file names but just use first element #
		$reads1_file_array[1]=$reads;
		$no_of_files = 1;
	}
	
	if ($data eq "PE")
	{
		until (-e $reads)
		{
			print "\nPlease input the name of your FIRST fastq sequence file (eg reads_01):      ";
			$reads = <STDIN>;
			chomp $reads;
			
			if ($reads eq "ls"){print "\n";system ("ls *.fastq")}
			
			if ($reads ne "ls")
			{
				if (! -e $reads){print "\n\n>>>>>>>>  File $reads not found.  Try again.  <<<<<<<<\n\n";}
			}
		}
		
		until (-e $reads2)
		{
			print "\nPlease input the name of your SECOND fastq sequence file (eg reads_02):      ";
			$reads2 = <STDIN>;
			chomp $reads2;
			
			
			if ($reads2 eq "ls"){print "\n";system ("ls *.fastq")}
			
			if ($reads2 ne "ls")
			{
				if (! -e $reads2){print "\n\n>>>>>>>>  File $reads2 not found.  Try again.  <<<<<<<<\n\n";}
			}
		}
		
		# Add .fastq suffixes if necessary
		#if (index($reads,".fastq") == -1){$reads = $reads.".fastq"}
		#if (index($reads2,".fastq") == -1){$reads2 = $reads2.".fastq"}
		
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
			print "\nPlease input the name of your file with a list of file names of the 'reads' files:      ";
		}
		if ($data eq "PE")
		{
			print "\nPlease input the name of your file with a list of file names of the 'reads' files\n\n";
			print "(The two paired-end 'reads' file names must be on the same line separated by a TAB)::      ";
		}
		
		
		$list_file = <STDIN>;
		chomp $list_file;
		
		if ($list_file eq "ls"){print "\n";system ("ls *.txt"); print "\n";}
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

		if ($data eq "SE")
		{
			@item=split(/\t/,$single_line);
			
			$array_size = scalar @item;
			
			if ($array_size == 1)
			{
				$reads = $single_line;
				
				$reads1_file_array[$list_count]=$reads;
			} # array size 1
			
			if ($array_size == 2)
			{
				$reads = $item[0];
				$new_sample_name = $item[1];
				
				
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
				
				$reads1_file_array[$list_count] = $reads;
				$reads2_file_array[$list_count] = $reads2;
				$sample_name_array[$list_count] = $new_sample_name;
				$third_column_found = "true";
			}
			

		}

		$list_count=$list_count + 1;
	}

	close LIST;

	$no_of_files=$list_count - 1;
	
	
	########################
	# Check if files exist #
	########################
	
	print "\n\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "Checking input files...\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~\n";
	
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
			
			###########################################
			# Check the two fastq files are different #
			###########################################
			if ($reads1_file_array[$list_count] eq $reads2_file_array[$list_count])
			{
				print "\n\n\nThe two FASTQ files cannot be the same!\n\n";
				print "$reads1_file_array[$list_count]  is the same as $reads2_file_array[$list_count]\n\n";
				print "PLEASE FIX THIS\n\n";
				exit;
			}
		}
	}
	

	
	###################
	# List file names #
	###################
	print "\nThere are $no_of_files pairs of 'reads' files in this file of file names.\n\n";
	
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
	 
	if (lc $proceed eq "q")
	{
		exit;
	} 

} # End of if ($read_file_method eq "multiple")

	

print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Please enter the memory setting for this analysis\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

print "Memory setting (default is -Xmx4g):      ";
$mem = <STDIN>;
chomp $mem;
 
if ($mem eq "")
{
 $mem = "-Xmx4g"
}


#################################
# Annotate variants preferences #
#################################

print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " Would you like to annotate your SNP and INDEL calls using Ensembl Variant Effect Predictor?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

print "   Enter 1 for YES\n";
print "   Enter 2 for NO\n\n";

$answer = <STDIN>;
chomp $answer;


if (substr($answer,0,1) eq "1"){$annotate_variants = "yes"}
if (substr($answer,0,1) eq "2"){$annotate_variants = "no"}
if ($answer eq ""){$annotate_variants = "yes"}





$command =  "clear";
#system("$command");   TEMPORARY

print "\n\n\n";
print color 'bold magenta';

print "================================================\n";
print "SUMMARY - PLEASE CHECK THESE DETAILS CAREFULLY!!\n";
print "================================================\n\n";

print color "bold green";
print "YOUR DETAILS :\n\n";

print color "bold white";
print "Name of this analysis:\t$run_title\n\n"; 

print "Your email address:\t$email_address\n\n";

print color "bold green";
print "\nSETTINGS :\n\n";
print color "bold white";

if ($data eq "SE")        {print "Single-end analysis    \t\t\tYES\n\n";}
if ($data eq "PE")        {print "Paired-end analysis    \t\t\tYES\n\n";}

if ($remove_duplicates eq "yes"){print "Remove duplicates\t\t\tYES\n\n";}
if ($remove_duplicates eq "no") {print "Do not remove duplicates\t\t\tNO\n\n";}

#if ($use_base_recalibrator eq "no"){print "Use CountCovariates from old GATK\n\n";}
#if ($use_base_recalibrator eq "yes"){print "Use BaseRecalibrator from new GATK\n\n";}

if ($annotate_variants eq "yes")        {print "Annotate variants     \t\t\tYES\n\n";}
if ($annotate_variants eq "no")         {print "Annotate variants     \t\t\tNO\n\n";}

print "Memory setting: \t\t\t$mem\n\n";

print "UnifiedGenotyper constants:\n\n";

print "    --stand_emit_conf:\t\t\t$stand_emit_conf\n";
print "    --stand_call_conf:\t\t\t$stand_call_conf\n";
print "    --max_alt_alleles:\t\t\t$max_alt_alleles\n\n";


print color "bold green";
print "\nDATA FILES :\n\n";
print color "bold white";

print "Reference sequence: \t\t\t$ref\n\n";

print "Bowtie index files: \t\t\t$bowtie_index\n\n";

if ($use_defined_region eq "yes")
{
	print "Defined region:  \t\t\t$region\n\n";
	print "  GATK region string:      \t\t\t$GATK_region_string\n";
	print "  UG region string:        \t\t\t$UG_region_string\n";
	print "  Samtools region string:  \t\t\t$samtools_region_string\n\n";
}
if ($use_defined_region eq "no")
{
	print "Defined region:\t\t\t\tNO\n\n";
}

print "SNP file for BaseRecalibrator:\t\t$dbsnp_file\n\n";

print "Number of FASTQ files to analyse:\t$no_of_files\n\n";
	

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	if ($data eq "SE")
	{
		print "  Reads file $list_count	\t$reads1_file_array[$list_count]\n";
	}
	if ($data eq "PE")
	{
			if ($third_column_found eq "true")
			{
				print "Pair of FASTQ files:  $list_count	\t$reads1_file_array[$list_count]  \t$reads2_file_array[$list_count]\t\t$sample_name_array[$list_count]\n";
			}
			if ($third_column_found eq "false")
			{
				print "Pair of FASTQ files:  $list_count	\t$reads1_file_array[$list_count]  \t$reads2_file_array[$list_count]\n";
			}
	}
	
			
}



print "\n\n\n";
print color 'bold white';
print "======================================================================\n";
print "Please press 'ENTER' to proceed with the analysis run (or 'Q' to quit)\n";
print "======================================================================\n\n";
print color 'reset';


$proceed = <STDIN>;
chomp $proceed; 
 
if (lc $proceed eq "q")
{
	exit;
} 
	



#########################
# open Command Log file #
#########################
$command_log = "$run_title"."_fastq2bam_RNA_command_log.out";
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

print COMMAND_LOG "COMMAND LOG for fastq2bam_RNA version $version\n\n";
print COMMAND_LOG "Analysis name:         \t$run_title\n\n";
print COMMAND_LOG "Reference sequence:    \t$ref\n";
print COMMAND_LOG "Ref seq name:          \t$ref_seq_name\n\n";
print COMMAND_LOG "Bowtie2 index files:   \t$bowtie_index\n\n";

if ($use_defined_region	eq "yes")
{
	print COMMAND_LOG "Defined region:  \t\t\t$region\n\n";
	print COMMAND_LOG "  GATK region string:     \t\t\t$GATK_region_string\n";
	print COMMAND_LOG "  UG region string:       \t\t\t$UG_region_string\n";
	print COMMAND_LOG "  Samtools region string: \t\t\t$samtools_region_string\n\n";
}
if ($use_defined_region	eq "no"){print COMMAND_LOG "No genome region defined\n\n";}

print COMMAND_LOG "Remove duplicates:       \t$remove_duplicates\n\n";

print COMMAND_LOG "SNP file for BaseRecalibrator:\t\t$dbsnp_file\n\n";

if ($use_base_recalibrator eq "yes"){print COMMAND_LOG "Used BaseRecalibrator\n\n";}
if ($use_base_recalibrator eq "no"){print COMMAND_LOG "Used CountCovariates\n\n";}

print COMMAND_LOG "List of input files\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	if ($data eq "SE")
	{
		print COMMAND_LOG "FASTQ FILE $list_count	\t$reads1_file_array[$list_count]\n";
	} # SE
	
	if ($data eq "PE")
	{
		if ($third_column_found eq "true")
		{
			print COMMAND_LOG "FASTQ files $list_count:	\t$reads1_file_array[$list_count]  \t$reads2_file_array[$list_count]\t\t$sample_name_array[$list_count]\n";
		}
		if ($third_column_found eq "false")
		{
			print COMMAND_LOG "FASTQ files $list_count:	\t$reads1_file_array[$list_count]  \t$reads2_file_array[$list_count]\n";
		}
	} # PE
}

print COMMAND_LOG "\n";


########################################
# Create results_all folder            #
# This holds the files from fastq2bam  #
# and the individual bam2vcf folders   #
########################################
&run_unix_command("mkdir $results_all_folder","Make results folder");	


##############################################################################
##############################################################################
# Start MAIN LOOP.  This runs through each of the files in the list of files #
##############################################################################
##############################################################################

$last_current_time = time();



for ($loop_count=1;$loop_count <=$no_of_files;$loop_count++)
{

	$no_of_steps=5; # fastq2bam
	
	######################################################
	# Assign correct name to $reads or $bam_file         #
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
	
	print "Loop count: $loop_count\n";
	#print "sample_name_array: $sample_name_array[$loop_count]\n";
	print "Sample name: $sample_name\n\n";
	
	
	#######################################################################
	# Set up various File names which vary for each loop  from FASTQ2BAM  #
	#######################################################################

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
	
	$final_tophat_bam = "$run_title"."_"."$sample_name"."_tophat.bam"; # Final BAM of TopHat stage
	$final_tophat_bai = "$run_title"."_"."$sample_name"."_tophat.bai";
	
	$validate_tophat_out = "$run_title"."_"."$sample_name"."_validate_summary_tophat.out";
	$flagstat_tophat_out = "$run_title"."_"."$sample_name"."_flagstat_tophat.out";
	
	
	######################################################################
	# Set up various File names which vary for each loop from BAM2VCF    #
	######################################################################
	
	$bai_file= "$prefix_name".".bai";
	$bam_bai_file = "$prefix_name".".bam.bai";
	
	# After AddReadGroups
	$rg_bam  = "$run_title"."_"."$sample_name"."_rg.bam";
	$rg_bai  = "$run_title"."_"."$sample_name"."_rg.bai";
	
	# After MarkDuplicates (deduped)
	$runtitle_dedup_bam = "$run_title"."_"."$sample_name".".dedup.bam";
	$runtitle_dedup_bai = "$run_title"."_"."$sample_name".".dedup.bai";
	
	# After RealignerTargetCreator
	$runtitle_dedup_bam_intervals = "$run_title"."_"."$sample_name".".dedup.bam.intervals";
	
	# After IndelRealigner (cleaned)
	$runtitle_clean_dedup_bam = "$run_title"."_"."$sample_name".".clean.dedup.bam";
	$runtitle_clean_dedup_bai = "$run_title"."_"."$sample_name".".clean.dedup.bai";
	
	# After CountCovariates/BaseRecalibrator
	$runtitle_sample_recal_csv  = "$run_title"."_"."$sample_name".".recal.csv";
	$runtitle_sample_recal_grp  = "$run_title"."_"."$sample_name".".recal.grp"; # Use for BaseRecalibrator
	$recal_file_count_covariates = "$run_title"."_"."$sample_name"."recal.csv"; # Use for CountCovariates
	
	# After TableRecalibration (recal)
	$runtitle_clean_dedup_recal_bam =  "$run_title"."_"."$sample_name".".clean.dedup.recal.bam";
	$runtitle_clean_dedup_recal_bai =  "$run_title"."_"."$sample_name".".clean.dedup.recal.bai";
	
	# After reducing BAM to a region
	$region_only_bam = "$run_title"."_"."$sample_name"."_region.bam";
	$region_only_bai = "$run_title"."_"."$sample_name"."_region.bai";
	
	# After all stages have finished, the BAM file is renamed to final_bam
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
	$annotated_indels_vcf	= "annotated_indels_"."$sample_name.vcf";
	$annotated_snps_vcf	= "annotated_snps_"."$sample_name.vcf";
	$insert_size_pdf = "$sample_name"."_insert_size.pdf";
	$gc_bias_pdf = "$sample_name"."_gc_bias.pdf";
	$validate_out = "$run_title"."_"."$sample_name"."_validate.out";
	$flagstat_out = "$run_title"."_"."$sample_name"."_flagstat.out";
	
	#file names
	
	$runtitle_sample_metrics = "$run_title"."_"."$sample_name".".metrics";

	#RNA file names
	$accepted_hits_bam = "$sample_name"."_accepted_hits.bam";
	$accepted_hits_bai = "$sample_name"."_accepted_hits.bai";
	$accepted_hits_rg_bam = "$sample_name"."_accepted_hits_rg.bam";
	$accepted_hits_rg_bai = "$sample_name"."_accepted_hits_rg.bai";
	$deletions_bed = "$sample_name"."_deletions.bed";
	$insertions_bed = "$sample_name"."_insertions.bed";
	$junctions_bed = "$sample_name"."_junctions.bed";
	$unmapped_bam = "$sample_name"."_unmapped.bam";
	
	$tophat_output_dir = "tophat_out_"."$sample_name";
	
	
	####################################################################
	# Stage 1 Uses TopHat to produce aligned SAM file                  #
	#                                                                  #
	# Makes files                                         #
	####################################################################

	#####################
	# SINGLE END DATA ! #
	#####################

	if  ($data eq "SE")
	{		

		print"\n============================================================================================================";
		print "\nFile $loop_count/$no_of_files   fastq2bam_RNA Step 1 Running TopHat to produce aligned BAM file - SE";
		print "\n===========================================================================================================\n\n\n";
		
		
		&run_unix_command("$tophat_path $mate_inner_string --output-dir $tophat_output_dir --no-coverage-search -x $bowtie_index -U $reads -S $aligned_sam","fastq2bam_RNA 1 SE");

		#$BT2_HOME/bowtie2 -x lambda_virus -U $BT2_HOME/example/reads/reads_1.fq -S eg1.sam
		
		&record_output_file_size("$tophat_output_dir/$accepted_hits_bam");
		
		
		print "\n==============================================================================================================";
		print "\nFile $loop_count/$no_of_files   fastq2bam_RNA Step 1 Running TopHat to produce aligned BAM file - SE";
		print "\n==============================================================================================================\n\n\n";

		&test_mode_subroutine;
		
	}

	#####################
	# PAIRED END DATA ! #
	#####################

	if  ($data eq "PE")
	{
	
		#####################################################################
		# Read 1 and 2 (unlike bwa bowtie does both reads at the same time) #
		#####################################################################
		
		print "\n===========================================================================================================";
		print "\nFile $loop_count/$no_of_files   fastq2bam_RNA Step 1a Running TopHat to produce aligned BAM file - PE";
		print "\n===========================================================================================================\n\n\n";
		
		
		&run_unix_command("$tophat_path $mate_inner_string --output-dir $tophat_output_dir --no-coverage-search $bowtie_index $reads $reads2","fastq2bam_RNA 1 PE");

		###############################################
		# Rename output files with prefix from sample #
		###############################################
		&run_unix_command("mv $tophat_output_dir/accepted_hits.bam $tophat_output_dir/$accepted_hits_bam","Rename accepted_hits file");
		&run_unix_command("mv $tophat_output_dir/deletions.bed $tophat_output_dir/$deletions_bed","Rename deletions_bed file");
		&run_unix_command("mv $tophat_output_dir/insertions.bed $tophat_output_dir/$insertions_bed","Rename insertions_bed file");
		&run_unix_command("mv $tophat_output_dir/junctions.bed $tophat_output_dir/$junctions_bed","Rename junctions_bed file");
		&run_unix_command("mv $tophat_output_dir/unmapped.bam $tophat_output_dir/$unmapped_bam","Rename unmapped_bam file");

		
		&record_output_file_size("$tophat_output_dir/$accepted_hits_bam");
		
		
		print "\n==============================================================================================================";
		print "\nFile $loop_count/$no_of_files   fastq2bam_RNA Step 1a Running TopHat to produce aligned BAM file";
		print "\n==============================================================================================================\n\n\n";

		&test_mode_subroutine;
		
	}

	&problem_check;

	########################################################
	# Stage 2 Make an index file for this BAM file         #
	########################################################
	
	print "\n\n============================================================";
	print "\nFile $loop_count/$no_of_files   fastq2bam_RNA stage 2/5:  New Index for BAM file being created";
	print "\n============================================================\n\n";


	&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/BuildBamIndex.jar I=$tophat_output_dir/$accepted_hits_bam  VALIDATION_STRINGENCY=LENIENT","fastq2bam_RNA 2");

	&record_output_file_size ("$tophat_output_dir/$accepted_hits_bai");
	
	print "\n============================================================";
	print "\nFile $loop_count/$no_of_files   fastq2bam_RNA stage 2/5:  New Index for BAM file created";
	print "\n============================================================\n\n\n";
	
	
	#MIGHT NEED TO REORDER THE BAM FILE USING PICARD/REORDERSAM.....
	
	###########################################################################################
	# Stage 3: Use picard/AddOrReplaceReadGroups to update ReadGroup info in BAM file         #
	###########################################################################################

	print	"\n======================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   fastq2bam_RNA stage 3/5: Using picard/AddOrReplaceReadGroups to update ReadGroup info in BAM file";
	print 	"\n======================================================================================================\n\n\n";

	&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/AddOrReplaceReadGroups.jar I=$tophat_output_dir/$accepted_hits_bam O=$tophat_output_dir/$accepted_hits_rg_bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true RGID=$sample_name  RGLB=$lib  RGPL=illumina  RGPU=$sample_name  RGSM=$sample_name","fastq2bam_RNA 2 1");

	&record_output_file_size ("$tophat_output_dir/$accepted_hits_rg_bam");

	print	"\n=======================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   fastq2bam_RNA stage 3/5 COMPLETED picard/AddOrReplaceReadGroups updated ReadGroup info in BAM file";
	print 	"\n=======================================================================================================\n\n\n";


	print "\n* * * \n";


	#######################################################
	# Stage 4  Use picard ValidateSamFile Summary         #
	#######################################################
	
	print "\n\n==============================================================";
	print "\nFile $loop_count/$no_of_files   Step 3/5 Picard ValidateSamFile to be carried out";
	print "\n==============================================================\n\n";
	
	
	&run_unix_command("java $mem -jar /opt/picard/ValidateSamFile.jar I=$tophat_output_dir/$accepted_hits_rg_bam O=$tophat_output_dir/$validate_tophat_out MODE=SUMMARY MAX_OPEN_TEMP_FILES=$max_open_temp_files","fastq2bam_RNA 3");

	
	&record_output_file_size("$tophat_output_dir/$validate_tophat_out");

	
	print "\n======================================================";
	print "\nFile $loop_count/$no_of_files   Step 3/5 Picard ValidateSamFile carried out";
	print "\n======================================================\n\n\n";
	
	
	################################################
	# Stage 5 Use samtools flagstat                #
	################################################
	
	print "\n\n=======================================================";
	print "\nFile $loop_count/$no_of_files   Step 5/5 samtools flagstat to be carried out";
	print "\n=======================================================\n\n";
	
	
	&run_unix_command("/opt/samtools/samtools flagstat $tophat_output_dir/$accepted_hits_rg_bam > $tophat_output_dir/$flagstat_tophat_out","fastq2bam_RNA 4");

	&record_output_file_size("$tophat_output_dir/$flagstat_tophat_out");

	
	print "\n==================================================";
	print "\nFile $loop_count/$no_of_files   Step 5/5 samtools flagstat carried out";
	print "\n==================================================\n\n\n";
	
	
	##########################
	# Final TopHat BAM files #
	# (maybe rename later)
	##########################
	
	&run_unix_command("mv $tophat_output_dir/$accepted_hits_rg_bam $tophat_output_dir/$final_tophat_bam","Rename final BAM file");
	&run_unix_command("mv $tophat_output_dir/$accepted_hits_rg_bai $tophat_output_dir/$final_tophat_bai","Rename final BAMIfile");
	
	$final_tophat_bam = "$tophat_output_dir/$accepted_hits_rg_bam";
	$final_tophat_bai = "$tophat_output_dir/$accepted_hits_rg_bai";
	
	
	##########
	# README #
	##########

	open (READMEFILE, '>>README_tophat.out'); 

	print READMEFILE "###########################################################\n";
	print READMEFILE "Summary of fastq2bam_RNA (fastq2bam section) results files\n";
	print READMEFILE "###########################################################\n\n";

	print READMEFILE "Program:\t\t\tfastq2bam_RNA (fastq2bam section)\tVersion: $version\n\n";
		
	print READMEFILE "Run title:\t\t\t$run_title\n\n";
	
	print READMEFILE "Reference sequence:\t$ref_seq_name\n\n";
	
	print READMEFILE "\nBAM ALIGNMENT FILE\n\n";

	print READMEFILE "\t$final_tophat_bam   (Note: tophat means this BAM file comes from the TopHat stage)\n\n\n";
	
	print READMEFILE "To check this BAM file look at these two files:\n\n";
	
	print READMEFILE "\t$validate_tophat_out\n\n";
	print READMEFILE "\t$flagstat_tophat_out\n\n\n";
	
	print READMEFILE "This BAM file may need to go through the NGS pipeline using the bam2vcf section\n\n";

	close (READMEFILE);

	&move_to_results_all_folder ("README_tophat.out");

	
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
	print MAIL "Subject: FASTQ2BAM_RNA: Sample $sample_name in run $run_title has finished\n\n";
	## Mail Body
	print MAIL "fastq2bam_RNA PERL script version $version\n\n";
	print MAIL "The processing of FASTQ files into BAM files, of $sample_name in run $run_title is complete\n\n";
	
	print MAIL "Use bam2vcf_RNA for the next stage of the processing.\n\n";
	
	print MAIL "Run time so far : $run_time seconds\n";
	printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
	close(MAIL);

	print COMMAND_LOG "\n\n";
	print COMMAND_LOG "========================================================================\n";
	print COMMAND_LOG "                       End of loop $loop_count                          \n";
	print COMMAND_LOG "========================================================================\n\n\n";
} # End of MAIN LOOP




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
print MAIL "Subject: FASTQ2BAMA: Run $run_title has finished\n\n";
## Mail Body
print MAIL "fastq2bam_RNA PERL script version $version\n\n";
print MAIL "The processing of FASTQ files into BAM files, in run $run_title is complete\n\n";
print MAIL "Use bam2vcf_RNA for the next stage of the processing.\n\n";
print MAIL "For each sample file information see README.out\n\n";
print MAIL "For a list of commands see $command_log and $log_file\n\n";
print MAIL "Also check $validate_out and $flagstat_out for checks on the BAM file\n\n";

printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
close(MAIL);

print "\n\n";
print "####################################################################\n";
print "# FASTQ2BAM_RNA ANALYSIS COMPLETE! YOU HAVE BEEN NOTIFIED BY EMAIL #\n";
print "####################################################################\n\n";

print "Run title: $run_title\n\n";
print "Run time: ";
printf "%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];

print "\n\nResults folder:                 \t$results_all_folder";

print "\n\nFor details of the analysis run:\tcheck $command_log and $log_file for full details\n\n";

print "For checks on the BAM file produced:\tcheck $validate_out and $flagstat_out\n\n";


#############################
# Turn off logging          #
#############################

close(STDERR);
close (COMMAND_LOG);

&move_to_results_all_folder("$log_file");
&move_to_results_all_folder("$command_log");



#############################################
#                                           #
# Subroutine to move file to results folder #
#                                           #
#############################################

sub move_to_results_folder
{

	my $file_to_be_moved = "";	

	$file_to_be_moved = $_[0];
	$command = "mv  $file_to_be_moved $results_folder_bam2vcf/$file_to_be_moved";
	print("\n$command\n");
	print COMMAND_LOG ("$command\n");
	system("$command");

}

#################################################
# Subroutine to move file to results_all folder #
#################################################

sub move_to_results_all_folder
{

	my $file_to_be_moved = "";	

	$file_to_be_moved = $_[0];
	$command = "mv  $file_to_be_moved $results_all_folder/$file_to_be_moved";
	print("\n$command\n");
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
	
	if ($loop_count > 0 ) {print COMMAND_LOG "File: $loop_count/$no_of_files \tStep: $step/$no_of_steps\n";}
	if ($loop_count == 0 ) {print COMMAND_LOG "Step: $step/$no_of_steps\n";}
	
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
			print MAIL "Subject: NGS ANALYSIS fastq2bam_RNA ZERO FILE SIZE: Run $run_title.  Sample $sample_name. File $outputfile.\n\n";
			## Mail Body
			print MAIL "fastq2bam_RNA script version $version\n\n";
			print MAIL "Run:    \t$run_title\n";
			print MAIL "Sample: \t$sample_name\n";
			print MAIL "File:   \t$outputfile\n\n";
			print MAIL "The output file $outputfile has zero file size or can't be found.\n\n";
			print MAIL "Something may be wrong.  Please check\n\n";
			close(MAIL);
		}
		
	} # if ($zero_size_reported eq "no";)
	
}
	
	

##############################################
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
# Subroutine to move log     			    #
#############################################

sub move_log
{
	my $file_to_be_moved = "";	

	$file_to_be_moved = $_[0];
	$command = "mv  $file_to_be_moved $results_folder_bam2vcf/$file_to_be_moved";
	system("$command");

}

#############################################
# Subroutine to rename files   			    #
#############################################

sub rename_file
{
	my $old_name = "";
	my $new_name = "";	

	$old_name = $_[0];
	$new_name = $_[1];
	
	if (-e "$old_name")
	{
		$command = "mv  $old_name $new_name";
		system("$command");
		print COMMAND_LOG "Rename file using command $command\n";
	}
	if (!-e "$old_name")
	{
		print COMMAND_LOG "Failed to find $old_name to rename as $new_name\n";
	}
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


##############################################
# Carries out error check and e-mails user   #
##############################################
sub problem_check
{
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

sub show_start_of_step
{
	my $program_details 	= "";
	my $info_string	 		= "";
	my $print_string		= "";
	my $count				= 0;
	
	$program_details = $_[0];
	$info_string = $_[1];
	
	$print_string = "File: $loop_count/$no_of_files  $program_details  $info_string";
	
	for ($count=1;$count < length ($print_string);$count++){print "=";}
	print "\n$print_string\n";
	for ($count=1;$count < length ($print_string);$count++){print "=";}
	print "\n\n";

}	

