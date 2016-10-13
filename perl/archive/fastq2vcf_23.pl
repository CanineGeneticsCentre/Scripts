#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	FASTQ2VCF           						                        #     
#									                                    #
#	PROCESS FASTQ FILES TO VCF FILES           	                        #
#									                                    #
#########################################################################

##############################
# Mike Boursnell Nov 2013    #
# Animal Health Trust        #
# Newmarket                  #
# UK                         #
# mike.boursnell@aht.org.uk  #
#                            #
# from an original script by #
# Oliver Forman              #
# oliver.forman@aht.org.uk   #
##############################

use strict;
use Getopt::Std ;
use File::Basename;
use Term::ANSIColor;
use Cwd;

# VERSION OF SOFTWARE #
my $version						= "23";


####################
# Define constants #
####################
#paths for software
my $bwa_path					="/opt/bwa/bwa"; # This is the current version of bwa. Older versions are in folders labeled /opt/bwa062 or whatever...

my $no_of_threads				= "8";
my $gatk_directory				= "gatk"; # "gatk_v1" is the previous GATK version and "gatk" is the current GATK version 2.0
my $testing_mode				= "off";
my $use_samtools_rmdup			= "no"; # This means use samtools rmdup rather than picard/MarkDuplicates
my $from 						= 'NGS_analysis@samba64.aht.org.uk'; # Who e-mails come from
my $bwa_alignment_method		= "mem"; # replaces aln
my $variant_caller				= "UnifiedGenotyper"; # could be HaplotypeCaller
my $ValidationStringency		= "LENIENT";
my $picard_validation_stringency = "SILENT";
my $GATK_validation_stringency	= "LENIENT";
my $ensembl_directory			= "/opt/ensembl.new";
my $validate_bam_files			= "yes";
my $check_depth					= "no";
my $max_open_temp_files			= 7900;
my $remove_duplicates			= "yes";
my $keep_all_for_first_file		= "no"; # This keeps all intermediate files for first file, for error-checking
my $delete_intermediate_files	= "yes"; # set to 'no' for de-bugging

my $temp_dir_string				= " -Djava.io.tmpdir=javatempdir";
my $tempdir						= "javatempdir";

#Various Parameters (e.g. for the GATK Haplotype Caller)
my $stand_emit_conf				= 30;
my $stand_call_conf				= 30;
my $max_alt_alleles				= 6; # This is how many alleles the INDEL caller in UnifiedGenotyper can allow


####################
# Define variables #
####################
my $list_count					= 0; #Counter for the list of input files
my $no_of_files					= 0; #No of files in the list of input files (or no of lines if Paired Ends)
my $loop_count					= 0; #Loops round once for each of the files in the list of files
my $array_size					= 0; # Size of item array to check there are two columns in the file
my $bam_missing_for_variant_calling	= 0;
my $no_of_steps					= 0; # Number of steps in a stage
my $zero_error_count			= 0;
my $file_count					= 0;
my $check_count					= 0;
my $pos							= 0;
my $no_of_files_in_dir			= 0;
my $filesize_main				= 0;

#Boolean
my $use_default_stand_values	= ""; # 'yes' or 'no'
my $fix_misencoded_qual_scores	= ""; # 'yes' or 'no'
my $use_defined_region			= ""; # yes or no
my $annotate_variants			= ""; # yes or no
my $use_other_ref_sequence 		= "no"; # yes or no
my $third_column_found			= "false"; # columns in put file of file names used as new sample names
my $second_column_found			= "false"; # columns in put file of file names

my $GATK_fix_quals_string		= ""; # GATK string to deal with old style quality scores
my $zero_size_reported			= ""; # To make sure program only e-mails about zero file size once per loop
my $variant_effect_predictor	= ""; # snpEff or variant_effect_predictor
my $command						= "";
my $ref							= "";
my $ref_prefix					= "";
my $fastq_file_1				= ""; # first FASTQ file
my $fastq_file_2				= ""; # second FASTQ file
my $fastq_file_2_guess			= "";
my $mem							= ""; ## memory setting for Java
my $paired_ends_type			= "";
my $email_address				= "";
my $run_title					= "";
my $list_file					= "";
my $read_file_method			= "";
my $answer						= "";
my $single_line					= "";
my $results_all_folder			= "";  # Folder for all output files
my $ref_seq_name				= "";  # Name of reference sequence for pindel analysis
my $log_file					= "";
my $sample_name					= "";  # Name of each sample. Used to give the .vcf file an individual sample name
my $new_sample_name				= "";  # sample name given in input file (optional)
my $start_time					= "";
my $start_clean_time			= "";
my $start_vcf_time				= "";
my $run_time					= "";
my $end_time					= "";
my $current_time				= "";
my $make_raw_bam_time			= "";
my $make_raw_bam_time_total		= "";
my $make_clean_bam_time			= "";
my $make_clean_bam_time_total	= "";
my $total_time					= ""; # sum of raw bam, clean bam and VCF
my $make_vcf_time				= "";
my $start_loop_time				= "";
my $loop_time					= "";
my $last_current_time			= "";
my $stage_time					= ""; # Time for each stage of the pipeline
my $lib							= ""; #This is used by the ReadGroups section and i have made it let us record the reference sequence
my $running_input_string		= ""; # string of BAM file names for the HaplotypeCaller
my $species						= "";
my $char_1						= "";
my $char_2						= "";
my $back_to_one					= "";
my $prefix_name					= "";  # Name of bam_file ommitting the suffix .bam
my $file_to_find				= ""; 
my $current_directory			= "";

# Region strings
my $region						= "";
my $GATK_region_string			= ""; # String to be used with the -L option of GATK tools such as UnifiedGenotyper
my $GATK_chr_only_region_string 			= ""; # String to be used with Unified Genotyper (just narrows down to chromosome e.g. -L chr15)
my $GATK_known_sites_string		= ""; # String for known SNPs for BaseRecalibrator/CountCovariates
my $chromosome					= "";

#### File names FASTQ2BAM ####
my $aln_sa_sai					= "";
my $aln_sa1_sai					= "";
my $aln_sa2_sai					= "";
my $aligned_sam					= "";
my $aligned_sorted_rg_bam		= "";
my $aligned_sorted_rg_bai		= "";
my $command_log					= "";


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

my $final_f2b_bam					= ""; # Final BAM file after FASTQ2BAM pipeline has run
my $final_f2b_bai					= "";

#VCF files for SNPs and Indels
my $final_VCF_file				= "";
my $final_VCF_idx_file			= "";
my $snps_final_vcf				= "";

#LOG files
my $validate_f2b_out			= "";
my $flagstat_f2b_out			= "";
my $validate_final_out			= "";
my $flagstat_final_out			= "";

#Other files
my $insert_size_pdf				= ""; # Used by CollectInsertSizeMetrics
my $gc_bias_pdf					= ""; # Used by CollectGcBiasMetrics
my $dummy_dbsnp_file			= ""; # The dummy dbSNP file used by....
my $actual_dbsnp_file			= ""; # Correct dbSNP file for this chromosome (if it exists)
my $dbsnp_file					= ""; # dbsnp file used by GATK
my $readme_file					= "";
my $readme_overall_file			= "";

my $runtitle_sample_recal_grp 	= ""; # For BaseRealigner

my @bam_file_array				= ();
my @reads1_file_array			= ();
my @reads2_file_array			= ();
my @sample_name_array			= ();
my @file_name_array				= (); # to store list of FASTQ files in the DIR when entering single file
my @item						= ();


################
# START TIMERS #
################
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


# Warn if testing mode is on #
if ($testing_mode eq "on"){print color 'bold cyan';print"\n\nTESTING MODE ON\n\n";print color 'reset';}

print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      FASTQ2VCF       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";

print "  - This program processes FASTQ files into VCF files.\n\n";

print "    This new version puts all the output files into a single folder.\n\n";

print "    The 'f2b' BAM files with 'raw alignments' are deleted and the 'final'\n";
print "    BAM files after processing through the GATK stages are kept.\n\n\n";


print color 'bold cyan';

print "  - When you get your final VCF file you need to do two things:\n\n";
print "    1: Use run_variant_effect_predictor (or run_snpEff) to annotate your variants.\n\n";
print "    2: Use vcf2excel to prepare your VCF file for NGS SNP Viewer.\n\n";

print color 'reset';

#############################
# Name the analysis run     #
#############################

&print_message("Please enter a name for this analysis","input");

print "(Keep if fairly short as it becomes part of the file name - and with no spaces):\n\n";
$run_title = <STDIN>;
chomp $run_title;


#################################################################
# Check if overall results folder with this name exists already #
#################################################################

$results_all_folder = "results_fastq2vcf_"."$run_title";
	
if (-e $results_all_folder)
{
	&print_message("WARNING: A results folder called $results_all_folder already exists","warning");
	print "Choose a new run title or delete the existing folder.\n\n";
	exit;
}

$log_file = "$run_title"."_fastq2vcf_log.out";

$| = 1;

open(STDERR, "| tee $log_file");


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

#####################################
# Ask if the data is SE or PE data? #
#####################################
&print_message("Do you have a Single-end or Paired-end dataset?","input");
print "  <1>  Single-end\n";
print "  <2>  Paired-end\n\n";

$answer = <STDIN>;
chomp $answer;
if ($answer eq ""){$answer = "2"}

if (substr($answer,0,1) eq "1" ){$paired_ends_type = "SE"}
if (substr($answer,0,1) eq "2" ){$paired_ends_type = "PE"}


##################################
# Define data files              #
##################################

&print_message("Which reference sequence do you want to use?","input");

print "   <1> CanFam3\n";
print "   <2> CanFam3nu (unknown chromosomes removed)\n";
print "   <3> CanFam2\n";
print "   <4> EquCab2\n";
print "   <5> Human\n\n";

print "   <6> Strep. equi\n";
print "   <7> Strep. zoo\n\n";

print "   <9> other\n\n";

$answer = <STDIN>; chomp $answer;

if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3/canfam3_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/canfam3nu/canfam3nu.fasta"; $ref_seq_name = "canfam3nu"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3nu/canfam3_dummy_DBSNP.vcf"}
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

#########################################################
# Check if REF file exists (if a REF file is specified) #
########################################################

if (! -e "$ref")
{ 
	&print_message("File $ref does not exist","warning");
	print "  You might need to make a new_bwa version of the file.\n\n";
	exit;
}

$lib = $ref_seq_name; # for use in adding Read Group information to BAM files
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

print "   <1> mem - this is a new method available with the latest version of bwa (for reads over 70bp)\n";
print "   <2> aln - this is the method for old Illumina reads of less than 70bp\n\n";

$answer = <STDIN>;
chomp $answer;
if ($answer eq ""){$answer = "1"}; #default is mem

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

if (substr($answer,0,1) eq "1" ){$fix_misencoded_qual_scores = "no"; $GATK_fix_quals_string = ""}
if (substr($answer,0,1) eq "2" ){$fix_misencoded_qual_scores = "yes"; $GATK_fix_quals_string = "-fixMisencodedQuals "}


#########################################################
# NOTE fix quality scores when you get to the BAM files #
#########################################################


############################################################################
# Remove duplicates                                                        #
# Would you like to remove duplicates - noramlly YES but might not want to #
############################################################################
&print_message("Would you like to remove duplicate reads?","input");
print "(This is normally strongly advised, but for special purposes you might not want to)\n\n";

print "   <1> YES - remove duplicates\n";
print "   <2> NO - don't remove duplicates\n\n";

$answer = <STDIN>;
chomp $answer;

if ($answer eq ""){$answer = "1"} # default

if (substr($answer,0,1) eq "1"){$remove_duplicates = "yes"}
if (substr($answer,0,1) eq "2"){$remove_duplicates = "no"}


###########################################################################
# Region preferences                                                      #
# Would you like to focus the analysis on a specific region of the genome #
###########################################################################


&print_message("Would you like to focus on alignments to specific region of the genome?","input");
print "(e.g. a sequence-captured region of a larger genome)\n";
print color "bold white";
print "(This is strongly advised as it speeds up the whole process)\n\n";
print color "reset";
print "  NB: For strep_equi and strep_zoo genomes use 'chr1'\n\n";

print "   <1> YES - focus in on region\n";
print "   <2> NO - use whole genome\n\n";

$answer = <STDIN>;
chomp $answer;

if ($answer eq ""){$answer = "1"} # default

if (substr($answer,0,1) eq "1"){$use_defined_region = "yes"}
if (substr($answer,0,1) eq "2"){$use_defined_region = "no"}



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
	&print_message("Please define your region of interest (eg 'chr5:21000000-23000000')","input");
	$region = <STDIN>;
	chomp $region;
	
	###################################################
	# If only chromosome  is given (e.g. 15 or chr15) #
	###################################################
	if (index($region,":") == -1)
	{
		if (index($region,"chr") == -1){$region = "chr"."$region";}

		&print_message("If you know the chromosomal coordinates of your region it is much better to put them in","message");

		print "Adding the chromosomal coordinates will speed up the process and make smaller files, so you should do it.\n\n";

		&print_message("Please define your region of interest (eg 'chr5:21000000-23000000' or just 'chr5' if you have to)","input");
		$region = <STDIN>;
		chomp $region;

	} # if no colon
	
	
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
	$GATK_chr_only_region_string = "-L $chromosome";
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
			&print_message("A SNP file for this chromosome has been found","message");

			print "SNP file $actual_dbsnp_file exists\n";
			
			$dbsnp_file = $actual_dbsnp_file;
			$GATK_known_sites_string = " -knownSites $dbsnp_file";
		}
		
		# File not found
		if (!-e $actual_dbsnp_file)
		{
			&print_message("A SNP file for this chromosome has not been found","message");

			print "SNP file $actual_dbsnp_file doesn't exist              \n";
			print "Dummy SNP file $dummy_dbsnp_file will be used instead  \n";
			
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
			
		} # actual DBSNP not found
	}
} # if ($use_defined_region eq "yes")


###########################################################################
# Which variant caller do you want?                                       #
###########################################################################

&print_message("Which variant caller would you like to use?","input");

print "   <1> GATK Unified Genotyper (faster - this is now the default until HaplotypeCaller is sorted out)\n";
print "   <2> GATK Haplotype Caller  (newer - recommended by GATK, but 10x slower and can give spurious calls)\n\n";

$answer = <STDIN>; chomp $answer;

if (substr($answer,0,1) eq "1" ){$variant_caller = "UnifiedGenotyper"}
if (substr($answer,0,1) eq "2" ){$variant_caller = "HaplotypeCaller"}
if ($answer eq ""){$variant_caller = "UnifiedGenotyper"}


##################################################################################
# Filtering preferences                                                          #
# If you want to chnge the default values of stand_call_conf and stand_emit_conf #
##################################################################################


&print_message("Would you like to use the default values of stand_emit_conf and stand_call_conf?","input");

print "(These are GATK variant calling parameters that control quality thresholds)\n\n";

print "   <1>  Yes - use default values of $stand_emit_conf and $stand_call_conf\n";
print "   <2>  No  - enter my own values\n\n";

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



#####################################################################
# ASK IF YOU WANT TO READ THE FILENAMES FROM A "FILE OF FILE NAMES" #
#####################################################################

$answer = "";

until ($read_file_method eq "multiple" || $read_file_method eq "single")
{
		&print_message("How do you want to read the input files?","input");

	print "   <1> Using a file of file names for your FASTQ files\n";
	if ($paired_ends_type eq "SE"){print "   <2> Enter a single FASTQ file\n\n";}
	if ($paired_ends_type eq "PE"){print "   <2> Enter a pair of FASTQ files individually\n\n";}

	$answer = <STDIN>;chomp $answer;
	if (substr($answer,0,1) eq "1"){$read_file_method = "multiple"}
	if (substr($answer,0,1) eq "2"){$read_file_method = "single"}
}


#####################################
# Input file names for fastq files  #
#####################################

if ($read_file_method eq "single")
{
	if ($paired_ends_type eq "SE")
	{
		until (-e $fastq_file_1)
		{
			print "\nPlease input the name of your FASTQ sequence file (type 'ls' to see a list of fastq files)\n";
			print "(type 'ls' to see a list of text files):      ";

			$fastq_file_1 = <STDIN>;chomp $fastq_file_1;
			
			if ($fastq_file_1 eq "ls"){print "\n";system ("ls *.fastq")}
			
			if ($fastq_file_1 ne "ls")
			{
				if (! -e $fastq_file_1){print "\n\n>>>>>>>>  File $fastq_file_1 not found.  Try again.  <<<<<<<<\n\n";}
			}
		}

		# Assign to array of file names but just use first element #
		$reads1_file_array[1]=$fastq_file_1;
		$no_of_files = 1;
	}
	
	if ($paired_ends_type eq "PE")
	{
		#############################################
		# Get list of FASTQ files in this directory #
		#############################################
		$current_directory=getcwd;
		opendir (DIR, $current_directory);
		while (my $file_to_find = readdir(DIR)) 
		{
			if ((substr($file_to_find, -6) eq ".fastq") && !($file_to_find =~ m/^\./))
			{
				$file_count = $file_count + 1;
				$file_name_array[$file_count] = $file_to_find;
			}
		}
		$no_of_files_in_dir = $file_count;
		close DIR;

		until (-e $fastq_file_1)
		{
			print "\nPlease input the name of your 1st FASTQ sequence file (type 'ls' to see a list of fastq files):      ";
			$fastq_file_1 = <STDIN>;chomp $fastq_file_1;
			
			if ($fastq_file_1 eq "ls"){print "\n";system ("ls *.fastq")}
			
			if ($fastq_file_1 ne "ls")
			{
				if (! -e $fastq_file_1){print "\n\n>>>>>>>>  File $fastq_file_1 not found.  Try again.  <<<<<<<<\n\n";}
			}
		}
		
		# Guess second FASTQ file name #
		for ($check_count =1; $check_count <= $no_of_files_in_dir; $check_count++)
		{
			$fastq_file_2 = $file_name_array[$check_count];
			for ($pos = 0; $pos <= length ($fastq_file_1); $pos++)
			{
				$char_1 = substr ($fastq_file_1,$pos,2);
				if ($char_1 eq "_1")
				{
					# Does same position in second file have a "_2"?
					$char_2 = substr ($fastq_file_2,$pos,2);
					if ($char_2 eq "_2")
					{
						$back_to_one = substr($fastq_file_2,0,$pos)."_1".substr($fastq_file_2,$pos+2,99);

						if ($fastq_file_1 eq $back_to_one)
						{
							$fastq_file_2_guess = $fastq_file_2;
						}
					}
				}
			} 
		} # check_count loop to guess 2nd FASTQ file name

		$fastq_file_2 = "";

		until (-e $fastq_file_2)
		{
			print "\nPlease input the name of your 2nd FASTQ sequence file: (default=$fastq_file_2_guess)      ";
			$fastq_file_2 = <STDIN>;chomp $fastq_file_2;
			if ($fastq_file_2 eq ""){$fastq_file_2 = $fastq_file_2_guess}
			if ($fastq_file_2 eq "ls"){print "\n";system ("ls *.fastq")}
			
			if ($fastq_file_2 ne "ls")
			{
				if (! -e $fastq_file_2){print "\n\n>>>>>>>>  File $fastq_file_2 not found.  Try again.  <<<<<<<<\n\n";}
			}
		}
		
		# Assign to array of file names but just use first element #
		$reads1_file_array[1]=$fastq_file_1;
		$reads2_file_array[1]=$fastq_file_2;
		$no_of_files = 1;
	}
} #  if ($read_file_method eq "single")


if ($read_file_method eq "multiple")
{
	until (-e "$list_file")
	{
		if ($paired_ends_type eq "SE")
		{
			&print_message("Please input the name of your file with a list of file names of the FASTQ files","input");
		}
		if ($paired_ends_type eq "PE")
		{
			&print_message("Please input the name of your file with a list of file names of the paired FASTQ files","input");
			print "(The two paired-end 'reads' file names must be on the same line separated by a TAB):      ";
		}
		
		$list_file = <STDIN>;chomp $list_file;
		
		if ($list_file eq "ls"){print "\n";system ("ls *.txt"); print "\n";}
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
	$no_of_files=0;
	while ($single_line = <LIST> ) 
	{
		chomp $single_line;

		if (length $single_line  > 1) # so if user has left blank lines at the end of the file of file names it won't matter
		{
			if ($paired_ends_type eq "SE")
			{
				@item=split(/\t/,$single_line);
				
				$array_size = scalar @item;
				
				if ($array_size == 1)
				{
					$fastq_file_1 = $single_line;
					$reads1_file_array[$list_count]=$fastq_file_1;
				} # array size 1, SE
				
				if ($array_size == 2)
				{
					$fastq_file_1 = $item[0];
					$new_sample_name = $item[1];
					
					$reads1_file_array[$list_count]=$fastq_file_1;
					
					$sample_name_array[$list_count] = $new_sample_name;
					$second_column_found = "true";
					
				} # array size 2, SE
			} # SE
			
			if ($paired_ends_type eq "PE") # Read two columns separated by a TAB for the two files reads1 and reads2
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
				}# array_size = 1
				
				if ($array_size == 2)
				{
					$fastq_file_1 = $item[0];
					$fastq_file_2 = $item[1];
					
					$reads1_file_array[$list_count]=$fastq_file_1;
					$reads2_file_array[$list_count]=$fastq_file_2;
					$sample_name_array[$list_count] = "";
				}# array_size = 2
				
				########################################################
				# If there is a third column to use as the sample name #
				########################################################
				if ($array_size == 3)
				{
					$fastq_file_1 = $item[0];
					$fastq_file_2 = $item[1];
					$new_sample_name = $item[2];
					
					$reads1_file_array[$list_count] = $fastq_file_1;
					$reads2_file_array[$list_count] = $fastq_file_2;
					$sample_name_array[$list_count] = $new_sample_name;
					$third_column_found = "true";
				} # array_size = 3
			} # PE

			$list_count=$list_count + 1;

		} # if length $single_line > 1

	} # while loop

	close LIST;

	$no_of_files=$list_count - 1;
	
	
	########################
	# Check if files exist #
	########################
	
	&print_message("Checking input files","message");
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		if ($paired_ends_type eq "SE")
		{
			if (! -e "$reads1_file_array[$list_count]")
			{ 
				print "\nFile $reads1_file_array[$list_count] does not exist.\n\n";
				exit;
			} 
		}
		if ($paired_ends_type eq "PE")
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
	} # list count loop
	

	
	###################
	# List file names #
	###################
	print "\nThere are $no_of_files pairs of FASTQ 'reads' files in this file of file names.\n\n";
	
	if ($third_column_found eq "true")
	{
		print color "bold white";
		print "There is also a third column in the input file which will be used for renaming the files.\n\n";
		print color "reset";
	}
	if ($second_column_found eq "true")
	{
		print color "bold white";
		print "There is also a second column in the input file which will be used for renaming the files.\n\n";
		print color "reset";
	}
	
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		if ($paired_ends_type eq "SE")
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
		if ($paired_ends_type eq "PE")
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

	&print_message("If these file names look OK, press enter to proceed (or 'Q' to quit)","input");

	$answer = <STDIN>;chomp $answer; 
	if (lc $answer eq "q"){exit;} 

} # End of if ($read_file_method eq "multiple")

	
####################################
# Memory setting for java programs #
####################################
&print_message("Please enter the memory setting for this analysis","input");

print "Memory setting (default is -Xmx4g but you could use -Xmx8g if not too many others are using the server):      ";
$mem = <STDIN>;chomp $mem;
if ($mem eq ""){$mem = "-Xmx4g"}


####################################
# Number of CPU threads to use     #
####################################
&print_message("How many CPU threads would you like to use?","input");

print "Number of threads (you had better use 1 at the moment because we have had errors using larger numbers):      ";
$no_of_threads = <STDIN>;chomp $no_of_threads;
if ($no_of_threads eq ""){$no_of_threads = "1"}
if ($no_of_threads > 8){$no_of_threads = "8"}


if ($no_of_threads > 1)
{
	&print_message("You have chosen $no_of_threads threads.  OK, but I did warn you - you might get errors","warning");
	print "\n  >> Press return to continue\n\n";
	$answer=<STDIN>;
}

#################################
# Annotate variants preferences #
#################################

&print_message("Would you like to annotate your SNP and INDEL calls (using final VCF file)?","input");

print "   <1> Yes using variant effect predictor (recommended)\n";
print "   <2> Yes using snpEff\n";
print "   <3> No - I will annotate the VCF file later\n\n";

$answer = <STDIN>;chomp $answer;
if ($answer eq ""){$answer = "3"}

if (substr($answer,0,1) eq "1"){$annotate_variants = "yes"; $variant_effect_predictor = "vep"}
if (substr($answer,0,1) eq "2"){$annotate_variants = "yes"; $variant_effect_predictor = "snpEff"}
if (substr($answer,0,1) eq "3"){$annotate_variants = "no"; $variant_effect_predictor = "none"}

&print_message("Variant annotation is now done after this pipeline has run","message");

if ($variant_effect_predictor eq "snpEff" )
{
	print "  - When you get your final VCF file you need to do two things:\n\n";
	print "    1: Use run_snpEff to annotate your variants.\n\n";
	print "    2: Use vcf2excel to prepare your VCF file for NGS SNP Viewer (dog only).\n\n";
}
if ($variant_effect_predictor eq "none")
{
	print "  - When you get your final VCF file you need to do two things:\n\n";
	print "    1: Use run_variant_effect_predictor (or run_snpEff) to annotate your variants.\n\n";
	print "    2: Use vcf2excel to prepare your VCF file for NGS SNP Viewer (dog only).\n\n";
}

print "  >> Press return to continue";
$answer=<STDIN>;

$command =  "clear";
#system("$command");   TEMPORARY

&print_message("SUMMARY DETAILS - PLEASE CHECK CAREFULLY!!","message");

print "  Name of this analysis:         \t$run_title\n"; 
print "  Your email address:            \t$email_address\n";

&print_message("SETTINGS","message");

if ($paired_ends_type eq "SE")                 {print "  Single-end analysis            \tYES\n";}
if ($paired_ends_type eq "PE")                 {print "  Paired-end analysis            \tYES\n";}

print "  Memory setting:                \t$mem\n";
print "  No of CPU threads:             \t$no_of_threads\n";
print "  BWA alignment option:          \t$bwa_alignment_method\n";

if ($remove_duplicates eq "yes")   {print "  Remove duplicates              \tYES\n";}
if ($remove_duplicates eq "no")    {print "  Remove duplicates              \tNO\n";}

if ($annotate_variants eq "yes")   {print "  Annotate variants              \tYES using $variant_effect_predictor\n";}
if ($annotate_variants eq "no")    {print "  Annotate variants              \tNO\n";}


print "  Variant caller:                 \tGATK $variant_caller\n\n";
print "  GATK variant caller constants:\n\n";

print "    --stand_emit_conf:\t\t\t$stand_emit_conf\n";
print "    --stand_call_conf:\t\t\t$stand_call_conf\n";
print "    --max_alt_alleles:\t\t\t$max_alt_alleles\n";

&print_message("DATA FILES","message");

print "  Reference sequence:            \t$ref\n\n";

if ($use_defined_region eq "yes")
{
	print "  Defined region:                  \t$region\n\n";
	print "  GATK region string:              \t$GATK_region_string\n";
}
if ($use_defined_region eq "no")
{
	print "Defined region:\t\t\t\tNO\n\n";
}

print "  SNP file for BaseRecalibrator:   \t$dbsnp_file\n\n";

print "  Number of FASTQ files to analyse:\t$no_of_files\n\n";
	
for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	if ($paired_ends_type eq "SE")
	{
		print "    FASTQ file $list_count	\t$reads1_file_array[$list_count]\n";
	}
	if ($paired_ends_type eq "PE")
	{
		print "    Pair of FASTQ files $list_count	\t$reads1_file_array[$list_count]\t$reads2_file_array[$list_count]\n";
	}
}

&print_message("  >> CHECK THESE SETTINGS.  If OK press ENTER to proceed with the analysis run (or Q to quit)","message");

$answer = <STDIN>;chomp $answer; 
if (lc $answer eq "q"){exit;} 
	

#########################
# open Command Log file #
#########################
$command_log = "$run_title"."_fastq2vcf_command_log.out";
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

print COMMAND_LOG "COMMAND LOG for fastq2vcf version $version\n\n";
print COMMAND_LOG "Analysis name:           \t$run_title\n";
print COMMAND_LOG "Reference sequence:      \t$ref\n";
print COMMAND_LOG "Ref seq name:            \t$ref_seq_name\n\n";

print COMMAND_LOG "Memory setting:          \t$mem\n";
print COMMAND_LOG "No of CPU threads:       \t$no_of_threads\n";
print COMMAND_LOG "BWA alignment option:    \t$bwa_alignment_method\n\n";

if ($remove_duplicates eq "yes")   {print COMMAND_LOG "Remove duplicates        \tYES\n\n";}
if ($remove_duplicates eq "no")    {print COMMAND_LOG "Remove duplicates        \tNO\n\n";}

if ($use_defined_region	eq "yes")
{
	print COMMAND_LOG "Defined region:          \t$region\n\n";
	print COMMAND_LOG "  GATK region string:    \t$GATK_region_string\n";
}
if ($use_defined_region	eq "no"){print COMMAND_LOG "No genome region defined\n\n";}

print COMMAND_LOG "Variant caller:           \tGATK $variant_caller\n\n";
print COMMAND_LOG "GATK variant caller constants:\n\n";

print COMMAND_LOG "  --stand_emit_conf:\t\t\t$stand_emit_conf\n";
print COMMAND_LOG "  --stand_call_conf:\t\t\t$stand_call_conf\n";
print COMMAND_LOG "  --max_alt_alleles:\t\t\t$max_alt_alleles\n\n";

print COMMAND_LOG "Fix quality scores:      \t$fix_misencoded_qual_scores\n\n";
print COMMAND_LOG "SNP file for BaseRecalibrator: \t$dbsnp_file\n\n";

print COMMAND_LOG "List of input files\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	if ($paired_ends_type eq "SE")
	{
		print COMMAND_LOG "FASTQ FILE $list_count	\t$reads1_file_array[$list_count]\n";
	}
	if ($paired_ends_type eq "PE")
	{
		print COMMAND_LOG "FASTQ FILES $list_count	\t$reads1_file_array[$list_count]\t$reads2_file_array[$list_count]\n";
	}
}

$final_VCF_file = "$run_title"."_variants.vcf";

print COMMAND_LOG "\n\n";
print COMMAND_LOG "Output VCF file:         \t$final_VCF_file\n\n";

########################################
# Create results_all folder            #
# This holds the files from fastq2bam  #
# and the bam2vcf folder(s)?           #
########################################
&run_unix_command("mkdir $results_all_folder","Make results folder");	


##############################################################################
##############################################################################
# Start MAIN LOOP.  This runs through each of the files in the list of files #
##############################################################################
##############################################################################

$last_current_time = time();

&log_time($last_current_time,"Current time logged");


for ($loop_count=1;$loop_count <=$no_of_files;$loop_count++)
{

	#################################################
	# Record elapsed time at the start of each loop #
	#################################################
	$start_loop_time = time();

	$no_of_steps=18; # fastq2vcf
	
	######################################################
	# Get names for FASTQ reads files                    #
	######################################################
	$fastq_file_1 = $reads1_file_array[$loop_count];
	$fastq_file_2 = $reads2_file_array[$loop_count];
	
	
	######################################################
	# Create sample name (to use to name BAM files etc)  #
	######################################################
	$sample_name = &get_prefix ($fastq_file_1);

	
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
	
	&print_message("Loop count: $loop_count  Sample name: $sample_name","message");
	

	#######################################################################
	# Set up various File names which vary for each loop  from FASTQ2BAM  #
	#######################################################################
	$aln_sa_sai = "$run_title"."_"."$sample_name"."_aln_sa.sai";
	$aln_sa1_sai = "$run_title"."_"."$sample_name"."_aln_sa1.sai";
	$aln_sa2_sai = "$run_title"."_"."$sample_name"."_aln_sa2.sai";
	$aligned_sam = "$run_title"."_"."$sample_name"."_aligned.sam";
	$aligned_sorted_rg_bam = "$run_title"."_"."$sample_name"."_aligned_sorted_rg.bam";
	$aligned_sorted_rg_bai = "$run_title"."_"."$sample_name"."_aligned_sorted_rg.bai";
	
	$final_f2b_bam = "$run_title"."_"."$sample_name"."_f2b.bam"; # Final BAM of FASTQ2BAM
	$final_f2b_bai = "$run_title"."_"."$sample_name"."_f2b.bai";
	
	$validate_f2b_out = "$run_title"."_"."$sample_name"."_validate_f2b.out";
	$flagstat_f2b_out = "$run_title"."_"."$sample_name"."_flagstat_f2b.out";
	

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
	
	# After BaseRecalibrator
	$runtitle_sample_recal_csv  = "$run_title"."_"."$sample_name".".recal.csv";
	$runtitle_sample_recal_grp  = "$run_title"."_"."$sample_name".".recal.grp"; # Use for BaseRecalibrator
	
	# After TableRecalibration (recal)
	$runtitle_clean_dedup_recal_bam =  "$run_title"."_"."$sample_name".".clean.dedup.recal.bam";
	$runtitle_clean_dedup_recal_bai =  "$run_title"."_"."$sample_name".".clean.dedup.recal.bai";
	
	# After reducing BAM to a region
	$region_only_bam = "$run_title"."_"."$sample_name"."_region.bam";
	$region_only_bai = "$run_title"."_"."$sample_name"."_region.bai";
	
	# After all stages have finished, the BAM file is renamed to final_bam
	$final_bam = "$run_title"."_"."$sample_name"."_final.bam";
	$final_bai = "$run_title"."_"."$sample_name"."_final.bai";
	
	# VCF files
	$final_VCF_file = "$run_title"."_variants.vcf";
	$final_VCF_idx_file = "$run_title"."_variants.vcf.idx";
	
	#Other files
	$insert_size_pdf = "$sample_name"."_insert_size.pdf";
	$gc_bias_pdf = "$sample_name"."_gc_bias.pdf";
	$validate_final_out = "$run_title"."_"."$sample_name"."_validate_final.out";
	$flagstat_final_out = "$run_title"."_"."$sample_name"."_flagstat_final.out";
	
	#file names
	
	$runtitle_sample_metrics = "$run_title"."_"."$sample_name"."_MarkDuplicates_metrics.out";


	####################################################################
	# 1 Run bwa                                                        #
	####################################################################

	#####################
	# SINGLE END DATA ! #
	#####################

	if  ($paired_ends_type eq "SE")
	{

		if ($bwa_alignment_method eq "aln")
		{

			&print_message("File $loop_count/$no_of_files   Step 1 Running BWA to produce aligned BAM file - SE","message");

			&run_unix_command("$bwa_path $bwa_alignment_method $ref $fastq_file_1 > $aln_sa_sai","1 Run bwa to make SAI file");

			&record_input_file_size("$fastq_file_1");
			&record_output_file_size("$aln_sa_sai");

			&print_message("File $loop_count/$no_of_files   Step 1 COMPLETED Running BWA to produce aligned BAM file - SE","message");

			&test_mode_subroutine;

		} # bwa aln

		if ($bwa_alignment_method eq "mem")
		{

			&print_message("File $loop_count/$no_of_files   Step 1 Running BWA mem to produce aligned SAM file - SE","message");

			&run_unix_command("$bwa_path mem -M -t $no_of_threads $ref $fastq_file_1 > $aligned_sam","Converting FASTQ to SAM using bwa mem");


			&record_input_file_size("$fastq_file_1");
			&record_output_file_size("$aligned_sam");

			&print_message("File $loop_count/$no_of_files   Step 1 COMPLETED Running BWA to produce aligned SAM file - SE","message");

			&test_mode_subroutine;

		} # bwa mem

	} # SE


	#####################
	# PAIRED END DATA ! #
	#####################

	if  ($paired_ends_type eq "PE")
	{
	
		if ($bwa_alignment_method eq "aln")
		{

			##########
			# Read 1 #
			##########
			
			&print_message("File $loop_count/$no_of_files   Step 1a Running BWA ALN to produce aligned BAM file - PE1","message");

			&run_unix_command("$bwa_path aln -t $no_of_threads $ref $fastq_file_1 > $aln_sa1_sai","1a Run bwa to make SAI file 1");

			&record_input_file_size("$fastq_file_1");
			&record_output_file_size("$aln_sa1_sai");
			
			&print_message("File $loop_count/$no_of_files   Step 1a COMPLETED Running BWA ALN to produce aligned BAM file - PE1","message");

			&test_mode_subroutine;
			
			
			##########
			# Read 2 #
			##########
			
			&print_message("File $loop_count/$no_of_files   Step 1b Running BWA ALN to produce aligned BAM file - PE2","message");
			
			&run_unix_command("$bwa_path aln -t $no_of_threads $ref $fastq_file_2 > $aln_sa2_sai","1b Run bwa to make SAI file 2");
			
			&record_input_file_size("$fastq_file_2");
			&record_output_file_size("$aln_sa2_sai");

			&print_message("File $loop_count/$no_of_files   Step 1b COMPLETED Running BWA ALN to produce aligned BAM file - PE2","message");

			&test_mode_subroutine;

		} #bwa aln

		if ($bwa_alignment_method eq "mem")
		{

			################################################
			# Read 1 and Read 2 are done together with mem #
			################################################
			
			&print_message("File $loop_count/$no_of_files   Step 1 Running BWA MEM to produce SAM file - PE","message");
			
			&run_unix_command("$bwa_path mem -M -t $no_of_threads $ref $fastq_file_1 $fastq_file_2 > $aligned_sam","1 Converting FASTQ to SAM using bwa mem");

			&record_input_file_size("$fastq_file_1");
			&record_input_file_size("$fastq_file_2");
			&record_output_file_size("$aligned_sam");
			
			&print_message("File $loop_count/$no_of_files   Step 1 COMPLETED Running BWA MEM to produce SAM file - PE","message");

			&test_mode_subroutine;
			
		} #bwa mem

	} # PE

	&problem_check;


	##########################################################################
	# 2 Convert SAI files into SAM files  (aln only, not required for mem)   #
	##########################################################################
	if ($bwa_alignment_method eq "aln")
	{
		if  ($paired_ends_type eq "SE")
		{
			&print_message("File $loop_count/$no_of_files   Step 2 Using SAI file to create SAM file","message");

			&run_unix_command("$bwa_path samse $ref $aln_sa_sai $fastq_file_1 > $aligned_sam","2 Converting SAI to SAM");

			&record_input_file_size("$aln_sa_sai");
			&record_output_file_size("$aligned_sam");

			&test_mode_subroutine;
		
			&problem_check;
		
			&print_message("File $loop_count/$no_of_files  Step 2 COMPLETED Using SAI file to create SAM file","message");
		}

		if  ($paired_ends_type eq "PE")
		{
			&print_message("File $loop_count/$no_of_files   Step 2 Using SAI files to create SAM file","message");

			&run_unix_command("$bwa_path sampe $ref $aln_sa1_sai $aln_sa2_sai $fastq_file_1 $fastq_file_2 > $aligned_sam","2 Converting SAI to SAM");

			&record_input_file_size("$aln_sa1_sai");
			&record_input_file_size("$aln_sa2_sai");
			&record_output_file_size("$aligned_sam");
			
			&test_mode_subroutine;
		
			&problem_check;
		
			&print_message("File $loop_count/$no_of_files   Step 2 COMPLETED Using SAI files to create SAM file","message");
		}

	} # aln only

    print " * * * \n";



	##########################################################
	# FASTQ2BAM 3                                             #
	# Run all the picard stages in one                       #
	#  - Sort the SAM file generated by BWA                  #
	#  - Convert the SAM file (from BWA) to a BAM file       #
	#  - Add Read Group information to the BAM file          #	
	##########################################################

	&print_message("File $loop_count/$no_of_files   Step 3 Picard stages","message");
	print "  - ReadGroup info to be added to the BAM file\n";
	print "  - The aligned reads to be sorted\n";
	print "  - The file to be converted to BAM format\n\n";

	&run_unix_command("java $mem -jar /opt/picard/AddOrReplaceReadGroups.jar I=$aligned_sam O=$aligned_sorted_rg_bam rgID=$run_title LB=$lib PL='ILLUMINA' PU=$sample_name SM=$sample_name SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=$picard_validation_stringency","3 Picard stages");	
	
	&record_input_file_size("$aligned_sam");
	&record_output_file_size("$aligned_sorted_rg_bam");

	&problem_check;
	
	&print_message("File $loop_count/$no_of_files   Step 3 Picard stages COMPLETED","message");
	
	
	################################################
	# FASTQ2BAM 4: Use picard ValidateSamFile      #
	################################################
	&print_message("File $loop_count/$no_of_files   Step 4 Picard ValidateSamFile","message");
	
	&run_unix_command("java $mem -jar /opt/picard/ValidateSamFile.jar I=$aligned_sorted_rg_bam O=$validate_f2b_out MODE=SUMMARY MAX_OPEN_TEMP_FILES=$max_open_temp_files VALIDATION_STRINGENCY=SILENT","4 Picard ValidateSamFile");

	&record_input_file_size("$aligned_sorted_rg_bam");
	&record_output_file_size("$validate_f2b_out");

	&print_message("File $loop_count/$no_of_files   Step 4 Picard ValidateSamFile COMPLETED","message");

	
	################################################
	# FASTQ2VCF: 5 Use samtools flagstat           #
	################################################
	&print_message("File $loop_count/$no_of_files   Step 5 samtools flagstat","message");
	
	&run_unix_command("/opt/samtools/samtools flagstat $aligned_sorted_rg_bam > $flagstat_f2b_out","5 samtools flagsta");

	&record_input_file_size("$aligned_sorted_rg_bam");
	&record_output_file_size("$flagstat_f2b_out");

	&print_message("File $loop_count/$no_of_files   Step 5 samtools flagstat COMPLETED","message");

	

	#############################################
	# FASTQ2VCF                                 #
	# Rename as final_f2b.bam  (for next stage) #
	#############################################

	&record_input_file_size("$aligned_sorted_rg_bam");
	&rename_file("$aligned_sorted_rg_bam","$final_f2b_bam");

	&record_input_file_size("$aligned_sorted_rg_bai");
	&run_unix_command("mv $aligned_sorted_rg_bai $final_f2b_bai","Change name to final_f2b_bai");

	&record_output_file_size("$final_f2b_bam");
	&record_output_file_size("$final_f2b_bai");
	
	print "\n\n";
	
	##########
	# README #
	##########

	$readme_file =  "$run_title"."_"."$sample_name"."_f2b_readme.out";

	open (READMEFILE, ">$readme_file"); 

	print READMEFILE "#######################################################\n";
	print READMEFILE "Summary of fastq2vcf results files for this sample     \n";
	print READMEFILE "#######################################################\n\n";

	print READMEFILE "Program:\t\t\tfastq2vcf (fastq2bam section)\tVersion: $version\n\n";
		
	print READMEFILE "Run title:\t\t\t$run_title\n\n";
	
	print READMEFILE "Reference sequence:\t$ref_seq_name\n\n";
	
	print READMEFILE "FASTQ file:        \t$fastq_file_1\n\n\n";

	print READMEFILE "BAM ALIGNMENT FILE (Raw alignments to the reference just using bwa)\n\n";

	print READMEFILE "\t$final_f2b_bam   (Note: f2b means this BAM file comes from the fastq2bam stage)\n\n\n";
	
	print READMEFILE "\tThis BAM file has been deleted as an intermediate file, but to check if it was OK look at these two files:\n\n";
	
	print READMEFILE "\t\t$validate_f2b_out\n\n";
	print READMEFILE "\t\t$flagstat_f2b_out\n\n\n";
	
	print READMEFILE "This BAM file then went through the rest of the NGS pipeline using the GATK pipeline\n\n";

	print READMEFILE "FINAL BAM ALIGNMENT FILES AFTER PROCESSING BY GATK\n\n";

	print READMEFILE "\t$final_bam		     \tBest alignments to the reference after processing by GATK\n\n\n";
	
	print READMEFILE "\tLook at these files to check the quality of your final BAM file:\n\n";

	print READMEFILE "\t\t$validate_final_out\n\n";
	print READMEFILE "\t\t$flagstat_final_out\n\n\n";

	close (READMEFILE);

	&move_to_results_all_folder ("$readme_file");

	#########################################################################################

	##############################
	# Set up bam_file to be used #
	##############################

	$bam_file = $final_f2b_bam;
	$bai_file = $final_f2b_bai;
		
	print "\n* * * \n";

	###########################################################################################
	# FASTQ2VCF 6: Use picard/AddOrReplaceReadGroups to update ReadGroup info in BAM file     #
	###########################################################################################


	&print_message("File $loop_count/$no_of_files   Step 6  Using picard/AddOrReplaceReadGroups to update ReadGroup info in BAM file","message");

	&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/AddOrReplaceReadGroups.jar I=$bam_file O=$rg_bam SO=coordinate VALIDATION_STRINGENCY=LENIENT CREATE_INDEX=true RGID=$sample_name  RGLB=$lib  RGPL=illumina  RGPU=$sample_name  RGSM=$sample_name","6 Add ReadGroups information");

	&record_input_file_size ("$bam_file");
	&record_output_file_size ("$rg_bam");

	&print_message("File $loop_count/$no_of_files   Step 6  COMPLETE Using picard/AddOrReplaceReadGroups to update ReadGroup info in BAM file","message");

	print "\n* * * \n";



	############################################
	# Record time for making the raw BAM files #
	############################################
	$current_time = time();
	$make_raw_bam_time = $current_time - $start_loop_time;

	$make_raw_bam_time_total = $make_raw_bam_time_total + $make_raw_bam_time;

print COMMAND_LOG "start_loop_time:         \t$start_loop_time\n";
print COMMAND_LOG "current_time:            \t$current_time\n";
print COMMAND_LOG "make_raw_bam_time:       \t$make_raw_bam_time\n";
print COMMAND_LOG "make_raw_bam_time_total: \t$make_raw_bam_time_total\n";

&log_time($start_loop_time,"Start loop $loop_count time");
&log_time($current_time,"Current time");
&log_time($make_raw_bam_time,"Time to make raw bam file $loop_count");
&log_time($make_raw_bam_time_total,"Total time to make raw BAM files");


	##########################################################
	# Record elapsed time at start of making clean BAM files #
	##########################################################
	$start_clean_time = time();

	########################################################################################
	# Now bam_file from fastq2bam has been used, it can be moved to the results_all_folder #
	########################################################################################
	&move_to_results_all_folder ("$bam_file");
	&move_to_results_all_folder ("$bai_file");


	########################################################################################
	# FASTQ2VCF 7: Use picard/MarkDuplicates to mark duplicates                            #
	########################################################################################
	if ($remove_duplicates eq "yes")
	{

		&print_message("File $loop_count/$no_of_files   Step 7  Use picard/MarkDuplicates to mark duplicates","message");

		&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/MarkDuplicates.jar I=$rg_bam O=$runtitle_dedup_bam  VALIDATION_STRINGENCY=$picard_validation_stringency CREATE_INDEX=true M=$runtitle_sample_metrics","7 Mark Duplicates");
		&record_input_file_size ("$rg_bam");
		&record_output_file_size ("$runtitle_dedup_bam");
		
		&print_message("File $loop_count/$no_of_files   COMPLETED Step 7  Use picard/MarkDuplicates to mark duplicates","message");
	}

	if ($remove_duplicates eq "no")
	{
		$runtitle_dedup_bam = $rg_bam;
		&record_output_file_size ("$runtitle_dedup_bam");
	}

	print "\n* * * \n";


	##############################################################################################################
	# FASTQ2VCF Step 8: Use GenomeAnalysisTK/RealignerTargetCreator to identify intervals in need of realignment #
	##############################################################################################################
	&print_message("File $loop_count/$no_of_files   Step 8:  Use GenomeAnalysisTK/RealignerTargetCreator to identify intervals in need of realignment","message");
	&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T RealignerTargetCreator $GATK_region_string -R $ref -I $runtitle_dedup_bam -o $runtitle_dedup_bam_intervals -mismatch 0.0 -S $GATK_validation_stringency $GATK_fix_quals_string","8 RealignerTargetCreator");
	&record_input_file_size ("$runtitle_dedup_bam");
	&record_output_file_size ("$runtitle_dedup_bam_intervals");

	&print_message("File $loop_count/$no_of_files   COMPLETED Step 8:  Use GenomeAnalysisTK/RealignerTargetCreator to identify intervals in need of realignment","message");

	print "\n* * * \n";
		

	##########################################################################################################
	# FASTQ2VCF Step 9:  Use GenomeAnalysisTK/IndelRealigner to Realign Indels                                 #
	##########################################################################################################
	&print_message("File $loop_count/$no_of_files   Step 9: Using GenomeAnalysisTK/IndelRealigner to Realign Indels","message");
	&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T IndelRealigner $GATK_region_string -R $ref -I $runtitle_dedup_bam -targetIntervals $runtitle_dedup_bam_intervals -o $runtitle_clean_dedup_bam -S $GATK_validation_stringency $GATK_fix_quals_string","9 IndelRealigner");

	&record_input_file_size ("$runtitle_dedup_bam");
	&record_input_file_size ("$runtitle_dedup_bam_intervals");
	&record_output_file_size ("$runtitle_clean_dedup_bam");

	&print_message("File $loop_count/$no_of_files   COMPLETED Step 9: Using GenomeAnalysisTK/IndelRealigner to Realign Indels","message");

	print "\n* * * \n";

	##############################################
	# Record that quality scores have been fixed #
	##############################################
	if ($GATK_fix_quals_string ne "")
	{
		print COMMAND_LOG "\n##########################################\n";
		print COMMAND_LOG "# IMPORTANT MESSAGE ABOUT QUALITY SCORES #\n";
		print COMMAND_LOG "##########################################\n\n";
		print COMMAND_LOG "Because you stated that the quality scores were of the old style\n";
		print COMMAND_LOG "the GATK RealignerTargetCreator and the IndelRealigner have been run using $GATK_fix_quals_string\n\n";

		print COMMAND_LOG "NB: If the scores were really new style, then this stage will not work properly (and may have crashed)\n\n";
	}

	##########################################################################################################
	# FASTQ2VCF Step 10:  NEW: Use GenomeAnalysisTK/BaseRecalibrator to generate a recalibration table       #
	##########################################################################################################
	&print_message("File $loop_count/$no_of_files   Step 10: Using GenomeAnalysisTK/BaseRecalibrator to generate a Recalibration Table","message");


	&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T BaseRecalibrator $GATK_region_string -R $ref -I $runtitle_clean_dedup_bam $GATK_known_sites_string -o $runtitle_sample_recal_grp","10 BaseRecalibrator");

	&record_input_file_size ("$runtitle_clean_dedup_bam");
	&record_output_file_size ("$runtitle_sample_recal_grp");

	&print_message("File $loop_count/$no_of_files   COMPLETED Step 10: Using GenomeAnalysisTK/BaseRecalibrator to generate a Recalibration Table","message");

	print "\n* * * \n";


	##########################################################################################################
	# FASTQ2VCF Step 11:  Use GenomeAnalysisTK/PrintReads to update the base quality scores                  #
	##########################################################################################################
	&print_message("File $loop_count/$no_of_files   Step 11:  Using GenomeAnalysisTK/PrintReads to update the base quality scores","message");

	&run_unix_command ("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T PrintReads $GATK_region_string -I $runtitle_clean_dedup_bam  -R $ref  -BQSR $runtitle_sample_recal_grp  -o $runtitle_clean_dedup_recal_bam -S $GATK_validation_stringency","11 PrintReads");

	&record_input_file_size ("$runtitle_clean_dedup_bam");
	&record_output_file_size ("$runtitle_clean_dedup_recal_bam");

	&print_message("File $loop_count/$no_of_files   COMPLETED Step 11:  Using GenomeAnalysisTK/PrintReads to update the base quality scores","message");


	###########################################################################################################
	# FASTQ2VCF Step 12: Use picard/ValidateSamFile to validate the cleaned BAM file (Summary only)             #
	###########################################################################################################

	if ($validate_bam_files eq "yes")
	{
		&print_message("File $loop_count/$no_of_files   Step 12:  Using picard/ValidateSamFile to validate the cleaned BAM file","message");

		&run_unix_command ("java $mem $temp_dir_string -jar /opt/picard/ValidateSamFile.jar INPUT=$runtitle_clean_dedup_recal_bam OUTPUT=$validate_final_out VALIDATION_STRINGENCY=$picard_validation_stringency MODE=SUMMARY MAX_OUTPUT=100 MAX_OPEN_TEMP_FILES=$max_open_temp_files","12 Picard ValidateSamFile");

		&record_input_file_size ("$runtitle_clean_dedup_recal_bam");
		&record_output_file_size ("$validate_final_out");
		
		&move_to_results_all_folder ("$validate_final_out");
		
		&print_message("File $loop_count/$no_of_files   COMPLETED Step 12: Using picard/ValidateSamFile to validate the cleaned BAM file","message");

		print "\n* * * \n";	


		#########################################################################################################
		# FASTQ2VCF Step 13: Use samtools flagstat to get simple stats on bam file                              #
		#########################################################################################################
		&print_message("File $loop_count/$no_of_files   Step 13: samtools flagstat to be carried out","message");

		&run_unix_command("/opt/samtools/samtools flagstat $runtitle_clean_dedup_recal_bam > $flagstat_final_out","13 samtools flagstat");

		&record_input_file_size ("$runtitle_clean_dedup_recal_bam");
		&record_output_file_size("$flagstat_final_out");

		&move_to_results_all_folder ("$flagstat_final_out");
			
		&print_message("File $loop_count/$no_of_files   COMPLETED Step 13: samtools flagstat to be carried out","message");

	} # if validate_bam_files

	print "\n* * * \n";

	########################################################################################################################
	# FASTQ2VCF Step 14: If you are using a defined region then create a new smaller BAM file with only this region in it  #
	########################################################################################################################
	if ($use_defined_region eq "yes")
	{

		&print_message("File $loop_count/$no_of_files   Step 14:  Making Bam file for region $region","message");

		&run_unix_command("/opt/samtools/samtools view $runtitle_clean_dedup_recal_bam $region -b -o $region_only_bam","14 Make smaller BAM file to specified region");

		&record_input_file_size ("$runtitle_clean_dedup_recal_bam");
		&record_output_file_size ("$region_only_bam");	
		
		&print_message("File $loop_count/$no_of_files   COMPLETED Step 14:  Making Bam file for region $region","message");
		
		############################################################################
		# FASTQ2VCF Step 15:  Now make an index file for this new smaller BAM file #
		############################################################################
		&print_message("File $loop_count/$no_of_files   Step 15:  New Index for region-only BAM file being created","message");

		&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/BuildBamIndex.jar I=$region_only_bam O=$region_only_bai VALIDATION_STRINGENCY=LENIENT","15 Make BAM index");

		&record_input_file_size ("$runtitle_clean_dedup_recal_bam");
		&record_output_file_size ("$region_only_bai");
		
		&print_message("File $loop_count/$no_of_files   COMPLETED Step 15:  New Index for region-only BAM file being created","message");
		
	} # if ($use_defined_region eq "yes"}


	##############################################################
	# If you didn't use a defined region, then change the final  #
	# BAM file from runtitle_clean_dedup_recal_bam to final_bam  #
	##############################################################
	if ($use_defined_region eq "no")
	{
		&record_input_file_size ("$runtitle_clean_dedup_recal_bam");
		&run_unix_command("mv $runtitle_clean_dedup_recal_bam $final_bam","bam2vcf Change name: $runtitle_clean_dedup_recal_bam to $final_bam");
		&record_output_file_size ("$final_bam");
		
		&record_input_file_size ("$runtitle_clean_dedup_recal_bai");
		&run_unix_command("mv $runtitle_clean_dedup_recal_bai $final_bai","bam2vcf Change name: $runtitle_clean_dedup_recal_bai to $final_bai");
		&record_output_file_size ("$final_bai");
		
	}


	##############################################################
	# If you DID use a defined region, then change the final     #
	# BAM file from region_only_bam to final_bam                 #
	##############################################################
	if ($use_defined_region eq "yes")
	{
		&record_input_file_size ("$region_only_bam");
		&run_unix_command("mv $region_only_bam $final_bam","bam2vcf Change name: $region_only_bam to $final_bam");
		&record_output_file_size ("$final_bam");
		
		&record_input_file_size ("$region_only_bai");
		&run_unix_command("mv $region_only_bai $final_bai","bam2vcf Change name: $region_only_bai to $final_bai");
		&record_output_file_size ("$final_bai");
	}

	print "\n* * * \n";


	##############################################################################################
	# Create input string for parallel use of the UnifiedGenotyper or HaplotypeCaller at the end #
	##############################################################################################

	if (-e $final_bam)
	{
		$running_input_string = $running_input_string." -I $results_all_folder/$final_bam";
	}
	if (! -e $final_bam)
	{
		$bam_missing_for_variant_calling = $bam_missing_for_variant_calling + 1;
	}


	########################################
	# FASTQ2VCF Step 16: Calculate GC bias #
	########################################
	&print_message("File $loop_count/$no_of_files   Step 16: Calculating GC Bias....","message");

	&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/CollectGcBiasMetrics.jar REFERENCE_SEQUENCE=$ref INPUT=$final_bam O=out1.junk CHART_OUTPUT=$gc_bias_pdf VALIDATION_STRINGENCY=LENIENT","16 calculate GC bias");

	&record_input_file_size ("$final_bam");
	&record_output_file_size("$gc_bias_pdf");

	&print_message("File $loop_count/$no_of_files   Step 16: GC bias has been calculated","message");


	##################################################
	# FASTQ2VCF Step 17:  Plot Insert Size Histogram #
	##################################################
	&print_message("File $loop_count/$no_of_files   Step 17: Plotting insert size histogram...","message");

	&run_unix_command("java $mem $temp_dir_string -jar /opt/picard/CollectInsertSizeMetrics.jar INPUT=$final_bam O=out2.junk HISTOGRAM_FILE=$insert_size_pdf VALIDATION_STRINGENCY=LENIENT","17 Plot insert size");

	&record_input_file_size ("$final_bam");
	&record_output_file_size_no_mail("$insert_size_pdf");

	&print_message("File $loop_count/$no_of_files   Step 17: Insert size histogram plotted","message");

	&print_message("Moving various files","message");

	
	##########################################
	# FASTQ2BAM                              #
	# Move various files to a Results Folder #
	##########################################
	
	&move_to_results_all_folder ("$validate_f2b_out");
	&move_to_results_all_folder ("$flagstat_f2b_out");
	

	##########################################
	# BAM2VCF                                #
	# Move various files to a Results Folder #
	##########################################
	
	&move_to_results_all_folder ("$final_bam");
	&move_to_results_all_folder ("$final_bai");
	&move_to_results_all_folder ("$runtitle_sample_metrics"); # MarkDuplicates info
	
	&move_to_results_all_folder ("$gc_bias_pdf");
	&move_to_results_all_folder ("$insert_size_pdf");
	

	##############################################################
	# For first loop copy over all files (for checking purposes) #
	# (only if $keep_all_for_first_file = "yes")                 #
	##############################################################
	if (($loop_count == 1) && ($keep_all_for_first_file eq "yes"))
	{
		&move_to_results_all_folder ("$rg_bam");
		&move_to_results_all_folder ("$rg_bai");
		&move_to_results_all_folder ("$runtitle_dedup_bam");
		&move_to_results_all_folder ("$runtitle_dedup_bai");
		&move_to_results_all_folder ("$runtitle_clean_dedup_bam");
		&move_to_results_all_folder ("$runtitle_clean_dedup_bai");
		&move_to_results_all_folder ("$runtitle_clean_dedup_recal_bam");
		&move_to_results_all_folder ("$runtitle_clean_dedup_recal_bai");
		#&move_to_results_all_folder ("$region_only_bam");
		#&move_to_results_all_folder ("$region_only_bai");
		&move_to_results_all_folder("$runtitle_dedup_bam_intervals");
	}
	
	########################################
	# Some files don't need to be kept     #
	# (keep files for first loop to check) #
	########################################
	
	if ($loop_count == 1){$delete_intermediate_files = "no"}else{$delete_intermediate_files = "yes"}
	if ($keep_all_for_first_file eq "no"){$delete_intermediate_files = "yes"}
	
	if ($delete_intermediate_files eq "yes")
	{
		#####################################
		# Delete files from fastq2bam stage #
		#####################################
		if ($bwa_alignment_method eq "aln")
		{
			if ($paired_ends_type eq "PE")
			{
				&delete_file("$aln_sa1_sai");
				&delete_file("$aln_sa2_sai");
			}
			if ($paired_ends_type eq "SE")
			{
				&delete_file("$aln_sa_sai");
			}
		}

		&delete_file("$aligned_sam");
		


 		###################################
		# Delete files from bam2vcf stage #
		###################################
		&delete_file("$rg_bam");
		&delete_file("$rg_bai");
		&delete_file("$runtitle_dedup_bam");
		&delete_file("$runtitle_dedup_bai");
		&delete_file("$runtitle_dedup_bam_intervals");
		&delete_file("$runtitle_clean_dedup_bam");
		&delete_file("$runtitle_clean_dedup_bai");
		&delete_file("$runtitle_clean_dedup_recal_bam");
		&delete_file("$runtitle_clean_dedup_recal_bai");
		&delete_file("out1.junk");
		&delete_file("out2.junk");


		##################################################
		# If final BAM has been produced, delete f2b BAM #
		##################################################
		if (-e "$results_all_folder/$final_bam")
		{
			$filesize_main = -s "$results_all_folder/$final_bam";

			if ($filesize_main < 100)
			{
				&delete_file("$results_all_folder/$final_f2b_bam");
				&delete_file("$results_all_folder/$final_f2b_bai");
			}
		}

		# Recal files
		&delete_file ("$runtitle_sample_recal_grp");
	}
	 # if delete files is yes


	
	####################################################################
	# Calculate the elapsed total run time at this point (end of loop) #
	####################################################################
	$current_time = time();
	$run_time = $current_time - $start_time;


	##########################################
	# Send an e-mail at the end of each loop #
	##########################################
	open(MAIL, "|/usr/sbin/sendmail -t");

	## Mail Header
	print MAIL "To: $email_address\n";
	print MAIL "From: $from\n";
	print MAIL "Subject: FASTQ2VCF: Sample $sample_name in run $run_title has finished\n\n";
	## Mail Body
	print MAIL "Run title:     \t$run_title\n\n";
	print MAIL "Sample name:   \t$sample_name  (sample number $loop_count)\n\n";
	print MAIL "Script:        \tfastq2vcf PERL script version $version\n\n";

	print MAIL "The processing of this FASTQ file into a final BAM file is complete.\n\n";
	
	print MAIL "Run time so far : $run_time seconds\n";
	printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
	close(MAIL);

	print COMMAND_LOG "\n\n";
	print COMMAND_LOG "========================================================================\n";
	print COMMAND_LOG "                       End of loop $loop_count                          \n";
	print COMMAND_LOG "========================================================================\n\n\n";

	###################################
	# Record time at end of each loop #
	###################################
	$current_time = time();
	$loop_time = $current_time - $start_loop_time;

	printf COMMAND_LOG "  Total time at the end of loop $loop_count:    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $end_time)[7,2,1,0];
	printf COMMAND_LOG "  Time for loop $loop_count:                    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $end_time)[7,2,1,0];


	##############################################
	# Record time for making the clean BAM files #
	##############################################
	$current_time = time();
	$make_clean_bam_time = $current_time - $start_clean_time;
	$make_clean_bam_time_total = $make_clean_bam_time_total + $make_clean_bam_time;
	
} # End of MA IN LOOP

###########################################################################################################
######################################### END OF MAIN LOOP ################################################
###########################################################################################################


##########################################################
#  Make a VCF file using all the BAM files in parallel   #
##########################################################


&print_message("Making VCF with all BAM files in parallel.  File of file name is $list_file","message");

print COMMAND_LOG "\n================================================================================\n";
print COMMAND_LOG   "Making VCF with all BAM files in parallel.  File of file name is $list_file\n";
print COMMAND_LOG "\n================================================================================\n\n";

################################
# Record time before VCF stage #
################################
$start_vcf_time = time();

$running_input_string = $running_input_string." ";


#####################################################################################
# Variant calling to make VCF file using either HaplotypeCaller or UnifiedGenotyper #
#####################################################################################

if ($variant_caller eq "UnifiedGenotyper")
{
	if ($use_default_stand_values eq "false")
	{
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm BOTH $running_input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o $final_VCF_file -S $GATK_validation_stringency -nct $no_of_threads","Run GATK HaplotypeCaller to call variants");	
	}
	if ($use_default_stand_values eq "true")
	{
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm BOTH $running_input_string -o $final_VCF_file -S $GATK_validation_stringency -nct $no_of_threads","Run GATK HaplotypeCaller to call SNP variants");	
	}
	&record_output_file_size ("$final_VCF_file");	
} # UG

if ($variant_caller eq "HaplotypeCaller")
{
	if ($use_default_stand_values eq "false")
	{
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller $GATK_region_string  $running_input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o $final_VCF_file -S $GATK_validation_stringency -nct $no_of_threads","Run GATK HaplotypeCaller to call variants");	
	}
	if ($use_default_stand_values eq "true")
	{
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller $GATK_region_string  $running_input_string -o $final_VCF_file -S $GATK_validation_stringency -nct $no_of_threads","Run GATK HaplotypeCaller to call variants");	
	}
	&record_output_file_size ("$final_VCF_file");	
} # HC


if ($bam_missing_for_variant_calling == 0)
{
	&print_message("A VCF file using all the BAM files in parallel has been created using GATK $variant_caller","message");
}

if ($bam_missing_for_variant_calling > 0)
{
	&print_message("VCF file was not created as not all BAM files were found","warning");

	print "Please check their location and use run_variant_caller later\n\n";
}

################################
# Record time after VCF stage  #
################################
$current_time = time();
$end_time = $current_time - $start_time;
$make_vcf_time = $current_time - $start_vcf_time;

$total_time = $make_raw_bam_time_total + $make_clean_bam_time_total + $make_vcf_time;

print COMMAND_LOG "\n";

&log_time($make_raw_bam_time_total,"Total time for making raw BAM files");
&log_time($make_clean_bam_time_total,"Total time for making clean BAM files");
&log_time($make_vcf_time,"Total time for making VCF files");
&log_time($total_time,"Total time for all stages - added up");
&log_time($end_time,"Total time for all stages");


##############################################################
# Annotate SNPs and Indels with the variant effect predictor #
##############################################################

if ($annotate_variants eq "yes")
{
	print COMMAND_LOG "\n#######################\n";
	print COMMAND_LOG "# Annotating variants #\n";
	print COMMAND_LOG "#######################\n\n";

	if ($variant_effect_predictor eq "vep")
	{
		&run_unix_command("perl /home/genetics/scripts/run_variant_effect_predictor.pl --file $final_VCF_file --ref $species --chr $chromosome","");
	}
	if ($variant_effect_predictor ne "vep")
	{
		print COMMAND_LOG "\n\nTo annotate your variants in the final VCF file $final_VCF_file, use run_variant_effect_predictor or run_snpEff\n\n";
	}
} # if annotate_variants eq "yes"

if ($annotate_variants eq "no")
{
	print COMMAND_LOG "\n\nTo annotate your variants in the final VCF file $final_VCF_file, use run_variant_effect_predictor or run_snpEff\n\n";
} # if annotate_variants eq "no"

#####################################################
# Move final VCF file to the results folder         #
#####################################################
&move_to_results_all_folder ("$final_VCF_file");
&move_to_results_all_folder ("$final_VCF_idx_file");



$current_time = time();
$run_time = $current_time - $start_time;

open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $email_address\n";
print MAIL "From: $from\n";
print MAIL "Subject: FASTQ2VCF: Run $run_title has finished\n\n";
## Mail Body
print MAIL "Run title:  $run_title\n\n";
print MAIL "Script:     fastq2vcf version $version\n\n";

print MAIL "Your next generation sequence analysis is complete\n\n";
print MAIL "For a list of commands see $command_log\n\n";
print MAIL "For more pipeline details see $log_file and README files\n\n";
print MAIL "For further details on intermediate 'raw' BAM files (f2b fast2bam files - now deleted):\n\n";
print MAIL "  $validate_f2b_out\n";
print MAIL "  $flagstat_f2b_out\n\n";

print MAIL "For further details on final BAM files:\n\n";
print MAIL "  $validate_final_out\n";
print MAIL "  $flagstat_final_out\n\n";

printf MAIL"Total run time: %d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
close(MAIL);


#############################
# Turn off logging          #
#############################

close(STDERR);
close (COMMAND_LOG);


###################################
# Write overall final README file #
###################################
$readme_overall_file = "$run_title"."_final_readme.out";

open (READMEOVERALL, ">$readme_overall_file"); 

print	READMEOVERALL	"==================================\n";
print	READMEOVERALL	"Summary of fastq2vcf results files\n";
print	READMEOVERALL	"==================================\n\n\n";


print	READMEOVERALL	"RAW BAM ALIGNMENT FILES FROM BWA\n\n";

print	READMEOVERALL	"\t$final_f2b_bam	     \tRaw alignments to the reference.\n\n";
print	READMEOVERALL	"\t(These are deleted but flagstat and validate files still exist for checking)\n\n\n";

print	READMEOVERALL	"FINAL BAM ALIGNMENT FILES AFTER PROCESSING BY GATK\n\n";

print	READMEOVERALL	"\t$final_bam		     \tBest alignments to the reference after processing by GATK\n\n\n";

print	READMEOVERALL	"\tLook at these files to check the quality of your final BAM files:\n\n";

print   READMEOVERALL   "\t\t$validate_final_out\n\n";
print   READMEOVERALL   "\t\t$flagstat_final_out\n\n\n";

print	READMEOVERALL	"PDF FILES\n\n";

print	READMEOVERALL	"\t$gc_bias_pdf          \tHistogram of GC content of aligned reads\n\n";
print	READMEOVERALL	"\t$insert_size_pdf      \tInsert size histogram\n\n\n";

print	READMEOVERALL	"LOG FILES\n\n";

print	READMEOVERALL	"\t$command_log	         \tCommand log for the NGS pipeline (very useful for tracking errors)\n\n";
print	READMEOVERALL	"\t$log_file             \tRun log for the NGS pipeline\n\n\n";


print	READMEOVERALL	"NOTES\n\n"; 
print	READMEOVERALL	"\t.bam files must have an associated .bai index file for loading into IGV\n\n\n";

print	READMEOVERALL	"NEXT STEPS\n\n"; 
print   READMEOVERALL	"\tWhen you get your final VCF file you need to do two things:\n\n";
print   READMEOVERALL	"\t\t1: Use run_variant_effect_predictor (or run_snpEff) to annotate your variants.\n\n";
print   READMEOVERALL	"\t\t2: Use vcf2excel to prepare your VCF file for NGS SNP Viewer.\n\n";

close (READMEOVERALL);

&move_to_results_all_folder ("$readme_overall_file");

#########################################################################################

######################################################
# Copy command log to command log folder in genetics #
######################################################

&run_unix_command("cp $command_log /home/genetics/command_logs/$command_log","Copy command log to genetics folder");

&move_to_results_all_folder("$log_file");
&move_to_results_all_folder("$command_log");




&print_message("FASTQ2VCF ANALYSIS COMPLETE! YOU HAVE BEEN SENT AN EMAIL WITH DETAILS","message");

print "  Run title: $run_title\n\n";
print "  Run time: ";
printf "%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];

print "\n\n  Results folder:                 \t$results_all_folder";

print "\n\n  For details of the analysis run:\tcheck $command_log and $log_file for full details\n\n";

print "  For checks on the BAM file produced:\tcheck $validate_final_out and $flagstat_final_out\n\n";

&print_message("Your final VCF file is called $final_VCF_file","message");

if ($ref_seq_name eq "canfam3")
{
	print "  - When you get your final VCF file you need to do two things:\n\n";
	print "    1: Use run_variant_effect_predictor (or run_snpEff) to annotate your variants.\n\n";
	print "    2: Use vcf2excel to prepare your VCF file for NGS SNP Viewer.\n\n";
}
exit;

#################################################################################################################
#################################################################################################################


#################################################
# Subroutine to move file to results_all folder #
# (Top folder conatining other sub folders)     #
#################################################

sub move_to_results_all_folder
{

	my $file_to_be_moved = "";	

	$file_to_be_moved = $_[0];

	if (! -e "$file_to_be_moved")
	{
		print "\n$file_to_be_moved could not be found to be moved to $results_all_folder\n";
		print COMMAND_LOG "\n$file_to_be_moved could not be found to be moved to $results_all_folder\n";
	}
	
	if (-e "$file_to_be_moved")
	{
		$command = "mv  $file_to_be_moved $results_all_folder/$file_to_be_moved";
		print "\n$file_to_be_moved was moved to $results_all_folder\n";
		system("$command");
		print COMMAND_LOG "\n$file_to_be_moved was moved to $results_all_folder\n";
	}

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

sub run_unix_command($;$)
{
	my $unix_command = "";
	my $step = "";	
	$unix_command = $_[0];
	$step = $_[1];
	print "\n";
	print("$unix_command\n");
	system("$unix_command");
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	
	if (($loop_count > 0 ) && ($loop_count <= $no_of_files)) {print COMMAND_LOG "File: $loop_count/$no_of_files \tStep: $step\n";}
	else {print COMMAND_LOG "Step: $step\n";}
	
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
				print MAIL "Subject: NGS ANALYSIS fastq2vcf ZERO FILE SIZE: Run $run_title.  Sample $sample_name. File $outputfile.\n\n";
				## Mail Body
				print MAIL "fastq2vcf script version $version\n\n";
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

sub log_time
{
	my $time_to_log = "";
	my $description = "";
	$time_to_log = $_[0];
	$description = $_[1];

	printf " $description: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $time_to_log)[7,2,1,0];
	printf COMMAND_LOG " $description: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $time_to_log)[7,2,1,0];
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