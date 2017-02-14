#!/usr/bin/perl -w

#################################################################################
#									                                            #      
#	fastq2vcf: processes FASTQ files through a pipeline to end up as VCF files  #
#									                                            #
#################################################################################

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
use Getopt::Long;
use File::Basename;
use Term::ANSIColor;
use Cwd;
use Sys::Hostname;
use POSIX qw(strftime);

# VERSION OF SOFTWARE #
#my $version						= "c71";  # 'a' means Workstation 1, "b" is Workstation 2 and "c" is Workstation 3
my $warning						= "on";

my %hostnames = (
	'gen-x1404-ws01'	=> 'ws1',	
	'GEN-X1404-WS02'	=> 'ws2',
	'ws3'				=> 'ws3',
	'workstation3'		=> 'ws3',
	'NGSBACK3'			=> 'samba64',
);

############################################
# For testing new versions of the software #
############################################
my $bwa_path					="/opt/bwa/bwa"; # This is the current version of bwa. Older versions are in folders labeled /opt/bwa062 or whatever...
my $gatk_directory				= "gatk-nightly"; # "gatk_v1" is the previous GATK version and "gatk" is the current GATK version 3.5
my $samtools_path				= "/opt/samtools/samtools";
my $picard_path					= "/opt/picard"; # slash added in command string

###############################################################################
#            GATK memory suggestions from the Broad website                   #
#-----------------------------------------------------------------------------#
#Tool				RTC		IR		BR		PR		RR		   UG             #
#Available modes 	NT 		SG 	  NCT,SG  	NCT	  	SG      NT,NCT,SG         #
#Cluster nodes 		1 		4 		4 		1 		4 	    4 / 4 / 4         #
#CPU threads (-nct) 1 		1 		8 		4-8 	1       3 / 6 / 24        #
#Data threads (-nt) 24 		1 		1 		1 		1       8 / 4 / 1         #
#Memory (Gb) 		48 		4 		4 		4 		4     32 / 16 / 4         #
###############################################################################


############################################################
# Memory and threading settings (can be altered later on)  #
#----------------------------------------------------------#
#                              Workstation      samba64    #
#                                                          #
# Java memory                       60             4       #
# $no_of_threads_bwa                8              2       #
#                                                          #
# $no_of_threads_gatk_ug_nct        1              1       #
# $no_of_threads_gatk_ug_nt         16             2       #
#                                                          #
# $no_of_threads_gatk_pr_nct        8              2       #
# $no_of_threads_gatk_br_nct        8              2       #
# $no_of_threads_gatk_rtc_nt        16             2       #
#                                                          #
############################################################
#if (index($region,":") == -1)


my $memory						= "4"; # memory setting for Java in gigabytes
my $no_of_threads_bwa			= "2";

my $no_of_threads_gatk_ug_nct	= "1";  # UnifiedGenotyper -nct (or HaplotypeCaller) (UG has to be 1)
my $no_of_threads_gatk_ug_nt	= "2";  # UnifiedGenotyper -nt (or HaplotypeCaller)

my $no_of_threads_gatk_pr_nct	= "2";  # PrintReads -nct                <== This one is used for PrintReads
my $no_of_threads_gatk_br_nct	= "2";  # BaseRecalibrator -nct          <== This one is used for BaseRecalibrator
my $no_of_threads_gatk_rtc_nt	= "2";
my $workstation 				= "unknown";
my $e_mail_from 				= ""; # Who e-mails come from
my $host	 					= hostname();
my $hostname 					= $hostnames{$host};

# Workstations 1 2 and 3 are denoted by 'a', 'b' and 'c' here (in the version name)
# Workstation 1/Workstation 2
if ($hostname eq 'ws1' || $hostname eq 'ws2')
{
	$memory						= "45"; # memory setting for Java in gigabytes
	$no_of_threads_bwa			= "8";
	$no_of_threads_gatk_ug_nct	= "8";  # UnifiedGenotyper -nct (or HaplotypeCaller) (UG has to be 1)
	$no_of_threads_gatk_ug_nt	= "16";  # UnifiedGenotyper -nt (but NOT HaplotypeCaller)
	$no_of_threads_gatk_pr_nct	= "8";  # PrintReads -nct                <== This one is used for PrintReads
	$no_of_threads_gatk_br_nct	= "8";  # BaseRecalibrator -nct          <== This one is used for BaseRecalibrator
	$no_of_threads_gatk_rtc_nt	= "16"; # RealignerTargetCreator         <== This one is used for RealignerTargetCreator
	$workstation 				= "true";
	$e_mail_from 				= 'NGS_analysis@'.$hostname.'.aht.org.uk'; # Who e-mails come from
}
#Workstation 3
elsif  ($hostname eq 'ws3')
{
	$memory						= "50"; # memory setting for Java in gigabytes
	$no_of_threads_bwa			= "24";
	$no_of_threads_gatk_ug_nct	= "24";  # UnifiedGenotyper -nct (or HaplotypeCaller) (UG has to be 1)
	$no_of_threads_gatk_ug_nt	= "24";  # UnifiedGenotyper -nt (but NOT HaplotypeCaller)
	$no_of_threads_gatk_pr_nct	= "24";  # PrintReads -nct                <== This one is used for PrintReads
	$no_of_threads_gatk_br_nct	= "24";  # BaseRecalibrator -nct          <== This one is used for BaseRecalibrator
	$no_of_threads_gatk_rtc_nt	= "24"; # RealignerTargetCreator         <== This one is used for RealignerTargetCreator
	$workstation 				= "true";
	$e_mail_from 				= 'NGS_analysis@'.$hostname.'.aht.org.uk'; # Who e-mails come from
}
elsif  ($hostname eq 'samba64')
# Samba64 version (if no 'w' in version name)
{
	$memory						= "4"; # memory setting for Java in gigabytes
	$no_of_threads_bwa			= "2";
	$no_of_threads_gatk_ug_nct	= "1";  # UnifiedGenotyper -nct (or HaplotypeCaller) (UG has to be 1)
	$no_of_threads_gatk_ug_nt	= "2";  # UnifiedGenotyper -nt (but NOT HaplotypeCaller)
	$no_of_threads_gatk_pr_nct	= "2";  # PrintReads -nct                <== This one is used for PrintReads
	$no_of_threads_gatk_br_nct	= "2";  # BaseRecalibrator -nct          <== This one is used for BaseRecalibrator
	$no_of_threads_gatk_rtc_nt	= "2"; # RealignerTargetCreator         <== This one is used for RealignerTargetCreator
	$workstation 				= "false";
	$e_mail_from 				= 'NGS_analysis@'.$hostname.'.aht.org.uk'; # Who e-mails come from
}
else{
	
}

#Other constants
my $testing_mode				= "off";
my $bwa_alignment_method		= "mem"; # replaces aln
my $variant_caller				= "UnifiedGenotyper"; # could be HaplotypeCaller
my $picard_validation_stringency = "SILENT";
my $GATK_validation_stringency	= "LENIENT";
my $filter_mismatching_string	= " --filter_mismatching_base_and_quals"; # Fixes a problem with mismatching qualities
my $validate_bam_files			= "yes";
my $max_open_temp_files			= 7900;
my $remove_duplicates			= "yes";
my $keep_all_for_first_file		= "no"; # This keeps all intermediate files for first file, for error-checking
my $delete_intermediate_files	= "yes"; # set to 'no' for de-bugging
my $delete_intermediates_option = "";

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
my $zero_error_count			= 0;
my $file_count					= 0;
my $check_count					= 0;
my $pos							= 0;
my $no_of_files_in_dir			= 0;
my $filesize_main				= 0;
my $final_VCF_size				= 0;
my $final_gVCF_size				= 0;
my $final_vep_size				= 0;
my $final_snpEff_vcf_size		= 0;

#Boolean
my $use_default_stand_values	= ""; # 'yes' or 'no'
my $fix_misencoded_qual_scores	= ""; # 'yes' or 'no'
my $use_defined_region			= ""; # yes or no
my $annotate_variants			= ""; # yes or no
my $use_other_ref_sequence 		= "no"; # yes or no
my $third_column_found			= "false"; # columns in put file of file names used as new sample names
my $second_column_found			= "false"; # columns in put file of file names
my $whole_genome_sequence		= ""; # true or false
my $choice_ok					= "false";
my $make_VCF_file				= "true";  # true or false. Do you want VCF file as well as gVCF file
my $bgzip_output_gVCFs			= ""; # True or false
my $pcr_free					= "false"; # true or false
my $make_gVCF_file_outside_loop	= "false"; # Under development. Leave as "false" for now

#strings
my $quality_scores				= ""; # Can be "New" or "Old". Assigned by user during program
my $memory_string				= ""; # e.g. -Xmx4g made up from 4 Gigabytes
my $pcr_free_string				= ""; # New setting for PCR free data
my $input						= "";
my $GATK_fix_quals_string		= ""; # GATK string to deal with old style quality scores
my $zero_size_reported			= "no"; # To make sure program only e-mails about zero file size once per loop
my $unix_error_reported			= "no";
my $effect_predictor			= ""; # snpEff or variant_effect_predictor
my $command						= "";
my $ref							= "";
my $ref_prefix					= "";
my $fastq_file_1				= ""; # first FASTQ file
my $fastq_file_2				= ""; # second FASTQ file
my $fastq_file_1_gzipped		= "";
my $fastq_file_2_gzipped		= "";
my $fastq_file_1_prefix			= "";
my $fastq_file_2_prefix			= "";
my $fastq_file_2_guess			= "";
my $paired_ends_type			= "";
my $email_address				= "";
my $run_title					= "";
my $list_file					= "";
my $read_file_method			= "";
my $answer						= "";
my $single_line					= "";
my $results_all_folder			= "";  # Folder for all output files
my $ref_seq_name				= "";  # Name of reference sequence for pindel analysis
my $screen_log_file				= "";  # STDERR screen log of run
my $sample_name					= "";  # Name of each sample. Used to give the .vcf file an individual sample name
my $new_sample_name				= "";  # sample name given in input file (optional)

my $lib							= ""; #This is used by the ReadGroups section and i have made it let us record the reference sequence
my $variant_caller_input_string		= ""; # string of BAM file names for the HaplotypeCaller
my $GVCF_input_string			= ""; # string of gVCF file names for the gVCF
my $species						= "";
my $char_1						= "";
my $char_2						= "";
my $back_to_one					= "";
my $file_to_find				= ""; 
my $current_directory			= "";
my $javabin						= ""; # Place where java is

#Timers
my $start_time					= ""; # The start of the whole pipeline
my $start_raw_time				= ""; # The start of making the raw BAM file
my $start_clean_time			= ""; # The start of "cleaning" the raw (raw) BAM file to the "final" BAM file
my $start_VCF_time				= ""; # The start of making the VCF files from the BAM files (variant calling)
my $start_gVCF_time				= ""; # The start of making the VCF files from the BAM files (variant calling)
my $run_time_so_far				= ""; # Running total of how much time has elapsed
my $total_run_time				= "";
my $current_time				= "";
my $total_time					= ""; # sum of raw bam, clean bam and VCF
my $start_loop_time				= "";
my $loop_time					= "";
my $last_current_time			= "";
my $stage_time					= ""; # Time for each stage of the pipeline
my $make_raw_bam_time			= ""; # Total time to make the raw BAM file
my $make_clean_bam_time			= ""; # Total time to make the clean BAM file
my $make_VCF_time				= ""; # Total time to make the VCF file from the BAM files
my $make_gVCF_time				= ""; # Total time to make the gVCF file from the BAM files
my $make_raw_bam_time_total		= 0;
my $make_clean_bam_time_total	= 0;
my $make_gVCF_time_total		= 0;
my $make_VCF_time_total			= 0;

# Region strings
my $region						= "";
my $GATK_region_string			= ""; # String to be used with the -L option of GATK tools such as UnifiedGenotyper
my $GATK_chr_only_region_string = ""; # String to be used with Unified Genotyper (just narrows down to chromosome e.g. -L chr15)
my $GATK_known_sites_string		= ""; # String for known SNPs for BaseRecalibrator/CountCovariates
my $chromosome					= "";

#### File names ####
my $aln_sa_sai						= "";
my $aln_sa1_sai						= "";
my $aln_sa2_sai						= "";
my $aligned_sam						= "";
my $aligned_sorted_rg_bam			= "";
my $aligned_sorted_rg_bai			= "";
my $command_log						= "";
my $command_times_log				= "";
my $output_gVCF_file				= "";
my $output_gVCF_file_idx			= "";

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

my $final_raw_bam					= ""; # Final BAM file after FASTQ2BAM pipeline has run
my $final_raw_bai					= "";

#VCF files for SNPs and Indels
my $final_VCF_file				= "";
my $final_VCF_idx_file			= "";
my $final_gVCF_file				= ""; # final multi-column gVCF file made outside the loop

#LOG files
my $validate_raw_out			= "";
my $flagstat_raw_out			= "";
my $validate_final_out			= "";
my $flagstat_final_out			= "";

#Other files
my $insert_size_pdf				= ""; # Used by picard CollectInsertSizeMetrics
my $insert_size_out				= ""; # Used by picard CollectInsertSizeMetrics
my $WGS_metrics_file			= ""; # Used by picard CollectWgsMetrics
my $dummy_dbsnp_file			= ""; # The dummy dbSNP file used by....
my $all_snps_file				= ""; # SNPs file for whole genome
my $actual_dbsnp_file			= ""; # Correct dbSNP file for this chromosome (if it exists)
my $dbsnp_file					= ""; # dbsnp file used by GATK
my $readme_file					= "";
my $readme_overall_file			= "";
my $vep_out_file				= "";
my $vep_html_file				= "";
my $vep_command_log				= "variant_effect_predictor_command_log.out";
my $final_snpEff_VCF_file		= "";

my $runtitle_sample_recal_grp 	= ""; # For BaseRealigner

my @fastq1_file_array			= ();
my @fastq2_file_array			= ();
my @sample_name_array			= ();
my @file_name_array				= (); # to store list of FASTQ files in the DIR when entering single file
my @item						= ();

# Set java...(depends on whether this is Workstation or samba64)
$javabin="java";

##############################################################
# Get command line options --delete_intermediate_files "yes" #
##############################################################
GetOptions("delete_intermediates:s"=>\$delete_intermediates_option);

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

print "Running for ".$hostname."\n\n";

print "  - This program processes FASTQ files into VCF files.\n\n";

print "    All the output files are put into a single folder.\n\n";

print "    The intermediate BAM files with 'raw alignments' are deleted and the\n";
print "    'final' BAM files after processing through the GATK stages are kept.\n\n\n";


print color 'bold cyan';

print "  - When you get your final VCF file you might want to do two things:\n\n";
print "    1: Use run_vep or run_snpEff to annotate your variants.\n\n";
print "    2: Use vcf2excel to prepare your VCF file for NGS SNP Viewer.\n\n";

print color 'reset';

if ($warning eq "on")
{
	&print_message("Using fastq2vcf on the Workstation","warning");
	print "Please check with other users of this Workstation before running on large files.\n\n";
}


if ($delete_intermediates_option eq "no")
{$delete_intermediate_files = "no";print "N.B. Delete intermediate files is set to 'no'\n\n";}



#########################################
# Warn if any debugging settings are on #
#########################################
if ($keep_all_for_first_file eq "yes")
{
	&print_message("Debugging is set to keep all intermediate files for the first BAM file.","warning");
	print "This means that very large intermediate files will not be deleted\n\n";
	&pause;
}
if ($delete_intermediate_files eq "no")
{
	&print_message("Debugging is set to keep all intermediate files.","warning");
	&pause;
}






#################################################################
# Name of run - this is attached to all the output files        #
#################################################################
$run_title = "ls";
until ($run_title ne "ls")
{
	&print_message("Type a name for this run (fairly short - no spaces)","input");

	$run_title = <STDIN>; chomp $run_title;

	if ($run_title eq "ls")
	{
		print "\n";
		system ("ls *_fastq2vcf_command_log.out");
		print "\n";
	}

	#####################################
	# Warn user if run title has spaces #
	#####################################
	if (index ($run_title," ") > -1)
	{
		$run_title =~ s/\ /_/g;

		&print_message("You can't have spaces in the run title $run_title","warning");

		print "Run title converted to:   $run_title\n\n";
		print "  >>> Press 'return' if this is OK   ";
		$answer=<STDIN>;
	}
}


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

$screen_log_file = "$run_title"."_fastq2vcf_screen_log.out";

$| = 1;

open(STDERR, "| tee $screen_log_file");


#################################################################
# E-mail address for notifications                              #
#################################################################

&print_message("Please enter your email address","input");
$email_address = <STDIN>;
chomp $email_address;

# Some short cuts
if ($email_address eq "e"){$email_address = 'ellen.schofield@aht.org.uk';}
if ($email_address eq "m"){$email_address = 'mike.boursnell@aht.org.uk';}
if ($email_address eq "s"){$email_address = 'sally.ricketts@aht.org.uk';}
if ($email_address eq "b"){$email_address = 'rebekkah.hitti@aht.org.uk';}
if ($email_address eq "g"){$email_address = 'graham.newland@aht.org.uk';}
if ($email_address eq "l"){$email_address = 'louise.pettitt@aht.org.uk';}

#################################################################
# Ask if the data is SE or PE data? (single or paired-end)      #
#################################################################
&print_message("Do you have a Single-end or Paired-end dataset?","input");
print "  <1>  Single-end\n";
print "  <2>  Paired-end [DEFAULT]\n\n";

$answer = <STDIN>;
chomp $answer;
if ($answer eq ""){$answer = "2"}

if (substr($answer,0,1) eq "1" ){$paired_ends_type = "SE"}
if (substr($answer,0,1) eq "2" ){$paired_ends_type = "PE"}



#####################################################################
# ASK IF YOU WANT TO READ THE FILENAMES FROM A "FILE OF FILE NAMES" #
#####################################################################

$answer = "";

until ($read_file_method eq "multiple" || $read_file_method eq "single")
{
	&print_message("How do you want to read in the input FASTQ files?","input");

	print "   <1> Using a file of file names for your FASTQ files [DEFAULT]\n";

	if ($paired_ends_type eq "SE"){print "   <2> Enter a single FASTQ file\n\n";}
	if ($paired_ends_type eq "PE"){print "   <2> Enter a pair of FASTQ files individually\n\n";}

	$answer = <STDIN>;chomp $answer;
	if (substr($answer,0,1) eq "2"){$read_file_method = "single"} else {$read_file_method = "multiple"}
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

		############################################################
		# Assign to array of file names but just use first element #
		############################################################
		$fastq1_file_array[1]=$fastq_file_1;
		$no_of_files = 1;
	}
	
	if ($paired_ends_type eq "PE")
	{
		############################################################
		# Get list of FASTQ files in this directory                #
		############################################################
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
		
		############################################################
		# Guess second FASTQ file name                             #
		############################################################
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
		
		############################################################
		# Assign to array of file names but just use first element #
		############################################################
		$fastq1_file_array[1]=$fastq_file_1;
		$fastq2_file_array[1]=$fastq_file_2;
		$sample_name_array[1] = "";
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
			print "(The two paired-end FASTQ file names must be on the same line separated by a TAB):      ";
		}
		
		$list_file = <STDIN>;chomp $list_file;
		
		if (($list_file eq "ls") || ($list_file eq "")){print "\n";system ("ls *.txt"); print "\n";}
	}
	
	############################################################
	# Make sure the list file is in Unix format                #
	############################################################
	$command = "dos2unix $list_file";
	system("$command");
	

	############################################################
	# Open the list file to get the list of file names         #
	############################################################
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
				@item=split(/\s/,$single_line);
				$array_size = scalar @item;
				
				if ($array_size == 1)
				{
					$fastq_file_1 = $single_line;
					$fastq1_file_array[$list_count]=$fastq_file_1;
				}
				
				if ($array_size == 2)
				{
					$fastq_file_1 = $item[0];
					$new_sample_name = $item[1];
					
					$fastq1_file_array[$list_count]=$fastq_file_1;
					$sample_name_array[$list_count] = $new_sample_name;
					$second_column_found = "true";
				} # array size 2, SE
			} # SE
			
			if ($paired_ends_type eq "PE") # Read two columns separated by a TAB for the two files FASTQ_1 and FASTQ_2
			{
				@item=split(/\s/,$single_line); # Can be TAB or white space (in case user makes a mistake)
				
				$array_size = scalar @item;
				
				if ($array_size == 1)
				{
					&print_message("FILE FORMAT ERROR!","warning");
					print "This file only appears to have a single column\n\n";
					print "For paired-end (PE) files you need two columns separated by a tab.\n\n";
					close LIST;
					exit;
				}# array_size = 1
				

				##########################################################
				# If there is NO third column for the sample name        #
				##########################################################
				if ($array_size == 2)
				{
					$fastq_file_1 = $item[0];
					$fastq_file_2 = $item[1];
					

					#############################################
					# Warn user if both file names are the same #
					#############################################
					if ($fastq_file_1 eq $fastq_file_2)
					{
						&print_message("Both FASTQ files are the same","warning");

						print "Both FASTQ files are called $fastq_file_1\n\n";
						print "Please fix this and start again\n\n";
						exit;
					}

					##########################################################
					# Create a sample name based on the two FASTQ file names #
					##########################################################
					for ($pos = 0; $pos <= length ($fastq_file_1); $pos++)
					{
						$char_1 = substr ($fastq_file_1,$pos,1);
						if ($char_1 eq "1")
						{
							$char_2 = substr ($fastq_file_2,$pos,1);
							if ($char_2 eq "2")
							{
								$back_to_one = substr($fastq_file_2,0,$pos)."1".substr($fastq_file_2,$pos+1,99);
								
								if ($fastq_file_1 eq $back_to_one)
								{
									$sample_name = substr($fastq_file_2,0,$pos).substr($fastq_file_2,$pos+1,99);
									$sample_name = &get_prefix($sample_name);

									if (rindex($sample_name,"_") == length($sample_name) - 1)
									{
										$sample_name = substr($sample_name,0,length($sample_name) - 1);
									}
								}
							} # if
						}
					} # guess 2nd file name

					$sample_name_array[$list_count] = $sample_name;
						

					###################################################
					# If it is not gzipped, make name of gzipped file #
					###################################################
					if (index($fastq_file_1,".gz") == -1)
					{
						$fastq_file_1_gzipped = "$fastq_file_1".".gz";
					}if (index($fastq_file_2,".gz") == -1)
					{
						$fastq_file_2_gzipped = "$fastq_file_2".".gz";
					}

					if ((!-e $fastq_file_1) && (-e $fastq_file_1_gzipped)){ $fastq_file_1 = $fastq_file_1_gzipped }
					if ((!-e $fastq_file_2) && (-e $fastq_file_2_gzipped)){ $fastq_file_2 = $fastq_file_2_gzipped }

					$fastq1_file_array[$list_count]=$fastq_file_1;
					$fastq2_file_array[$list_count]=$fastq_file_2;
					$third_column_found = "false";
				}# array_size = 2
				

				##########################################################
				# If there is a third column use this as the sample name #
				##########################################################
				if ($array_size == 3)
				{
					$fastq_file_1 = $item[0];
					$fastq_file_2 = $item[1];
					$new_sample_name = $item[2];
					
					#############################################
					# Warn user if both file names are the same #
					#############################################
					if ($fastq_file_1 eq $fastq_file_2)
					{
						&print_message("Both FASTQ files are the same","warning");

						print "Both FASTQ files are called $fastq_file_1\n\n";
						print "Please fix this and start again\n\n";
						exit;
					}

					###################################################
					# If it is not gzipped, make name of gzipped file #
					###################################################
					if (index($fastq_file_1,".gz") == -1)
					{
						$fastq_file_1_gzipped = "$fastq_file_1".".gz";
					}if (index($fastq_file_2,".gz") == -1)
					{
						$fastq_file_2_gzipped = "$fastq_file_2".".gz";
					}

					if ((!-e $fastq_file_1) && (-e $fastq_file_1_gzipped)){ $fastq_file_1 = $fastq_file_1_gzipped }
					if ((!-e $fastq_file_2) && (-e $fastq_file_2_gzipped)){ $fastq_file_2 = $fastq_file_2_gzipped }


					$fastq1_file_array[$list_count] = $fastq_file_1;
					$fastq2_file_array[$list_count] = $fastq_file_2;
					$sample_name_array[$list_count] = $new_sample_name;
					$third_column_found = "true";
				} # array_size = 3
			} # PE

			$list_count=$list_count + 1;

		} # if length $single_line > 1

	} # while ($single_line = <LIST> ) 

	close LIST;

	$no_of_files=$list_count - 1;
	

	##########################################################
	# Check if FASTQ files exist                             #
	##########################################################
	&print_message("Checking input files","message");
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		if ($paired_ends_type eq "SE")
		{
			if (! -e "$fastq1_file_array[$list_count]")
			{ 
				print "\n  >>>>>>>>>> File $fastq1_file_array[$list_count] does not exist.\n\n";
				exit;
			} 
		}
		if ($paired_ends_type eq "PE")
		{
			if (! -e "$fastq1_file_array[$list_count]")
			{ 
				print "\n >>>>>>>>>> File $fastq1_file_array[$list_count] does not exist..\n\n";
				exit;
			} 
			
			if (! -e "$fastq2_file_array[$list_count]")
			{ 
				print "\n >>>>>>>>>> File $fastq2_file_array[$list_count] does not exist...\n\n";
				exit;
			} 
		}
	} # list count loop for checking if FASTQ files exist
	

	##########################################################
	# List FASTQ file names for user                         #
	##########################################################
	if ($paired_ends_type eq "SE"){print "\nThere are $no_of_files FASTQ files in this file of file names.\n\n";}
	if ($paired_ends_type eq "PE"){print "\nThere are $no_of_files pairs of FASTQ files in this file of file names.\n\n";}
	
	if ($third_column_found eq "true")
	{
		print color "bold white";
		print "There is also a third column in the input file which will be used as sample prefixes for naming the output files.\n\n";
		print color "reset";
	}
	
	if ($second_column_found eq "true")
	{
		print color "bold white";
		print "There is also a second column in the input file which will be used as sample prefixes for naming the output files.\n\n";
		print color "reset";
	}
	
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		if ($paired_ends_type eq "SE")
		{
			if ($second_column_found eq "true")
			{
				print "File $list_count	\t$fastq1_file_array[$list_count]\t\t$sample_name_array[$list_count]\n";
			}
			if ($second_column_found eq "false")
			{
				print "File $list_count	\t$fastq1_file_array[$list_count]\n";
			}
		} # SE
		if ($paired_ends_type eq "PE")
		{
			if ($third_column_found eq "true")
			{
				print "Pair $list_count	\t$fastq1_file_array[$list_count]  \t$fastq2_file_array[$list_count]\t\t$sample_name_array[$list_count]\n";
			}
			if ($third_column_found eq "false")
			{
				print "Pair $list_count	\t$fastq1_file_array[$list_count]  \t$fastq2_file_array[$list_count]\t\t$sample_name_array[$list_count]\n";
			}
		} # PE
	}

	if ($third_column_found eq "false")
	{
		print color "bold white";
		print "\nIn the third column the script has created sample prefixes which will be used for naming the output files.\n\n";
		print color "reset";
	}

	&print_message("If these file names and sample prefixes look OK, press enter to proceed (or 'Q' to quit)","message");

	$answer = <STDIN>;chomp $answer; 
	if (lc $answer eq "q"){exit;} 

} # End of if ($read_file_method eq "multiple")


#################################################################
# Choose which reference sequene you want                       #
#################################################################

&print_message("Which reference sequence do you want to use?","input");

print "   <1> CanFam3\n";
print "   <2> CanFam2\n\n";

print "   <4> EquCab2\n";
print "   <5> Human\n\n";

print "   <6> Strep. equi\n";
print "   <7> Strep. zoo\n\n";

print "   <9> other\n\n";

$answer = <STDIN>; chomp $answer;

#if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3/canfam3_dummy_DBSNP.vcf"; $all_snps_file = "/home/genetics/canfam3/canfam3_snps_all.vcf" }
if (substr($answer,0,1) eq "1" ){$ref = "/geneticsdata/alces-flight-image/canfam3/ensembl/canfam3.fasta"; $ref_seq_name = "canfam3"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3/canfam3_dummy_DBSNP.vcf"; $all_snps_file = "/geneticsdata/alces-flight-image/canfam3/canfam3_snps_all.vcf" }
if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/canfam2/canfam2.fasta"; $ref_seq_name = "canfam2"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam2/canfam2_dummy_DBSNP.vcf"; $all_snps_file = ""}

if (substr($answer,0,1) eq "4" ){$ref = "/home/genetics/equcab2/equcab2.fasta"; $ref_seq_name = "equcab2"; $species = "equus_caballus";$dummy_dbsnp_file = "/home/genetics/equcab2/equcab2_dummy_DBSNP.vcf"; $all_snps_file = "/home/genetics/equcab2/equcab2_snps_all.vcf"}
if (substr($answer,0,1) eq "5" ){$ref = "/home/genetics/human/human.fasta"; $ref_seq_name = "human"; $species = "homo_sapiens";$dummy_dbsnp_file = "/home/genetics/human/human_dummy_DBSNP.vcf"; $all_snps_file = ""}

if (substr($answer,0,1) eq "6" ){$ref = "/home/genetics/strep_equi/strep_equi.fasta"; $ref_seq_name = "s_equi"; $species = "streptococcus_equi";$dummy_dbsnp_file = "/home/genetics/strep_equi/strep_equi_dummy_DBSNP.vcf"; $all_snps_file = ""}
if (substr($answer,0,1) eq "7" ){$ref = "/home/genetics/strep_zoo/strep_zoo.fasta"; $ref_seq_name = "s_zoo"; $species = "streptococcus_zoo";$dummy_dbsnp_file = "/home/genetics/strep_zoo/strep_zoo_dummy_DBSNP.vcf"; $all_snps_file = ""}


#################################################################
# Choose your own reference sequence (if not on the list above) #
#################################################################
if (substr($answer,0,1) eq "9" )
{
	$use_other_ref_sequence = "yes";
	 
	until (-e "$ref")
	{
		print "What is the full path of your reference sequence?\n";
		print "e.g. /geneticsdata/alces-flight-image/canfam3/ensembl/canfam3.fasta\n\n";
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


#################################################################
# Check if REF file exists (if a new REF file is specified)     #
#################################################################

if (! -e "$ref")
{ 
	&print_message("File $ref does not exist","warning");
	print "  (You might need to make a new_bwa version of the file?)\n\n";
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


#################################################################
# Ask which bwa method you want to use: aln or mem              #
#################################################################

&print_message("Which bwa alignment method do you want to use?","input");

print "   <1> mem - use this if the reads are over 70bp                [DEFAULT]\n";
print "   <2> aln - this is the method for old Illumina reads of less than 70bp\n\n";

$answer = <STDIN>;
chomp $answer;
if ($answer eq ""){$answer = "1"}; #default is mem

if (substr($answer,0,1) eq "1" ){$bwa_alignment_method = "mem";}
if (substr($answer,0,1) eq "2" ){$bwa_alignment_method = "aln";}



#################################################################
# Ask if the fastq file has old Illumina quality scores         #
#################################################################

&print_message("Are your Quality Scores of the correct Sanger type (i.e. new style)?","input");

print "  Illumina scores prior to Illumina 1.8 were 'old-style'              \n";
print "  All current FASTQ files should be OK.                               \n\n";

print "   <1> New style     [DEFAULT]\n";
print "   <2> Old style\n\n";

$answer = <STDIN>;
chomp $answer;

if ($answer eq ""){$answer = "1"} # default is New style

if (substr($answer,0,1) eq "1" ){$fix_misencoded_qual_scores = "no"; $GATK_fix_quals_string = ""; $quality_scores = "New"}
if (substr($answer,0,1) eq "2" ){$fix_misencoded_qual_scores = "yes"; $GATK_fix_quals_string = "-fixMisencodedQuals "; $quality_scores = "Old"}


#################################################################
# NOTE fix quality scores when you get to the BAM files         #
#################################################################


#################################################################
# Would you like to remove duplicates?                          #
# normally YES but might not want to                            #
#################################################################
&print_message("Would you like to remove duplicate reads?","input");
print "(This is normally advised, but for special purposes you might not want to)\n\n";

print "   <1> YES - remove duplicates                [DEFAULT]\n";
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
print "(For targeted sequencing this is strongly advised if possible as it speeds up the whole process)\n\n";

print color "reset";

if (($ref_seq_name eq "s_equi") || ($ref_seq_name eq "s_zoo"))
{
	print "  NB: For strep_equi and strep_zoo genomes use 'chr1'\n\n";
}

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
	
	#################################################################
	# If you are using canfam3 or equcab2  reference                #
	#################################################################
	if (($use_other_ref_sequence eq "no")  && (($ref_seq_name eq "canfam3") || ($ref_seq_name eq "equcab2")))
	{
		###########################################################################
		# If $use_defined_region eq "no" then ask if it is whole genome           #
		###########################################################################
		&print_message("Is this whole genome sequence?","input");

		print "   <1> YES - whole genome sequence\n";
		print "   <2> NO  - single chromosome\n\n";

		$answer = <STDIN>; chomp $answer;

		if (substr($answer,0,1) eq "1"){$whole_genome_sequence = "yes"}
		if (substr($answer,0,1) eq "2"){$whole_genome_sequence = "no"}

		if ($whole_genome_sequence eq "yes")
		{
			if ($all_snps_file ne "")
			{
				$dbsnp_file = $all_snps_file;
				print "\n>> File with all SNPs for $ref_seq_name is $all_snps_file <<\n\n";
				&pause;
			}
			else
			{
				&print_message("There is no file with all the SNPs for $ref_seq_name","warning");
				print "You will have to create one\n\n";
				exit;
			}
			$GATK_known_sites_string = " -knownSites $dbsnp_file";
			
		} # if ($whole_genome_sequence eq "yes")


		################################################################################################
		# If you say it is just one chromosome then user should have chosen use_defined_region earlier #
		################################################################################################
		if ($whole_genome_sequence eq "no")
		{
			&print_message("Earlier on you said you didn't want to focus on a region","message");

			print "I suggest you start again and then say you DO want to focus on a region.\n\n";
			exit;
		} #if ($whole_genome_sequence eq "no")

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

	} # if no colon (e.g. chr15)
	
	
	####################################################
	# If full region given chr15:34000000-390000000    #
	####################################################
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

	#####################################################
	# If no colon or chr then get the chromosome number #
	#####################################################
	if ((index($region,":") == -1) && (index($region,"chr") == 0))
	{
		$chromosome = substr($region,3,99);
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
		
		# SNP file exists
		if (-e $actual_dbsnp_file)
		{
			&print_message("A SNP file for this chromosome has been found","message");

			print "SNP file $actual_dbsnp_file exists\n";
			
			$dbsnp_file = $actual_dbsnp_file;
			$GATK_known_sites_string = " -knownSites $dbsnp_file";
		}
		
		# SNP file not found
		if (! -e $actual_dbsnp_file)
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

print "   <1> GATK Haplotype Caller using gVCF            [DEFAULT]\n";
print "   <2> GATK Haplotype Caller              \n";
print "   <3> GATK Unified Genotyper  (no longer supported by GATK)\n";


$answer = <STDIN>; chomp $answer;

if (substr($answer,0,1) eq "1" ){$variant_caller = "HaplotypeCallerGVCF"}
if (substr($answer,0,1) eq "2" ){$variant_caller = "HaplotypeCaller"}
if (substr($answer,0,1) eq "3" ){$variant_caller = "UnifiedGenotyper"}

if ($answer eq ""){$variant_caller = "HaplotypeCallerGVCF"}


###########################################################################
# Do you want to specify PCR-free?                                        #
###########################################################################
$pcr_free = "false";

if (($variant_caller eq "HaplotypeCaller") || ($variant_caller eq "HaplotypeCallerGVCF"))
{
	&print_message("Is your sequencing data PCR-free?","input");

	print "There is a new option to specify PCR-free data.\n\n";

	print "GATK says 'When working with PCR-free data, be sure to set `-pcr_indel_model NONE`'\n\n";

	print "   <1> -pcr_indel_model = NONE            (for PCR-free data)       [DEFAULT]\n";
	print "   <2> -pcr_indel_model = CONSERVATIVE    (previous default)            \n";


	$answer = <STDIN>; chomp $answer;

	if (substr($answer,0,1) eq "1" ){$pcr_free = "true"}
	if (substr($answer,0,1) eq "2" ){$pcr_free = "false"}

	if ($answer eq ""){$pcr_free = "true"}

	if ($pcr_free eq "true"){$pcr_free_string = " --pcr_indel_model NONE"} else{$pcr_free_string = ""}
}



###########################################################################
# Do you want to compress the gVCF files                                  #
###########################################################################

&print_message("Would you like to compress your output gVCF files?","input");

print "This saves space on the server and still works for GenotypeGVCFs\n\n";

print "   <1> Yes - compress the gVCF files using bgzip and tabix [DEFAULT]\n";
print "   <2> No - don't compress the gVCF files.              \n\n";

$answer = <STDIN>; chomp $answer;

if (substr($answer,0,1) eq "1" ){$bgzip_output_gVCFs = "true"}
if (substr($answer,0,1) eq "2" ){$bgzip_output_gVCFs = "false"}

if ($answer eq ""){$bgzip_output_gVCFs = "true"}



###########################################################################
# If you opt for gVCF do you want VCF as well?                            #
# And do you want to use GenotypeGVCFs to combine any gVCFs made?         #
###########################################################################
$make_VCF_file = "true"; 


if ($variant_caller eq "HaplotypeCaller")
{

	##########################################################################
	# If you are making several BAMs, do you want them combined into a VCF? #
	##########################################################################
	if ($no_of_files > 1)
	{
		&print_message("As you are making $no_of_files BAM files, do you want to combine these into a multi-column VCF file","input");

		print "   <1> YES - make a multi-column VCF file             [DEFAULT]\n";
		print "   <2> NO  - just a BAM file is fine             \n";


		$answer = <STDIN>; chomp $answer;

		if (substr($answer,0,1) eq "1" ){$make_VCF_file = "true"}
		if (substr($answer,0,1) eq "2" ){$make_VCF_file = "false"}
		if ($answer eq ""){$make_VCF_file = "true"}
	}


	##########################################################################
	# If you are making only 1 gVCF file, do you want a VCF as well?         #
	##########################################################################
	if ($no_of_files == 1)
	{
		&print_message("As you have only one sample, do you want a VCF made from this?","input");

		print "Normally we make multi-column VCF files which need multiple samples\n\n";

		print "   <1> YES - make a single column VCF file  \n";
		print "   <2> NO  - just a BAM file is fine              [DEFAULT]\n";

		$answer = <STDIN>; chomp $answer;

		if (substr($answer,0,1) eq "1" ){$make_VCF_file = "false"}
		if (substr($answer,0,1) eq "2" ){$make_VCF_file = "true"}
		if ($answer eq ""){$make_VCF_file = "false"}
	}
} # variant caller is HaplotypeCaller


if ($variant_caller eq "HaplotypeCallerGVCF")
{

	##########################################################################
	# If you are making several gVCFs, do you want them combined into a VCF? #
	##########################################################################
	if ($no_of_files > 1)
	{
		&print_message("As you are making $no_of_files gVCF files, do you want to combine these into a multi-column VCF file","input");

		print "Normally you would combine them later, possibly with several more gVCF files (using gvcf2vcf)\n\n";

		print "   <1> NO  - leave the gVCFs uncombined             [DEFAULT]\n";
		print "   <2> YES - combine these gVCFs into a VCF              \n";


		$answer = <STDIN>; chomp $answer;

		if (substr($answer,0,1) eq "1" ){$make_VCF_file = "false"}
		if (substr($answer,0,1) eq "2" ){$make_VCF_file = "true"}
		if ($answer eq ""){$make_VCF_file = "false"}
	}


	##########################################################################
	# If you are making only 1 gVCF file, do you want a VCF as well?         #
	##########################################################################
	if ($no_of_files == 1)
	{
		&print_message("As you are making only one gVCF file, do you want a VCF version of this as well?","input");

		print "Normally you would just use the gVCF to make a multi-column VCF file later, possibly with several more gVCF files (using gvcf2vcf)\n\n";

		print "   <1> NO  - gVCF only is fine                [DEFAULT]\n";
		print "   <2> YES - make a single VCF file as well              \n";


		$answer = <STDIN>; chomp $answer;

		if (substr($answer,0,1) eq "1" ){$make_VCF_file = "false"}
		if (substr($answer,0,1) eq "2" ){$make_VCF_file = "true"}
		if ($answer eq ""){$make_VCF_file = "false"}
	}
}



##################################################################################
# Filtering preferences for GATK variant callers.                                #
# If you want to chnge the default values of stand_call_conf and stand_emit_conf #
##################################################################################
&print_message("Would you like to use the default values of stand_emit_conf and stand_call_conf?","input");

print "(These are GATK variant calling parameters that control quality thresholds)\n\n";

print "   <1>  Yes - use default values of $stand_emit_conf and $stand_call_conf [DEFAULT]\n";
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

	print "New value for stand_call_conf (default is 30 - a lower value produces more variants):    ";
	$stand_call_conf = <STDIN>;
	chomp $stand_call_conf;
	
	print "New value for stand_emit_conf (default is 30 - a lower value produces more variants):    ";
	$stand_emit_conf = <STDIN>;
	chomp $stand_emit_conf;

	if ($stand_call_conf eq ""){ $stand_call_conf = 30; }
	if ($stand_emit_conf eq ""){ $stand_emit_conf = 30; }
} # New values





##########################################################
# Get user input on all memory and threads options       #
##########################################################
$memory_string = "-Xmx".$memory."g";

 while ($choice_ok eq "false")
{
		&print_message("Current memory and CPU thread settings","message");
		print "You may as well use these settings unless you have good reason to change them)\n\n";

		print "Workstation: $workstation\n\n";

		print "Memory setting for java steps (Gigabytes):                     \t$memory\n";
		print "Number of CPU threads for bwa:                                 \t$no_of_threads_bwa\n";
		print "Number of CPU threads (-nct) for GATK variant calling:         \t$no_of_threads_gatk_ug_nct\n";

		if ($variant_caller eq "UnifiedGenotyper")
		{
			print "Number of data threads (-nt) for GATK variant calling:         \t$no_of_threads_gatk_ug_nt\n";
		}

		print "Number of data threads (-nt) for GATK RealignerTargetCreator:  \t$no_of_threads_gatk_rtc_nt\n";
		print "Number of CPU threads (-nct) for GATK BaseRecalibrator:        \t$no_of_threads_gatk_br_nct\n";
		print "Number of CPU threads (-nct) for GATK PrintReads:              \t$no_of_threads_gatk_pr_nct\n";

		print "\nWould you like to change these? (y/n)    ";
		$answer=<STDIN>;chomp $answer;$answer = lc $answer;

		if ($answer eq "y")
		{
			&print_message("Enter new values (press 'return' to keep existing value)","input");

			$choice_ok = "false";

			print "  Memory setting for java steps (in Gigabytes)                  [current value = $memory]:      ";
			$input = <STDIN>;chomp $input;
			if ($input ne "" && $input <= $memory){$memory = $input}
			#if ($memory > 60){$memory = "60"}

			print "  Number of CPU threads for bwa                                 [current value = $no_of_threads_bwa]:      ";
			$input = <STDIN>;chomp $input;
			if ($input ne "" && $input <= $no_of_threads_bwa){$no_of_threads_bwa = $input}
			#if ($no_of_threads_bwa > 16){$no_of_threads_bwa = "16"}

			print "  Number of CPU threads (-nct) for GATK variant calling?        [current value = $no_of_threads_gatk_ug_nct]:      ";
			$input = <STDIN>;chomp $input;
			if ($input ne "" && $input <= $no_of_threads_gatk_ug_nct){$no_of_threads_gatk_ug_nct = $input}
			#if ($no_of_threads_gatk_ug_nct > 16){$no_of_threads_gatk_ug_nct = "16"}

			if ($variant_caller eq "UnifiedGenotyper")
			{
				print "  Number of data threads (-nt) for GATK variant calling?        [current value = $no_of_threads_gatk_ug_nt]:      ";
				$input = <STDIN>;chomp $input;
				if ($input ne "" && $input <= $no_of_threads_gatk_ug_nt){$no_of_threads_gatk_ug_nt = $input}
				#if ($no_of_threads_gatk_ug_nt > 16){$no_of_threads_gatk_ug_nt = "16"}
			}

			print "  Number of data threads (-nt) for GATK RealignerTargetCreator? [current value = $no_of_threads_gatk_rtc_nt]:      ";
			$input = <STDIN>;chomp $input;
			if ($input ne "" && $input <= $no_of_threads_gatk_rtc_nt){$no_of_threads_gatk_rtc_nt = $input}
			#if ($no_of_threads_gatk_rtc_nt > 16){$no_of_threads_gatk_rtc_nt = "16"}

			print "  Number of CPU threads (-nt) for GATK BaseRecalibrator?        [current value = $no_of_threads_gatk_br_nct]:      ";
			$input = <STDIN>;chomp $input;
			if ($input ne "" && $input <= $no_of_threads_gatk_br_nct){$no_of_threads_gatk_br_nct = $input}
			#if ($no_of_threads_gatk_br_nct > 16){$no_of_threads_gatk_br_nct = "16"}

			print "  Number of CPU threads (-nt) for GATK PrintReads?              [current value = $no_of_threads_gatk_pr_nct]:      ";
			$input = <STDIN>;chomp $input;
			if ($input ne "" && $input <= $no_of_threads_gatk_pr_nct){$no_of_threads_gatk_pr_nct = $input}
			#if ($no_of_threads_gatk_pr_nct > 16){$no_of_threads_gatk_pr_nct = "16"}

		}# Yes to change memory options
		else
		{
			$choice_ok = "true";
		}

} # while choice_ok eq "false"

# Make up java memory_string
$memory_string = "-Xmx".$memory."g";


###################################################################
# Annotate variants preferences                                   #
# (Leave this out. Do this separately with run_snpEff or run_VEP) #
###################################################################

$annotate_variants = "no";
$effect_predictor = "none";


##########################################################
# open Command Log file to record all commands           #
##########################################################
$command_log = "$run_title"."_fastq2vcf_command_log.out";
$command_times_log = "$run_title"."_fastq2vcf_times_log.out";
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
open (COMMAND_TIMES_LOG, ">$command_times_log") || die "Cannot create output file: $command_times_log";

my $version = strftime "%F", localtime($start_time);
print COMMAND_LOG "Command log of perl script fastq2vcf (started $version)\n\n";

&print_message("SUMMARY DETAILS - PLEASE CHECK CAREFULLY!!","message");
print COMMAND_LOG "SUMMARY DETAILS\n\n";

&print_both("  Name of this analysis:         \t$run_title\n"); 
&print_both("  Your email address:            \t$email_address\n");

&print_message("SETTINGS","message");

print COMMAND_LOG "\nSETTINGS\n\n";

if ($paired_ends_type eq "SE")     {&print_both("  Single-end analysis            \tYES\n");}
if ($paired_ends_type eq "PE")     {&print_both("  Paired-end analysis            \tYES\n");}
if ($remove_duplicates eq "yes")   {&print_both("  Remove duplicates              \tYES\n");}
if ($remove_duplicates eq "no")    {&print_both("  Remove duplicates              \tNO\n");}

if ($annotate_variants eq "yes")   {&print_both("  Annotate variants              \tYES using $effect_predictor\n");}
if ($annotate_variants eq "no")    {&print_both("  Annotate variants              \tNO\n");}

									&print_both("  BWA alignment option:          \t$bwa_alignment_method\n");
									&print_both("  Quality scores:                \t$quality_scores\n");
									&print_both("  Delete intermediate files:     \t$delete_intermediate_files\n\n");


&print_both("  Memory setting for Java:                                   \t$memory_string\n");
&print_both("  No of CPU threads bwa:                                     \t$no_of_threads_bwa\n");
&print_both("  No of CPU threads (-nct) GATK Variant calling:             \t$no_of_threads_gatk_ug_nct\n");
if ($variant_caller eq "UnifiedGenotyper")
{
	&print_both("  No of data threads (-nt) GATK Variant calling:             \t$no_of_threads_gatk_ug_nt\n");
}
&print_both("  No of data threads (-nt) GATK RealignerTargetCreator:      \t$no_of_threads_gatk_rtc_nt\n");
&print_both("  No of CPU threads (-nct) GATK BaseRecalibrator:            \t$no_of_threads_gatk_br_nct\n");
&print_both("  No of CPU threads (-nct) GATK PrintReads:  	              \t$no_of_threads_gatk_pr_nct\n\n");


&print_both("  Variant caller:                 \tGATK $variant_caller\n\n");
&print_both("  GATK variant caller constants:\n\n");

&print_both("    --stand_emit_conf:\t\t\t$stand_emit_conf\n");
&print_both("    --stand_call_conf:\t\t\t$stand_call_conf\n");
&print_both("    --max_alt_alleles:\t\t\t$max_alt_alleles\n");

&print_both("  Make VCF file:                  \t$make_VCF_file\n\n");

&print_message("DATA FILES","message");
print COMMAND_LOG "\nDATA FILES\n\n";

&print_both("  Reference sequence:            \t$ref\n\n");

if ($use_defined_region eq "yes")
{
	&print_both("  Defined region:                  \t$region\n\n");
	&print_both("  GATK region string:              \t$GATK_region_string\n");
}
if ($use_defined_region eq "no")
{
	&print_both("  Defined region:\t\t\t\tNO\n\n");
}

&print_both("  SNP file for BaseRecalibrator:   \t$dbsnp_file\n\n");

&print_both("  Number of FASTQ files to analyse:\t$no_of_files\n\n");
	
for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	if ($paired_ends_type eq "SE")
	{
		&print_both("    FASTQ file $list_count	\t$fastq1_file_array[$list_count]\n");
	}
	if ($paired_ends_type eq "PE")
	{
		&print_both("    Pair of FASTQ files $list_count	\t$fastq1_file_array[$list_count]\t$fastq2_file_array[$list_count]\t$sample_name_array[$list_count]\n");
	}
}

&print_message("  >> CHECK THESE SETTINGS.  If OK press ENTER to proceed with the analysis run (or Q to quit)","message");

$answer = <STDIN>;chomp $answer; 
if (lc $answer eq "q"){exit;} 
	
##########################################################
# Create some file names                                 # 
##########################################################
$final_VCF_file = "$run_title"."_variants.vcf";
$vep_out_file = "$run_title"."_variants_vep.out";
$vep_html_file = "$run_title"."_variants_vep.out_summary.html";
$final_snpEff_VCF_file = "$run_title"."_variants_snpEff.vcf";



print COMMAND_LOG "\n\n";
print COMMAND_LOG "Output VCF file:         \t$final_VCF_file\n\n";

##########################################################
# Create results_all folder (holds the results files)    #
##########################################################
&run_unix_command("mkdir $results_all_folder","Make results folder");	


##############################################################################
##############################################################################
# Start MAIN LOOP.  This runs through each of the files in the list of files #
##############################################################################
##############################################################################

$last_current_time = time();
&log_time($last_current_time,"Start of whole fastq2vcf run");

for ($loop_count=1;$loop_count <=$no_of_files;$loop_count++)
{
	##########################################################
	# Record elapsed time at the start of each loop          #
	##########################################################
	$start_loop_time = time();
	$start_raw_time = time();


	############################################################################
	# Reset error warnings - so that e-mail warnings only happen once per loop #
	############################################################################
	$unix_error_reported = "no";
	$zero_size_reported = "no";


	##########################################################
	# Get names for FASTQ reads files from arrays            #
	##########################################################
	$fastq_file_1 = $fastq1_file_array[$loop_count];
	$fastq_file_2 = $fastq2_file_array[$loop_count];
	
	

	##########################################################
	# Create a sample name based on the two FASTQ file names #
	##########################################################
	if ($sample_name_array[$loop_count] eq "")
	{
		for ($pos = 0; $pos <= length ($fastq_file_1); $pos++)
		{
			$char_1 = substr ($fastq_file_1,$pos,1);
			
			if ($char_1 eq "1")
			{
				# Does same position in second file have a "2"?
				
				$char_2 = substr ($fastq_file_2,$pos,1);
				
				if ($char_2 eq "2")
				{
					$back_to_one = substr($fastq_file_2,0,$pos)."1".substr($fastq_file_2,$pos+1,99);
					
					if ($fastq_file_1 eq $back_to_one)
					{
						$sample_name = substr($fastq_file_2,0,$pos).substr($fastq_file_2,$pos+1,99);
						
						print "FASTQ files: $fastq_file_1    \t $fastq_file_2\t";
						$sample_name = &get_prefix($sample_name);
						$back_to_one = length($sample_name);

						if (index($sample_name,"_") == length($sample_name) - 1)
						{
							$sample_name = substr($sample_name,0,length($sample_name) - 1)
						}
						print "$sample_name\n\n";
					}
				} # if
			}
		} # guess 2nd file name

		$sample_name_array[$loop_count] = $sample_name;

	} # If no third column for the sample name



	#############################################################################
	# If a sample name has been specified in the third column of the input file #	
	#############################################################################
	$sample_name = &get_prefix ($fastq_file_1); # default if none has been specified

	if ($sample_name_array[$loop_count] ne "")
	{
		$sample_name = $sample_name_array[$loop_count];
	}
	
	
	&print_message("Loop count: $loop_count  Sample name: $sample_name","message");
	

	######################################################################
	# Set up various File names which vary for each loop                 #
	######################################################################
	$aln_sa_sai = "$run_title"."_"."$sample_name"."_aln_sa.sai";
	$aln_sa1_sai = "$run_title"."_"."$sample_name"."_aln_sa1.sai";
	$aln_sa2_sai = "$run_title"."_"."$sample_name"."_aln_sa2.sai";

	$aligned_sam = "$run_title"."_"."$sample_name"."_aligned.sam";

	
	
	$final_raw_bam = "$run_title"."_"."$sample_name"."_raw.bam"; # Final BAM of FASTQ2BAM (change to 'raw'?)
	$final_raw_bai = "$run_title"."_"."$sample_name"."_raw.bai";
	
	$validate_raw_out = "$run_title"."_"."$sample_name"."_validate_raw.out";
	$flagstat_raw_out = "$run_title"."_"."$sample_name"."_flagstat_raw.out";

	# After AddReadGroups
	$aligned_sorted_rg_bam = "$run_title"."_"."$sample_name"."_aligned_sorted_rg.bam";
	$aligned_sorted_rg_bai = "$run_title"."_"."$sample_name"."_aligned_sorted_rg.bai";
	
	# After MarkDuplicates (deduped)
	$runtitle_dedup_bam = "$run_title"."_"."$sample_name".".dedup.bam";
	$runtitle_dedup_bai = "$run_title"."_"."$sample_name".".dedup.bai";
	
	# After RealignerTargetCreator (output file from this stage)
	$runtitle_dedup_bam_intervals = "$run_title"."_"."$sample_name".".dedup.bam.intervals";
	
	# After IndelRealigner (cleaned)
	$runtitle_clean_dedup_bam = "$run_title"."_"."$sample_name".".clean.dedup.bam";
	$runtitle_clean_dedup_bai = "$run_title"."_"."$sample_name".".clean.dedup.bai";
	
	# After BaseRecalibrator (output file from this stage)
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
	$insert_size_out = "$sample_name"."_insert_size.out";
	$WGS_metrics_file = "$sample_name"."_WGS_metrics.out";
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

			&run_unix_command("$bwa_path $bwa_alignment_method $ref $fastq_file_1 > $aln_sa_sai","Step 1 Run bwa to make SAI file");

			&record_input_file_size("$fastq_file_1");
			&record_output_file_size("$aln_sa_sai");

			&print_message("File $loop_count/$no_of_files   Step 1 COMPLETED Running BWA to produce aligned BAM file - SE","message");

			&test_mode_subroutine;

		} # bwa aln

		if ($bwa_alignment_method eq "mem")
		{

			&print_message("File $loop_count/$no_of_files   Step 1 Running BWA mem to produce aligned SAM file - SE","message");

			&run_unix_command("$bwa_path mem -M -t $no_of_threads_bwa $ref $fastq_file_1 > $aligned_sam","Step 1 Converting FASTQ to SAM using bwa mem");


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
			################
			# FASTQ file 1 #
			################
			
			&print_message("File $loop_count/$no_of_files   Step 1.1 Running BWA ALN to produce aligned BAM file - PE1","message");

			&run_unix_command("$bwa_path aln -t $no_of_threads_bwa $ref $fastq_file_1 > $aln_sa1_sai","Step 1.1 Run bwa to make SAI file 1");

			&record_input_file_size("$fastq_file_1");
			&record_output_file_size("$aln_sa1_sai");
			
			&print_message("File $loop_count/$no_of_files   Step 1.1 COMPLETED Running BWA ALN to produce aligned BAM file - PE1","message");

			&test_mode_subroutine;
			
			
			################
			# FASTQ file 2 #
			################
			&print_message("File $loop_count/$no_of_files   Step 1.2 Running BWA ALN to produce aligned BAM file - PE2","message");
			
			&run_unix_command("$bwa_path aln -t $no_of_threads_bwa $ref $fastq_file_2 > $aln_sa2_sai","Step 1.2 Run bwa to make SAI file 2");
			
			&record_input_file_size("$fastq_file_2");
			&record_output_file_size("$aln_sa2_sai");

			&print_message("File $loop_count/$no_of_files   Step 1.2 COMPLETED Running BWA ALN to produce aligned BAM file - PE2","message");

			&test_mode_subroutine;
		} #bwa aln

		if ($bwa_alignment_method eq "mem")
		{
			################################################
			# Read 1 and Read 2 are done together with mem #
			################################################
			&print_message("File $loop_count/$no_of_files   Step 1 Running BWA MEM to produce SAM file - PE","message");
			
			&run_unix_command("$bwa_path mem -M -t $no_of_threads_bwa $ref $fastq_file_1 $fastq_file_2 > $aligned_sam","Step 1 Converting FASTQ to SAM using bwa mem");

			&record_input_file_size("$fastq_file_1");
			&record_input_file_size("$fastq_file_2");
			&record_output_file_size("$aligned_sam");
			
			&print_message("File $loop_count/$no_of_files   Step 1 COMPLETED Running BWA MEM to produce SAM file - PE","message");

			&test_mode_subroutine;	
		} #bwa mem
	} # PE


	##########################################################################
	# 2 Convert SAI files into SAM files  (aln only, not required for mem)   #
	##########################################################################
	if ($bwa_alignment_method eq "aln")
	{
		if  ($paired_ends_type eq "SE")
		{
			&print_message("File $loop_count/$no_of_files   Step 2 Using SAI file to create SAM file","message");

			&run_unix_command("$bwa_path samse $ref $aln_sa_sai $fastq_file_1 > $aligned_sam","Step 2 bwa: Converting SAI to SAM");

			&record_input_file_size("$aln_sa_sai");
			&record_output_file_size("$aligned_sam");

			&test_mode_subroutine;
		
			&print_message("File $loop_count/$no_of_files  Step 2 COMPLETED Using SAI file to create SAM file","message");
		}

		if  ($paired_ends_type eq "PE")
		{
			&print_message("File $loop_count/$no_of_files   Step 2 Using SAI files to create SAM file","message");

			&run_unix_command("$bwa_path sampe $ref $aln_sa1_sai $aln_sa2_sai $fastq_file_1 $fastq_file_2 > $aligned_sam","Step 2 bwa: Converting SAI to SAM");

			&record_input_file_size("$aln_sa1_sai");
			&record_input_file_size("$aln_sa2_sai");
			&record_output_file_size("$aligned_sam");
			
			&test_mode_subroutine;
		
			&print_message("File $loop_count/$no_of_files   Step 2 COMPLETED Using SAI files to create SAM file","message");
		}

	} # aln only

    print " * * * \n";


	##########################################################
	# Run all the picard stages in one                       #
	#  - Sort the SAM file generated by BWA                  #
	#  - Convert the SAM file (from BWA) to a BAM file       #
	#  - Add Read Group information to the BAM file          #	
	##########################################################

	&print_message("File $loop_count/$no_of_files   Step 3 Picard stages","message");
	print "  - ReadGroup info to be added to the BAM file\n";
	print "  - The aligned reads to be sorted\n";
	print "  - The file to be converted to BAM format\n\n";

	&run_unix_command("java $memory_string $temp_dir_string -jar $picard_path/AddOrReplaceReadGroups.jar I=$aligned_sam O=$aligned_sorted_rg_bam rgID=$run_title LB=$lib PL='ILLUMINA' PU=$sample_name SM=$sample_name SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=$picard_validation_stringency","Step 3 Picard stages to convert SAM to BAM file");	
	
	&record_input_file_size("$aligned_sam");
	&record_output_file_size("$aligned_sorted_rg_bam");

	&print_message("File $loop_count/$no_of_files   Step 3 Picard stages COMPLETED","message");
	

	################################################
	# Only delete SAM file if next file is made OK #
	################################################
	&delete_first_file_if_second_ok("$aligned_sam","$aligned_sorted_rg_bam");


	################################################
	# Step 4: Use picard ValidateSamFile           #
	################################################
	&print_message("File $loop_count/$no_of_files   Step 4 Picard ValidateSamFile","message");
	
	&run_unix_command("java $memory_string $temp_dir_string -jar $picard_path/ValidateSamFile.jar I=$aligned_sorted_rg_bam O=$validate_raw_out MODE=SUMMARY MAX_OPEN_TEMP_FILES=$max_open_temp_files VALIDATION_STRINGENCY=SILENT","Step 4 Picard ValidateSamFile");

	&record_input_file_size("$aligned_sorted_rg_bam");
	&record_output_file_size("$validate_raw_out");

	&print_message("File $loop_count/$no_of_files   Step 4 Picard ValidateSamFile COMPLETED","message");

	
	################################################
	# Step 5: Use samtools flagstat                #
	################################################
	&print_message("File $loop_count/$no_of_files   Step 5 samtools flagstat","message");
	
	&run_unix_command("$samtools_path flagstat $aligned_sorted_rg_bam > $flagstat_raw_out","Step 5 samtools flagstat");

	&record_input_file_size("$aligned_sorted_rg_bam");
	&record_output_file_size("$flagstat_raw_out");

	&print_message("File $loop_count/$no_of_files   Step 5 samtools flagstat COMPLETED","message");

	
	################################################
	# Rename as final_raw.bam  (for next stage)    #
	################################################
	&record_input_file_size("$aligned_sorted_rg_bam");
	&rename_file("$aligned_sorted_rg_bam","$final_raw_bam");
	&record_output_file_size("$final_raw_bam");

	&record_input_file_size("$aligned_sorted_rg_bai");
	&rename_file("$aligned_sorted_rg_bai","$final_raw_bai");
	&record_output_file_size("$final_raw_bai");
	
	print "\n\n";
	
	################################################
	# Write info to the README file for FASTQ2BAM  #
	################################################

	$readme_file =  "$run_title"."_"."$sample_name"."_raw_readme.out";

	open (READMEFILE, ">$readme_file"); 

	print READMEFILE "#######################################################\n";
	print READMEFILE "Summary of fastq2vcf results files for this sample     \n";
	print READMEFILE "#######################################################\n\n";

	print READMEFILE "Program:\t\t\tfastq2vcf (fastq2bam section)\tStared: $version\n\n";
		
	print READMEFILE "Run title:\t\t\t$run_title\n\n";
	
	print READMEFILE "Reference sequence:\t$ref_seq_name\n\n";
	
	print READMEFILE "FASTQ file:        \t$fastq_file_1\n\n\n";

	print READMEFILE "BAM ALIGNMENT FILE (Raw alignments to the reference just using bwa)\n\n";

	print READMEFILE "\t$final_raw_bam   (Note: raw means this BAM file comes from the fastq2bam stage)\n\n\n";
	
	print READMEFILE "\tThis BAM file has been deleted as an intermediate file, but to check if it was OK look at these two files:\n\n";
	
	print READMEFILE "\t\t$validate_raw_out\n\n";
	print READMEFILE "\t\t$flagstat_raw_out\n\n\n";
	
	print READMEFILE "This BAM file then went through the rest of the NGS pipeline using the GATK pipeline\n\n";

	print READMEFILE "FINAL BAM ALIGNMENT FILES AFTER PROCESSING BY GATK\n\n";

	print READMEFILE "\t$final_bam		     \tBest alignments to the reference after processing by GATK\n\n\n";
	
	print READMEFILE "\tLook at these files to check the quality of your final BAM file:\n\n";

	print READMEFILE "\t\t$validate_final_out\n\n";
	print READMEFILE "\t\t$flagstat_final_out\n\n\n";

	close (READMEFILE);

	&move_to_results_all_folder ("$readme_file");


	##########################################################
	# Record time for making the raw BAM files               #
	##########################################################
	$current_time = time();
	$make_raw_bam_time = $current_time - $start_raw_time;

	$make_raw_bam_time_total = $make_raw_bam_time_total + $make_raw_bam_time;

	print COMMAND_LOG "\n  Record time to make raw BAM file\n\n";
	print COMMAND_TIMES_LOG "\n  Record time to make raw BAM file\n\n";

	&log_time($start_loop_time,"Start time for loop $loop_count");
	&log_time($make_raw_bam_time,"Time to make single raw BAM file $loop_count");
	&log_time($make_raw_bam_time_total,"Total time to make all raw BAM files");
	&log_time($current_time,"Current time");


	##########################################################
	# Record time at start of making clean BAM files         #
	##########################################################
	$start_clean_time = time();

	##########################################################
	# Step 6: Use picard/MarkDuplicates to mark duplicates   #
	##########################################################
	if ($remove_duplicates eq "yes")
	{
		&print_message("File $loop_count/$no_of_files   Step 6  Use picard/MarkDuplicates to mark duplicates","message");

		&run_unix_command ("$javabin $memory_string $temp_dir_string -jar $picard_path/MarkDuplicates.jar I=$final_raw_bam O=$runtitle_dedup_bam  VALIDATION_STRINGENCY=$picard_validation_stringency CREATE_INDEX=true M=$runtitle_sample_metrics","Step 6 picard: Mark Duplicates");
		&record_input_file_size ("$final_raw_bam");
		&record_output_file_size ("$runtitle_dedup_bam");
		
		&delete_first_file_if_second_ok("$final_raw_bam","$runtitle_dedup_bam");

		&print_message("File $loop_count/$no_of_files   COMPLETED Step 6  Use picard/MarkDuplicates to mark duplicates","message");
	}


	#&delete_first_file_if_second_ok("$aligned_sam","$aligned_sorted_rg_bam"); This happens earlier


	##############################################################################################
	# If you didn't have the remove duplicates stage, the output file won't be called            #
	# runtitle_dedup.bam so we have to put the name of the final_raw_bam into runtitle_dedup.bam #
	##############################################################################################                                 
	if ($remove_duplicates eq "no")
	{
		$runtitle_dedup_bam = $final_raw_bam;

		&print_both("Setting run_title_dedup_bam as final_raw_bam\n");
		&record_output_file_size ("$runtitle_dedup_bam");
	}
	print "\n* * * \n";


	####################################################################################################
	# Step 7: Use GenomeAnalysisTK/RealignerTargetCreator to identify intervals in need of realignment #
	####################################################################################################
	&print_message("File $loop_count/$no_of_files   Step 7:  Use GenomeAnalysisTK/RealignerTargetCreator to identify intervals in need of realignment","message");
	&run_unix_command ("$javabin $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T RealignerTargetCreator -I $runtitle_dedup_bam $GATK_region_string -nt $no_of_threads_gatk_rtc_nt -R $ref -o $runtitle_dedup_bam_intervals -mismatch 0.0 -S $GATK_validation_stringency $GATK_fix_quals_string $filter_mismatching_string","Step 7 RealignerTargetCreator");
	&record_input_file_size ("$runtitle_dedup_bam");
	&record_output_file_size ("$runtitle_dedup_bam_intervals");

	&print_message("File $loop_count/$no_of_files   COMPLETED Step 7:  Use GenomeAnalysisTK/RealignerTargetCreator to identify intervals in need of realignment","message");

	print "\n* * * \n";
		

	####################################################################################################
	# Step 8:  Use GenomeAnalysisTK/IndelRealigner to Realign Indels                                   #
	####################################################################################################
	&print_message("File $loop_count/$no_of_files   Step 8: Using GenomeAnalysisTK/IndelRealigner to Realign Indels","message");
	&run_unix_command ("$javabin $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T IndelRealigner $GATK_region_string -R $ref -I $runtitle_dedup_bam -targetIntervals $runtitle_dedup_bam_intervals -o $runtitle_clean_dedup_bam -S $GATK_validation_stringency $GATK_fix_quals_string $filter_mismatching_string","Step 8 IndelRealigner");

	&record_input_file_size ("$runtitle_dedup_bam");
	&record_input_file_size ("$runtitle_dedup_bam_intervals");
	&record_output_file_size ("$runtitle_clean_dedup_bam");

	&delete_first_file_if_second_ok("$runtitle_dedup_bam","$runtitle_clean_dedup_bam");
	&delete_first_file_if_second_ok("$runtitle_dedup_bam_intervals","$runtitle_clean_dedup_bam");

	&print_message("File $loop_count/$no_of_files   COMPLETED Step 8: Using GenomeAnalysisTK/IndelRealigner to Realign Indels","message");

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


	####################################################################################################
	# Step 9:  NEW: Use GenomeAnalysisTK/BaseRecalibrator to generate a recalibration table            #
	####################################################################################################
	&print_message("File $loop_count/$no_of_files   Step 9: Using GenomeAnalysisTK/BaseRecalibrator to generate a Recalibration Table","message");

	&run_unix_command ("$javabin $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T BaseRecalibrator $GATK_region_string -nct $no_of_threads_gatk_br_nct -R $ref -I $runtitle_clean_dedup_bam $GATK_known_sites_string -o $runtitle_sample_recal_grp -S $GATK_validation_stringency","Step 9 BaseRecalibrator");

	&record_input_file_size ("$runtitle_clean_dedup_bam");
	&record_output_file_size ("$runtitle_sample_recal_grp");

	&print_message("File $loop_count/$no_of_files   COMPLETED Step 9: Using GenomeAnalysisTK/BaseRecalibrator to generate a Recalibration Table","message");

	print "\n* * * \n";


	####################################################################################################
	# Step 10:  Use GenomeAnalysisTK/PrintReads to update the base quality scores                      #
	####################################################################################################
	&print_message("File $loop_count/$no_of_files   Step 10:  Using GenomeAnalysisTK/PrintReads to update the base quality scores","message");

	&run_unix_command ("$javabin $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T PrintReads $GATK_region_string -I $runtitle_clean_dedup_bam  -nct $no_of_threads_gatk_pr_nct -R $ref  -BQSR $runtitle_sample_recal_grp  -o $runtitle_clean_dedup_recal_bam -S $GATK_validation_stringency","Step 10 PrintReads");

	&record_input_file_size ("$runtitle_clean_dedup_bam");
	&record_output_file_size ("$runtitle_clean_dedup_recal_bam");

	&delete_first_file_if_second_ok("$runtitle_clean_dedup_bam","$runtitle_clean_dedup_recal_bam");
	&delete_first_file_if_second_ok("$runtitle_sample_recal_grp","$runtitle_clean_dedup_recal_bam"); # recal_grp file not needed any more if recal_bam has been made

	&print_message("File $loop_count/$no_of_files   COMPLETED Step 10:  Using GenomeAnalysisTK/PrintReads to update the base quality scores","message");


	####################################################################################################
	# Step 11: Use picard/ValidateSamFile to validate the cleaned BAM file (Summary only)              #
	####################################################################################################
	if ($validate_bam_files eq "yes")
	{
		&print_message("File $loop_count/$no_of_files   Step 11:  Using picard/ValidateSamFile to validate the cleaned BAM file","message");

		&run_unix_command ("$javabin $memory_string $temp_dir_string -jar $picard_path/ValidateSamFile.jar INPUT=$runtitle_clean_dedup_recal_bam OUTPUT=$validate_final_out VALIDATION_STRINGENCY=$picard_validation_stringency MODE=SUMMARY MAX_OUTPUT=100 MAX_OPEN_TEMP_FILES=$max_open_temp_files","Step 11 Picard ValidateSamFile on final.bam");

		&record_input_file_size ("$runtitle_clean_dedup_recal_bam");
		&record_output_file_size_no_mail ("$validate_final_out");
		
		&move_to_results_all_folder ("$validate_final_out");
		
		&print_message("File $loop_count/$no_of_files   COMPLETED Step 11: Using picard/ValidateSamFile to validate the cleaned BAM file","message");

		print "\n* * * \n";	


		####################################################################################################
		# Step 12: Use samtools flagstat to get simple stats on bam file                                   #
		####################################################################################################
		&print_message("File $loop_count/$no_of_files   Step 12: samtools flagstat to be carried out","message");

		&run_unix_command("$samtools_path flagstat $runtitle_clean_dedup_recal_bam > $flagstat_final_out","Step 12 samtools flagstat on final.bam");

		&record_input_file_size ("$runtitle_clean_dedup_recal_bam");
		&record_output_file_size_no_mail("$flagstat_final_out");

		&move_to_results_all_folder ("$flagstat_final_out");
			
		&print_message("File $loop_count/$no_of_files   COMPLETED Step 12: samtools flagstat to be carried out","message");

	} # if validate_bam_files

	print "\n* * * \n";


	##############################################################################################################
	# Step 13: If you are using a defined region then create a new smaller BAM file with only this region in it  #
	##############################################################################################################
	if ($use_defined_region eq "yes")
	{
		&print_message("File $loop_count/$no_of_files   Step 13:  Making Bam file for region $region","message");

		&run_unix_command("$samtools_path view $runtitle_clean_dedup_recal_bam $region -b -o $region_only_bam","Step 13 Make smaller BAM file to specified region");

		&record_input_file_size ("$runtitle_clean_dedup_recal_bam");
		&record_output_file_size ("$region_only_bam");	
		
		&print_message("File $loop_count/$no_of_files   COMPLETED Step 13:  Making Bam file for region $region","message");
		
		############################################################################
		# Step 14:  Now make an index file for this new smaller BAM file           #
		############################################################################
		&print_message("File $loop_count/$no_of_files   Step 14:  New Index for region-only BAM file being created","message");

		&run_unix_command("java $memory_string $temp_dir_string -jar $picard_path/BuildBamIndex.jar I=$region_only_bam O=$region_only_bai VALIDATION_STRINGENCY=LENIENT","Step 14 Make BAM index of smaller BAM file");

		&record_input_file_size ("$runtitle_clean_dedup_recal_bam");
		&record_output_file_size ("$region_only_bai");
		
		&delete_first_file_if_second_ok("$runtitle_clean_dedup_recal_bam","$region_only_bam");

		&print_message("File $loop_count/$no_of_files   COMPLETED Step 14:  New Index for region-only BAM file being created","message");
		
	} # if ($use_defined_region eq "yes")


	##############################################################
	# If you didn't use a defined region, then change the final  #
	# BAM file from runtitle_clean_dedup_recal_bam to final_bam  #
	##############################################################
	if ($use_defined_region eq "no")
	{
		&record_input_file_size ("$runtitle_clean_dedup_recal_bam");
		&run_unix_command("mv $runtitle_clean_dedup_recal_bam $final_bam","Change file name: $runtitle_clean_dedup_recal_bam to $final_bam");
		&record_output_file_size ("$final_bam");
		
		&record_input_file_size ("$runtitle_clean_dedup_recal_bai");
		&run_unix_command("mv $runtitle_clean_dedup_recal_bai $final_bai","Change file name: $runtitle_clean_dedup_recal_bai to $final_bai");
		&record_output_file_size ("$final_bai");

	} # if ($use_defined_region eq "no")


	##############################################################
	# If you DID use a defined region, then change the final     #
	# BAM file from region_only_bam to final_bam                 #
	##############################################################
	if ($use_defined_region eq "yes")
	{
		&record_input_file_size ("$region_only_bam");
		&run_unix_command("mv $region_only_bam $final_bam","Change file name: $region_only_bam to $final_bam");
		&record_output_file_size ("$final_bam");
		
		&record_input_file_size ("$region_only_bai");
		&run_unix_command("mv $region_only_bai $final_bai","Change file name: $region_only_bai to $final_bai");
		&record_output_file_size ("$final_bai");
	} # if ($use_defined_region eq "yes")

	print "\n* * * \n";


	####################################################
	# Record time for making the final clean BAM files #
	####################################################
	$current_time = time();
	$make_clean_bam_time = $current_time - $start_clean_time;
	$make_clean_bam_time_total = $make_clean_bam_time_total + $make_clean_bam_time; # if there are many BAM files

	print COMMAND_LOG "\n  Record time to make final clean BAM file\n\n";
	print COMMAND_TIMES_LOG "\n  Record time to make final clean BAM file\n\n";

	&log_time($make_clean_bam_time,"Time to clean BAM file $loop_count  (from 'raw' to 'final')");
	&log_time($make_clean_bam_time_total,"Total time to clean all BAM files (from 'raw' to 'final')");
	&log_time($current_time,"Current time");


	##############################################################################################
	# Create input string for parallel use of the UnifiedGenotyper or HaplotypeCaller at the end #
	##############################################################################################
	if (-e $final_bam)
	{
		$variant_caller_input_string = $variant_caller_input_string." -I $results_all_folder/$final_bam";
	}
	if (! -e $final_bam)
	{
		$bam_missing_for_variant_calling = $bam_missing_for_variant_calling + 1;
		# Maybe have gVCF missing at some stage?
	}


	#########################################################################################
	# If the variant caller is the gVCF option then you can do this now (inside main looop) #
	# (if you use HaplotypeCaller (not gVCF) then you only do variant calling at the end)   #
	#########################################################################################
	if ($variant_caller eq "HaplotypeCallerGVCF")
	{
		$start_gVCF_time= time();

		##########################################################################
		# First make the input string for later                                  #
		##########################################################################
		$output_gVCF_file = "$run_title"."_"."$sample_name".".gVCF";
		$output_gVCF_file_idx= "$run_title"."_"."$sample_name".".gVCF.tbi";
		if ($bgzip_output_gVCFs eq "true")
		{
			$output_gVCF_file = $output_gVCF_file.".gz";
			$output_gVCF_file_idx = $output_gVCF_file.".tbi";
		}

		$GVCF_input_string = $GVCF_input_string." -V $results_all_folder/$output_gVCF_file"; # for use outside the main loop

		##########################################################################
		# Now make the actual gVCF file (which doesn't have to be done using all #
		# the BAM files together, so can be done inside the Ma in Loop)          #
		##########################################################################

		&print_message("File $loop_count/$no_of_files   Step 14a: Making gVCF file....","message");


		##########################################################################
		# stand_emit_conf and stand_call_conf are not used when HaplotypeCaller  #
		# is used in emitRefConfidence mode. You use them when recombining with  #
		# GenotypeGVCFs                                                          #
		##########################################################################
		&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 $pcr_free_string $GATK_region_string  -I $final_bam -o $output_gVCF_file -S $GATK_validation_stringency -nct $no_of_threads_gatk_ug_nct","Step 14a Run HaplotypeCaller in GVCF mode");	


		&record_input_file_size ("$final_bam"); # At this stage final_bam isn't in results_all_folder
		&record_output_file_size_no_mail("$output_gVCF_file");

		##############################################
		# Move gVCF file to results folder           #
		##############################################
		&move_to_results_all_folder ("$output_gVCF_file");
		&move_to_results_all_folder ("$output_gVCF_file_idx");

		&print_message("File $loop_count/$no_of_files   Step 14a: gVCF file has been made....","message");

		#################################
		# Record time after gVCF stage  #
		#################################
		$current_time = time();
		$make_gVCF_time = $current_time - $start_gVCF_time;

		$make_gVCF_time_total = $make_gVCF_time_total + $make_gVCF_time;

		print COMMAND_LOG "\n  Record time to make gVCF file\n\n";
		print COMMAND_TIMES_LOG "\n  Record time to make gVCF file\n\n";

		&log_time($make_gVCF_time,"Time to make gVCF file $loop_count");
		&log_time($make_gVCF_time_total,"Total time to make all gVCF files");
		&log_time($current_time,"Current time");


	} # if ($variant_caller eq "HaplotypeCallerGVCF")


	##################################################
	# Step 15:  Plot Insert Size Histogram           #
	##################################################
	&print_message("File $loop_count/$no_of_files   Step 15: Plotting insert size histogram...","message");

	&run_unix_command("java $memory_string $temp_dir_string -jar $picard_path/CollectInsertSizeMetrics.jar INPUT=$final_bam OUTPUT=$insert_size_out HISTOGRAM_FILE=$insert_size_pdf VALIDATION_STRINGENCY=LENIENT","Step 15 Plot insert size");

	&record_input_file_size ("$final_bam");
	&record_output_file_size_no_mail("$insert_size_pdf");
	&record_output_file_size_no_mail("$insert_size_out");

	&print_message("File $loop_count/$no_of_files   Step 15: Insert size histogram plotted","message");

	
#java -Xmx45g  -Djava.io.tmpdir=/home/LANPARK/mboursnell/javatempdir -jar $picard_path/CollectWgsMetrics.jar INPUT=ESS_ESS_26602_case_final.bam  OUTPUT=WGS_metrics.out  REFERENCE_SEQUENCE=/home/genetics/canfam3/canfam3.fasta STOP_AFTER=100000000   VALIDATION_STRINGENCY=LENIENT

	##################################################
	# Step 16: WGS metrics using CollectWgsMetrics   #
	##################################################
	&print_message("File $loop_count/$no_of_files   Step 16: Calculating WGS Metrics...","message");

	&run_unix_command("java $memory_string $temp_dir_string -jar $picard_path/CollectWgsMetrics.jar INPUT=$final_bam OUTPUT=$WGS_metrics_file REFERENCE_SEQUENCE=$ref STOP_AFTER=100000000 VALIDATION_STRINGENCY=LENIENT","Step 15 Plot insert size");

	&record_input_file_size ("$final_bam");
	&record_output_file_size_no_mail("$WGS_metrics_file");

	&print_message("File $loop_count/$no_of_files   Step 16: IWGS metrics calculated","message");


	##################################################
	# Move various files to the Results_All Folder   #
	##################################################
	&print_message("Moving various files","message");

	&print_both("------------------------------------------------------------------------------------------------------------------------------------\n");
	&print_both("Moving various files\n");


	&move_to_results_all_folder ("$validate_raw_out");
	&move_to_results_all_folder ("$flagstat_raw_out");
	
	&move_to_results_all_folder ("$final_raw_bam"); # will be deleted later on
	&move_to_results_all_folder ("$final_raw_bai"); # will be deleted later on

	&move_to_results_all_folder ("$final_bam");
	&move_to_results_all_folder ("$final_bai");
	&move_to_results_all_folder ("$runtitle_sample_metrics"); # MarkDuplicates info
	
	&move_to_results_all_folder ("$insert_size_pdf");
	&move_to_results_all_folder ("$insert_size_out");
	#&move_to_results_all_folder ("$WGS_metrics_file");
	

	###############################################################
	# For first loop only, save all files (for checking purposes) #
	# (only if $keep_all_for_first_file = "yes")                  #
	###############################################################
	if (($loop_count == 1) && ($keep_all_for_first_file eq "yes"))
	{
		&move_to_results_all_folder ("$rg_bai");
		&move_to_results_all_folder ("$runtitle_dedup_bam");
		&move_to_results_all_folder ("$runtitle_dedup_bai");
		&move_to_results_all_folder ("$runtitle_clean_dedup_bam");
		&move_to_results_all_folder ("$runtitle_clean_dedup_bai");
		&move_to_results_all_folder ("$runtitle_clean_dedup_recal_bam");
		&move_to_results_all_folder ("$runtitle_clean_dedup_recal_bai");
		&move_to_results_all_folder ("$runtitle_dedup_bam_intervals");
	}
	

	##################################################
	# Some files don't need to be kept               #
	# (keep files for first loop to check)           #
	##################################################
	
	if ($loop_count == 1){$delete_intermediate_files = "no"} else {$delete_intermediate_files = "yes"}

	if ($delete_intermediates_option eq "no"){$delete_intermediate_files = "no"}

	if ($keep_all_for_first_file eq "no"){$delete_intermediate_files = "yes"}
	
	if ($delete_intermediate_files eq "no")
	{
		print COMMAND_LOG "   Delete intermediate files is turned off\n\n";
	}

	if ($delete_intermediate_files eq "yes")
	{
		#####################################
		# Delete files from bwa aln stage   #
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

	
 		##################################################
		# Delete files from bam2vcf stage                #
		# If they haven't been deleted already           #
		##################################################

		#&delete_file("$aligned_sam"); # should have been deleted already
		#&delete_file("$rg_bai");
		#&delete_file("$runtitle_dedup_bam");
		#&delete_file("$runtitle_dedup_bai");
		#&delete_file("$runtitle_dedup_bam_intervals");
		#&delete_file("$runtitle_clean_dedup_bam");
		#&delete_file("$runtitle_clean_dedup_bai");
		#&delete_file("$runtitle_clean_dedup_recal_bam");
		#&delete_file("$runtitle_clean_dedup_recal_bai");
		&delete_file("out1.junk");
		&delete_file("out2.junk");
		#&delete_file ("$runtitle_sample_recal_grp");


		##################################################
		# If final BAM has been produced, delete raw BAM #
		##################################################

		&delete_first_file_if_second_ok("$results_all_folder/$final_raw_bam","$results_all_folder/$final_bam");
		&delete_first_file_if_second_ok("$results_all_folder/$final_raw_bai","$results_all_folder/$final_bai");

	} #if ($delete_intermediate_files eq "yes")


	
	####################################################################
	# Calculate the elapsed total run time at this point (end of loop) #
	####################################################################
	$current_time = time();
	$run_time_so_far = $current_time - $start_time;

	##################################################
	# Check if final BAM is made                     #
	##################################################
	if (-e "$results_all_folder/$final_bam") {$filesize_main = -s "$results_all_folder/$final_bam";} else {$filesize_main = "Not found"}


	##################################################
	# Send an e-mail at the end of each loop         #
	##################################################
	open(MAIL, "|/usr/sbin/sendmail -t");

	## Mail Header
	print MAIL "To: $email_address\n";
	print MAIL "From: $e_mail_from\n";

	if ($filesize_main eq "Not found")
	{
		print MAIL "Subject: FASTQ2VCF: Sample $sample_name in run $run_title has FAILED to reach the final BAM file stage\n\n";
	}
	else
	{
		print MAIL "Subject: FASTQ2VCF: Sample $sample_name in run $run_title has reached the final BAM file stage\n\n";
	}
	
	## Mail Body
	print MAIL "Run title:     \t$run_title\n\n";
	print MAIL "Sample name:   \t$sample_name  (sample number $loop_count)\n\n";
	print MAIL "Script:        \tfastq2vcf PERL script from GIT started $version\n\n";

	print MAIL "FASTQ files:\t$fastq_file_1\t$fastq_file_2\n";
	print MAIL "BAM file:   \t$final_bam\n\n";

	if ($filesize_main eq "Not found")
	{
		print MAIL "The processing of these FASTQ files into a final BAM file has finished running,\n\n";
		print MAIL "but the final BAM file has NOT BEEN CREATED.  Check the command log and the screen log.\n\n";;
	}
	else
	{
		print MAIL "The processing of these FASTQ files into a final BAM file is complete.\n\n";
		print MAIL "The script may still be running for the variant calling stage (making VCF files).\n\n";
	}

	print MAIL  "\n  Output file: $final_bam\t\tSize: $filesize_main\n\n";
	&print_both("\n  Output file: $results_all_folder/$final_bam\t\tSize: $filesize_main\n\n");

	print MAIL "Run time so far : $run_time_so_far seconds\n";
	printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time_so_far)[7,2,1,0];
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
	$run_time_so_far = $current_time - $start_time;
	printf COMMAND_LOG "  Total time at the end of loop $loop_count:    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time_so_far)[7,2,1,0];
	printf COMMAND_LOG "  Time for loop $loop_count:                    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $loop_time)[7,2,1,0];


} # End of MA IN LOOP

###########################################################################################################
######################################### END OF MAIN LOOP ################################################
###########################################################################################################


####################################################################
# If there is only one FASTQ file then there is only one gVCF file #
# This is now the $final_gVCF_file                                 #
####################################################################
if ($no_of_files == 1)
{
	$final_gVCF_file = $output_gVCF_file;
}
##########################################################
#  Make a VCF file using all the BAM files in parallel   #
##########################################################

if ($make_VCF_file eq "true")
{
	if ($no_of_files > 1)
	{
		&print_message("Making VCF with $no_of_files BAM files in parallel.  File of file names is $list_file","message");

		print COMMAND_LOG "\n================================================================================\n";
		print COMMAND_LOG   "Making VCF with all BAM files in parallel.  File of file names is $list_file";
		print COMMAND_LOG "\n================================================================================\n\n";
	}
	if ($no_of_files == 1)
	{
		&print_message("Making VCF with single BAM file.  File of file names is $list_file","message");

		print COMMAND_LOG "\n================================================================================\n";
		print COMMAND_LOG   "Making VCF with single BAM file.  File of file names is $list_file";
		print COMMAND_LOG "\n================================================================================\n\n";
	}


	################################
	# Record time before VCF stage #
	################################
	$start_VCF_time = time();

	$variant_caller_input_string = $variant_caller_input_string." ";
	$GVCF_input_string = $GVCF_input_string." ";


	#####################################################################################
	# Variant calling to make VCF file using either HaplotypeCaller or UnifiedGenotyper #
	#####################################################################################

	if ($variant_caller eq "UnifiedGenotyper")
	{
		&print_message("Running UnifiedGenotyper on all $no_of_files files... ","message");

		if ($use_default_stand_values eq "false")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm BOTH $variant_caller_input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o $final_VCF_file -S $GATK_validation_stringency -nt $no_of_threads_gatk_ug_nt -nct $no_of_threads_gatk_ug_nct","18 Run GATK UnifiedGenotyper to call variants");	
		}
		if ($use_default_stand_values eq "true")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm BOTH $variant_caller_input_string -o $final_VCF_file -S $GATK_validation_stringency -nt $no_of_threads_gatk_ug_nt -nct $no_of_threads_gatk_ug_nct","18 Run GATK UnifiedGenotyper to call variants");	
		}
		&record_output_file_size ("$final_VCF_file");

		&print_message("Finished running UnifiedGenotyper on all $no_of_files files... ","message");

	} # UG

	if ($variant_caller eq "HaplotypeCaller")
	{

		&print_message("Running HaplotypeCaller on all $no_of_files files... ","message");

		if ($use_default_stand_values eq "false")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -pairHMM VECTOR_LOGLESS_CACHING -R $ref -T HaplotypeCaller $pcr_free_string $GATK_region_string  $variant_caller_input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o $final_VCF_file -S $GATK_validation_stringency -nct $no_of_threads_gatk_ug_nct","18 Run GATK HaplotypeCaller to call variants");	
		}
		if ($use_default_stand_values eq "true")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -pairHMM VECTOR_LOGLESS_CACHING -R $ref -T HaplotypeCaller $pcr_free_string $GATK_region_string  $variant_caller_input_string -o $final_VCF_file -S $GATK_validation_stringency  -nct $no_of_threads_gatk_ug_nct","18 Run GATK HaplotypeCaller to call variants");	
		}
		&record_output_file_size ("$final_VCF_file");

		&print_message("Finished running HaplotypeCaller on all $no_of_files files... ","message");	
	} # HC


	if ($variant_caller eq "HaplotypeCallerGVCF")
	{
		####################################################################
		# GenotypeGVCFs - run on all the files at once (outside main loop) #
		####################################################################

		&print_message("Running GenotypeGVCFs on all $no_of_files files... ","message");

		&print_both("Input string: $GVCF_input_string\n\n");

		&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T GenotypeGVCFs -nt $no_of_threads_gatk_ug_nt $GATK_region_string  $GVCF_input_string -o $final_VCF_file -S $GATK_validation_stringency","Running GenotypeGVCFs on all gVCF files to make a single VCF file");	

		&record_output_file_size("$final_VCF_file");

		&print_message("The VCF files from Haplotype_Caller gVCF mode have been merged into a single VCF file $final_VCF_file","message");
		
	} # gVCF

	$current_time = time();
	$make_VCF_time = $current_time - $start_VCF_time;
	$make_VCF_time_total = $make_VCF_time_total + $make_VCF_time;

	print COMMAND_LOG "\n  Record time to make VCF file\n\n";
	print COMMAND_TIMES_LOG "\n  Record time to make VCF file\n\n";

	&log_time($make_VCF_time,"Time to make VCF file $loop_count");
	&log_time($make_VCF_time_total,"Total time to make all VCF files");
	&log_time($current_time,"Current time");

} # if $make_VCF is true


#############################################################
# If the FASTQ files are all from the same sample then      #
# you might opt to make a GVCF file from all the BAM files  #
# together (rather than one from each BAM file)             #                          
#############################################################
if ($make_gVCF_file_outside_loop eq "true")
{
	if ($variant_caller eq "HaplotypeCallerGVCF")
	{
		####################################################################
		# GenotypeGVCFs - run on all the files at once (outside main loop) #
		####################################################################
		#need to find out how to combine BAM files into a SINGEL-COLUMN gVCF

		&print_message("Running HaplotypeCallerGVCF on all $no_of_files BAM files to make a single gVCF... ","message");

		&print_both("Input string: $GVCF_input_string\n\n");

		&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T GenotypeGVCFs -nt $no_of_threads_gatk_ug_nt $GATK_region_string  $GVCF_input_string -o $final_VCF_file -S $GATK_validation_stringency","Running GenotypeGVCFs on all gVCF files to make a single VCF file");	

		&record_output_file_size("$final_VCF_file");

		&print_message("The VCF files from Haplotype_Caller gVCF mode have been merged into a single VCF file $final_VCF_file","message");	
	}
}

#########################################################################
# Annotate SNPs and Indels with the finished VCF file (if there is one) #
#########################################################################

if ($annotate_variants eq "yes")
{
	print COMMAND_LOG "\n#######################\n";
	print COMMAND_LOG "# Annotating variants #\n";
	print COMMAND_LOG "#######################\n\n";

	if ($effect_predictor eq "vep")
	{
		if ($use_defined_region eq "yes"){&run_unix_command("perl /home/genetics/scripts/run_variant_effect_predictor.pl --file $final_VCF_file --ref $ref_seq_name --chr $chromosome","19 Run variant_effect_predictor on chr $chromosome");}

		if ($use_defined_region eq "no"){&run_unix_command("perl /home/genetics/scripts/run_variant_effect_predictor.pl --file $final_VCF_file --ref $ref_seq_name --chr none","19 Run variant_effect_predictor");}
	}
	if ($effect_predictor eq "snpEff")
	{
		if ($use_defined_region eq "yes"){&run_unix_command("perl /home/genetics/scripts/run_snpEff.pl --file $final_VCF_file --ref $ref_seq_name --chr $chromosome","19 Run snpEff on chr $chromosome");}
		if ($use_defined_region eq "no"){&run_unix_command("perl /home/genetics/scripts/run_snpEff.pl --file $final_VCF_file --ref $ref_seq_name","19 Run snpEff");}
	}
	if ($effect_predictor eq "none")
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
if ($make_VCF_file eq "true")
{
	&move_to_results_all_folder ("$final_VCF_file");
	&move_to_results_all_folder ("$final_VCF_idx_file");
	$final_VCF_size = -s "$results_all_folder/$final_VCF_file";
}

if ($effect_predictor eq "vep")
{
	&move_to_results_all_folder ("$vep_out_file");
	&move_to_results_all_folder ("$vep_html_file");
	&move_to_results_all_folder ("$vep_command_log"); # Add run_title later
	$final_vep_size = -s "$results_all_folder/$vep_out_file";
}
if ($effect_predictor eq "snpEff")
{
	&move_to_results_all_folder ("$final_snpEff_VCF_file");
	$final_snpEff_vcf_size = -s "$results_all_folder/$final_snpEff_VCF_file";
}

if ($variant_caller eq "HaplotypeCallerGVCF")
{
	$final_gVCF_size = -s "$results_all_folder/$output_gVCF_file";
}



################################################
# Record time after VCF stage outside the loop #
################################################

$current_time = time();
$run_time_so_far = $current_time - $start_time;
$total_run_time = $current_time - $start_time;


#########################################################
# Calculate total time by adding up stages (as a check) #
#########################################################
if ($make_VCF_file eq "true"){ $total_time = $make_raw_bam_time_total + $make_clean_bam_time_total + $make_VCF_time; }
if ($variant_caller eq "HaplotypeCallerGVCF"){ $total_time = $make_raw_bam_time_total + $make_clean_bam_time_total + $make_gVCF_time;}


#################################################
# Add times to the time log and the command log #
#################################################
print COMMAND_LOG "\n";
print COMMAND_TIMES_LOG "\nFinal times:\n\n";


&log_time($make_raw_bam_time_total,"Total time for making all raw BAM files");

&log_time($make_clean_bam_time_total,"Total time for making all clean BAM files");


if ($make_VCF_file eq "true")
{
	&log_time($make_VCF_time_total,"Total time for making all VCF files");
}

if ($variant_caller eq "HaplotypeCallerGVCF")
{
	&log_time($make_gVCF_time_total,"Total time for making all gVCF files");
}

&log_time($total_time,"Total time for all stages - added up");
&log_time($total_run_time,"Total time for all stages");


#############################
# Turn off logging          #
#############################

close(STDERR);

print COMMAND_LOG "\n\n";
print COMMAND_LOG "###################################################################\n";
print COMMAND_LOG "#                        End of fastq2vcf run                     #\n";
print COMMAND_LOG "###################################################################\n";

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

print	READMEOVERALL	"\t$final_raw_bam	     \tRaw alignments to the reference.\n\n";
print	READMEOVERALL	"\t(These are deleted but flagstat and validate files still exist for checking)\n\n\n";

print	READMEOVERALL	"FINAL BAM ALIGNMENT FILES AFTER PROCESSING BY GATK\n\n";

print	READMEOVERALL	"\t$final_bam		     \tBest alignments to the reference after processing by GATK\n\n\n";

print	READMEOVERALL	"\tLook at these files to check the quality of your final BAM files:\n\n";

print   READMEOVERALL   "\t\t$validate_final_out\n\n";
print   READMEOVERALL   "\t\t$flagstat_final_out\n\n\n";

print	READMEOVERALL	"PDF FILES\n\n";

print	READMEOVERALL	"\t$insert_size_pdf      \tInsert size histogram\n\n";
print	READMEOVERALL	"\t$insert_size_out      \tInsert size metrics file\n\n";

print	READMEOVERALL	"\t$WGS_metrics_file     \tWGS metrics file\n\n\n";


print	READMEOVERALL	"LOG FILES\n\n";

print	READMEOVERALL	"\t$command_log	         \tCommand log for the NGS pipeline (very useful for tracking errors)\n\n";
print	READMEOVERALL	"\t$screen_log_file      \tScreen log for the NGS pipeline\n\n\n";


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

&run_unix_command_no_error_reporting("cp $command_log /home/genetics/command_logs/$command_log","Copy command log to genetics folder");
&run_unix_command_no_error_reporting("cp $WGS_metrics_file /home/genetics/command_logs/$WGS_metrics_file","Copy WGS metrics file to genetics folder");

&move_to_results_all_folder("$WGS_metrics_file");
&move_to_results_all_folder("$screen_log_file");
&move_to_results_all_folder("$command_log");
&move_to_results_all_folder("$command_times_log");



#############################
# Send final e-mail to user #
#############################
open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $email_address\n";
print MAIL "From: $e_mail_from\n";
print MAIL "Subject: FASTQ2VCF: Run $run_title has finished\n\n";
## Mail Body
print MAIL "Run title:  $run_title\n\n";
print MAIL "Script:     fastq2vcf (GIT) started $version\n\n";

print MAIL "Your next generation sequence pipeline is complete\n\n";
print MAIL "For a list of commands see $command_log\n\n";
print MAIL "For more pipeline details see $screen_log_file and README files\n\n";
print MAIL "For further details on intermediate 'raw' BAM files (now deleted):\n\n";
print MAIL "  $validate_raw_out\n";
print MAIL "  $flagstat_raw_out\n\n";

print MAIL "For further details on final BAM files:\n\n";
print MAIL "  $validate_final_out\n";
print MAIL "  $flagstat_final_out\n\n";

if ($make_VCF_file eq "true")
{
	print MAIL "Size of final multi-column VCF file: \t$final_VCF_size\n\n";
}
if ($make_VCF_file eq "true")
{
	print MAIL "Size of final multi-column VCF file: \t$final_VCF_size\n\n";
}

if ($effect_predictor eq "vep")
{
	print MAIL "Size of VEP file made from VCF file: \t$final_vep_size\n\n";
}
if ($effect_predictor eq "snpEff")
{
	print MAIL "Size of snpEff-processed multi-column VCF file: \t$final_snpEff_vcf_size\n\n";
}

printf MAIL"Total run time: %d days, %d hours, %d minutes and %d seconds\n",(gmtime $total_run_time)[7,2,1,0];
close(MAIL);


&print_message("FASTQ2VCF ANALYSIS COMPLETE! YOU HAVE BEEN SENT AN EMAIL WITH DETAILS","message");

print "  Run title: $run_title\n\n";
print "  Run time: ";
printf "%d days, %d hours, %d minutes and %d seconds\n",(gmtime $total_run_time)[7,2,1,0];

print "\n\n  Results folder:                 \t$results_all_folder";

print "\n\n  For details of the analysis run:\tcheck $command_log and $screen_log_file for full details\n\n";

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
		print "\n  $file_to_be_moved could not be found to be moved to $results_all_folder\n";
		print COMMAND_LOG "\n  $file_to_be_moved could not be found to be moved to $results_all_folder\n";
	}
	
	if (-e "$file_to_be_moved")
	{
		$command = "mv  $file_to_be_moved $results_all_folder/$file_to_be_moved";
		print "\n  $file_to_be_moved was moved to $results_all_folder\n";
		system("$command");
		print COMMAND_LOG "\n  Move to results folder:  $file_to_be_moved was moved to $results_all_folder\n";
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
		print COMMAND_LOG "\n  $file_to_be_deleted could not be found to be deleted\n";
	}
	
	if (-e "$file_to_be_deleted")
	{
		$command = "rm  $file_to_be_deleted";
		system("$command");
		print COMMAND_LOG "\n  Deleted File: $file_to_be_deleted was deleted\n";
	}
}



#############################################
# Subroutine to execute unix command        #
#############################################

sub run_unix_command($;$)
{
	my $_unix_command 	= "";
	my $_step 			= "";	

	$_unix_command 		= $_[0];
	$_step 				= $_[1];

	print("\n$_unix_command\n");

	###########################
	# Run the command in unix #
	###########################
	my $returnCode = system("$_unix_command");
	
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	
	if (($loop_count > 0 ) && ($loop_count <= $no_of_files)) {print COMMAND_LOG "File: $loop_count/$no_of_files \t$_step\n";}
	else {print COMMAND_LOG "Step: $_step\n";}
	
	print COMMAND_LOG "\n$_unix_command\n";

	###############################################
	# Warn user if return code is not zero or one #
	###############################################
	if (($returnCode > 1) && (index ($_unix_command,"ValidateSamFile") == -1))
	{
		print COMMAND_LOG "\nWARNING: run_unix_command detected an ERROR ($returnCode) when running command: $_unix_command\n";

		# Only send unix error e-mail once per loop.
		if ($unix_error_reported eq "no") {&send_email_err_alert($returnCode,$_unix_command);}
		$unix_error_reported = "yes";
	}
}

sub run_unix_command_no_error_reporting($;$)
{
	my $_unix_command 	= "";
	my $_step 			= "";	

	$_unix_command 		= $_[0];
	$_step 				= $_[1];

	print("\n$_unix_command\n");

	###########################
	# Run the command in unix #
	###########################
	my $returnCode = system("$_unix_command");
	
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	
	if (($loop_count > 0 ) && ($loop_count <= $no_of_files)) {print COMMAND_LOG "File: $loop_count/$no_of_files \t$_step\n";}
	else {print COMMAND_LOG "Step: $_step\n";}
	
	print COMMAND_LOG "\n$_unix_command\n";
}

##############################################
# Send an email alert on error               #
##############################################
sub send_email_err_alert($;$) 
{
	my $_returnCode = "";
	my $_unix_command = "";
		
	$_returnCode = $_[0];
	$_unix_command = $_[1];

	open(MAIL, "|/usr/sbin/sendmail -t");

	## Mail Header
	print MAIL "To: $email_address\n";
	print MAIL "From: $e_mail_from\n";
	print MAIL "Subject: FASTQ2VCF: Error detected in run $run_title\n\n";

	## Mail Body
	print MAIL "Script:    \tfastq2vcf (GIT) started $version\n";
	print MAIL "Run title: \t$run_title\n\n";
	
	print MAIL "A unix error was detected during script processing (error number $_returnCode)\n\n";
	print MAIL "The script will continue to run, but you may wish to review the output so far to decide whether you wish to manually abort processing or allow the script to continue.\n\n";
	print MAIL "(Note: if the command is a 'mv' command, sometimes it seems to give an error when it shouldn't be there.).\n\n";

	print MAIL "Command that generated the error:\n\n";
	print MAIL "$_unix_command\n\n";

	close(MAIL);
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
		print COMMAND_LOG "\n  Input file:  $inputfile\t\tSize: $filesize\n";
		print "\n  Input file: $inputfile\t\tSize: $filesize\n\n";
	}
	else
	{
		$filesize=0;
		print COMMAND_LOG "\n  Input file:  $inputfile\t\tSize: Not found\n";
		print "\n  Input file: $inputfile\t\tSize: Not found\n\n";
	}
} # record_input_file_size


##############################################
# Subroutine to record size of output file   #
##############################################

sub record_output_file_size
{
	my $_outputfile = "";
	my $filesize = "";	

	$_outputfile = $_[0];
	
	if (-e $_outputfile)
	{
		$filesize = -s "$_outputfile";
		print COMMAND_LOG "\n  Output file: $_outputfile\t\tSize: $filesize\n\n";
		print "\n  Output file: $_outputfile\t\tSize: $filesize\n\n";
	}
	else
	{
		$filesize = 0;
		print COMMAND_LOG "\n  Output file: $_outputfile\t\tSize: Not found\n\n";
		print "\n  Output file: $_outputfile\t\tSize: Not found\n\n";
	}

	$current_time = time();
	$total_run_time = $current_time - $start_time;
	$stage_time	= $current_time - $last_current_time;
	printf COMMAND_LOG "  Time for stage:\t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $stage_time)[7,2,1,0];
	printf COMMAND_LOG "  Total time:    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $total_run_time)[7,2,1,0];
	
	
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
				print MAIL "From: $e_mail_from\n";
				print MAIL "Subject: NGS ANALYSIS fastq2vcf ZERO FILE SIZE: Run $run_title.  Sample $sample_name. File $_outputfile.\n\n";
				## Mail Body
				print MAIL "fastq2vcf script (GIT) started $version\n\n";
				print MAIL "Run:    \t$run_title\n";
				print MAIL "Sample: \t$sample_name\n";
				print MAIL "File:   \t$_outputfile\n\n";
				print MAIL "The output file $_outputfile has zero file size or can't be found.\n\n";
				print MAIL "Something may be wrong.  Please check\n\n";
			close(MAIL);
		}
		
	} # if ($zero_size_reported eq "no";)
	
} # record_output_file_size
	
	

##############################################
# Subroutine to record size of output file   #
# (without sending an e-mail if it is zero)  #
##############################################

sub record_output_file_size_no_mail
{
	my $_outputfile = "";
	my $filesize = "";	

	$_outputfile = $_[0];
	
	if (-e $_outputfile)
	{
		$filesize = -s "$_outputfile";
		print COMMAND_LOG "\n  Output file: $_outputfile\t\tSize: $filesize\n";
		print "\n  Output file: $_outputfile\t\tSize: $filesize\n\n";
	}
	else
	{
		$filesize=0;
		print COMMAND_LOG "\n  Output file: $_outputfile\t\tSize: Not found\n";
		print "\n  Output file: $_outputfile\t\tSize: Not found\n\n";
	}
	
	$current_time = time();
	$total_run_time = $current_time - $start_time;
	$stage_time	= $current_time - $last_current_time;
	printf COMMAND_LOG "  Time for stage:\t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $stage_time)[7,2,1,0];
	printf COMMAND_LOG "  Total time:    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $total_run_time)[7,2,1,0];
	
	$last_current_time = $current_time;

	if ($testing_mode eq "on")
	{
		print "\nPRESS 'RETURN' TO CONTINUE\n";
		$answer = <STDIN>;
	}

} # record_output_file_size_no_mail

	

#############################################
# Subroutine to rename files   			    #
#############################################

sub rename_file
{
	my $_old_name = "";
	my $_new_name = "";	

	$_old_name = $_[0];
	$_new_name = $_[1];
	
	if (-e "$_old_name")
	{
		$command = "mv $_old_name $_new_name";
		system("$command");
		print COMMAND_LOG "  Renamed File:  $_old_name ==> $_new_name\n";
	}
	else
	{
		print COMMAND_LOG "  Failed to find $_old_name to rename as $_new_name\n";
	}
} # rename_file



####################################################################
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test.seq" from "test.seq.txt")                     #
####################################################################

sub get_prefix_right
{
	my $_filename = "";
	$_filename = $_[0];
	
	if (rindex($_filename,".") > 0)
	{
		$_filename = substr($_filename, 0, rindex($_filename,"."));
	}
	else
	{
		$_filename = $_filename;
	}
} # get_prefix_right


####################################################################
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test" from "test.seq.txt")                         #
####################################################################
sub get_prefix
{
	my $_filename = "";
	$_filename = $_[0];
	
	if (index($_filename,".") > 0)
	{
		$_filename = substr($_filename, 0, index($_filename,"."));
	}
	else
	{
		$_filename = $_filename;
	}
} # get_prefix




	
########################################################################
# If testing is set to ON then this pauses the script after each stage #
########################################################################
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

	printf "   $description: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $time_to_log)[7,2,1,0];
	printf COMMAND_LOG "   $description: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $time_to_log)[7,2,1,0];
	printf COMMAND_TIMES_LOG "   $description: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $time_to_log)[7,2,1,0];
}


####################################################
# Subroutine to print to screen and to COMMAND_LOG #
####################################################
sub print_both
{
	my $_message = $_[0];

	print "$_message";
	print COMMAND_LOG "$_message";
}

##########################################################
# Subroutine to pause until user hits 'return'           #
##########################################################
sub pause
{
	my $_answer = "";
	print "\n Press RETURN to continue\n";
	$_answer=<STDIN>;
}

######################################################################################
# Delete previous file in the pipeline if output file has been created successfully  #
######################################################################################
sub delete_first_file_if_second_ok
{
	my $_input_file = $_[0];
	my $_output_file = $_[1];

	my $_output_file_size 	= 0;
	my $_bai_file			= "";
	my $_prefix				= "";

	if (-e $_output_file)
	{
		$_output_file_size = -s "$_input_file";
	}
	else
	{
		$_output_file_size = 0;
	}

	if (-e $_output_file)
	{
		if ($_output_file_size > 2000)
		{
			&print_both("\n  Output file $_output_file was greater than 2K \n  so previous file $_input_file can be deleted\n\n");
			&delete_file("$_input_file");

			# If it is a BAM file then delete the BAI file as well.
			if (substr($_input_file, -4) eq ".bam")
			{
				$_prefix = &get_prefix_right($_input_file);
				$_bai_file = $_prefix.".bai";
				&delete_file("$_bai_file");
			}
		} 
		else 
		{
			&print_both("  Output file $_output_file was not greater than 2K so previous file $_input_file was not deleted\n\n");
		}
	} # if outfile exists
	else
	{
		&print_both("  Output file $_output_file was not found so previous file $_input_file was not deleted\n\n");
	}
}# delete_first_file_if_second_ok



#############################################
# Subroutine to delete files                #
#############################################

sub delete_file_1
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
