#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	run_variant_caller 						                            #     
#									                                    #
#	THIS PERL SCRIPT WILL RUN various variant callers on all BAM files  #
#									                                    #
#########################################################################

#############################
# Mike Boursnell July 2012  #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;
use Cwd;

my $version							= "w23";
my $testing_mode					= "off";

###############################################################################
#            GATK memory suggestions                                          #
#-----------------------------------------------------------------------------#
#Tool				RTC		IR		BR		PR		RR		   UG             #
#Available modes 	NT 		SG 	  NCT,SG  	NCT	  	SG      NT,NCT,SG         #
#Cluster nodes 		1 		4 		4 		1 		4 	    4 / 4 / 4         #
#CPU threads (-nct) 1 		1 		8 		4-8 	1       3 / 6 / 24        #
#Data threads (-nt) 24 		1 		1 		1 		1       8 / 4 / 1         #
#Memory (Gb) 		48 		4 		4 		4 		4      32 / 16 / 4        #
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


my $memory_in_Gb				= "4"; # memory setting for Java in gigabytes
my $no_of_threads_bwa			= "2";

my $no_of_threads_gatk_ug_nct	= "1";  # UnifiedGenotyper -nct (or HaplotypeCaller) (UG has to be 1)
my $no_of_threads_gatk_ug_nt	= "2";  # UnifiedGenotyper -nt (or HaplotypeCaller)

my $no_of_threads_gatk_pr_nct	= "2";  # PrintReads -nct                <== This one is used for PrintReads
my $no_of_threads_gatk_br_nct	= "2";  # BaseRecalibrator -nct          <== This one is used for BaseRecalibrator
my $no_of_threads_gatk_rtc_nt	= "2";
my $workstation 				= "unknown";
my $e_mail_from 				= ""; # Who e-mails come from

# Workstation version - (if 'w' in version name)
if (index($version,"w") > -1)
{
	$memory_in_Gb				= "45"; # memory setting for Java in gigabytes
	$no_of_threads_bwa			= "8";
	$no_of_threads_gatk_ug_nct	= "1";  # UnifiedGenotyper -nct (or HaplotypeCaller) (UG has to be 1)
	$no_of_threads_gatk_ug_nt	= "16";  # UnifiedGenotyper -nt (but NOT HaplotypeCaller)
	$no_of_threads_gatk_pr_nct	= "8";  # PrintReads -nct                <== This one is used for PrintReads
	$no_of_threads_gatk_br_nct	= "8";  # BaseRecalibrator -nct          <== This one is used for BaseRecalibrator
	$no_of_threads_gatk_rtc_nt	= "16"; # RealignerTargetCreator         <== This one is used for RealignerTargetCreator
	$workstation 				= "true";
	$e_mail_from 				= 'NGS_analysis@gen-x1404-ws01.aht.org.uk'; # Who e-mails come from
}
else
# Samba64 version (if no 'w' in version name)
{
	$memory_in_Gb				= "4"; # memory setting for Java in gigabytes
	$no_of_threads_bwa			= "2";
	$no_of_threads_gatk_ug_nct	= "1";  # UnifiedGenotyper -nct (or HaplotypeCaller) (UG has to be 1)
	$no_of_threads_gatk_ug_nt	= "2";  # UnifiedGenotyper -nt (but NOT HaplotypeCaller)
	$no_of_threads_gatk_pr_nct	= "2";  # PrintReads -nct                <== This one is used for PrintReads
	$no_of_threads_gatk_br_nct	= "2";  # BaseRecalibrator -nct          <== This one is used for BaseRecalibrator
	$no_of_threads_gatk_rtc_nt	= "2"; # RealignerTargetCreator         <== This one is used for RealignerTargetCreator
	$workstation 				= "false";
	$e_mail_from 				= 'NGS_analysis@samba64.org.uk'; # Who e-mails come from
}



#Various Parameters (e.g. for the Unified Genotyper)
my $gatk_directory					= "gatk";
my $default_mem						= "-Xmx60g";
my $GATK_validation_stringency		= "LENIENT";
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
my $has_threading_option			= ""; # true or false. Can it use multiple CPU threads
my $has_memory_option				= ""; # true or false. Can it use memory setting like -Xmx4g
my $choice_ok						= "false";
my $fix_misencoded_qual_scores		= ""; # 'yes' or 'no'

# File names
my $bcf_file						= ""; # for samtools mpileup
my $list_file						= ""; # file of file names
my $bam_file						= "";
my $bai_file						= "";
my $bam_file_in_results				= "";
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
my $bam_file_for_unified_genotyper 	= "";
my $species							= "";
my $ref_seq_name					= "";
my $other_ref_sequence 				= "false";
my $current_directory 				= "";
my $folder 							= "";
my $bam_file_in_folder				= "";
my $variant_caller					= ""; # UnifiedGenotyper, HaplotypeCaller, mpileup or platypus
my $search_string					= "";

#Timers
my $start_time						= "";
my $end_time						= "";
my $current_time					= "";

my @bam_file_array					= ();
my @chromosome_end_array			= (); # for running all chromosomes separately

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "              run_variant_caller   \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';

if ( index ($version,"w") == -1 ){ print "Version $version\n\n"; } else { print "Version $version (for Workstation)\n\n"; }

print "    - This program runs a Variant Caller to convert multiple BAM files to a VCF file\n\n";
print "    - Note: this is also done in fastq2vcf but this program can be used for just the last stage, i.e. to \n";
print "      process the final BAM files from fastq2vcf into VCF files\n\n";
print "      There are several reasons why you might want to do this:\n\n";

print "      - You might want to try a different variant caller\n";
print "      - You might want to try different parameters for the variant caller\n";
print "      - The BAM files may have been made in two runs of fastq2vcf, and you want to combine them all for the VCF file\n\n";

print "    - It allows a choice of Variant Callers: \n\n";

print "              - GATK HaplotypeCaller [GATK preferred method]\n";
print "              - GATK UnifiedGenotyper\n";
print "              - GATK HaplotypeCaller (new GVCF method)\n";
print "              - samtools mpileup\n";
print "              - platypus\n";
print "              - freeBayes\n\n";

print color 'reset';


################################
# Make the temp java directory #
# Check if staticdata exists   #
################################
$staticdata = "$ENV{HOME}/staticdata";
if (! stat($staticdata))
{
            $tempdir = "$ENV{HOME}/javatempdir"; # Individual's HOME space (if staticdata isn't online)
            $temp_dir_string = " -Djava.io.tmpdir=$tempdir";
            if (! -e $tempdir)
            {
                        unless(mkdir $tempdir){die "Unable to create temporary Java directory $tempdir";}
                        $temp_dir_string = " -Djava.io.tmpdir=$tempdir";            
            }
} else {
            $tempdir = "$ENV{HOME}/staticdata/javatempdir"; # Moved from individual's HOME space into staticdata 02/05/14
            $temp_dir_string = " -Djava.io.tmpdir=$tempdir";
            if (! -e $tempdir)
            {
                        unless(mkdir $tempdir){die "Unable to create temporary Java directory $tempdir";}
                        $temp_dir_string = " -Djava.io.tmpdir=$tempdir";            
            }
} 

&print_message("The input is a file with a list of the original BAM file names.","input");


##########################################################
# Get the file that contains a list of the BAM filenames #
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
# Open the list file to get the list of BAM file names #
########################################################
open (LIST, "$list_file") || die "Cannot open $list_file";
$list_count=1;

while ($bam_file = <LIST> ) 
{
	chomp $bam_file;
	
	$prefix = &get_prefix("$bam_file");
	
	print "$list_count\t$bam_file\t";
	
	if (-e $bam_file){print "File found\n"}

	if (! -e $bam_file)
	{
		print "File not found\n";
		#print "\n\n";
		#print "###########################################################################################\n";
		#print "$bam_file cannot be found\n\n";
		#print "Correct the input file of file names, or the locations of the input BAM files and try again\n";
		#print "###########################################################################################\n\n";
		#close LIST;
		#exit;
	}
	
	$bam_file_array[$list_count]=$bam_file;
	$list_count=$list_count + 1;
}

close LIST;


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

if (substr($answer,0,1) eq "1" ){$fix_misencoded_qual_scores = "no"; $GATK_fix_quals_string = ""; $quality_scores = "New"}
if (substr($answer,0,1) eq "2" ){$fix_misencoded_qual_scores = "yes"; $GATK_fix_quals_string = "-fixMisencodedQuals "; $quality_scores = "Old"}




###############
# Name of run #
###############
$run_title = "ls";
until ($run_title ne "ls")
{
	&print_message("Type a name for this run (for output file prefix):","input");

	$run_title = <STDIN>; chomp $run_title;

	if ($run_title eq "ls")
	{
		print "\n";
		system ("ls *_run_variant_caller_command_log.out");
		print "\n";
	}
}


############################
# Open the screen log file #
############################
$screen_log_file = "$run_title"."_run_variant_caller_screen_log.out";
$| = 1;
open(STDERR, "| tee $screen_log_file");


##################################
# Define data files              #
##################################

&print_message("Which reference sequence do you want to use?","input");

print "   <1>  CanFam3\n";
print "   <2>  EquCab2\n";
print "   <3>  Human\n\n";

print "   <4>  Strep. equi\n";
print "   <5>  Strep. zoo\n\n";

print "   <9>  Other\n\n";


$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3";}

if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/equcab2/equcab2.fasta"; $ref_seq_name = "equcab2";}
if (substr($answer,0,1) eq "3" ){$ref = "/home/genetics/human/human.fasta"; $ref_seq_name = "human";}

if (substr($answer,0,1) eq "4" ){$ref = "/home/genetics/strep_equi/strep_equi.fasta"; $ref_seq_name = "s_equi";}
if (substr($answer,0,1) eq "5" ){$ref = "/home/genetics/strep_zoo/strep_zoo.fasta"; $ref_seq_name = "s_zoo";}

if ($ref eq ""){print "\n\nYou have to choose a reference sequence\n\n";exit;}


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


###########################################################################
# Which variant caller do you want?                                       #
###########################################################################

&print_message("Which variant caller would you like to use?","input");

print "   <1> GATK Haplotype Caller              [this is GATK's preferred caller]\n\n";
print "   <2> GATK Unified Genotyper             [this is what we used to use]\n\n";
print "   <3> GATK Haplotype Caller using gVCF   [this is GATK's most up to date method]\n\n";
print "   <4> samtools mpileup       \n\n";
print "   <5> platypus               \n\n";
print "   <6> freeBayes              \n\n";

$answer = <STDIN>;
chomp $answer;

# Default is now Haplotype Caller
if ($answer eq ""){$answer = "1"}

if (substr($answer,0,1) eq "1" ){$variant_caller = "HaplotypeCaller";$has_threading_option = "true";$has_memory_option = "true"}
if (substr($answer,0,1) eq "2" ){$variant_caller = "UnifiedGenotyper";$has_threading_option = "true";$has_memory_option = "true"}
if (substr($answer,0,1) eq "3" ){$variant_caller = "HaplotypeCallerGVCF";$has_threading_option = "true";$has_memory_option = "true"}
if (substr($answer,0,1) eq "4" ){$variant_caller = "mpileup";$has_threading_option = "false";$has_memory_option = "false"}
if (substr($answer,0,1) eq "5" ){$variant_caller = "platypus";$has_threading_option = "false";$has_memory_option = "false"}
if (substr($answer,0,1) eq "6" ){$variant_caller = "freeBayes";$has_threading_option = "false";$has_memory_option = "false"}


if ($variant_caller eq "")
{
	&print_message("No variant caller chosen","warning");
	exit;
}

$calls = "BOTH_TOGETHER"; # default

if ($variant_caller eq "UnifiedGenotyper") 
{
	##################################
	# Ask if you want SNPs or INDELS #
	##################################

	&print_message("Do you want to call SNPs, Indels or both?","input");

	print "  <1> Both together      [default]\n";
	print "  <2> Both separately\n";
	print "  <3> SNPs\n";
	print "  <4> Indels\n\n";
	
	$answer = <STDIN>;
	chomp $answer;

	if ($answer eq "1"){$calls = "BOTH_TOGETHER"}
	if ($answer eq "2"){$calls = "BOTH"}
	if ($answer eq "3"){$calls = "SNPS"}
	if ($answer eq "4"){$calls = "INDELS"}
	if ($answer eq ""){$calls = "BOTH_TOGETHER"}
}


###########################################################################
# Region preferences                                                      #
# Would you like to focus the analysis on a specific region of the genome #
###########################################################################

&print_message("Would you like to focus on a specific region of the genome?","input");

print "   <1>  YES\n";
print "   <2>  NO\n\n";

$answer = <STDIN>;
chomp $answer;

if ($answer eq ""){$answer = "1"} # default

if (substr($answer,0,1) eq "1"){$use_defined_region = "yes"}
if (substr($answer,0,1) eq "2"){$use_defined_region = "no"}

if ($use_defined_region eq "yes")
{
	print "\nPlease define your region of interest (eg 'chr5:21000000-23000000'):      ";
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
	$samtools_region_string = "-r $region";
	$GATK_chromosome_only_string = "-L $chromosome";
	$platypus_region_string = "--regions $region";
	$freeBayes_region_string = " --region $region";

} # If use defined region

$no_of_files=$list_count - 1;


###########################################################################
# Would you like to run all the chromosomes separately                    #
# (only for Unified Genotype currently)                                   #
###########################################################################
if (($variant_caller eq "UnifiedGenotyper") && (($ref_seq_name eq "canfam3")  || ($ref_seq_name eq "canfam3ens")))
{
	&print_message("Would you like to run each of the chromosomes separately? (i.e. make a separate VCF for each chromosome)","input");

	print "   <1>  YES\n";
	print "   <2>  NO  [default]\n\n";

	$answer = <STDIN>;
	chomp $answer;

	if ($answer eq ""){$answer = "2"} # default

	if (substr($answer,0,1) eq "1"){$run_chromosomes_separately = "yes"}
	if (substr($answer,0,1) eq "2"){$run_chromosomes_separately = "no"}
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


	print "   <1>  YES [default]\n";
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



##################################################################
# Number of CPU threads to use (but only if program uses them)   #
##################################################################
if ($has_threading_option eq "true")
{
	
	####################################################
	# Get user input on all memory and threads options #
	####################################################
	$memory_string = "-Xmx".$memory_in_Gb."g";

	 while ($choice_ok eq "false")
	{
			&print_message("Current memory and CPU thread settings","message");

			print "Workstation: $workstation\n\n";

			print "Memory setting for java steps (Gigabytes):                     \t$memory_in_Gb\n";
			print "Number of CPU threads for bwa:                                 \t$no_of_threads_bwa\n";
			print "Number of CPU threads (-nct) for GATK variant calling:         \t$no_of_threads_gatk_ug_nct\n";

			if ($variant_caller eq "UnifiedGenotyper")
			{
				print "Number of data threads (-nt) for GATK variant calling:         \t$no_of_threads_gatk_ug_nt\n";
			}

			print "Number of data threads (-nt) for GATK RealignerTargetCreator:  \t$no_of_threads_gatk_rtc_nt\n";
			print "Number of CPU threads (-nct) for GATK BaseRecalibrator:        \t$no_of_threads_gatk_br_nct\n";
			print "Number of CPU threads (-nct) for GATK PrintReads:              \t$no_of_threads_gatk_pr_nct\n";

			print "\nWould you like to change these? (y/n)   [default = 'n'] ";
			$answer=<STDIN>;chomp $answer;$answer = lc $answer;

			if ($answer eq "y")
			{
				&print_message("Enter new values (press 'return' to keep existing value)","input");

				$choice_ok = "false";

				print "  Memory setting for java steps (in Gigabytes)                  [current value = $memory_in_Gb]:      ";
				$input = <STDIN>;chomp $input;
				if ($input ne ""){$memory_string = $input}
				if ($memory_string > 60){$memory_string = "60"}

				print "  Number of CPU threads for bwa                                 [current value = $no_of_threads_bwa]:      ";
				$input = <STDIN>;chomp $input;
				if ($input ne ""){$no_of_threads_bwa = $input}
				if ($no_of_threads_bwa > 16){$no_of_threads_bwa = "16"}

				print "  Number of CPU threads (-nct) for GATK variant calling?        [current value = $no_of_threads_gatk_ug_nct]:      ";
				$input = <STDIN>;chomp $input;
				if ($input ne ""){$no_of_threads_gatk_ug_nct = $input}
				if ($no_of_threads_gatk_ug_nct > 16){$no_of_threads_gatk_ug_nct = "16"}

				if ($variant_caller eq "UnifiedGenotyper")
				{
					print "  Number of data threads (-nt) for GATK variant calling?        [current value = $no_of_threads_gatk_ug_nt]:      ";
					$input = <STDIN>;chomp $input;
					if ($input ne ""){$no_of_threads_gatk_ug_nt = $input}
					if ($no_of_threads_gatk_ug_nt > 16){$no_of_threads_gatk_ug_nt = "16"}
				}

				print "  Number of data threads (-nt) for GATK RealignerTargetCreator? [current value = $no_of_threads_gatk_rtc_nt]:      ";
				$input = <STDIN>;chomp $input;
				if ($input ne ""){$no_of_threads_gatk_rtc_nt = $input}
				if ($no_of_threads_gatk_rtc_nt > 16){$no_of_threads_gatk_rtc_nt = "16"}

				print "  Number of CPU threads (-nt) for GATK BaseRecalibrator?        [current value = $no_of_threads_gatk_br_nct]:      ";
				$input = <STDIN>;chomp $input;
				if ($input ne ""){$no_of_threads_gatk_br_nct = $input}
				if ($no_of_threads_gatk_br_nct > 16){$no_of_threads_gatk_br_nct = "16"}

				print "  Number of CPU threads (-nt) for GATK PrintReads?              [current value = $no_of_threads_gatk_pr_nct]:      ";
				$input = <STDIN>;chomp $input;
				if ($input ne ""){$no_of_threads_gatk_pr_nct = $input}
				if ($no_of_threads_gatk_pr_nct > 8){$no_of_threads_gatk_pr_nct = "8"}

			}# Yes to change memory options
			else
			{
				$choice_ok = "true";
			}

	} # while choice_ok eq "false"

	# Make up java memory_string
	$memory_string = "-Xmx".$memory_in_Gb."g";

} # has_threading_option



#########################
# open Command Log file #
#########################
$command_log_file = "$run_title"."_run_variant_caller_command_log.out";
open (COMMAND_LOG, ">$command_log_file")|| die "Cannot create output file: $command_log_file";
print COMMAND_LOG "COMMAND LOG for run_variant_caller version $version\n";


################################
# Show all the summary details #
################################
&print_message_both("SUMMARY details","message");

&print_both("Run title:                  \t$run_title\n");
&print_both("Variant caller:             \t$variant_caller\n");
&print_both("Reference sequence:         \t$ref\n\n");


if (($variant_caller eq "UnifiedGenotyper") || ($variant_caller eq "HaplotypeCaller") || ($variant_caller eq "HaplotypeCallerGVCF"))
{
	&print_both("GATK region string:         \t$GATK_region_string\n\n");

	&print_both("GATK constants:\n\n");

	&print_both("    --stand_emit_conf:      \t$stand_emit_conf\n");
	&print_both("    --stand_call_conf:      \t$stand_call_conf\n");
	&print_both("    --max_alt_alleles:      \t$max_alt_alleles\n\n");
} # GATK

if ($has_memory_option eq "true"){&print_both("Memory setting:             \t$memory_string\n");}

if ($has_threading_option eq "true")
{
	&print_both("No. of CPU threads:         \t$no_of_threads_gatk_ug_nct\n");
	&print_both("No. of data threads:        \t$no_of_threads_gatk_ug_nt\n\n");
}

if ($use_defined_region eq "yes")
{	
	&print_both("Chromosome:                 \t$chromosome\n");

	if (($variant_caller eq "UnifiedGenotyper") || ($variant_caller eq "HaplotypeCaller") || ($variant_caller eq "HaplotypeCallerGVCF"))
	{&print_both("GATK region string:         \t$GATK_region_string\n")}

	if ($variant_caller eq "mpileup"){&print_both("samtools region string:         \t$samtools_region_string\n")}
	if ($variant_caller eq "platypus"){&print_both("platypus region string:         \t$platypus_region_string\n")}
	if ($variant_caller eq "freeBayes"){&print_both("freeBayes region string:         \t$freeBayes_region_string\n")}
	if ($variant_caller eq "platypus"){&print_both("platypus region string:         \t$platypus_region_string\n")}
}

&print_both("\nThere are $no_of_files BAM files\n\n");

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	&print_both("File $list_count	\t$bam_file_array[$list_count]\n");
}

print "\n\n  >> Press return to continue\n\n";
$answer=<STDIN>;
chomp $answer;
if ($answer eq "q" || $answer eq "Q"){exit;}



#####################################
# First get a list of the BAM files #
#####################################
for ($list_count=1;$list_count <=$no_of_files;$list_count++)
{
	$bam_file = $bam_file_array[$list_count];
	$sample_name = &get_prefix ($bam_file);

	print "File: $list_count \t$bam_file\n";
		
	$GATK_input_string = $GATK_input_string." -I $bam_file";

} # list_count loop

$GATK_input_string = $GATK_input_string." ";


##########################################################
#  Make a VCF file using all the BAM files in parallel   #
##########################################################

&print_message_both("Making VCF using $variant_caller","message");
$start_time = time();
printf COMMAND_LOG "  Start time:    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $start_time)[7,2,1,0];


#####################
# Unified Genotyper #
#####################
if (($variant_caller eq "UnifiedGenotyper")  && ($run_chromosomes_separately eq "no"))
{
	if (($calls eq "SNPS") || ($calls eq "BOTH"))
	{
		&print_message ("Running UnifiedGenotyper for SNPs... ","message");

		$output_vcf = "$run_title"."_SNPs_UG.vcf";
			
		if ($use_default_stand_values eq "false")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_fix_quals_string $GATK_region_string -glm SNP $GATK_input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o $output_vcf -S $GATK_validation_stringency -nt $no_of_threads_gatk_ug_nt -nct $no_of_threads_gatk_ug_nct");	
		}
		if ($use_default_stand_values eq "true")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_fix_quals_string $GATK_region_string -glm SNP $GATK_input_string -o $output_vcf -S $GATK_validation_stringency -nt $no_of_threads_gatk_ug_nt -nct $no_of_threads_gatk_ug_nct");	
		}
	
		&print_message("A SNP VCF file using all the BAM files in parallel has been created","message");

		&record_output_file_size("$output_vcf");
	} # SNPS


	if (($calls eq "INDELS") || ($calls eq "BOTH"))
	{
		&print_message ("Running UnifiedGenotyper for INDELs... ","message");
		
		$output_vcf = "$run_title"."_Indels_UG.vcf";
		
		if ($use_default_stand_values eq "false")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_fix_quals_string $GATK_region_string -glm INDEL $GATK_input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency -nt $no_of_threads_gatk_ug_nt -nct $no_of_threads_gatk_ug_nct");
		}
		if ($use_default_stand_values eq "true")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_fix_quals_string $GATK_region_string -glm INDEL $GATK_input_string --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency -nt $no_of_threads_gatk_ug_nt -nct $no_of_threads_gatk_ug_nct");
		}
				
		&print_message("An Indel VCF file using all the BAM files in parallel has been created","message");

		&record_output_file_size("$output_vcf");
	} # INDELS


	if ($calls eq "BOTH_TOGETHER")
	{
		&print_message ("Running UnifiedGenotyper for SNPs and INDELS together... ","message");
		
		$output_vcf = "$run_title"."_variants_UG.vcf";
		
		if ($use_default_stand_values eq "false")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_fix_quals_string $GATK_region_string -glm BOTH $GATK_input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency-nt $no_of_threads_gatk_ug_nt -nct $no_of_threads_gatk_ug_nct");
		}
		if ($use_default_stand_values eq "true")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_fix_quals_string $GATK_region_string -glm BOTH $GATK_input_string --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency -nt $no_of_threads_gatk_ug_nt -nct $no_of_threads_gatk_ug_nct");
		}

		&print_message("A SNP and Indel VCF file using all the BAM files in parallel has been created","message");

		&record_output_file_size("$output_vcf");
	} # BOTH

} # variant_caller is Unified Genotyper


#############################################
# Unified Genotyper on separate chromosomes #
#############################################
if (($variant_caller eq "UnifiedGenotyper")  && ($run_chromosomes_separately eq "yes"))
{
	$calls = "BOTH_TOGETHER"; # fix this for calling SNPs and Indels together

	for ($chromosome_count = 1; $chromosome_count <= 39; $chromosome_count++)
	{
		#####################################################
		# Change region and output file for each chromosome #
		#####################################################
		$GATK_region_string = "-L chr".$chromosome_count;

		# Deal with chromosome X
		if (($ref_seq_name eq "canfam3") || ($ref_seq_name eq "canfam3ens"))
		{
			if ($chromosome_count == 39){$GATK_region_string = "-L chrX"}
		}

		$output_vcf = "$run_title"."_variants_UG_chr".$chromosome_count.".vcf";

		if ($use_default_stand_values eq "false")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_fix_quals_string $GATK_region_string -glm BOTH $GATK_input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency -nt $no_of_threads_gatk_ug_nt -nct $no_of_threads_gatk_ug_nct");
		}
		if ($use_default_stand_values eq "true")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_fix_quals_string $GATK_region_string -glm BOTH $GATK_input_string --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency -nt $no_of_threads_gatk_ug_nt -nct $no_of_threads_gatk_ug_nct");
		}

		&record_output_file_size("$output_vcf");

	} # chromosome count

	&print_message("A SNP and Indel VCF file using all the BAM files in parallel has been created [separate chromosomes]","message");
} # UnifiedGenotyper on chromosomes separately


#####################
# Haplotype Caller  #
#####################
if ($variant_caller eq "HaplotypeCaller")
{
	
	&print_message("Running HaplotypeCaller for SNPs and Indels... ","message");
	
	$output_vcf = "$run_title"."_variants_HaplotypeCaller.vcf";
	
	if ($use_default_stand_values eq "false")
	{
		&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller $GATK_fix_quals_string $GATK_region_string  $GATK_input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads_gatk_ug_nct");	
	}
	if ($use_default_stand_values eq "true")
	{
		&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller $GATK_fix_quals_string $GATK_region_string  $GATK_input_string -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads_gatk_ug_nct");	
	}

	&print_message("A SNP and Indel VCF file using GATK HaplotypeCaller has been created","message");

	&record_output_file_size("$output_vcf");

} # variant_caller is Haplotype Caller



#########################
# Haplotype Caller gVCF #
#########################
if ($variant_caller eq "HaplotypeCallerGVCF")
{

	&print_message("Running HaplotypeCaller with the gVCF option for SNPs and Indels... ","message");
	

	#############################################################################
	# Loop through the BAM files - creating a gVCF file for each one separately #
	#############################################################################
	
	for ($list_count=1;$list_count <=$no_of_files;$list_count++)
	{
		$bam_file = $bam_file_array[$list_count];
		$sample_name = &get_prefix ($bam_file);

		$output_gvcf = "$run_title"."_"."$sample_name".".gVCF";

		$GVCF_input_string = $GVCF_input_string." -V $output_gvcf";

		&print_message("HaplotypeCaller gVCF method, stage 1. File $list_count.  Making gVCF file from $bam_file","message");

		if ($use_default_stand_values eq "false")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 $GATK_fix_quals_string $GATK_region_string  -I $bam_file -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o $output_gvcf -S $GATK_validation_stringency -nct $no_of_threads_gatk_ug_nct");	
		}
		if ($use_default_stand_values eq "true")
		{
			&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 $GATK_fix_quals_string $GATK_region_string  -I $bam_file -o $output_gvcf -S $GATK_validation_stringency -nct $no_of_threads_gatk_ug_nct");	
		}
		
		&record_output_file_size("$output_gvcf");

	} # list_count loop

	$GVCF_input_string = $GVCF_input_string." ";

	&print_message("HaplotypeCaller GVCF method, stage 1: A separate gVCF file has been created from each of the BAM files","message");


	########################################################
	# GenotypeGVCFs - run on all the files at once         #
	########################################################

	&print_message("Running GenotypeGVCFs on all $no_of_files files... ","message");

	&print_both("Input string: $GVCF_input_string\n\n");

	$output_vcf = "$run_title"."_variants_gVCF.vcf";

	&run_unix_command("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T GenotypeGVCFs -nt $no_of_threads_gatk_ug_nt $GATK_region_string  $GVCF_input_string -o $output_vcf -S $GATK_validation_stringency","Running GenotypeGVCFs on all VCF files");	

	&record_output_file_size("$output_vcf");

	&print_message("The VCF files from Haplotype_Caller gVCF mode have been merged into a single VCF file $output_vcf","message");

} # variant_caller is Haplotype Caller GVCF



#####################
# samtools mpileup  #
#####################
if ($variant_caller eq "mpileup")
{
	$bcf_file = "SNPS_and_Indels_MP_"."$run_title.bcf";

	$output_vcf = "$run_title"."_variants_mpileup.vcf";

	&print_message("Stage 1: Running samtools mpileup to create BCF file $bcf_file","message");
	&run_unix_command("samtools mpileup -b $list_file -uf $ref $samtools_region_string |bcftools view -bvcg - > $bcf_file");

	&print_message("Stage 2: Running bcftools view to create VCF file $output_vcf","message");	
	&run_unix_command("bcftools view $bcf_file > $output_vcf");

	&record_output_file_size("$output_vcf");

	&print_message("A SNP and Indel VCF file using samtools mpileup has been created","message");

} # variant_caller is samtools mpileup
#samtools mpileup -uf ref.fa aln1.bam aln2.bam | bcftools view -bvcg - > var.raw.bcf  
#bcftools view var.raw.bcf | vcfutils.pl varFilter -D100 > var.flt.vcf 


#####################
# platypus          #
#####################
if ($variant_caller eq "platypus")
{
	$output_vcf = "$run_title"."_variants_platypus.vcf";

	&print_message("Running platypus to create VCF file $output_vcf","message");
		
	&run_unix_command("python2.6 /opt/platypus/Platypus.py callVariants --refFile $ref --bamFiles $list_file --output $output_vcf --bufferSize=10000 --minReads=$min_reads_platypus  $platypus_region_string");	

	
	&record_output_file_size("$output_vcf");

	&print_message("A SNP and Indel VCF file using platypus has been created","message");

} # variant_caller is platypus
# python /opt/platypus/Platypus.py callVariants --refFile /home/genetics/canfam3/canfam3.fasta --bamFiles SHY8_input_bam.txt
	

#####################
# freeBayes         #
#####################
if ($variant_caller eq "freeBayes")
{
	print "\n\n\nNeed to rename BAI files\n\n";
	#exit;

	$output_vcf = "$run_title"."_variants_freeBayes.vcf";

	&print_message("Running freeBayes to create VCF file $output_vcf","message");
		
	&run_unix_command("/opt/freebayes/bin/freebayes -f $ref $freeBayes_region_string --bam-list $list_file --vcf $output_vcf")

	&record_output_file_size("$output_vcf");
	
	&print_message("A SNP and Indel VCF file using freeBayes has been created: $output_vcf","message");

} # variant_caller is freeBayes



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
&print_message("FINISHED RUNNING run_variant_caller","message");

&print_both("Variant caller used: \t$variant_caller\n\n");

if (($variant_caller  eq "UnifiedGenotyper") && ($calls eq "BOTH"))
{
		if (-e $output_SNPs_vcf){$size = -s "$output_SNPs_vcf"}else{$size=0}
		print "Output file SNPs:    \t$output_SNPs_vcf\tFile size: $size\n";

		if (-e $output_Indels_vcf){$size = -s "$output_Indels_vcf"}else{$size=0}
		print "Output file Indels:  \t$output_Indels_vcf\tFile size: $size\n\n";
}
else
{
	if (-e $output_vcf){$size = -s "$output_vcf"}else{$size=0}
	print "Output VCF file: \t$output_vcf\tFile size: $size\n\n";
}

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
