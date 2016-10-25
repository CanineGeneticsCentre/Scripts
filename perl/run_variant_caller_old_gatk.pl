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

my $version							= "12";


my $testing_mode					= "off";
my $include_reduce_reads			= "false";

#Various Parameters (e.g. for the Unified Genotyper)
my $gatk_directory					= "gatk_v2";
my $ref								= "";
my $make_parallel_VCF				= "yes";
my $mem								= "-Xmx8g";
my $no_of_threads					= "1"; # No of server CPUs to use (maximum is 8)
my $GATK_validation_stringency		= "SILENT";
my $temp_dir_string					= " -Djava.io.tmpdir=javatempdir";
my $tempdir							= "";
my $stand_emit_conf					= 30;
my $stand_call_conf					= 30;
my $min_reads_platypus				= 2; # Minimum number of reads required for a call in platypus (2 is the default)
my $max_alt_alleles					= 6; # This is how many alleles the INDEL caller in UnifiedGenotyper can allow

my $list_count						= 0;
my $no_of_files						= 0;
my $dot_pos							= 0;
my $size							= 0;
my $chromosome_count				= 0;

# Boolean
my $use_default_stand_values		= ""; # default values of stand_emit_conf and stand_call_conf
my $use_defined_region				= "";
my $reduce_reads					= ""; # decide whether to use GATK ReduceReads
my $all_reduced_files_found			= "true"; # check if all BAM files reduced with ReduceReads are found
my $use_existing_reduced			= "no"; # if the reduced BAM files are already there, use them
my $run_chromosomes_separately		= "";

# File names
my $bcf_file						= ""; # for samtools mpileup
my $list_file						= ""; # file of file names
my $bam_file						= "";
my $bai_file						= "";
my $bam_file_in_results				= "";
my $bai_file_in_results				= "";
my $output_vcf						= ""; # Output VCF file
my $output_SNPs_vcf					= ""; # Output VCF in separate SNPs
my $output_Indels_vcf				= ""; # Output VCF if separate Indels
my $command_log_file				= "";
my $screen_log_file					= ""; # log file for screen output
my $bam_file_reduced				= ""; # After ReduceReads



#Other variables
my $region							= "";
my $command							= "";
my $vcf_merge_string				= "";
my $prefix							= ""; # Part of sfile name before the dot
my $input_string					= "";
my $input_string_reduced			= "";
my $run_title						= "";
my $chromosome						= "";
my $GATK_chromosome_only_string				= ""; #This is the chr only string for the UnifiedGenotyper command line -L chr12
my $GATK_region_string				= ""; # -L chr12:50000000-60000000
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
print "Version $version\n\n";
print "    - This program runs a Variant Caller to convert multiple BAM files to a VCF file\n\n";
print "    - Note: this is also done in fastq2vcf but this program can be used for just the last stage, i.e. to \n";
print "      process the final BAM files from fastq2vcf into VCF files\n\n";
print "      There are several reasons why you might want to do this:\n\n";

print "      - You might want to try a different variant caller\n";
print "      - You might want to try different parameters for the variant caller\n";
print "      - The BAM files may have been made in two runs of fastq2vcf, and you want to combine them all for the VCF file\n\n";

print "    - It allows a choice of Variant Callers: \n\n";
print "              - GATK UnifiedGenotyper\n";
print "              - GATK HaplotypeCaller\n";
print "              - samtools mpileup\n";
print "              - platypus\n\n";

print color 'reset';





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

&print_message("The input is a file with a list of the original BAM file names.","input");

until (-e $list_file)
{
	print "Name of the file of file names:      ";
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
print "\n";


####################################################
# Open the list file to get the list of file names #
####################################################
open (LIST, "$list_file") || die "Cannot open $list_file";
$list_count=1;

while ($bam_file = <LIST> ) 
{
	chomp $bam_file;
	
	$prefix = &get_prefix("$bam_file");
	$bam_file_reduced = "$prefix"."_reduced.bam";
	
	if (! -e $bam_file_reduced){$all_reduced_files_found = "false"}
	
	print "$list_count\t$bam_file\t";
	
	if (-e $bam_file_reduced){print "\t$list_count\t$bam_file_reduced\t";}
	
	#if (-e $bam_file){print "$bam_file exists\n"}
	#if (! -e $bam_file){print "$bam_file DOESN'T exist\n"}
	
	if (-e $bam_file)
	{
		print "File found\n";
	}
	if (! -e $bam_file)
	{
		$bam_file_in_folder = &find_file($bam_file);
		
		if (-e $bam_file_in_folder)
		{
			$bam_file = $bam_file_in_folder;
			print "File found in sub directory\n";
		}
		
	}
	if (! -e $bam_file)
	{
		print "File not found\n";
		
		print "\n\n";
		print "###########################################################################################\n";
		print "$bam_file cannot be found\n\n";
		print "Correct the input file of file names, or the locations of the input BAM files and try again\n";
		print "###########################################################################################\n\n";
		close LIST;
		exit;

	}
	
	$bam_file_array[$list_count]=$bam_file;
	
	$list_count=$list_count + 1;
}

close LIST;


###############
# Name of run #
###############
&print_message("Type a name for this run (for output file prefix):","input");

$run_title = <STDIN>;
chomp $run_title;


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
print "   <2>  CanFam3nu (no unknown chromosomes)\n\n";

print "   <3>  EquCab2\n";
print "   <4>  Human\n\n";

print "   <5>  Strep. equi\n";
print "   <6>  Strep. zoo\n\n";

print "   <9>  Other\n\n";


$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3";}
if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/canfam3nu/canfam3nu.fasta"; $ref_seq_name = "canfam3nu";}
if (substr($answer,0,1) eq "3" ){$ref = "/home/genetics/equcab2/equcab2.fasta"; $ref_seq_name = "equcab2";}
if (substr($answer,0,1) eq "4" ){$ref = "/home/genetics/human/human.fasta"; $ref_seq_name = "human";}

if (substr($answer,0,1) eq "5" ){$ref = "/home/genetics/strep_equi/strep_equi.fasta"; $ref_seq_name = "s_equi";}
if (substr($answer,0,1) eq "6" ){$ref = "/home/genetics/strep_zoo/strep_zoo.fasta"; $ref_seq_name = "s_zoo";}


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

print "   <1> GATK Unified Genotyper [recommended]\n\n";
print "   <2> GATK Haplotype Caller  [still under test - not yet reliable]\n";
print "   <3> samtools mpileup       [still under test - feel free to try]\n";
print "   <4> platypus               [still under test - feel free to try]\n\n";


$answer = <STDIN>;
chomp $answer;

# Default is Unified Genotyper
if ($answer eq ""){$answer = "1"}

if (substr($answer,0,1) eq "1" ){$variant_caller = "UnifiedGenotyper"}
if (substr($answer,0,1) eq "2" ){$variant_caller = "HaplotypeCaller"}
if (substr($answer,0,1) eq "3" ){$variant_caller = "mpileup"}
if (substr($answer,0,1) eq "4" ){$variant_caller = "platypus"}

if ($variant_caller eq "")
{
	&print_message("No variant caller chosen","warning");
	exit;
}

$calls = "BOTH_TOGETHER"; # default


####################################################################
# If using GATK then ask about compressing reads using ReduceReads #
####################################################################

$reduce_reads = "no";

if ($include_reduce_reads eq "true")
{
	if (($variant_caller eq "UnifiedGenotyper") || ($variant_caller eq "HaplotypeCaller"))
	{

		&print_message("Do you want to try the experimental 'Reduce Reads' option which makes special reduced size BAM files first?","input");

		print "Do you want to reduce your BAM reads with GATK ReduceReads?";

		print "<1> Yes [ not fully tested yet ]\n";
		print "<2> No\n\n";

		$answer=<STDIN>;
		if (substr($answer,0,1) eq "1" ){$reduce_reads = "yes"} else {$reduce_reads = "no"}
		
		if ($reduce_reads eq "yes")
		{
			if ($all_reduced_files_found eq "true")
			{
				&print_message("Reduced BAM files found for each BAM file.  Do you want to use these existing reduced files?","input");

				print "<1> Yes - use existing reduced BAMs\n";
				print "<2> No  - create them again\n\n";
				
				$answer=<STDIN>;
				if (substr($answer,0,1) eq "1" ){$use_existing_reduced = "yes"} else {$use_existing_reduced = "no"}
			}
			if ($all_reduced_files_found eq "false")
			{
				&print_message("Reduced BAM files NOT found for each BAM file.  Are you happy to make these from scratch?","input");

				print "<1> Yes - make new reduced BAMs\n";
				print "<2> No  - forget about reducing them\n";
				print "<3> No  - quit\n";
				
				$answer=<STDIN>;
				if (substr($answer,0,1) eq "2" ){$reduce_reads = "no"}
				if (substr($answer,0,1) eq "3" ){exit;}
			}
		} # reduce reads yes
	}
} # include reduce reads option


if ($variant_caller eq "UnifiedGenotyper") 
{
	##################################
	# Ask if you want SNPs or INDELS #
	##################################

	&print_message("Do you want to call SNPs, Indels or both?","input");

	print "  <1> Both together [recommended]\n";
	print "  <2> Both separately\n";
	print "  <3> SNPS\n";
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

} # If use defined region

$no_of_files=$list_count - 1;


###########################################################################
# Would you like to run all the chromosomes separately                    #
# (only for Unified Genotype currently)                                   #
###########################################################################
if (($variant_caller eq "UnifiedGenotyper") && (($ref_seq_name eq "canfam3")  || ($ref_seq_name eq "canfam3nu")))
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



if (($variant_caller eq "UnifiedGenotyper") || ($variant_caller eq "HaplotypeCaller"))
{
	##################################################################################
	# Filtering preferences                                                          #
	# If you want to chnge the default values of stand_call_conf and stand_emit_conf #
	##################################################################################


	print "\n\n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
	print "Would you like to use the default values of stand_emit_conf and stand_call_conf?\n";
	print "(These are GATK parameters that control quality thresholds)                     \n";
	print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

	print "   <1>  YES\n";
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



###############################
# Ask user for memory setting #
###############################
&print_message("Please enter the memory setting for this analysis","input");

print "Memory setting (default is -Xmx4g):      ";
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


#####################
# Create file names #
#####################
$command_log_file = "$run_title"."_run_variant_caller_command_log.out";


#########################
# open Command Log file #
#########################
open (COMMAND_LOG, ">$command_log_file")|| die "Cannot create output file: $command_log_file";
print COMMAND_LOG "COMMAND LOG for run_variant_caller version $version\n";



&print_message_both("SUMMARY details","message");

&print_both("Run title:                  \t$run_title\n");
&print_both("Variant caller:             \t$variant_caller\n");
&print_both("Reference sequence:         \t$ref\n");
&print_both("GATK region string:         \t$GATK_region_string\n\n");

&print_both("GATK constants:\n\n");

&print_both("    --stand_emit_conf:      \t$stand_emit_conf\n");
&print_both("    --stand_call_conf:      \t$stand_call_conf\n");
&print_both("    --max_alt_alleles:      \t$max_alt_alleles\n\n");

&print_both("Memory setting:             \t$mem\n");
&print_both("No. of CPU threads:         \t$no_of_threads\n\n");

if ($use_defined_region eq "true")
{	
&print_both("Chromosome:                 \t$chromosome\n");
&print_both("GATK region string:         \t$GATK_region_string\n");
&print_both("UG region string:           \t$GATK_chromosome_only_string\n");
&print_both("samtools region string:     \t$samtools_region_string\n");
}

if ($reduce_reads eq "yes")
{
	&print_both("Reduce reads:                 \tYes\n");
	if ($use_existing_reduced eq "yes"){&print_both("Use existing reduced BAM files:\t\tYes\n");}
}
if ($reduce_reads ne "yes"){&print_both("Reduce reads:                 \tNo\n")}

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
	$bam_file_reduced = "$sample_name"."_reduced.bam";

	#######################
	# Show user the files #
	#######################
	print "FILE: $list_count \tBAM FILE: $bam_file\n";
		
	$input_string = $input_string." -I $bam_file";
	$input_string_reduced = $input_string_reduced." -I $bam_file_reduced";

} # list_count loop

if ($reduce_reads eq "yes"){$input_string = $input_string_reduced}
$input_string = $input_string." ";

if (($reduce_reads eq "yes") && ($use_existing_reduced eq "no"))
{
	#############################################
	# Reduce the reads in each of the BAM files #
	#############################################
	&print_message_both("Reducing the BAM files with GATK ReduceReads","message");
	printf COMMAND_LOG "  Start time:    \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $start_time)[7,2,1,0];
	
	for ($list_count=1;$list_count <=$no_of_files;$list_count++)
	{
		$bam_file = $bam_file_array[$list_count];
		$sample_name = &get_prefix ($bam_file);
		$bam_file_reduced = "$sample_name"."_reduced.bam";
		
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T ReduceReads -I $bam_file -o $bam_file_reduced  -S $GATK_validation_stringency","Reduce Reads for file $list_count/$no_of_files");	
		
		&record_input_and_output_file_size("$bam_file","$bam_file_reduced");
	} # list_count loop
}

##########################################################
#  Make a VCF file using all the BAM files in parallel   #
##########################################################

&print_message_both("Making VCF with all BAM files in parallel using $variant_caller","message");
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
			&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm SNP $input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads");	
		}
		if ($use_default_stand_values eq "true")
		{
			&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm SNP $input_string -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads");	
		}
	
		print "\n===================================================================";
		print "\nA SNP VCF file using all the BAM files in parallel has been created";
		print "\n===================================================================\n\n\n";
	} # SNPS


	if (($calls eq "INDELS") || ($calls eq "BOTH"))
	{
		&print_message ("Running UnifiedGenotyper for INDELs... ","message");
		
		$output_vcf = "$run_title"."_Indels_UG.vcf";
		
		if ($use_default_stand_values eq "false")
		{
			&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm INDEL $input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads");
		}
		if ($use_default_stand_values eq "true")
		{
			&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm INDEL $input_string --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads");
		}
				
		print "\n======================================================================";
		print "\nAn Indel VCF file using all the BAM files in parallel has been created";
		print "\n======================================================================\n\n\n";
	} # INDELS


	if ($calls eq "BOTH_TOGETHER")
	{
		&print_message ("Running UnifiedGenotyper for SNPs and INDELS together... ","message");
		
		$output_vcf = "$run_title"."_variants_UG.vcf";
		
		if ($use_default_stand_values eq "false")
		{
			&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm BOTH $input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads");
		}
		if ($use_default_stand_values eq "true")
		{
			&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm BOTH $input_string --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads");
		}

		print "\n=============================================================================";
		print "\nA SNP and Indel VCF file using all the BAM files in parallel has been created";
		print "\n=============================================================================\n\n\n";
	} # INDELS


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
		if (($ref_seq_name eq "canfam3") || ($ref_seq_name eq "canfam3nu"))
		{
			if ($chromosome_count == 39){$GATK_region_string = "-L chrX"}
		}

		
		$output_vcf = "$run_title"."_variants_UG_chr".$chromosome_count.".vcf";

		if ($use_default_stand_values eq "false")
		{
			&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm BOTH $input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads");
		}
		if ($use_default_stand_values eq "true")
		{
			&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper $GATK_region_string -glm BOTH $input_string --max_alternate_alleles $max_alt_alleles -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads");
		}
	}

}

#####################
# Haplotype Caller  #
#####################
if ($variant_caller eq "HaplotypeCaller")
{
	
	&print_message("Running HaplotypeCaller for SNPs and Indels... ","message");
	
	$output_vcf = "$run_title"."_variants_HaplotypeCaller.vcf";
	
	if ($use_default_stand_values eq "false")
	{
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller $GATK_region_string  $input_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads");	
	}
	if ($use_default_stand_values eq "true")
	{
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller $GATK_region_string  $input_string -o $output_vcf -S $GATK_validation_stringency -nct $no_of_threads");	
	}

	&print_message("A SNP and Indel VCF file using GATK HaplotypeCaller has been created","message");

} # variant_caller is Haplotype Caller


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

	print "\n\nRunning platypus in parallel for SNPs and Indels...\n\n";
		
	&run_unix_command("python2.6 /opt/platypus/Platypus.py callVariants --refFile $ref --bamFiles $list_file --output $output_vcf --bufferSize=10000 --minReads=$min_reads_platypus  $platypus_region_string");	

# python /opt/platypus/Platypus.py 
#callVariants --refFile /home/genetics/canfam3/canfam3.fasta --bamFiles SHY8_input_bam.txt


	&print_message("A SNP and Indel VCF file using platypus has been created","message");

} # variant_caller is platypus




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

&print_message("FINISHED RUNNING run_variant_caller","message");


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
	my $step = "";	
	$unix_command = $_[0];
		
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	print COMMAND_LOG "$unix_command\n";

	print("$unix_command\n\n");
	system("$unix_command");
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
