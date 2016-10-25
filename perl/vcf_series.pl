#!/usr/bin/perl -w

#############################################################################
#									                                        #      
#	VCF SERIES       						                                #     
#									                                        #
#	THIS PERL SCRIPT WILL RUN Unified Genotype on all BAM files in series   #
#									                                        #
#############################################################################

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

my $version							= "4";

#Various Parameters (e.g. for the Unified Genotyper)
my $stand_emit_conf					= 30;
my $stand_call_conf					= 30;
my $max_alt_alleles					= 6; # This is how many alleles the INDEL caller in UnifiedGenotyper can allow

my $mem								= "-Xmx4g";
my $GATK_validation_stringency		= "LENIENT";
my $temp_dir_string					= " -Djava.io.tmpdir=javatempdir";
my $tempdir							= "";
my $list_count						= 0;
my $no_of_files						= 0;

my $list_file				= "";
my $output_file_SNPS		= "";
my $output_file_INDELS		= "";
my $command					= "";
my $bam_file				= "";
my $bam_file_in_results		= "";
my $vcf_merge_string		= "";
my $results_folder			= "";
my $input_string			= "";
my $run_name				= "";
my $chromosome						= "";
my $UG_region_string				= ""; #This is the string for the UnifiedGenotyper command line -L chr12
my $sample_name				= "";
my $cleaned_sorted_bam		= "";
my $answer					= "";
my $calls					= "";
my $species					= "";
my $gatk_directory			= "gatk";
my $ref						= "";
my $ref_seq_name			= "";
my $other_ref_sequence		= "false";

my @bam_file_array			= ();


print "\n\n\n\n\n";

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "#######################################################################################\n";
print color 'bold white';
print "   PERL script to run Unified Genotyper on all the BAM files consecutively (in series)   \n";
print color 'bold magenta';
print "#######################################################################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This program produces one SNP VCF file and one Indel VCF file for each BAM file\n";

print color 'reset';
print "\n\n";

print "The input is a file with a list of the BAM file names.\n\n";


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
		if (! -e "$list_file"){print "\nFile doesn't exist. Try again...    \n\n";}
	}
		
}


#############################################
# Make sure the list file is in Unix format #
#############################################

$command = "dos2unix $list_file";
print("\n$command\n");
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
	
	#########################################################
	# Add file type suffix .bam if user hasn't added it     #
	#########################################################

	if (index($bam_file,".bam") == -1 ){$bam_file = $bam_file.".bam"}
	
	$bam_file_array[$list_count]=$bam_file;
	
	#print "$list_count\t$bam_file\n";
	
	#if (-e $bam_file){print "$bam_file exists\n"}
	#if (! -e $bam_file){print "$bam_file DOESN'T exist\n"}
	
	if (! -e $bam_file)
	{
		
		print "\n\n";
		print "###########################################################################################\n";
		print "$bam_file cannot be found\n\n";
		print "Correct the input file of file names, or the locations of the input BAM files and try again\n";
		print "###########################################################################################\n\n";
		close LIST;
		exit;

	}
	$list_count=$list_count + 1;
}

close LIST;

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

if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3";}
if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/canfam2/canfam2.fasta"; $ref_seq_name = "canfam2";}
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

##################################
# Ask if you want SNPs or INDELS #
##################################

print "\nWhich do you want to call:      \n\n";

print "  <1> SNPS\n";
print "  <2> Indels\n";
print "  <3> Both\n\n";

$answer = <STDIN>;
chomp $answer;

if ($answer eq "1"){$calls = "SNPS"}
if ($answer eq "2"){$calls = "INDELS"}
if ($answer eq "3"){$calls = "BOTH"}
if ($answer eq ""){$calls = "BOTH"}



print "\nChromosome:      ('return' for whole genome)";
$chromosome = <STDIN>;
chomp $chromosome;

if (index($chromosome,"chr") == -1 ){$chromosome = "chr".$chromosome}

$UG_region_string = "-L $chromosome";

###################################################################
# If user has not specified a chromosome, don't add the -L option #
###################################################################
if ($chromosome eq "chr")
{
	print "\nNo chromosome specified. Use whole genome\n\n";
	$UG_region_string = "";
	$chromosome = "None";
}

$no_of_files=$list_count - 1;

###################
# SUMMARY DATA    #
###################

print "==================\n";
print "VCF_SERIES SUMMARY \n";
print "==================\n\n";

print "UnifiedGenotyper constants:\n\n";

print "    --stand_emit_conf:\t\t\t$stand_emit_conf\n";
print "    --stand_call_conf:\t\t\t$stand_call_conf\n";
print "    --max_alt_alleles:\t\t\t$max_alt_alleles\n\n";

print "Memory setting: \t\t\t$mem\n\n";
print "GATK directory: \t\t\t$gatk_directory\n";
print "Chromosome:     \t\t\t$chromosome\n";
print "UG region:      \t\t\t$UG_region_string\n\n";

print "There are $no_of_files BAM files:\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print "   File $list_count	\t$bam_file_array[$list_count]\n";
}

print "\n\n";
print "############################\n";
print "# Press return to continue #\n";
print "############################\n\n";
$answer=<STDIN>;
chomp $answer;
if ($answer eq "q" || $answer eq "Q"){exit;}


#########################
# open Command Log file #
#########################
open (COMMAND_LOG, ">vcf_series_command_log.out")|| die "Cannot create output file: vcf_series_command_log.out";
print COMMAND_LOG "COMMAND LOG for VCF_PARALLEL\n\n";


##########################################################
#  Make VCF files using all the BAM files in series      #
##########################################################

print "\n\n";
print "#------------------------------------------#\n";
print "#  Making VCF with all BAM files in series #\n";
print "#------------------------------------------#\n\n";


############################################################
# Run through all the BAM files in order to make SNP files #
############################################################

print "no_of_files: $no_of_files\n\n";

for ($list_count=1;$list_count <=$no_of_files;$list_count++)
{

	$bam_file = $bam_file_array[$list_count];
	
	print "BAM file: $bam_file\n\n";
	
	if (-e $bam_file)
	{
		$bam_file = $bam_file;
		$sample_name = &get_prefix ($bam_file);
		$output_file_SNPS = "SNPS"."_"."$sample_name".".vcf";
		$output_file_INDELS = "INDELS"."_"."$sample_name".".vcf";
		
	}
	
	if (($calls eq "SNPS") || ($calls eq "BOTH"))
	{
		print "\n\n";
		print "========================\n";
		print "SNPS:  File $list_count/$no_of_files\n";
		print "========================\n\n";
		
		print "BAM file: \t$bam_file\n";
		print "OUTPUT file: \t$output_file_SNPS\n\n";
		if (! -e $bam_file){print "File $bam_file can't be found\n\n";}
		
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper -glm SNP -o $output_file_SNPS -I $bam_file $UG_region_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf -S $GATK_validation_stringency","$list_count");

		&record_output_file_size("$output_file_SNPS");
		
		print "\n===================================================================";
		print "\nA SNP VCF file using BAM file $bam_file has been created";
		print "\n===================================================================\n\n\n";
		
	} # SNPS
	
	if (($calls eq "INDELS") || ($calls eq "BOTH"))
	{
		print "\n\n";
		print "========================\n";
		print "INDELS:  File $list_count/$no_of_files\n";
		print "========================\n\n";
		
		print "BAM file: \t$bam_file\n";
		print "OUTPUT file: \t$output_file_INDELS\n\n";
		if (! -e $bam_file){print "File $bam_file can't be found\n\n";}
		
		
		&run_unix_command("java $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper -glm INDEL -o $output_file_INDELS -I $bam_file $UG_region_string -stand_emit_conf $stand_emit_conf -stand_call_conf $stand_call_conf --max_alternate_alleles $max_alt_alleles -S $GATK_validation_stringency","$list_count");
		
		&record_output_file_size("$output_file_INDELS");
		
		print "\n===================================================================";
		print "\nAn INDEL VCF file using BAM file $bam_file has been created";
		print "\n===================================================================\n\n\n";

	} # INDELS
	
	
} # End of SNP loop

close (COMMAND_LOG);
		
print "\n===========================================================================";
print "\nSNP and INDEL VCF files using all the BAM files in series have been created";
print "\n===========================================================================\n\n\n";

	
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
	$step = $_[1];
		
	print("$unix_command\n\n");
	system("$unix_command");
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	print COMMAND_LOG "File: $step/$no_of_files\n";
	print COMMAND_LOG "$unix_command\n";

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
