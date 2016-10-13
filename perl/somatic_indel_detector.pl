 #!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	SOMATIC INDEL DETECTOR  				                            #     
#									                                    #
#	THIS PERL SCRIPT WILL RUN SomaticIndelDetector on all BAM files     #
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

#Various Parameters (e.g. for the Unified Genotyper)
my $stand_emit_conf				= 30;
my $stand_call_conf				= 30;

my $list_count				= 0;
my $no_of_files				= 0;

my $list_file				= "";
my $command					= "";
my $bam_file				= "";
my $vcf_merge_string		= "";
my $results_folder			= "";
my $input_string			= "";
my $run_name				= "";
my $sample_name				= "";
my $cleaned_sorted_bam		= "";
my $bam_file_in_results		= "";
my $answer					= "";
my $indels_raw_vcf			= "";
my $gatk_directory			= "newgatk";
my $ref						= "/home/genetics/canfam2/canfam2.fasta";
my $make_parallel_VCF		= "yes";
my $mem						= "-Xmx4g";
my $ValidationStringency	= "LENIENT";

my @bam_file_array			= ();

print "\n\n";
print "#######################################################################################\n";
print "#  PERL script to run the SomaticIndelDetector on all the BAM files in a list         #\n";
print "#######################################################################################\n\n";

print "The input is a file with a list of the BAM file names.\n\n";

print "Name of the file:      ";
$list_file = <STDIN>;
chomp $list_file;


print "Name of the run:      ";
$run_name = <STDIN>;
chomp $run_name;


#############################################
# Make sure the list file is in Unix format #
#############################################

$command = "dos2unix $list_file";
print("\n$command\n");
system("$command");



####################################################
# Open the list file to get the list of file names #
####################################################
open (LIST, "$list_file") || die "Cannot open $list_file";
$list_count=1;

while ($bam_file = <LIST> ) 
{
	chomp $bam_file;
	
	#########################################################
	# Add file type suffix .vcf if user hasn't added it     #
	#########################################################

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


print "\n\n";
print "Press return to continue\n\n";
$answer=<STDIN>;
if ($answer eq "q" || $answer eq "Q"){exit;}

##########################################################
#  Make a VCF file using all the BAM files in parallel   #
##########################################################

print "\n\n";
print "#---------------------------------------------------#\n";
print "#  Making Indel VCF files with SomaticIndelDetector #\n";
print "#---------------------------------------------------#\n\n";



# First get a list of the BAM files #
for ($list_count=1;$list_count <=$no_of_files;$list_count++)
{

	print "LIST COUNT: $list_count\n";
	
	$results_folder = "results_"."$run_name"."_"."$list_count";
	print "RESULTS FOLDER: $results_folder\n";
	
	$bam_file = $bam_file_array[$list_count];
	print "BAM FILE: $bam_file\n";
	
	$sample_name = &get_prefix ($bam_file);
	print "SAMPLE NAME: $sample_name\n";
	
	$cleaned_sorted_bam = "$results_folder/best_align_"."$sample_name".".bam";
	$indels_raw_vcf = "INDELS_"."$sample_name".".vcf";
	
	##############################################
	# Check if BAM file in results folder exists #
	# If not use file outside folder             #
	##############################################
	
	if (! -e "$cleaned_sorted_bam")
	{
		$cleaned_sorted_bam = $bam_file;
	}
	
	if (-e "$cleaned_sorted_bam")
	{
		$command =  "/usr/bin/java16 $mem -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T SomaticIndelDetector -R $ref -I $cleaned_sorted_bam -o $indels_raw_vcf  -ws 300 --unpaired";
		print("\n$command\n");
		system("$command");
	}
	if (! -e "$cleaned_sorted_bam")
	{
		print "\n\nFILE CAN'T BE FOUND: $cleaned_sorted_bam\n\n"; 
	}
}


print "\n\n";
print "##############\n";
print "#  FINISHED   #\n";
print "###############\n\n";

print "Output Indel VCF file have been created\n\n";
	
exit;

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

