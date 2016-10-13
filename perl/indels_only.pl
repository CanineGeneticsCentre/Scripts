#!/usr/bin/perl -w


use strict;
use Getopt::Std ;
use File::Basename ;

my $indels_caller 		= "";
my $indels_raw_vcf 		= "";
my $indels_raw_bed 		= "";
my $detailed_output_bed = "";
my $command 			= "";
my $sample_name 		= "test";
my $mem 				= "-Xmx4g";
my $ref 				= "/home/genetics/canfam2/canfam2.fasta";
my $gatk_directory 		= "newgatk";
my $cleaned_sorted_bam 	= "best_align_S1_1.bam";
my $ValidationStringency	= "STRICT";
my $file_count			= 0;


print "\n\nINDEL CALLING...\n\n";

for ($file_count = 4; $file_count <=10; $file_count++)

{

	########################################################################################
	# 16 Indel genotyper   
	########################################################################################

		
		
	$indels_caller="UnifiedGenotyper";
	$sample_name = "S"."$file_count"."_1";
	$cleaned_sorted_bam= "best_align_"."$sample_name".".bam";
	
	print "\n\n>>>>>>>>>>>>   FILE COUNT: $file_count\n\n";
	print "                   Sample name: $sample_name\n";
	print "                   Bam file: $cleaned_sorted_bam\n";
	
	
	if ($indels_caller eq "UnifiedGenotyper")
	{
	
		$indels_raw_vcf = "$sample_name"."_indels_raw.vcf";

		$command =  "/usr/bin/java16 $mem -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T UnifiedGenotyper -glm INDEL -R $ref -I $cleaned_sorted_bam -o $indels_raw_vcf  -S $ValidationStringency";

		print("$command\n");
		system("$command");

	}

}
	
print "\n\nFINISHED\n";
	
	