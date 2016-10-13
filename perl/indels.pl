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

print "\n\nINDEL CALLING...\n\n";

	########################################################################################
	# 16 Indel genotyper   
	########################################################################################

		$indels_caller="UnifiedGenotyper";
		
	if ($indels_caller eq "SomaticIndelDetector")
	{
	
		$indels_raw_vcf = "$sample_name"."_indels_raw.vcf";
		$indels_raw_bed = "$sample_name"."_indels_raw.bed";
		$detailed_output_bed = "$sample_name"."_detailed_output.bed";
		
		$command =  "/usr/bin/java16 $mem -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T SomaticIndelDetector -R $ref -I $cleaned_sorted_bam -O $indels_raw_bed  -o $detailed_output_bed -ws 300 -verbose -S $ValidationStringency";

		print("$command\n");
		system("$command");

	}
	
	if ($indels_caller eq "UnifiedGenotyper")
	{
	
		$indels_raw_vcf = "$sample_name"."_indels_raw.vcf";

		$command =  "/usr/bin/java16 $mem -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T UnifiedGenotyper -glm INDEL -R $ref -I $cleaned_sorted_bam -o $indels_raw_vcf  -S $ValidationStringency";

		print("$command\n");
		system("$command");

	}
	
print "\n\nFINISHED\n";
	
	