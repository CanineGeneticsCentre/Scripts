#!/usr/bin/perl -w

##################################################################################
#									                                             #      
#	vcf2excel			                                                         #     
#									                                             #
#	This PERL script converts VCF files to a format ready for NGS SNP Handler	 #
#									                                             #
##################################################################################

use strict;
use File::Basename;
use Term::ANSIColor;
use List::Util qw[min max];

# VERSION OF SOFTWARE #
my $version							= "2";

my $show							= "false"; # for debugging
my $show_seg						= "";

# Constants
my $merge_threshold					= 10; # maximum distance apart for merging nearby indels

# File names
my $vcf_file						= "";
my $command_log						= ""; 
my $output_file						= "";
my $output_file_filtered			= ""; # output file filtered for best segregations and effects
my $indels_merged_file				= ""; # output file indels only
my $alleles_file					= ""; # Equivalent to 'SNPs' on NGS SNP Handler
my $disease_status_file				= ""; 
my $vep_file						= ""; 
my $simplified_indels_output_file	= ""; # Indels output file for use as input for merging indels

my $temp_indels_file				= ""; # TEMPORARY!!

#Boolean - 'true' or 'false'
my $snpEff_data_found				= "false";
my $VEP_data_found					= "false";
my $start_storing 					= "false";
my $write_to_filtered				= "false"; # decides whether each line is written to the filtered output file
my $include_extra_columns			= "false"; # This is fixed as false (but could be changed here to 'true' to get extra output columns)
my $vep_id_found					= "false"; # whether a VEP position is found at any VCF position

#Strings
my $single_line						= ""; # variable for readsing single line of text files
my $check_line						= "";
my $variant_caller					= ""; # e.g. 'UnifiedGenotyper', 'HaplotypeCaller etc'
my $answer							= "";
my $default_vep_file				= ""; 
my $char							= "";
my $files_string					= "";
my $chrom_string					= "";
my $prefix							= "";
my $array_string					= ""; # for debugging
my $sample_name						= "";
my $chromosome						= "";
my $position						= "";
my $position_plus_one				= "";
my $position_check					= "";
my $REF_base						= "";
my $ALT_base_overall				= "";
my $QUAL							= "";
my $FILTER							= "";
my $GT_string 						= "";
my $GT_score						= "";
my $AD_string 						= "";
my $DP_genotype 					= "";
my $GQ_string 						= "";
my $PL_string 						= "";
my $myString						= "";
my $allele							= "";
my $allele_check					= "";
my $merged_allele					= "";
my $allele_number_1					= ""; # as in the genotype 0/1, 1/1 etc
my $allele_number_2					= "";
my $allele_letter_1					= ""; # as in the genotype AA, AB etc
my $allele_letter_2					= "";
my $allele_base_1					= ""; # actual base of the genotype, A,C, G or T
my $allele_base_2					= "";
my $allele_base_check_1				= ""; # actual base of the genotype, A,C, G or T when looking for merged indels
my $allele_base_check_2				= "";
my $allele_base_merged_1			= ""; # actual base of the genotype, A,C, G or T when looking for merged indels
my $allele_base_merged_2			= "";
my $mark_missing_as_X				= ""; # Controls when genotype has not been called
my $score_X_as_REF					= ""; # If you do use X, do you count these as REF for scoring?
my $num								= "";
my $comment							= "";
my $minor_allele					= "";
my $allele_A_original				= "";
my $allele_B_original				= "";
my $allele_C_original				= "";
my $allele_A						= "";
my $allele_B						= "";
my $allele_C						= ""; # more than 3 alleles isn't dealt with.
my $main_affected_allele			= ""; # if variant is SNP then this is the actual base
my $main_affected_allele_check		= ""; 
my $second_affected_allele			= "";
my $main_normal_allele				= ""; # IS THIS NEEDED?
my $genotype						= "";
my $genotype_check					= "";
my $original_genotype				= "";
my $merged_genotype					= "";
my $swapped							= "";
my $reassigned						= ""; # flag to check whether alleles A and B have been reassigned
my $disease_status_line				= "";
my $disease_status					= "";
my $sample_name_from_status_file 	= "";
my $effect_predictor				= ""; # variant_effect_predictor or snpEff
my $variant_type					= ""; # snp, indel (or both?)
my $variation_id					= ""; # VEP variant ID
my $vep_location					= ""; # Second column in VEP output 18:25134510
my $vep_chr							= ""; # CHR split from VEP vep_location in the form 7:21821942
my $vep_pos							= ""; # POS split similarly
my $vep_pos_minus_one				= ""; # To cover VEP renumbering position

#snpEff and VEP variables
my $snpEff_effect					= ""; # The main name for the snpEff effect
my $vep_effect						= ""; # VEP effect or consequence
my $snpEff_effect_full_string		= ""; # holds string of effects from snpEff
my $vep_effect_full_string			= ""; # holds string of effects from VEP
my $snpEff_effect_string			= ""; # holds individual effects from snpEff
my $snpEff_effect_impact 			= "";
my $snpEff_effect_codon_change		= "";
my $snpEff_effect_gene_name 		= "";
my $snpEff_effect_transcript 		= "";
my $rest_of_string					= "";
my $snpEff_results_string			= ""; # String which is written to output file for consequence column in Excel
my $vep_results_string				= ""; # String which is written to output file for consequence column in Excel
my $snpEff_gene_string				= ""; # String which is written to output file for gene column in Excel
my $vep_id							= "";
my $vep_id_minus_one					= "";
my $vcf_id							= "";
my $vcf_id_plus_one					= "";
my $consequence						= "";
my $consequence_check				= "";
my $consequence_merged				= "";

#Other
my $segregation_score_threshold		= 0;
my $line_count						= 0;
my $col_count						= 0;
my $check_count						= 0;
my $vep_line_count					= 0;
my $pos								= 0;
my $open_bracket_pos				= 0;
my $close_bracket_pos				= 0;
my $sample_count					= 0;
my $new_sample_count				= 0;
my $sample_test_count				= 0;
my $chrom_line_array_size			= 0;
my $array_count						= 0;
my $no_of_samples_UG_line			= 0;
my $no_of_samples_chrom_line		= 0;
my $no_of_samples					= 0; # final number of samples
my $file_count						= 0;
my $no_of_data_columns				= 0;
my $no_alleles						= 0;
my $array_size_1					= 0;
my $array_size_7					= 0;
my $array_size_8					= 0;
my $array_size_9					= 0;
my $slash_pos						= 0;
my $homozygosity_ratio				= 0;
my $homozygosity_ratio_check		= 0;
my $homozygosity_ratio_merged		= 0;
my $no_of_alt_alleles				= 0;
my $no_of_alt_alleles_check			= 0;
my $no_of_alt_alleles_merged		= 0;
my $base_count						= 0;
my $base_count_1					= 0;
my $base_count_2					= 0;
my $col_1							= 0;
my $col_2							= 0;
my $status_line_array_size			= 0;
my $disease_status_count			= 0;
my $name_match_count				= 0; # checks number of times each name in VCF matches with name in Disease Status files
my $total_match_count				= 0; # checks how many names in VCF matche with name in Disease Status files
my $allele_count					= 0;
my $max_allele_count				= 0;
my $second_max_allele_count			= 0;
my $no_of_snpEff_effects			= 0;
my $no_of_vep_effects				= 0;
my $effect_count					= 0;
my $effect_score					= 0;
my $effect_score_check				= 0; # When checking for merged indels
my $effect_score_merged				= 0; # 
my $vep_effect_score				= 0;
my $snpEff_effect_score				= 0;
my $max_effect_score   				= 0; # If VEP and snpEff are both used
my $max_vep_effect_score   			= 0;
my $max_snpEff_effect_score   		= 0;
my $missing_genotypes_total			= 0;
my $missing_genotypes_percent		= 0;
my $genotypes_counted_total			= 0;
my $genotypes_counted_per_sample	= 0;
my $vep_id_count					= 0;
my $vcf_id_count					= 0;
my $filter_count_no_of_alleles		= 0; # Counts the number of SNPs added to the filtered results because of no of alleles
my $filter_count_seg_score			= 0; # Counts the number of SNPs added to the filtered results because of segreagation score best
my $filter_count_effect_score		= 0; # Counts the number of SNPs added to the filtered results because of effect score
my $vep_variant_count				= 0;
my $vcf_variant_count				= 0;

my $no_of_alleles					= 0;
my $no_of_alleles_check				= 0;
my $max_no_of_alleles				= 0;
my $second_no_of_alleles			= 0;
my $no_of_indels					= 0; # number of indels in the separate indel file (or array)
my $vep_id_found_count				= 0; # count of when a VEP position is found at any VCF position
my $vep_id_not_found_count			= 0; # count of when a VEP position is NOT found at any VCF position
my $vep_position_matches_count		= 0; # count of when VCF position matches exactly to a VEP position
my $vep_position_minus_one_count	= 0; # count of times you have to use minus_one to get VEP to match to VCF
my $vcf_id_1_count					= 0; # Count number of times we have to use the 18456646_1 system for two positions the same
my $vcf_id_2_count					= 0; # Count number of times we have to use the 18456646_2 system for two positions the same


# AA_AB_BB         	Recessive:  Affecteds AA, carriers AB, normals BB
# AA_AB_ABorBB     	Recessive, some normals may be carriers:  Affecteds AA, carriers AB, normals AB or BB
# AA_ABorBB_BB      	Recessive, some carriers may be normals:  Affecteds AA, carriers AB or BB, normals BB
# AA_ABorBB_ABorBB 	Recessive:  Affecteds AA, carriers AB or BB, normals AB or BB
# AAorAB_BB        	Dominant:   Affecteds AA or AB, normals BB
# A_B              	Additive:   Affecteds Awith As, Normals with Bs
# AA_BB             Strict:   Affecteds AA, Normals BB (no score for AB)


#Segregation scores
my $segregation_score_1				= 0; # AA_AB_BB         Recessive:  Affecteds AA, carriers AB, normals BB
my $segregation_score_2				= 0; # AA_AB_ABorBB     Recessive, some normals may be carriers:  Affecteds AA, carriers AB, normals AB or BB
my $segregation_score_3				= 0; # AA_ABorBB_BB     Recessive, some carriers may be normals:  Affecteds AA, carriers AB or BB, normals BB
my $segregation_score_4				= 0; # AA_ABorBB_ABorBB Recessive:  Affecteds AA, carriers AB or BB, normals AB or BB
my $segregation_score_5				= 0; # AAorAB_BB        Dominant:   Affecteds AA or AB, normals BB
my $segregation_score_6				= 0; # A_B              Additive:   Affecteds Awith As, Normals with Bs
my $segregation_score_7				= 0; # AA_BB            Strict:   Affecteds AA, Normals BB (no score for AB)
my $segregation_score_1B			= 0; 
my $segregation_score_2B			= 0; 
my $segregation_score_3B			= 0; 
my $segregation_score_4B			= 0; 
my $segregation_score_5B			= 0; 
my $segregation_score_6B			= 0; 
my $segregation_score_7B			= 0; 
my $segregation_score_best_1		= 0; 
my $segregation_score_best_2		= 0; 
my $segregation_score_best_3		= 0; 
my $segregation_score_best_4		= 0; 
my $segregation_score_best_5		= 0; 
my $segregation_score_best_6		= 0; 
my $segregation_score_best_7		= 0; 
my $segregation_score_best			= 0; # Best overall
	
#Statistics	
my $A_count_all						= 0;
my $B_count_all						= 0;
my $MAF								= 0;
my $A_count_affected				= 0;
my $B_count_affected				= 0;
my $C_count_affected				= 0;
my $A_count_normal					= 0;
my $B_count_normal					= 0;
my $C_count_normal					= 0;
my $AA_count_affected				= 0;
my $AB_count_affected				= 0;
my $BB_count_affected				= 0;
my $AA_AB_count_affected			= 0;
my $AB_BB_count_affected			= 0;
my $AA_count_normal					= 0;
my $AB_count_normal					= 0;
my $BB_count_normal					= 0;
my $AA_AB_count_normal				= 0;
my $AB_BB_count_normal				= 0;
my $affected_dom_count				= 0;
my $normal_dom_count				= 0;
my $HH_count_affected				= 0;
my $HH_count_normal					= 0;
my $affected_homozygosity			= 0;
my $normal_homozygosity				= 0;
my $no_affected_samples				= 0;
my $no_normal_samples				= 0;
my $no_carrier_samples				= 0;
my $no_samples_in_status_file		= 0;
my $source_col						= 0;
my $three_alleles_count_snps		= 0;
my $three_alleles_count_indels		= 0;
my $four_alleles_count_snps			= 0;
my $four_alleles_count_indels		= 0;
my $no_of_lines_vep					= 0;
my $main_affected_homozygous_genotype = "";


# New for vcf2excel_single
my $genotype_string					= ();
my $allele_1_for_writing			= (); # The base to write for the first allele: A, C, G or T
my $allele_2_for_writing			= (); # The base to write for the second allele: A, C, G or T
my $missing_genotypes_count			= (); # no of missing genotypes for THE sample
my $genotypes_counted				= (); # no of counted genotypes for each sample

# Arrays
my @input_file_array				= (); # names of input files (on UnifiedGenotyper line)
my @sample_name_array				= (); # names of the samples (from the #CHROM line)
my @chrom_line_array				= (); # array of items on line with #CHROM in it.
my @disease_status_sample_array 	= (); # array of sample names from disease status file
my @disease_status_array 			= (); # array of disease_statuses from disease status file
my @myArray1						= ();
my @myArray7						= (); # contains snpEff results (if they are in the VCF file)
my @myArray8						= ();
my @myArray9						= ();
my @ALT_allele_array				= ();

my @base_array						= (); # array of single bases (twice as many as samples)
my @base_array_check				= ();
my @base_orig_array					= (); # array of single bases original (in case they are 'X')
my @affection_array					= (); # 'affected', 'carrier' or 'normal'


my @vep_array						= (); # Array to whole elements of vVEP variation_id e.g. 7_21821942_G/A
my @allele_count_array				= ();
my @indel_file_array				= ();
my @merged_allele_array				= ();
my @merged_genotype_string			= ();
my @vep_file_array					= ();


#Arrays in original input order
my @genotype_string_input_order		= (); # genotypes in original input order (@genotype_string is in output order or disease_status order)
my @sample_status_array_input_order	= (); # disease statuses of samples in original input order

#New arrays for PERL script
#my @QUAL_array						= ();
my @DP_array						= ();
my @GT_array						= ();
my @PL_array						= ();
my @item							= ();
my @item_check						= ();
my @status_line_array				= ();

my @source_column_array				= (); # stores destination column for each sample (to make order Affected -> Carrier -> Normal)
my @sample_status_array				= (); # Disease statuses attached to samples in VCF file (may be different order than Disease Status file)
my @snpEff_effects_array			= (); # Stores effects from snpEff
my @vep_effects_array				= (); # Stores effects from VEP
my @snpEff_items_array				= (); # Holds the different items in a snpEff effect, e.g. gene name

# Hash table for indel letters "a", "b", "c" etc (looked up with actual sequence e.g. might be TAAA)
my %indel_hash						= ();
my %vep_hash						= ();
my %vep_hash_minus_one				= ();
my %merged_hash						= ();

print color 'reset';


print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      vcf2excel_single            \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This PERL script converts a single-column VCF file to a format ready for NGS SNP VIEWER\n\n";
print "    It creates the final alleles page directly from the VCF file.\n";
print "    It also add the variant annotations (effects) but it requires that\n";
print "    variant_effect_predictor or snpEff has been run on this file before you run vcf2excel.\n\n";
print color 'reset';


if (index($version,"t") > -1)
{
	&print_message("THS IS ONLY A TEST VERSION","warning");
}


 ######   ######## ########    ######## #### ##       ########  ######  
##    ##  ##          ##       ##        ##  ##       ##       ##    ## 
##        ##          ##       ##        ##  ##       ##       ##       
##   #### ######      ##       ######    ##  ##       ######    ######  
##    ##  ##          ##       ##        ##  ##       ##             ## 
##    ##  ##          ##       ##        ##  ##       ##       ##    ## 
 ######   ########    ##       ##       #### ######## ########  ######  


&print_message("Please input the name of the SINGLE-SAMPLE VCF file","input");

until (-e $vcf_file)
{
	print "   VCF file:  ";
	$vcf_file = <STDIN>;
	chomp $vcf_file;
	if ($vcf_file eq ""){$vcf_file = "test.vcf"}  # TEMP!!!
	if ($vcf_file eq "ls"){print "\n";system ("ls *.vcf");print "\n"}
	if ($vcf_file ne "ls"){if (! -e $vcf_file){print "\n\n>>>>>>>>  File $vcf_file not found.  Try again.  <<<<<<<<\n\n";}}
}

print "\n\n";


################################################
# Ask user how the SNP Effects are to be added #
################################################
until ($effect_predictor ne "")
{
	&print_message("How are the consequences/effects of the mutations to be added?","input");

	print "   <1> Variant Effect Predictor (VEP) (missing consequences should be fixed, but please check)\n";
	print "   <2> snpEff\n\n";

	$answer = <STDIN>;
	chomp $answer;

	if ($answer eq "" ){$effect_predictor = "VEP"}
	if (substr($answer,0,1) eq "1" ){$effect_predictor = "VEP"}
	if (substr($answer,0,1) eq "2" ){$effect_predictor = "snpEff"}

	print "\n\n";
}


#################################
# Create names for output files #
#################################
$prefix = &get_prefix ($vcf_file);

$output_file = $prefix."_for_excel.txt";
$output_file_filtered = $prefix."_filtered_for_excel.txt";

$simplified_indels_output_file = $prefix."_indels_only.txt";
$indels_merged_file = $prefix."_indels_merged.txt";
$temp_indels_file = "temp_indels_only_file.txt"; # temporary

$alleles_file = $prefix."_alleles.txt";
$command_log = $prefix."_log.out";
$default_vep_file = $prefix."_vep.out";


#########################
# Open command log file #
#########################
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";


####################
# Get VEP filename #
####################
if ($effect_predictor eq "VEP")
{
	&print_message("The VEP file from the variant_effect_predictor can be added during the running of this program.","message");
	&print_message("Please input the name of the VEP file (variant effect predictor)","input");

	until (-e $vep_file)
	{
		print "   VEP file (default = $default_vep_file):  ";
		$answer = <STDIN>;
		chomp $answer;

		if ($answer eq ""){$vep_file = $default_vep_file} else {$vep_file = $answer}

		if ($vep_file eq ""){$vep_file = "snps_vep.out"}  # TEMP!!!
		if ($vep_file eq "ls"){print "\n";system ("ls *vep*");print "\n"}
		if ($vep_file ne "ls"){if (! -e $vep_file){print "\n\n>>>>>>>>  File $vep_file not found.  Try again.  <<<<<<<<\n\n";}}
	}
} # VEP


##########
# snpEff #
##########
print "\n\n";
if ($effect_predictor eq "snpEff")
{
	&print_message("snpEff adds its results to the VCF file so to use the data you should have already processed this VCF file.","message");

	&print_message("Has this VCF file already been processed by snpEff?","input");

	print "   <1> Yes\n";
	print "   <2> No\n\n";

	$answer = <STDIN>;
	chomp $answer;

	if (substr($answer,0,1) eq "1")
	{
		print "\n\tOK. This program will use the snpEff data\n";
	}
	if (substr($answer,0,1) eq "2")
	{
		&print_message ("OK. Process it with snpEff and start again","warning");
		print "(The easiest way is to use the perl script 'run_snpEff.pl')\n\n";
		exit;
	}
} # snpEff

print "\n";



#########################################################
# Ask about how you want to deal with missing genotypes #
#########################################################
$mark_missing_as_X = "true";

&print_message("How do you want to mark the allele if a genotype has not been called?","input");

print "   <1> Mark with X (recommended: you can still choose to score it as REF)\n";
print "   <2> Use Reference base\n\n";

$answer = <STDIN>;
chomp $answer;

if ($answer eq ""){$mark_missing_as_X = "true"} # Default
if (substr($answer,0,1) eq "1" ){$mark_missing_as_X = "true"}
if (substr($answer,0,1) eq "2" ){$mark_missing_as_X = "false"}

$score_X_as_REF = "true";

if ($mark_missing_as_X eq "true")
{
	&print_message("When doing the scoring, do you want to count Xs as the REF allele?","input");

	print "   <1> Count X as REF allele (recommended)\n";
	print "   <2> Don't count Xs in scoring at all.\n\n";

	$answer = <STDIN>;
	chomp $answer;
	
	if ($answer eq ""){$score_X_as_REF = "true"} # Default
	if (substr($answer,0,1) eq "1" ){$score_X_as_REF = "true"}
	if (substr($answer,0,1) eq "2" ){$score_X_as_REF = "false"}
}

########  ########    ###    ########     ##     ## ######## ########     ######## #### ##       ######## 
##     ## ##         ## ##   ##     ##    ##     ## ##       ##     ##    ##        ##  ##       ##       
##     ## ##        ##   ##  ##     ##    ##     ## ##       ##     ##    ##        ##  ##       ##       
########  ######   ##     ## ##     ##    ##     ## ######   ########     ######    ##  ##       ######   
##   ##   ##       ######### ##     ##     ##   ##  ##       ##           ##        ##  ##       ##       
##    ##  ##       ##     ## ##     ##      ## ##   ##       ##           ##        ##  ##       ##       
##     ## ######## ##     ## ########        ###    ######## ##           ##       #### ######## ########

#################################################################################
# If variants are from the variant_effect_predictor VEP then read this file now #
#################################################################################

$line_count = 0;

if ($effect_predictor eq "VEP")
{
	&print_message("Reading VEP file $vep_file","message");

	open (VEP, "$vep_file") || die "Cannot open $vep_file";

	@vep_file_array = <VEP>;
	$no_of_lines_vep = scalar @vep_file_array;
	close VEP;

	$line_count = 0;
	$start_storing = "false";
	$array_count = 0;

	for ($line_count = 0; $line_count < $no_of_lines_vep; $line_count++)
	{
		$single_line = $vep_file_array[$line_count];

		chomp $single_line;
		&chomp_all ($single_line);

		if (($line_count % 10000) == 0 ){print "Line: $line_count/$no_of_lines_vep\n";}

		#Split line at tabs
		@item = split (/\t/,$single_line);
		$array_size_1 = scalar @item;
		
		if ($start_storing eq "true")
		{
			$VEP_data_found	= "true";

			$variation_id 	= $item[0]; # This is the 'Uploaded_variation' column in the VEP output file' e.g. 13_52901557_G/-

			$vep_location   = $item[1]; # This is the 'Location' column.      e.g. 39009362 or 39073821-39073822

			$vep_effect 	= $item[6]; # This is the 'Consequence' column 

			$vep_variant_count = $vep_variant_count + 1;

			##############################################################
			# Deal with VEP positions that are in the form 185355-185356 #
			##############################################################
			if (index($vep_location,"-") > -1)
			{

				$vep_location = substr ($vep_location, 0,index($vep_location,"-"));

			}


			#######################################################################
			# Get CHR and POS from vep_location                                   #
			# NOTE: for Indels the Position can be different than in the VCF file #
			#######################################################################
			@vep_array = split(":",$vep_location);
			$vep_chr = $vep_array[0];
			$vep_pos = $vep_array[1];


			#######################
			# convert X to 39 etc #
			#######################
			if ($vep_chr eq "X"){$vep_chr = 39}

			#if ($vep_pos eq "38900305"){$show = "true"}
			if ($vep_pos eq "38900306xx"){$show = "true"}else{$show = "false"}
	
			###########################################################
			# Check if chr is numeric and only do this bit if it is   #
			# This excludes the unknown extra chromosomes and X and M #
			###########################################################
			if ($vep_chr =~ /^\d+?$/) 
			{

				# convert chr1 to chr01
				if ($vep_chr < 10){$vep_chr = "0"."$vep_chr";}

				################################################################################
				# Make vep_pos_minus_one to cover cases where VEP increases the position by one #
				# (more likely to need VCF position plus one in practice)                      #
				################################################################################
				$vep_pos_minus_one = $vep_pos - 1;

				$vep_id = $vep_chr."_".$vep_pos;
				$vep_id_minus_one = $vep_chr."_".$vep_pos_minus_one;
				
				if ($show eq "true")
				{
					print "\nvep_id:           $vep_id\n";
					print "vep_id_minus_one: $vep_id_minus_one\n";
				}
				
				###########################################################
				# Store in hash table (faster than array)                 #
				# If there are two effects for one position, store the    #
				# subsequent ones one as 12_12134562_1, 12_12134562_2 etc #
				#                                                         #
				# This will be a SNP and an Indel at the same position    #
				###########################################################
				$vep_id_count = 0;
				
				if ($show eq "true")
				{
					if (defined $vep_hash{$vep_id})
					{
						print "vep_hash($vep_id) defined as $vep_hash{$vep_id}\n";
					}
					if (not defined $vep_hash{$vep_id})
					{
						print "vep_hash($vep_id) not defined\n";
					}
					if (defined $vep_hash_minus_one{$vep_id_minus_one})
					{
						print "vep_hash_minus_one($vep_id) defined as $vep_hash{$vep_id}\n";
					}
					if (not defined $vep_hash_minus_one{$vep_id_minus_one})
					{
						print "vep_hash_minus_one($vep_id_minus_one) not defined\n";
					}
				}
				
				if (not defined $vep_hash{$vep_id})
				{
					$vep_hash{$vep_id} = $vep_effect;
				}
				else
				{
					while(defined $vep_hash{$vep_id})
					{
						$vep_id_count = $vep_id_count + 1;
						$vep_id = $vep_id."_".$vep_id_count;
					}
					$vep_hash{$vep_id} = $vep_effect;
				}

				###################################################
				# Back up hash in case first one doesn't match up #
				###################################################
				$vep_id_count = 0;
				if (not defined $vep_hash_minus_one{$vep_id_minus_one})
				{
					$vep_hash_minus_one{$vep_id_minus_one} = $vep_effect;
				}
				else
				{
					while(defined $vep_hash_minus_one{$vep_id_minus_one})
					{
						$vep_id_count = $vep_id_count + 1;
						$vep_id_minus_one = $vep_id_minus_one."_".$vep_id_count;
					}
					$vep_hash_minus_one{$vep_id_minus_one} = $vep_effect;
				}

				
				if ($show eq "true")
				{
					if (defined $vep_hash{$vep_id})
					{
						print "  vep_hash($vep_id) defined as $vep_hash{$vep_id}\n";
					}
					if (not defined $vep_hash{$vep_id})
					{
						print "  vep_hash($vep_id) not defined\n";
					}
					if (defined $vep_hash_minus_one{$vep_id_minus_one})
					{
						print "  vep_hash_minus_one($vep_id_minus_one) defined as $vep_hash{$vep_id}\n";
					}
					if (not defined $vep_hash_minus_one{$vep_id_minus_one})
					{
						print "  vep_hash_minus_one($vep_id_minus_one) not defined\n";
					}
				}
				
			}#if chr is numeric
			
		}

		if ($array_size_1 > 12){$start_storing = "true"}

	} # end of reading VEP file

	close VEP;

	&print_message("VEP file $vep_file read successfully","message");

} # Effect predictor is VEP


###########################
# Open the input VCF file #
###########################
open (VCF, "$vcf_file") || die "Cannot open $vcf_file";
$line_count = 0;
$vep_line_count = 1;
$start_storing = "false";

while ($single_line = <VCF>) 
{
	chomp $single_line;
	&chomp_all ($single_line);

	$line_count = $line_count + 1;
	
	###################################################################
	# This checks the first line to see that it looks like a VCF file #
	###################################################################
	if($line_count == 1)
	{
		if (index($single_line,"VCF") > -1)
		{
			&print_message("File format check","message");
			print "This looks like a VCF file\n\n";
			
			&pause;
		}		
	}
	

	 ######  ########    ###    ########  ########     ######  ########  #######  ########  #### ##    ##  ######      ##     ##  ######  ######## 
	##    ##    ##      ## ##   ##     ##    ##       ##    ##    ##    ##     ## ##     ##  ##  ###   ## ##    ##     ##     ## ##    ## ##       
	##          ##     ##   ##  ##     ##    ##       ##          ##    ##     ## ##     ##  ##  ####  ## ##           ##     ## ##       ##       
	 ######     ##    ##     ## ########     ##        ######     ##    ##     ## ########   ##  ## ## ## ##   ####    ##     ## ##       ######   
	      ##    ##    ######### ##   ##      ##             ##    ##    ##     ## ##   ##    ##  ##  #### ##    ##      ##   ##  ##       ##       
	##    ##    ##    ##     ## ##    ##     ##       ##    ##    ##    ##     ## ##    ##   ##  ##   ### ##    ##       ## ##   ##    ## ##       
	 ######     ##    ##     ## ##     ##    ##        ######     ##     #######  ##     ## #### ##    ##  ######         ###     ######  ## 

	######################################################
	# Only do this next bit if 'start_storing is True    #   
	# i.e. if the line containing #CHROM has been found  #
	######################################################
	if ($start_storing eq "true")
	{
		
		#######################################
		# Set various things to blank or zero #
		#######################################
		
		# Alleles are set to Q so we can easily see if they have not been assigned correctly
		$allele_A = "Q";
		$allele_B = "Q";
		$allele_C = "Q";
		
		$main_affected_allele = "";
		$main_normal_allele = "";
		
		# Allele counts
		$A_count_affected = 0;
		$B_count_affected = 0;
		$C_count_affected = 0;
		
		$A_count_normal = 0;
		$B_count_normal = 0;
		$C_count_normal = 0;
		
		# Genotype counts
		$AA_count_affected = 0;
		$AB_count_affected = 0;
		$BB_count_affected = 0;
		$HH_count_affected = 0;
		
		$AA_count_normal = 0;
		$AB_count_normal = 0;
		$BB_count_normal = 0;
		$HH_count_normal = 0;
		
		# Homozygosities
		$affected_homozygosity = 0;
		$normal_homozygosity = 0;
		
		# Dominant counts
		$AA_AB_count_affected = 0;
		$AB_BB_count_affected = 0;
		$AA_AB_count_normal = 0;
		$AB_BB_count_normal = 0;
		$affected_dom_count = 0;
		
		$homozygosity_ratio = 0;
		
		# Segregation scores
		$segregation_score_1 = 0;
		$segregation_score_2 = 0;
		$segregation_score_3 = 0;
		$segregation_score_4 = 0;
		$segregation_score_5 = 0;
		$segregation_score_6 = 0;
		$segregation_score_7 = 0;
		
		$segregation_score_1B = 0;
		$segregation_score_2B = 0;
		$segregation_score_3B = 0;
		$segregation_score_4B = 0;
		$segregation_score_5B = 0;
		$segregation_score_6B = 0;
		$segregation_score_7B = 0;
		
		$segregation_score_best_1 = 0;
		$segregation_score_best_2 = 0;
		$segregation_score_best_3 = 0;
		$segregation_score_best_4 = 0;
		$segregation_score_best_5 = 0;
		$segregation_score_best_6 = 0;
		$segregation_score_best_7 = 0;
		

		#######################################################################
		# Initially we assume this is not to be written to the filtered file  #
		# The criteria are:                                                   #
		#                                                                     #
		#  1.  If there are more than one ALT allele                          #
		#                                                                     #
		#  2.  If there is a snpEff or VEP effect score > 1 (i.e. any effect) #
		#                                                                     #
		#######################################################################

		$write_to_filtered = "false";


		#########################################################
        # Split whole line at TABs into an array myArray1 (1-9) #
        #########################################################
		
		@myArray1 = split (/\t/,$single_line);
		$array_size_1 = scalar @myArray1;
		$no_of_data_columns = $array_size_1 - 9;
		
	
		##################################################
        # Warn user if the number of samples seems wrong #
        ##################################################
		
		if ($no_of_data_columns != $no_of_samples)
		{
			&print_message("POSSIBLE FATAL ERROR on line $line_count of $vcf_file","warning");

			print "\tThe number of samples is unclear.\n\n";
			print "\tNumber of samples listed in the VCF file is $no_of_samples\n";
			print "\tNumber of samples in this data line is $no_of_data_columns.\n\n";
			
			print "Maybe the line is truncated?\n\n";
			print "$single_line\n\n";

			print "You should check the VCF file.\n\n";

			print "If you want to continue anyway, press 'return'   \n\n";

			$answer = <STDIN>;
		}
		
		
		#############################################
        # Get first fixed columns of the VCF file   #
        # CHR, POS, REF, ALT, QUAL, FILTER          #
        #############################################
		 
		$chromosome = $myArray1[0];
        $position = $myArray1[1];
        $REF_base = $myArray1[3];
        $ALT_base_overall = $myArray1[4];
        $QUAL = $myArray1[5];
        $FILTER = $myArray1[6];
		 
		if ($position =~ /^\d+?$/) {$vcf_variant_count = $vcf_variant_count + 1}
		
		# Convert chrX to chr39 for dog 
		if ($chromosome eq "chrX"){$chromosome = "chr39"}
         
		# Remove string 'chr'
        if (substr($chromosome,0,3) eq "chr"){$chromosome = substr($chromosome,3,99)}
		
		# If chromosome is an integer and less than 10 then add a leading zero
		if ($chromosome =~ /^\d+?$/){if ($chromosome < 10){$chromosome = "0".$chromosome}}




##     ## ######## ########  
##     ## ##       ##     ## 
##     ## ##       ##     ## 
##     ## ######   ########  
 ##   ##  ##       ##        
  ## ##   ##       ##        
   ###    ######## ##        

		####################################
		# If variant predictor is VEP then #
		####################################
		if ($effect_predictor eq "VEP")
		{
			if ($position eq "38939806xxx") {$show = "true"}


			#Set these to zero to start with
			$vep_results_string = "";
			$max_vep_effect_score = 0;
			$vep_effect_score = 0;
			$vep_effect = "";

			$vcf_id = $chromosome."_".$position;
			$vcf_id_count = 0;

			$position_plus_one = $position + 1;
			$vcf_id_plus_one = $chromosome."_".$position_plus_one;


			########################################################
			# If the vcf_id from the VCF file (chr_pos) is defined #
			########################################################
			$vep_id_found = "false";
			
			if ($show eq "true")
			{
				print "\n\nChecking VCF ID ($vcf_id) against vep_hash:\n\n";
				if (defined $vep_hash{$vcf_id}){print "vep_hash($vcf_id) is defined\n"} else {print "vep_hash($vcf_id) is NOT defined\n"}
			}
				
			if (defined $vep_hash{$vcf_id})
			{
				$vep_id_found = "true";
				
				$vep_position_matches_count = $vep_position_matches_count + 1;

				if ($show eq "true"){print "Position from VCF file is $position and was found in VEP files.\n\n";}
				
				$vep_effect_full_string = $vep_hash{$vcf_id};

				# Split effect string at commas (if there is more than one effect)
				@vep_effects_array = split (/,/,$vep_effect_full_string);

				$no_of_vep_effects = scalar @vep_effects_array;

				# Loop through the VEP effects to find the maximum effect #
				for ($effect_count = 1; $effect_count <= $no_of_vep_effects; $effect_count++)
				{
					$vep_effect_score = 0;
					$vep_effect = $vep_effects_array[$effect_count-1];

					###########################################
					# Get vep effect_score from subroutine    #
					###########################################
					&get_effect_score_VEP($vep_effect);

					# Remember maximum effect score for this position
					if ($vep_effect_score > $max_vep_effect_score){$max_vep_effect_score = $vep_effect_score}

					############################################################
					# Add vep result to string if not already in the string    #
					############################################################
					if (index($vep_results_string,$vep_effect)  == -1)
					{
						$vep_results_string = $vep_results_string."+".$vep_effect;
					}

				} # VEP effects loop

			} # if vep_hash is defined at this vcf_id
			

			############################################################
			# If the vcf_id from the VCF file (chr_pos) is NOT defined #
			# first look to see if it found at position minus one      #
			############################################################
			if (not defined $vep_hash{$vcf_id}) 
			{
					
				######################################################
				# If vep_id is not found for this VCF position, try  #
				# looking at one position less in vep_hash_minus_one #
				# (VEP can add one to the position)                  #
				######################################################
				if (defined $vep_hash_minus_one{$vcf_id})
				{
					$vep_effect_full_string = $vep_hash_minus_one{$vcf_id};
					$vep_id_found = "true";
					$vep_position_minus_one_count = $vep_position_minus_one_count + 1;

					# Split effect string at commas (if there is more than one effect)
					@vep_effects_array = split (/,/,$vep_effect_full_string);

					$no_of_vep_effects = scalar @vep_effects_array;

					# Loop through the VEP effects to find the maximum effect #
					for ($effect_count = 1; $effect_count <= $no_of_vep_effects; $effect_count++)
					{
						$vep_effect_score = 0;
						$vep_effect = $vep_effects_array[$effect_count-1];

						###########################################
						# Get vep effect_score from subroutine    #
						###########################################
						&get_effect_score_VEP($vep_effect);

						# Remember maximum effect score for this position
						if ($vep_effect_score > $max_vep_effect_score){$max_vep_effect_score = $vep_effect_score}

						############################################################
						# Add vep result to string if not already in the string    #
						############################################################
						if (index($vep_results_string,$vep_effect)  == -1)
						{
							$vep_results_string = $vep_results_string."+".$vep_effect;
						}

					} # VEP effects loop (minus_one)

				} # if (defined $vep_hash_minus_one{$vcf_id})


				###############################################################
				# If not defined vep_hash_minus_one (and vep_hash from above) #
				# then mark the effect as NOT FOUND so we can see it clearly  #
				###############################################################
				if (not defined $vep_hash_minus_one{$vcf_id})
				{
					$vep_effect_full_string = "NOT FOUND in VCF";
					$vep_results_string = "NOT FOUND in VCF.";
				}
					
			} # if not defined vep_hash(vcf_id)


			##################################################################
			# This section deals with the vep_hashes which were not defined  #
			# In other words there is no match for chr_pos from the VCF file #
			# This may be because it is a SNP and Indel at the same position #
			# so the second one will have had "_1" added to it,              #
			##################################################################

			do   # until (not defined $vep_hash{$vcf_id});
			{
				$vcf_id_count = $vcf_id_count + 1;
				$vcf_id = $vcf_id."_".$vcf_id_count;

				# Look at _1s where positions match exactly
				if (defined $vep_hash{$vcf_id})
				{
					$vep_id_found = "true";

					# If new effect (on second line) is not already in the string, add it
					if (index($vep_results_string,$vep_hash{$vcf_id}) == -1)
					{
						$vep_effect_full_string = $vep_hash{$vcf_id};

						
						######################################################
						# Count the number of times we have to use           #
						# the _1 system for two positions which are the same #
						######################################################
						if ($vcf_id_count == 1) {$vcf_id_1_count = $vcf_id_1_count + 1;}
						if ($vcf_id_count == 2) {$vcf_id_2_count = $vcf_id_2_count + 1;}

						# Split effect string at commas (if there is more than one effect)
						@vep_effects_array = split (/,/,$vep_effect_full_string);
						$no_of_vep_effects = scalar @vep_effects_array;

						# Loop through the VEP effects again
						for ($effect_count = 1; $effect_count <= $no_of_vep_effects; $effect_count++)
						{
							$vep_effect_score = 0;

							$vep_effect = $vep_effects_array[$effect_count-1];

							###########################################
							# Get vep effect_score from subroutine    #
							###########################################
							&get_effect_score_VEP($vep_effect);


							# Remember maximum effect score for this position
							if ($vep_effect_score > $max_vep_effect_score){$max_vep_effect_score = $vep_effect_score}

							############################################################
							# Add vep result to string if not already in the string    #
							############################################################
							if (index($vep_results_string,$vep_effect)  == -1)
							{
								$vep_results_string = $vep_results_string."+".$vep_effect;
							}

						} # VEP effects loop 3

					} # add effect to string

				} # if if (defined $vep_hash{$vcf_id}) when you add _1

				# Look at _1s where positions might match with minus_one (because of VEP)
				if (not defined $vep_hash{$vcf_id})
				{
					if (defined $vep_hash_minus_one{$vcf_id})
					{
						$vep_effect_full_string = $vep_hash_minus_one{$vcf_id};

						if ($vcf_id_count == 1) {$vcf_id_1_count = $vcf_id_1_count + 1;}
						if ($vcf_id_count == 2) {$vcf_id_2_count = $vcf_id_2_count + 1;}

						# Split effect string at commas (if there is more than one effect)
						@vep_effects_array = split (/,/,$vep_effect_full_string);
						$no_of_vep_effects = scalar @vep_effects_array;

						# Loop through the VEP effects again
						for ($effect_count = 1; $effect_count <= $no_of_vep_effects; $effect_count++)
						{
							$vep_effect_score = 0;

							$vep_effect = $vep_effects_array[$effect_count-1];

							###########################################
							# Get vep effect_score from subroutine    #
							###########################################
							&get_effect_score_VEP($vep_effect);


							# Remember maximum effect score for this position
							if ($vep_effect_score > $max_vep_effect_score){$max_vep_effect_score = $vep_effect_score}

							############################################################
							# Add vep result to string if not already in the string    #
							############################################################
							if (index($vep_results_string,$vep_effect)  == -1)
							{
								$vep_results_string = $vep_results_string."+".$vep_effect;
							}

						} # VEP effects loop 4
					}
				}

			}
			until (not defined $vep_hash{$vcf_id});


			######################################################################################
			# Are there any times we need to use the 1856675_1 system with the minus_one system? #
			# carrots
			######################################################################################


			###################################################
			# Count those VCF positions found in the VEP hash #
			###################################################
			
			if ($vep_id_found eq "true"){$vep_id_found_count = $vep_id_found_count + 1} else {$vep_id_not_found_count = $vep_id_not_found_count + 1}
			
			
			# Remove leading plus sign
			if (index($vep_results_string,"+") == 0)
			{
				$vep_results_string = substr($vep_results_string,1,999);
			} 

			#99999999
			if ($position eq "99999999")
			{
				print "\tCHECK 4 - final listing\n\n";
				print "\tPosition:              \t$position\n";
				print "\tmax_vep_effect_score:  \t$max_vep_effect_score\n";
				print "\tVEP_results_string:    \t$vep_results_string\n\n\n";

				$answer=<STDIN>;
			}


		} # effect predictor is VEP


		#################################################
        # Split the 7th item in the array at semicolons #
        # The 7th item is INFO (various):               #
        #                                               #
        # It also contains the snpEff results if there! #  <-- snpEff Results
        #                                               #
        #  AC=2;AF=1.00;AN=2;DP=27;Dels=0.00;FS=0.000;  #
        #  HRun=0;HaplotypeScore=0.6651;MQ=60.00;MQ0=0; #
        #  QD=37.00;SB=-376.46                          #
        #################################################
         
		
		#Effect ( Effect_Impact | Codon_Change | Amino_Acid_change | Gene_Name | Gene_BioType | Coding | Transcript | Rank [ | ERRORS | WARNINGS ] ) 


		 ######  ##    ## ########  ######## ######## ######## 
		##    ## ###   ## ##     ## ##       ##       ##       
		##       ####  ## ##     ## ##       ##       ##       
		 ######  ## ## ## ########  ######   ######   ######   
		      ## ##  #### ##        ##       ##       ##       
		##    ## ##   ### ##        ##       ##       ##       
		 ######  ##    ## ##        ######## ##       ##      


		##############################
		# If the 7th item is present #
		##############################
        if ($array_size_1 > 6)
		{ 
				@myArray7 = split(/;/, $myArray1[7]);
				$array_size_7 = scalar @myArray7;

				for ($array_count = 0; $array_count < $array_size_7; $array_count++)
				{

					#########################################################
					# snpEff records the effect in a string starting "EFF=" #
					#########################################################

					#EFF=NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|aCg/aTg|T80M|404|BEND4|protein_coding|CODING|ENSCAFT00000002832|2|1)	
					#EFF=DOWNSTREAM(MODIFIER||||970|ATP8A1|protein_coding|CODING|ENSCAFT00000043202||1),SYNONYMOUS_CODING(LOW|SILENT|ctG/ctA|L47|146|SHISA3|protein_coding|CODING|ENSCAFT00000002835|1|1|WARNING_TRANSCRIPT_NO_START_CODON)	

					if (index ($myArray7[$array_count],"EFF=") > -1)
					{
						#################################################
						# Record that this VCF fle contains snpEff data #
						#################################################
						$snpEff_data_found = "true";


						#Set these to zero to start with
						$snpEff_results_string 		= "";
						$snpEff_gene_string 		= "";
						$max_snpEff_effect_score 	= 0;
						$snpEff_effect_score 		= 0;
						$snpEff_effect 				= "";

						$snpEff_effect_full_string = $myArray7[$array_count];

						# Remove EFF= from the start using regex
						$snpEff_effect_full_string =~ s/EFF=//g;

						# Split effect string at commas (if there is more than one effect)
						@snpEff_effects_array = split (/,/,$snpEff_effect_full_string);

						$no_of_snpEff_effects = scalar @snpEff_effects_array;

						# Loop through the snpEff effects
						for ($effect_count = 1; $effect_count <= $no_of_snpEff_effects; $effect_count++)
						{
							$snpEff_effect_score = 0;

							$snpEff_effect_string = $snpEff_effects_array[$effect_count-1];

							############################################################################
							# The snpEff info has the Effect first, then the extra data in brackets () #
							# NON_SYNONYMOUS_CODING(MODERATE|MISSENSE|.... etc )                       #
							############################################################################
							$open_bracket_pos = index($snpEff_effect_string,"(");
							if ($open_bracket_pos > -1)
							{
								$snpEff_effect = substr($snpEff_effect_string,0,$open_bracket_pos);
								$rest_of_string = substr($snpEff_effect_string, $open_bracket_pos + 1,length($snpEff_effect_string) - $open_bracket_pos - 2);
							}
							
							#############################
							# Split the extra data at | #
							#############################
							@snpEff_items_array = split (/\|/,$rest_of_string);

							$snpEff_effect_impact = $snpEff_items_array[0];
							$snpEff_effect_codon_change = $snpEff_items_array[1];
							$snpEff_effect_gene_name = $snpEff_items_array[5];
							$snpEff_effect_transcript = $snpEff_items_array[8];


							#####################################################################
							# Set the snpEff effect score according to the impact to start with #
							# (This is overruled later on so maybe omit)                        #
							# Perhaps use as a safety net after assigning the Effect scaore. <<< carrots
							#####################################################################
							if ($snpEff_effect_impact eq "HIGH") {$snpEff_effect_score = 5}
							if ($snpEff_effect_impact eq "MODERATE") {$snpEff_effect_score = 4}
							if ($snpEff_effect_impact eq "LOW") {$snpEff_effect_score = 3}
							if ($snpEff_effect_impact eq "MODIFIER") {$snpEff_effect_score = 1}


							###########################################
							# Get snpEff effect_score from subroutine #
							###########################################

							&get_effect_score_snpEff($snpEff_effect);


							# Remember maximum effect score for this position #
							if ($snpEff_effect_score > $max_snpEff_effect_score){$max_snpEff_effect_score = $snpEff_effect_score}



							############################################################
							# Add snpEff result to string if not already in the string #
							############################################################
							if (index($snpEff_results_string,$snpEff_effect)  == -1)
							{
								$snpEff_results_string = $snpEff_results_string."+".$snpEff_effect;
							}


							############################################################
							# Add snpEff result to string if not already in the string #
							############################################################
							if (index($snpEff_gene_string,$snpEff_effect_gene_name)  == -1)
							{
								$snpEff_gene_string = $snpEff_gene_string."+".$snpEff_effect_gene_name;
							}



							if ($snpEff_effect_score > 22)
							{
								print "\tsnpEff CHECK 3   In loop. Effect count: $effect_count\n\n";
								print "\tPosition:                 \t$position\n";
								print "\tVEP Effect:               \t$snpEff_effect\n";
								print "\tsnpEff_effect_score:      \t$snpEff_effect_score\n";
								print "\tmax_snpEff_effect_score:  \t$max_snpEff_effect_score\n";
								print "\tsnpEff_results_string:    \t$snpEff_results_string\n";
								print "\tsnpEff_effect_gene_name:  \t$snpEff_effect_gene_name\n";
								print "\tsnpEff_gene_string:       \t$snpEff_gene_string\n";

								$answer=<STDIN>;
							}

						} # snpEff effects loop

						# Remove leading plus sign
						if (index($snpEff_results_string,"+") == 0)
						{
							$snpEff_results_string = substr($snpEff_results_string,1,999);
						} 
						# Remove leading plus sign
						if (index($snpEff_gene_string,"+") == 0)
						{
							$snpEff_gene_string = substr($snpEff_gene_string,1,999);
						} 


						#99999999
						if ($position eq "99999999")
						{
							print "\tsnpEff CHECK 4 - final listing\n\n";
							print "\tPosition:                 \t$position\n";
							print "\tmax_snpEff_effect_score:  \t$max_snpEff_effect_score\n";
							print "\tsnpEff_results_string:    \t$snpEff_results_string\n\n\n";

							$answer=<STDIN>;
						}

					} # snpEff data found in VCF file
				}

		} # if 7th item is present


		###############################################################
		# Get highest effect score from VEP and snpEff (if both used) #
		###############################################################
		if ($max_vep_effect_score > $max_snpEff_effect_score)
		{
			$max_effect_score = $max_vep_effect_score;
		}
		else
		{
			$max_effect_score = $max_snpEff_effect_score;
		}


		#####################################################################################
		# Parse myArray1(9)                                                                 #
        # These are the genotype columns (one for each sample - so there should be only 1)  # 
        # Split the 9th item in the array at COLONS                                         #
        #                                                                                   #
        # The 8th column should be GT:AD:DP:GQ:PL but it could be in any order so this      #
		# filed is used to determine the order of the data in the 9th item                  #
        #                                                                                   #
        # Note DP is twice in the row, at the start and in this genotype section.           #
        # They are not necessarily the same.  We use the genotype one.                      #
        #####################################################################################
        
         
		##############################
		# If the 9th item is present #
		##############################
        if ($array_size_1 > 8)
		{  
		 
            ################################################
            # Read the FORMAT in myArray1(8)               #
            #                                              #
            # This looks like this:  GT:AD:DP:GQ:PL        #
            # The genotype data is then read in this order #
            ################################################
            @myArray8 = split(":", $myArray1[8]);
			
            $array_size_8 = scalar @myArray8;
            
            ##########################################################
            # Get the genotype into a string                         #
            ##########################################################
            $genotype_string = $myArray1[9];

		
			################################################
			# Split ALT_base_overall at commas             #
			# to get all the different alleles in an array #
			################################################
			@ALT_allele_array = split(",",$ALT_base_overall);
			$no_of_alt_alleles = (scalar @ALT_allele_array);
			

			#####################################################
			# Make last item in ALT_allele_array the REF_base   #
			# This can be accessed if you want all the possible #
			# alleles, but looping with <= $no_of_alt_alleles   #
			# rather than < $no_of_alt_alleles                	#
			#####################################################
			$ALT_allele_array[$no_of_alt_alleles] = $REF_base;


			#############################################################################################
			# Make hash table for indels with REF_base as 'a' and remaining ALT alleles as 'b', 'c' etc #
			#############################################################################################
			$indel_hash{$REF_base} = "a";

			for ($array_count = 0;$array_count < $no_of_alt_alleles; $array_count++)
			{
				$indel_hash{$ALT_allele_array[$array_count]} = chr($array_count + 98);
			}

			
			$GT_string = "";
			$AD_string = "";
			$DP_genotype = "";
			$GQ_string = "";
			$PL_string = "";
		 	
				
			############################################################
            # Split myArray9 (genotype data string) using colons       #
			# Remember if genotype is ./. then there may be no colons  #
			# but splitting should still work if GT is there           #
			############################################################
           				
            @myArray9 = split(":",$genotype_string);
            $array_size_9 = scalar @myArray9;
            
		
            ################################################################
            # Now use order of fields in 'FORMAT' (myArray8)               #
            # to parse the genotype data fields of                         # 
            # each sample (order may not be the same for every VCF file)   #
            ################################################################
            for($array_count = 0; $array_count < $array_size_8; $array_count++)
			{
				$myString = $myArray8[$array_count];
                
				if ($myString eq "GT"){$GT_string = $myArray9[$array_count]}
				if ($myString eq "AD"){$AD_string = $myArray9[$array_count]}
				if ($myString eq "DP"){$DP_genotype = $myArray9[$array_count]}
				if ($myString eq "GQ"){$GQ_string = $myArray9[$array_count]}
				if ($myString eq "PL"){$PL_string = $myArray9[$array_count]}

            } # next array_count
            


			##############################################################################
			# Decide whether the variant type is a SNP or an Indel                       # 
			# If any variant is a different length than the REF base then it is an indel #
			##############################################################################
			$variant_type = "snp";
			for ($allele_count = 0; $allele_count <$no_of_alt_alleles; $allele_count++)
			{
				if (length($REF_base) != length($ALT_allele_array[$allele_count]))
				{
					$variant_type = "indel";
				}
			}

			 ######   #######  ##    ## ##     ## ######## ########  ########    ########  #######        ###    ########  
			##    ## ##     ## ###   ## ##     ## ##       ##     ##    ##          ##    ##     ##      ## ##   ##     ## 
			##       ##     ## ####  ## ##     ## ##       ##     ##    ##          ##    ##     ##     ##   ##  ##     ## 
			##       ##     ## ## ## ## ##     ## ######   ########     ##          ##    ##     ##    ##     ## ########  
			##       ##     ## ##  ####  ##   ##  ##       ##   ##      ##          ##    ##     ##    ######### ##     ## 
			##    ## ##     ## ##   ###   ## ##   ##       ##    ##     ##          ##    ##     ##    ##     ## ##     ## 
			 ######   #######  ##    ##    ###    ######## ##     ##    ##          ##     #######     ##     ## ########  


			#########################################################
			# Now convert the GT string from 0/1 format to AB format #
			#                                                       #
			#   - allele_number is the number from the string '0/1' #
			#   - allele_letter is the assigned alleles A,B,C,D etc #
			#   - allele_base is the actual base (for SNPs)         #
			#########################################################
            $slash_pos = index($GT_string, "/");
            
			$allele_number_1 = "";
			$allele_number_2 = "";

			$allele_base_1 = "";
			$allele_base_2 = "";
				
			# If there is a slash in the genotype e.g. 0/1, 2/2 etc		
            if ($slash_pos > 0)
			{
				$allele_number_1 = substr($GT_string,$slash_pos-1,1);  # The allele number comes from 0/1
                $allele_number_2 = substr($GT_string,$slash_pos+1,1);  
                

				#####################################################
				# We now want to convert this into actual bases     #
				# (allele_base_1 and allele_base_2)                 #
				# First we need to know how many alleles there are  #
				# by looking at ALT_base_overall.  This has all the #
				# ALT alleles found separated by commas. e.g. A,G,T #
				#                                                   #
				# Remember the 0 in 0/1 means the REF base          #
				#####################################################
				
				##############################################################
				# If genotype is NOT ./. (i.e. there is a valid genotype)    #
				##############################################################
				if ($allele_number_1 ne ".")
                {
                	if ($allele_number_1 == 0){$allele_base_1 = $REF_base}
					elsif ($allele_number_1 > 0){$allele_base_1 = $ALT_allele_array[$allele_number_1 - 1]}
					
					if ($allele_number_2 == 0){$allele_base_2 = $REF_base}
					elsif ($allele_number_2 > 0){$allele_base_2 = $ALT_allele_array[$allele_number_2 - 1]}

				} # if not ./.


				#############################################################
				# If genotype IS ./.  (this means no genotype was called)   #
				#############################################################
				if ($allele_number_1 eq ".")
                {

					if ($mark_missing_as_X eq "true")
					{
						$allele_base_1 = "X"; $allele_base_2 = "X";
					}
					else
					{
						$allele_base_1 = $REF_base;
						$allele_base_2 = $REF_base;
					}	

					#####################################
					# Count how many of these there are #
					# in total and for each sample 	    #
					#####################################
					$missing_genotypes_total = $missing_genotypes_total + 1;
					$missing_genotypes_count = $missing_genotypes_count + 1;

            	} # ./.

            } 
			#End If 'If slash_pos > 0  i.e. if the genotype is in the form 0/1 (with a forward slash)
				
			

			if ($show eq "true")
			{
				print ("At 1693\n");
				print "$single_line\n";
				print "Genotype_string:\t>$genotype_string<\n";
				print "Pos:\t$position\n";
				print "GT_string:       \t$GT_string\n";
				print "allele_number:   \t$allele_number_1\t$allele_number_2\n";
				print "Genotype:        \t$allele_base_1\t$allele_base_2\n\n";
				print "\n>\n";
				$answer = <STDIN>;
			}
			#show
			
			

			#######################################################################
			# If score_X_as_REF then treat Xs as REF_base for the sake of scoring #
			# Do this by putting the REF_base into $base_array 
			#
			# (the actual bases  #
			# have been written to the output file already so will appear as 'X'  #
			#######################################################################
			$col_1 = 1;
			$col_2 = 2;

			# First store original base in case it is X
			$base_orig_array[$col_1] = $allele_base_1;
			$base_orig_array[$col_2] = $allele_base_2;

			# Then replace it by X if required
			if ($score_X_as_REF eq "true")
			{
				if ($allele_base_1 eq "X") {$allele_base_1 = $REF_base}
				if ($allele_base_2 eq "X") {$allele_base_2 = $REF_base}
			}

			###############################################
			# Store the two bases in an array of bases    #
			# for later use in statistics                 #
			# (Note these are the actual bases, ACGT etc) #
			###############################################
			$base_array[$col_1] = $allele_base_1;
			$base_array[$col_2] = $allele_base_2;


			# Update some counters of the genotypes counted
			$genotypes_counted_total = $genotypes_counted_total + 1;
			$genotypes_counted = $genotypes_counted + 1;
			

			
			###########################################################################
			#                                Statistics                               #
			###########################################################################

			#########################################################################
			# New way to get main affected allele (works for any number of alleles) #
			# This doesn't need allele_A and allele_B  .                            #
			#                                                                       #
			# When analysing for segregation score:                                 #
			#    Main affected allele is called 'A'.  All the rest are called 'B'   #
			#########################################################################
			
			

			 ######  ##     ##  #######   #######   ######  ########    ########     ###     ######  ########  ######  
			##    ## ##     ## ##     ## ##     ## ##    ## ##          ##     ##   ## ##   ##    ## ##       ##    ## 
			##       ##     ## ##     ## ##     ## ##       ##          ##     ##  ##   ##  ##       ##       ##       
			##       ######### ##     ## ##     ##  ######  ######      ########  ##     ##  ######  ######    ######  
			##       ##     ## ##     ## ##     ##       ## ##          ##     ## #########       ## ##             ## 
			##    ## ##     ## ##     ## ##     ## ##    ## ##          ##     ## ##     ## ##    ## ##       ##    ## 
			 ######  ##     ##  #######   #######   ######  ########    ########  ##     ##  ######  ########  ######  
 			

 			########################################################################################
 			#  Choose the alleles for all the samples and save them in allele_for_writing_array    #
 			#                                                                                      #
 			# For the SNPs store the actual base (ACGT) and for the indels use a code a,b,c,d etc  #
 			########################################################################################

			$col_1 = 1;
			$col_2 = 2;

			if ($variant_type eq "snp")
				{	
				$allele_base_1 = $base_array[$col_1];
				$allele_base_2 = $base_array[$col_2];
			}

			if ($variant_type eq "indel")
				{
				$allele_base_1 = $indel_hash{$base_array[$col_1]};
				$allele_base_2 = $indel_hash{$base_array[$col_2]};
			}


			#########################################################
			# Put alleles in alphabetical order for Excel           #
			# (just because this is the way NGS SNP Handler did it) #
			#########################################################

			if (ord $allele_base_1 > ord $allele_base_2)
			{
				$myString =$allele_base_1;
				$allele_base_1 = $allele_base_2;
				$allele_base_2 = $myString;
			}


			##########################################
			# Save alleles for writing to file later #
			##########################################
			$allele_1_for_writing = $allele_base_1;
			$allele_2_for_writing = $allele_base_2;



			################################################
			# If the original base was X save this instead #
			################################################
			if ($base_orig_array[$col_1] eq "X")
			{
				$allele_1_for_writing = $base_orig_array[$col_1];
				$allele_2_for_writing = $base_orig_array[$col_1];
			}


		

			##      ## ########  #### ######## ########     #######  ##     ## ######## ########  ##     ## ######## 
			##  ##  ## ##     ##  ##     ##    ##          ##     ## ##     ##    ##    ##     ## ##     ##    ##    
			##  ##  ## ##     ##  ##     ##    ##          ##     ## ##     ##    ##    ##     ## ##     ##    ##    
			##  ##  ## ########   ##     ##    ######      ##     ## ##     ##    ##    ########  ##     ##    ##    
			##  ##  ## ##   ##    ##     ##    ##          ##     ## ##     ##    ##    ##        ##     ##    ##    
			##  ##  ## ##    ##   ##     ##    ##          ##     ## ##     ##    ##    ##        ##     ##    ##    
			 ###  ###  ##     ## ####    ##    ########     #######   #######     ##    ##         #######     ##  


			######################################################### 
			# Now write the first three chr, pos and ref columns    # 
			# followed by the columns with the actual alleles       #  
			#########################################################   
			
			print OUT "$chromosome\t$position\t$REF_base\t$ALT_base_overall";

			print OUT "\t$allele_1_for_writing";
			print OUT "\t$allele_2_for_writing";
			

			#####################################################################
			# Choose the first ALT allele as the "main affected allele" for now #  TEMP!!
			#####################################################################
			$main_affected_allele = $ALT_allele_array[0]; 
			
			
			print OUT "\t$max_effect_score"; # Effect score
			
			if ($effect_predictor eq "snpEff") {print OUT "\t$snpEff_results_string"} # Consequence
			if ($effect_predictor eq "VEP") {print OUT "\t$vep_results_string"} # Consequence
			
			if ($show eq "true")
			{
				print "Variant type: $variant_type\n";
				print "Main affected allele: $main_affected_allele\n";
				print "indel_hash{$main_affected_allele}: $indel_hash{$main_affected_allele}\n";
				print "main_affected_allele: $main_affected_allele\n";

				&pause;
			}

			######################################################### 
			# Now write to the OUT_INDELS file, but only if the     # 
			# variant is an Indel                                   #
			#########################################################  
			
			if ($variant_type eq "indel")
			{
				print INDELS_OUTPUT "$chromosome\t$position\t$REF_base\t$ALT_base_overall";

				print INDELS_OUTPUT "\t$allele_1_for_writing";
				print INDELS_OUTPUT "\t$allele_2_for_writing";
				
				
				if ($main_affected_allele ne $REF_base){print INDELS_OUTPUT "\tN"} else {print INDELS_OUTPUT "\tR";}
				
				print INDELS_OUTPUT "\t$max_effect_score"; # Effect score
				
				if ($effect_predictor eq "snpEff") {print INDELS_OUTPUT "\t$snpEff_results_string"} # Consequence
				if ($effect_predictor eq "VEP") {print INDELS_OUTPUT "\t$vep_results_string"} # Consequence

				

				# No option for extra columns 
			
				print INDELS_OUTPUT "\t$no_of_alt_alleles";

				print INDELS_OUTPUT "\n";

			} # If variant is Indel

			
			print OUT "\t$no_of_alt_alleles";

			print OUT "\t$snpEff_gene_string";
			

			#############################
			# Now the final end of line #
			#############################
			print OUT "\n";


			######################################################### 
			# Now write to the OUT_FILTERED file, but only if the   # 
			# criteria have been met                                #
			#########################################################   

			# Number of ALT alleles criterion
			if ($no_of_alt_alleles > 1)
			{
				$write_to_filtered = "true";
				$filter_count_no_of_alleles = $filter_count_no_of_alleles + 1;
			}

			# Varaiant effect criterion
			if ($max_effect_score > 1)
			{
				$write_to_filtered = "true";
				$filter_count_effect_score = $filter_count_effect_score + 1;
			}
			
			if ($write_to_filtered eq "true")
			{
				print OUT_FILTERED "$chromosome\t$position\t$REF_base\t$ALT_base_overall";

				print OUT_FILTERED "\t$allele_1_for_writing";
				print OUT_FILTERED "\t$allele_2_for_writing";
				
				print OUT_FILTERED "\t$max_effect_score"; # Effect score
				
				if ($effect_predictor eq "snpEff") {print OUT_FILTERED "\t$snpEff_results_string"} # Consequence
				if ($effect_predictor eq "VEP") {print OUT_FILTERED "\t$vep_results_string"} # Consequence
				
				print OUT_FILTERED "\t$no_of_alt_alleles";
				print OUT_FILTERED "\t$snpEff_gene_string";
				
			
				#############################
				# Now the final end of line #
				#############################
				print OUT_FILTERED "\n";
				
			} # if write_to_filtered eq "true"

        } #End If 'If Ubound > 9
 
	}
	# if start_storing = true


	
  ## ##    ######  ##     ## ########   #######  ##     ##    ##       #### ##    ## ######## 
  ## ##   ##    ## ##     ## ##     ## ##     ## ###   ###    ##        ##  ###   ## ##       
######### ##       ##     ## ##     ## ##     ## #### ####    ##        ##  ####  ## ##       
  ## ##   ##       ######### ########  ##     ## ## ### ##    ##        ##  ## ## ## ######   
######### ##       ##     ## ##   ##   ##     ## ##     ##    ##        ##  ##  #### ##       
  ## ##   ##    ## ##     ## ##    ##  ##     ## ##     ##    ##        ##  ##   ### ##       
  ## ##    ######  ##     ## ##     ##  #######  ##     ##    ######## #### ##    ## ######## 


	#########################################################
    # Has the start of the CHR data been reached yet?       #
	# The line before the actual data starts with '#CHROM'  #
	# This part of the code checks that line, which has the #
	# (very important) headers for the sample columns.      #
    #########################################################
	
	if (index($single_line,"#CHROM") > -1)
	{
		$start_storing = "true";
		
		############################################################
		# Parsing of #CHROM data line to get a list of input files #
		############################################################

		&print_message("Parsing of #CHROM data line to get list of input sample names","message");

		@chrom_line_array = split(/\s+/,$single_line);
		$chrom_line_array_size = (scalar @chrom_line_array) - 1;
		$no_of_samples_chrom_line = $chrom_line_array_size - 8; # ignore first columns 0-9
		
		for ($array_count = 9; $array_count <= $chrom_line_array_size; $array_count++)
		{
			# Sample names go into sample_name_array
			$sample_name_array[$array_count - 8] = $chrom_line_array[$array_count];
			$file_count = $array_count - 8;
			print "$chrom_line_array[$array_count]\n";
		}
	

         
		$no_of_samples = $no_of_samples_chrom_line;



		if ($no_of_samples > 1)
		{
			&print_message("This program only works on single column VCF files","warning");

			print "\nThere are $no_of_samples samples in this multi-column VCF file.\n\n";
			print "\tYou need to run vcf2excel for this file\n\n";

			print "\tIf it really is all one sample then you may need to change the \n";
			print "\tread group headings so that it looks like data from a single sample.\n\n";
		
        	exit;
    	}
    	if ($no_of_samples == 1)
		{
			&print_message("VCF file confirmed as single sample files","message");

			print "\nThere is $no_of_samples sample in this multi-column VCF file.\n\n";
		
        	&pause;
    	}


		#################################
		# Show user a list of the files #
		#################################
		
		&print_message("Sample columns found in the VCF file (from the '#CHROM' line)","message");
		

		for ($sample_count = 1; $sample_count <= $no_of_samples; $sample_count++)
		{
			print "$sample_count \t$sample_name_array[$sample_count]\n";
		}
		
		
		
		$no_samples_in_status_file = 1;
		
		
		
		print "\nConverting VCF file to for_excel file...\n\n";

		##############################
		# Set various things to zero #
		##############################

		$missing_genotypes_count = 0;
		$genotypes_counted = 0;
		
        ######################################
        # Open the output files              #
        ######################################
        # Out for all variants (unfiltered)
		open (OUT, ">$output_file")|| die "Cannot create output file: $output_file";

		# Output for selected (filtered) variants
		open (OUT_FILTERED, ">$output_file_filtered")|| die "Cannot create output file: $output_file_filtered";
		
		# Output for indels
		open (INDELS_OUTPUT, ">$simplified_indels_output_file")|| die "Cannot create output file: $simplified_indels_output_file";


		############################################################
		# Write headers for the first 3 columns in the output file #
		############################################################
		print OUT "Chr\tPos\tRef\tAlt";
		print OUT_FILTERED "Chr\tPos\tRef\tAlt";
		print INDELS_OUTPUT "Chr\tPos\tRef\tAlt";
		

		########################################
        # Now add headers for the output file  #
        ########################################

        ##########################################################
        # Add name of sample as header in first batch of columns #
        # (ie those after the first three CHR,POS,REF columns    #
        # with the ALT base for each sample                      #
        ##########################################################
        $sample_count=1;

        $sample_name = $sample_name_array[1];
		
        print OUT "\t$sample_name_array[1]"; # Once for allele_1 
		print OUT "\t$sample_name_array[1]"; # Once for allele_2

		print INDELS_OUTPUT "\t$sample_name_array[1]"; # Once for allele_1 
		print INDELS_OUTPUT "\t$sample_name_array[1]"; # Once for allele_2

		print OUT_FILTERED "\t$sample_name_array[1]"; # Once for allele_1 
		print OUT_FILTERED "\t$sample_name_array[1]"; # Once for allele_2

	
		#Main file
		print OUT "\tEffect score";
		print OUT "\tConsequence";

		#Indels file
		print INDELS_OUTPUT "\tEffect score";
		print INDELS_OUTPUT "\tConsequence";

		#Filtered output
		print OUT_FILTERED "\tEffect score";
		print OUT_FILTERED "\tConsequence";

		# No of ALT alleles (useful for sorting)
		print OUT "\tNo of ALT alleles";
		print INDELS_OUTPUT "\tNo of ALT alleles";
		print OUT_FILTERED "\tNo of ALT alleles";

		print OUT "\tGene";
		print INDELS_OUTPUT "\tGene";
		print OUT_FILTERED "\tGene";

		
		#############################
		# Now the final end of line #
		#############################
		print OUT "\n";
		print INDELS_OUTPUT "\n";
		print OUT_FILTERED "\n";
		

		&print_message("Processing VCF file...","message");

	} # End of: if (index($single_line,"#CHROM") > -1)
	
	
	if (($line_count % 10000) == 0)
	{
		if ($start_storing eq "true"){print "Reading VCF file.  Line: $line_count\n"}
	}
	 
} # End of: while ($single_line = <VCF>) 

close OUT; # close out file for excel
close OUT_FILTERED; # close out file for excel
close INDELS_OUTPUT;


print "\n\n";


#############################
# Write details to log file #
#############################
print COMMAND_LOG "\n\n####################\n";
print COMMAND_LOG "# Analysis details #\n";
print COMMAND_LOG "####################\n\n";

&print_message("FINISHED ANALYSIS","message");

# Missing genotypes
&print_both("\nMissing genotypes:  \n\n");

	&print_both("\tSample\tNo missing\t% missing\n");
	&print_both("\t======\t==========\t=========\n");

	
	#Calculate precentage missing for THE sample
	$missing_genotypes_percent = $missing_genotypes_count/$genotypes_counted * 100;
	$missing_genotypes_percent = sprintf("%.3f", $missing_genotypes_percent);

	&print_both("\t$sample_name_array[1]\t$missing_genotypes_count   \t$missing_genotypes_percent%\n");
	

	# Some calculations
	$missing_genotypes_percent = $missing_genotypes_total/$genotypes_counted_total * 100;
	$missing_genotypes_percent = sprintf("%.3f", $missing_genotypes_percent);
	$genotypes_counted_per_sample = $genotypes_counted_total;
	$genotypes_counted_per_sample = sprintf("%.0f", $genotypes_counted_per_sample);

	&print_both("\n\tTotal number of genotypes counted:  \t$genotypes_counted_total\n");
	&print_both("\tTotal number of missing genotypes:  \t$missing_genotypes_total\n");
	&print_both("\tTotal percentage missing genotypes: \t$missing_genotypes_percent%\n");


&print_both ("\nVariant annotations: \n\n");

	&print_both ("\tVariant effect predictor selected: $effect_predictor\n");

	if (($VEP_data_found eq "false") && ($snpEff_data_found eq "false"))
	{
		&print_both ("\tNo variant annotations found\n");
	}

	if (($effect_predictor eq "snpEff") && ($snpEff_data_found eq "true"))
	{
		&print_both ("\tVariant effect predictor used:     $effect_predictor\n");
	}
	if (($effect_predictor eq "VEP") && ($VEP_data_found eq "true"))
	{
		&print_both ("\tVariant effect predictor used:     $effect_predictor\n");
	}

	if (($effect_predictor eq "snpEff") && ($snpEff_data_found eq "false"))
	{
		&print_message("WARNING. snpEff data not found","warning");

		&print_both ("\tYou said the file was annotated by snpEff but snpEff data was not found\n");
	}
	if (($effect_predictor eq "VEP") && ($snpEff_data_found eq "true") && ($VEP_data_found eq "true"))
	{
		&print_both ("\tYou said you wanted to use VEP annotations but snpEff data was found as well (VEP only was used)\n");
	}


# SNPs written to filtered file
&print_both("\n\nSNPs retained and written to filtered file:\n\n");

	&print_both("\tSNPs retained because of effect score:         \t$filter_count_effect_score\n");
	&print_both("\tSNPs retained because more than 2 alleles:     \t$filter_count_no_of_alleles\n");



&print_both("\n\nCheck on numbers found:\n\n");

    &print_both("\tNo. of variants in VCF file:                             \t$vcf_variant_count\n");

if ($effect_predictor eq "VEP") 
{
	&print_both("\tNo. of variants in VEP file:                             \t$vep_variant_count\n");

	&print_both("\tNo of VCF positions found in VEP:                        \t$vep_id_found_count\n");
	&print_both("\tNo of VCF positions not found in VEP:                    \t$vep_id_not_found_count   \t<<<< This should be zero!!\n");

	&print_both("\tNo of VCF positions directly matching VEP position:      \t$vep_position_matches_count\n");
	&print_both("\tNo of VCF positions matching at VEP position minus one:  \t$vep_position_minus_one_count\n");
}



&print_both("\n\nFiles used\n\n");

&print_both ("\tInput VCF file:           \t$vcf_file\n");
&print_both ("\tOutput file (filtered):   \t$output_file_filtered\n");
&print_both ("\tOutput file (full):       \t$output_file\n");
&print_both ("\tMerged Indels file:       \t$indels_merged_file\n");

&print_message("The output files can now be viewed in NGS SNP Viewer","message");

close COMMAND_LOG;

exit;

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


##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

foreach (@_) {s/\n//g}  
foreach (@_) {s/\r//g}  

}


sub list_genotypes
{
	print "\tAffecteds: ";
	for ($sample_test_count = 1; $sample_test_count <= $no_affected_samples; $sample_test_count++)
	{
		$col_1 =  1;
		$col_2 =  2;
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		 print "\t$allele_base_1$allele_base_2";
	}
	print "\n\tCarriers: ";
	for ($sample_test_count = $no_affected_samples + 1; $sample_test_count <= ($no_affected_samples + $no_carrier_samples); $sample_test_count++)
	{
		$col_1 = 1;
		$col_2 = 2;
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		print "\t$allele_base_1$allele_base_2";
	}
	print "\n\tNormals: ";
	for ($sample_test_count = $no_affected_samples + $no_carrier_samples + 1; $sample_test_count <= $no_of_samples; $sample_test_count++)
	{
		$col_1 = 1;
		$col_2 = 2;
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		print "\t$allele_base_1$allele_base_2";
	}
	print"\n\n";
}


sub list_genotypes_brief
{
	print "\nGenotypes affected: ";
	for ($sample_test_count = 1; $sample_test_count <= $no_affected_samples; $sample_test_count++)
	{
		$col_1 = 1;
		$col_2 = 2;
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		 print "\t$allele_base_1-$allele_base_2 ";
	}
	if ($no_carrier_samples > 0)
	{
		print "\nGenotypes carrier: ";
		for ($sample_test_count = $no_affected_samples + 1; $sample_test_count <= ($no_affected_samples + $no_carrier_samples); $sample_test_count++)
		{
			$col_1 = 1;
			$col_2 = 2;
			$allele_base_1 = $base_array[$col_1];
			$allele_base_2 = $base_array[$col_2];
			print "\t$allele_base_1-$allele_base_2 ";
		}
	}
	print "\nGenotypes normal: ";
	for ($sample_test_count = $no_affected_samples + $no_carrier_samples + 1; $sample_test_count <= $no_of_samples; $sample_test_count++)
	{
		$col_1 = 1;
		$col_2 = 2;
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		print "\t$allele_base_1-$allele_base_2 ";
	}
	print"\n\n";
}

sub list_genotypes_simple
{
	my $check_count = 0;
	
	print "\n\tAffecteds: ";
	for ($check_count = 1; $check_count <= $no_affected_samples; $check_count++)
	{
		$col_1 = 1;
		$col_2 = 2;
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		 print "\t$allele_base_1-$allele_base_2";
	}
	print "\n\tNormals: ";
	for ($check_count = $no_affected_samples + 1; $check_count <= $no_of_samples; $check_count++)
	{
		$col_1 = 1;
		$col_2 = 2;
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		print "\t$allele_base_1-$allele_base_2";
	}
}
sub list_genotypes_AB
{
	my $check_count = 0;
	
	print "\n\tAffecteds: ";
	for ($check_count = 1; $check_count <= $no_affected_samples; $check_count++)
	{
		$col_1 = 1;
		$col_2 = 2;
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];

		if    (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype = "AA"}
	 	elsif (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype = "AB"}
	 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype = "AB"}
	 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype = "BB"}

		 print "\t$genotype";
	}
	print "\n\tNormals: ";
	for ($check_count = $no_affected_samples + 1; $check_count <= $no_of_samples; $check_count++)
	{
		$col_1 = 1;
		$col_2 = 2;
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];

		if    (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype = "AA"}
	 	elsif (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype = "AB"}
	 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype = "AB"}
	 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype = "BB"}

		 print "\t$genotype";
	}
}





sub show_counters
{
			print "\nREF: $REF_base\n";
			
			print "\tAA_count_affected:    \t$AA_count_affected\tAA_count_normal:      \t$AA_count_normal\n";
			print "\tAB_count_affected:    \t$AB_count_affected\tAB_count_normal:      \t$AB_count_normal\n";
			print "\tBB_count_affected:    \t$BB_count_affected\tBB_count_normal:      \t$BB_count_normal\n";
			print "\tAA_AB_count_affected: \t$AA_AB_count_affected\tAA_AB_count_normal:   \t$AA_AB_count_normal\n";
			print "\tAB_BB_count_affected: \t$AB_BB_count_affected\tAB_BB_count_normal:   \t$AB_BB_count_normal\n";
			print "\tHH_count_affected:    \t$HH_count_affected\tHH_count_normal:      \t$HH_count_normal\n\n";
			
}

sub show_counters_allele
{
			print "\tREF: $REF_base\n";
			
			print "\tA_count_affected:    \t$A_count_affected\n";
			print "\tB_count_affected:    \t$B_count_affected\n";
			print "\tC_count_affected:    \t$C_count_affected\n";
			print "\tA_count_normal:      \t$A_count_normal\n";
			print "\tB_count_normal:      \t$B_count_normal\n";
			print "\tC_count_normal:      \t$C_count_normal\n";
			print "\tMain affected:     \t$main_affected_allele\n";
			print "\tMain normal:       \t$main_normal_allele\n";
			print "--------------------------------------------\n";

}

sub show_segregation_scores
{
# 1 AA_AB_BB         	Recessive:  								Affecteds AA, carriers AB, normals BB
# 2 AA_AB_ABorBB     	Recessive, some normals may be carriers:  	Affecteds AA, carriers AB, normals AB or BB
# 3 AA_ABorBB_BB      	Recessive, some carriers may be normals:  	Affecteds AA, carriers AB or BB, normals BB
# 4 AA_ABorBB_ABorBB 	Recessive:  								Affecteds AA, carriers AB or BB, normals AB or BB
# 5 AAorAB_BB        	Dominant:   								Affecteds AA or AB, normals BB
# 6 A_B              	Additive:   								Affecteds A with As, Normals with Bs
# 7 AA_BB               Strict:                                     Affecteds AA, Normals BB
	
	print "Main affected allele: $main_affected_allele\n\n";
	
	print "1: AA_AB_BB         \t$segregation_score_1\t$segregation_score_1B\t$segregation_score_best_1\n";
	print "2: AA_AB_ABorBB     \t$segregation_score_2\t$segregation_score_2B\t$segregation_score_best_2\n";
	print "3: AA_ABorBB_BB     \t$segregation_score_3\t$segregation_score_3B\t$segregation_score_best_3\n";
	print "4: AA_ABorBB_ABorBB \t$segregation_score_4\t$segregation_score_4B\t$segregation_score_best_4\n";
	print "5: AAorAB_BB        \t$segregation_score_5\t$segregation_score_5B\t$segregation_score_best_5\n";
	print "6: A_B              \t$segregation_score_6\t$segregation_score_6B\t$segregation_score_best_6\n";
	print "7: AA_BB            \t$segregation_score_7\t$segregation_score_7B\t$segregation_score_best_7\n";

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
	if ($style eq "help"){$char = "+"}
	
	print "\n\n";
	print color ' bold yellow';
	if ($style eq "warning"){print color ' bold red'}
	if ($style eq "input"){print color ' bold white'}
	if ($style eq "help"){print color ' bold green'}
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n$char    $message    $char\n";
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n\n";
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

sub get_effect_score_snpEff
{
	my $effect = $_[0];

#frameshift_variant+splice_region_variant+intron_variant+feature_elongation =5
#INTRON+SPLICE_SITE_ACCEPTOR+SPLICE_SITE_DONOR = 4

	#Effect score 5
	if ($effect eq "START_LOST"){$snpEff_effect_score = 5}
	if ($effect eq "EXON_DELETED"){$snpEff_effect_score = 5}
	if ($effect eq "FRAME_SHIFT"){$snpEff_effect_score = 5}
	if ($effect eq "STOP_GAINED"){$snpEff_effect_score = 5}
	if ($effect eq "STOP_LOST"){$snpEff_effect_score = 5}
	if ($effect eq "RARE_AMINO_ACID"){$snpEff_effect_score = 5}

	#Effect score promoted to 5 from snpEff 'Moderate'
	if ($effect eq "NON_SYNONYMOUS_CODING"){$snpEff_effect_score = 5}
	if ($effect eq "CODON_INSERTION"){$snpEff_effect_score = 5}
	if ($effect eq "CODON_CHANGE_PLUS_CODON_INSERTION"){$snpEff_effect_score = 5}
	if ($effect eq "CODON_DELETION"){$snpEff_effect_score = 5}
	if ($effect eq "CODON_CHANGE_PLUS_CODON_DELETION"){$snpEff_effect_score = 5}

	# Effect score 4 - might easily have an effect
	if ($effect eq "SPLICE_SITE_ACCEPTOR"){$snpEff_effect_score = 4}
	if ($effect eq "SPLICE_SITE_DONOR"){$snpEff_effect_score = 4}

	#Can't work out if the change is synonymous or non-synonymous
	if ($effect eq "CODON_CHANGE"){$snpEff_effect_score = 4}

	# Effect score 3 - could have some effect
	if ($effect eq "UTR_5_DELETED"){$snpEff_effect_score = 3}
	if ($effect eq "UTR_3_DELETED"){$snpEff_effect_score = 3}
	if ($effect eq "INTRON_CONSERVED"){$snpEff_effect_score = 3}
	if ($effect eq "INTERGENIC_CONSERVED"){$snpEff_effect_score = 3}
	if ($effect eq "UTR_3_PRIME"){$snpEff_effect_score = 3}
	if ($effect eq "UTR_5_PRIME"){$snpEff_effect_score = 3}
	if ($effect eq "TRANSCRIPT"){$snpEff_effect_score = 3}
	if ($effect eq "SPLICE_SITE_BRANCH"){$snpEff_effect_score = 3}
	if ($effect eq "SPLICE_SITE_BRANCH_U12"){$snpEff_effect_score = 3}
	if ($effect eq "EXON"){$snpEff_effect_score = 3} # These could be non-coding exons


	# Effect score 2 - upstream or downstream so could effect regulatory elements
	if ($effect eq "DOWNSTREAM"){$snpEff_effect_score = 2}
	if ($effect eq "UPSTREAM"){$snpEff_effect_score = 2}


	# Effect score 1 - unlikely to have an effect 
	if ($effect eq "INTRON"){$snpEff_effect_score = 1}
	if ($effect eq "SYNONYMOUS_CODING"){$snpEff_effect_score = 1}
	if ($effect eq "INTRAGENIC"){$snpEff_effect_score = 1}
	if ($effect eq "INTERGENIC"){$snpEff_effect_score = 1}
	if ($effect eq "CDS"){$snpEff_effect_score = 1}

	if ($snpEff_effect_score == 0){$snpEff_effect_score = 4.4} # Have this as default so none slip through

}

sub get_effect_score_VEP
{
	my $effect = $_[0];

#intron_variant+non_coding_exon_variant+nc_transcript_variant = 3
	#Effect score 5
	if ($effect eq "transcript_ablation"){$vep_effect_score = 5}
	if ($effect eq "stop_gained"){$vep_effect_score = 5}
	if ($effect eq "frameshift_variant"){$vep_effect_score = 5}
	if ($effect eq "stop_lost"){$vep_effect_score = 5}
	if ($effect eq "initiator_codon_variant"){$vep_effect_score = 5} # Same as start_lost?
	if ($effect eq "inframe_insertion"){$vep_effect_score = 5}
	if ($effect eq "inframe_deletion"){$vep_effect_score = 5}
	if ($effect eq "missense_variant"){$vep_effect_score = 5}

	# Effect score 4 - might easily have an effect
	if ($effect eq "splice_donor_variant"){$vep_effect_score = 4}
	if ($effect eq "splice_acceptor_variant"){$vep_effect_score = 4}
	if ($effect eq "incomplete_terminal_codon_variant"){$vep_effect_score = 4}

	#Can't work out if the change is synonymous or non-synonymous
	if ($effect eq "coding_sequence_variant"){$vep_effect_score = 4}
	

	# Effect score 3 - could have some effect
	if ($effect eq "splice_region_variant"){$vep_effect_score = 3}
	if ($effect eq "5_prime_UTR_variant"){$vep_effect_score = 3}
	if ($effect eq "3_prime_UTR_variant"){$vep_effect_score = 3}
	if ($effect eq "non_coding_exon_variant"){$vep_effect_score = 3}
	if ($effect eq "regulatory_region_ablation"){$vep_effect_score = 3}
	if ($effect eq "regulatory_region_amplification"){$vep_effect_score = 3}
	if ($effect eq "nc_transcript_variant"){$vep_effect_score = 3}

	# Effect score 2 - upstream or downstream so could effect regulatory elements
	if ($effect eq "upstream_gene_variant"){$vep_effect_score = 2} 
	if ($effect eq "downstream_gene_variant"){$vep_effect_score = 2}

	# Effect score 1 - unlikely to have an effect
	if ($effect eq "intron_variant"){$vep_effect_score = 1}
	if ($effect eq "feature_truncation"){$vep_effect_score = 1}
	if ($effect eq "feature_elongation"){$vep_effect_score = 1}
	if ($effect eq "intergenic_variant"){$vep_effect_score = 1}
	if ($effect eq "synonymous_variant"){$vep_effect_score = 1}
	if ($effect eq "stop_retained_variant"){$vep_effect_score = 1}
	if ($effect eq "feature_elongation"){$vep_effect_score = 1}

	if ($vep_effect_score == 0){$vep_effect_score = 4.5} # Have this as default so none slip through
}


sub choose_best_segregation_scores
{
	###########################################################
	# Choose best segregation score out of segregation_score  #
	# and segregation_score_B.  This removes the problem that #
	# the scores depends on how the "main affected allele" is #
	# chosen.                                                 #
	###########################################################
	$swapped = "false";
	
	if($segregation_score_1 >= $segregation_score_1B)
	{$segregation_score_best_1 = $segregation_score_1}
	else
	{$segregation_score_best_1 = $segregation_score_1B; $swapped = "true"}
	
	if($segregation_score_2 >= $segregation_score_2B)
	{$segregation_score_best_2 = $segregation_score_2}
	else
	{$segregation_score_best_2 = $segregation_score_2B; $swapped = "true"}
	
	if($segregation_score_3 >= $segregation_score_3B)
	{$segregation_score_best_3 = $segregation_score_3}
	else
	{$segregation_score_best_3 = $segregation_score_3B; $swapped = "true"}
	
	if($segregation_score_4 >= $segregation_score_4B)
	{$segregation_score_best_4 = $segregation_score_4}
	else
	{$segregation_score_best_4 = $segregation_score_4B; $swapped = "true"}
	
	if($segregation_score_5 >= $segregation_score_5B)
	{$segregation_score_best_5 = $segregation_score_5}
	else
	{$segregation_score_best_5 = $segregation_score_5B; $swapped = "true"}
	
	if($segregation_score_6 >= $segregation_score_6B)
	{$segregation_score_best_6 = $segregation_score_6}
	else
	{$segregation_score_best_6 = $segregation_score_6B; $swapped = "true"}
	 
		if($segregation_score_7 >= $segregation_score_7B)
	{$segregation_score_best_7 = $segregation_score_7}
	else
	{$segregation_score_best_7 = $segregation_score_7B; $swapped = "true"}


	############################################################################################
	# Normalise additive score (segregation_score_best_6) so it has the same maximum as others #
	############################################################################################

	$segregation_score_best_6 = ($segregation_score_best_6 * $no_of_samples)/ ($no_affected_samples * 2 + $no_normal_samples * 2 + $no_carrier_samples);

	$segregation_score_best_6 = sprintf("%.1f", $segregation_score_best_6);

	$segregation_score_best = max($segregation_score_best_1, $segregation_score_best_2, $segregation_score_best_3, $segregation_score_best_4, $segregation_score_best_5, $segregation_score_best_6, $segregation_score_best_7);
			
}

sub zero_segregation_scores
{
	# Segregation scores
	$segregation_score_1 = 0;
	$segregation_score_2 = 0;
	$segregation_score_3 = 0;
	$segregation_score_4 = 0;
	$segregation_score_5 = 0;
	$segregation_score_6 = 0;
	$segregation_score_7 = 0;
	
	$segregation_score_1B = 0;
	$segregation_score_2B = 0;
	$segregation_score_3B = 0;
	$segregation_score_4B = 0;
	$segregation_score_5B = 0;
	$segregation_score_6B = 0;
	$segregation_score_7B = 0;
	
	$segregation_score_best_1 = 0;
	$segregation_score_best_2 = 0;
	$segregation_score_best_3 = 0;
	$segregation_score_best_4 = 0;
	$segregation_score_best_5 = 0;
	$segregation_score_best_6 = 0;
	$segregation_score_best_7 = 0;
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

sub pause
{
	print "\n >>> Press RETURN to continue\n";
	$answer=<STDIN>;
}