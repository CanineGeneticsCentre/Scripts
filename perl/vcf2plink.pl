#!/usr/bin/perl -w

##################################################################################
#									                                             #      
#	vcf2plink			                                                         #     
#									                                             #
#	This PERL script converts VCF files to a format ready for NGS SNP Handler	 #
#									                                             #
##################################################################################

use strict;
use File::Basename;
use Term::ANSIColor;
use List::Util qw[min max];

# VERSION OF SOFTWARE #
my $version							= "10";

my $show							= "false"; # for debugging
my $show_seg						= "";

# Constants
my $merge_threshold					= 10; # maximum distance apart for merging nearby indels
my $delimiter						= "   ";
my $flanking_rows					= "750"; # FOr the smaller Haploview file. Only 1500 rows altogether, so 750 falnking before and after

# File names
my $vcf_file						= "";
my $command_log						= ""; 
my $tfam_file						= "";
my $tped_file						= "";
my $ped_file						= "";
my $map_file						= "";
my $info_file						= "";
my $disease_status_file				= ""; 
my $tfam_file_haploview				= "";
my $tped_file_haploview				= "";
my $ped_file_haploview				= "";


#Boolean - 'true' or 'false'
my $base_format						= "ACGT"; # or could be ACGT
my $base_format_plink_string		= ""; # string used in plink either --alleleACGT or --alleleACGT
my $indel_format					= "12";   # "12" means REF is 1, rest are 2.  "1234" means REF is 1, rest are 2,3,4 etc etc
my $snpEff_data_found				= "false";
my $VEP_data_found					= "false";
my $start_storing 					= "false";
my $write_to_filtered				= "false"; # decides whether each line is written to the filtered output file
my $include_extra_columns			= "false"; # This is fixed as false (but could be changed here to 'true' to get extra output columns)
my $vep_id_found					= "false"; # whether a VEP position is found at any VCF position
my $use_disease_status_file			= "no"; # 'yes' or 'no'
my $set_carrier_as_normal			= "no"; # set later on
my $miss_this_one					= ""; # 'true' or 'false'
my $handle_triallelic				= ""; # 'remove' or 'convert'
my $make_haploview_file				= ""; # yes or no
my $haploview_centre_not_found		= "false";

#Strings
my $status							= "";
my $family_id						= "";
my $individual_id					= "";
my $snp_name						= "";
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
my $allele_number_A					= "";
my $allele_number_B					= "";
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
my $command							= "";

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
my $vep_id							= "";
my $vep_id_minus_one					= "";
my $vcf_id							= "";
my $vcf_id_plus_one					= "";
my $consequence						= "";
my $consequence_check				= "";
my $consequence_merged				= "";

#Other numbers#
my $first_row						= 0;
my $last_row						= 0;	
my $haploview_centre_position		= 0;
my $miss_this_one_count				= 0; # count if loci are removed because they are tri-allelic
my $count							= 0; # vcf2plink
my $segregation_score_threshold		= 0;
my $line_count						= 0;
my $start_storing_line				= 0; # first line of reall data (after the headers)
my $haploview_centre_row			= 0; #
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
my $second_max_no_of_alleles		= 0;
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
my @genotype_array					= ();
my @base_array						= (); # array of single bases (twice as many as samples)
my @base_number_array				= (); # array of single bases (as 0,1,2 etc) (twice as many as samples)
my @base_array_check				= ();
my @base_orig_array					= (); # array of single bases original (in case they are 'X')
my @affection_array					= (); # 'affected', 'carrier' or 'normal'
my @missing_genotypes_array			= (); # no of missing genotypes for each sample
my @genotypes_counted_array			= (); # no of counted genotypes for each sample
my @vep_array						= (); # Array to whole elements of vVEP variation_id e.g. 7_21821942_G/A
my @indel_file_array				= ();
my @merged_allele_array				= ();
my @merged_genotype_array			= ();
my @vep_file_array					= ();
my @allele_count_array				= (); # for vcf2plink
my @allele_count_array_sorted		= (); # for vcf2plink
my @allele_number_array_1			= (); # for vcf2plink
my @allele_number_array_2			= (); # for vcf2plink



#New arrays for PERL script
#my @QUAL_array						= ();
my @DP_array						= ();
my @GT_array						= ();
my @PL_array						= ();
my @item							= ();
my @item_check						= ();
my @status_line_array				= ();
my @allele_1_for_writing_array		= (); # The base to write for the first allele: A, C, G or T
my @allele_2_for_writing_array		= (); # The base to write for the second allele: A, C, G or T
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
print "      vcf2plink            \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This PERL script converts VCF files to PLINK files\n\n";

print "      You can choose 'ACGT' format or '1234' format\n\n\n";

print color 'bold magenta';
print "  - THERE ARE SOME IMPORTANT THINGS TO NOTE:   \n\n";
print color 'yellow';

print "  - Tri-allelic loci (or more):\n\n";

print "      PLINK can only handle bi-allelic loci, so you have the option of making loci bi-allelic\n";
print "      by excluding any genotypes with the third most common allele, or by excluding tri-allelic loci altogether\n\n\n";

print "  - Indels:\n\n";

print "      Indels are entered as 'A' and 'C' (or '1' and '2' if you choose the '1234' format option)\n";
print "      Again only the two most common alleles are used.\n\n\n";


print "  - Phenotypes:\n\n";

print "      There is an option to add the phenotypes from an external file.\n\n\n";

print "  - The script makes a TPED and TFAM file and then uses 'PLINK --recode' to make a PED and MAP file.\n\n";



 ######   ######## ########    ######## #### ##       ########  ######  
##    ##  ##          ##       ##        ##  ##       ##       ##    ## 
##        ##          ##       ##        ##  ##       ##       ##       
##   #### ######      ##       ######    ##  ##       ######    ######  
##    ##  ##          ##       ##        ##  ##       ##             ## 
##    ##  ##          ##       ##        ##  ##       ##       ##    ## 
 ######   ########    ##       ##       #### ######## ########  ######  


&print_message("Please input the name of the VCF file","input");

until (-e $vcf_file)
{
	print "   VCF file:  ";
	$vcf_file = <STDIN>;
	chomp $vcf_file;
	if ($vcf_file eq ""){$vcf_file = "AS_HC_UG.vcf"}  # TEMP!!!
	if ($vcf_file eq "ls"){print "\n";system ("ls *.vcf");print "\n"}
	if ($vcf_file ne "ls"){if (! -e $vcf_file){print "\n\n>>>>>>>>  File $vcf_file not found.  Try again.  <<<<<<<<\n\n";}}
}

print "\n\n";

$answer="";
until ((substr($answer,0,1) eq "1" ) || (substr($answer,0,1) eq "2" ))
{
	&print_message("Do you want alleles written as A,C,G,T or 1,2,3,4?","input");

	print "   <1> ACGT\n";
	print "   <2> 1234\n\n";

	$answer = <STDIN>;
	chomp $answer;

	if (substr($answer,0,1) eq "1" ){$base_format = "ACGT"; $base_format_plink_string = "--alleleACGT"}
	if (substr($answer,0,1) eq "2" ){$base_format = "1234"; $base_format_plink_string = "--allele1234"}
}



$answer="";
until ((substr($answer,0,1) eq "1" ) || (substr($answer,0,1) eq "2" ))
{
	&print_message("How do you want to handle loci with more than two alleles (PLINK only handles bi-allelic loci)?","input");

	print "   <1> Make them bi-allelic by using the most common two alleles\n";
	print "   <2> Remove tri-allelic alleles altogether\n\n";

	$answer = <STDIN>;
	chomp $answer;

	if (substr($answer,0,1) eq "1" ){$handle_triallelic = "convert"}
	if (substr($answer,0,1) eq "2" ){$handle_triallelic = "remove"}
}


###############################################################
# Ask about whether you want a Haploview compatible file made #
###############################################################
$answer="";
until ((substr($answer,0,1) eq "1" ) || (substr($answer,0,1) eq "2" ))
{
	&print_message("Do you want make a smaller version of the file for Haploview?","input");

	print "   Haploview can only handle about 1500 variants, for the LD mapping, so you\n";
	print "   opt to make a smaller ped file, centred round your SNP of interest.\n\n";

	print "   <1> Yes\n";
	print "   <2> No\n\n";

	$answer = <STDIN>;
	chomp $answer;

	if (substr($answer,0,1) eq "1" ){$make_haploview_file = "yes"}
	if (substr($answer,0,1) eq "2" ){$make_haploview_file = "no"}
}

if ($make_haploview_file eq "yes")
{
	&print_message("What is the position that you want the Haploview file centred round?","input");

	print "(It will have $flanking_rows variants either side of this position)\n\n";

	$haploview_centre_position = <STDIN>;
	chomp $haploview_centre_position;

	unless( $haploview_centre_position =~ /^\d+?$/) # numeric
	{
		print "\n>>>>>> This is not numeric.  No Haploview file will be made.    <<<<<<<<\n\n";

		$make_haploview_file = "no";
	}
}


#################################
# Create names for output files #
#################################
$prefix = &get_prefix ($vcf_file);

$tfam_file = $prefix."_vcf2plink.tfam";
$tped_file = $prefix."_vcf2plink.tped";
$ped_file = $prefix."_vcf2plink.ped";
$map_file = $prefix."_vcf2plink.map";
$info_file = $prefix."_vcf2plink_haploview.info";
$tped_file_haploview = $prefix."_vcf2plink_haploview.tped";
$tfam_file_haploview = $prefix."_vcf2plink_haploview.tfam";
$ped_file_haploview = $prefix."_vcf2plink_haploview.ped";

$command_log = "vcf2plink_"."$prefix"."_command_log.out";


#########################
# Open command log file #
#########################
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output command log file: $command_log";

print COMMAND_LOG "Program vcf2plink\tversion $version\n\n";

if ($handle_triallelic eq "convert")
{
	print COMMAND_LOG "If there are more than 2 alleles at a locus, only genotypes with the most common two alleles will be used\n\n";
}
if ($handle_triallelic eq "remove")
{
	print COMMAND_LOG "If there are more than 2 alleles at a locus, the locus will be removed from the output PLINK files\n\n";
}

print  COMMAND_LOG "Base format chosen:\t$base_format\n\n";


#################################
# Initialise allele_count_array #
#################################
for ($count = 0; $count < 20; $count++)
{
	$allele_count_array[$count] = 0;
}





###########################
# Open the input VCF file MAIN #
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
			
			print "\n\t>>Press RETURN to continue ";
			$answer = <STDIN>;
		}		
	}
	
	
	#########################################################################################
	# This looks at the data in the VCF file from GATK UnifiedGenotyper or HaplotypeCaller  #
	# (if the VCF file was not generated by UnifiedGenotyper or HaplotypeCaller this may    #
	# have to be changed. Basically it gets a list of the files that were                   #
	# used by the UnifiedGenotyper or HaplotypeCaller so it can check against the columns   #
	#########################################################################################
	

	if ((index(lc($single_line),"unifiedgenotyper") > -1) || (index(lc($single_line),"haplotypecaller") > -1))
	{
		#########################################################
        # Check it is a SNP VCF file from whatever...           #
        #########################################################
		if (index(lc($single_line),"genotype_likelihoods_model=snp") > -1)
		{
			$variant_type = "snp";
			$variant_caller = "UnifiedGenotyper";

			&print_message("Variant caller check","message");
			print "The file looks as if comes from the use of the GATK Unified Genotyper ";
			print "to call SNPs (with -glm as SNP).\n\n";
			
			print "\n\t>>Press RETURN if you still want to continue ";
			$answer = <STDIN>;
		}
		if (index(lc($single_line),"genotype_likelihoods_model=indel") > -1)
		{
			$variant_type = "indel";
			$variant_caller = "UnifiedGenotyper";

			&print_message("Variant caller check","message");
			print "The file looks as if comes from the use of the GATK Unified Genotyper ";
			print "to call Indels (with -glm as Indel).\n\n";
			
			print "\n\t>>Press RETURN if you still want to continue ";
			$answer = <STDIN>;
		}

		if (index(lc($single_line),"genotype_likelihoods_model=both") > -1)
		{
			$variant_caller = "UnifiedGenotyper";

			&print_message("Variant caller check","message");
			print "The file looks as if comes from the use of the GATK Unified Genotyper ";
			print "to call SNPs and Indels (with -glm as BOTH).\n\n";
			
			print "\n\t>>Press RETURN if you still want to continue ";
			$answer = <STDIN>;
		}

		if (index(lc($single_line),"haplotypecaller") > -1)
		{
			$variant_caller = "HaplotypeCaller";

			&print_message("Variant caller check","message");
			print "The file looks as if comes from the use of the GATK Haplotype Caller";
			print "to call SNPs and Indels.\n\n";
			
			print "\n\t>>Press RETURN if you still want to continue ";
			$answer = <STDIN>;
		}
		
		################################################
        # Count the number of samples from the UG line #
        ################################################
		for ($pos=0;$pos<=length($single_line);$pos++)
		{
			$char = substr($single_line,$pos,1);
			
			if ($char eq "["){$open_bracket_pos = $pos}
			if ($char eq "]"){$close_bracket_pos = $pos; last}
		}
		
		#####################################################
        # Warn if the list of file names doesn't look right #
        # [file1,file2,file3 etc....]                       #
        #####################################################
		if (($open_bracket_pos == 0 ) || ($close_bracket_pos == 0))
		{
			print "The line in this file containing the words\n";
			print "'UnifiedGenotyper' or 'HaplotypeCaller' should have a list of file names\n";
			print "between square brackets and separated by commas.\n\n";
			
			print "Please check the VCF file\n\n";
			close VCF;
			exit;
		}
		
		
		#############################################################
		# Now get the list of files from this UnifiedGenotyper line #
		# (This is only used to check that it is the same number    #
		# as the number of samples from the #CHROM line)            #
		#############################################################

		$files_string = substr($single_line, $open_bracket_pos + 1, $close_bracket_pos - $open_bracket_pos - 1);
		
		@input_file_array = split (", ",$files_string);
		
		$no_of_samples_UG_line = scalar @input_file_array;
		
		
		
	} # End of: if (index(lc($single_line),"unifiedgenotyper")> - 1)  OR (index(lc($single_line),"haplotypecaller")> - 1)
	

	
	

	 ######  ########    ###    ########  ######## 
	##    ##    ##      ## ##   ##     ##    ##    
	##          ##     ##   ##  ##     ##    ##    
	 ######     ##    ##     ## ########     ##    
	      ##    ##    ######### ##   ##      ##    
	##    ##    ##    ##     ## ##    ##     ##    
	 ######     ##    ##     ## ##     ##    ##    
	 ######  ########  #######  ########  #### ##    ##  ######   
	##    ##    ##    ##     ## ##     ##  ##  ###   ## ##    ##  
	##          ##    ##     ## ##     ##  ##  ####  ## ##        
	 ######     ##    ##     ## ########   ##  ## ## ## ##   #### 
	      ##    ##    ##     ## ##   ##    ##  ##  #### ##    ##  
	##    ##    ##    ##     ## ##    ##   ##  ##   ### ##    ##  
	 ######     ##     #######  ##     ## #### ##    ##  ###### 

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
		

		#########################################################
        # Split whole line at TABs into an array myArray1 (1-9) #
        #########################################################
		
		@myArray1 = split (/\t/,$single_line);
		$array_size_1 = scalar @myArray1;
		$no_of_data_columns = $array_size_1 - 9;
		
		#################################
		# Initialise allele_count_array #
		#################################
		for ($count = 0; $count < 20; $count++)
		{
			$allele_count_array[$count] = 0;
		}

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
		
		##################################
		# Look for central Haploview row #
		##################################
		if ($make_haploview_file eq "yes")
		{
			if ($position =~ /^\d+?$/)
			{
				if ($position == $haploview_centre_position)
				{
					$haploview_centre_row= $line_count;
				}
			}
		}
		# Convert chrX to chr39 for dog 
		if ($chromosome eq "chrX"){$chromosome = "chr39"}
         
		# Remove string 'chr'
        if (substr($chromosome,0,3) eq "chr"){$chromosome = substr($chromosome,3,99)}
		
		# If chromosome is an integer and less than 10 then add a leading zero
		if ($chromosome =~ /^\d+?$/){if ($chromosome < 10){$chromosome = "0".$chromosome}}


		################################################################################
		# Parse myArray1(9)                                                            #
        # These are the genotype columns (one for each sample)                         # 
        # Split the 9th item in the array at COLONS                                    #
        #                                                                              #
        # The 8th column should be GT:AD:DP:GQ:PL but it could be in any order so this #
		# filed is used to determine the order of the data in the 9th item             #
        #                                                                              #
        # Note DP is twice in the row, at the start and in this genotype section.      #
        # They are not necessarily the same.  We use the genotype one.                 #
        ################################################################################
        
         
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
            

            ################################################
            # Get the genotype from myArray[8]             #
            ################################################
            if ($array_size_8 > 0)
            {
                for($sample_count = 1; $sample_count <= $no_of_samples; $sample_count++)
				{
                    $genotype_array[$sample_count] = $myArray1[$sample_count + 8];
                }
            }
			
			################################################
			# Split ALT_base_overall at commas             #
			# to get all the different alleles in an array #
			#                                              #
			# The first ALT allele is not the REF allele   #
			################################################
			@ALT_allele_array = split(",",$ALT_base_overall);
			$no_of_alt_alleles = (scalar @ALT_allele_array);


			##############################################################################
			# Decide whether the variant type is a SNP or an Indel                       # 
			# If any variant is a different length than the REF base then it is an indel #
			##############################################################################
			
			$variant_type = "snp";
			for ($allele_count = 0; $allele_count < $no_of_alt_alleles; $allele_count++)
			{
				if (length($REF_base) != length($ALT_allele_array[$allele_count]))
				{
					$variant_type = "indel";
				}
			}

			#################################
			# Make up name for SNP or Indel #
			#################################
			if ($variant_type eq "snp"){$snp_name = "$chromosome"."_".$position."_SNP";}
			if ($variant_type eq "indel"){$snp_name = "$chromosome"."_".$position."_IND";}


			######################################
			# Write first 4 columns of TPED file #
			######################################

			#######################################################
			# Only do this bit if we are keeping tri-allelic loci #
			#######################################################
			if (($no_of_alt_alleles > 1) && ($handle_triallelic eq "remove"))
			{
				$miss_this_one = "true";
				$miss_this_one_count = $miss_this_one_count + 1;
			}
			else{$miss_this_one = "false"}

			if ($miss_this_one eq "false")
			{
				print TPED "$chromosome$delimiter$snp_name$delimiter"."0"."$delimiter$position";
			}

			$ALT_allele_array[$no_of_alt_alleles] = "Z"; # set extra element put on the end of the array to Z.  This should disappear completely

			#if (($no_of_alt_alleles == 1) && ($variant_type eq "indel")){$show = "true"} else {$show = "false"}

			if ($show eq "true")
			{
				print "Before moving array elements...\n\n";
				for ($count = 0; $count <= $no_of_alt_alleles; $count++)
				{
					print "$count: $ALT_allele_array[$count]\n";
				}
			}

			###############################################
			# Move ALT_alleles along so REF_base is first #
			###############################################
			for ($count = 0; $count <= $no_of_alt_alleles; $count++)
			{
				$array_count = $no_of_alt_alleles - $count;
				$ALT_allele_array[$array_count] = $ALT_allele_array[$array_count -1];
			}

			$ALT_allele_array[0] = $REF_base;
			# NOTE if no_of_ALT_alleles > 1 then PLINK won't handle it.

			if ($show eq "true")
			{
				print "\nAfter moving array elements...\n\n";
				for ($count = 0; $count <=$no_of_alt_alleles; $count++)
				{
					print "$count: $ALT_allele_array[$count]\n";
				}
				$answer=<STDIN>;
			}


			######################################################
			# Make first item in ALT_allele_array the REF_base   #
			# This matches the numbers in the VCF file, i.e. "0" #
			# means the REF base                                 #
			######################################################
			$ALT_allele_array[0] = $REF_base;

			
			$GT_string = "";
			$AD_string = "";
			$DP_genotype = "";
			$GQ_string = "";
			$PL_string = "";
		 
			
            #######################################################################################
            # Now work through all the samples, splitting each separate genotype_array at colons  #
            #######################################################################################
            
            for($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
			{
				
				############################################################
                # Split myArray9 (genotype data) using colons              #
				# Remember if genotype is ./. then there may be no colons  #
				# but splitting should still work if GT is there           #
				############################################################
               				
                @myArray9 = split(":",$genotype_array[$sample_count]);
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
                

 				#####################################################
				# FreeBayes has a single dot . for no genotype      #
				# so there is no slash.  Convert "." to "./." first #
				# so it behaves as GATK VCF files                   #
				#####################################################	

				if ($GT_string eq "."){$GT_string = "./."}


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
				
				# all genotypes should have a slash though....
                if ($slash_pos > 0)
				{
					$allele_number_1 = substr($GT_string,$slash_pos-1,1);  # The allele number comes from 0/1
                    $allele_number_2 = substr($GT_string,$slash_pos+1,1);  

                    #############################################################
                    # Store this to use later when we have decided which to use #
                    #############################################################
                    $allele_number_array_1[$sample_count] = $allele_number_1;
                    $allele_number_array_2[$sample_count] = $allele_number_2;
 

 					###############################################################################
 					# Count how many of each number we get, as we want to use the two most common #
 					###############################################################################
 					if ($allele_number_1 =~ /^\d+?$/) # numeric
					{
	 					$allele_count_array[$allele_number_1] = $allele_count_array[$allele_number_1] + 1;
	 					$allele_count_array[$allele_number_2] = $allele_count_array[$allele_number_2] + 1;
					}
 				} # if slash_pos > 0

 			} # sample_count loop to see which alleles are most common and store all


			###############################################
 			# Now we decide which are the top two alleles #
			###############################################
			$max_no_of_alleles= 0;
			$second_max_no_of_alleles = 0;
			$max_allele_count = 0;
			$second_max_allele_count = 0;

			# First find max allele #
			for ($allele_count = 0; $allele_count <= $no_of_alt_alleles; $allele_count++)
			{
				if ($allele_count_array[$allele_count] > $max_no_of_alleles)
				{
					$max_no_of_alleles = $allele_count_array[$allele_count];
					$max_allele_count = $allele_count;
				}
			}


			# Then find second max allele #
			for ($allele_count = 0; $allele_count <= $no_of_alt_alleles; $allele_count++)
			{
				if ($allele_count != $max_allele_count)
				{	
					if ($allele_count_array[$allele_count] > $second_max_no_of_alleles)
					{
						$second_max_no_of_alleles = $allele_count_array[$allele_count];
						$second_max_allele_count = $allele_count;
						
					}
				}
			}

			############################################
			# Specify the main two alleles             #
			# (what if there is only one ALT allele? ) #
			# REF is other allele                      #
			############################################
			$allele_number_A = $max_allele_count;
			$allele_number_B = $second_max_allele_count;
			$allele_A_original = $ALT_allele_array[$max_allele_count];
			$allele_B_original = $ALT_allele_array[$second_max_allele_count];


			######################################
			# Set what allele_A and allele_B are #
			######################################
			if ($variant_type eq "snp")
			{
				$allele_A = $ALT_allele_array[$max_allele_count];
				$allele_B = $ALT_allele_array[$second_max_allele_count];
			}#snp

			if ($variant_type eq "indel")
			{
				$allele_A = "A";
				$allele_B = "C";
				
			}#indel

			#if (($no_of_alt_alleles == 1) && ($variant_type eq "indel")){$show = "true"} else {$show = "false"} # TEMP!!!!


			########################
			# Show results TEMP!!! #
			########################
			if ($show eq "true")
			{
				print "========================================================\n";
				print "Line $line_count\tvariant_type: $variant_type\n";
				print "No alt alleles: $no_of_alt_alleles\t";

				for ($allele_count = 0; $allele_count <= $no_of_alt_alleles; $allele_count++)
				{
					print "$allele_count: $ALT_allele_array[$allele_count]   ";
				}
				print "REF_base: $REF_base\n\n";
				for($count = 1;$count <=$no_of_samples; $count++)
				{
					print "$count\t$genotype_array[$count]\n";
				}
				print "\n";
				for ($allele_count = 0; $allele_count <= $no_of_alt_alleles; $allele_count++)
				{
					print "Allele $allele_count\tallele_count_array[$allele_count]: $allele_count_array[$allele_count]\n";
				}
				print "\n";
				print "Max no of alleles:    \t$max_no_of_alleles\n";
				print "Second no of alleles: \t$second_max_no_of_alleles\n\n";
				print "Max allele count:     \t$max_allele_count\n";
				print "Second allele count:  \t$second_max_allele_count\n\n";

				print "allele_A_original:    \t$allele_A_original\n";
				print "allele_B_original:    \t$allele_B_original\n";
				print "allele_A:             \t$allele_A\n";
				print "allele_B:             \t$allele_B\n";
				print "allele_number_A:      \t$allele_number_A\n";
				print "allele_number_B:      \t$allele_number_B\n\n";

				#$answer=<STDIN>;
			}

			#############################################################
			# We now want to convert this into actual data to write     #
			# We know what allele_A and allele_B are so we only want to #
			# write these (as PLINK is bi-alleleic only) so others will #
			# be written as -9                                          #
			#############################################################


			#######################################################
			# Only do this bit if we are keeping tri-allelic loci #
			#######################################################
			if (($no_of_alt_alleles > 1) && ($handle_triallelic eq "remove"))
			{
				$miss_this_one = "true";
				$miss_this_one_count = $miss_this_one_count + 1;
			}
			else{$miss_this_one = "false"}

			if ($miss_this_one eq "false")
			{
	 			for($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
				{
					
					if ($variant_type eq "snp") 
					{
						$allele_number_1 = $allele_number_array_1[$sample_count];
		                $allele_number_2 = $allele_number_array_2[$sample_count];

		                # If either are not A or B then both must be.
		                if ((($allele_number_1 ne $allele_number_A)  && ($allele_number_1 ne $allele_number_B)) || (($allele_number_2 ne $allele_number_A)  && ($allele_number_2 ne $allele_number_B)))
		                {
		                	$allele_number_1 = ".";
		                	$allele_number_2 = "."
		                }

						if ($allele_number_1 eq $allele_number_A)
						{
							print TPED "$delimiter"."$allele_A"; 
							if ($show eq "true"){print " $allele_A"}
						}
						elsif
						($allele_number_1 eq $allele_number_B)
						{
							print TPED "$delimiter$allele_B"; 
							if ($show eq "true"){print " $allele_B"}
						}
						else {print TPED "    0";if ($show eq "true"){print " 0"}}


						if ($allele_number_2 eq $allele_number_A)
						{
							print TPED "    $allele_A"; 
							if ($show eq "true"){print "/$allele_A "}
						}
						elsif
						($allele_number_2 eq $allele_number_B)
						{
							print TPED "    $allele_B"; 
							if ($show eq "true"){print "/$allele_B "}
						}
						else {print TPED "    0";if ($show eq "true"){print "/0 "}}
						
					} # SNP


					if ($variant_type eq "indel") 
					{
						$allele_number_1 = $allele_number_array_1[$sample_count];
		                $allele_number_2 = $allele_number_array_2[$sample_count];

						#If either are not A or B then both must be.
		                if ((($allele_number_1 ne $allele_number_A)  && ($allele_number_1 ne $allele_number_B)) || (($allele_number_2 ne $allele_number_A)  && ($allele_number_2 ne $allele_number_B)))
		                {
		                	$allele_number_1 = ".";
		                	$allele_number_2 = "."
		                }

						if ($allele_number_1 eq $allele_number_A)
						{
							print TPED "    $allele_A"; 
							if ($show eq "true"){print " $allele_A"}
						}
						elsif
						($allele_number_1 eq $allele_number_B)
						{
							print TPED "    $allele_B"; 
							if ($show eq "true"){print " $allele_B"}
						}
						else {print TPED "    0";if ($show eq "true"){print " 0"}}


						if ($allele_number_2 eq $allele_number_A)
						{
							print TPED "    $allele_A"; 
							if ($show eq "true"){print "/$allele_A "}
						}
						elsif
						($allele_number_2 eq $allele_number_B)
						{
							print TPED "    $allele_B"; 
							if ($show eq "true"){print "/$allele_B "}
						}
						else {print TPED "    0";if ($show eq "true"){print "/0 "}}

					} # indel

				} 
				#Next sample_count - loop to process genotypes
				
				#print "=================================\n";
				

				#############################
				# Now the final end of line #
				#############################
				print TPED "\n";
				if ($show eq "true"){print "\n"}

			} # if miss_this_one eq "false"

        } #End If 'If Ubound > 9


 		if ($show eq "true"){$answer=<STDIN>}
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

		$start_storing_line = $line_count;
		
		############################################################
		# Parsing of #CHROM data line to get a list of input files #
		############################################################
		@chrom_line_array = split(/\s+/,$single_line);
		$chrom_line_array_size = (scalar @chrom_line_array) - 1;
		$no_of_samples_chrom_line = $chrom_line_array_size - 8; # ignore first columns 0-9
		
		for ($array_count = 9; $array_count <= $chrom_line_array_size; $array_count++)
		{
			# Sample names go into sample_name_array
			$sample_name_array[$array_count - 8] = $chrom_line_array[$array_count];
			$file_count = $array_count - 8;
		}
		
		

		##############################################################
        # The 'unified_genotyper' line has a list of the input files #
        # The '#CHROM' line has a list of the column headings        #
        # This second is the order to use (they MAY BE DIFFERENT)    #
        ##############################################################
        if (($variant_caller eq "UnifiedGenotyper") || ($variant_caller eq "HaplotypeCaller"))
        {
        	if ($no_of_samples_UG_line != $no_of_samples_chrom_line)
			{
				&print_message("POSSIBLE ERROR!","warning");
				print "There is a problem with this VCF file.\n\n";
				print "There are $no_of_samples_UG_line samples in the list of input file names.\n";
				print "There are $no_of_samples_chrom_line columns in the file (from #CHROM line).\n\n";
				print "Something is wrong. Please check\n\n";

				print "\t>>If you want to carry on anyway, press 'return' to continue .\n\n";
		
        		$answer = <STDIN>;
			}
        
	        if ($no_of_samples_UG_line == $no_of_samples_chrom_line)
			{
				$no_of_samples = $no_of_samples_chrom_line;
			}
		} # 
         


		###################################################################
		# If it is not a GATK program UnifiedGenotyper or HaplotypeCaller #
		# then you have to ignore the no_of_samples_UG_line variable      #
		###################################################################
		if ($no_of_samples_UG_line == 0)
		{
			$no_of_samples = $no_of_samples_chrom_line;
		}



		#################################
		# Show user a list of the files #
		#################################
		
		&print_message("Sample columns found in the VCF file (from the '#CHROM' line)","message");
		
		for ($sample_count = 1; $sample_count <= $no_of_samples; $sample_count++)
		{
			print "$sample_count \t$sample_name_array[$sample_count]\n";
		}
		
		print "\nThere are $no_of_samples samples in this multi-column VCF file.\n\n";
		print "\t>>If these columns look OK, press 'return' to continue .\n\n";
		
        $answer = <STDIN>;
		
		
		#######################################################
		# Ask if you want to use disease statuses from a file #
		#######################################################

		&print_message("Do you want to add disease statuses from a file","input");

		print "   <1> Yes\n";
		print "   <2> No\n\n";

		$answer = <STDIN>;
		chomp $answer;

		if (substr($answer,0,1) eq "1" ){$use_disease_status_file = "yes"}
		if (substr($answer,0,1) eq "2" ){$use_disease_status_file = "no"}


		if ($use_disease_status_file eq "no")
		{
			print "\n\t>>>>>>>>  All disease statuses will be set to zero (unknown)\n\n";
			print COMMAND_LOG "Diseae statuses not used. All disease statuses set to zero\n\n";
		}


		########################################
		# Get the disease statuses from a file #
		########################################
		if ($use_disease_status_file eq "yes")
		{
			&print_message("Please input the name of the file with the disease statuses of the samples","input");

			print "   This has two columns separated by TAB. The first is the Sample Name as in the VCF file, the \n";
			print "   second column has 'affected' or 'normal' (or 'case' or 'control')\n\n";
		
			until (-e $disease_status_file)
			{
				print "   Disease status file: ";
				$disease_status_file = <STDIN>;
				chomp $disease_status_file;

				if ($disease_status_file eq ""){$disease_status_file = "AS_HC_disease_statuses.txt"} # TEMP!!!
				if ($disease_status_file eq "ls"){print "\n";system ("ls *.txt");print "\n"}


				if (($disease_status_file ne "ls") && ($disease_status_file ne "help")){if (! -e $disease_status_file){print "\n\n>>>>>>>>  File $disease_status_file not found.  Try again.  <<<<<<<<\n\n";}}

				if ($disease_status_file eq "help")
				{
					&print_message("HELP ON DISEASE STATUS FILE","help");
					print "\tThe disease status file has two columns, separated by a TAB.\n\n";
					print "\tThe first column has the name of the sample, as used as the column header in the VCF file.\n";
					print "\tThe second column must be 'Affected', 'Carrier' or 'Normal'\n";
					print "\t------------------------------------------------------------------------------------------\n\n";
				}
			}

		} # use_disease_status_file

	



		############################
		# Read Disease Status file #
		############################

		if ($use_disease_status_file eq "yes")
		{

			open (STATUS_FILE, "$disease_status_file") || die "Cannot open $vcf_file";
			$disease_status_count = 0;
			
			while ($disease_status_line = <STATUS_FILE>) 
			{
				chomp $disease_status_line;
				$disease_status_count = $disease_status_count + 1;
				
				@status_line_array = split(/\s+/,$disease_status_line);
				$status_line_array_size = scalar @status_line_array;

				if ($status_line_array_size == 1)
				{
					&print_message("There should be two columns in this file","warning");
					close STATUS_FILE;
					exit;
				}

				if ($status_line_array_size == 2)
				{
					$sample_name_from_status_file = $status_line_array[0];
					$disease_status = lc $status_line_array[1];
					
					
					##########################################
					# Warn if disease status is unrecognised #
					##########################################
					if (($disease_status ne "affected") &&
						($disease_status ne "case") &&
						($disease_status ne "carrier") &&
						($disease_status ne "control") &&
						($disease_status ne "clear") &&
						($disease_status ne "normal"))
					{
						&print_message("Warning. Disease status not recognised","warning");
						
						print "Disease status $disease_status is not a recognised type\n\n";
						print "It must be 'Affected', 'Carrier' or 'Normal'\n\n\n";
						close STATUS_FILE;
						exit;
					}
						
					###########################################
					# Rationalise disease status descriptions #
					###########################################
					if (($disease_status eq "affected") || ($disease_status eq "case"))
					{
						$disease_status = "affected";
					}
					if (($disease_status eq "clear") || ($disease_status eq "normal") || ($disease_status eq "control"))
					{
						$disease_status = "normal";
					}
					

					if (($disease_status eq "carrier") && ($set_carrier_as_normal ne "yes"))
					{
						&print_message("For PLINK, the status must be 'affected' or 'normal'","warning");

						print "Sample $sample_name_from_status_file  is set as $disease_status\n\n";

						print "Do you want to set the carriers as normals?     \n\n";

						print "   <1> YES\n";
						print "   <2> NO\n\n";

						$answer = <STDIN>;
						chomp $answer;

						$set_carrier_as_normal = "yes";

						if (substr($answer,0,1) eq "1" ){$set_carrier_as_normal = "yes"}
						if (substr($answer,0,1) eq "2" ){$set_carrier_as_normal = "no"}

						if ($set_carrier_as_normal eq "no")
						{
							print ">>>>>>>>>  Quitting program\n\n";
							exit;
						}
						if ($set_carrier_as_normal eq "yes")
						{
							print "\n\t>>>>>>>> Setting disease status for sample $sample_name_from_status_file from $disease_status to ";
							$disease_status = "normal";
							print "$disease_status\n";
						}
					} # carrier

					if (($disease_status eq "carrier") && ($set_carrier_as_normal eq "yes"))
					{
						print "\n\t>>>>>>>> Setting disease status for sample $sample_name_from_status_file from $disease_status to ";

						$disease_status = "normal";

						print "$disease_status\n";
						
					} # carrier


					$disease_status_array[$disease_status_count] = $disease_status;
					$disease_status_sample_array[$disease_status_count] = $sample_name_from_status_file;
				
				} # if array size = 2

			}	# reading STATUS_FILE
			
			$no_samples_in_status_file = $disease_status_count;
			
		

		
			############################################
			# Show user a list of the disease statuses #
			############################################
			&print_message ("List of samples and disease statuses","message");

			print COMMAND_LOG "List of samples and disease statuses:\n\n";

			print "\n";
			for ($disease_status_count = 1;$disease_status_count <=$no_samples_in_status_file; $disease_status_count++)
			{
				print "\t$disease_status_count\t$disease_status_sample_array[$disease_status_count]\t$disease_status_array[$disease_status_count]\n";
				print COMMAND_LOG "\t$disease_status_count\t$disease_status_sample_array[$disease_status_count]\t$disease_status_array[$disease_status_count]\n";
			}
			
			if ($no_samples_in_status_file != $no_of_samples)
			{
				&print_message("FATAL ERROR!  Numbers of samples don't match.","warning");
				print "There are $no_of_samples samples in the VCF file and $no_samples_in_status_file samples in the Disease Status file\n\n";
				exit;
			}
			
			
			############################################################################################
			# Check each name in the VCF file only matches once to the list in the Disease Status file #
			############################################################################################
			for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
			{
				$name_match_count = 0;
				
				for ($disease_status_count = 1;$disease_status_count <=$no_samples_in_status_file; $disease_status_count++)
				{
					if ($disease_status_sample_array[$disease_status_count] eq $sample_name_array[$sample_count])
					{
						$name_match_count = $name_match_count + 1; # Checks each name only matches once
					}
				}
				# Check each name only matches once
				if ($name_match_count == 1){$total_match_count = $total_match_count + 1}
				if ($name_match_count > 1){&print_message("Too many names match!  Check file","warning");exit}
			} # sample_count loop
			
			if ($total_match_count != $no_of_samples)
			{
				&print_message("ERROR!  All names in Disease status file do not match those in VCF file.","warning");
				
				print "Sample columns found in the VCF file (from the '#CHROM' line)\n\n";
				
				for ($sample_count = 1; $sample_count <= $no_of_samples; $sample_count++)
				{
					print "\t$sample_count \t$sample_name_array[$sample_count]\n";
				}

				
				print "\n\nSample names in the disease_status file\n\n";

				for ($disease_status_count = 1;$disease_status_count <=$no_samples_in_status_file; $disease_status_count++)
				{
					print "\t$disease_status_count\t$disease_status_sample_array[$disease_status_count]\t$disease_status_array[$disease_status_count]\n";
				}

				print "\nOnly $total_match_count names out of $no_of_samples match\n\n";
				exit;
			}
			if ($total_match_count == $no_of_samples)
			{
				print "\n\t>>If these look correct, press 'RETURN' ";
				$answer=<STDIN>;
				
				&print_message("All names in Disease status file match those in VCF file.","message");
				
				print "\t($total_match_count names out of $no_of_samples match)\n\n";
			}
			
		
		


			######  #######       #     #    #    ######  
			#     # #             ##   ##   # #   #     # 
			#     # #             # # # #  #   #  #     # 
			######  #####   ##### #  #  # #     # ######  
			#   #   #             #     # ####### #       
			#    #  #             #     # #     # #       
			#     # #######       #     # #     # #   

			#############################################################################################
			#  Start to 'Re-map' the columns so that the order is Affected -> Carrier -> Normal         #
			#                                                                                           #
			# This is then used by starting with the new order and looking up a "source column         #
			# from the VCF file to get the data from.                                                   #
			#                                                                                           #
			# So if your second affected sample was in column 5 of the VCF file then the new order      #  
			# count would be 2 and the $source_col would be 5.  Data would be transferred from column 5 #   
			# and put into new column 2.                                                                #
			#############################################################################################
			$new_sample_count = 0;
			
			# Get affecteds		
			for ($disease_status_count = 1;$disease_status_count <=$no_samples_in_status_file; $disease_status_count++)
			{
				if ($disease_status_array[$disease_status_count] eq "affected")
				{
					$no_affected_samples = $no_affected_samples + 1;
					
					# Find this name in VCF file array (sample_name_array)
					for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
					{
						if ($disease_status_sample_array[$disease_status_count] eq $sample_name_array[$sample_count])
						{
							$new_sample_count = $new_sample_count + 1;

							$source_column_array[$new_sample_count] = $sample_count; # Stores array of source columns
							$sample_status_array[$sample_count] = "affected";
						} # if names match
					}
				} # If affected
			} # Affecteds
			
			# Get carriers		
			for ($disease_status_count = 1;$disease_status_count <=$no_samples_in_status_file; $disease_status_count++)
			{
				if ($disease_status_array[$disease_status_count] eq "carrier")
				{
					$no_carrier_samples = $no_carrier_samples + 1;
					
					# Find this name in VCF file array (sample_name_array)
					for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
					{
						if ($disease_status_sample_array[$disease_status_count] eq $sample_name_array[$sample_count])
						{
							$new_sample_count = $new_sample_count + 1;

							$source_column_array[$new_sample_count] = $sample_count; # Stores array of source columns
							$sample_status_array[$sample_count] = "carrier";
						} # if names match
					}
				} # If carrier
			} # Carriers
			
			# Get Normals		
			for ($disease_status_count = 1;$disease_status_count <=$no_samples_in_status_file; $disease_status_count++)
			{
				if ($disease_status_array[$disease_status_count] eq "normal")
				{
					$no_normal_samples = $no_normal_samples + 1;
					
					# Find this name in VCF file array (sample_name_array)
					for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
					{
						if ($disease_status_sample_array[$disease_status_count] eq $sample_name_array[$sample_count])
						{
							$new_sample_count = $new_sample_count + 1;

							$source_column_array[$new_sample_count] = $sample_count; # Stores array of source columns
							$sample_status_array[$sample_count] = "normal";			
						} # if names match
					}
				} # If normal
			} # Normals

			
			print COMMAND_LOG "\nDisease status file: \t$disease_status_file\n\n";

			&print_both ("\tNumber of affecteds: \t$no_affected_samples\n");
			&print_both ("\tNumber of carriers:  \t$no_carrier_samples\n");
			&print_both ("\tNumber of normals:   \t$no_normal_samples\n\n");

		} #if ($use_disease_status_file eq "yes")

		

		#
		# If there are no disease statuses #
		####################################
		if ($use_disease_status_file eq "no")
		{
			for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
			{
				$sample_status_array[$sample_count] = "unknown";
			}
		}



		##############################
		# Set various things to zero #
		##############################
		for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
		{
			$missing_genotypes_array[$sample_count] = 0;
			$genotypes_counted_array[$sample_count] = 0;
		}
		
        ######################################
        # Open the output files              #
        ######################################

		

		open (TPED, ">$tped_file")|| die "Cannot create output file: $tped_file";


		
		
		################################################################
		# Write headers for the SIX columns in the output TFAM file    #
		################################################################
		#print TPED "CHR\tSNP\tcM\tPOS";


		############################################
        # Now add headers for the output TPED file  #   TEMP!!!
        ############################################
        for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
		{

            ##########################################################
            # Add name of sample as header in first batch of columns #
            # (ie those after the first six pedfile columns          #
            # with the ALT base for each sample                      #
            ##########################################################
            
            $sample_name = &get_prefix($sample_name_array[$sample_count]);
			
            #print TPED "\t$sample_name_array[$sample_count]"; # Once for allele_1 
			#print TPED "\t$sample_name_array[$sample_count]"; # Once for allele_2

		}
		
	
		
		#############################
		# Now the final end of line #
		#############################
		
		#print TPED "\n";


		&print_message("Processing VCF file...","message");

	} # End of: if (index($single_line,"#CHROM") > -1)
	
	
	if (($line_count % 10000) == 0)
	{
		if ($start_storing eq "true"){print "Reading VCF file.  Line: $line_count\n"}
	}
	 
} # End of: while ($single_line = <VCF>) 


close TPED;


###########################################################
# Tell user of haploview centre row has been found or not #
###########################################################

if ($haploview_centre_row > 0)
{
	&print_both("\nHaploview centre variant found at row $haploview_centre_row\n\n");
}
else
{
	&print_both("\nHaploview centre variant position not found in TPED file\n\n");
	&print_both("Haploview-compatible files will not be created.\n\n");
	$make_haploview_file = "no";
	$haploview_centre_not_found = "true";
}




################################################################
# Write headers for the SIX columns in the output TFAM file    #
# These are the usual six columns in a PED file
# FID, IID, PA, MA, Sex, STATUS
################################################################

print "\nWriting to TFAM file $tfam_file...\n\n";

open (TFAM, ">$tfam_file")|| die "Cannot create output file: $tfam_file";

#print TFAM "TFAM1\tTFAM2\tTFAM3\tTFAM4\tTFAM5\tTFAM6\n";


for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
{
    
    $sample_name = &get_prefix($sample_name_array[$sample_count]);
	
	$family_id = "F1";
	$individual_id = $sample_name_array[$sample_count];

	if ($sample_status_array[$sample_count] eq "affected"){$status = "2"} 
	elsif ($sample_status_array[$sample_count] eq "normal"){$status = "1"} 
	elsif ($sample_status_array[$sample_count] eq "unknown"){$status = "0"} 
	else {$status = "0"}

    print TFAM "$family_id\t$individual_id\t0\t0\t0\t$status\n";
}

close TFAM;


##################################################################
# If the option of making a small Haploview file has been chosen #
##################################################################

if ($make_haploview_file eq "yes")
{
	# Open the TPED file
	open (TPED, "$tped_file")|| die "Cannot open file: $tped_file";
	open (HAP, ">$tped_file_haploview")|| die "Cannot create output file: $tped_file_haploview";

	$first_row = $haploview_centre_row - $flanking_rows;
	$last_row = $haploview_centre_row + $flanking_rows;
	$line_count = 0;

	&print_message("Making file compatible with Haploview","message");

	print "   First_row:   \t$first_row\n";
	print "   Centre row:  \t$haploview_centre_row\n";
	print "   Last row:    \t$last_row\n\n";


	$line_count = 0;
	while ($single_line = <TPED>)
	{
		$line_count = $line_count + 1;
		#chomp $single_line;
		#&chomp_all ($single_line);

		if (($line_count > $first_row)  && ($line_count < $last_row))
		{
			print HAP $single_line;
		}

	}

	close TPED;
	close HAP;
}


##################################################
# Convert transposed files to normal PLINK files #
##################################################
$prefix=&get_prefix($tped_file);

$command = "plink --noweb --dog --allow-no-sex --make-founders --nonfounders --tfile $prefix --recode $base_format_plink_string --out $prefix";
&print_message("Converting TPED and TFAM files to PED and MAP...","message");
print "$command\n";
print COMMAND_LOG "$command\n";

system("$command");


####################################################
# Make files for Haploview - smaller TPED and INFO #
####################################################
if ($make_haploview_file eq "yes")
{
	# Make PED and INFO files

	$prefix=&get_prefix($tped_file_haploview);

	$command = "plink --noweb --dog --allow-no-sex --make-founders --nonfounders --alleleACGT --tped $tped_file_haploview --tfam $tfam_file --recodeHV $base_format_plink_string --out $prefix";
	&print_message("Converting TPED and TFAM files to PED and MAP for Haploview compatibility...","message");
	print "$command\n";
	print COMMAND_LOG "$command\n";

	system("$command");

} # Make Haploview files


#############################
# Write details to log file #
#############################
print COMMAND_LOG "\n\n####################\n";
print COMMAND_LOG "# Analysis details #\n";
print COMMAND_LOG "####################\n\n";

&print_message("FINISHED ANALYSIS","message");


&print_both("Input file:                        \t$vcf_file\n\n");

&print_both("Output TPED file:                  \t$tped_file\n");
&print_both("Output TFAM file:                  \t$tfam_file\n\n");

&print_both("Output PED file:                   \t$ped_file\n");
&print_both("Output MAP file:                   \t$map_file\n\n");


if ($make_haploview_file eq "yes")
{
	&print_both("\nFiles for Haploview:\n\n");

	&print_both("Output PED file:                   \t$ped_file_haploview\n");
	&print_both("Output INFO file:                  \t$info_file\n\n");
}
if ($haploview_centre_not_found eq "true")
{
	&print_both("Haploview centre position not found\n");
	&print_both("so no Haploview output files were created.\n\n");
}
#print "Output INFO file (for Haploview):  \t$info_file\n\n";

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
		$col_1 = ($sample_test_count * 2 ) - 1;
		$col_2 = ($sample_test_count * 2 );
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		 print "\t$allele_base_1$allele_base_2";
	}
	print "\n\tCarriers: ";
	for ($sample_test_count = $no_affected_samples + 1; $sample_test_count <= ($no_affected_samples + $no_carrier_samples); $sample_test_count++)
	{
		$col_1 = ($sample_test_count * 2 ) - 1;
		$col_2 = ($sample_test_count * 2 );
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		print "\t$allele_base_1$allele_base_2";
	}
	print "\n\tNormals: ";
	for ($sample_test_count = $no_affected_samples + $no_carrier_samples + 1; $sample_test_count <= $no_of_samples; $sample_test_count++)
	{
		$col_1 = ($sample_test_count * 2 ) - 1;
		$col_2 = ($sample_test_count * 2 );
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
		$col_1 = ($sample_test_count * 2 ) - 1;
		$col_2 = ($sample_test_count * 2 );
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		 print "\t$allele_base_1-$allele_base_2 ";
	}
	if ($no_carrier_samples > 0)
	{
		print "\nGenotypes carrier: ";
		for ($sample_test_count = $no_affected_samples + 1; $sample_test_count <= ($no_affected_samples + $no_carrier_samples); $sample_test_count++)
		{
			$col_1 = ($sample_test_count * 2 ) - 1;
			$col_2 = ($sample_test_count * 2 );
			$allele_base_1 = $base_array[$col_1];
			$allele_base_2 = $base_array[$col_2];
			print "\t$allele_base_1-$allele_base_2 ";
		}
	}
	print "\nGenotypes normal: ";
	for ($sample_test_count = $no_affected_samples + $no_carrier_samples + 1; $sample_test_count <= $no_of_samples; $sample_test_count++)
	{
		$col_1 = ($sample_test_count * 2 ) - 1;
		$col_2 = ($sample_test_count * 2 );
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
		$col_1 = ($check_count * 2 ) - 1;
		$col_2 = ($check_count * 2 );
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		 print "\t$allele_base_1-$allele_base_2";
	}
	print "\n\tNormals: ";
	for ($check_count = $no_affected_samples + 1; $check_count <= $no_of_samples; $check_count++)
	{
		$col_1 = ($check_count * 2 ) - 1;
		$col_2 = ($check_count * 2 );
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
		$col_1 = ($check_count * 2 ) - 1;
		$col_2 = ($check_count * 2 );
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
		$col_1 = ($check_count * 2 ) - 1;
		$col_2 = ($check_count * 2 );
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


sub get_sample_cols
{
	my $sample_count_sub = $_[0];

	$col_1 = ($sample_count_sub * 2) - 1;
	$col_2 = ($sample_count_sub * 2);
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
	print "\n Press RETURN to continue\n";
	$answer=<STDIN>;
}