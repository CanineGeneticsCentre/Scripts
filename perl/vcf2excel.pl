#!/usr/bin/perl -w

#######################################################
#									                  #      
#	vcf2excel			                              #     
#									                  #
#	This PERL script converts VCF files to a format   #
#   ready for the Excel sheet 'NGS SNP Viewer'	      #
#									                  #
#     Mike Boursnell: mike.boursnell@aht.org.uk       #    
#######################################################

use strict;
#use File::Basename; # Is this used?
use Term::ANSIColor;
use List::Util qw[min max];


# VERSION OF SOFTWARE #
my $version							= "80";

my $show							= "false"; # for debugging

#DEBUGGING
my $position_to_debug				= ""; # Can add a genome position here for debugging, e.g. set as "12013456"

#Constants
my $merge_threshold					= 10; # maximum distance apart for merging nearby indels
my $single_case_scoring				= "true"; # In cases where there is one case and many controls ($only_one_case is true), use a special scoring system
my $call_frequency_threshold		= 30; # Call frequencies below this can be filtered out.

my $include_extra_columns_to_debug	= "false"; # This is fixed as false (but could be changed here to 'true' to get extra output columns for debugging)
my $mark_missing_as_X				= "true"; # Controls when genotype has not been called. This is now fixed as 'true'
my $score_X_as_REF					= "true"; # If you do use X, do you count these as REF for scoring?  This is now fixed as 'true'
my $score_X_method					= "best";  # 'reference', best' or 'missing'
my $filter_multiple_alleles			= "false"; # if a variant has more than one ALT allele, handle like any other variant. This is now fixed as 'false'
my $keep_unfiltered_output			= "false"; # Whether or not to keep the very large unfiltered output file (if this is set to "" then the user has a choice)
my $keep_indels_only_file			= "false"; # Whether or not to keep the indels_only output file  (if this is set to "" then the user has a choice)
my $make_merged_indels_file			= "false"; # Experimental merged indels output file  (if this is set to "" then the user has a choice)
my $combine_filters					= "AND"; # AND or OR  (if this is set to "" then the user has a choice)

#File names
my $vcf_file						= ""; # A VCF file, processed by VEP or snpEff
my $disease_status_file				= ""; # The second most important file. A list of which samples are Affected, Carriers, Normals (or Omit)
my $command_log						= "";
my $missing_effects					= ""; # File to record any efects which are not assigned a score
my $output_file						= "";
my $output_file_filtered			= ""; # output file filtered for best segregations and effects
my $output_file_filtered_and_sorted	= ""; # output file filtered for best segregations and effects, and then sorted in unix
my $indels_merged_file				= ""; # output file indels only

my $simplified_indels_output_file	= ""; # Indels output file for use as input for merging indels
my $snp_list_file					= ""; # file with a lists of known SNPs
my $ensembl_names_file				= ""; # file with a list of ensembl names and equivalent gene names"

my $temp_indels_file				= ""; # TEMPORARY!!
my $check_output_file				= "check_output.out";

# Debugging and checking
my %check_output_genotype_hash		= ();
my $check_output_string				= "";

#Boolean - 'true' or 'false' (or 'yes' or 'no')
my $snpEff_data_found				= "false";
my $LOF_data_found					= "false";
my $NMD_data_found					= "false";
my $vep_data_found					= "false";
my $passed_header_lines 			= "false"; # Used to be start storing
my $write_to_filtered				= "false"; # decides whether each line is written to the filtered output file
my $filter_unknown_chromosomes		= "false"; # Could make this a user choice
my $include_base_change_info		= "false"; # Include extra details in the output fil, base/amino changes, CDS/protein positions
my $on_snp_list						= "true"; # "true" or "false"
my $pause_cancelled					= "false"; # works with the 'pause_with_cancel' subroutine
my $vcf_file_closed					= "false"; # true or false
my $miss_this_one					= "false"; # Boolean to decide if next bit of code is missed
my $sift_deleterious_found			= "false"; # Another criterion for inclusion in the filtered file

my $monomorphic						= "false"; # Are all the genotypes for the samples in this VCF file the same?
my $omit_low_call_frequencies		= "false"; # Omit variants with low caal frequencies from the filtered output (see $call_frequency_threshold)
my $only_one_case					= "false"; # If thee is only one case and many controls, the seg scores are calculated differently

my $check_snp_list					= "yes"; # Yes or No.  Check each SNP to filtered file is on known SNP list. Dog only
my $convert_ensembl_names			= "yes"; # Yes or No.  Convert ensembl names to genes names (if snpEff or VEP doesn't) Dog only

#Strings
my $run_title						= ""; # this becomes part of the name of the output file
my $vcf_line						= ""; # variable for readsing single line of text files
my $check_line						= "";
my $snp_list_line					= "";
my $ensembl_names_line				= "";
my $variant_caller					= ""; # e.g. 'UnifiedGenotyper', 'HaplotypeCaller etc', freeBayes
my $variant_caller_from_info		= ""; # variant caller read from the INFO field (if put there by pool_vcf_files)
my $depth_coverage_from_info		= "";
my $ref_seq_name					= "";
my $answer							= "";
my $default_vep_file				= ""; 
my $char							= "";
my $prefix							= "";
my $prefix_disease_status_file		= "";
my $sample_name						= "";
my $chromosome						= "";
my $position						= "";
my $snp_name						= "";
my $position_check					= "";
my $REF_base						= "";
my $ALT_base_overall				= "";
my $QUAL							= "";
my $FILTER							= "";
my $GT_string 						= "";
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
my $allele_base_1					= ""; # actual base of the genotype, A,C, G or T
my $allele_base_2					= "";
my $allele_base_check_1				= ""; # actual base of the genotype, A,C, G or T when looking for merged indels
my $allele_base_check_2				= "";
my $allele_base_merged_1			= ""; # actual base of the genotype, A,C, G or T when looking for merged indels
my $allele_base_merged_2			= "";
my $allele_A						= "";
my $allele_B						= "";
my $allele_C						= ""; # more than 3 alleles isn't dealt with.
my $main_affected_allele			= ""; # if variant is SNP then this is the actual base
my $main_affected_allele_check		= ""; 
my $second_affected_allele			= "";
my $main_normal_allele				= ""; # IS THIS NEEDED?
my $genotype_AB						= ""; # Genotype in the form 'AB' rathet than using the actual bases e.g. 'CT' 'AT'  'CG' etc
my $main_affected_genotype			= "";
my $single_case_genotype			= ""; # If there is only one case, this is the genotype of that one 
my $genotype_check					= "";
my $original_genotype				= "";
my $genotype_sample_1				= "";
my $merged_genotype					= "";
my $swapped							= "";
my $disease_status_line				= "";
my $disease_status					= "";
my $sample_name_from_status_file 	= "";
my $effect_predictor				= ""; # variant_effect_predictor (VEP) or snpEff
my $variant_type					= ""; # snp, indel (or both?)
my $vep_gene						= ""; # VEP ensembl gene name
my $command							= "";
my $strand							= "1"; # This is 1 or -1 for the SIFT output file
my $chr_pos							= ""; # chromsome and position separated by an underscore (e.g. 12_459022)
my $ensembl_gene_name				= "";
my $real_gene_name					= "";
my $search_string					= ""; #
my $VEP_fields_string				= ""; # String holding the list of fields used by VEP (in VCF header)
my $snpEff_fields_string			= ""; # String holding the list of fields used by snpEff (in VCF header)
my $vep_check_string				= ""; # Used to check whether a set of effect data is merely a dulication of previous data onthat line
my $vep_check_string_previous		= ""; # This is a string made up of a previous effect to check with current string
my $snpEff_check_string				= ""; # Used to check whether a set of effect data is merely a dulication of previous data onthat line
my $snpEff_check_string_previous	= ""; # This is a string made up of a previous effect to check with current string
my $start_time						= "";
my $end_time						= "";
my $run_time						= "";
my $seg_score_filter_choice			= "best"; # Later to be used to choose which score you filter by.  Is this worthwhile?
my $x_chromosome_number				= ""; # This is set to "chr39" for dog (for example)

# snpEff and VEP variables
my $snpEff_effect					= ""; # The main name for the snpEff effect
my $vep_effect						= ""; # VEP effect or consequence
my $snpEff_effect_full_string		= ""; # holds string of effects from snpEff
my $LOF_full_string					= ""; # holds LOF data from snpEff
my $NMD_full_string					= ""; # holds NMD data from snpEff
my $vep_effect_full_string			= ""; # holds string of effects from VEP
my $snpEff_effect_string			= ""; # holds individual effects from snpEff
my $vep_effect_string				= ""; # holds individual effects from VEP
my $snpEff_effect_impact 			= "";
my $vep_effect_impact	 			= "";
my $snpEff_effect_codon_change		= "";
my $vep_effect_codon_change			= "";
my $snpEff_effect_gene_name 		= "";
my $vep_effect_gene_name 			= "";
my $snpEff_effect_gene_symbol		= ""; #
my $vep_effect_gene_symbol			= ""; #
my $snpEff_effect_real_gene_name	= ""; # Gene name looked up in ensembl genes names file (Dog only)
my $vep_effect_real_gene_name		= ""; # Gene name looked up in ensembl genes names file (Dog only)
my $snpEff_effect_base_change		= "";
my $vep_effect_base_change			= "";
my $snpEff_effect_amino_acid_change	= "";
my $vep_effect_amino_acid_change	= "";
my $snpEff_CDS_position				= "";
my $vep_CDS_position				= "";
my $snpEff_protein_position			= "";
my $vep_protein_position			= "";
my $snpEff_effect_transcript 		= "";
my $vep_effect_transcript 			= "";
my $vep_sift_prediction				= ""; # SIFT prediction is an option for VEP
my $snpEff_sift_prediction			= ""; # not currently an option for snpEff
my $LOF_string						= "";
my $LOF_gene						= "";
my $LOF_gene_string					= "";
my $NMD_string						= "";
my $NMD_gene						= "";
my $NMD_gene_string					= "";
my $snpEff_results_string			= ""; # String which is written to output file for consequence column in Excel
my $vep_results_string				= ""; # String which is written to output file for consequence column in Excel
my $snpEff_gene_string				= ""; # String which is written to output file for gene column in Excel
my $snpEff_base_change_string		= ""; # String which is written to output file for base change column in Excel
my $vep_base_change_string			= ""; # String which is written to output file for base change column in Excel
my $snpEff_amino_acid_change_string	= ""; # String which is written to output file for amino acid change column in Excel
my $vep_amino_acid_change_string	= ""; # String which is written to output file for amino acid change column in Excel
my $snpEff_CDS_position_string		= "";
my $vep_CDS_position_string			= ""; 
my $snpEff_protein_position_string	= "";
my $vep_protein_position_string		= ""; 
my $vep_gene_string					= "";
my $vep_sift_prediction_string		= "";
my $consequence						= "";
my $consequence_check				= "";
my $consequence_merged				= "";
my $main_affected_homozygous_genotype = "";

#Numbers
# Counters
my $line_count						= 0;
my $col_count						= 0;
my $check_count						= 0;
my $vep_line_count					= 0;
my $sample_count					= 0;
my $new_sample_count				= 0;
my $sample_test_count				= 0;
my $array_count						= 0;
my $ampersand_array_count			= 0;
my $file_count						= 0;
my $disease_status_count			= 0;
my $name_match_count				= 0; # checks number of times each name in VCF matches with name in Disease Status files
my $total_match_count				= 0; # checks how many names in VCF matche with name in Disease Status files
my $allele_count					= 0;
my $max_allele_count				= 0;
my $second_max_allele_count			= 0;
my $NMD_count						= 0;
my $LOF_count						= 0;
my $max_seg_score_count				= 0;
my $second_seg_score_count			= 0;
my $third_seg_score_count			= 0;
my $effect_count					= 0;
my $filter_count_no_of_alleles		= 0; # Counts the number of SNPs added to the filtered results because of no of alleles
my $filter_count_seg_score			= 0; # Counts the number of SNPs added to the filtered results because of segregation score best
my $filter_count_effect_score		= 0; # Counts the number of SNPs added to the filtered results because of effect score
my $filter_count_combined_score		= 0; # Counts the number of SNPs added to the filtered results for seg AND effect score
my $filter_count_all				= 0; # Sum of all the above counts
my $effect_score_of_5_count			= 0;
my $effect_score_of_4_count			= 0;
my $effect_score_of_3_count			= 0;
my $deleterious_less_than_5_count	= 0;
my $vcf_variant_count				= 0;
my $omit_count						= 0; # How many samples are labelled "omit" in the disease status file
my $monomorphic_count				= 0;
my $call_freq_omit_count			= 0; # Count of variants omitted from filtered output because of low call frequency
my $read_depth_omit_count			= 0; # Count of variants omitted from filtered output because of low cread depth
my $affected_count_ds_file			= 0;
my $normal_count_ds_file			= 0;
my $missing_effects_count			= 0;
my $max_both_scores_count			= 0;
my $ensembl_gene_replaced_count		= 0; # Counts when an ensembl gene name is replaced by the real name

#Array sizes
my $array_size						= 0;
my $snpEff_items_array_size			= 0;
my $vep_items_array_size			= 0;
my $split_ampersand_array_size		= 0;
my $LOF_items_array_size			= 0;
my $NMD_items_array_size			= 0;
my $chrom_line_array_size			= 0;
my $INFO_field_array_size			= 0;
my $array_size_1					= 0;
my $array_size_8					= 0;
my $array_size_9					= 0;
my $status_line_array_size			= 0;

#Total numbers etc
my $no_of_samples_chrom_line		= 0;
my $no_of_samples					= 0; # final number of samples
my $no_of_samples_without_omits		= 0; # final number of samples minus the number omitted
my $no_of_data_columns				= 0;
my $no_of_snpEff_effects			= 0;
my $no_of_snpEff_effects_full		= 0; # counting those linked by '&' as well
my $no_of_LOFs						= 0;
my $no_of_NMDs						= 0;
my $no_of_vep_effects				= 0;
my $no_of_alt_alleles				= 0;
my $no_of_alleles					= 0;
my $no_of_alleles_check				= 0;
my $max_no_of_alleles				= 0;
my $second_no_of_alleles			= 0;
my $no_affected_samples				= 0;
my $no_normal_samples				= 0;
my $no_carrier_samples				= 0;
my $no_omit_samples					= 0;
my $no_samples_in_status_file		= 0;

#Scores
my $minimum_segregation_score		= 0;
my $max_seg_score					= 0;
my $second_seg_score				= 0;
my $third_seg_score					= 0;
my $minimum_effect_score			= 3;
my $effect_score					= 0;
my $effect_score_check				= 0; # When checking for merged indels
my $effect_score_merged				= 0; #
my $vep_effect_score				= 0;
my $vep_effect_score_ampersand_max	= 0;
my $snpEff_effect_score				= 0;
my $snpEff_effect_score_ampersand_max= 0;
my $max_effect_score   				= 0; # If VEP and snpEff are both used
my $max_vep_effect_score   			= 0;
my $max_snpEff_effect_score   		= 0;
my $highest_possible_segregation_score		= 0;

#Column positions
my $col_1							= 0;
my $col_2							= 0;
my $effect_score_col 				= 0;
my $seg_score_best_col 				= 0;
my $source_col						= 0;

#Other numbers
my $pos								= 0;
my $max_sample_name_length			= 0;
my $slash_pos						= 0;
my $homozygosity_ratio				= 0;
my $homozygosity_ratio_check		= 0;
my $missing_genotypes_total			= 0;
my $missing_genotypes_percent		= 0;
my $genotypes_counted_total			= 0;
my $genotypes_counted_per_sample	= 0;
my $no_of_indels					= 0; # number of indels in the separate indel file (or array)
my $missing_genotypes_at_this_position = 0; # This is used to determine the call frequency at this position
my $call_frequency_at_this_position = 0;    # This is the call frequency at this position (percentage of good calls)
my $IMPACT_field_order				= 0;  # These are the order of fields specified by snpEff or VEP in the VCF file header
my $SYMBOL_field_order				= 0;  # These are the order of fields specified by snpEff or VEP in the VCF file header
my $Gene_field_order				= 0;  # These are the order of fields specified by snpEff or VEP in the VCF file header
my $CDS_position_field_order		= 0;  # These are the order of fields specified by snpEff or VEP in the VCF file header
my $Protein_position_field_order	= 0;  # These are the order of fields specified by snpEff or VEP in the VCF file header
my $Amino_acids_field_order			= 0;  # These are the order of fields specified by snpEff or VEP in the VCF file header
my $Codons_field_order				= 0;  # These are the order of fields specified by snpEff or VEP in the VCF file header
my $SIFT_field_order				= 0;  # These are the order of fields specified by snpEff or VEP in the VCF file header

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
my $segregation_score_best_overall	= 0; # Best overall
my $segregation_score_to_filter_by	= 0; # Default is best, but in this version you can choose
	
#Statistics	
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
my $MM_count						= 0; #Count of missing genotypes
my $normalise_missing_factor		= 0;

#Arrays
my @sample_name_array				= (); # names of the samples (from the #CHROM line)
my @chrom_line_array				= (); # array of items on line with #CHROM in it.
my @disease_status_sample_array 	= (); # array of sample names from disease status file
my @disease_status_array 			= (); # array of disease_statuses from disease status file
my @myArray1						= ();
my @INFO_field_array				= (); # contains snpEff results (if they are in the VCF file)
my @myArray8						= ();
my @myArray9						= ();
my @ALT_allele_array				= ();
my @genotype_array					= ();
my @base_array						= (); # array of single bases (twice as many as samples)
my @base_array_check				= ();
my @base_orig_array					= (); # array of single bases original (in case they are 'X') # Is this still used?
my @missing_genotypes_array			= (); # no of missing genotypes for each sample
my @genotypes_counted_array			= (); # no of counted genotypes for each sample
my @indel_file_array				= ();
my @merged_allele_array				= ();
my @merged_genotype_array			= ();
my @split_ampersand_array			= (); # when snpEff effect is like this: splice_acceptor_variant&splice_donor_variant
my @GT_string_array					= (); # Stores the GT_string for each sample for checking output
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
my @vep_items_array					= (); # Holds the different items in a snpEff effect, e.g. gene name
my @LOF_items_array					= (); # Holds the different items in a LOF field, e.g. gene name
my @LOF_array						= (); # Holds LOF (Loss Of Function) data
my @NMD_items_array					= (); # Holds the different items in a NMD field, e.g. gene name
my @NMD_array						= (); # Holds NMD (Nonsense Mediated Decay) data
my @VEP_fields_array				= (); # Holds the list and order of fields used by VEP
my @snpEff_fields_array				= (); # Holds the list and order of fields used by VEP

# Arrays in original input order
# Note.  The order of samples in the original VCF file may not be the same
# as the order required in the output file for excel.  In the output file
# we will want the cases in the first columns and the controls in the last columns
#
my @genotype_array_input_order		= (); # genotypes in original input order (@genotype_array is in output order or disease_status order)
my @sample_name_array_input_order	= (); # names of the samples (from the #CHROM line in original input order)
my @sample_status_array_input_order	= (); # disease statuses of samples in original input order


# New approach where each gene is linked to all other factors
# For each gene, the effect predictors (VEP and snpEff) give a list of data
# We want to store all the odata for eacg gene before deciding which bits
# to print to the output file
#
# snpEff arrays
my @snpEff_single_genes_array				= ();
my @snpEff_single_effects_array				= ();
my @snpEff_single_effect_scores_array		= ();
my @snpEff_single_base_change_array			= ();
my @snpEff_single_amino_acid_change_array	= ();
my @snpEff_single_CDS_position_array		= ();
my @snpEff_single_protein_position_array	= ();
my @snpEff_single_LOF_array					= ();
my @snpEff_single_NMD_array					= ();

# VEP arrays
my @vep_single_genes_array					= ();
my @vep_single_effects_array				= ();
my @vep_single_effect_scores_array			= ();
my @vep_single_base_change_array			= ();
my @vep_single_amino_acid_change_array		= ();
my @vep_single_CDS_position_array			= ();
my @vep_single_protein_position_array		= ();
my @vep_single_sift_prediction_array		= ();
my @vep_single_LOF_array					= ();
my @vep_single_NMD_array					= ();

# Hash tables
my %indel_hash								= (); # Hash table for indel letters "a", "b", "c" etc (looked up with actual sequence e.g. might be TAAA)
my %merged_hash								= (); # merged indels
my %snp_list_hash							= (); # List of SNPs in known SNP list file
my %ensembl_names_hash						= (); # List of ensembl gene names with real gene name (dog only)

print color 'reset';
print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      vcf2excel             \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This PERL script converts VCF files to a format ready for the Excel sheet NGS SNP VIEWER\n\n";
print "    It creates the final alleles page directly from the VCF file, and also\n";
print "    calculates all the possible segregation scores.\n\n";
print "    It also adds the variant annotations (effects) from snpEff or variant_effect_predictor (VEP).\n\n";
print "    If the VCF file has been processed by VEP it will extract the SIFT predictions as well.\n\n";

print color 'reset';

print color 'bold green';

print "    The main affected allele is called 'A' and all other alleles are called 'B',\n";
print "    and the segregation score is worked out using these 'AA,'AB', and 'BB' genotypes..\n\n";

if ($make_merged_indels_file eq "true")
{
	print "    It also makes a 'merged indels' file which merges adjacent indels if they are within $merge_threshold bases of each other.\n\n";
}

print color 'yellow';
print "Files required:\n\n";
print "\t1.  Multi-column VCF file\n";
print "\t2.  Disease status file (sample names with disease statuses - 'Affected', 'Carrier', 'Normal' or 'Omit')\n";
print "\t         [use 'make_file_of_file_names' for this]\n\n";

if ($position_to_debug ne "")
{
	&print_message("This is set to stop at position $position_to_debug for debugging","warning");
	$answer=<STDIN>;
}

###########################################################
# Get the name of the input VCF file                      #
###########################################################
$vcf_file = &get_input_file ("What is the name of the input VCF file","vcf");


###########################################################
# Open the input VCF file and check the file (including   #
# headers) then close once checks have been made          #
###########################################################
&check_VCF_file;

if ($vcf_file_closed eq "false"){print "CLOSE!!!!";close VCF_TEMPORARY;}


###########################################################
# Warn if no effect predictor annotations have been found #
###########################################################
if ($effect_predictor eq "")
{
	&print_message("Neither snpEff nor VEP annotations have been detected in this VCF file","warning");
	print "Make sure the original VCF is processed by run_snpEff or run_VEP\n\n";
	exit;
}


###########################################################
# Ask user if they are using VEP  or snpEff               #
###########################################################
if ($effect_predictor eq "both"){&get_effect_predictor;}


###########################################################
# Ask about how you want to deal with missing genotypes   #
# This is now removed as an option (2/9/2015)             #
# There is no reason no to use the defaults.              #
###########################################################
#&ask_about_dealing_with_missing_genotypes;


###########################################################
# Ask about how you want to deal with missing genotypes   #
###########################################################
&ask_about_extra_columns;


###########################################################
# Load in lists of existing SNPs and ensembl gene symbols #
###########################################################
&load_existing_snp_list;




###########################################################
# Open the input VCF file                                 #
# This is the MAIN LOOP of the script                     #
###########################################################
open (VCF, "$vcf_file") || die "Cannot open $vcf_file";
$line_count = 0;
$vep_line_count = 1;
$passed_header_lines = "false";


while ($vcf_line = <VCF>) 
{
	chomp $vcf_line;
	&chomp_all ($vcf_line);

	$line_count = $line_count + 1;

	
	#######################################################################################
	# Only do this next bit if 'passed_header_lines' is True, i.e. if the line containing #   
	# #CHROM has been found and the real variant data starts (i.e. past the headers)      #
	#######################################################################################
	if ($passed_header_lines eq "true")
	{
		
		#############################################################
		# Set various things to blank or zero for this variant line #
		#############################################################
		&set_all_to_zero;
		

		#########################################################################
		# Initially we assume this is not to be written to the filtered file    #
		# The criteria are:                                                     #
		#                                                                       #
		#  1.  If segregation_score_best is >= a user set threshold             #
		#                                                                       #
		#  2.  If there is a snpEff or VEP effect score > a user set threshold  #
		#                                                                       #
		#  3.  If there are more than one ALT allele                            #
		#                                                                       #
		#########################################################################
		$write_to_filtered = "false";


		###########################################################################
		# Read part of the line from the VCF file into these global variables:    #
		# Column	Variable                                                      #
		# 1   		$chromosome                                                   #
		# 2   		$position                                                     #
		# 4   		$REF_base                                                     #
		# 5   		$ALT_base_overall                                             #
		# 6   		$QUAL                                                         #
		# 7   		$FILTER                                                       #
		#$#########################################################################
		&read_vcf_first_seven_columns;
    

		############################################################################################
		# myArray1 contains the whole VCF line split at TABs                                       #
		# The INFO field is column 8 [array element 7]                                             #
		# If the INFO field is present get data from there: contains genotypes and snpEff/VEP info #
		############################################################################################
		&get_data_from_VCF_INFO_field (\@myArray1);

	
		###############################################################
		# Get highest effect score from VEP and snpEff (if both used) #
		###############################################################
		if ($max_vep_effect_score > $max_snpEff_effect_score)
		{$max_effect_score = $max_vep_effect_score} else {$max_effect_score = $max_snpEff_effect_score}



		################################################################################
		# Parse myArray1(9)   These are the genotype columns (one for each sample)     # 
        # Split the 9th item in the array at COLONS                                    #
        #                                                                              #
        # The 8th column should be GT:AD:DP:GQ:PL but it could be in any order so this #
		# field is used to determine the order of the data in the 9th item             #
        #                                                                              #
        # Note DP (depth) occurs twice in the row, at the start and in this genotype   #
        # section. They are not necessarily the same.  We use the genotype one.        #
        ################################################################################
        
         
		##########################################################
		# If the 9th item is present - get the genotypes from it #
		##########################################################
        if ($array_size_1 > 8)
		{  
            #######################################################################################
			# If the 9th item is present get data from there (this includes the GENOTYPES)        #
			# The data goes into an array called '@genotype_array[$sample_count]'                 #
			#######################################################################################
			&get_data_from_VCF_GENOTYPES_field (\@myArray1); 

			
			#####################################################
			# Split ALT_base_overall at commas                  #
			# to put all the different alleles in an array      #
			#####################################################
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
				
				#############################################################################
				# Not sure why this is here. I assume that the ALT_allele can't be . or ./. #
				#############################################################################
				if ($ALT_allele_array[$array_count] eq ".")
				{
					print " >ALT_allele_array[$array_count]: >$ALT_allele_array[$array_count]<\n";
					$answer=<STDIN>;
				}
				elsif ($ALT_allele_array[$array_count] eq "./.")
				{
					print "  >>ALT_allele_array[$array_count]: >$ALT_allele_array[$array_count]<\n";
					$answer=<STDIN>;
				}

				$indel_hash{$ALT_allele_array[$array_count]} = chr($array_count + 98);

				if ($position eq $position_to_debug)
				{
					print "==============================================================================================\n";
					print "This is just after where the indel_hash gets assigned\n\n";

					print "Array count: $array_count out of $no_of_alt_alleles\n";
					print "   ALT_allele_array[$array_count]: $ALT_allele_array[$array_count]\n\n";

					print "indel_hash{$ALT_allele_array[$array_count]}: $indel_hash{$ALT_allele_array[$array_count]}\n";

					#print " indel_hash{$base_array[$col_1]}: $indel_hash{$base_array[$col_1]}\n";
					#print " indel_hash{$base_array[$col_2]}: $indel_hash{$base_array[$col_2]}\n\n";

					print " base_array[$col_1]: $base_array[$col_1]\n";
					print " base_array[$col_2]: $base_array[$col_2]\n\n";

					print "  base_orig_array[$col_1]: $base_orig_array[$col_1]\n";
					print "  base_orig_array[$col_2]: $base_orig_array[$col_2]\n\n";

					print "   allele_base_1: $allele_base_1\n";
					print "   allele_base_2: $allele_base_2\n\n";

					for($sample_count = 1; $sample_count <= $no_of_samples; $sample_count++)
					{
				        print "  genotype_array[$sample_count]: $genotype_array[$sample_count]\n";
				    }

					print "==============================================================================================\n";
					$answer=<STDIN>;
				}
			} # array_count loop

			
			#############################################################################################
		    # Now work through all the samples, splitting each separate genotype_array at colons        #
		    # The key bit of data is put into GT_string and then into two base_arrays:                  #
		    # $base_array[$col_1] = $allele_base_1;                                                     #
			# $base_array[$col_2] = $allele_base_2;                                                     #
		    #############################################################################################
         	&get_genotypes_from_all_samples (\@genotype_array);
			
			
			
			#########################################################################
			# New way to get main affected allele (works for any number of alleles) #
			# This doesn't need allele_A and allele_B  .                            #
			#                                                                       #
			# When analysing for segregation score:                                 #
			#    Main affected allele is called 'A'.  All the rest are called 'B'   #
			#    The variable with these in is $genotype_AB                         #                  
			#########################################################################
			&get_main_affected_allele; # for this VCF row
 			

 			########################################################################################
 			#  Choose the alleles for all the samples and save them in allele_for_writing_array    #
 			# For the SNPs store the actual base (ACGT) and for the indels use a code a,b,c,d etc  #
 			# Note: the segregation score comparison are done on the genotypes in the form 'AA',   #
 			# 'AB' etc and the 'allele_for_writing' is the actual base (A, C, G or T) for SNPs     #
 			# and 'a', 'b', 'c' etc for indels.                                                    #
 			########################################################################################
 			&choose_alleles_for_all_samples;

 			
			##########################################################################
			# Now score affected GENOTYPES.                                          #
			# HH_count is the most common homozygous allele in the affected.         #
			# At this point the counts of single alleles such as A_count_affected    #
			# are not used any more.                                                 #
			##########################################################################
			&score_affected_genotypes;


			##########################################################################
			# Now score CARRIER genotypes.                                           #
			# (if there are any)                                                     #
			##########################################################################
			if ($no_carrier_samples > 0)
			{
				&score_carrier_genotypes;
			}


			##########################################################################
			# Now score NORMAL genotypes. (i.e. non-affecteds including carriers)    #
			# HH_count is the most common homozygous allele in the normal.           #
			##########################################################################
			&score_normal_genotypes;

		 
		 
			##########################################################################
			# Format numbers before writing to output file                           #
			##########################################################################
			$affected_homozygosity = sprintf("%.1f", $affected_homozygosity);
			$normal_homozygosity = sprintf("%.1f", $normal_homozygosity);
			$homozygosity_ratio = sprintf("%.1f", $homozygosity_ratio);
			
			
			##########################################################################################
			# The segregation score" depends on which allele is chosen as the "main affected allele" #
			#                                                                                        #
			# So two different segregation scores are chosen, which are called "segregation_score"   #
			# and "segregation_score_B".                                                             #
			#                                                                                        #
			# The highest of these two segregation scores is chosen.                                 #
			##########################################################################################
			&choose_best_segregation_score_A_or_B;
			
			
			##################################################################
			# Other sums (only needed if extra de-bugging columns are used   #
			# so could be cut for speed                                      #
			##################################################################
			$AA_AB_count_affected = $AA_count_affected + $AB_count_affected;
			$AA_AB_count_normal = $AA_count_normal + $AB_count_normal;
			
			$AB_BB_count_affected = $AB_count_affected + $BB_count_affected;
			$AB_BB_count_normal = $AB_count_normal + $BB_count_normal;
	  


			########################################################################################## 
			# Check if all the genotypes are the same (if so don't write to filtered file)           # 
			##########################################################################################   
			$genotype_sample_1 = "$allele_1_for_writing_array[1]$allele_2_for_writing_array[1]";
			$monomorphic = "true";

			for ($sample_count = 2; $sample_count <= $no_of_samples_without_omits; $sample_count++)  # Started at 2. New 14/10/2015
			{
				$genotype_check = "$allele_1_for_writing_array[$sample_count]$allele_2_for_writing_array[$sample_count]";
				if ($genotype_check ne $genotype_sample_1)
				{
					$monomorphic = "false";
					last;
				}
			} # sample_count loop


			########################################################################################## 
			# Now write the first three chr, pos and ref columns to the OUT file (unfiltered)        # 
			# followed by the columns with the actual alleles                                        #  
			########################################################################################## 
			if ($keep_unfiltered_output eq "true")  
			{
				print "keep_unfiltered_output: $keep_unfiltered_output\n";
				$answer=<STDIN>;

				print OUT "$chromosome\t$position\t$REF_base";

				for ($sample_count = 1; $sample_count <= $no_of_samples_without_omits; $sample_count++)
				{
					print OUT "\t$allele_1_for_writing_array[$sample_count]\t$allele_2_for_writing_array[$sample_count]"; # New 14/10/2015
				} # sample_count loop
			

				#####################################################################
				# Print out all the segregation scores                              #
				#####################################################################
				print OUT "\t$segregation_score_best_1";
				print OUT "\t$segregation_score_best_2";
				print OUT "\t$segregation_score_best_3";
				print OUT "\t$segregation_score_best_4";
				print OUT "\t$segregation_score_best_5";
				print OUT "\t$segregation_score_best_6";
				print OUT "\t$segregation_score_best_7";
				print OUT "\t$segregation_score_best_overall"; # Best overall
				

				#####################################################################
				# Print "R" if the main affected allele is the same as the REF base #
				#####################################################################
				if ($main_affected_allele ne $REF_base){print OUT "\tN"} else {print OUT "\tR";}
				

				#####################################################################
				# Record whether SNP is on the SNP list                             #
				#####################################################################
				if ($on_snp_list eq "false")
				{ print OUT "\tN"; }
				else
				{ print OUT "\t$on_snp_list"; }


				#####################################################################
				# Print the maximum effect score                                    #
				#####################################################################
				print OUT "\t$max_effect_score"; # Effect score
				

				#####################################################################
				# Print VEP or snpEff results                                       #
				#####################################################################
				if ($effect_predictor eq "VEP") 
				{
					print OUT "\t$vep_results_string"; # Consequence
					print OUT "\t$vep_gene_string"; # Gene
				} 
				elsif ($effect_predictor eq "snpEff") 
				{
					print OUT "\t$snpEff_results_string"; # Consequence
					print OUT "\t$snpEff_gene_string"; # Gene
				} 
				


				#####################################################################
				# Print the LOF and NMD from snpEff or the SIFT results from VEP    #
				#####################################################################
				if ($effect_predictor eq "VEP")
				{
					print OUT "\t$vep_sift_prediction_string";
					print OUT "\t";	
				}
				elsif ($effect_predictor eq "snpEff")
				{
					print OUT "\t$LOF_gene_string";
					print OUT "\t$NMD_gene_string";						
				}
				


				#####################################################################
				# Put extra columns in here with base and amino acid change         #
				# details (but if the user specifies this)                          #
				# If not specified then the column are omitted altogether)          #
				#####################################################################
				if ($include_base_change_info eq "true")
				{ 
					if ($effect_predictor eq "VEP")
					{
						print OUT "\t$vep_base_change_string"; # Shows which bases were changed
						print OUT "\t$vep_amino_acid_change_string"; # Shows which amino acids were changed
						print OUT "\t$vep_CDS_position_string"; # Shows base position in CDS (coding sequence)
						print OUT "\t$vep_protein_position_string"; # Shows amino acid psoition in protein
					}
					elsif ($effect_predictor eq "snpEff")
					{
						print OUT "\t$snpEff_base_change_string"; # Shows which bases were changed
						print OUT "\t$snpEff_amino_acid_change_string"; # Shows which amino acids were changed
						print OUT "\t$snpEff_CDS_position_string"; # Shows base position in CDS (coding sequence)
						print OUT "\t$snpEff_protein_position_string"; # Shows amino acid psoition in protein
					}
					
				} # if ($include_base_change_info eq "true")


				#####################################################################
				# Print the rest of the columns                                     #
				#####################################################################
				print OUT "\t$depth_coverage_from_info";
				print OUT "\t$call_frequency_at_this_position";
				print OUT "\t$no_of_alt_alleles";
				print OUT "\t$ALT_base_overall";
				print OUT "\t$homozygosity_ratio";
				print OUT "\t$variant_caller_from_info";
				print OUT "\t$variant_type";
				
				if ($variant_type eq "indel"){print OUT "\t$indel_hash{$main_affected_allele}"} else {print OUT "\t$main_affected_allele"}

			} # if ($keep_unfiltered_output eq "true") 


			##################################################################### 
			# Now write to the OUT_INDELS file, but only if the                 # 
			# variant is an Indel                                               #
			##################################################################### 
			if ($make_merged_indels_file eq "true")
			{			
				if ($variant_type eq "indel")
				{
					print INDELS_OUTPUT "$chromosome\t$position\t$REF_base";

					for ($sample_count = 1; $sample_count <= $no_of_samples_without_omits; $sample_count++)
					{
						print INDELS_OUTPUT "\t$allele_1_for_writing_array[$sample_count]";
						print INDELS_OUTPUT "\t$allele_2_for_writing_array[$sample_count]";
					}

					#####################################################################
					# Print out all the segregation scores                              #
					#####################################################################
					print INDELS_OUTPUT "\t$segregation_score_best_1";
					print INDELS_OUTPUT "\t$segregation_score_best_2";
					print INDELS_OUTPUT "\t$segregation_score_best_3";
					print INDELS_OUTPUT "\t$segregation_score_best_4";
					print INDELS_OUTPUT "\t$segregation_score_best_5";
					print INDELS_OUTPUT "\t$segregation_score_best_6";
					print INDELS_OUTPUT "\t$segregation_score_best_7";		 
				 
					print INDELS_OUTPUT "\t$segregation_score_best_overall"; # Best overall
					
					if ($main_affected_allele ne $REF_base){print INDELS_OUTPUT "\tN"} else {print INDELS_OUTPUT "\tR";}
					
					print INDELS_OUTPUT "\t$max_effect_score"; # Effect score
					
					if ($effect_predictor eq "VEP") 
					{
						print INDELS_OUTPUT "\t$vep_results_string"; # Consequence
						print INDELS_OUTPUT "\t$vep_gene_string"; # Consequence
					} 
					elsif ($effect_predictor eq "snpEff") 
					{
						print INDELS_OUTPUT "\t$snpEff_results_string"; # Consequence
						print INDELS_OUTPUT "\t$snpEff_gene_string"; # Consequence
					}

					

					print INDELS_OUTPUT "\t$ALT_base_overall";
					
					print INDELS_OUTPUT "\t$allele_A";
					print INDELS_OUTPUT "\t$allele_B";
					
					print INDELS_OUTPUT "\t$indel_hash{$main_affected_allele}";

					print INDELS_OUTPUT "\t$homozygosity_ratio";
					print INDELS_OUTPUT "\t$no_of_alt_alleles";
					print INDELS_OUTPUT "\t$variant_caller_from_info";
					print INDELS_OUTPUT "\t$depth_coverage_from_info";
					print INDELS_OUTPUT "\t$call_frequency_at_this_position";
					print INDELS_OUTPUT "\n";

				} # If variant is Indel

			} # if ($make_merged_indels_file eq "true")


			if ($keep_unfiltered_output eq "true")
			{ 
				####################################################################
				# These extra columns can be included in the output for de-bugging #
				####################################################################
				if ($include_extra_columns_to_debug eq "true")
				{
					print OUT "\t$main_affected_homozygous_genotype";
					print OUT "\t$AA_count_affected";
					print OUT "\t$AB_count_affected";
					print OUT "\t$BB_count_affected";
					print OUT "\t$AA_AB_count_affected";
					print OUT "\t$AB_BB_count_affected";
					print OUT "\t$HH_count_affected";
					
					print OUT "\t$AA_count_normal";
					print OUT "\t$AB_count_normal";
					print OUT "\t$BB_count_normal";
					print OUT "\t$AA_AB_count_normal";
					print OUT "\t$AB_BB_count_normal";
					print OUT "\t$HH_count_normal";

					print OUT "\t$affected_homozygosity";
					print OUT "\t$normal_homozygosity";

				} # Include extra columns for de-bugging

				
				####################################################################
				# Now the final end of line                                        #
				####################################################################
				print OUT "\n";

			} #if ($keep_unfiltered_output eq "true")


			######################################################################## 
			# Check the criteria for writing to the FILTERED file                  #
			########################################################################   
		
			if ($position eq $position_to_debug)
			{
				print "Line 1069\n";
				print "\ncombine_filters: $combine_filters\n";
				print "sift_deleterious_found: $sift_deleterious_found\n";
				print "write_to_filtered: $write_to_filtered\n\n";
			}


			if ($position eq $position_to_debug)
			{
				print "77777777777777777777777777777777777777777777777777777777777777777\n\n";
				print "max_effect_score:            \t$max_effect_score\n";
				print "minimum_effect_score:        \t$minimum_effect_score\n";
				print "seg_score_to_filter_by:      \t$segregation_score_to_filter_by\n";
				print "minimum_segregation_score:   \t$minimum_segregation_score\n";
				print "combine_filters:             \t$combine_filters\n";
				print "sift_deleterious_found:      \t$sift_deleterious_found\n";
				print "write_to_filtered:           \t$write_to_filtered\n\n";

				&show_segregation_scores;
				print "77777777777777777777777777777777777777777777777777777777777777777\n\n";
				&pause;
			}

			
			###########################################################################
			# If filtered variants must meet seg_score AND effect score criteria      #
			###########################################################################
			if ($combine_filters eq "AND")
			{
				# Segregation score criterion
				if (($segregation_score_to_filter_by >= $minimum_segregation_score) && (($max_effect_score >= $minimum_effect_score)  || ($sift_deleterious_found eq "true")))
				{
					$write_to_filtered = "true";
					$filter_count_combined_score = $filter_count_combined_score + 1;
				}
			} # combine_filters eq "AND" (--> smaller output files)
			
			###########################################################################
			# If filtered variants can meet either seg_score OR effect score criteria #
			###########################################################################
			elsif ($combine_filters eq "OR")
			{
				# Segregation score criterion
				if ($segregation_score_to_filter_by >= $minimum_segregation_score)
				{
					$write_to_filtered = "true";
					$filter_count_seg_score = $filter_count_seg_score + 1;
				}

				# Variant effect criterion
				if ($max_effect_score >= $minimum_effect_score) 
				{
					$write_to_filtered = "true";
					$filter_count_effect_score = $filter_count_effect_score + 1;
				}

				# SIFT deleterious criterion
				if ($sift_deleterious_found eq "true")
				{
					$write_to_filtered = "true";
				}
			} # combine_filters eq "OR" (--> large output files)
			
				 

			###########################################################################
			# Check the number of ALT alleles criterion                               #
			# $filter_multiple_alleles is now fixed as "false"                        #
			###########################################################################
			#if ($filter_multiple_alleles eq "true")
			#{
			#	if ($no_of_alt_alleles > 1)
			#	{
			#		$write_to_filtered = "true";
			#		$filter_count_no_of_alleles = $filter_count_no_of_alleles + 1;
			#	}
			#}

			#if ($filter_multiple_alleles eq "snps_only")
			#{
			#	if (($no_of_alt_alleles > 1) && ($variant_type eq "snp"))
			#	{
			#		$write_to_filtered = "true";
			#		$filter_count_no_of_alleles = $filter_count_no_of_alleles + 1;
			#	}
			#}


			###########################################################################
			# Unknown Chromosomes                                                     #
			# This is currently fixed at "false" and any Unknown chromosomes still    #
			# get written to the filtered file                                        # 
			###########################################################################
			#if ($filter_unknown_chromosomes eq "true")
			#{
			#	if (index($chromosome,"Un") > -1)
			#	{
			#		$write_to_filtered = "false";
			#	}
			#}


			###########################################################################
			# Record more details of numbers of variants passing each criterion       #
			# ($effect_score_of_4_count means "at least" a score of 4)                #
			###########################################################################
			#Effect score criterion
			if ($max_effect_score == 5)
			{
				$effect_score_of_5_count = $effect_score_of_5_count + 1;
				$effect_score_of_4_count = $effect_score_of_4_count + 1;
				$effect_score_of_3_count = $effect_score_of_3_count + 1;
				
				if ($segregation_score_best_1 == $no_of_samples_without_omits)
				{
					$max_both_scores_count = $max_both_scores_count + 1;
				}
			}
			elsif ($max_effect_score == 4)
			{
				$effect_score_of_4_count = $effect_score_of_4_count + 1;
				$effect_score_of_3_count = $effect_score_of_3_count + 1;
			}
			elsif ($max_effect_score == 3)
			{
				$effect_score_of_3_count = $effect_score_of_3_count + 1;
			}

			# Segregation score criterion
			if ($segregation_score_best_overall == $no_of_samples_without_omits)
			{
				$max_seg_score_count = $max_seg_score_count + 1;
			}
			elsif ($segregation_score_best_overall == $no_of_samples_without_omits - 1)
			{
				$second_seg_score_count = $second_seg_score_count + 1;
			}
			elsif ($segregation_score_best_overall == $no_of_samples_without_omits - 2)
			{
				$third_seg_score_count = $third_seg_score_count + 1;
			}

			
			###########################################################################
			# Remove lines in which Call Frequency is very low (if user has opted)    #
			###########################################################################
			if ($omit_low_call_frequencies eq "true")
			{
				if ($call_frequency_at_this_position < $call_frequency_threshold)
				{
					$write_to_filtered = "false"; $call_freq_omit_count = $call_freq_omit_count + 1;
				}
			}


			###########################################################################
			# If all the genotypes in the samples in the VCF line are the same        #
			# then we don't want it in the filtered file.                             #
			###########################################################################
			if ($monomorphic eq "true"){$write_to_filtered = "false"; $monomorphic_count = $monomorphic_count + 1;}


			###########################################################################
			# Now if "write_to_filtered" is true                                      #
			# write to the filtered output file                                       #
			###########################################################################
			if ($write_to_filtered eq "true")
			{
				$filter_count_all = $filter_count_all + 1;
				
				print OUT_FILTERED "$chromosome\t$position\t$REF_base";

				for ($sample_count = 1; $sample_count <= $no_of_samples_without_omits; $sample_count++)
				{
					print OUT_FILTERED "\t$allele_1_for_writing_array[$sample_count]";
					print OUT_FILTERED "\t$allele_2_for_writing_array[$sample_count]";
				}

				#####################################################################
				# Print out all the segregation scores                              #
				#####################################################################
				#print "  >> Just before writing to OUT_FILTERED\t";
				#print "segregation_score_best_1: $segregation_score_best_1\n\n";

				#$answer = <STDIN>;

				print OUT_FILTERED "\t$segregation_score_best_1";
				print OUT_FILTERED "\t$segregation_score_best_2";
				print OUT_FILTERED "\t$segregation_score_best_3";
				print OUT_FILTERED "\t$segregation_score_best_4";
				print OUT_FILTERED "\t$segregation_score_best_5";
				print OUT_FILTERED "\t$segregation_score_best_6";
				print OUT_FILTERED "\t$segregation_score_best_7";
			
				print OUT_FILTERED "\t$segregation_score_best_overall"; # Best overall




				#####################################################################
				# Print "R" if the main affected allele is the same as the REF base #
				#####################################################################
				if ($main_affected_allele ne $REF_base){print OUT_FILTERED "\tN"} else {print OUT_FILTERED "\tR";}
				

				#########################################
				# Record whether SNP is on the SNP list #
				#########################################
				if ($on_snp_list eq "false")
				{ print OUT_FILTERED "\tN"; }
				else
				{ print OUT_FILTERED "\t$on_snp_list"; }


				#####################################################################
				# Print the maximum effect score                                    #
				#####################################################################
				print OUT_FILTERED "\t$max_effect_score"; # Effect score
				

				#####################################################################
				# Print the results of VEP or snpEff  (whichever is used)           #
				#####################################################################
				if ($effect_predictor eq "VEP") 
				{
					print OUT_FILTERED "\t$vep_results_string"; # Consequence
					print OUT_FILTERED "\t$vep_gene_string"; # Gene
				} 
				elsif ($effect_predictor eq "snpEff") 
				{
					print OUT_FILTERED "\t$snpEff_results_string"; # Consequence
					print OUT_FILTERED "\t$snpEff_gene_string"; # Gene
				} 

				

				#####################################################################
				# Print the LOF and NMD from snpEff or the SIFT results from VEP    #
				#####################################################################
				if ($effect_predictor eq "VEP")
				{
					print OUT_FILTERED "\t$vep_sift_prediction_string";
					print OUT_FILTERED "\t";	
				}
				elsif ($effect_predictor eq "snpEff")
				{
					print OUT_FILTERED "\t$LOF_gene_string";
					print OUT_FILTERED "\t$NMD_gene_string";						
				}
				
				
				
				#####################################################################
				# Put extra columns in here with base and amino acid change         #
				# details (but if the user specifies this)                          #
				# If not specified then the column are omitted altogether)          #
				#####################################################################
				if ($include_base_change_info eq "true")
				{ 
					if ($effect_predictor eq "VEP")
					{
						print OUT_FILTERED "\t$vep_base_change_string"; # Shows which bases were changed
						print OUT_FILTERED "\t$vep_amino_acid_change_string"; # Shows which amino acids were changed
						print OUT_FILTERED "\t$vep_CDS_position_string"; # Shows base position in CDS (coding sequence)
						print OUT_FILTERED "\t$vep_protein_position_string"; # Shows amino acid psoition in protein
					}
					elsif ($effect_predictor eq "snpEff")
					{
						print OUT_FILTERED "\t$snpEff_base_change_string"; # Shows which bases were changed
						print OUT_FILTERED "\t$snpEff_amino_acid_change_string"; # Shows which amino acids were changed
						print OUT_FILTERED "\t$snpEff_CDS_position_string"; # Shows base position in CDS (coding sequence)
						print OUT_FILTERED "\t$snpEff_protein_position_string"; # Shows amino acid psoition in protein
					}
				}


				print OUT_FILTERED "\t$depth_coverage_from_info";
				print OUT_FILTERED "\t$call_frequency_at_this_position";
				print OUT_FILTERED "\t$no_of_alt_alleles";
				print OUT_FILTERED "\t$ALT_base_overall";
				
				print OUT_FILTERED "\t$homozygosity_ratio";
				print OUT_FILTERED "\t$variant_caller_from_info";
				print OUT_FILTERED "\t$variant_type";
				
				if ($variant_type eq "indel"){print OUT_FILTERED "\t$indel_hash{$main_affected_allele}"} else {print OUT_FILTERED "\t$main_affected_allele"}


				####################################################################
				# These extra columns can be included in the output for de-bugging #
				####################################################################
				if ($include_extra_columns_to_debug eq "true")
				{
					print OUT_FILTERED "\t$main_affected_homozygous_genotype";
					print OUT_FILTERED "\t$AA_count_affected";
					print OUT_FILTERED "\t$AB_count_affected";
					print OUT_FILTERED "\t$BB_count_affected";
					print OUT_FILTERED "\t$AA_AB_count_affected";
					print OUT_FILTERED "\t$AB_BB_count_affected";
					print OUT_FILTERED "\t$HH_count_affected";
					
					print OUT_FILTERED "\t$AA_count_normal";
					print OUT_FILTERED "\t$AB_count_normal";
					print OUT_FILTERED "\t$BB_count_normal";
					print OUT_FILTERED "\t$AA_AB_count_normal";
					print OUT_FILTERED "\t$AB_BB_count_normal";
					print OUT_FILTERED "\t$HH_count_normal";

					print OUT_FILTERED "\t$affected_homozygosity";
					print OUT_FILTERED "\t$normal_homozygosity";
				}

				
				#############################
				# Now the final end of line #
				#############################
				print OUT_FILTERED "\n";
				
			} # if write_to_filtered eq "true"

        } #End If 'If array_size_1 > 8
 
	}
	# if passed_header_lines = true (VCF file)


	#########################################################
    # Has the start of the CHR data been reached yet?       #
	# The line before the actual data starts with '#CHROM'  #
	# This part of the code checks that line, which has the #
	# (very important) headers for the sample columns.      #
    #########################################################
	
	elsif (index($vcf_line,"#CHROM") > -1)  # New elsif to speed things up
	{
		$passed_header_lines = "true";
		
		############################################
		# Report which variant caller was detected #
		############################################
		&print_message("Variant caller check...","message");
		print "Variant caller detected from header lines: $variant_caller\n\n";
		print "\n\t>>Press RETURN to continue ";   $answer = <STDIN>;


		############################################################
		# Parsing of #CHROM data line to get a list of input files #
		############################################################
		@chrom_line_array = split(/\s+/,$vcf_line);
		$chrom_line_array_size = (scalar @chrom_line_array) - 1;
		$no_of_samples_chrom_line = $chrom_line_array_size - 8; # ignore first columns 0-9
		
		for ($array_count = 9; $array_count <= $chrom_line_array_size; $array_count++)
		{
			# Sample names go into sample_name_array_input_order
			$sample_name_array_input_order[$array_count - 8] = $chrom_line_array[$array_count];
			$file_count = $array_count - 8;
		}
		
		
        ############################################################
        # Total number of samples is defined by the number of      #
        # samples on the #CHROM line, $no_of_samples_chrom_line    #
        ############################################################
		$no_of_samples = $no_of_samples_chrom_line;


		############################################################
		# Show user a list of the files                            #
		# (and get maximum sample name length for later padding)   #
		############################################################
		&print_message("Sample columns found in the VCF file (from the '#CHROM' line)","message");
		$max_sample_name_length = 0;
		
		for ($sample_count = 1; $sample_count <= $no_of_samples; $sample_count++)
		{
			print "$sample_count \t$sample_name_array_input_order[$sample_count]\n";
			
			if (length($sample_name_array_input_order[$sample_count]) > $max_sample_name_length)
			{
				$max_sample_name_length = length($sample_name_array_input_order[$sample_count]);
			}
		}
		
		print "\nThere are $no_of_samples samples in this multi-column VCF file.\n\n";
		print "\t>>If these columns look OK, press 'return' to continue .\n\n";
		
        $answer = <STDIN>;
		
		
		############################################################
		# Get the disease statuses from a file                     #
		############################################################
		&print_message("Please input the name of the file with the disease statuses of the samples","input");
		
		until (-e $disease_status_file)
		{
			print "   Disease status file: ";
			$disease_status_file = <STDIN>;
			chomp $disease_status_file;

			if ($disease_status_file eq ""){$disease_status_file = "ls"}
			if ($disease_status_file eq "ls"){print "\n";system ("ls *.txt");print "\n"}

			# Starts with 'ls' followed by search string
			if (($disease_status_file ne "ls") && (index ($disease_status_file,"ls") == 0))
			{
				$search_string = substr($disease_status_file,3,99);
				print "\n";
				system ("ls *$search_string*");
				$disease_status_file = "ls";
				print "\n";
			}

			#DEBUG
			if (! -e $disease_status_file){print "File $disease_status_file does not exist\n";}
			if (-e $disease_status_file){print "File $disease_status_file exists\n";}

			if (($disease_status_file ne "ls") && ($disease_status_file ne "help")){if (! -e $disease_status_file){print "\n\n>>>>>>>>  File $disease_status_file not found.  Try again.  <<<<<<<<\n\n";}}

			if ($disease_status_file eq "help")
			{
				&print_message("HELP ON DISEASE STATUS FILE","help");
				print "\tThe disease status file has two columns, separated by a TAB.\n\n";
				print "\tThe first column has the name of the sample, as used as the column header in the VCF file.\n";
				print "\tThe second column must be 'Affected', 'Carrier' or 'Normal'\n";
				print "\t------------------------------------------------------------------------------------------\n\n";
			}
		} #until (-e $disease_status_file)


		############################################################
		# Make sure the list file is in Unix format                #
		############################################################
		print "\n";
		$command = "dos2unix $disease_status_file";
		system("$command");
		print "\n";


		############################################################
		# Reads the disease status file to get which samples are   #
		# cases and which are controls.                            #
		############################################################
		&read_disease_status_file;



		##########################################################################
		# Special case when there is only one case and many controls             #
		# Set $only_one_case as "true" and change the segregation scoring system #
		# so that to get a score, as well as matching the correct genotype       #
		# expected, you also have to be DIFFERENT FROM THE SINGLE CASE           #
		# (The single case gets a point for matching its required genotype too)  #
		##########################################################################
		

		if (($affected_count_ds_file == 1) && ($normal_count_ds_file > 1) && ($single_case_scoring eq "true"))
		{
			&print_message("Single case segregation scores selected","message");

			print "Number of affected samples:    \t$affected_count_ds_file\n";
			print "Number of unaffected samples:  \t$normal_count_ds_file\n\n";

			print "Because there is a single case and multiple controls, a special\n";
			print "segregation scoring system will be used.  Controls will only\n";
			print "score if their genotype is DIFFERENT to the single case\n\n";

			&pause;
			$only_one_case = "true";
		}


		############################################################
		# Get run title to add to name of output file              #
		############################################################
		&get_name_of_output_files;


		############################################################
		# Open command log file                                    #
		############################################################
		open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
		print COMMAND_LOG "Command log of vcf2excel, version $version\n\n";
		print COMMAND_LOG "VCF file:            \t$vcf_file\n";
		print COMMAND_LOG "Disease status file: \t$disease_status_file\n";
		print COMMAND_LOG "Reference genome:    \t$ref_seq_name\n";

		
		####################################################################
		# Open file to record effects which have not been assigned a score #
		####################################################################
		open (MISSING_EFFECTS, ">$missing_effects")|| die "Cannot create output file: $missing_effects";


		#############################################################################################
		#  Start to 'Re-map' the columns so that the order is Affected -> Carrier -> Normal         #
		#                                                                                           #
		# This is then used by starting with the new order and looking up a "source column"         #
		# from the VCF file to get the data from.                                                   #
		#                                                                                           #
		# So if your second affected sample was in column 5 of the VCF file then the new order      #  
		# count would be 2 and the $source_col would be 5.  Data would be transferred from column 5 #   
		# and put into new column 2.                                                                #
		#############################################################################################
		$new_sample_count = 0;

		
		########################
		# Get affected samples #
		########################		
		for ($disease_status_count = 1;$disease_status_count <=$no_samples_in_status_file; $disease_status_count++)
		{
			if ($disease_status_array[$disease_status_count] eq "affected")
			{
				$no_affected_samples = $no_affected_samples + 1;
				
				# Find this name in VCF file array (sample_name_array)
				for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
				{
					if ($disease_status_sample_array[$disease_status_count] eq $sample_name_array_input_order[$sample_count])
					{
						$new_sample_count = $new_sample_count + 1;

						$source_column_array[$new_sample_count] = $sample_count; # Stores array of source columns
						$sample_status_array_input_order[$sample_count] = "affected";
					} # if names match
				}
			} # If affected
		} # Affecteds
		

		#######################
		# Get carrier samples #
		#######################		
		for ($disease_status_count = 1;$disease_status_count <=$no_samples_in_status_file; $disease_status_count++)
		{
			if ($disease_status_array[$disease_status_count] eq "carrier")
			{
				$no_carrier_samples = $no_carrier_samples + 1;
				
				# Find this name in VCF file array (sample_name_array)
				for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
				{
					if ($disease_status_sample_array[$disease_status_count] eq $sample_name_array_input_order[$sample_count])
					{
						$new_sample_count = $new_sample_count + 1;

						$source_column_array[$new_sample_count] = $sample_count; # Stores array of source columns
						$sample_status_array_input_order[$sample_count] = "carrier";
					} # if names match
				}
			} # If carrier
		} # Carriers
		

		######################
		# Get Normal samples #
		######################		
		for ($disease_status_count = 1;$disease_status_count <=$no_samples_in_status_file; $disease_status_count++)
		{
			if ($disease_status_array[$disease_status_count] eq "normal")
			{
				$no_normal_samples = $no_normal_samples + 1;
				
				# Find this name in VCF file array (sample_name_array)
				for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
				{
					if ($disease_status_sample_array[$disease_status_count] eq $sample_name_array_input_order[$sample_count])
					{
						$new_sample_count = $new_sample_count + 1;

						$source_column_array[$new_sample_count] = $sample_count; # Stores array of source columns
						$sample_status_array_input_order[$sample_count] = "normal";			
					} # if names match
				}
			} # If normal
		} # Normals


		#######################
		# Get Omitted samples #
		#######################	
		for ($disease_status_count = 1;$disease_status_count <=$no_samples_in_status_file; $disease_status_count++)
		{
			if ($disease_status_array[$disease_status_count] eq "omit")
			{
				$no_omit_samples = $no_omit_samples + 1;
				
				# Find this name in VCF file array (sample_name_array)
				for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
				{
					if ($disease_status_sample_array[$disease_status_count] eq $sample_name_array_input_order[$sample_count])
					{
						$new_sample_count = $new_sample_count + 1;

						$source_column_array[$new_sample_count] = $sample_count; # Stores array of source columns
						$sample_status_array_input_order[$sample_count] = "omit";			
					} # if names match
				}
			} # If omit
		} # Omit

	

		##########################################################################################################
		# Re-map sample names and sample statuses to create $sample_name_array in output or disease status order #
		##########################################################################################################
		for($sample_count = 1; $sample_count <= $no_of_samples; $sample_count++)
		{
            #Get the source column
            $source_col = $source_column_array[$sample_count];

            #Fill two new arrays in output order
            $sample_name_array[$sample_count] = $sample_name_array_input_order[$source_col];

            $sample_status_array[$sample_count] = $sample_status_array_input_order[$source_col];
        }


		################################################################
		# Show new column positions                                    #
		################################################################	
		&print_message("New destination columns for samples (so they are Affected -> Carrier -> Normal)","message");

		&print_both ("\n\n\n\tOriginal arrays (this is the order in the VCF file):");
		&print_both ("\n\n\tOld column\tSample\tStatus\t\tSource\n");
		&print_both ("\t----------\t------\t------\t\t------\n");

		for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
		{
			# pad out disease status to 12 characters
			$disease_status = sprintf("%-*s", 12, $sample_status_array_input_order[$sample_count]);

			
			# pad out sample name to max_sample_name_length characters
			$sample_name = sprintf("%-*s", $max_sample_name_length + 2, $sample_name_array_input_order[$sample_count]);
			
			&print_both("\t$sample_count\t\t$sample_name\t\t$disease_status\t\t$sample_count\n");
		}

		&print_both ("\n\n\tRe-mapped arrays (this is the order for the Excel file):");
		&print_both ("\n\n\tNew column\tSample\tStatus\t\tSource\n");
		&print_both ("\t----------\t------\t------\t\t------\n");

		$omit_count = 0;

		for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
		{
			# pad out disease status to 12 characters			
			$disease_status = sprintf("%-*s", 12, $sample_status_array[$sample_count]);
			
			# pad out sample name to $max_sample_name_length + 2 characters
			$sample_name = sprintf("%-*s", $max_sample_name_length + 2, $sample_name_array[$sample_count]);
			
			#Get the source column
            $source_col = $source_column_array[$sample_count];

			&print_both ("\t$sample_count\t\t$sample_name\t\t$disease_status\t\t$source_col\n");


			############################################
			# Count how many samples are to be omitted #
			############################################
			if ($sample_status_array[$sample_count] eq "omit")
			{
				$omit_count = $omit_count + 1;
			}
		} # sample_count loop
		

		################################################################
		# Recalculate the number of samples                            #
		################################################################
		$no_of_samples_without_omits = $no_of_samples - $omit_count;

		if ($omit_count == 1)
		{
			print "\nNo of samples to be used: $no_of_samples_without_omits  ($omit_count sample omitted)\n\n"
		}
		if ($omit_count > 1)
		{
			print "\nNo of samples to be used: $no_of_samples_without_omits  ($omit_count samples omitted)\n\n"
		}

		print "\n\t>>PLEASE CHECK THESE CAREFULLY.  If they look correct, press 'RETURN'\n";
		$answer=<STDIN>;
		

		################################################################
		# Ask about what the segregation score for filtering should be #
		################################################################
		&print_message("Choose minimum segregation score for filtering","input");

		$minimum_segregation_score = $no_of_samples_without_omits - 1;
		$max_seg_score = $no_of_samples_without_omits;
		$second_seg_score = $no_of_samples_without_omits - 1;
		$third_seg_score = $no_of_samples_without_omits - 2;

		print "There are $no_of_samples_without_omits samples.";
		print "\tThe maximum segregation score is therefore $no_of_samples_without_omits\n\n";

		print "Enter the value for the minimum segregation score allowed into the filtered output file? (DEFAULT = $second_seg_score)\n\n";
		print "(i.e. segregation scores below this are not written to the filtered file):  ";

		$answer=<STDIN>;
		chomp $answer;

		if  ($answer ne "")
		{
			if (($answer <=$no_of_samples_without_omits) && ($answer > 0)){$minimum_segregation_score = $answer}
		}
		if  ($answer eq "")
		{
			$minimum_segregation_score = $second_seg_score;
		}

		print "\n\nMinimum segregation score = $minimum_segregation_score\n\n";

	
		
		################################################################
		# Ask about what the effect score for filtering should be      #
		################################################################
		&print_message("Choose minimum effect score for the filtered file","input");

		$minimum_effect_score = 3;

		print "Effect scores below this minimum are NOT written to the filtered output file\n\n";
		print "The default is to keep all SNPs with an effect score of 2 or greater.";
		print "You can increase this if you want to make the output filtered file smaller.\n\n";
		print "The maximum possible effect score is 5 so the highest minimum effect score is 5.\n\n";

		print "Enter the value for the minimum effect score allowed into the filtered output file? (DEFAULT = $minimum_effect_score)\n\n";
		
		$answer=<STDIN>;
		chomp $answer;

		if  ($answer ne "")
		{
			if (($answer <=5) && ($answer > 0)){$minimum_effect_score = $answer}
		}
		if  ($answer eq "")
		{
			$minimum_effect_score = $minimum_effect_score; #default
		}

		print "\n\nMinimum effect score = $minimum_effect_score\n\n";



		################################################################
		# Ask what you want to do with unscored genotypes ('XX')       #
		################################################################
		&print_message("How do you want score ungenotyped variants ('XX')?","input");

		print "There are three methods to do this.\n\n";
		
		print "The old method was to score them all as the Reference allele.\n\n";
		print "This is good because that is the allele they are most likely to be, but is bad because\n";
		print "you might miss a really good segregation because of one un-genotyped variant in one sample.\n\n";
		print "I suggest you start with method 1, and then maybe try the others as an alternative if you like.\n\n";

		print "   <1> Score 'XX' genotypes as two Reference alleles.                                [DEFAULT]\n";
		print "   <2> Score 'XX' genotypes as whatever will give the best segregation score.   \n";
		print "   <3> Don't score 'XX' genotypes at all, but then normalise the segregation score to compensate for the missing data.\n\n";


		$answer=<STDIN>; chomp $answer;

		# Default is 1
		if ($answer eq ""){$answer = "1"}

		if ($answer eq "1" ){$score_X_method = "reference"}
		if ($answer eq "2" ){$score_X_method = "best"}
		if ($answer eq "3" ){$score_X_method = "missing"}



		################################################################
		# Ask how you want to combine these filters. AND or OR?        #
		################################################################
		if ($combine_filters eq "")
		{
			&print_message("How do you want to combine the segregation and effect score filters?","input");

			print "You can combine these using AND or using OR.\n\n";
			
			print "If you use 'AND' you will get a much smaller output file,\n";
			print "which is probably the best way to handle Whole Genome Sequencing data.\n\n";
			print "Note: with 'AND' you may miss out on variants that segregate perfectly but have a low effect score.\n\n";

			print "   <1> Combine with 'AND'. Filtered SNPs must have a seg score >= $minimum_segregation_score AND an effect score >= $minimum_effect_score   [DEFAULT]\n";
			print "   <2> Combine with 'OR'.  Filtered SNPs must have a seg score >= $minimum_segregation_score OR an effect score >= $minimum_effect_score\n\n";


			$answer=<STDIN>; chomp $answer;

			# Default is 1
			if ($answer eq ""){$answer = "1"}

			if (substr($answer,0,1) eq "2" ){$combine_filters = "OR"} else {$combine_filters = "AND"}
		}
				
		
		################################################################################
		# Ask if you want to keep the LARGE unfiltered file                            #
		################################################################################
		if ($keep_unfiltered_output ne "false")
		{
			&print_message("Do you want to keep the large unfiltered output file?","input");

			print "These are usually too large to be read into Excel anyway!\n\n";

			print "   <1>  No -  delete the unfiltered output file                        [DEFAULT]\n";
			print "   <2>  Yes - keep the unfiltered and the filtered output files\n";
			
			$answer=<STDIN>; chomp $answer;
			if ($answer eq ""){$answer = "1"} # default

			if (substr($answer,0,1) eq "1" ){$keep_unfiltered_output = "false"};
			if (substr($answer,0,1) eq "2" ){$keep_unfiltered_output = "true"} 

		} # if ($keep_unfiltered_output ne "false")


		################################################################################
		# Ask if you want to make a merged indels file                                 #
		# (if this is set to false in the declarations at the top then don't even ask) #
		################################################################################

		if ($make_merged_indels_file ne "false")
		{
			&print_message("Do you want to make a merged indels output file?","input");

			print "This is an experimental analysis, in which very close indels are merged together.\n";
			print "The idea is to overcome situations where indels are treated as different variants.\n";
			print "just because of a small alignment error.\n\n";

			print "   <1>  No  - don't make a merged indels file                     [DEFAULT]\n";
			print "   <2>  Yes - make a merged indels file\n";
			
			$answer=<STDIN>; chomp $answer;
			if ($answer eq ""){$answer = "1"} # default

			if (substr($answer,0,1) eq "1" ){$make_merged_indels_file = "false"};
			if (substr($answer,0,1) eq "2" ){$make_merged_indels_file = "true"} 
		}


		################################################################################
		# Ask if you want to keep the indels_only output file                          #
		# (if this is set to false in the declarations at the top then don't even ask) #
		################################################################################
		if ($keep_indels_only_file ne "false")
		{
			&print_message("Do you want to keep the 'indels_only' output file?","input");

			print "The filtered output contains indels as well as SNPs so this file is not required\n";
			print "unless you want an Indels Only file for some reason\n\n";

			print "   <1>  No  - delete the 'indels_only' output file                     [DEFAULT]\n";
			print "   <2>  Yes - keep the 'indels_only' output file\n";
			
			$answer=<STDIN>; chomp $answer;
			if ($answer eq ""){$answer = "1"} # default

			if (substr($answer,0,1) eq "1" ){$keep_indels_only_file = "false"};
			if (substr($answer,0,1) eq "2" ){$keep_indels_only_file = "true"} 
		}


		################################################################################
		# Ask if you want to filter on call frequency                                  #
		# (this is only offered if you are using one of the two new score_X_method     #
		# options, i.e. @best' and 'missing'.  These can allow a large number of       #
		# low call frequency variants into the filetered file, so the user is          #
		# offered the chance to add an extra filter here                               #
		################################################################################
		if ($score_X_method ne "reference")
		{
			&print_message("Do you want to filter off variants with low Call Frequency?","input");

			print "Call frequency is the percentage of samples with valid genotypes (current threshold: $call_frequency_threshold%) \n\n";
			print "so if 90% of the genotypes are scored as 'XX' then the Call Frequency is 10%\n\n";

			print "   <1>  Yes - omit variants with low call frequencies                   [DEFAULT]\n";
			print "   <2>  No -  keep samples whatever the call frequency\n";
			
			$answer=<STDIN>; chomp $answer;
			if ($answer eq ""){$answer = "1"} # default

			if (substr($answer,0,1) eq "1" ){$omit_low_call_frequencies = "true"};
			if (substr($answer,0,1) eq "2" ){$omit_low_call_frequencies = "false"} 

			if ($omit_low_call_frequencies eq "true")
			{
				&print_message("What cut-off frequency would you like to use?","input");

				print "If the Call Frequency cut-off is 20% then any variants with less than 20% Call Frequency are omitted\n\n";
				print "The default here is $call_frequency_threshold\n\n";

				print "  Enter a number between 1-100  >  ";

				$answer=<STDIN>; chomp $answer;
				if ($answer ne ""){$call_frequency_threshold = $answer} # default
			}

		} # if ($score_X_method ne "reference")
		


		######################################################################################
		# START CONVERTING VCF FILE                                                          #
		######################################################################################
		&print_message("Converting VCF file to text file for Excel sheet NGS SNP Viewer...","message");


		################
		# START TIMERS #
		################
		$start_time = time();


		############################################################
		# Set various things to zero                               #
		############################################################
		for ($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
		{
			$missing_genotypes_array[$sample_count] = 0;
			$genotypes_counted_array[$sample_count] = 0;
		}
		
        ############################################################
        # Open the output files                                    #
        ############################################################
        # Out for all variants (unfiltered)
		if ($keep_unfiltered_output eq "true")
        {
			open (OUT, ">$output_file")|| die "Cannot create output file: $output_file";
		}

		# Output for selected (filtered) variants
		open (OUT_FILTERED, ">$output_file_filtered")|| die "Cannot create output file: $output_file_filtered";
		
		# Output for indels
		if ($make_merged_indels_file eq "true")
		{
			open (INDELS_OUTPUT, ">$simplified_indels_output_file")|| die "Cannot create output file: $simplified_indels_output_file";
		}




		###################################################################
		# Write headers for the first 3 columns in the output file        #
		###################################################################
		if ($keep_unfiltered_output eq "true") { print OUT "Chr\tPos\tRef"; }

		print OUT_FILTERED "Chr\tPos\tRef";
		
		if ($make_merged_indels_file eq "true"){ print INDELS_OUTPUT "Chr\tPos\tRef"; }
		

        ###################################################################
        # Add name of sample as header in first batch of columns          #
        # (ie those after the first three CHR,POS,REF columns             #
    	#                                                                 #
    	# Add "_A" to Affected samples                                    #
        # Add "_C" to Carrier samples                                     #
        # Add "_N" to Normal samples                                      #
        #                                                                 #
        # This means the disease statuses are hard coded into the output  #
        # file, rather than just in the file title                        #
        # Also Excel NGS SNP VIEWER won't need the affected samples to    #
        # be coloured red.                                                #
        ###################################################################

        ####################
        # Affected samples #
        ####################
        for ($sample_count = 1;$sample_count <=$no_of_samples_without_omits; $sample_count++)
		{
            $sample_name = &get_prefix($sample_name_array[$sample_count]);
			
			#########################
			# Affected suffix added #
			#########################
			if ($sample_count < $no_affected_samples + 1){ $sample_name = $sample_name."_A"; }

			#########################
			# Carrier suffix added  #
			#########################
			if ($no_carrier_samples > 0)
			{
				if (($sample_count > $no_affected_samples) && ($sample_count <= ($no_carrier_samples + $no_affected_samples)))
			    {
			    	$sample_name = $sample_name."_C";
			    }
			}	

			#########################
			# Normal suffix added   #
			#########################
			if ($sample_count > ($no_affected_samples + $no_carrier_samples)){ $sample_name = $sample_name."_N"; }

			if ($keep_unfiltered_output eq "true") 
			{
	            print OUT "\t$sample_name"; # Once for allele_1 
				print OUT "\t$sample_name"; # Once for allele_2
			}

			print OUT_FILTERED "\t$sample_name"; # Once for allele_1 
			print OUT_FILTERED "\t$sample_name"; # Once for allele_2

			if ($make_merged_indels_file eq "true")
			{
				print INDELS_OUTPUT "\t$sample_name"; # Once for allele_1 
				print INDELS_OUTPUT "\t$sample_name"; # Once for allele_2
			}

		} # Sample_count loop to add sample names as headers


		#########################################################################
		# Print headers for the group of columns to the right of the alleles    #
		# This is used by Excel                                                 #
		#########################################################################

		# 1 AA_AB_BB         	Recessive:  								Affecteds AA, carriers AB, normals BB
		# 2 AA_AB_ABorBB     	Recessive, some normals may be carriers:  	Affecteds AA, carriers AB, normals AB or BB
		# 3 AA_ABorBB_BB      	Recessive, some carriers may be normals:  	Affecteds AA, carriers AB or BB, normals BB
		# 4 AA_ABorBB_ABorBB 	Recessive:  								Affecteds AA, carriers AB or BB, normals AB or BB
		# 5 AAorAB_BB        	Dominant:   								Affecteds AA or AB, normals BB
		# 6 A_B              	Additive:   								Affecteds A with As, Normals with Bs
		# 7 AA_BB               Difference:                                 Normals different to main affected homozygous genotypeB
		
		#######################
		# Unfiltered OUT file                                                    #
		# Don't write to this file unless the unfiltered output file is required #
		##########################################################################
		if ($keep_unfiltered_output eq "true")
		{
			print OUT "\tSeg score 1 AA/AB/BB";
			print OUT "\tSeg score 2 AA/AB/ABorBB";
			print OUT "\tSeg score 3 AA/ABorBB/BB";
			print OUT "\tSeg score 4 AA/ABorBB/ABorBB";
			print OUT "\tSeg score 5 AAorAB/BB dom";
			print OUT "\tSeg score 6 A/B add";
			print OUT "\tSeg score 7 Aff <> Norm";
			print OUT "\tSeg score BEST";
			print OUT "\tMain affected allele is REF";
			print OUT "\tOn SNP list?";
			print OUT "\tEffect score";
			print OUT "\tConsequence";
			print OUT "\tGene";
			
			if ($effect_predictor eq "snpEff"){print OUT "\tLOF"} else {print OUT "\tSIFT"}

			print OUT "\tNMD";
					

			#####################################################################
			# Put extra columns in here with base and amino acid change         #
			# details (but if the user specifies this)                          #
			# If not specified then the columns are omitted altogether)         #
			#####################################################################
			if ($include_base_change_info eq "true")
			{
				print OUT "\tBase change"; # Shows which bases were changed
				print OUT "\tAmino acid change"; # Shows which amino acids were changed
				print OUT "\tCDS_pos"; # Shows base position in CDS (coding sequence)
				print OUT "\tProtein_pos"; # Shows amino acid psoition in protein
			}

			print OUT "\tDepth";
			print OUT "\tCall_Freq";
			print OUT "\tNo of ALT alleles";
			print OUT "\tAlt";

			print OUT "\tHomozygosity Ratio";
			print OUT "\tVariant Caller";
			print OUT "\tVariant type";
			print OUT "\tMain_affected_allele";

		} # if ($keep_unfiltered_output eq "true")
		

		if ($make_merged_indels_file eq "true")
		{
			#######################
			# INDELS OUT file     #
			#######################
			print INDELS_OUTPUT "\tSeg score 1 AA/AB/BB";
			print INDELS_OUTPUT "\tSeg score 2 AA/AB/ABorBB";
			print INDELS_OUTPUT "\tSeg score 3 AA/ABorBB/BB";
			print INDELS_OUTPUT "\tSeg score 4 AA/ABorBB/ABorBB";
			print INDELS_OUTPUT "\tSeg score 5 AAorAB/BB dom";
			print INDELS_OUTPUT "\tSeg score 6 A/B add";
			print INDELS_OUTPUT "\tSeg score 7 Aff <> Normt";
			print INDELS_OUTPUT "\tSeg score BEST";
			print INDELS_OUTPUT "\tMain affected allele is REF";
			print INDELS_OUTPUT "\tEffect score";
			print INDELS_OUTPUT "\tConsequence";
			print INDELS_OUTPUT "\tGene";

			print INDELS_OUTPUT "\tAlt";
			print INDELS_OUTPUT "\tallele_A";
			print INDELS_OUTPUT "\tallele_B";
			print INDELS_OUTPUT "\tMain_affected_allele";

			print INDELS_OUTPUT "\tHomozygosity Ratio";
			print INDELS_OUTPUT "\tNo of ALT alleles";
			print INDELS_OUTPUT "\tVariant Caller";
			print INDELS_OUTPUT "\tDepth";
			print INDELS_OUTPUT "\tVariant type";
			print INDELS_OUTPUT "\tLOF";
			print INDELS_OUTPUT "\tNMD";
		} #if ($make_merged_indels_file eq "true")
		

		#Filtered output
		print OUT_FILTERED "\tSeg score 1 AA/AB/BB";
		print OUT_FILTERED "\tSeg score 2 AA/AB/ABorBB";
		print OUT_FILTERED "\tSeg score 3 AA/ABorBB/BB";
		print OUT_FILTERED "\tSeg score 4 AA/ABorBB/ABorBB";
		print OUT_FILTERED "\tSeg score 5 AAorAB/BB dom";
		print OUT_FILTERED "\tSeg score 6 A/B add";
		print OUT_FILTERED "\tSeg score 7 Aff <> Norm";
		print OUT_FILTERED "\tSeg score BEST";
		print OUT_FILTERED "\tMain affected allele is REF";
		print OUT_FILTERED "\tOn SNP list?";
		print OUT_FILTERED "\tEffect score";
		print OUT_FILTERED "\tConsequence";
		print OUT_FILTERED "\tGene";
		

		if ($effect_predictor eq "snpEff"){print OUT_FILTERED "\tLOF"} else {print OUT_FILTERED "\tSIFT"}

		print OUT_FILTERED "\tNMD";

		#####################################################################
		# Put extra columns in here with base and amino acid change         #
		# details (but if the user specifies this)                          #
		# If not specified then the column are omitted altogether)          #
		#####################################################################
		if ($include_base_change_info eq "true")
		{
			print OUT_FILTERED "\tBase change"; # Shows which bases were changed
			print OUT_FILTERED "\tAmino acid change"; # Shows which amino acids were changed
			print OUT_FILTERED "\tCDS_pos"; # Shows base position in CDS (coding sequence)
			print OUT_FILTERED "\tProtein_pos"; # Shows amino acid psoition in protein
		}

		print OUT_FILTERED "\tDepth";
		print OUT_FILTERED "\tCall_Freq";
		print OUT_FILTERED "\tNo of ALT alleles";
		print OUT_FILTERED "\tAlt";

		print OUT_FILTERED "\tHomozygosity Ratio";
		print OUT_FILTERED "\tVariant Caller";
		print OUT_FILTERED "\tVariant type";
		print OUT_FILTERED "\tMain_affected_allele";
		
		
		####################################################################
		# These extra columns can be included in the output for de-bugging #
		####################################################################
		if ($include_extra_columns_to_debug eq "true")
		{
			print OUT_FILTERED "\tMain_affected_homozygous_genotype";
			print OUT_FILTERED "\tAA_count_affected";
			print OUT_FILTERED "\tAB_count_affected";
			print OUT_FILTERED "\tBB_count_affected";
			print OUT_FILTERED "\tAA_AB_count_affected";
			print OUT_FILTERED "\tAB_BB_count_affected";
			print OUT_FILTERED "\tHH_count_affected";
			print OUT_FILTERED "\tAA_count_normal";
			print OUT_FILTERED "\tAB_count_normal";
			print OUT_FILTERED "\tBB_count_normal";
			print OUT_FILTERED "\tAA_AB_count_normal";
			print OUT_FILTERED "\tAB_BB_count_normal";
			print OUT_FILTERED "\tHH_count_normal";
			print OUT_FILTERED "\tAffected Homozygosity";
			print OUT_FILTERED "\tNormal Homozygosity";
		} # extra columns


		if ($keep_unfiltered_output eq "true")
		{
			####################################################################
			# These extra columns can be included in the output for de-bugging #
			####################################################################
			if ($include_extra_columns_to_debug eq "true")
			{
				print OUT "\tMain_affected_homozygous_genotype";
				print OUT "\tAA_count_affected";
				print OUT "\tAB_count_affected";
				print OUT "\tBB_count_affected";
				print OUT "\tAA_AB_count_affected";
				print OUT "\tAB_BB_count_affected";
				print OUT "\tHH_count_affected";
				print OUT "\tAA_count_normal";
				print OUT "\tAB_count_normal";
				print OUT "\tBB_count_normal";
				print OUT "\tAA_AB_count_normal";
				print OUT "\tAB_BB_count_normal";
				print OUT "\tHH_count_normal";
				print OUT "\tAffected Homozygosity";
				print OUT "\tNormal Homozygosity";
			} # extra columns

		} # if ($keep_unfiltered_output eq "true")


		#############################
		# Now the final end of line #
		#############################
		if ($keep_unfiltered_output eq "true"){ print OUT "\n"; }
		if ($make_merged_indels_file eq "true"){ print INDELS_OUTPUT "\n"; }

		print OUT_FILTERED "\n";
		

	} # End of: if (index($vcf_line,"#CHROM") > -1)

	if (($line_count % 20000) == 0)
	{
		if ($passed_header_lines eq "true")
		{
			print "Reading VCF file.  Line: $line_count\tChromosome: $chromosome\tNumber in filtered file: $filter_count_all\n";
		}
	}
	 
} # End of: while ($vcf_line = <VCF>)  top of loop is at 587

if ($keep_unfiltered_output eq "true"){close OUT;} # close out file for excel}

close OUT_FILTERED; # close filtered out file for excel

if ($make_merged_indels_file eq "true"){ close INDELS_OUTPUT; }

print "\n\n";


############################################################################
# Merge Indels                                                             #
# This procedure merges nearby indels into one, in case they are actually  #
# the same indels by the variant caller has decided thay are in a slightly #
# different position.                                                      #
############################################################################

if ($make_merged_indels_file eq "true")
{
	&merge_indels;
}

##########################################################################
# Sort output filtered file so that highest effect scores are at the top #
##########################################################################

#  (head -n 1 test.txt; tail -n +3 test.txt | sort -k33,33rn -k31,31rn) > sorted_2.out


print "no_of_samples: $no_of_samples_without_omits\n\n";

$effect_score_col = ($no_of_samples_without_omits * 2) + 14;
$seg_score_best_col = ($no_of_samples_without_omits * 2) + 11;

&print_message("Sorting output file so that highest effect scores are at the top","message");



$command = "(head -n 1 $output_file_filtered; tail -n +2 $output_file_filtered | sort -k$effect_score_col,$effect_score_col"."rn -k$seg_score_best_col,$seg_score_best_col"."rn) > $output_file_filtered_and_sorted";
&print_both("\n");
&print_both("$command");
system("$command");


################################################################
# If sorted file is same size as unsorted then delete unsorted #
################################################################
if (-s $output_file_filtered == -s $output_file_filtered_and_sorted)
{
	$command = "rm $output_file_filtered";
	&print_both("$command");
	system("$command");
}


################################################################
# If user has chosen to remove large unfiltered output file    #
################################################################
if ($keep_unfiltered_output eq "false")
{
	$command = "rm $output_file";
	&print_both("$command");
	system("$command");
}


################################################################
# If user has chosen to remove the 'indels_only' output file   #
# (or hasn't opted to make the merged indels file, which the   #
# indels only file is needed for)                              #
################################################################
if (($keep_indels_only_file eq "false")  && ($make_merged_indels_file eq "true"))
{
	$command = "rm $simplified_indels_output_file";
	&print_both("$command");
	system("$command");
}



#############################
# Write details to log file #
#############################
print COMMAND_LOG "\n\n####################\n";
print COMMAND_LOG "# Analysis details #\n";
print COMMAND_LOG "####################\n\n";

&print_message("FINISHED ANALYSIS","message");


#Samples
&print_both ("\tNumber of affecteds: \t$no_affected_samples\n");
&print_both ("\tNumber of carriers:  \t$no_carrier_samples\n");
&print_both ("\tNumber of normals:   \t$no_normal_samples\n\n");

# Missing genotypes
&print_both("\nMissing genotypes:  \n\n");

&print_both("\tSample\tNo missing\t% missing\n");
&print_both("\t======\t==========\t=========\n");

for ($sample_count = 1;$sample_count <=$no_of_samples_without_omits; $sample_count++)
{
		#Calculate precentage missing for each sample
		$missing_genotypes_percent = $missing_genotypes_array[$sample_count]/$genotypes_counted_array[$sample_count] * 100;
		$missing_genotypes_percent = sprintf("%.3f", $missing_genotypes_percent);
		
		# pad out sample name to $max_sample_name_length + 2 characters
		$sample_name = sprintf("%-*s", $max_sample_name_length + 2, $sample_name_array[$sample_count]);
		
		&print_both("\t$sample_name\t$missing_genotypes_array[$sample_count]   \t$missing_genotypes_percent%\n")
}

# Some calculations
$missing_genotypes_percent = $missing_genotypes_total/$genotypes_counted_total * 100;
$missing_genotypes_percent = sprintf("%.3f", $missing_genotypes_percent);
$genotypes_counted_per_sample = $genotypes_counted_total/$no_of_samples;
$genotypes_counted_per_sample = sprintf("%.0f", $genotypes_counted_per_sample);

&print_both("\n\tTotal number of genotypes counted:  \t$genotypes_counted_total\n");
&print_both("\tTotal number of missing genotypes:  \t$missing_genotypes_total\n");
&print_both("\tTotal percentage missing genotypes: \t$missing_genotypes_percent%\n");


&print_both ("\nVariant annotations: \n\n");

&print_both ("\tVariant effect predictor selected: $effect_predictor\n");

if (($vep_data_found eq "false") && ($snpEff_data_found eq "false"))
{
	&print_both ("\tNo variant annotations found\n");
}

if (($effect_predictor eq "snpEff") && ($snpEff_data_found eq "true"))
{
	&print_both ("\tVariant effect predictor used:     $effect_predictor\n");
}
if (($effect_predictor eq "VEP") && ($vep_data_found eq "true"))
{
	&print_both ("\tVariant effect predictor used:     $effect_predictor\n");
}

if (($effect_predictor eq "snpEff") && ($snpEff_data_found eq "false"))
{
	&print_message("WARNING. snpEff data not found","warning");

	&print_both ("\tYou said the file was annotated by snpEff but snpEff data was not found\n");
}
if (($effect_predictor eq "VEP") && ($snpEff_data_found eq "true") && ($vep_data_found eq "true"))
{
	&print_both ("\tYou said you wanted to use VEP annotations but snpEff data was found as well (VEP only was used)\n");
}


&print_both("\tNumber of ensembl gene names replaced by real names: \t$ensembl_gene_replaced_count\n");


#####################################################
# Details of minimum scores kept for filtered file  #
#####################################################
&print_both("\n\nMinimum score thresholds for retaining in filtered file:\n\n");

&print_both("\tMinimum segregation score:                                    \t$minimum_segregation_score\n");
&print_both("\tMinimum effect score:                                         \t$minimum_effect_score\n\n");

&print_both("\tMethod of combining these:                                    \t$combine_filters\n");
&print_both("\tBase change info included:                                    \t$include_base_change_info\n"); # true or false
&print_both("\tMethod of using ungenotyped Xs in scoring:                    \t$score_X_method\n\n"); # reference, best or missing

if(($only_one_case eq "true")  && ($single_case_scoring eq "true")){
&print_both("\tSingle case scoring was used.  This is a variant of the scoring when there is only a single case.\n\n");}


#####################################################
# Details of variants written to filtered file      #
#####################################################
&print_both("\nVariants retained and written to filtered file:\n\n");
&print_both("\tTotal number of variants in VCF file:                         \t$vcf_variant_count\n");

if ($combine_filters eq "OR")
{
	&print_both("\tVariants retained because of segregation score:              \t$filter_count_seg_score\n");
	&print_both("\tVariants retained because of effect score:                   \t$filter_count_effect_score\n");
} # OR

if ($combine_filters eq "AND")
{
	&print_both("\tVariants retained because of segregation AND effect score:    \t$filter_count_combined_score\n");
} # AND

&print_both("\tVariants retained because more than 2 alleles:                \t$filter_count_no_of_alleles\n");

&print_both("\tVariants omitted because they were monomorphic:               \t$monomorphic_count\n");

&print_both("\tVariants omitted because they had low call frequency          \t$call_freq_omit_count\n");

&print_both("\tTotal number of variants in the filtered file:                \t$filter_count_all\n\n");


#####################################################
# More details of variants written to filtered file #
#####################################################
&print_both("\nMore details of numbers of variants retained for filtered file:\n\n");

&print_both("\tVariants with effect score 5 and max segregation score ($max_seg_score):\t$max_both_scores_count\n\n");

&print_both("\tVariants in which effect score is = 5:                        \t$effect_score_of_5_count\n");
&print_both("\tVariants in which effect score is >= 4:                       \t$effect_score_of_4_count\n");
&print_both("\tVariants in which effect score is >= 3:                       \t$effect_score_of_3_count\n\n");

if ($effect_predictor eq "VEP")
{
&print_both("\tSIFT 'deleterious' positions with effect score < 5:           \t$deleterious_less_than_5_count\n\n");
}

&print_both("\tVariants in which segregation score is = $max_seg_score:                \t$max_seg_score_count\n");
&print_both("\tVariants in which segregation score is >= $second_seg_score:               \t$second_seg_score_count\n");
&print_both("\tVariants in which segregation score is >= $third_seg_score:               \t$third_seg_score_count\n\n");

&print_both("\nFiles used\n\n");

&print_both ("\tInput VCF file:                      \t$vcf_file\n");
&print_both ("\tDisease status file:                 \t$disease_status_file\n");
&print_both ("\tReference genome:                    \t$ref_seq_name\n\n");

&print_both ("\tOutput file (filtered and sorted):   \t$output_file_filtered_and_sorted\n");

if ($keep_unfiltered_output eq "true") {&print_both ("\tOutput file (full):                  \t$output_file\n");}

if ($make_merged_indels_file eq "true") {&print_both ("\n\tMerged Indels file:                  \t$indels_merged_file\n");}

&print_message("The output files can now be viewed in the Excel sheet 'NGS SNP Viewer'","message");

$end_time = time();
$run_time = $end_time - $start_time;
&print_both("\n");

printf "Total run time: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
printf COMMAND_LOG "Total run time: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];

close COMMAND_LOG;
#close CHECK_OUTPUT;
close MISSING_EFFECTS;


######################################################
# Copy command log to command log folder in genetics #
######################################################
$command = "cp $command_log /home/genetics/command_logs/$command_log";
system("$command");


if ($missing_effects_count > 0)
{
	&print_message("Some effects in the VCF file were not recognised by this program. Please check!","warning");
	print "Have a look in the missing effects file: $missing_effects\n\n";
}
exit; # end of program #
################################################################################################################################
################################################################################################################################
################################################################################################################################


############################################################
# Get the name of the effect predictor file: VEP or snpEff #
############################################################
sub get_effect_predictor_filename
{

	########
	# VEP  #
	########
	if ($effect_predictor eq "VEP")
	{
		&print_message("VEP now adds its results to the VCF file so to use the data you should have already processed this VCF file.","message");

		&print_message("Has this VCF file already been processed by VEP?","input");

		print "   <1> Yes\n";
		print "   <2> No\n\n";

		$answer = <STDIN>;
		chomp $answer;

		if (substr($answer,0,1) eq "1")
		{
			print "\n\tOK. This program will use the VEP data\n";
		}
		if (substr($answer,0,1) eq "2")
		{
			&print_message ("OK. Process it with VEP and start again","warning");
			print "(The easiest way is to use the perl script 'run_variant_effect_predictor')\n\n";
			exit;
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
			print "(The easiest way is to use the perl script 'run_snpEff')\n\n";
			exit;
		}
	} # snpEff


	
} # get_effect_predictor_filename



###########################################################
# Ask about how you want to deal with missing genotypes   #
###########################################################
sub ask_about_dealing_with_missing_genotypes
{
	$mark_missing_as_X = "true";

	&print_message("How do you want to mark the allele if a genotype has not been called?","input");

	print "   <1> Mark with X (recommended: you can still choose to score it as REF)       [DEFAULT]\n";
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

		print "   <1> Count X as REF allele                 [DEFAULT]\n";
		print "   <2> Don't count Xs in scoring at all.\n\n";

		$answer = <STDIN>;
		chomp $answer;
		
		if ($answer eq ""){$score_X_as_REF = "true"} # Default
		if (substr($answer,0,1) eq "1" ){$score_X_as_REF = "true"}
		if (substr($answer,0,1) eq "2" ){$score_X_as_REF = "false"}
	}

} # ask_about_dealing_with_missing_genotypes


#############################################
# Ask if you want extra snpEff data columns #
#############################################
sub ask_about_extra_columns
{
	if ($effect_predictor eq "snpEff")
	{
		&print_message("Do you want these extra columns of snpEff data?","input");

		print "\tsnpEff base change data       \te.g. c.2022C>A\n";
		print "\tsnpEff amino acid change data \te.g. p.Asp674Glu\n";
		print "\tsnpEff CDS position data      \te.g. 2022/2922\n";
		print "\tsnpEff protein position data  \te.g. 674/973\n\n";

		print "   <1> Yes. Include extra snpEff data. \n";
		print "   <2> No.  Don't include extra data.             [DEFAULT]\n\n";

	}

	if ($effect_predictor eq "VEP")
	{
		&print_message("Do you want these extra columns of VEP data?","input");

		print "\tVEP base change data        \te.g. cCc/cAc\n";
		print "\tVEP amino acid change data  \te.g. P/H\n";
		print "\tVEP CDS position data       \te.g. 2022\n";
		print "\tVEP protein position data   \te.g. 674\n\n";

		print "   <1> Yes. Include extra VEP data. \n";
		print "   <2> No.  Don't include extra data.             [DEFAULT]\n\n";

	}

	$answer = <STDIN>;
	chomp $answer;
	
	if ($answer eq ""){$answer = "2"} # Default

	if (substr($answer,0,1) eq "1" )
	{
		$include_base_change_info = "true";  # Include  base change position details in the output file
	}
	if (substr($answer,0,1) eq "2" )
	{
		$include_base_change_info = "false";  # Include  base change position details in the output file
	}
} # ask_about_extra_columns


######################################################
# Read in a list of existing SNPs so any SNPs in the #
# output file can be marked as exisiting or new      #
#                                                    #
# This procedure also reads in the list of ensembl   #
# gene names alongside theer gene symbols so the     #
# ensembl name can be converted if required.         #
######################################################
sub load_existing_snp_list
{

	if ($check_snp_list eq "yes")
	{
		&print_message("Which reference sequence do you want to use?","input");

		print "   <1> CanFam3\n";
		print "   <2> EquCab2\n";
		print "   <3> Cat\n";
		print "   <4> Human\n\n";

		$answer=<STDIN>;
		chomp $answer;

		if ($answer eq ""){$answer="1"} # DEFAULT

		if (substr($answer,0,1) eq "1" )
		{
			#$snp_list_file = "/home/genetics/canfam3/canfam3_snps_all.vcf";
			#$ensembl_names_file = "/home/genetics/canfam3/canfam3_ensembl_gene_names.txt"; 
			$snp_list_file = "/geneticsdata/alces-flight-image/canfam3/canfam3_snps_all.vcf";
			$ensembl_names_file = "/geneticsdata/alces-flight-image/canfam3/ensembl_gene_names.txt";
			$ref_seq_name = "canfam3";
			$x_chromosome_number = "39"; 
		}

		if (substr($answer,0,1) eq "2" )
		{
			$snp_list_file = "/home/genetics/equcab2/equcab2_snps_all.vcf";
			$ensembl_names_file = "/home/genetics/equcab2/equcab2_ensembl_gene_names.txt"; 
			$ref_seq_name = "equcab2";
			$x_chromosome_number = "32";
		}

		if (substr($answer,0,1) eq "3" )
		{
			print "\nNo SNP file for cat has been set up in /home/genetics/\n\n";
			$check_snp_list = "no";
			$convert_ensembl_names = "no";
			$ref_seq_name = "felcat5";
			$x_chromosome_number = "19";
		}

		if (substr($answer,0,1) eq "4" )
		{
			print "\nNo SNP file for human has been set up in /home/genetics/\n\n";
			$check_snp_list = "no";
			$convert_ensembl_names = "no";
			$ref_seq_name = "human";
			$x_chromosome_number = "23";
		}


		##################################################
		# Open the input text file with the list of SNPs #
		##################################################
		if ($check_snp_list eq "yes")
		{
			open (SNP_LIST, "$snp_list_file") || die "Cannot open $snp_list_file";

			print "Reading list of known SNPs in file $snp_list_file....\n\n";

			$line_count = 0;

			while ($snp_list_line = <SNP_LIST>) 
			{
					$line_count = $line_count + 1;
					next if $snp_list_line =~ /^#/;

					chomp $snp_list_line;
					@item = split(/\s+/,$snp_list_line);

					if ($line_count > 2)
					{
						$chromosome = $item[0];
						$position   = $item[1];
						$snp_name   = $item[2];

						$chr_pos = $chromosome."_".$position;

						$snp_list_hash{$chr_pos}=$snp_name;

					}
					
					if ($line_count % 500000 == 0)
					{print "  Reading SNP list file.  Line: $line_count\n";}

			}

			close SNP_LIST;

		} # if $check_snp_list eq "yes"


		###################################################
		# Open the file with a list of ensembl gene names #
		###################################################
		if ($convert_ensembl_names eq "yes")
		{
			open (ENSEMBL_NAMES, "$ensembl_names_file") || die "Cannot open $ensembl_names_file";

			print "\nReading ensembl names file $ensembl_names_file....\n\n";

			$line_count = 0;

			print "  Ensembl name\tGene name\n";
			print "  ============\t=========\n";

			while ($ensembl_names_line = <ENSEMBL_NAMES>) 
			{
					$line_count = $line_count + 1;

					chomp $ensembl_names_line;
					@item = split(/\s+/,$ensembl_names_line);

					if ($line_count > 1)
					{
						$ensembl_gene_name = $item[0];
						$real_gene_name   = $item[1];
						$ensembl_names_hash{$ensembl_gene_name}=$real_gene_name;
					}

					if (($line_count < 10) && ($line_count > 1))
					{
						print "  $ensembl_gene_name\t$ensembl_names_hash{$ensembl_gene_name}\n";
					}
					
					if ($line_count % 50000 == 0)
					{print "  Reading ensembl names file.  Line: $line_count\n";}

			}

			print "\n etc etc ...\n\n";
			close ENSEMBL_NAMES;

		} # if $convert_ensembl_names eq "yes"

	} # if $convert_ensembl_names eq "yes"

}


#################################################################
# Check the VCF header lines (before starting on the real data) #
# Check for which variant caller was used.                      #
# Also get information on the order of fields in the INFO field #
# used by VEP or snpEff                                         #
#################################################################
sub check_vcf_header_lines
{
	if ($vcf_file_closed eq "false")
	{
		############################
		# Variant Caller detection #
		############################
		if ($variant_caller eq "")
		{
			if (index(lc($vcf_line),"id=unifiedgenotyper") > -1)
			{
				$variant_caller= "UnifiedGenotyper";
			}
			if (index(lc($vcf_line),"id=haplotypecaller") > -1)
			{
				$variant_caller= "HaplotypeCaller";
			}
			if (index(lc($vcf_line),"samtools") > -1)
			{
				$variant_caller= "mpileup";
			}
			if (index(lc($vcf_line),"freebayes") > -1)
			{
				$variant_caller= "freeBayes";
			}
			if (index(lc($vcf_line),"genotypegvcfs") > -1)
			{
				$variant_caller= "HaplotypeCallergVCF";
			}
		}

		

		################################
		# VEP field order detection    #
		################################
		# ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. 
		# Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|CANONICAL|TSL|CCDS|ENSP|SWISSPROT|TREMBL|UNIPARC|SIFT|DOMAINS|GMAF|AFR_MAF|AMR_MAF|ASN_MAF|EAS_MAF|EUR_MAF|SAS_MAF|AA_MAF|EA_MAF|CLIN_SIG|SOMATIC|PHENO|PUBMED">

		if (index($vcf_line,"#INFO=<ID=CSQ") > -1)
		{
			&print_message("Checking which effect predictor has been used...","message");

			$effect_predictor = "VEP";

			$VEP_fields_string = substr($vcf_line,index($vcf_line,"Format:") + 8);

			###########################################
			# Split VEP fields string at pipes "|"    #
			###########################################
			@VEP_fields_array=split(/\|/, $VEP_fields_string);

			$array_size = scalar @VEP_fields_array;

			for ($array_count = 0; $array_count < $array_size; $array_count++)
			{
				######################################################
				# Remove white space from ends of individual strings #
				######################################################
				$VEP_fields_array[$array_count] =~ s/^\s+|\s+$//g;

				if ($VEP_fields_array[$array_count] eq "IMPACT"){$IMPACT_field_order = $array_count}
				if ($VEP_fields_array[$array_count] eq "SYMBOL"){$SYMBOL_field_order = $array_count}
				if ($VEP_fields_array[$array_count] eq "Gene"){$Gene_field_order = $array_count}
				if ($VEP_fields_array[$array_count] eq "CDS_position"){$CDS_position_field_order = $array_count}
				if ($VEP_fields_array[$array_count] eq "Protein_position"){$Protein_position_field_order = $array_count}
				if ($VEP_fields_array[$array_count] eq "Amino_acids"){$Amino_acids_field_order = $array_count}
				if ($VEP_fields_array[$array_count] eq "Codons"){$Codons_field_order = $array_count}
				if ($VEP_fields_array[$array_count] eq "SIFT"){$SIFT_field_order = $array_count}
			}

			###########################################################
			# Check that all field order positions have been detected #
			###########################################################
			&print_message("VEP annotations detected","message");

			print "   This is the detected order of the data in the VEP 'CSQ' field\n\n";
			print "   IMPACT_field_order:           \t$IMPACT_field_order\n";
			print "   SYMBOL_field_order:           \t$SYMBOL_field_order\n";
			print "   Gene_field_order:             \t$Gene_field_order\n";
			print "   CDS_position_field_order:     \t$CDS_position_field_order\n";
			print "   Protein_position_field_order: \t$Protein_position_field_order\n";
			print "   Amino_acids_field_order:      \t$Amino_acids_field_order\n";
			print "   Codons_field_order:           \t$Codons_field_order\n";
			print "   SIFT_field_order:             \t$SIFT_field_order  \n\n";

			print "   >> NOTE:  If any of these are zero then the data will not be detected and transferred to the output file.\n\n";

			&pause;

		} # if (index($vcf_line,"#INFO=<ID=CSQ") > -1)


		################################
		# snpEff field order detection #
		################################
		## INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 
		# 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">

		if (index($vcf_line,"#INFO=<ID=ANN") > -1)
		{
			&print_message("Checking which effect predictor has been used...","message");

			$effect_predictor = "snpEff";

			$snpEff_fields_string = substr($vcf_line,index($vcf_line,"Format:") + 8);

			###########################################
			# Split snpEff fields string at pipes "|" #
			###########################################
			@snpEff_fields_array=split(/\|/, $snpEff_fields_string);

			$array_size = scalar @snpEff_fields_array;

			for ($array_count = 0; $array_count < $array_size; $array_count++)
			{
				######################################################
				# Remove white space from ends of individual strings #
				######################################################
				$snpEff_fields_array[$array_count] =~ s/^\s+|\s+$//g;

				if ($snpEff_fields_array[$array_count] eq "Annotation_Impact"){$IMPACT_field_order = $array_count}
				if ($snpEff_fields_array[$array_count] eq "Gene_Name"){$SYMBOL_field_order = $array_count}
				if ($snpEff_fields_array[$array_count] eq "Gene_ID"){$Gene_field_order = $array_count}
				if ($snpEff_fields_array[$array_count] eq "cDNA.pos / cDNA.length"){$CDS_position_field_order = $array_count}
				if ($snpEff_fields_array[$array_count] eq "CDS.pos / CDS.length"){$Protein_position_field_order = $array_count}
				if ($snpEff_fields_array[$array_count] eq "HGVS.p"){$Amino_acids_field_order = $array_count}
				if ($snpEff_fields_array[$array_count] eq "HGVS.c"){$Codons_field_order = $array_count}
				if ($snpEff_fields_array[$array_count] eq "SIFT"){$SIFT_field_order = $array_count} # Should be zero as thee is no SIFT field
			}

			###########################################################
			# Check that all field order positions have been detected #
			###########################################################

			&print_message("snpEff annotations detected","message");

			print "   This is the detected order of the data in the snpEff 'ANN' field\n\n";
			print "   Annotation Impact:           \t$IMPACT_field_order\n";
			print "   Gene Symbol:                 \t$SYMBOL_field_order\n";
			print "   Gene Name:                   \t$Gene_field_order\n";
			print "   Codons change:               \t$Codons_field_order\n";
			print "   Amino acids change:          \t$Amino_acids_field_order\n";
			print "   DNA position:                \t$CDS_position_field_order\n";
			print "   Protein position:            \t$Protein_position_field_order\n\n";

			print "   >> NOTE: If any of these are zero then the data will not be detected and transferred to the output file.\n\n";

			&pause;

		} # if (index($vcf_line,"#INFO=<ID=ANN") > -1)

	} # file close eq false

	if (index($vcf_line,"#CHROM") > -1)
	{
		close VCF_TEMPORARY;
		$vcf_file_closed = "true";
	}

} # check_vcf_header_lines



###############################################################
# The original array is passed to the subroutine by reference #
# using &get_data_from_VCF_INFO_field(\@myArray1)             #
# and then dereferenced in the array using @{$_[0]};          #
#                                                             #
# INFO_field_array is the INFO field array which is the 7th element   #
# in vcf_line_array and then is split at semicolons           #
###############################################################
sub get_data_from_VCF_INFO_field
{
	my @vcf_line_array = @{$_[0]};

	#############################################
	# Set the key strings to zero               #
	#############################################
	$snpEff_results_string 					= "";
	$snpEff_gene_string 					= "";
	$snpEff_base_change_string 				= "";
	$snpEff_amino_acid_change_string 		= "";
	$snpEff_CDS_position_string				= "";
	$snpEff_protein_position_string			= "";
	$LOF_string 							= "";
	$LOF_gene_string 						= "";
	$NMD_string 							= "";
	$NMD_gene_string 						= "";
	

	if (scalar @vcf_line_array > 6)
	{ 
		if ($position eq $position_to_debug)
		{
			print "===================================================\n";
			print "Line: $line_count \tPosition: $position\t Whole VCF INFO FIELD\n";
			print "$vcf_line_array[7]\n\n";
			
			$answer=<STDIN>;
		}


		##############################################################
		# Split the INFO fields at semi-colons into INFO_field_array #
		##############################################################		
		@INFO_field_array = split(/;/, $vcf_line_array[7]);
		$INFO_field_array_size = scalar @INFO_field_array;

		for ($array_count = 0; $array_count < $INFO_field_array_size; $array_count++)
		{
			################################################
			# Get which Variant CALLER it is               #
			# Only if put there by pool_vcf_files)         #
			################################################
			if (index ($INFO_field_array[$array_count],"CALLER=") > -1)
			{
				$variant_caller_from_info = $INFO_field_array[$array_count];

				# Remove CALLER= from the start using regex
				$variant_caller_from_info =~ s/CALLER=//g;
			}
			else
			{$variant_caller_from_info = $variant_caller} # if not in INFO get it from header lines

			# get depth of coverage
			if (index ($INFO_field_array[$array_count],"DP=") > -1)
			{
				$depth_coverage_from_info = $INFO_field_array[$array_count];

				# Remove DP= from the start using regex
				$depth_coverage_from_info =~ s/DP=//g;
			}


			###########################################################################
			# Get VEP EFFECT info                                                     #
			# This is marked with the string 'CSQ' this may change update_required    #
			###########################################################################

			if (index ($INFO_field_array[$array_count],"CSQ=") > -1)
			{
				#################################################
				# Record that this VCF fle contains VEP data    #
				#################################################
				$vep_data_found = "true";

				#Set these VEP variables to zero to start with
				$vep_results_string 				= "";
				$vep_gene_string 					= "";
				$vep_base_change_string 			= "";
				$vep_amino_acid_change_string 		= "";
				$vep_CDS_position_string			= "";
				$vep_protein_position_string		= "";
				$vep_sift_prediction_string			= "";
				$max_vep_effect_score 				= 0;
				$vep_effect_score 					= 0;
				$vep_effect 						= "";
				$vep_sift_prediction				= "";
				$sift_deleterious_found				= "false";

				$vep_effect_full_string = $INFO_field_array[$array_count];

				# VEP output
				# snpEff output

				# T|missense_variant|MODERATE|ENPP1|ENSCAFG00000000001|Transcript|ENSCAFT00000000001|protein_coding|18/25||||1919|1802|601|P/H|cCc/cAc|||-1|HGNC|HGNC:3356|tolerated(0.34)|||
				# T|missense_variant|MODERATE|ENPP1|ENSCAFG00000000001|transcript|ENSCAFT00000000001|protein_coding|18/25|c.1802C>A|p.Pro601His|1919/7431|1802/2751|601/916||,

				############################################
				# Remove CSQ= from the start using regex   #
				############################################
				$vep_effect_full_string =~ s/CSQ=//g;


				#####################################################################
				# Split effect string at commas (if there is more than one effect)  #
				# If there is ony one effect then add to first element in the array #
				#####################################################################
				if (index($vep_effect_full_string,",") > -1)
				{
					@vep_effects_array = split (/,/,$vep_effect_full_string);
					$no_of_vep_effects = scalar @vep_effects_array;
				}
				else
				{
					$no_of_vep_effects = 1;
					$vep_effects_array[0] = $vep_effect_full_string;
				}

				## DEBUGGING ##
				if ($position eq $position_to_debug)
				{
					print ">>>>>>>> There are $no_of_vep_effects effects on this line <<<<<<<<\n\n";
				}



				##########################################################
				# Loop through the VEP effects                           #
				# (This string can have several effects on the one line) #
				##########################################################

				for ($effect_count = 1; $effect_count <= $no_of_vep_effects; $effect_count++)
				{		
					$vep_effect_score = 0;
					$vep_effect_string = $vep_effects_array[$effect_count-1];


					###############################################################
					# Split vep_effect_string at pipes into 25 separate fields    #
					# (This includes SIFT data which snpEff doesn't)              #
					###############################################################
	
					# T|missense_variant|MODERATE|ENPP1|ENSCAFG00000000001|Transcript|ENSCAFT00000000001|protein_coding|
					# 18/25||||1919|1802|601|P/H|cCc/cAc|||-1|HGNC|HGNC:3356|tolerated(0.34)|||

					@vep_items_array = split (/\|/,$vep_effect_string);

					$vep_items_array_size = (scalar @vep_items_array) - 1;


					##################################################
					# Assign items in the array to various variables #
					##################################################
					$vep_effect = $vep_items_array[1];

					if ($vep_items_array_size > 1) {$vep_effect_impact = $vep_items_array[$IMPACT_field_order];} else {$vep_effect_impact = ""}

					if ($vep_items_array_size > 2) {$vep_effect_gene_symbol = $vep_items_array[$SYMBOL_field_order];} else {$vep_effect_gene_name = ""} # real name like SRF1
					if ($vep_items_array_size > 3) {$vep_effect_gene_name = $vep_items_array[$Gene_field_order];} else {$vep_effect_gene_name = ""}   # ensembl name like ENSCAFG00000321

					if ($vep_items_array_size > 15) {$vep_effect_base_change = $vep_items_array[$Codons_field_order];} else {$vep_effect_base_change = ""} 
					if ($vep_items_array_size > 14) {$vep_effect_amino_acid_change = $vep_items_array[$Amino_acids_field_order];} else {$vep_effect_amino_acid_change = ""}

					if ($vep_items_array_size > 12) {$vep_CDS_position = $vep_items_array[$CDS_position_field_order];} else {$vep_CDS_position = ""}
					if ($vep_items_array_size > 13) {$vep_protein_position = $vep_items_array[$Protein_position_field_order];} else {$vep_protein_position = ""}

					$vep_sift_prediction = "";
					if ($vep_items_array_size > 29) {$vep_sift_prediction = $vep_items_array[$SIFT_field_order]} else {$vep_sift_prediction = ""}



					###################################################
					# Try to look up real name from ensembl gene name #
					# (Dog only)                                      #
					###################################################
					$vep_effect_real_gene_name = "";
					if ($convert_ensembl_names eq "yes")
					{
						##############################################
						# If gene name is only ensembl name repeated #
						##############################################
						if ($vep_effect_gene_symbol eq $vep_effect_gene_name)
						{
							if (defined $ensembl_names_hash{$vep_effect_gene_name})
							{
								$vep_effect_gene_symbol = $ensembl_names_hash{$vep_effect_gene_name};
								$ensembl_gene_replaced_count = $ensembl_gene_replaced_count + 1;
							}
						}	
					}


					###################################################
					# Choose best gene name                           #
					# If there is a non-ensembl name, use it          #
					###################################################
					if ($vep_effect_gene_symbol ne "") 
					{
						$vep_effect_gene_name = $vep_effect_gene_symbol;

						if ($convert_ensembl_names eq "yes")
						{
							if (($vep_effect_gene_symbol ne $vep_effect_real_gene_name)  && ($vep_effect_real_gene_name ne ""))
							{
								$vep_effect_gene_symbol = $vep_effect_gene_symbol." or ".$vep_effect_real_gene_name;
							}
						}
					}


					## DEBUGGING ##
					if ($position eq $position_to_debug)
					{
						print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nEffect count $effect_count/$no_of_vep_effects\n\n";

						print "$vep_effect_string\n\n";

						print " vep_effect:             \t$vep_effect\n";
						print " vep_effect_impact:      \t$vep_effect_impact\n";
						print " vep_effect_gene_symbol: \t$vep_effect_gene_symbol\n";
						print " vep_effect_gene_name:   \t$vep_effect_gene_name\n";
						print " vep_sift_prediction:    \t$vep_sift_prediction\n";
						$answer=<STDIN>;
					}

					
				
					#####################################################################
					# Set the VEP effect score according to the impact to start with    #
					# (Later this will be looked up properly in a full table)           #
					#####################################################################
					if ($vep_effect_impact eq "MODIFIER") {$vep_effect_score = 1.1}
					elsif ($vep_effect_impact eq "LOW") {$vep_effect_score = 3.1}
					elsif ($vep_effect_impact eq "MODERATE") {$vep_effect_score = 4.1}
					elsif ($vep_effect_impact eq "HIGH") {$vep_effect_score = 5.1}


					##############################################
					# Get VEP effect_score from subroutine       # 
					##############################################

					##############################################
					# If there is NO ampersand in the effect     #
					##############################################
					if (index($vep_effect,"&") == -1)
					{
						$vep_effect_score = &get_effect_score_vep_and_snpEff($vep_effect);

						if ($position eq $position_to_debug)
						{
								print "No ampersand. Effect score from $vep_effect = $vep_effect_score\n";
						}
					}

					##############################################
					# If there IS AN ampersand in the effect     #
					# get the highest effect score in the string #
					# e.g. missense_variant&splice_variant       #
					##############################################
					if (index($vep_effect,"&") > -1)
					{
						# Split list of effects at "&"
						@split_ampersand_array = split ("&",$vep_effect);
						$split_ampersand_array_size = scalar @split_ampersand_array;

						$vep_effect_score_ampersand_max	= 0;

						for ($ampersand_array_count = 0; $ampersand_array_count < $split_ampersand_array_size; $ampersand_array_count++)
						{
							$vep_effect_score = &get_effect_score_vep_and_snpEff($split_ampersand_array[$ampersand_array_count]);

							if ($position eq $position_to_debug)
							{
								print "WITH ampersand. Effect score $ampersand_array_count from $vep_effect = $vep_effect_score\n";
							}

							##########################
							# Keep the highest score #
							##########################
							if ($vep_effect_score > $vep_effect_score_ampersand_max)
							{
								$vep_effect_score_ampersand_max = $vep_effect_score;
							}

							if ($position eq $position_to_debug)
							{
								print "  Max score WITH ampersand. vep_effect_score = $vep_effect_score\tvep_effect_score_ampersand_max = $vep_effect_score_ampersand_max\n";
							}
						}

						$vep_effect_score = $vep_effect_score_ampersand_max;

					} # if there is an ampersand


					if ($position eq $position_to_debug)
					{
						print "  >>> vep_effect_score = $vep_effect_score\n";
					}

					#######################################################
					# Now store the separate data from this effect in a   #
					# set of arrays, e.g. gene name, SIFT predcition etc. #                                   
					#######################################################
					$vep_single_genes_array[$effect_count] = $vep_effect_gene_name;
					$vep_single_effects_array[$effect_count] = $vep_effect;
					$vep_single_effect_scores_array[$effect_count] = $vep_effect_score;
					$vep_single_protein_position_array[$effect_count] = $vep_protein_position;
					$vep_single_CDS_position_array[$effect_count] = $vep_CDS_position;
					$vep_single_base_change_array[$effect_count] = $vep_effect_base_change;
					$vep_single_amino_acid_change_array[$effect_count] = $vep_effect_amino_acid_change;
					$vep_single_sift_prediction_array[$effect_count] = $vep_sift_prediction;


					# Remember maximum effect score for this position
					if ($vep_effect_score > $max_vep_effect_score){$max_vep_effect_score = $vep_effect_score}

					
					#############
					# DEBUGGING #
					#############
					if ($position eq $position_to_debug)
					{
						print "\nLine: $line_count \t CSQ FIELD SECTION\n\n";
						print "$INFO_field_array[$array_count]\n\n";
						print "no_of_vep_effects:       \t$no_of_vep_effects\n";
						print "Effect count:            \t$effect_count\n";
						print "vep_effect:              \t$vep_effect\n";
						print "vep_effect_string:       \t$vep_effect_string\n\n";
						print "vep_effect_gene_name:    \t$vep_effect_gene_name\n";
						print "vep_gene_string:         \t$vep_gene_string\n";
						print "vep_results_string:      \t$vep_results_string\n";

						$answer=<STDIN>;
					}	# debugging

					if ($position eq $position_to_debug){print "\nEnd of effect count $effect_count\n";print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";$answer=<STDIN>;}

				} # VEP effects loop


				if ($position eq $position_to_debug)
				{

					print ">>>>>>>>>>>>   THERE ARE $no_of_vep_effects VEP EFFECTS  <<<<<<<<<<<<<<<<\n\n";

					for ($effect_count = 1; $effect_count <= $no_of_vep_effects; $effect_count++)
					{
						print "\t======\nEffect count: $effect_count\n";
						print "\t\tvep_single_genes_array[$effect_count]:             \t$vep_single_genes_array[$effect_count]\n";
						print "\t\tvep_single_effects_array[$effect_count]:           \t$vep_single_effects_array[$effect_count]\n";
						print "\t\tvep_single_effect_scores_array[$effect_count]:     \t$vep_single_effect_scores_array[$effect_count]\n";
						print "\t\tvep_single_sift_prediction_array[$effect_count]:   \t$vep_single_sift_prediction_array[$effect_count]\n\n\n";
					}

					$answer=<STDIN>;
				}

				############################################################
				# Create results strings based on the list of effects      #
				############################################################
				$max_effect_score = 0;
				for ($effect_count = 1; $effect_count <= $no_of_vep_effects; $effect_count++)
				{
					$effect_score = $vep_single_effect_scores_array[$effect_count];

					if ($effect_score >= $minimum_effect_score)
					{
						
						###################################################################################
						# Now check if this data has been added before by making up a check string of the #
						# four key bits of data, and using it to compare to previous check strings        #
						###################################################################################
						$vep_check_string = $vep_single_effect_scores_array[$effect_count]."_".$vep_single_effects_array[$effect_count]."_".$vep_single_genes_array[$effect_count]."_".$vep_single_sift_prediction_array[$effect_count];

						###################################################
						# Check through all previous effect check_strings #
						###################################################
						$miss_this_one = "false";
						for ($check_count = 1; $check_count < $effect_count; $check_count++)
						{
							$vep_check_string_previous = $vep_single_effect_scores_array[$check_count]."_".$vep_single_effects_array[$check_count]."_".$vep_single_genes_array[$check_count]."_".$vep_single_sift_prediction_array[$check_count];

							if ($vep_check_string_previous eq $vep_check_string)
							{
								$miss_this_one = "true";
								last; # New Oct 2015
							}

							if ($position eq $position_to_debug)
							{
								print "\nCompare with previous strings:\n\n";
								print "This:     $vep_check_string\n";
								print "Previous $check_count: $vep_check_string_previous\n\n";
							}
						}

						########################################################################
						# Only do this bit if the data are not the same as any previous string #
						#######################################################################
						if ($miss_this_one eq "false")
						{
							$vep_results_string = $vep_results_string."+".$vep_single_effects_array[$effect_count];
							$vep_gene_string = $vep_gene_string."+".$vep_single_genes_array[$effect_count];

							########################################
							# Only include these if chosen by user #
							########################################
							if ($include_base_change_info eq "true")
							{
								$vep_base_change_string = $vep_base_change_string." ".$vep_single_base_change_array[$effect_count];
								$vep_amino_acid_change_string = $vep_amino_acid_change_string." ".$vep_single_amino_acid_change_array[$effect_count];
								$vep_CDS_position_string = $vep_CDS_position_string."+".$vep_single_CDS_position_array[$effect_count];
								$vep_protein_position_string = $vep_protein_position_string."+".$vep_single_protein_position_array[$effect_count];
							}
							
							# Add SIFT precition to SIFT prediction string
							if ($vep_single_sift_prediction_array[$effect_count] ne "")
							{
								$vep_sift_prediction_string = $vep_sift_prediction_string."+".$vep_single_sift_prediction_array[$effect_count];
							}
							else
							{
								$vep_sift_prediction_string = $vep_sift_prediction_string."+None";
							}

						} # miss_this_one eq false

					} # if Effect score > minimum threshold effect score


					##########################################################
					# Get the maximum effect score from this loop of effects #
					# (This is written inhe Effect score column)             #
					##########################################################
					if ($effect_score > $max_effect_score)
					{
						$max_effect_score = $effect_score;
					}

					if ($position eq $position_to_debug)
					{
						print "Effect_count: $effect_count/$no_of_vep_effects\n\n";
						print "  vep_results_string:             \t$vep_results_string\n";
						print "  vep_gene_string:                \t$vep_gene_string\n";
						print "  vep_sift_prediction_string:     \t$vep_sift_prediction_string\n";
						print "  max_effect_score:               \t$max_effect_score\n\n";

						$answer=<STDIN>;
					}
					
				} # effect_count loop


				##################################################
				# Remove leading plus sign from combined strings #
				##################################################
				if (index($vep_results_string,"+") == 0)
				{
					$vep_results_string = substr($vep_results_string,1);
				} 
				if (index($vep_gene_string,"+") == 0)
				{
					$vep_gene_string = substr($vep_gene_string,1);
				}
				if (index($vep_base_change_string," ") == 0)
				{
					$vep_base_change_string = substr($vep_base_change_string,1);
				} 
				if (index($vep_amino_acid_change_string," ") == 0)
				{
					$vep_amino_acid_change_string = substr($vep_amino_acid_change_string,1);
				}
				if (index($vep_CDS_position_string,"+") == 0)
				{
					$vep_CDS_position_string = substr($vep_CDS_position_string,1);
				}
				if (index($vep_protein_position_string,"+") == 0)
				{
					$vep_protein_position_string = substr($vep_protein_position_string,1);
				}
				if (index($vep_sift_prediction_string,"+") == 0)
				{
					$vep_sift_prediction_string = substr($vep_sift_prediction_string,1);
				}

				###########################################################
				# Record if SIFT prediction string mentions "deleterious" #
				###########################################################
				if (index($vep_sift_prediction_string,"deleterious") > -1)
				{
					$sift_deleterious_found	= "true";
				}


				####################################################################################
				# Add prefix to position strings so Excel doesn't think they are division formulae #
				####################################################################################
				if ($vep_CDS_position_string ne ""){$vep_CDS_position_string = "CDS: ".$vep_CDS_position_string;}

				if ($vep_protein_position_string ne ""){$vep_protein_position_string = "PROT: ".$vep_protein_position_string;}


				if ($position eq $position_to_debug)
				{
					print "After leading pluses have been removed etc\n";print "no_of_vep_effects: $no_of_vep_effects\n\n";
					print "  vep_results_string:             \t$vep_results_string\n";
					print "  vep_gene_string:                \t$vep_gene_string\n";
					print "  vep_sift_prediction_string:     \t$vep_sift_prediction_string\n";
					print "  max_effect_score:               \t$max_effect_score\n\n";

					$answer=<STDIN>;
				}

				
			} # VEP data found in VCF file


			###########################################################################
			# Get snpEff EFFECT info                                                  #
			# This is marked with the string 'ANN' this may change update_required    #
			###########################################################################

			if (index ($INFO_field_array[$array_count],"ANN=") > -1)
			{
				#################################################
				# Record that this VCF fle contains snpEff data #
				#################################################
				$snpEff_data_found = "true";

				#Set these snpEff variables to zero to start with
				$snpEff_results_string 				= "";
				$snpEff_gene_string 				= "";
				$snpEff_base_change_string 			= "";
				$snpEff_amino_acid_change_string 	= "";
				$snpEff_CDS_position_string			= "";
				$snpEff_protein_position_string		= "";
				$max_snpEff_effect_score 			= 0;
				$snpEff_effect_score 				= 0;
				$snpEff_effect 						= "";

				$snpEff_effect_full_string = $INFO_field_array[$array_count];

				############################################
				# Remove CSQ= from the start using regex   #
				############################################
				$snpEff_effect_full_string =~ s/ANN=//g;


				#####################################################################
				# Split effect string at commas (if there is more than one effect)  #
				# If there is ony one effect then add to first element in the array #
				#####################################################################
				if (index($snpEff_effect_full_string,",") > -1)
				{
					@snpEff_effects_array = split (/,/,$snpEff_effect_full_string);
					$no_of_snpEff_effects = scalar @snpEff_effects_array;
				}
				else
				{
					$no_of_snpEff_effects = 1;
					$snpEff_effects_array[0] = $snpEff_effect_full_string;
				}

				## DEBUGGING ##
				if ($position eq $position_to_debug)
				{
					print ">>>>>>>> There are $no_of_snpEff_effects effects on this line <<<<<<<<\n\n";
				}

				##########################################################
				# Loop through the snpEff effects                        #
				# (This string can have several effects on the one line) #
				##########################################################

				for ($effect_count = 1; $effect_count <= $no_of_snpEff_effects; $effect_count++)
				{		
					$snpEff_effect_score = 0;
					$snpEff_effect_string = $snpEff_effects_array[$effect_count-1];


					###############################################################
					# Split svep_effect_string at pipes into 25 separate fields    #
					# (This includes SIFT data which snpEff doesn't)              #
					###############################################################

					@snpEff_items_array = split (/\|/,$snpEff_effect_string);

					$snpEff_items_array_size = (scalar @snpEff_items_array) - 1;


					##################################################
					# Assign items in the array to various variables #
					##################################################
					$snpEff_effect = $snpEff_items_array[1];

					if ($snpEff_items_array_size > 1) {$snpEff_effect_impact = $snpEff_items_array[$IMPACT_field_order];} else {$snpEff_effect_impact = ""}

					if ($snpEff_items_array_size > 2) {$snpEff_effect_gene_symbol = $snpEff_items_array[$SYMBOL_field_order];} else {$snpEff_effect_gene_name = ""} # real name like SRF1
					if ($snpEff_items_array_size > 3) {$snpEff_effect_gene_name = $snpEff_items_array[$Gene_field_order];} else {$snpEff_effect_gene_name = ""}   # ensembl name like ENSCAFG00000321

					if ($snpEff_items_array_size > 15) {$snpEff_effect_base_change = $snpEff_items_array[$Codons_field_order];} else {$snpEff_effect_base_change = ""} 
					if ($snpEff_items_array_size > 14) {$snpEff_effect_amino_acid_change = $snpEff_items_array[$Amino_acids_field_order];} else {$snpEff_effect_amino_acid_change = ""}

					if ($snpEff_items_array_size > 12) {$snpEff_CDS_position = $snpEff_items_array[$CDS_position_field_order];} else {$snpEff_CDS_position = ""}
					if ($snpEff_items_array_size > 13) {$snpEff_protein_position = $snpEff_items_array[$Protein_position_field_order];} else {$snpEff_protein_position = ""}

					$snpEff_sift_prediction = "";
					if ($snpEff_items_array_size > 29) {$snpEff_sift_prediction = $snpEff_items_array[$SIFT_field_order]} else {$snpEff_sift_prediction = ""}



					###################################################
					# Try to look up real name from ensembl gene name #
					# (Dog only)                                      #
					###################################################
					$snpEff_effect_real_gene_name = "";
					if ($convert_ensembl_names eq "yes")
					{
						##############################################
						# If gene name is only ensembl name repeated #
						##############################################
						if ($snpEff_effect_gene_symbol eq $snpEff_effect_gene_name)
						{
							if (defined $ensembl_names_hash{$snpEff_effect_gene_name})
							{
								$snpEff_effect_gene_symbol = $ensembl_names_hash{$snpEff_effect_gene_name};

								#print "\t>> snpEff lookup. Gene name: $snpEff_effect_gene_name\tGene symbol: $snpEff_effect_gene_symbol\n";

							}
						}	
					}

					###################################################
					# Choose best gene name                           #
					# If there is a non-ensembl name, use it          #
					###################################################
					if ($snpEff_effect_gene_symbol ne "") 
					{
						$snpEff_effect_gene_name = $snpEff_effect_gene_symbol;

						if ($convert_ensembl_names eq "yes")
						{
							if (($snpEff_effect_gene_symbol ne $snpEff_effect_real_gene_name)  && ($snpEff_effect_real_gene_name ne ""))
							{
								$snpEff_effect_gene_symbol = $snpEff_effect_gene_symbol." or ".$snpEff_effect_real_gene_name;
							}
						}
					}


					## DEBUGGING ##
					if ($position eq $position_to_debug)
					{
						print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\nEffect count $effect_count/$no_of_snpEff_effects\n\n";

						print "$snpEff_effect_string\n\n";

						print " snpEff_effect:             \t$snpEff_effect\n";
						print " snpEff_effect_impact:      \t$snpEff_effect_impact\n";
						print " snpEff_effect_gene_symbol: \t$snpEff_effect_gene_symbol\n";
						print " snpEff_effect_gene_name:   \t$snpEff_effect_gene_name\n";
						$answer=<STDIN>;
					}

					
				
					#####################################################################
					# Set the snpEff effect score according to the impact to start with    #
					# (Later this will be looked up properly in a full table)           #
					#####################################################################
					if ($snpEff_effect_impact eq "HIGH") {$snpEff_effect_score = 5.1}
					elsif ($snpEff_effect_impact eq "MODERATE") {$snpEff_effect_score = 4.1}
					elsif ($snpEff_effect_impact eq "LOW") {$snpEff_effect_score = 3.1}
					elsif ($snpEff_effect_impact eq "MODIFIER") {$snpEff_effect_score = 1.1}


					##############################################
					# Get snpEff effect_score from subroutine    # 
					##############################################

					##############################################
					# If there is NO ampersand in the effect     #
					##############################################
					if (index($snpEff_effect,"&") == -1)
					{
						$snpEff_effect_score = &get_effect_score_vep_and_snpEff($snpEff_effect);

						if ($position eq $position_to_debug)
						{
								print "No ampersand. Effect score from $snpEff_effect = $snpEff_effect_score\n";
						}
					}

					##############################################
					# If there IS AN ampersand in the effect     #
					# get the highest effect score in the string #
					# e.g. missense_variant&splice_variant       #
					##############################################
					if (index($snpEff_effect,"&") > -1)
					{
						# Split list of effects at "&"
						@split_ampersand_array = split ("&",$snpEff_effect);
						$split_ampersand_array_size = scalar @split_ampersand_array;

						$snpEff_effect_score_ampersand_max	= 0;

						for ($ampersand_array_count = 0; $ampersand_array_count < $split_ampersand_array_size; $ampersand_array_count++)
						{
							$snpEff_effect_score = &get_effect_score_vep_and_snpEff($split_ampersand_array[$ampersand_array_count]);

							if ($position eq $position_to_debug)
							{
								print "WITH ampersand. Effect score $ampersand_array_count from $snpEff_effect = $snpEff_effect_score\n";
							}

							##########################
							# Keep the highest score #
							##########################
							if ($snpEff_effect_score > $snpEff_effect_score_ampersand_max)
							{
								$snpEff_effect_score_ampersand_max = $snpEff_effect_score;
							}

							if ($position eq $position_to_debug)
							{
								print "  Max score WITH ampersand. snpEff_effect_score = $snpEff_effect_score\tsnpEff_effect_score_ampersand_max = $snpEff_effect_score_ampersand_max\n";
							}
						}

						$snpEff_effect_score = $snpEff_effect_score_ampersand_max;

					} # if there is an ampersand


					if ($position eq $position_to_debug)
					{
						print "  >>> snpEff_effect_score = $snpEff_effect_score\n";
					}

					#######################################################
					# Now store the separate data from this effect in a   #
					# set of arrays, e.g. gene name, SIFT predcition etc. #                                   
					#######################################################
					$snpEff_single_genes_array[$effect_count] = $snpEff_effect_gene_name;
					$snpEff_single_effects_array[$effect_count] = $snpEff_effect;
					$snpEff_single_effect_scores_array[$effect_count] = $snpEff_effect_score;
					$snpEff_single_protein_position_array[$effect_count] = $snpEff_protein_position;
					$snpEff_single_CDS_position_array[$effect_count] = $snpEff_CDS_position;
					$snpEff_single_base_change_array[$effect_count] = $snpEff_effect_base_change;
					$snpEff_single_amino_acid_change_array[$effect_count] = $snpEff_effect_amino_acid_change;


					# Remember maximum effect score for this position
					if ($snpEff_effect_score > $max_snpEff_effect_score){$max_snpEff_effect_score = $snpEff_effect_score}

					
					#############
					# DEBUGGING #
					#############
					if ($position eq $position_to_debug)
					{
						print "\nLine: $line_count \t ANN FIELD SECTION\n\n";
						print "$INFO_field_array[$array_count]\n\n";
						print "no_of_snpEff_effects:       \t$no_of_snpEff_effects\n";
						print "Effect count:            \t$effect_count\n";
						print "snpEff_effect:              \t$snpEff_effect\n";
						print "snpEff_effect_string:       \t$snpEff_effect_string\n\n";
						print "snpEff_effect_gene_name:    \t$snpEff_effect_gene_name\n";
						print "snpEff_gene_string:         \t$snpEff_gene_string\n";
						print "snpEff_results_string:      \t$snpEff_results_string\n";

						$answer=<STDIN>;
					}	# debugging

					if ($position eq $position_to_debug){print "\nEnd of effect count $effect_count\n";print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";$answer=<STDIN>;}

				} # snpEff effects loop


				if ($position eq $position_to_debug)
				{

					print ">>>>>>>>>>>>   THERE ARE $no_of_snpEff_effects snpEff EFFECTS  <<<<<<<<<<<<<<<<\n\n";

					for ($effect_count = 1; $effect_count <= $no_of_snpEff_effects; $effect_count++)
					{
						print "\t======\nEffect count: $effect_count\n";
						print "\t\tsnpEff_single_genes_array[$effect_count]:             \t$snpEff_single_genes_array[$effect_count]\n";
						print "\t\tsnpEff_single_effects_array[$effect_count]:           \t$snpEff_single_effects_array[$effect_count]\n";
						print "\t\tsnpEff_single_effect_scores_array[$effect_count]:     \t$snpEff_single_effect_scores_array[$effect_count]\n";
					}

					$answer=<STDIN>;
				}

				############################################################
				# Create results strings based on the list of effects      #
				############################################################
				$max_effect_score = 0;
				for ($effect_count = 1; $effect_count <= $no_of_snpEff_effects; $effect_count++)
				{
					$effect_score = $snpEff_single_effect_scores_array[$effect_count];

					if ($effect_score >= $minimum_effect_score)
					{
						
						###################################################################################
						# Now check if this data has been added before by making up a check string of the #
						# four key bits of data, and using it to compare to previous check strings        #
						###################################################################################
						$snpEff_check_string = $snpEff_single_effect_scores_array[$effect_count]."_".$snpEff_single_effects_array[$effect_count]."_".$snpEff_single_genes_array[$effect_count];

						###################################################
						# Check through all previous effect check_strings #
						###################################################
						$miss_this_one = "false";
						for ($check_count = 1; $check_count < $effect_count; $check_count++)
						{
							$snpEff_check_string_previous = $snpEff_single_effect_scores_array[$check_count]."_".$snpEff_single_effects_array[$check_count]."_".$snpEff_single_genes_array[$check_count];

							if ($snpEff_check_string_previous eq $snpEff_check_string)
							{
								$miss_this_one = "true";
								last; # new Oct 2015
							}

							if ($position eq $position_to_debug)
							{
								print "\nCompare with previous strings:\n\n";
								print "This:     $snpEff_check_string\n";
								print "Previous $check_count: $snpEff_check_string_previous\n\n";
							}
						}

						########################################################################
						# Only do this bit if the data are not the same as any previous string #
						#######################################################################
						if ($miss_this_one eq "false")
						{
							$snpEff_results_string = $snpEff_results_string."+".$snpEff_single_effects_array[$effect_count];
							$snpEff_gene_string = $snpEff_gene_string."+".$snpEff_single_genes_array[$effect_count];

							########################################
							# Only include these if chosen by user #
							########################################
							if ($include_base_change_info eq "true")
							{
								$snpEff_base_change_string = $snpEff_base_change_string." ".$snpEff_single_base_change_array[$effect_count];
								$snpEff_amino_acid_change_string = $snpEff_amino_acid_change_string." ".$snpEff_single_amino_acid_change_array[$effect_count];
								$snpEff_CDS_position_string = $snpEff_CDS_position_string."+".$snpEff_single_CDS_position_array[$effect_count];
								$snpEff_protein_position_string = $snpEff_protein_position_string."+".$snpEff_single_protein_position_array[$effect_count];
							}

						} # miss_this_one eq false

					} # if Effect score > minimum threshold effect score


					##########################################################
					# Get the maximum effect score from this loop of effects #
					# (This is written inhe Effect score column)             #
					##########################################################
					if ($effect_score > $max_effect_score)
					{
						$max_effect_score = $effect_score;
					}

					if ($position eq $position_to_debug)
					{
						print "Effect_count: $effect_count/$no_of_snpEff_effects\n\n";
						print "  snpEff_results_string:             \t$snpEff_results_string\n";
						print "  snpEff_gene_string:                \t$snpEff_gene_string\n";
						print "  max_effect_score:               \t$max_effect_score\n\n";

						$answer=<STDIN>;
					}
					
				} # effect_count loop


				##################################################
				# Remove leading plus sign from combined strings #
				##################################################
				if (index($snpEff_results_string,"+") == 0)
				{
					$snpEff_results_string = substr($snpEff_results_string,1);
				} 
				if (index($snpEff_gene_string,"+") == 0)
				{
					$snpEff_gene_string = substr($snpEff_gene_string,1);
				}
				if (index($snpEff_base_change_string," ") == 0)
				{
					$snpEff_base_change_string = substr($snpEff_base_change_string,1);
				} 
				if (index($snpEff_amino_acid_change_string," ") == 0)
				{
					$snpEff_amino_acid_change_string = substr($snpEff_amino_acid_change_string,1);
				}
				if (index($snpEff_CDS_position_string,"+") == 0)
				{
					$snpEff_CDS_position_string = substr($snpEff_CDS_position_string,1);
				}
				if (index($snpEff_protein_position_string,"+") == 0)
				{
					$snpEff_protein_position_string = substr($snpEff_protein_position_string,1);
				}


				####################################################################################
				# Add prefix to position strings so Excel doesn't think they are division formulae #
				####################################################################################
				if ($snpEff_CDS_position_string ne ""){$snpEff_CDS_position_string = "CDS: ".$snpEff_CDS_position_string;}

				if ($snpEff_protein_position_string ne ""){$snpEff_protein_position_string = "PROT: ".$snpEff_protein_position_string;}


				if ($position eq $position_to_debug)
				{
					print "After leading pluses have been removed etc\n\n";
					print "  snpEff_results_string:             \t$snpEff_results_string\n";
					print "  snpEff_gene_string:                \t$snpEff_gene_string\n";
					print "  max_effect_score:                  \t$max_effect_score\n\n";

					$answer=<STDIN>;
				}

				
			} # snpEff data found in VCF file 2



			#########################################################################
			# Look for the LOF tag to see if snpEff has added Loss Of Function data #
			#########################################################################
			if (index ($INFO_field_array[$array_count],"LOF=") > -1)
			{
				#################################################
				# Record that this VCF file contains LOF data   #
				#################################################
				$LOF_data_found = "true";

				$LOF_full_string = $INFO_field_array[$array_count];

				$LOF_string 		= "";
				$LOF_gene_string 	= "";

				# Remove LOF= from the start using regex
				$LOF_full_string =~ s/LOF=//g;


				# Split LOF string at commas (if there is more than one effect) #
				# If there is ony one effect then add to first element in the array #`
				if (index($LOF_full_string,",") > -1)
				{
					@LOF_array = split (/,/,$LOF_full_string);
					$no_of_LOFs = scalar @LOF_array;
				}
				else
				{
					$no_of_LOFs = 1;
					$LOF_array[0] = $LOF_full_string;
				}


				##########################################################
				# Loop through the LOF   (Loss Of Function)              #
				# (This string can have several LOFs on the one line)    #
				########################################################## 111
				for ($LOF_count = 1; $LOF_count <= $no_of_LOFs; $LOF_count++)
				{

					$LOF_string = $LOF_array[$LOF_count-1];

					# Remove brackets
					$LOF_string = substr($LOF_string,1,length($LOF_string)-2);
					#$LOF_full_string = substr($LOF_full_string,1,99);

					###############################################################
					# Split LOF_string at pipes into separate fields              #
					###############################################################
	
					# ANN=T|upstream_gene_variant|MODIFIER|HSBP1L1|ENSCAFG00000031133|transcript|ENSCAFT00000045122|protein_coding||c.-1_-1insA|||||103|WARNING_TRANSCRIPT_NO_START_CODON,
				
					@LOF_items_array = split (/\|/,$LOF_string);

					$LOF_items_array_size = (scalar @LOF_items_array) - 1;

					###############################################################
					# Assign items in the array to various variables              #
					###############################################################
					$LOF_gene = $LOF_items_array[0];


					###############################################################
					# Add snpEff result to string if not already in the string    #
					###############################################################
					if (index($LOF_gene_string,$LOF_gene)  == -1)
					{
						$LOF_gene_string = $LOF_gene_string."+".$LOF_gene;
					}
					
				} # LOF loop

				# Remove leading plus sign
				if (index($LOF_gene_string,"+") == 0)
				{
					$LOF_gene_string = substr($LOF_gene_string,1);
				} 


			}# If LOF data found in VCF file (loss of function)


			######################################################################
			# If effect predictor is VEP then put SIFT predictions in LOF column #
			######################################################################


			################################################################################
			# Look for the NMD tag to see if snpEff has added Nonsense Mediated Decay data #
			################################################################################
			if (index ($INFO_field_array[$array_count],"NMD=") > -1)
			{
				#################################################
				# Record that this VCF fle contains NMD data    #
				#################################################
				$NMD_data_found = "true";

				$NMD_full_string = $INFO_field_array[$array_count];

				$NMD_string 		= "";
				$NMD_gene_string 	= "";

				# Remove NMD= from the start using regex
				$NMD_full_string =~ s/NMD=//g;


				# Split NMD string at commas (if there is more than one effect) #
				# If there is ony one effect then add to first element in the array #`
				if (index($NMD_full_string,",") > -1)
				{
					@NMD_array = split (/,/,$NMD_full_string);
					$no_of_NMDs = scalar @NMD_array;
				}
				else
				{
					$no_of_NMDs = 1;
					$NMD_array[0] = $NMD_full_string;
				}


				##########################################################
				# Loop through the NMD   (Loss Of Function)              #
				# (This string can have several NMDs on the one line)    #
				########################################################## 111
				for ($NMD_count = 1; $NMD_count <= $no_of_NMDs; $NMD_count++)
				{

					$NMD_string = $NMD_array[$NMD_count-1];

					# Remove brackets
					$NMD_string = substr($NMD_string,1,length($NMD_string)-2);
					#$NMD_full_string = substr($NMD_full_string,1,99);

					###############################################################
					# Split NMD_string at pipes into separate fields              #
					###############################################################
	
					# ANN=T|upstream_gene_variant|MODIFIER|HSBP1L1|ENSCAFG00000031133|transcript|ENSCAFT00000045122|protein_coding||c.-1_-1insA|||||103|WARNING_TRANSCRIPT_NO_START_CODON,
				
					@NMD_items_array = split (/\|/,$NMD_string);

					$NMD_items_array_size = (scalar @NMD_items_array) - 1;

					##################################################
					# Assign items in the array to various variables #
					##################################################
					$NMD_gene = $NMD_items_array[0];


					###########################################
					# Get snpEff effect_score from subroutine #
					###########################################
					#if (index($NMD_string,"&") > -1)
					#{
					#	print "NMD_string: $NMD_string has ampersand\n";
					#	$answer=<STDIN>;
					#}


					###############################################################
					# Add snpEff result to string if not already in the string    #
					###############################################################
					if (index($NMD_gene_string,$NMD_gene)  == -1)
					{
						$NMD_gene_string = $NMD_gene_string."+".$NMD_gene;
					}

					
				} # NMD loop

				# Remove leading plus sign
				if (index($NMD_gene_string,"+") == 0)
				{
					$NMD_gene_string = substr($NMD_gene_string,1);
				} 


			}# If NMD data found in VCF file (loss of function)

			#################################################################################
			# If INFO_field_array contains TYPE=snp, TYPE=del etc then it is from freeBayes #
			# If there are Indels where the REF_base is the same as the ALT then we         #
			# have to set the variant_type as "indel" (or maybe indel_FB" ??)               #
			#################################################################################

			# e.g. at position 50693881

			#########################################################################
			# If there is a TYPE= in the INFO line then it is called with freeBayes #
			#########################################################################
			if (index ($vcf_line_array[7],"TYPE=") > -1)
			{
				$variant_caller="freeBayes";
			}
		}
	} # if ($array_size_1 > 6)


} # sub get_data_from_VCF_INFO_field



sub get_data_from_VCF_GENOTYPES_field
{
	my @vcf_line_array = @{$_[0]};

	################################################
    # Read the FORMAT in myArray1(8)               #
    #                                              #
    # This looks like this:  GT:AD:DP:GQ:PL        #
    # The genotype data is then read in this order #
    ################################################
    @myArray8 = split(":", $vcf_line_array[8]);
    $array_size_8 = scalar @myArray8;
    

    ##########################################################
    # Get the genotypes (i.e. all fields in myArray1(9)      #
    ##########################################################
    if ($array_size_8 > 0)
    {
        for($sample_count = 1; $sample_count <= $no_of_samples; $sample_count++)
		{
            $genotype_array_input_order[$sample_count] = $myArray1[$sample_count + 8];
        }
    }       


	############################################################################
	# RE_MAP the genotype_array data as this is basically the only data we use #
	# (the sample name array and the disease status array were re-mapped       #
	# earlier when processed the #CHROM line)                                  #                           
	############################################################################
	for($sample_count = 1; $sample_count <= $no_of_samples; $sample_count++)
	{
        # Get original position (source column)
        $source_col = $source_column_array[$sample_count];
        $genotype_array[$sample_count] = $genotype_array_input_order[$source_col];
    }

} # get_data_from_VCF_GENOTYPES_field


############################################################################################
# This procedure looks at the FORAMT field to determine where the genotype data is exactly #
# Then it gets the genotype in the form '0/1'  '1/1'  1/2' etc                             #
# It also converts the '0/1' format into actual bases (e.g. 'C/G'  or 'A/C' etc)           #
# The original 'bases' go into the array 'base_array_original' (includes Xs)               #            
# The actual bases go into the array 'base_array'                                          #
############################################################################################
sub get_genotypes_from_all_samples
{
	my @genotype_array_local = @{$_[0]};

	#################################################################################################
    # Now work through all the re-mapped samples, splitting each separate genotype_array at colons  #
    #################################################################################################
    
    $GT_string = "";
	$AD_string = "";
	$DP_genotype = "";
	$GQ_string = "";
	$PL_string = "";
	$missing_genotypes_at_this_position = 0; # This is used to determine the call frequency at this position
	$call_frequency_at_this_position = 0;    # This is the call frequency at this position (percentage of good calls)

    for($sample_count = 1;$sample_count <=$no_of_samples; $sample_count++)
	{
		############################################################
        # Split myArray9 (genotype data) using colons              #
		# Remember if genotype is ./. then there may be no colons  #
		# but splitting should still work if GT is there           #
		#                                                          #
		# It is possible to have a single dot (as well as ./.) for #
		# the whole field so we will have to deal with that too.   #
		############################################################
       				
        @myArray9 = split(":",$genotype_array_local[$sample_count]);
        $array_size_9 = scalar @myArray9;
        
	
        ################################################################
        # Now use order of fields in 'FORMAT' (myArray8)               #
        # to parse the genotype data fields of                         # 
        # each sample (order may not be the same for every VCF file)   #
        #                                                              #
        # The FORMAT field and the equivalent data field for each      #
        # sample might look like this:                                 #
        #                                                              #
        # FORMAT: GT:AD:DP:GQ:PL         <- this could be in any order #
        # SAMPLE: 0/1:0,251:251:99:10318,0   <- genotype is '0/1'      #
        #                                                              #
        # So in this case the GT part is 0/1, and the DP part is 251   #
        ################################################################
        for($array_count = 0; $array_count < $array_size_8; $array_count++)
		{
			$myString = $myArray8[$array_count];
            
			if ($myString eq "GT"){$GT_string = $myArray9[$array_count]}  # This is the field we are interested in - Genotype
			elsif ($myString eq "AD"){$AD_string = $myArray9[$array_count]}
			elsif ($myString eq "DP"){$DP_genotype = $myArray9[$array_count]}
			elsif ($myString eq "GQ"){$GQ_string = $myArray9[$array_count]}
			elsif ($myString eq "PL"){$PL_string = $myArray9[$array_count]}

        } # next array_count
        

		##############################################################################
		# Decide whether the variant type is a SNP or an Indel                       # 
		# If any variant is a different length than the REF base then it is an indel #
		##############################################################################
		$variant_type = "snp";
		for ($allele_count = 0; $allele_count <$no_of_alt_alleles; $allele_count++)
		{
			#freeBayes only
			if ($variant_caller eq "freeBayes")
			{
				if ((length($REF_base) > 1) || (length($ALT_allele_array[$allele_count]) > 1))
				{
					$variant_type = "indel";
					last;
				}
			} # caller is freeBayes

			#All
			if (length($REF_base) != length($ALT_allele_array[$allele_count]))
			{
				$variant_type = "indel";
				last;
			}
		}


		#####################################################
		# FreeBayes has a single dot . for 'no genotype'    #
		# so there is no slash.  Convert "." to "./." first #
		# so it behaves as GATK VCF files                   #
		#####################################################	
		#freeBayes only
		if ($variant_caller eq "freeBayes")
		{
			if ($GT_string eq "."){$GT_string = "./."}
		}
		else
		{
			if ($genotype_array_local[$sample_count] eq ".")
			{
				#################################################################
				# This dot is occurring when two VCF files are merged, and one  #
				# of the files has no variants at a particular position, so it  #
				# puts simply a dot '.' for the data field for that sample      #
				# This should be the REF_base, so we set GT_string as "0"       #
				#################################################################
				$GT_string = "0/0"; # Set GT_string as REF_base
			}
		} # not Freebayes


		##################################################################
		# Now convert the GT string from 0/1 format to ACGT format       #
		#                                                                #
		#   - allele_number: This is the number from the string '0/1'    #
		#   - allele_letter: This is the assigned alleles A,B,C,D etc    #
		#   - allele_base:   This is the actual base (for SNPs)          #
		##################################################################
        $slash_pos = index($GT_string, "/");
        
		$allele_number_1 = "";
		$allele_number_2 = "";

		$allele_base_1 = "";
		$allele_base_2 = "";
		

		####################################################################
		# If there is a slash (i.e. if genotype data is in the form '0/1') #
		####################################################################	
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
			
			###################################
			# If genotype is OK, i.e. NOT ./. #
			###################################
			if ($allele_number_1 ne ".")
            {
            	if ($allele_number_1 == 0){$allele_base_1 = $REF_base}
				elsif ($allele_number_1 > 0){$allele_base_1 = $ALT_allele_array[$allele_number_1 - 1]}
				
				if ($allele_number_2 == 0){$allele_base_2 = $REF_base}
				elsif ($allele_number_2 > 0){$allele_base_2 = $ALT_allele_array[$allele_number_2 - 1]}

			} # if not ./.


			###################################
			# If genotype IS ./.              #
			###################################
			if ($allele_number_1 eq ".")
            {

				##############################################
				# $mark_missing_as_X is set at true nowadays #
				# Anything else is not sensible              # 
				##############################################
				$allele_base_1 = "X"; 
				$allele_base_2 = "X";


				#####################################
				# Count how many of these there are #
				# in total and for each sample 	    #
				#####################################
				$missing_genotypes_total = $missing_genotypes_total + 1;
				$missing_genotypes_array[$sample_count] = $missing_genotypes_array[$sample_count] + 1;

				$missing_genotypes_at_this_position = $missing_genotypes_at_this_position + 1;

        	} # ./.

        } # if there is a slash

			
		###############################################
		# Store original base in case it is X         #
		###############################################
		$col_1 = ($sample_count * 2) - 1;
		$col_2 = ($sample_count * 2);
		$base_orig_array[$col_1] = $allele_base_1;
		$base_orig_array[$col_2] = $allele_base_2;


		###############################################
		# Store the two bases in an array of bases    #
		# for later use in statistics                 #
		# (Note these are the actual bases, ACGT etc) #
		###############################################
		$base_array[$col_1] = $allele_base_1;
		$base_array[$col_2] = $allele_base_2;


		###############################################
		# Update some counters of the genotype        #
		###############################################
		$genotypes_counted_total = $genotypes_counted_total + 1;
		$genotypes_counted_array[$sample_count] = $genotypes_counted_array[$sample_count] + 1;
		
		$GT_string_array[$sample_count] = $GT_string; # store for later output checking

	} # sample_count loop  to process genotypes


	###############################################################################
	# The call frequency is the percentage of samples having a good genotype call #
	###############################################################################
	$call_frequency_at_this_position = ($no_of_samples_without_omits - $missing_genotypes_at_this_position)/$no_of_samples_without_omits * 100;
	$call_frequency_at_this_position = sprintf("%.1f", $call_frequency_at_this_position);

} # get_genotypes_from_all_samples



#########################################################################
# New way to get main affected allele (works for any number of alleles) #
# This doesn't need allele_A and allele_B  .                            #
#                                                                       #
# The variable $main_affected_allele holds the actual base (A,C,G,T)    #
#                                                                       #
# When analysing for segregation score the main affected allele is 
# used to specify what the main affected genotype is, i.e. 'AA'
# (CAn it also be 'AB'? ) or 'AB'
#########################################################################
sub get_main_affected_allele
{
	$main_affected_allele   = "";
	$second_affected_allele = "";

	$max_no_of_alleles = 0;       # Maximum no of alleles for any of all the ALT alleles
	$second_no_of_alleles = 0;    # Second max no of alleles
	$max_allele_count = 0;        # allele_count at which maximum no of alleles occurs (i.e. array position)
	$second_max_allele_count = 0; # second allele_count (ditto)

	########################################################################
	# Loop through all the alleles to find the most common affected allele #
	########################################################################
	for ($allele_count = 0; $allele_count <= $no_of_alt_alleles; $allele_count++)
	{
		$no_of_alleles = 0;
		#########################################
		# Loop through all the Affected samples #
		#########################################
		for ($sample_count = 1; $sample_count <= $no_affected_samples; $sample_count++)
		{
			$col_1 = ($sample_count * 2) - 1; $col_2 = $sample_count * 2;

			if ($base_array[$col_1] eq $ALT_allele_array[$allele_count]){$no_of_alleles = $no_of_alleles + 1}
			if ($base_array[$col_2] eq $ALT_allele_array[$allele_count]){$no_of_alleles = $no_of_alleles + 1}
		} 

		######################################################################
		# Record maximum no of alleles and at which allele_count it occurred #
		######################################################################
		if ($no_of_alleles >= $max_no_of_alleles)
		{
			#Store second best first (before it changes)
			$second_max_allele_count = $max_allele_count;
			$second_no_of_alleles = $max_no_of_alleles;

			#Store best
			$max_no_of_alleles = $no_of_alleles;
			$max_allele_count = $allele_count;
		}
	}

	$second_affected_allele = $ALT_allele_array[$second_max_allele_count];
	$main_affected_allele = $ALT_allele_array[$max_allele_count];


	##########################################################################################
	# If top two alleles have the same numbers, use normals (and carriers, if any) to decide #
	##########################################################################################
	if ($max_no_of_alleles == $second_no_of_alleles)
	{
		$no_of_alleles = 0;
		$second_no_of_alleles = 0;

		########################################################
		# Look through normals and count top and second allele #
		# (This also counts through carriers if any)           #
		########################################################
		for ($sample_count = $no_affected_samples + 1; $sample_count <= $no_of_samples_without_omits; $sample_count++)
		{
			$col_1 = ($sample_count * 2) - 1; $col_2 = $sample_count * 2;

			if ($base_array[$col_1] eq $main_affected_allele){$no_of_alleles = $no_of_alleles + 1}
			if ($base_array[$col_1] eq $second_affected_allele){$second_no_of_alleles = $second_no_of_alleles + 1}
			if ($base_array[$col_2] eq $main_affected_allele){$no_of_alleles = $no_of_alleles + 1}
			if ($base_array[$col_2] eq $second_affected_allele){$second_no_of_alleles = $second_no_of_alleles + 1}
		} 

		##################################################################################
		# If there are more main than second in the NORMALS, then replace main by second #
		##################################################################################
		if ($no_of_alleles > $second_no_of_alleles)
		{
			$main_affected_allele = $second_affected_allele;
		}
	} # if ($max_no_of_alleles == $second_no_of_alleles)

	$allele_A = $main_affected_allele;

	if ($no_of_alleles == 2){$allele_B = $second_affected_allele;}
	elsif ($no_of_alleles > 2){$allele_B = "M";}

} # end get_main_affected_allele



##############################################################################
# Now count Affected GENOTYPES.                                              #
# This starts scoring segregation scores for each variant                    #
# It puts the genotype to be used for scoring into the variable $genotype_AB #
# Then it scores for segregation comparing 'AA', 'AB' and 'BB' etc           #
# I think this works for Indels as well                                    
##############################################################################
sub score_affected_genotypes
{	 
	 #SCORE AFFECTED GENOTYPES
	for ($sample_count = 1; $sample_count <= $no_affected_samples; $sample_count++)
	{		 
		$col_1 = ($sample_count * 2 ) - 1;
		$col_2 = ($sample_count * 2 );
		 
		$allele_base_1 = $base_array[$col_1]; # These are the actual bases A,C,G or T
		$allele_base_2 = $base_array[$col_2];

		$original_genotype = $allele_base_1.$allele_base_2;
		$genotype_AB = "";


		 ###############################################
		 # Affecteds:                                  #
		 # $score_X_method can be three things         #
		 # 1. "reference"  (old method)                #
		 # 2. "best"  (what gives highest seg score)   #
		 # 3. "missing" (ignores missing data)         #
		 #     (Then scaled up seg score in proportion #
		 #      to the number of samples used)         #
		 ###############################################
		 if ($original_genotype ne "XX")
		 {
			 #######################################################################################
			 # If allele is main affected allele call it 'A', otherwise call it 'B'                #
			 # This is only used for the segregation scores so only AA, AB and BB make sense       #
			 #######################################################################################
			
			if    (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AA"}
		 	elsif (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "AB"}
		 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AB"}
		 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "BB"}
		 }
		 else
		 {
		 	if ($score_X_method eq "reference")
			 {
				 if ($REF_base eq $main_affected_allele){$genotype_AB = "AA"}else{$genotype_AB = "BB"}
			 }
			 elsif ($score_X_method eq "best")
			 {
				 $genotype_AB = "AA"; # Because these are Affecteds the 'best' option would be AA
			 }
			 elsif ($score_X_method eq "missing")
			 {
				 $genotype_AB = "MM";
			 }
		 }


		 ################################################################################################################
		 # Counter increase for genotypes containing one OR two copies of the main affected allele (for Dominant score) #
		 ################################################################################################################
		 if (($allele_base_1 eq $main_affected_allele) || ($allele_base_2 eq $main_affected_allele))
		 {
			$affected_dom_count = $affected_dom_count + 1; # NOT NEEDED?
		 }
		
		 ############################################################################
		 # If you are comparing a single case against a whole loads of normals then #
		 # you need to know what the genotype of the single case is                 #
		 ############################################################################
		if ($sample_count == 1 ){$single_case_genotype = $genotype_AB} else {$single_case_genotype = "WW"};


		 #########################################
		 # Affecteds:                            #
		 # Count AAs, ABs and BBs                #
		 # These are used later in ????????????  #
		 #########################################
		 if    ($genotype_AB eq "AA"){$AA_count_affected = $AA_count_affected + 1}
		 elsif ($genotype_AB eq "AB"){$AB_count_affected = $AB_count_affected + 1}
		 elsif ($genotype_AB eq "BB"){$BB_count_affected = $BB_count_affected + 1}
		 elsif ($genotype_AB eq "MM"){$MM_count = $MM_count + 1}
		 

		 #############
		 # DEBUGGING #
		 #############
		 if ($position eq $position_to_debug)
		 {
		 	print "AFFECTED\n";
		 	print "Position      \t$position\tSample: $sample_count\tVariant type: \t$variant_type\n";
		 	print "Allele_base_1:\t$allele_base_1    \tAllele_base_2: $allele_base_2\n\n";

		 	print " allele_1_for_writing_array[$sample_count]: $allele_1_for_writing_array[$sample_count]\n";
		 	print " allele_1_for_writing_array[$sample_count]: $allele_1_for_writing_array[$sample_count]\n\n";

		 	print "score_X_method: $score_X_method\n\n";
		 	print "Original genotype: $original_genotype\tGenotype for scoring: $genotype_AB\n\n";

		 	if ($original_genotype eq "XX")
		 	{
		 		print "NOTE: the setting of 'score_X_method' decides what an 'XX' genotype will be scored as.\n\n";
		 		print "      'reference':  the genotype is set as the reference base ($REF_base) [could be AA or BB depending on the main main_affected_genotype which in this case is $main_affected_genotype]\n";
		 		print "      'best':       the genotype for affecteds is set as AA (which will give the highest segregation score))\n";
		 		print "      'missing':    the genotype is set as 'MM' and is not used in the scoring\n";
		 	}
		 	
		 	$answer=<STDIN>;
		 }


		#####################################
		# Segregation scoring for Affecteds #
		#####################################
		# 1 AA_AB_BB         	Recessive:  								Affecteds AA, carriers AB, normals BB
		# 2 AA_AB_ABorBB     	Recessive, some normals may be carriers:  	Affecteds AA, carriers AB, normals AB or BB
		# 3 AA_ABorBB_BB        Recessive, some carriers may be normals:  	Affecteds AA, carriers AB or BB, normals BB
		# 4 AA_ABorBB_ABorBB 	Recessive:  								Affecteds AA, carriers AB or BB, normals AB or BB
		# 5 AAorAB_BB        	Dominant:   								Affecteds AA or AB, normals BB
		# 6 A_B              	Additive:   								Affecteds with As, Normals with Bs
		# 7 AA_BB               Difference:                                 Normals different to main affected homozygous genotype
		
		######################################################
		# Affected: New meaning for seg score 7 - Difference #
		######################################################
		if ($genotype_AB ne $main_affected_homozygous_genotype)
		{
			$segregation_score_7 = $segregation_score_7 + 1;
			$segregation_score_7B = $segregation_score_7B + 1;
		}

		if ($genotype_AB eq "AA")
		{
			$segregation_score_1 = $segregation_score_1 + 1;
			$segregation_score_2 = $segregation_score_2 + 1;
			$segregation_score_3 = $segregation_score_3 + 1;
			$segregation_score_4 = $segregation_score_4 + 1;
			$segregation_score_5 = $segregation_score_5 + 1;
			$segregation_score_6 = $segregation_score_6 + 2;
		}
		elsif ($genotype_AB eq "AB")
		{
			$segregation_score_5 = $segregation_score_5 + 1;
			$segregation_score_6 = $segregation_score_6 + 1;

			$segregation_score_5B = $segregation_score_5B + 1;
			$segregation_score_6B = $segregation_score_6B + 1;
		}
		elsif ($genotype_AB eq "BB")
		{
			$segregation_score_1B = $segregation_score_1B + 1;
			$segregation_score_2B = $segregation_score_2B + 1;
			$segregation_score_3B = $segregation_score_3B + 1;
			$segregation_score_4B = $segregation_score_4B + 1;
			$segregation_score_5B = $segregation_score_5B + 1;
			$segregation_score_6B = $segregation_score_6B + 2;
		}
		
	} # sample_count loop for affecteds


	####################################################
	# Calculate Main homozygous AFFECTED genotype      #
	# HH_count is the count of this main homozygote    #
	####################################################
	$main_affected_homozygous_genotype = "";

	if($AA_count_affected >= $BB_count_affected)
	{
		$HH_count_affected = $AA_count_affected;
		$main_affected_homozygous_genotype = "AA"; # 
	}
	else
	{
		$HH_count_affected = $BB_count_affected;
		$main_affected_homozygous_genotype = "BB"; #
	}


	if (($AA_count_affected + $AB_count_affected + $BB_count_affected) > 0)
	{
		#$affected_homozygosity = $HH_count_affected / ($AA_count_affected + $AB_count_affected + $BB_count_affected);
		$affected_homozygosity = $HH_count_affected / $no_affected_samples; # is this used??
	}


	####################################
	# Calculate main affected genotype #
	####################################
	$main_affected_genotype = "";
	if (($AA_count_affected >= $BB_count_affected) && ($AA_count_affected >= $AB_count_affected))
	{
		$main_affected_genotype = "AA"; #
	}
	if (($AB_count_affected >= $AA_count_affected) && ($AB_count_affected >= $BB_count_affected))
	{
		$main_affected_genotype = "AB"; #
	}
	if (($BB_count_affected >= $AA_count_affected) && ($BB_count_affected >= $AB_count_affected))
	{
		$main_affected_genotype = "BB"; #
	}

} # score_affected_genotypes


#####################################################################
# Count CARRIER genotypes and continue to update segregation scores #
# (This is only called if there are some carriers)                  #
#####################################################################
sub score_carrier_genotypes
{

	for ($sample_count = $no_affected_samples + 1; $sample_count <= $no_affected_samples + $no_carrier_samples; $sample_count++)
	{
		$col_1 = ($sample_count * 2 ) - 1;
		$col_2 = ($sample_count * 2 );
		 
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		 
		$original_genotype = $allele_base_1.$allele_base_2;
		$genotype_AB = "";


		 ###############################################
		 # Carriers:                                   #
		 # $score_X_method can be three things         #
		 # 1. "reference"  (old method)                #
		 # 2. "best"  (what gives highest seg score)   #
		 # 3. "missing" (ignores missing data)         #
		 #     (Then scaled up seg score in proportion #
		 #      to the number of samples used)         #
		 ###############################################
		 if ($original_genotype ne "XX")
		 {
			 #######################################################################################
			 # If allele is main affected allele call it 'A', otherwise call it 'B'                #
			 # This is only used for the segregation scores so only AA, AB and BB make sense       #
			 #######################################################################################
			
			if    (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AA"}
		 	elsif (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "AB"}
		 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AB"}
		 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "BB"}
		 }
		 else
		 {
		 	if ($score_X_method eq "reference")
			 {
				 if ($REF_base eq $main_affected_allele){$genotype_AB = "AA"}else{$genotype_AB = "BB"}
			 }
			 elsif ($score_X_method eq "best")
			 {
				 $genotype_AB = "AB"; # Because these are Carriers the 'best' option would be AB
			 }
			 elsif ($score_X_method eq "missing")
			 {
				 $genotype_AB = "MM";
			 }
		 }
		

		 #######################################################
		 # Carriers:                                           #
		 # Count AAs, ABs and BBs                              #
		 #######################################################
		 if ($genotype_AB eq "AA"){$AA_count_normal = $AA_count_normal + 1}
		 elsif ($genotype_AB eq "AB"){$AB_count_normal = $AB_count_normal + 1}
		 elsif ($genotype_AB eq "BB"){$BB_count_normal = $BB_count_normal + 1}
		 elsif ($genotype_AB eq "MM"){$MM_count = $MM_count + 1}


		 #############
		 # DEBUGGING #
		 #############
		 if ($position eq $position_to_debug)
		 {
		 	print "CARRIER\n";
		 	print "Position      \t$position\tSample: $sample_count\tVariant type: \t$variant_type\n";
		 	print "Allele_base_1:\t$allele_base_1    \tAllele_base_2: $allele_base_2\n\n";

		 	print " allele_1_for_writing_array[$sample_count]: $allele_1_for_writing_array[$sample_count]\n";
		 	print " allele_1_for_writing_array[$sample_count]: $allele_1_for_writing_array[$sample_count]\n\n";

		 	print "score_X_method: $score_X_method\n\n";
		 	print "Original genotype: $original_genotype\tGenotype for scoring: $genotype_AB\n\n";

		 	if ($original_genotype eq "XX")
		 	{
		 		print "NOTE: the setting of 'score_X_method' decides what an 'XX' genotype will be scored as.\n\n";
		 		print "      'reference':  the genotype is set as the reference base ($REF_base) [could be AA or BB depending on the main main_affected_genotype which in this case is $main_affected_genotype]\n";
		 		print "      'best':       the genotype for carriers is set as AB (which will give the highest segregation score))\n";
		 		print "      'missing':    the genotype is set as 'MM' and is not used in the scoring\n";
		 	}
		 	
		 	$answer=<STDIN>;
		 }
		 

		 #######################################################
		 # Count if there is main_affected_homozygous_genotype #
		 #######################################################
		if ($genotype_AB eq $main_affected_homozygous_genotype)
		{
			$HH_count_normal = $HH_count_normal + 1;
		}


		######################################################
		# Carriers: New meaning for seg score 7 - Difference #
		######################################################
		if ($genotype_AB ne $main_affected_genotype)
		{
			$segregation_score_7 = $segregation_score_7 + 1;
			$segregation_score_7B = $segregation_score_7B + 1;
		}

		####################################
		# Segregation scoring for carriers #
		####################################
		# 1 AA_AB_BB         	Recessive:  								Affecteds AA, carriers AB, normals BB
		# 2 AA_AB_ABorBB     	Recessive, some normals may be carriers:  	Affecteds AA, carriers AB, normals AB or BB
		# 3 AA_ABorBB_BB        Recessive, some carriers may be normals:  	Affecteds AA, carriers AB or BB, normals BB
		# 4 AA_ABorBB_ABorBB 	Recessive:  								Affecteds AA, carriers AB or BB, normals AB or BB
		# 5 AAorAB_BB        	Dominant:   								Affecteds AA or AB, normals BB (no such thing as a carrier)
		# 6 A_B              	Additive:   								Affecteds with As, Normals with Bs (no such thing as a carrier)
		# 7 AA_BB               Difference:                                 Normals different to main affected homozygous genotype


		#########################################################################
		# Carriers:                                                             #
		# If there is more than one case OR single_case_scoring is switched OFF #
		#########################################################################
		if ($only_one_case eq "false")
		{
			if ($genotype_AB eq "AA")
			{
				$segregation_score_3B = $segregation_score_3B + 1;
				$segregation_score_4B = $segregation_score_4B + 1;
			}
			elsif ($genotype_AB eq "AB")
			{
				$segregation_score_1 = $segregation_score_1 + 1;
				$segregation_score_2 = $segregation_score_2 + 1;
				$segregation_score_1B = $segregation_score_1B + 1;
				$segregation_score_2B = $segregation_score_2B + 1;

				$segregation_score_3 = $segregation_score_3 + 1;
				$segregation_score_4 = $segregation_score_4 + 1;
				$segregation_score_3B = $segregation_score_3B + 1;
				$segregation_score_4B = $segregation_score_4B + 1;
			}
			elsif ($genotype_AB eq "BB")
			{
				
				######################################################
				# Carriers                                           #
				# An AA genotype could be originally XX              #
				# How should we deal with this?                      #
				# Previously it was set as AA but maybe AB is better #
				######################################################

				####################################
				# Segregation scoring for carriers #
				####################################
				# 1 AA_AB_BB         	Recessive:  								Affecteds AA, carriers AB, normals BB
				# 2 AA_AB_ABorBB     	Recessive, some normals may be carriers:  	Affecteds AA, carriers AB, normals AB or BB
				# 3 AA_ABorBB_BB        Recessive, some carriers may be normals:  	Affecteds AA, carriers AB or BB, normals BB
				# 4 AA_ABorBB_ABorBB 	Recessive:  								Affecteds AA, carriers AB or BB, normals AB or BB
				# 5 AAorAB_BB        	Dominant:   								Affecteds AA or AB, normals BB (no such thing as a carrier)
				# 6 A_B              	Additive:   								Affecteds with As, Normals with Bs (no such thing as a carrier)
				# 7 AA_BB               Difference:                                 Normals different to main affected homozygous genotype

				$segregation_score_3 = $segregation_score_3 + 1;
				$segregation_score_4 = $segregation_score_4 + 1;
			
			}

		} ## NOT single case only


		#########################################################################
		# Carriers:                                                             #
		# If there is only one case and single_case_scoring is switched on      #
		# (See box above for explanation)                                       #
		#########################################################################
		if ($only_one_case eq "true")
		{
			if ($genotype_AB ne $single_case_genotype)
			{
				if ($genotype_AB eq "AA")
				{
					$segregation_score_3B = $segregation_score_3B + 1;
					$segregation_score_4B = $segregation_score_4B + 1;
				}
				elsif ($genotype_AB eq "AB")
				{
					$segregation_score_1 = $segregation_score_1 + 1;
					$segregation_score_2 = $segregation_score_2 + 1;
					$segregation_score_1B = $segregation_score_1B + 1;
					$segregation_score_2B = $segregation_score_2B + 1;

					$segregation_score_3 = $segregation_score_3 + 1;
					$segregation_score_4 = $segregation_score_4 + 1;
					$segregation_score_3B = $segregation_score_3B + 1;
					$segregation_score_4B = $segregation_score_4B + 1;
				}
				elsif ($genotype_AB eq "BB")
				{
					$segregation_score_3 = $segregation_score_3 + 1;
					$segregation_score_4 = $segregation_score_4 + 1;
				}

			} # Only score if Carrier is different from single case

		} # only one case


		#####################################################
		# Seg score 5 (dominant) doesn't apply to carriers  #
		#####################################################
	
	
		################################################################
		# Segregation scoring for Carriers                             #
		################################################################
		# 1 AA_AB_BB         	Recessive:  								Affecteds AA, carriers AB, normals BB
		# 2 AA_AB_ABorBB     	Recessive, some normals may be carriers:  	Affecteds AA, carriers AB, normals AB or BB
		# 3 AA_ABorBB_BB      	Recessive, some carriers may be normals:  	Affecteds AA, carriers AB or BB, normals BB
		# 4 AA_ABorBB_ABorBB 	Recessive:  								Affecteds AA, carriers AB or BB, normals AB or BB
		# 5 AAorAB_BB        	Dominant:   								Affecteds AA or AB, normals BB
		# 6 A_B              	Additive:   								Affecteds A with As, Normals with Bs
		# 7 AA_BB               Difference:                                 Normals different to main affected homozygous genotype

	} # sample_count loop	

} # score_carrier_genotypes


#####################################################################
# Score normal genotypes and continue to update segregation scores  #
# This used to count CARRIERS as well.                              #
#####################################################################
sub score_normal_genotypes
{	 
	for ($sample_count = $no_affected_samples + $no_carrier_samples + 1; $sample_count <= $no_of_samples_without_omits; $sample_count++)
	{
		$col_1 = ($sample_count * 2 ) - 1;
		$col_2 = ($sample_count * 2 );
		 
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];
		 
		 $original_genotype = $allele_base_1.$allele_base_2;
		$genotype_AB = "";


		 ###############################################
		 # Normals:                                    #
		 # $score_X_method can be three things         #
		 # 1. "reference"  (old method)                #
		 # 2. "best"  (what gives highest seg score)   #
		 # 3. "missing" (ignores missing data)         #
		 #     (Then scaled up seg score in proportion #
		 #      to the number of samples used)         #
		 ###############################################
		 if ($original_genotype ne "XX")
		 {
			 #######################################################################################
			 # If allele is main affected allele call it 'A', otherwise call it 'B'                #
			 # This is only used for the segregation scores so only AA, AB and BB make sense       #
			 #######################################################################################
			
			if    (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AA"}
		 	elsif (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "AB"}
		 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AB"}
		 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "BB"}
		 }
		 else
		 {
		 	if ($score_X_method eq "reference")
			 {
				 if ($REF_base eq $main_affected_allele){$genotype_AB = "AA"}else{$genotype_AB = "BB"}
			 }
			 elsif ($score_X_method eq "best")
			 {
				 $genotype_AB = "BB"; # Because these are Carriers the 'best' option would be BB
			 }
			 elsif ($score_X_method eq "missing")
			 {
				 $genotype_AB = "MM";
			 }
		 }


		 ##########################
		 # Normals:               #
		 # Count AAs, ABs and BBs #
		 ##########################
		 if ($genotype_AB eq "AA"){$AA_count_normal = $AA_count_normal + 1}
		 elsif ($genotype_AB eq "AB"){$AB_count_normal = $AB_count_normal + 1}
		 elsif ($genotype_AB eq "BB"){$BB_count_normal = $BB_count_normal + 1}
		 elsif ($genotype_AB eq "MM"){$MM_count = $MM_count + 1}


		 #############
		 # DEBUGGING #
		 #############
		 if ($position eq $position_to_debug)
		 {
		    print "NORMAL\n";
		 	print "Position      \t$position\tSample: $sample_count\tVariant type: \t$variant_type\n";
		 	print "Allele_base_1:\t$allele_base_1    \tAllele_base_2: $allele_base_2\n\n";

		 	print " allele_1_for_writing_array[$sample_count]: $allele_1_for_writing_array[$sample_count]\n";
		 	print " allele_1_for_writing_array[$sample_count]: $allele_1_for_writing_array[$sample_count]\n\n";

		 	print "score_X_method: $score_X_method\n\n";
		 	print "Original genotype: $original_genotype\tGenotype for scoring: $genotype_AB\n\n";

		 	if ($original_genotype eq "XX")
		 	{
		 		print "NOTE: the setting of 'score_X_method' decides what an 'XX' genotype will be scored as.\n\n";
		 		print "      'reference':  the genotype is set as the reference base ($REF_base) [could be AA or BB depending on the main main_affected_genotype which in this case is $main_affected_genotype]\n";
		 		print "      'best':       the genotype for normals is set as BB (which will give the highest segregation score))\n";
		 		print "      'missing':    the genotype is set as 'MM' and is not used in the scoring\n";

		 		if ($score_X_method eq "missing")
		 		{
		 			print "MM_count: $MM_count\n";
		 		}
		 	}
		 	
		 	$answer=<STDIN>;
		 }

		 
		 #######################################################
		 # Count if there is main_affected_homozygous_genotype #
		 #######################################################
		if ($genotype_AB eq $main_affected_homozygous_genotype)
		{
			$HH_count_normal = $HH_count_normal + 1;
		}

	
		################################################################
		# Segregation scoring for Normals                              #
		################################################################
		# 1 AA_AB_BB         	Recessive:  								Affecteds AA, carriers AB, normals BB
		# 2 AA_AB_ABorBB     	Recessive, some normals may be carriers:  	Affecteds AA, carriers AB, normals AB or BB
		# 3 AA_ABorBB_BB      	Recessive, some carriers may be normals:  	Affecteds AA, carriers AB or BB, normals BB
		# 4 AA_ABorBB_ABorBB 	Recessive:  								Affecteds AA, carriers AB or BB, normals AB or BB
		# 5 AAorAB_BB        	Dominant:   								Affecteds AA or AB, normals BB
		# 6 A_B              	Additive:   								Affecteds A with As, Normals with Bs
		# 7 AA_BB               Difference:                                 Normals different to main affected homozygous genotype

		####################################################################################
		# Normals:                                                                        #
		# Scoring system for when there is only one case:                                 #
		#                                                                                 #
		# The single case scores 1 if it matches the required genotype (in this case AA)  #
		#                                                                                 #
		# The Normals score 1 if they meet two requirements:                              #
		#                                                                                 #
		#    1.  They match their required genotype as before                             #
		#    2.  They are different to the single Affected Case                           #
		###################################################################################

	
		
		#####################################################
		# Normals: New meaning for seg score 7 - Difference #
		#####################################################
		if ($genotype_AB ne $main_affected_genotype)
		{
			$segregation_score_7 = $segregation_score_7 + 1;
			$segregation_score_7B = $segregation_score_7B + 1;

			if ($position eq $position_to_debug)
			{
				print "Normals\n";print "Position: $position\n";
				print " Genotype: $genotype_AB\n";
				print " main_affected_homozygous_genotype: $main_affected_homozygous_genotype\n\n";
				print " main_affected_genotype: $main_affected_genotype\n\n";
				print " segregation_score_7: $segregation_score_7\n";
				print " segregation_score_7B: $segregation_score_7B\n";
				print "score_X_method: $score_X_method\n\n";
				$answer=<STDIN>;
			}
		}


		#########################################################################
		# Normals:                                                              #
		# If there is more than one case OR single_case_scoring is switched OFF #
		# Case is AA so Carriers might need to be AB                            #
		#########################################################################
		if ($only_one_case eq "false")
		{
			if ($genotype_AB eq "AA")
			{
				$segregation_score_1B = $segregation_score_1B + 1;
				$segregation_score_2B = $segregation_score_2B + 1;
				$segregation_score_3B = $segregation_score_3B + 1;
				$segregation_score_4B = $segregation_score_4B + 1;	
				$segregation_score_5B = $segregation_score_5B + 1;	
			}
			elsif ($genotype_AB eq "AB")
			{
				$segregation_score_2 = $segregation_score_2 + 1;
				$segregation_score_4 = $segregation_score_4 + 1;
				$segregation_score_2B = $segregation_score_2B + 1;
				$segregation_score_4B = $segregation_score_4B + 1;
			}
			elsif ($genotype_AB eq "BB")
			{
				$segregation_score_1 = $segregation_score_1 + 1;
				$segregation_score_2 = $segregation_score_2 + 1;
				$segregation_score_3 = $segregation_score_3 + 1;
				$segregation_score_4 = $segregation_score_4 + 1;
				$segregation_score_5 = $segregation_score_5 + 1;
			}

			# Normal segregation scoring (additive)
			if ($genotype_AB eq "BB"){$segregation_score_6 = $segregation_score_6 + 2}
			elsif ($genotype_AB eq "AA"){$segregation_score_6B = $segregation_score_6B + 2}
			elsif ($genotype_AB eq "AB"){$segregation_score_6 = $segregation_score_6 + 1; $segregation_score_6B = $segregation_score_6B + 1}

		} # Normal scoring
		
		
		###################################################################################
		# Normals:                                                                        #
		# Scoring system for when there is only one case:                                 #
		#                                                                                 #
		# The single case scores 1 if it matches the required genotype (in this case AA)  #
		#                                                                                 #
		# The Normals score 1 if they meet TWO requirements:                              #
		#                                                                                 #
		#    1.  They match their required genotype as before                             #
		#    2.  They are different to the single Affected Case                           #
		# TRY FOR SEG SCORE 7 FIRST. Make seg score 7 AA/AB/BB
		###################################################################################
		if ($only_one_case eq "true")
		{
			################################################################
			# Segregation scoring for Normals (new with elsif)             #
			################################################################
			# 1 AA_AB_BB         	Recessive:  								Affecteds AA, carriers AB, normals BB
			# 2 AA_AB_ABorBB     	Recessive, some normals may be carriers:  	Affecteds AA, carriers AB, normals AB or BB
			# 3 AA_ABorBB_BB      	Recessive, some carriers may be normals:  	Affecteds AA, carriers AB or BB, normals BB
			# 4 AA_ABorBB_ABorBB 	Recessive:  								Affecteds AA, carriers AB or BB, normals AB or BB
			# 5 AAorAB_BB        	Dominant:   								Affecteds AA or AB, normals BB
			# 6 A_B              	Additive:   								Affecteds A with As, Normals with Bs
			# 7 AA_BB               Strict:                                     Affecteds AA, Normals BB
			if ($genotype_AB ne $single_case_genotype)
			{
				if ($genotype_AB eq "BB")
				{
					$segregation_score_1 = $segregation_score_1 + 1;
					$segregation_score_2 = $segregation_score_2 + 1;
					$segregation_score_3 = $segregation_score_3 + 1;
					$segregation_score_4 = $segregation_score_4 + 1;
					$segregation_score_5 = $segregation_score_5 + 1;
				}
				elsif ($genotype_AB eq "AB")
				{
					$segregation_score_2 = $segregation_score_2 + 1;
					$segregation_score_4 = $segregation_score_4 + 1;
					$segregation_score_2B = $segregation_score_2B + 1;
					$segregation_score_4B = $segregation_score_4B + 1;
				}
				elsif ($genotype_AB eq "AA")
				{
					$segregation_score_1B = $segregation_score_1B + 1;
					$segregation_score_2B = $segregation_score_2B + 1;
					$segregation_score_3B = $segregation_score_3B + 1;
					$segregation_score_4B = $segregation_score_4B + 1;
					$segregation_score_5B = $segregation_score_5B + 1;
				}
				
				#################################################
				# Normals:                                      #
				# Additive segregation scoring                  #
				# (single case doesn't apply here maybe?)    Try with single case for now    #
				#################################################
				if ($genotype_AB eq "BB"){$segregation_score_6 = $segregation_score_6 + 2}
				elsif ($genotype_AB eq "AA"){$segregation_score_6B = $segregation_score_6B + 2}
				elsif ($genotype_AB eq "AB"){$segregation_score_6 = $segregation_score_6 + 1; $segregation_score_6B = $segregation_score_6B + 1}
	
			} # Only score if Normal is different from single case

		} # Single case scoring

	} # score_normal_genotypes


	#################################
	# Calculate normal homozygosity #
	#################################
	if (($AA_count_normal + $AB_count_normal + $BB_count_normal) > 0)
	{
		$normal_homozygosity = $HH_count_normal / ($no_normal_samples + $no_carrier_samples);
	}
	

	################################
	# Calculate homozygosity ratio #
	################################
	if ($normal_homozygosity > 0)
	{
		$homozygosity_ratio = $affected_homozygosity/$normal_homozygosity;
	}
	else
	{
		$homozygosity_ratio = $affected_homozygosity/0.5;
	}

	if ($position eq $position_to_debug)
	{
		&show_segregation_scores;
		&pause;

	}
} # score_normal_genotypes


#333
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

		if    (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AA"}
	 	elsif (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "AB"}
	 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AB"}
	 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "BB"}

		 print "\t$genotype_AB";
	}
	print "\n\tNormals: ";
	for ($check_count = $no_affected_samples + 1; $check_count <= $no_of_samples; $check_count++)
	{
		$col_1 = ($check_count * 2 ) - 1;
		$col_2 = ($check_count * 2 );
		$allele_base_1 = $base_array[$col_1];
		$allele_base_2 = $base_array[$col_2];

		if    (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AA"}
	 	elsif (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "AB"}
	 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AB"}
	 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "BB"}

		 print "\t$genotype_AB";
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
			print "\tMain affected:       \t$main_affected_allele\n";
			print "\tMain normal:         \t$main_normal_allele\n";
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
	
	print "Seg score type      \tSEG\tSEG B\tBEST\n";
	print "1: AA_AB_BB         \t$segregation_score_1\t$segregation_score_1B\t$segregation_score_best_1\n";
	print "2: AA_AB_ABorBB     \t$segregation_score_2\t$segregation_score_2B\t$segregation_score_best_2\n";
	print "3: AA_ABorBB_BB     \t$segregation_score_3\t$segregation_score_3B\t$segregation_score_best_3\n";
	print "4: AA_ABorBB_ABorBB \t$segregation_score_4\t$segregation_score_4B\t$segregation_score_best_4\n";
	print "5: AAorAB_BB        \t$segregation_score_5\t$segregation_score_5B\t$segregation_score_best_5\n";
	print "6: A_B              \t$segregation_score_6\t$segregation_score_6B\t$segregation_score_best_6\n";
	print "7: AA_BB            \t$segregation_score_7\t$segregation_score_7B\t$segregation_score_best_7\n";

	print "Best overall: $segregation_score_best_overall\n\n";

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



###############################################
# Both VEP and snpEff use the same scores now #
###############################################
sub get_effect_score_vep_and_snpEff
{
	my $_effect = $_[0];
	my $_effect_score = 0;

	# Effect score 5 - judged as maximum effect
	if($_effect eq "chromosome_number_variation"){$_effect_score = 5}
	elsif($_effect eq "exon_loss_variant"){$_effect_score = 5}
	elsif($_effect eq "frameshift_variant"){$_effect_score = 5}
	elsif($_effect eq "rare_amino_acid_variant"){$_effect_score = 5}
	elsif($_effect eq "start_lost"){$_effect_score = 5}
	elsif($_effect eq "stop_gained"){$_effect_score = 5}
	elsif($_effect eq "stop_lost"){$_effect_score = 5}
	elsif($_effect eq "exon_loss"){$_effect_score = 5}
	elsif($_effect eq "exon_loss_variant"){$_effect_score = 5}
	elsif($_effect eq "disruptive_inframe_deletion"){$_effect_score = 5}
	elsif($_effect eq "disruptive_inframe_insertion"){$_effect_score = 5}
	elsif($_effect eq "inframe_deletion"){$_effect_score = 5}
	elsif($_effect eq "inframe_insertion"){$_effect_score = 5}
	elsif($_effect eq "missense_variant"){$_effect_score = 5}
	elsif($_effect eq "protein_altering_variant"){$_effect_score = 5} # Not sure what this is exactly but it seems best to give it 5

	elsif($_effect eq "splice_acceptor_variant"){$_effect_score = 5}                        # Promoted from 4 to 5 by Oliver
	elsif($_effect eq "splice_donor_variant"){$_effect_score = 5}                           # Promoted from 4 to 5 by Oliver
	elsif($_effect eq "transcript_ablation"){$_effect_score = 5}                            # Promoted from 4 to 5 by Oliver
	elsif($_effect eq "5_prime_UTR_premature_start_codon_gain_variant"){$_effect_score = 5} # Promoted from 4 to 5 by Oliver
	elsif($_effect eq "initiator_codon_variant"){$_effect_score = 5}                        # Promoted from 4 to 5 by Oliver
	elsif($_effect eq "initiator_codon_variant+non_canonical_start_codon"){$_effect_score = 5}                        # 5?


	# Effect score 4 - might easily have an effect
	elsif($_effect eq "3_prime_UTR_truncation"){$_effect_score = 4}
	elsif($_effect eq "5_prime_UTR_truncation"){$_effect_score = 4}
	elsif($_effect eq "regulatory_region_ablation"){$_effect_score = 4}
	elsif($_effect eq "TFBS_ablation"){$_effect_score = 4}
	elsif($_effect eq "non_coding_exon_variant"){$_effect_score = 4}
	elsif($_effect eq "splice_region_variant"){$_effect_score = 4}                          # Promoted from 3 to 4 by Mike


	# Effect score 3 - could have some effect
	elsif($_effect eq "coding_sequence_variant"){$_effect_score = 3}
	elsif($_effect eq "conserved_intergenic_variant"){$_effect_score = 3}
	elsif($_effect eq "conserved_intron_variant"){$_effect_score = 3}
	elsif($_effect eq "nc_transcript_variant"){$_effect_score = 3}


	# Effect score 2 - upstream or downstream so could effect regulatory elements
	elsif($_effect eq "start_retained"){$_effect_score = 2}
	elsif($_effect eq "stop_retained_variant"){$_effect_score = 2}
	elsif($_effect eq "exon_variant"){$_effect_score = 2}
	elsif($_effect eq "feature_elongation"){$_effect_score = 2}
	elsif($_effect eq "feature_truncation"){$_effect_score = 2}
	elsif($_effect eq "gene_variant"){$_effect_score = 2}
	elsif($_effect eq "intron_variant"){$_effect_score = 2}
	elsif($_effect eq "mature_miRNA_variant"){$_effect_score = 2}
	elsif($_effect eq "miRNA"){$_effect_score = 2}
	elsif($_effect eq "NMD_transcript_variant"){$_effect_score = 2}
	elsif($_effect eq "non_coding_transcript_exon_variant"){$_effect_score = 2}
	elsif($_effect eq "non_coding_transcript_variant"){$_effect_score = 2}
	elsif($_effect eq "regulatory_region_amplelsification"){$_effect_score = 2}
	elsif($_effect eq "regulatory_region_variant"){$_effect_score = 2}
	elsif($_effect eq "TF_binding_site_variant"){$_effect_score = 2}
	elsif($_effect eq "TFBS_amplelsification"){$_effect_score = 2}
	elsif($_effect eq "transcript_amplelsification"){$_effect_score = 2}
	elsif($_effect eq "transcript_variant"){$_effect_score = 2}
	elsif($_effect eq "upstream_gene_variant"){$_effect_score = 2}

	elsif($_effect eq "3_prime_UTR_variant"){$_effect_score = 2}                          # Promoted from 1 to 2 by Mike
	elsif($_effect eq "5_prime_UTR_variant"){$_effect_score = 2}                          # Promoted from 1 to 2 by Mike

	# Effect score 1 - unlikely to have an effect
	elsif($_effect eq "synonymous_variant"){$_effect_score = 1}
	elsif($_effect eq "downstream_gene_variant"){$_effect_score = 1}
	elsif($_effect eq "intergenic_region"){$_effect_score = 1}
	elsif($_effect eq "intragenic_variant"){$_effect_score = 1}
	elsif($_effect eq "coding_sequence_variant"){$_effect_score = 1}
	elsif($_effect eq "intergenic_variant"){$_effect_score = 1}

	# Sometimes comes up (with indels?)
	elsif($_effect eq "transcript"){$_effect_score = 1}


	#########################################################
	# If no effects score is found then write to a file so  #
	# that user can check up on why this happened.          #
	# (If the chr is Unknown don't do this)                 #
	#########################################################
	if (($_effect_score == 0) && (index($chromosome,"Un") == -1) && ($chromosome ne "M"))
	{
		print MISSING_EFFECTS "Chr: $chromosome\tPos: $position\t$_effect\t$effect_predictor\n";
		$_effect_score = 4.4;

		$missing_effects_count = $missing_effects_count + 1;
	} #


	###############################################
	# Both VEP and snpEff use the same scores now #
	# but we will still use the same separate     #
	# variables as it may affect code furhter on. #
	###############################################
	$vep_effect_score = $_effect_score;
	$snpEff_effect_score = $_effect_score;

	# For use as a function
	$_effect_score = $_effect_score;

	######################################
	#  These are snpEff's impact ratings #
	######################################

	#HIGH	chromosome_number_variation
	#HIGH	exon_loss_variant
	#HIGH	frameshift_variant
	#HIGH	rare_amino_acid_variant
	#HIGH	splice_acceptor_variant
	#HIGH	splice_donor_variant
	#HIGH	start_lost
	#HIGH	stop_gained
	#HIGH	stop_lost
	#HIGH	transcript_ablation

	#MODERATE         3_prime_UTR_truncation+exon_loss
	#MODERATE         5_prime_UTR_truncation+exon_loss_variant
	#MODERATE         coding_sequence_variant
	#MODERATE         disruptive_inframe_deletion
	#MODERATE         disruptive_inframe_insertion
	#MODERATE         inframe_deletion
	#MODERATE         inframe_insertion
	#MODERATE         missense_variant                
	#MODERATE         regulatory_region_ablation
	#MODERATE         splice_region_variant
	#MODERATE         TFBS_ablation

	#LOW              5_prime_UTR_premature_start_codon_gain_variant
	#LOW              initiator_codon_variant
	#LOW              splice_region_variant
	#LOW              start_retained
	#LOW              stop_retained_variant
	#LOW              synonymous_variant

	#MODIFIER         3_prime_UTR_variant
	#MODIFIER         5_prime_UTR_variant
	#MODIFIER         coding_sequence_variant
	#MODIFIER         conserved_intergenic_variant
	#MODIFIER         conserved_intron_variant
	#MODIFIER         downstream_gene_variant
	#MODIFIER         exon_variant
	#MODIFIER         feature_elongation
	#MODIFIER         feature_truncation
	#MODIFIER         gene_variant
	#MODIFIER         intergenic_region
	#MODIFIER         intragenic_variant
	#MODIFIER         intron_variant
	#MODIFIER         mature_miRNA_variant
	#MODIFIER         miRNA
	#MODIFIER         NMD_transcript_variant
	#MODIFIER         non_coding_transcript_exon_variant
	#MODIFIER         non_coding_transcript_variant
	#MODIFIER         regulatory_region_amplification
	#MODIFIER         regulatory_region_variant
	#MODIFIER         TF_binding_site_variant
	#MODIFIER         TFBS_amplification
	#MODIFIER         transcript_amplification
	#MODIFIER         transcript_variant
	#MODIFIER         upstream_gene_variant
}



sub get_sample_cols
{
	my $sample_count_sub = $_[0];

	$col_1 = ($sample_count_sub * 2) - 1;
	$col_2 = ($sample_count_sub * 2);
}

sub choose_best_segregation_scores_indels_merged
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


	#############################################################
	# If segregation score is greater than number of samples    #
	# then cut it down to no of samples (e.g. 14.1 ==> 14)      #
	#############################################################
	$highest_possible_segregation_score = $no_of_samples_without_omits;


	############################################################################################
	# Normalise additive score (segregation_score_best_6) so it has the same maximum as others #
	############################################################################################

	$segregation_score_best_6 = ($segregation_score_best_6 * $no_of_samples_without_omits)/ ($no_affected_samples * 2 + $no_normal_samples * 2 + $no_carrier_samples);

	$segregation_score_best_6 = sprintf("%.1f", $segregation_score_best_6);

	################################################################################################
	# Find best overall segregation score (but miss out seg_score_7 as that is a bit experimental) #
	################################################################################################
	$segregation_score_best_overall= max($segregation_score_best_1, $segregation_score_best_2, $segregation_score_best_3, $segregation_score_best_4, $segregation_score_best_5, $segregation_score_best_6);

} # choose_best_segregation_scores_indels_merged

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

	$segregation_score_best_overall= 0;
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

sub pause_with_cancel
{
	if ($pause_cancelled eq "false")
	{
		print "\n Press RETURN to continue.\n\n";
		print "Type 'u' if you want the pause at this point to be ignored, and the programm to continue uninterrupted\n";
		$answer=<STDIN>;
		chomp $answer;
		if (lc $answer eq "u")
		{
			$pause_cancelled = "true";
		}
	}
} # pause_with_cancel


################################################
# Get name of input VCF file                   #
################################################
sub get_input_file
{
	my $question = "";
	my $inputfile = "";
	my $filetype = "";	

	$question = $_[0];
	$filetype = $_[1];

	&print_message("$question","input");

	until (-e $inputfile)
	{
		print "   File:  ";
		$inputfile = <STDIN>;
		chomp $inputfile;
		if ($inputfile eq ""){$inputfile = "ls"}
		if ($inputfile eq "ls"){print "\n";system ("ls *.vcf");print "\n"}
		if ($inputfile ne "ls"){if (! -e $inputfile){print "\n\n>>>>>>>>  File $inputfile not found.  Try again.  <<<<<<<<\n\n";}}
	}
	$inputfile = $inputfile;
}



 ######   ######## ########         ######## ######## ######## ########  ######  ########         ########  ########  ######## ########  ####  ######  ########  #######  ########  
##    ##  ##          ##            ##       ##       ##       ##       ##    ##    ##            ##     ## ##     ## ##       ##     ##  ##  ##    ##    ##    ##     ## ##     ## 
##        ##          ##            ##       ##       ##       ##       ##          ##            ##     ## ##     ## ##       ##     ##  ##  ##          ##    ##     ## ##     ## 
##   #### ######      ##            ######   ######   ######   ######   ##          ##            ########  ########  ######   ##     ##  ##  ##          ##    ##     ## ########  
##    ##  ##          ##            ##       ##       ##       ##       ##          ##            ##        ##   ##   ##       ##     ##  ##  ##          ##    ##     ## ##   ##   
##    ##  ##          ##            ##       ##       ##       ##       ##    ##    ##            ##        ##    ##  ##       ##     ##  ##  ##    ##    ##    ##     ## ##    ##  
 ######   ########    ##    ####### ######## ##       ##       ########  ######     ##    ####### ##        ##     ## ######## ########  ####  ######     ##     #######  ##     ## 

################################################
# Ask user how the SNP Effects are to be added #
################################################
sub get_effect_predictor
{

	until (($effect_predictor ne "")  && ($effect_predictor ne "both"))
	{
		&print_message("How are the consequences/effects of the mutations to be added?","input");

		if ($effect_predictor eq "VEP")
		{
			print "   <1> snpEff                                \n";
			print "   <2> Variant Effect Predictor (VEP)          [DEFAULT]\n\n";

			$answer = <STDIN>;
			chomp $answer;

			if ($answer eq "" ){$effect_predictor = "VEP"} # default

			if (substr($answer,0,1) eq "1" ){$effect_predictor = "snpEff"}
			if (substr($answer,0,1) eq "2" ){$effect_predictor = "VEP"}
		}
		else
		{
			print "   <1> snpEff                                  [DEFAULT]\n";
			print "   <2> Variant Effect Predictor (VEP)\n\n";

			$answer = <STDIN>;
			chomp $answer;

			if ($answer eq "" ){$effect_predictor = "snpEff"} # default

			if (substr($answer,0,1) eq "1" ){$effect_predictor = "snpEff"}
			if (substr($answer,0,1) eq "2" ){$effect_predictor = "VEP"}
		}
		
		print "\n\n";
	}
}




###################################################
# This sets a load of things to zero at the start #
###################################################
sub set_all_to_zero
{

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

	$MM_count = 0; # count of missing data
	
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
	
	$segregation_score_best_overall= 0;  # New Aug 22 2015

} # set_all_to_zero


###################################################################
# This checks the first line to see that it looks like a VCF file #
# (Not very sophisticated. May need improving)                    #
###################################################################
sub check_vcf_format
{
	if($line_count == 1)
	{
		&print_message(" File format check  ","message");
		if (index($vcf_line,"VCF") > -1)
		{
			print "   First line is $vcf_line.  This looks like a VCF file\n\n";
		}	
		else
		{
			print "   First line is $vcf_line.  This doesn't look like a VCF file\n\n";
			exit;
		}	
	}
	
} # check_vcf_format

#################################################
# Reads the first seven columns of the VCF file #
# myArray1 gets the whole line split at TABS    #
#################################################
sub read_vcf_first_seven_columns
{
	#########################################################
    # Split whole line at TABs into an array myArray1 (1-9) #
    #########################################################
	
	@myArray1 = split (/\t/,$vcf_line);
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
		print "$vcf_line\n\n";
		print "You should check the VCF file.\n\n";

		print "If you want to continue anyway, press 'return'   \n\n";  $answer = <STDIN>;
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

    # What is this?
	#if ($position =~ /^\d+?$/) {$vcf_variant_count = $vcf_variant_count + 1}
	
	$vcf_variant_count = $vcf_variant_count + 1;

    ####################################################################
	# Remove string 'chr'                                              #
	####################################################################
    if (substr($chromosome,0,3) eq "chr"){$chromosome = substr($chromosome,3,99)}


	#################################################################
	# Convert chrX to relevant number (e.g. chrX --> chr39 for dog) #
	#################################################################
	if ($chromosome eq "X")
	{
		$chromosome = $x_chromosome_number;
	}
	
	####################################################################
    # Make joint chromosome and position to check if it is on SNP list #
    ####################################################################
    $chr_pos = $chromosome."_".$position;

    if (defined $snp_list_hash{$chr_pos})
	{
		#$on_snp_list = "true";
		$on_snp_list = $snp_list_hash{$chr_pos};
	}
    elsif (defined $snp_list_hash{"chr".$chr_pos})
	{
		$on_snp_list = $snp_list_hash{"chr".$chr_pos};
	}
	else
	{
		$on_snp_list = "false";
	}

} # read_vcf_first_seven_columns





###########################################################
# Choose best segregation score out of segregation_score  #
# and segregation_score_B.  This removes the problem that #
# the scores depends on how the "main affected allele" is #
# chosen.                                                 #
###########################################################
sub choose_best_segregation_score_A_or_B
{
	$swapped = "false";
	$segregation_score_best_1 = 0;
	$segregation_score_best_2 = 0;
	$segregation_score_best_3 = 0;
	$segregation_score_best_4 = 0;
	$segregation_score_best_5 = 0;
	$segregation_score_best_6 = 0;
	$segregation_score_best_7 = 0;
	$segregation_score_best_overall  = 0;
	
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

	$segregation_score_best_6 = ($segregation_score_best_6 * $no_of_samples_without_omits)/ ($no_affected_samples * 2 + $no_normal_samples * 2 + $no_carrier_samples);

	$segregation_score_best_6 = sprintf("%.1f", $segregation_score_best_6);

	

	########################################################################################
	# If score_X_method = "missing" then you have to normalise the score for each variant. #
	# So for example if there are 14 samples, and 3 of them are counted as missing,        #
	# then the segregation score for that variant would be multiplied by 14/11             #
	# to bring it up to the value it would be if all samples were counted                  #
	########################################################################################

	if ($score_X_method eq "missing")
	{
		if ($MM_count < $no_of_samples_without_omits)
		{
			$normalise_missing_factor = $no_of_samples_without_omits / ($no_of_samples_without_omits - $MM_count);

			if ($position eq $position_to_debug)
			{
				print "no_of_samples_without_omits: $ no_of_samples_without_omits\n";
				print "MM_count: $ MM_count\n";
				print "normalise_missing_factor: $normalise_missing_factor\n\n";
				print "Before normalisation - segregation_score_best_1: $segregation_score_best_1\n\n";

			}
			
			# Segregation scores
			$segregation_score_best_1 =  $segregation_score_best_1 * $normalise_missing_factor;
			$segregation_score_best_2 =  $segregation_score_best_2 * $normalise_missing_factor;
			$segregation_score_best_3 =  $segregation_score_best_3 * $normalise_missing_factor;
			$segregation_score_best_4 =  $segregation_score_best_4 * $normalise_missing_factor;
			$segregation_score_best_5 =  $segregation_score_best_5 * $normalise_missing_factor;
			$segregation_score_best_6 =  $segregation_score_best_6 * $normalise_missing_factor;
			$segregation_score_best_7 =  $segregation_score_best_7 * $normalise_missing_factor;

			# Format for 2 decimal places
			$segregation_score_best_1 = sprintf("%.1f", $segregation_score_best_1);
			$segregation_score_best_2 = sprintf("%.1f", $segregation_score_best_2);
			$segregation_score_best_3 = sprintf("%.1f", $segregation_score_best_3);
			$segregation_score_best_4 = sprintf("%.1f", $segregation_score_best_4);
			$segregation_score_best_5 = sprintf("%.1f", $segregation_score_best_5);
			$segregation_score_best_6 = sprintf("%.1f", $segregation_score_best_6);
			$segregation_score_best_7 = sprintf("%.1f", $segregation_score_best_7);

			if ($position eq $position_to_debug)
			{
				print "After normalisation - segregation_score_best_1: $segregation_score_best_1\n\n";

				$answer=<STDIN>;
			}
			#carrot
		}
	} # score_X_method = missing

	################################################################################################
	# Find best overall segregation score (but miss out seg_score_7 as that is a bit experimental) #
	################################################################################################
	$segregation_score_best_overall= max($segregation_score_best_1, $segregation_score_best_2, $segregation_score_best_3, $segregation_score_best_4, $segregation_score_best_5, $segregation_score_best_6);
	
	$segregation_score_to_filter_by = $segregation_score_best_overall;

} # choose_best_segregation_score_A_or_B


############################################################################
# Merge Indels                                                             #
# This procedure merges nearby indels into one, in case they are actually  #
# the same indels by the variant caller has decided thay are in a slightly #
# different position.                                                      #
############################################################################
sub merge_indels
{
	########################################################################################
	# Open INDELS_ONLY ONLY file which was saved earlier as $simplified_indels_output_file #
	########################################################################################

	&print_message("Merging Indels into $indels_merged_file","message");

	open (INDELS_OUTPUT, "$simplified_indels_output_file") || die "Cannot open indels file: $simplified_indels_output_file";

	@indel_file_array= <INDELS_OUTPUT>;

	close INDELS_OUTPUT;

	$no_of_indels = scalar @indel_file_array;

	####################################
	# Open out file for merged indels  #
	####################################
	open (INDELS_MERGED, ">$indels_merged_file") || die "Cannot open indels file: $indels_merged_file";


	###############################################
	# Now add headers for the indels merged file  #
	###############################################
	print INDELS_MERGED "Chr\tPos\tRef";

	for ($sample_count = 1;$sample_count <=$no_of_samples_without_omits; $sample_count++)
	{
	    ##########################################################
	    # Add name of sample as header in first batch of columns #
	    # (ie those after the first three CHR,POS,REF columns    #
	    ##########################################################
	    
	    $sample_name = &get_prefix($sample_name_array[$sample_count]);

		print INDELS_MERGED "\t$sample_name_array[$sample_count]"; # Once for allele_1 
		print INDELS_MERGED "\t$sample_name_array[$sample_count]"; # Once for allele_2
	} # Sample_count loop

	print INDELS_MERGED "\tSeg score 1 AA/AB/BB";
	print INDELS_MERGED "\tSeg score 2 AA/AB/ABorBB";
	print INDELS_MERGED "\tSeg score 3 AA/ABorBB/BB";
	print INDELS_MERGED "\tSeg score 4 AA/ABorBB/ABorBB";
	print INDELS_MERGED "\tSeg score 5 AAorAB/BB dom";
	print INDELS_MERGED "\tSeg score 6 A/B add";
	print INDELS_MERGED "\tSeg score 7 Aff <> Norm";
	print INDELS_MERGED "\tSeg score BEST";
	print INDELS_MERGED "\tMain affected allele is REF";
	print INDELS_MERGED "\tEffect score";
	print INDELS_MERGED "\tConsequence";
	print INDELS_MERGED "\tAlt";
	print INDELS_MERGED "\tallele_A";
	print INDELS_MERGED "\tallele_B";
	print INDELS_MERGED "\tMain_affected_allele";
	print INDELS_MERGED "\tHomozygosity Ratio";
	print INDELS_MERGED "\tNo of ALT alleles";
	print INDELS_MERGED "\tVariant Caller";
	print INDELS_MERGED "\tDepth\n";


	for ($line_count = 1; $line_count < $no_of_indels; $line_count++)
	{

		if (($line_count % 20000) == 0)
		{
			if ($passed_header_lines eq "true"){print "Reading simplified indels file.  Line: $line_count\n"}
		}

		$vcf_line = $indel_file_array[$line_count];

		chomp $vcf_line;

		#Split line at tabs
		@item = split (/\t/,$vcf_line);
		$array_size_1 = scalar @item;

		#################################################################
		# Indels Merging                                                #
		# $array_size_1 is the number of data fields when split by tabs #
		# $no_of_data_columns is the no_of_samples x 2                  #
		#################################################################
		$no_of_data_columns = $no_of_samples_without_omits * 2;
		$chromosome = $item[0];
	    $position = $item[1];
	    $REF_base = $item[3];

	    $main_affected_allele =  $item[$no_of_data_columns + 17];

		for ($sample_count = 1; $sample_count <= $no_of_samples_without_omits; $sample_count++)
	    {
	    	$col_1 = ($sample_count * 2) - 1;
			$col_2 = ($sample_count * 2);
			$base_array[$col_1] = $item[$col_1 + 2];
			$base_array[$col_2] = $item[$col_2 + 2];
	    }


		###########################################################
	    # Read the columns to the right of the segregation scores #
	    ###########################################################
	    $effect_score = 0;
	    $consequence = "";

	    $effect_score = $item[$no_of_data_columns + 12]; #This use of 12 etc could be improved
	    $consequence = $item[$no_of_data_columns + 13];
	    $homozygosity_ratio = $item[$no_of_data_columns + 18];
	    $no_of_alleles = $item[$no_of_data_columns + 19];

		if ($effect_score eq "") {$effect_score = 0}


		####################################################################
		# Indels merging                                                   #
	    # Check forward several rows to see if there is an indel(s) nearby #
	    ####################################################################
	    for ($check_count = 1; $check_count <= $merge_threshold; $check_count++)
	    {
	    	#Only do this if you don't go beyond the end of the file #
	    	if (($line_count + $check_count) < $no_of_indels)
	    	{
		    	$check_line = $indel_file_array[$line_count + $check_count];
				@item_check = split (/\t/,$check_line);
				$position_check = $item_check[1];
			}

			
			################################################################
			# Indels merging                                               #
			# Only proceed if this line has not be used in a merge already #
			################################################################
			if (not defined $merged_hash{$line_count + $check_count})
			{
				#print "Position: $position\tPosition check: $position_check\n";

				if (($position_check - $position) < $merge_threshold)
				{
					##########################################
					# Merge this line with the original line #
					##########################################

					# Mark as merged 
					#$merged_hash{$line_count + $check_count} = "M";

					# Get the alleles from this check position
					for ($sample_count = 1; $sample_count <= $no_of_samples_without_omits; $sample_count++)
				    {
				    	$col_1 = ($sample_count * 2) - 1;
						$col_2 = ($sample_count * 2);
						$base_array_check[$col_1] = $item_check[$col_1 + 2];
						$base_array_check[$col_2] = $item_check[$col_2 + 2];
				    }

				    #################################################################
					# Indels merging                                                #
				    # Read the columns to the right of the segregation scores check #
				    #################################################################
				    $effect_score_check = 0;
				    $consequence_check = "";

					$effect_score_check = $item_check[$no_of_data_columns + 12];
				    $consequence_check = $item_check[$no_of_data_columns + 13];
				    $homozygosity_ratio_check = $item_check[$no_of_data_columns + 18];
				    $no_of_alleles_check = $item_check[$no_of_data_columns + 19];



				    if ($effect_score_check eq "") {$effect_score_check = 0}

			    	# Get the Main Affected Allele (check)
				    $main_affected_allele_check =  $item_check[$no_of_data_columns + 17];

				    ##########################################################################################################
					# Indels merging                                                                                         #
					# Check across columns for alleles and merge alleles according to the following rules:                   #
					# In the cases if one of the pair is the main affected allele (MAA) make merged the main affected allele #
					# In the controls if one of pair is NOT the MAA then make the merged allele NOT the MAA                  #
					##########################################################################################################


				    ##############################################
					# Indels merging                             #
					# Check through affected samples BY GENOTYPE #
					##############################################
					for ($sample_count = 1; $sample_count <= $no_affected_samples; $sample_count++)
				    {
				    	$col_1 = ($sample_count * 2) - 1;
						$col_2 = ($sample_count * 2);
						$allele_base_1 = $base_array[$col_1];
						$allele_base_2 = $base_array[$col_2];
						$allele_base_check_1 = $base_array_check[$col_1];
						$allele_base_check_2 = $base_array_check[$col_2];

						#######################################################################################
						# NEW METHOD: if allele is main affected allele call it 'A', otherwise call it 'B'    #
						# This is only used for the segregation scores so only AA, AB and BB make sense       #
						#######################################################################################
						
						if    (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AA"}
					 	elsif (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "AB"}
					 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AB"}
					 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "BB"}

					 	if    (($allele_base_check_1 eq $main_affected_allele_check) && ($allele_base_check_2 eq $main_affected_allele_check)){$genotype_check = "AA"}
					 	elsif (($allele_base_check_1 eq $main_affected_allele_check) && ($allele_base_check_2 ne $main_affected_allele_check)){$genotype_check = "AB"}
					 	elsif (($allele_base_check_1 ne $main_affected_allele_check) && ($allele_base_check_2 eq $main_affected_allele_check)){$genotype_check = "AB"}
					 	elsif (($allele_base_check_1 ne $main_affected_allele_check) && ($allele_base_check_2 ne $main_affected_allele_check)){$genotype_check = "BB"}

						#########################################################################
					 	# Merged the genotypes, choosing the one with most As for the Affecteds #
					 	#########################################################################
				    	if (($genotype_AB eq "AA") || ($genotype_check eq "AA")){$merged_genotype = "AA"}

				    	if (($genotype_AB eq "AB") && ($genotype_check eq "BB")){$merged_genotype = "AB"}
				    	if (($genotype_AB eq "BB") && ($genotype_check eq "AB")){$merged_genotype = "AB"}

				    	if (($genotype_AB eq "BB") && ($genotype_check eq "BB")){$merged_genotype = "BB"}

				    	# Store in merged genotype array
				    	$merged_genotype_array[$sample_count] = $merged_genotype;

				    	
				    	# Check whether genotype is 2 bases long
				    	if (length ($merged_genotype) == 2)
				    	{
				    		$allele_base_merged_1 = substr($merged_genotype,0,1);
					    	$allele_base_merged_2 = substr($merged_genotype,1,1);

					    	# Store in merged allele array
					    	$merged_allele_array[$col_1] = $allele_base_merged_1;
					    	$merged_allele_array[$col_2] = $allele_base_merged_2;

				    	}

				    	if (length ($merged_genotype) < 2)
				    	{
				    		print "Sample $sample_count\tLine: $line_count\tPos: $position\tmerged genotype: $merged_genotype\tLength: ".length ($merged_genotype)."\n";
				    		# Store in merged allele array
					    	$merged_allele_array[$col_1] = "y";
					    	$merged_allele_array[$col_2] = "y";
				    	}
				    	
				    }


				    ##############################################
					# Indels merging                             #
					# Check through Normal samples BY GENOTYPE   #
					##############################################
					for ($sample_count = $no_affected_samples + 1; $sample_count <= $no_of_samples_without_omits; $sample_count++)
				    {
				    	$col_1 = ($sample_count * 2) - 1;
						$col_2 = ($sample_count * 2);
						$allele_base_1 = $base_array[$col_1];
						$allele_base_2 = $base_array[$col_2];
						$allele_base_check_1 = $base_array_check[$col_1];
						$allele_base_check_2 = $base_array_check[$col_2];

						#######################################################################################
						# NEW METHOD: if allele is main affected allele call it 'A', otherwise call it 'B'    #
						# This is only used for the segregation scores so only AA, AB and BB make sense       #
						#######################################################################################
						
						if    (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AA"}
					 	elsif (($allele_base_1 eq $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "AB"}
					 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 eq $main_affected_allele)){$genotype_AB = "AB"}
					 	elsif (($allele_base_1 ne $main_affected_allele) && ($allele_base_2 ne $main_affected_allele)){$genotype_AB = "BB"}

					 	if    (($allele_base_check_1 eq $main_affected_allele_check) && ($allele_base_check_2 eq $main_affected_allele_check)){$genotype_check = "AA"}
					 	elsif (($allele_base_check_1 eq $main_affected_allele_check) && ($allele_base_check_2 ne $main_affected_allele_check)){$genotype_check = "AB"}
					 	elsif (($allele_base_check_1 ne $main_affected_allele_check) && ($allele_base_check_2 eq $main_affected_allele_check)){$genotype_check = "AB"}
					 	elsif (($allele_base_check_1 ne $main_affected_allele_check) && ($allele_base_check_2 ne $main_affected_allele_check)){$genotype_check = "BB"}


						#########################################################################
					 	# Merged the genotypes, choosing the one with most Bs for the Normals   #
					 	#########################################################################
				    	if (($genotype_AB eq "BB") || ($genotype_check eq "BB")){$merged_genotype = "BB"}

				    	if (($genotype_AB eq "AA") && ($genotype_check eq "AB")){$merged_genotype = "AB"}
				    	if (($genotype_AB eq "AB") && ($genotype_check eq "AA")){$merged_genotype = "AB"}

				    	if (($genotype_AB eq "AA") && ($genotype_check eq "AA")){$merged_genotype = "AA"}

				    	# Store in merged genotype array
				    	$merged_genotype_array[$sample_count] = $merged_genotype;

				    	$allele_base_merged_1 = substr($merged_genotype,0,1);
				    	$allele_base_merged_2 = substr($merged_genotype,1,1);

				    	# Store in merged allele array
				    	$merged_allele_array[$col_1] = $allele_base_merged_1;
				    	$merged_allele_array[$col_2] = $allele_base_merged_2;
				    }


				    #####################################################
					# Indels merging                                    #
					# Find main affected allele in the merged genotypes #
					#####################################################
					$A_count_affected = 0;
					$B_count_affected = 0;
					$main_affected_allele = "R";

					for ($sample_count = 1; $sample_count <= $no_affected_samples; $sample_count++)
				    {
				    	$col_1 = ($sample_count * 2) - 1;
						$col_2 = ($sample_count * 2);
						$allele_base_1 = $merged_allele_array[$col_1];
						$allele_base_2 = $merged_allele_array[$col_2];

						if ($allele_base_1 eq "A"){$A_count_affected = $A_count_affected + 1}
						if ($allele_base_2 eq "A"){$A_count_affected = $A_count_affected + 1}
						if ($allele_base_1 eq "B"){$B_count_affected = $B_count_affected + 1}
						if ($allele_base_2 eq "B"){$B_count_affected = $B_count_affected + 1}
				    }

				    if ($B_count_affected > $A_count_affected){$main_affected_allele = "B"} else {$main_affected_allele = "A"}

				   

				    ###################################################
				    # Write first three columns to Indels merged file #
				    ###################################################
					print INDELS_MERGED "$chromosome\t$position\t$REF_base";

				    for ($col_count = 1; $col_count <= $no_of_samples_without_omits * 2; $col_count++)
				    {
				    	print INDELS_MERGED "\t$merged_allele_array[$col_count]";
				    }
				    #print "\n";
				    


				    #######################################
				    # Calculate segregation scores Merged #
				    #######################################
				    &zero_segregation_scores;

				    # Affected
				    for ($sample_count = 1; $sample_count <= $no_affected_samples; $sample_count++)
				    {
				    	$col_1 = ($sample_count * 2) - 1;
						$col_2 = ($sample_count * 2);
						$allele =  $merged_allele_array[$col_1];
				    	$allele_check =  $merged_allele_array[$col_2];
				    	$genotype_AB = $allele.$allele_check;

						################################################################
						# Segregation scoring for Affecteds (Merged Indels)            #
						################################################################
						# 1 AA_AB_BB         	Recessive:  								Affecteds AA, carriers AB, normals BB
						# 2 AA_AB_ABorBB     	Recessive, some normals may be carriers:  	Affecteds AA, carriers AB, normals AB or BB
						# 3 AA_ABorBB_BB      	Recessive, some carriers may be normals:  	Affecteds AA, carriers AB or BB, normals BB
						# 4 AA_ABorBB_ABorBB 	Recessive:  								Affecteds AA, carriers AB or BB, normals AB or BB
						# 5 AAorAB_BB        	Dominant:   								Affecteds AA or AB, normals BB
						# 6 A_B              	Additive:   								Affecteds A with As, Normals with Bs
						# 7 AA_BB               Strict:                                     Affecteds AA, Normals BB

				    	if ($genotype_AB eq "AA") # Merged Indels
						{
							$segregation_score_1 = $segregation_score_1 + 1;
							$segregation_score_2 = $segregation_score_2 + 1;
							$segregation_score_3 = $segregation_score_3 + 1;
							$segregation_score_4 = $segregation_score_4 + 1;
							$segregation_score_6 = $segregation_score_6 + 2;
							$segregation_score_7 = $segregation_score_7 + 1;
						}
						
						if ($genotype_AB eq "BB") # Merged Indels
						{
							$segregation_score_1B = $segregation_score_1B + 1;
							$segregation_score_2B = $segregation_score_2B + 1;
							$segregation_score_3B = $segregation_score_3B + 1;
							$segregation_score_4B = $segregation_score_4B + 1;
							$segregation_score_6B = $segregation_score_6B + 2;
							$segregation_score_7B = $segregation_score_7B + 1;
						}
						
						
						# Affected segregation scoring (dominant) Merged Indels
						if (($genotype_AB eq "AA") || ($genotype_AB eq "AB") || ($genotype_AB eq "BA"))
						{
							$segregation_score_5 = $segregation_score_5 + 1;
						}
						
						# Affected segregation scoring (dominant) Merged Indels
						if (($genotype_AB eq "BB") || ($genotype_AB eq "AB") || ($genotype_AB eq "BA"))
						{
							$segregation_score_5B = $segregation_score_5B + 1;
						}
						# Merged Indels
						if (($genotype_AB eq "AB") || ($genotype_AB eq "BA"))
						{
							$segregation_score_6 = $segregation_score_6 + 1;
							$segregation_score_6B = $segregation_score_6B + 1;
						}
						
					 } # Affected seg scores


				    # Non-affected 
				    # (if any are omitted then they are at the end and not scored ) #
				    for ($sample_count = $no_affected_samples + 1; $sample_count <= $no_of_samples_without_omits; $sample_count++)
				    {
				    	$col_1 = ($sample_count * 2) - 1;
						$col_2 = ($sample_count * 2);
						$allele =  $merged_allele_array[$col_1];
				    	$allele_check =  $merged_allele_array[$col_2];
				    	$genotype_AB = $allele.$allele_check;


				    	################################################################
						# Segregation scoring for Non-Affecteds   (Merged Indels)      #
						################################################################
						# 1 AA_AB_BB         	Recessive:  								Affecteds AA, carriers AB, normals BB
						# 2 AA_AB_ABorBB     	Recessive, some normals may be carriers:  	Affecteds AA, carriers AB, normals AB or BB
						# 3 AA_ABorBB_BB      	Recessive, some carriers may be normals:  	Affecteds AA, carriers AB or BB, normals BB
						# 4 AA_ABorBB_ABorBB 	Recessive:  								Affecteds AA, carriers AB or BB, normals AB or BB
						# 5 AAorAB_BB        	Dominant:   								Affecteds AA or AB, normals BB
						# 6 A_B              	Additive:   								Affecteds A with As, Normals with Bs
						# 7 AA_BB               Strict:                                     Affecteds AA, Normals BB

				    	# Normals

						if ($genotype_AB eq "BB") # Merged Indels
						{
							$segregation_score_1 = $segregation_score_1 + 1;
							$segregation_score_3 = $segregation_score_3 + 1;
							$segregation_score_5 = $segregation_score_5 + 1;
							$segregation_score_7 = $segregation_score_7 + 1;
						}
						if (($genotype_AB eq "AB") || ($genotype_AB eq "BA") || ($genotype_AB eq "BB"))
						{
							$segregation_score_2 = $segregation_score_2 + 1;
							$segregation_score_4 = $segregation_score_4 + 1;
						}
						if ($genotype_AB eq "AA") # Merged Indels
						{
							$segregation_score_1B = $segregation_score_1B + 1;
							$segregation_score_3B = $segregation_score_3B + 1;
							$segregation_score_5B = $segregation_score_5B + 1;
							$segregation_score_7B = $segregation_score_7B + 1;
						}
						if (($genotype_AB eq "AA") || ($genotype_AB eq "AB") || ($genotype_AB eq "BA"))
						{
							$segregation_score_2B = $segregation_score_2B + 1;
							$segregation_score_4B = $segregation_score_4B + 1;
						}
						
						# Normal segregation scoring (additive)
						if ($genotype_AB eq "BB"){$segregation_score_6 = $segregation_score_6 + 2}
						elsif ($genotype_AB eq "AA"){$segregation_score_6B = $segregation_score_6B + 2}
						elsif (($genotype_AB eq "AB")  || ($genotype_AB eq "BA")){$segregation_score_6 = $segregation_score_6 + 1; $segregation_score_6B = $segregation_score_6B + 1}
					
				    } # Merged Indels

					###########################################################
					# Choose best segregation score out of segregation_score  #
					# and segregation_score_B.  This removes the problem that #
					# the scores depends on how the "main affected allele" is #
					# chosen.                                                 #
					###########################################################
				    &choose_best_segregation_scores_indels_merged;

				    print INDELS_MERGED "\t$segregation_score_best_1";
					print INDELS_MERGED "\t$segregation_score_best_2";
					print INDELS_MERGED "\t$segregation_score_best_3";
					print INDELS_MERGED "\t$segregation_score_best_4";
					print INDELS_MERGED "\t$segregation_score_best_5";
					print INDELS_MERGED "\t$segregation_score_best_6";
					print INDELS_MERGED "\t$segregation_score_best_7";
					
					print INDELS_MERGED "\t$segregation_score_best_overall"; # Best overall


					##############################
					# Work out the merged scores #
					##############################
					if ($effect_score >= $effect_score_check){$effect_score_merged = $effect_score} else{$effect_score_merged = $effect_score_check}
					$consequence_merged = $consequence." ".$consequence_check;


					#print "Effect score:        $effect_score\n";
					#print "Effect score check:  $effect_score_check\n";
					#print "Effect score merged: $effect_score_merged\n\n";

					#print "Consequence:         $consequence\n";
					#print "Consequence check:   $consequence_check\n";
					#print "Consequence merged:  $consequence_merged\n\n";


					#######################################################################
				    # Now write the columns to the right of the merged segregation scores #
				    #######################################################################
				    print INDELS_MERGED "\tNA"; # Main Affected Allele is REF - not applicable
				    print INDELS_MERGED "\t$effect_score_merged";
				    print INDELS_MERGED "\t$consequence_merged";
				    print INDELS_MERGED "\tNA"; # Alt
				    print INDELS_MERGED "\tNA"; # Allele_A
				    print INDELS_MERGED "\tNA"; # Allele_B
				    print INDELS_MERGED "\t$main_affected_allele"; # Main Affected Allele
				    print INDELS_MERGED "\tNA";
				    print INDELS_MERGED "\tNA";
				    print INDELS_MERGED "\tNA";

					print INDELS_MERGED "\n"; # end of line

					last;  # Don't look forward any more as we only want to merge adjacent indels

					#print "\n\n++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n";
					#if ($position eq "99999999"){$answer=<STDIN>;}

				} # If next indel is within merge threshold

			} # if not merged already "M" # not currently used at all

	    } # check_count loop to check forward for indels to merge  

	} # line_count loop

	close INDELS_MERGED;
} # merge_indels


###########################################################
# Open the input VCF file and check the file (including   #
# headers) then close once checks have been made          #
###########################################################
sub check_VCF_file
{
	open (VCF_TEMPORARY, "$vcf_file") || die "Cannot open $vcf_file";
	$line_count = 0;
	$vep_line_count = 1;
	$passed_header_lines = "false";
	$vcf_file_closed = "false";

	while ($vcf_file_closed eq "false")
	{
		$vcf_line = <VCF_TEMPORARY>;

		if ($vcf_file_closed eq "false")
		{
			chomp $vcf_line;
			&chomp_all ($vcf_line);

			$line_count = $line_count + 1;

			############################################################################
			# Check that first line looks like a VCF file (subroutine needs improving) #
			############################################################################
			&check_vcf_format;

			#################################################################
			# Check the VCF header lines (before starting on the real data) #
			# Check for which variant caller was used.                      #
			# Also get information on the order of fields in the INFO field #
			# used by VEP or snpEff                                         #
			#################################################################
			&check_vcf_header_lines;

		} # if ($vcf_file_closed eq "false")

	} # while loop

} # check_VCF_file

#######################################################################################
# Reads the disease status file to get which samples are cases and which are controls #
#######################################################################################
sub read_disease_status_file
{
	open (STATUS_FILE, "$disease_status_file") || die "Cannot open $vcf_file";
	$disease_status_count = 0;
	$affected_count_ds_file = 0;
	$normal_count_ds_file = 0;

	while ($disease_status_line = <STATUS_FILE>) 
	{
		chomp $disease_status_line;

		if ($disease_status_line ne "")
		{
			$disease_status_count = $disease_status_count + 1;
			
			#print "$disease_status_count \t$disease_status_line\n";

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
					($disease_status ne "normal") &&
					($disease_status ne "omit"))
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
					$affected_count_ds_file = $affected_count_ds_file + 1;
				}
				if (($disease_status eq "clear") || ($disease_status eq "normal") || ($disease_status eq "control"))
				{
					$disease_status = "normal";
					$normal_count_ds_file = $normal_count_ds_file + 1;
				}
				
				$disease_status_array[$disease_status_count] = $disease_status;
				$disease_status_sample_array[$disease_status_count] = $sample_name_from_status_file;
			
			} # if array size = 2

		} # if ($disease_status_line ne "")

	}	# reading STATUS_FILE
	
	$no_samples_in_status_file = $disease_status_count;
	
	
	############################################
	# Show user a list of the disease statuses #
	############################################
	&print_message ("List of samples and disease statuses","message");
	print "\n";

	for ($disease_status_count = 1;$disease_status_count <=$no_samples_in_status_file; $disease_status_count++)
	{
		# pad out sample name to max_sample_name_length characters + 2
		$sample_name = sprintf("%-*s", $max_sample_name_length + 2, $disease_status_sample_array[$disease_status_count]);

		print "\t$disease_status_count\t$sample_name\t$disease_status_array[$disease_status_count]\n";
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
			if ($disease_status_sample_array[$disease_status_count] eq $sample_name_array_input_order[$sample_count])
			{
				$name_match_count = $name_match_count + 1; # Checks each name only matches once
			}
		}
		# Check each name only matches once
		if ($name_match_count == 1){$total_match_count = $total_match_count + 1}
		if ($name_match_count > 1){&print_message("Too many names match!  Check file","warning");exit}

	} # sample_count loop
	
	
	################################################################
	# Warn if not all the sample names match                       #
	################################################################
	if ($total_match_count != $no_of_samples)
	{
		&print_message("ERROR!  All names in Disease status file do not match those in VCF file.","warning");
		
		print "Sample columns found in the VCF file (from the '#CHROM' line)\n\n";
		
		for ($sample_count = 1; $sample_count <= $no_of_samples; $sample_count++)
		{
			print "\t$sample_count \t$sample_name_array_input_order[$sample_count]\n";
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
		print "\n\t>>If these samples and disease statuses look correct, press 'RETURN' ";
		$answer=<STDIN>;
		
		&print_message("All names in Disease status file match those in VCF file.","message");
		
		print "\t($total_match_count names out of $no_of_samples match)\n\n";
	}

	close STATUS_FILE;

}# read_disease_status_file

############################################################
# Get run title to add to name of output file              #
############################################################
sub get_name_of_output_files
{
	$prefix = &get_prefix ($vcf_file);
	$prefix_disease_status_file = &get_prefix($disease_status_file);
	$answer = "n";

	until ($answer eq "y")
	{
		&print_message("Type a prefix that will identify the output files","message");

		print "The input VCF file is:      \t$vcf_file\n";
		print "The disease status file is: \t$disease_status_file\n\n";

		print "Choose file prefix:    [DEFAULT = $prefix"."_".$prefix_disease_status_file."]\n\n";

		print "> ";

		$run_title = <STDIN>;
		chomp $run_title;

		if ($run_title eq ""){$run_title = $prefix."_".$prefix_disease_status_file}


		################################################################
		# Create names for output files with run_title                 #
		################################################################
		$output_file = 						$run_title."_for_excel.txt";
		$output_file_filtered = 			$run_title."_filtered_for_excel.txt";
		$output_file_filtered_and_sorted = 	$run_title."_filtered_and_sorted_for_excel.txt";

		$simplified_indels_output_file = 	$run_title."_indels_only.txt";
		$indels_merged_file = 				$run_title."_indels_merged.txt";
		$temp_indels_file = 				"temp_indels_only_file.txt"; # temporary

		$command_log =      				$run_title."_vcf2excel_command_log.out";
		$default_vep_file = 				$prefix."_vep.out";			

		$missing_effects = 					$run_title."_vcf2excel_missing_effects.out";

		print "\nOutput file for Excel: \t$output_file\n\n";

		print "Is this output prefix OK? (y/n)         ";
		$answer=<STDIN>;
		chomp $answer;

		# Default is y
		if ($answer eq "") {$answer = "y"}
		$answer= lc($answer);


		#####################################################
		# Warn if output file with this name already exists #
		#####################################################
		if ((-e $output_file) && ($output_file ne "test_for_excel.txt"))
		{
			&print_message("File $output_file already exists!","warning");

			print "Do you want to overwrite this file? (y/n)   ";
			$answer=<STDIN>;
			$answer = lc $answer;
			if ($answer eq "") {$answer="n";} # default
		}

	} # until $answer eq "y"

}# get_name_of_output_files


########################################################################################
# Choose the alleles for all the samples and save them in allele_for_writing_array     #
# For the SNPs store the actual base (ACGT) and for the indels use a code a,b,c,d etc  #
#
# For indels, when the genotype is 'X' you can't use the indel_hash, so...
#
# Note: the segregation score comparison are done on the genotypes in the form 'AA',   #
# 'AB' etc and the 'allele_for_writing' is the actual base (A, C, G or T) for SNPs     #
# and 'a', 'b', 'c' etc for indels.                                                    #
########################################################################################
sub choose_alleles_for_all_samples
{
	for ($sample_count = 1; $sample_count <=$no_of_samples_without_omits; $sample_count++)
	{

		$col_1 = ($sample_count * 2) - 1;
		$col_2 = $sample_count * 2;

		if ($variant_type eq "snp")
		{	
			$allele_base_1 = $base_array[$col_1];
			$allele_base_2 = $base_array[$col_2];
		}
		elsif ($variant_type eq "indel")
		{
			if ($base_array[$col_1] ne "X")
			{
				$allele_base_1 = $indel_hash{$base_array[$col_1]};
				$allele_base_2 = $indel_hash{$base_array[$col_2]};
			}
			else
			{
				$allele_base_1 = "x";
				$allele_base_2 = "x";
			}
		}

		###########################################################################
		# OK. When the indel_hash is created it stores the ALT alleles listed     #
		# in the ALT_alleles column.  If there is no genotype the base_array is X #
		# but this isn't one of the ALT alleles so doesn't go into the indel_hash #
		###########################################################################
		if ($position eq $position_to_debug)
		{
			print "===========================================================================000===================\n";
			print "Sample count: $sample_count\n\n";

			#print "indel_hash{$base_array[$col_1]}: $indel_hash{$base_array[$col_1]}\n";
			#print "indel_hash{$base_array[$col_2]}: $indel_hash{$base_array[$col_2]}\n\n";

			print "base_array[$col_1]: $base_array[$col_1]\n";
			print "base_array[$col_2]: $base_array[$col_2]\n\n";

			print "base_orig_array[$col_1]: $base_orig_array[$col_1]\n";
			print "base_orig_array[$col_2]: $base_orig_array[$col_2]\n\n";

			print "allele_base_1: $allele_base_1\n";
			print "allele_base_2: $allele_base_2\n";

			print "===========================================================================000===================\n";
			$answer=<STDIN>;
		}

		#########################################################
		# Put alleles in alphabetical order for Excel           #
		# (just because this is the way NGS SNP Handler did it) #
		#########################################################
		if (ord $allele_base_1 > ord $allele_base_2)
		{
			$myString = $allele_base_1;
			$allele_base_1 = $allele_base_2;
			$allele_base_2 = $myString;
		}


		#########################################################
		# Save alleles for writing to file later                #
		#########################################################
		$allele_1_for_writing_array[$sample_count] = $allele_base_1;
		$allele_2_for_writing_array[$sample_count] = $allele_base_2;


		#########################################################
		# If the original base was X save this instead          #
		#########################################################
		if ($base_orig_array[$col_1] eq "X")
		{
			if ($variant_type eq "snp")
			{
				$allele_1_for_writing_array[$sample_count] = $base_orig_array[$col_1];
				$allele_2_for_writing_array[$sample_count] = $base_orig_array[$col_1];
			}
			else
			{
				$allele_1_for_writing_array[$sample_count] = "x";
				$allele_2_for_writing_array[$sample_count] = "x";
			}

			# CHECKING
			if (($base_array[$col_1]) ne ($base_orig_array[$col_1]))
			{
				print "\n\n>>>>>>>>>>>>> DIFFERENCE FOUND! <<<<<<<<<<<<<<<\n\n";
				print "base_array[$col_1]:      $base_array[$col_1]\n";
				print "base_orig_array[$col_1]: $base_orig_array[$col_1]\n\n";

				print "Please tell Mike this useful information\n";

				$answer=<STDIN>;
			}
		}

	} # sample_count loop for choosing the alleles

}# choose_alleles_for_all_samples