#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	pool_vcf_files  						                            #     
#									                                    #
#	THIS PERL SCRIPT POOLS VCF FILES                                    #
#									                                    #
#########################################################################

#############################
# Mike Boursnell Nov 2013   #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;
use Cwd;


my $version							= "11";
my $debug							= "false";
my $show							= "false";

my $answer							= "";
my $single_line						= "";
my $single_line_new					= ""; # For pooled output
my $start_storing					= "false";
my $chromosome						= "";
my $position						= "";
my $position_new					= "";
my $ID								= "";
my $INFO							= "";
my $FORMAT							= "";
my $chr_and_pos						= ""; #chromosome and position linked e.g. 12_6456378
my $chr_and_pos_original			= "";
my $chr_and_pos_duplicates			= ""; #
my $REF_base						= "";
my $ALT_base_overall				= "";
my $ALT_base_overall_2				= "";
my $QUAL							= "";
my $FILTER							= "";
my $run_title						= "";
my $genotypes_string				= "";
my $depth_2_string					= "";
my $depth							= "";
my $myString						= "";
my $GT_string						= "";
my $DP_string 						= "";
my $variant_type					= ""; # "snp" or "indel"
my $variant_caller					= "UNKNOWN";
my $default							= "";
my $padding							= "";
my $string_1						= "";
my $string_2						= "";
#files
my $vcf_file						= "";
my $pooled_vcf_file					= "";
my $command_log						= "";
my $discordant_variants_file		= "";
my $concordant_variants_file		= "";
my $only_in_1_variants_file			= "";
my $only_in_2_variants_file			= "";
my $in_both_variants_file			= "";
my $discordant_variants_details		= "";
my $venn_file						= "";
my $allele_1						= "";
my $allele_2						= "";
my $locus_genotype_string			= ""; # a string og CHR POS GENOTYPE_STRING (chr_pos-genotype_string)

#Boolean
my $all_genotypes_match				= ""; #"true" or "false"
my $sample_name_mismatch			= ""; #"true" or "false"  - checks if sample names are the same for the two VCF files
my $matching_genotype				= ""; #"true" or "false"
my $matching_position				= "";

#Counts of various types
my $defined_count_1					= 0; # no of counts in position_hash
my $position_hash_size				= 0; # no of keys in position_hash
my $compare_hash_size				= 0; # no of keys in compare_hash
my $concordant_position_count		= 0; # Each position will have a number of genotpyes (same as no of samples)
my $discordant_position_count		= 0;
my $discordant_snps_count			= 0; # This refers to position not individual genotype
my $discordant_indels_count			= 0; # This refers to position not individual genotype
my $only_in_file_1_count			= 0;
my $only_in_file_2_count			= 0;
my $variants_in_file_count			= 0;
my $variants_in_file_1_count		= 0;
my $discordant_one_null_vcount		= 0; # vcount means count of individual sample genotype, not just position. Null genotype is ./.
my $discordant_neither_null_vcount	= 0; # vcount means count of individual sample genotype, not just position.
my $concordant_null_vcount			= 0; # vcount means count of individual sample genotype, not just position. Null genotype is ./.
my $concordant_neither_null_vcount	= 0; # vcount means count of individual sample genotype, not just position.
my $snps_in_file_1_count			= 0;
my $snps_in_file_2_count			= 0;
my $indels_in_file_1_count			= 0;
my $indels_in_file_2_count			= 0;
my $unknowns_in_file_1_count		= 0;
my $unknowns_in_file_2_count		= 0;
my $total_variants_in_file			= 0;
my $total_variants_in_file_check	= 0;
my $total_variants_in_file_1		= 0;
my $file_1_and_2_count				= 0;
my $file_1_only_calc				= 0;
my $file_count						= 0;
my $check_genotypes_count			= 0; # All genotypes counted
my $pos								= 0;
my $pooled_snp_count				= 0;
my $total_variants_in_pooled		= 0;
my $check_count						= 0;
my $duplicates_count				= 0;
my $row_count						= 0;

#Other
my $no_of_vcf_files					= 0;
my $line_count						= 0;
my $chrom_line_array_size			= 0;
my $array_count						= 0;
my $no_of_samples_chrom_line		= 0;
my $vcf_line_array_size					= 0;
my $vcf_format_field_array_size					= 0;
my $vcf_data_field_array_size					= 0;
my $no_of_data_columns				= 0;
my $sample_count					= 0;
my $no_of_alt_alleles_1				= 0;
my $no_of_alt_alleles_2				= 0;
my $array_size_genotypes_1			= 0;
my $array_size_genotypes_2			= 0;
my $variant_at_same_pos_count		= 0; # if there are more than one variant for a single position
my $max_variants_at_same_pos		= 0; # store this to check if there are any with 3 at the same position
my $allele_count					= 0;
my $percent_discordant				= 0;

my @vcf_file_array					= ();
my @chrom_line_array				= ();
my @sample_name_array				= ();
my @sample_name_array_2				= ();
my @sample_name_string				= (); # holds a string of the samples names for each input VCF file
my @vcf_line_array						= ();
my @vcf_format_field_array						= ();
my @vcf_data_field_array						= ();
my @ALT_allele_array				= ();
my @ALT_allele_array_2				= ();
my @genotype_block_array			= ();
my @genotype_all_1_array			= ();
my @genotype_all_2_array			= ();
my @genotypes_1_array				= (); # split genotype string into individual genotypes whilst comparing
my @genotypes_2_array				= (); # split genotype string into individual genotypes whilst comparing
my @depth_2_array					= (); # split depth string into individual genotypes whilst comparing
my @variant_caller_array			= (); # The variant caller used by each VCF file
my @lines_added_array				= (); # records the number of lines added with each new file
my @variants_in_file_array			= ();
my @found_in_pooled					= (); # if variant in single VCF file is found in the pool of all variants
my @already_defined_array			= ();
my @not_defined_array				= ();
my @defined_but_different_array		= ();
my @total_rows_array				= (); # stores number of rows in each file

# 2 dimensional array for storing file number and locus_genotype string
my @file_number 					= ();
my @row_number 						= ();
my @locus_genotype_array_2d 		= (\@file_number, \@row_number);

my %position_hash					= (); # hash table for positions in all VCF files
my %compare_hash					= (); # hash table for comparing files
my %compare_hash_SNP				= (); # hash table for comparing files - SNPs only
my %compare_hash_IND				= (); # hash table for comparing files - Indels only
my %genotypes_hash					= (); # hash table for genotypes in VCF file 1
my %pooled_hash						= (); # hash table for all pooled locus_genotype strings
my %venn_hash						= (); # stores data in form XXXMXX for Venn diagram
my %venn_hash_SNP					= (); # stores data in form XXXMXX for Venn diagram - SNPs only
my %venn_hash_IND					= (); # stores data in form XXXMXX for Venn diagram - Indels only

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "              pool_vcf_files   \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "    - This program takes a series of VCF files produced by different variant callers\n";
print "      (for example GATK, samtools, platypus)\n";
print "      and combines the list of variants produced into a single file, with no duplicates.\n";

print color 'reset';

&print_message("Here is a list of VCF files in this directory","message");
print "\n";
system("ls -1 *.vcf");
print "\n";


&print_message("How many VCF files do you want to combine","input");

$answer=<STDIN>;
chomp $answer;
$no_of_vcf_files = $answer;


#################
# Get VCF files #
#################
for ($file_count = 1; $file_count <= $no_of_vcf_files; $file_count++)
{
	until (-e "$vcf_file")
	{
		print "Name of file $file_count:    ";

		$vcf_file = <STDIN>;
		chomp $vcf_file;

		if (($file_count ==1) && ($vcf_file eq "")){$vcf_file = "AS_HC_HC.vcf"}
		if (($file_count ==2) && ($vcf_file eq "")){$vcf_file = "AS_HC_UG.vcf"}
		if (($file_count ==3) && ($vcf_file eq "")){$vcf_file = "AS_HC_mpileup.vcf"}
		if (($file_count ==4) && ($vcf_file eq "")){$vcf_file = "AS_HC_platypus.vcf"}
		if (($file_count ==5) && ($vcf_file eq "")){$vcf_file = "AS_HC_gVCF.vcf"}
		if (($file_count ==6) && ($vcf_file eq "")){$vcf_file = "AS_HC_freebayes.vcf"}

		if (! -e $vcf_file){print "\n\n>>>>>>    $vcf_file not found!.  Try again   <<<<<<\n\n";}
	}
	$vcf_file_array[$file_count] = $vcf_file;
	$vcf_file = "";
}

###################################
# Choose name for pooled VCF file #
###################################

&print_message("Choose a name for the pooled VCF file","input");

print "\n";

$pooled_vcf_file = <STDIN>;
chomp $pooled_vcf_file;

if ($pooled_vcf_file eq ""){$pooled_vcf_file = "Pooled_VCF_".$no_of_vcf_files."_files.vcf"}


########################
# Open pooled VCF file #
########################

open (POOLED, ">$pooled_vcf_file")|| die "Cannot create output file: $pooled_vcf_file";


print "\nList of VCF files chosen:\n\n";


#################################
# Show user a list of the files #
#################################
for ($file_count = 1; $file_count <= $no_of_vcf_files; $file_count++)
{
	print "$file_count:\t$vcf_file_array[$file_count]\n";
}
$run_title = "run_title";

$command_log = "pool_vcf_files_"."$run_title"."_command_log.out";


##################################################
# Ask which variant caller is used for each file #
##################################################
&print_message("Which variant callers were used for these files?","input");

for ($file_count = 1; $file_count <= $no_of_vcf_files; $file_count++)
{
	print "\n\n$vcf_file_array[$file_count].  Which variant caller was used for this file?\n\n";

	# Try to guess variant caller #
	$default="";
	if (index($vcf_file_array[$file_count],"HC") > -1){$default = "HaplotypeCaller"}
	if (index($vcf_file_array[$file_count],"UG") > -1){$default = "UnifiedGenotyper"}
	if (index($vcf_file_array[$file_count],"gVCF") > -1){$default = "HaplotypeCallergVCF"}
	if (index($vcf_file_array[$file_count],"plat") > -1){$default = "platypus"}
	if (index($vcf_file_array[$file_count],"mpile") > -1){$default = "mpileup"}
	if (index($vcf_file_array[$file_count],"ayes") > -1){$default = "freeBayes"}

	print "  <1> Unified Genotyper\n";
	print "  <2> Haplotype Caller\n";
	print "  <3> Haplotype Caller gVCF mode\n";
	print "  <4> Platypus\n";
	print "  <5> samtools mpileup\n";
	print "  <6> freeBayes\n\n";

	if ($default eq ""){print "> "}
	if ($default ne ""){print "> (default = $default)   "}
	$answer=<STDIN>;
	chomp $answer;

	if ($answer eq ""){$variant_caller = $default}

	if (substr($answer,0,1) eq "1"){$variant_caller = "UnifiedGenotyper"}
	if (substr($answer,0,1) eq "2"){$variant_caller = "HaplotypeCaller"}
	if (substr($answer,0,1) eq "3"){$variant_caller = "HaplotypeCallergVCF"}
	if (substr($answer,0,1) eq "4"){$variant_caller = "platypus"}
	if (substr($answer,0,1) eq "5"){$variant_caller = "mpileup"}
	if (substr($answer,0,1) eq "6"){$variant_caller = "freeBayes"}

	$variant_caller_array[$file_count] = $variant_caller;
}

&print_message("Check on variant callers for the files","message");

for ($file_count = 1; $file_count <= $no_of_vcf_files; $file_count++)
{
		print"File $file_count\t$vcf_file_array[$file_count]\t$variant_caller_array[$file_count]\n";
}

print "\n  >> Press 'return' to continue      ";
$answer=<STDIN>;


#########################
# Open command log file #
#########################
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "pool_vcf_files  version $version\n\n";	
		


####################################################
# Print the headers for the POOLED output VCF file #
####################################################
print POOLED "##fileformat=VCFv4.0\n";
print POOLED "##source=pool_vcf_files_v$version\n";
for ($file_count = 1; $file_count <= $no_of_vcf_files; $file_count++)
{
	print POOLED "##input_vcf_file_$file_count $vcf_file_array[$file_count] $variant_caller_array[$file_count]\n";
}
print POOLED "##INFO=<ID=CALLER,Number=1,Type=String,Description=\"Which variant_caller was used on this line\">\n";
print POOLED "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
print POOLED "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">\n";
print POOLED "##INFO=<ID=TYPE,Number=A,Type=String,Description=\"The type of allele, either snp, mnp, ins, del, or complex.\">\n";



##################################
# Open the VCF files MAIN LOOP   #
##################################
$variants_in_file_1_count = 0;

&print_message("Reading and pooling VCF files...","message");

for ($file_count = 1; $file_count <= $no_of_vcf_files; $file_count++)
{
	$vcf_file = $vcf_file_array[$file_count];
	$start_storing = "false";
	$lines_added_array[$file_count] = 0;
	$already_defined_array[$file_count] = 0;
	$not_defined_array[$file_count] = 0;
	$defined_but_different_array[$file_count] = 0;
	$total_rows_array[$file_count] = 0;
	$row_count = 0;
	$variant_caller = $variant_caller_array[$file_count]; # set variant caller as chosen by user for each VCF file

	&print_both("Reading VCF file $file_count/$no_of_vcf_files  \t$vcf_file\n");

	open (VCF, "$vcf_file") || die "Cannot open $vcf_file";
	$line_count = 0;

	$line_count = 0;
	$variants_in_file_count = 0; # set to zero for each VCF file
	

	while ($single_line = <VCF> ) 
	{
		chomp $single_line;
		&chomp_all ($single_line);

		$line_count = $line_count + 1;

		if ($line_count % 20000 == 0) {print "  File $file_count/$no_of_vcf_files    \tLine $line_count\n"}

		######################################################
		# Only do this next bit if 'start_storing is True    #   
		# i.e. if the line containing #CHROM has been found  #
		######################################################
		if ($start_storing eq "true")
		{
			
			###############################################################
			# Read VCF line into the array vcf_line_array - split at tabs #
			###############################################################
			&read_vcf_line;


	        ############################################################################
	        # Initial reading VCF files..                                              #
	        # Now we have read all the elements we can write to the output pooled file #
	        # and store details in the hash arrays                                     #
	        ############################################################################


			###########################################################################
			# First VCF file....                                                      #
			# If this is the first file, write the line to the pooled VCF file        #
			#                                                                         #
			# Also store all the positions in the position_hash
			###########################################################################

			if ($file_count == 1)
			{   
				&read_line_if_first_vcf_file;   # i.e. all lines are retained.
			}

			
			###########################################################################
			# Subsequent VCF files...                                                 #    
			# If it is a subsequent file, check to see if there is a position there   #
			# and if not, add the new position into the hash, and also write the line #
			# at the bottom of the new VCF file.                                      #
			###########################################################################

			if ($file_count > 1)
			{   
				&read_line_if_subsequent_vcf_file;   # i.e. only unique lines are retained.
			}

		} # if start_storing is true


		#########################################################
	    # VCF:  Has the start of the CHR data been reached yet? #
		# The line before the actual data starts with '#CHROM'  #
		# This part of the code checks that line, which has the #
		# (very important) headers for the sample columns.      #
	    #########################################################
		
		if (index($single_line,"#CHROM") > -1)
		{
			if ($file_count==1)
			{
				print POOLED "$single_line\n";
			}
			$start_storing = "true";
			&get_samples_from_chrom_line
		}

	} # reading VCF file

	close VCF;
	print "\n";

	$variants_in_file_array[$file_count] = $variants_in_file_count;

	$total_rows_array[$file_count] = $row_count;

} # End of MAIN file_count loop

$total_variants_in_pooled = $pooled_snp_count;

close POOLED;
&print_both("Finished reading $no_of_vcf_files VCF files\n\n");


################################################################
# Compare all VCF files.  See which variants are common to all #
# and get data for a complete Venn diagram                     #
################################################################

&print_message("Comparing all $no_of_vcf_files VCF files to the pooled list of variants...","message");
print COMMAND_LOG "Comparing all $no_of_vcf_files VCF files...\n\n";


for ($file_count = 1; $file_count <= $no_of_vcf_files; $file_count++)
{
	$vcf_file = $vcf_file_array[$file_count];
	print "Comparing VCF file $file_count/$no_of_vcf_files      $vcf_file\n";

	#############################################
	# Loop the number of rows in each VCF file  #
	#############################################
	for ($row_count = 1; $row_count <= $total_rows_array[$file_count]; $row_count++)
	{
		$locus_genotype_string = $locus_genotype_array_2d[$file_count][$row_count];

		if (defined $pooled_hash{$locus_genotype_string})
		{
			if (index($locus_genotype_string,"_IND") > -1){$variant_type = "indel"}
			elsif (index($locus_genotype_string,"_SNP") > -1){$variant_type = "snp"}
			else {$variant_type = ""}


			#############################################################
			# If this is the first time that compare_hash has been used #
			# Create the Venn string (MMXXXM) from scratch              #
			#############################################################
			#All
			if (not defined $compare_hash{$locus_genotype_string})
			{
				$string_1 = substr("XXXXXX",0,$file_count -1);
				$string_2 = substr("XXXXXX",$file_count,99);
				$compare_hash{$locus_genotype_string} = "$string_1"."M"."$string_2";
			}

			#SNP
			if ($variant_type eq "snp")
			{
				if (not defined $compare_hash_SNP{$locus_genotype_string})
				{
					$string_1 = substr("XXXXXX",0,$file_count -1);
					$string_2 = substr("XXXXXX",$file_count,99);
					$compare_hash_SNP{$locus_genotype_string} = "$string_1"."M"."$string_2";
				}
			} # compare_hash_SNP

			#Indel
			if ($variant_type eq "indel")
			{
				if (not defined $compare_hash_IND{$locus_genotype_string})
				{
					$string_1 = substr("XXXXXX",0,$file_count -1);
					$string_2 = substr("XXXXXX",$file_count,99);
					$compare_hash_IND{$locus_genotype_string} = "$string_1"."M"."$string_2";
				}
			} # compare_hash_IND

			
			
			###################################################################
			# If compare_hash has already been defined by a previous VCF file #
			# we need to add data to the existing string                      #
			###################################################################
			#All
			if (defined $compare_hash{$locus_genotype_string})
			{
				$string_1 = substr($compare_hash{$locus_genotype_string},0,$file_count -1);
				$string_2 = substr($compare_hash{$locus_genotype_string},$file_count,99);
				$compare_hash{$locus_genotype_string} = "$string_1"."M"."$string_2";
			}

			#SNP
			if ($variant_type eq "snp")
			{
				if (defined $compare_hash_SNP{$locus_genotype_string})
				{
					$string_1 = substr($compare_hash_SNP{$locus_genotype_string},0,$file_count -1);
					$string_2 = substr($compare_hash_SNP{$locus_genotype_string},$file_count,99);
					$compare_hash_SNP{$locus_genotype_string} = "$string_1"."M"."$string_2";
				}
			} # compare_hash_SNP

			#Indel
			if ($variant_type eq "indel")
			{
				if (defined $compare_hash_IND{$locus_genotype_string})
				{
					$string_1 = substr($compare_hash_IND{$locus_genotype_string},0,$file_count -1);
					$string_2 = substr($compare_hash_IND{$locus_genotype_string},$file_count,99);
					$compare_hash_IND{$locus_genotype_string} = "$string_1"."M"."$string_2";
				}
			} # compare_hash_IND
		}

	}

} # second main loop to compare all against pooled. file_count


&print_message("Getting data for Venn diagrams...","message");
print COMMAND_LOG "Getting data for Venn diagrams...\n";


#########################################################################
# Now go through all the positions again making data for a Venn diagram #
#########################################################################

##########################
# Open Venn diagram file #
##########################
$venn_file = "Venn_".$no_of_vcf_files.".txt";

open (VENN, ">$venn_file")|| die "Cannot create output file: $venn_file";

print VENN "INDEX\tCOUNT\tSTRING\n";

$pooled_snp_count =0;
#All
while ( my ($key, $value) = each(%compare_hash) ) 
{
    $pooled_snp_count = $pooled_snp_count + 1;
    print VENN "$pooled_snp_count\t$key\t$value\n";

    if (not defined $venn_hash{$value}) {$venn_hash{$value} = 1}

    if (defined $venn_hash{$value}) {$venn_hash{$value} = $venn_hash{$value} + 1}
}

#SNPs
while ( my ($key, $value) = each(%compare_hash_SNP) ) 
{
    if (not defined $venn_hash_SNP{$value}) {$venn_hash_SNP{$value} = 1}
    if (defined $venn_hash_SNP{$value}) {$venn_hash_SNP{$value} = $venn_hash_SNP{$value} + 1}
}

#Indels
while ( my ($key, $value) = each(%compare_hash_IND) ) 
{
    if (not defined $venn_hash_IND{$value}) {$venn_hash_IND{$value} = 1}
    if (defined $venn_hash_IND{$value}) {$venn_hash_IND{$value} = $venn_hash_IND{$value} + 1}
}

#All
print VENN "STRING\tCOUNT\tVARIANT\n";
while ( my ($key, $value) = each(%venn_hash) ) 
{
		print VENN "$key\t$value\tAll\n";
}

#SNPs
print VENN "STRING\tCOUNT\tVARIANT\n";
while ( my ($key, $value) = each(%venn_hash_SNP) ) 
{
		print VENN "$key\t$value\tSNP\n";
}

#Indels
print VENN "STRING\tCOUNT\tVARIANT\n";
while ( my ($key, $value) = each(%venn_hash_IND) ) 
{
		print VENN "$key\t$value\tIND\n";
}

print "\nPooled SNP count: $pooled_snp_count\n\n";
close VENN;


##########################################################################
# Check on sample names to make sure they are the same for each VCF file #
##########################################################################
print COMMAND_LOG "\nCheck on sample names in the input VCF files:\n\n";
$sample_name_mismatch = "false";
for ($file_count = 1; $file_count <= $no_of_vcf_files; $file_count++)
{
	print COMMAND_LOG "VCF file $file_count:\t$sample_name_string[$file_count]\n";

	# if any sample name string do not match no. 1 then warn (TEMP!! all need checking against all)
	if ($sample_name_string[1] ne $sample_name_string[$file_count]){$sample_name_mismatch = "true"}
}


if ($sample_name_mismatch eq "true")
{
	&print_message("SAMPLE NAME MISMATCH","warning");

	&print_both("Sample names do not match.  You need to check the #CHROM lines of the VCF files\n\n");

	$answer=<STDIN>;

}


&print_message("FINISHED POOLING VCF FILES","message");

#print "Maximum no. of variants at a single position: $max_variants_at_same_pos\n\n";


	&print_both ("Number of new variants added by each VCF file:\n\n");
    &print_both ("File 1\t$vcf_file_array[1]\n");

    &print_both ("\tTotal variants in file:   \t$variants_in_file_array[1]\n\n");

for ($file_count = 2; $file_count <= $no_of_vcf_files; $file_count++)
{
	#$padding = &add_space($vcf_file_array[$file_count]);
	&print_both( "File $file_count\t$vcf_file_array[$file_count]\t$variant_caller\n");

	&print_both ("\tTotal variants in file:   \t$variants_in_file_array[$file_count]\n");
	&print_both ("\tVariants added:           \t$lines_added_array[$file_count]\n");
	&print_both ("\tAlready defined:          \t$already_defined_array[$file_count]\n");
	&print_both ("\tDefined but different:   \t$defined_but_different_array[$file_count]\n");
	&print_both ("\tNot defined:              \t$not_defined_array[$file_count]\n\n");
}

&print_both("Output pooled VCF: $pooled_vcf_file\n\n");

exit;
################################################################################################################################################


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
		
	print("$unix_command\n\n");
	system("$unix_command");
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	print COMMAND_LOG "$unix_command\n";

}

##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

foreach (@_) {s/\n//g}  
foreach (@_) {s/\r//g}  

}

####################################################################
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test" from "test.fasta")                           #
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

####################################################
# Subroutine to print to screen and to COMMAND_LOG #
####################################################

sub print_both
{
	my $message = $_[0];

	print "$message";
	print COMMAND_LOG "$message";
}

sub add_space
{
	my $string = $_[0];
	my $space = "";
	my $count = 0;

	for($count=1; $count < 40-length($string) ; $count++)
	{
		$space = $space." ";
	}
}


sub read_vcf_line
{
	#########################################################
    # Split whole line at TABs into an array vcf_line_array (1-9) #
    #########################################################
	@vcf_line_array = split (/\t/,$single_line);
	$vcf_line_array_size = scalar @vcf_line_array;
	$no_of_data_columns = $vcf_line_array_size - 9;


	###################################################
    # VCF: Get main columns of the VCF file           #
    # CHR, POS, REF, ALT, QUAL, FILTER, INFO, FORMAT  #
    ###################################################
	$chromosome = $vcf_line_array[0];
    $position = $vcf_line_array[1];
    $ID = $vcf_line_array[2];
    $REF_base = $vcf_line_array[3];
    $ALT_base_overall = $vcf_line_array[4];
    $QUAL = $vcf_line_array[5];
    $FILTER = $vcf_line_array[6];
    $INFO = $vcf_line_array[7];
    $FORMAT = $vcf_line_array[8];

    $variants_in_file_count = $variants_in_file_count + 1;

    if ($position eq "39013405xx"){$show = "true"}else{$show = "false"} ## current check TEMP!!!!!


    ############################################
    # VCF:  Various conversions on chromosome  #
    ############################################
    # Convert chrX to chr39 for dog 
	if ($chromosome eq "chrX"){$chromosome = "chr39"}
     
	# Remove string 'chr'
    if (substr($chromosome,0,3) eq "chr"){$chromosome = substr($chromosome,3,99)}
	
	# If chromosome is an integer and less than 10 then add a leading zero
	if ($chromosome =~ /^\d+?$/){if ($chromosome < 10){$chromosome = "0".$chromosome}}


	#######################################################################################
	# VCF: Split ALT_base_overall at commas to get all the different alleles in the array #
	#######################################################################################
	@ALT_allele_array = split(",",$ALT_base_overall);
	$no_of_alt_alleles_1 = (scalar @ALT_allele_array);


	###########################################################
	# VCF:  Decide whether it is a SNP or an Indel VCF1       #
	###########################################################
	$variant_type = "snp";
	for ($allele_count = 0; $allele_count <$no_of_alt_alleles_1; $allele_count++)
	{
		if (length($REF_base) != length($ALT_allele_array[$allele_count])){$variant_type = "indel"}
	}

	if ($show eq "true"){print "At $position, variant_type is $variant_type\n";}


	#########################
	# Count SNPs and Indels #
	#########################
	if ($variant_type eq "snp"){$snps_in_file_1_count = $snps_in_file_1_count + 1}
	if ($variant_type eq "indel"){$indels_in_file_1_count = $indels_in_file_1_count + 1}
	if (($variant_type ne "snp") && ($variant_type ne "indel")){$unknowns_in_file_1_count = $unknowns_in_file_1_count + 1}


	#################################################################
	# Any VCF file...                                               #
	# Make up a position which combines CHR and POS                 #
	# (in case two positions are the same on different chromosomes) #
	#################################################################
	$chr_and_pos = $chromosome."_".$position;
	
	###########################################################
	# If it is an indel add "_IND" to the $chr_and_pos string #
	# and if it is a SNP add "_SNP"                           #
	###########################################################
	if ($variant_type eq "indel"){$chr_and_pos = $chr_and_pos."_IND"}
	if ($variant_type eq "snp"){$chr_and_pos = $chr_and_pos."_SNP"}

	$chr_and_pos_original = $chr_and_pos; # store for when adding _1, _2 etc.  Put it here where it will retain the SNP and IND suffixes
	
	################################################################################
	# VCF:  Parse vcf_line_array(9)                                                      #
    # These are the genotype columns (one for each sample, starting at 9)          # 
    # Split the 9th item in the array at COLONS                                    #
    #                                                                              #
    # The 8th column should start with GT but it could be in any order so this     #
	# field is used to determine the order of the data in the 9th item             #
    #                                                                              #
    # Note DP occurs twice in the row, at the start and in this genotype section.  #
    # They are not necessarily the same.  We use the genotype one (DP_2)           #
    ################################################################################
    

     
	 #############################################################
	 # If the 9th item (or more for multiple samples) is present #
	 #############################################################
     if ($vcf_line_array_size > 8)
	{  
	 
        ################################################
        # VCF:  Read the FORMAT in vcf_line_array(8)   #
        #                                              #
        # This looks like this:  GT:AD:DP:GQ:PL        #
        # The genotype data is then read in this order #
        ################################################
        @vcf_format_field_array = split(":", $vcf_line_array[8]);
        $vcf_format_field_array_size = scalar @vcf_format_field_array;

			
		##########################################################
        # Get the genotype block (looks like 1/1:0,1:1:3:34,3,0) #
        ##########################################################
        if ($vcf_format_field_array_size > 0)
        {
            for($sample_count = 1; $sample_count <= $no_of_data_columns; $sample_count++)
			{
                $genotype_block_array[$sample_count] = $vcf_line_array[$sample_count + 8];
            }
        }

        ################################################################
        # VCF:  Now use order of fields in 'FORMAT'                    #
        # to parse the genotype data fields of                         # 
        # each sample (order may not be the same for every VCF file)   #
        ################################################################
        $genotypes_string = "";
        $GT_string = "";

        for($sample_count = 1; $sample_count <= $no_of_data_columns; $sample_count++)
		{
			@vcf_data_field_array = split(":",$genotype_block_array[$sample_count]);
        	$vcf_data_field_array_size = scalar @vcf_data_field_array;

            $GT_string = "";

            ################################################################
            # Now use order of fields in 'FORMAT' (vcf_format_field_array)               #
            # to parse the genotype data fields of                         # 
            # each sample (order may not be the same for every VCF file)   #
            ################################################################
            for($array_count = 0; $array_count < $vcf_format_field_array_size; $array_count++)
			{
				$myString = $vcf_format_field_array[$array_count];
                
				if ($myString eq "GT"){$GT_string = $vcf_data_field_array[$array_count]; last;}

            } # next array_count


			#####################################################
			# FreeBayes has a single dot . for no genotype      #
			# so there is no slash.  Convert "." to "./." first #
			# so it behaves as GATK VCF files                   #
			#####################################################	
			if ($GT_string eq "."){$GT_string = "./."}


			######################################################
			# platypus seems to have 1/0 as a possible genotype. #
			# Change these to 0/1                                #
			######################################################
			$allele_1 = substr($GT_string,0,1);
			$allele_2 = substr($GT_string,2,1);

			if (($GT_string ne "./.") && ($allele_1 > $allele_2))
			{
				$GT_string = $allele_2."/".$allele_1;
			}

			####################################################################
			# Make up genotypes_string which is a string of all the genotypes  #
			# at each position - across all samples                            #
			####################################################################
			if ($genotypes_string eq "") {$genotypes_string = $GT_string}
			else
			{$genotypes_string = $genotypes_string.",".$GT_string}

        } # next array_count

    } # if 9th element is present (data field with genotypes) VCF1


    #########################################################################
	# Store the locus_genotype string in the array locus_genotype_array_2d  #
	# This is for all variants, whether or not they go into the pooled file #
	#########################################################################

	$row_count = $row_count + 1;

	if ($variant_type eq "indel"){$locus_genotype_string = "$chr_and_pos"."_ind-"."$genotypes_string"}
	if ($variant_type eq "snp"){$locus_genotype_string = "$chr_and_pos"."_snp-"."$genotypes_string"}

	$locus_genotype_array_2d[$file_count][$row_count] = $locus_genotype_string;

}#read_vcf_line


sub read_line_if_first_vcf_file
{

	##############################################################################
	# First VCF file...                                                          #
	# Deal with situation where there is a SNP and an Indel at the same position #
	# (using the "_SNP" and "_IND" suffixes may sort this anyway)                #
	##############################################################################
	$variant_at_same_pos_count = 0;


	##########################################################################
	# First VCF file....                                                     #
	# If not defined then add it, and define genotypes_hash at the same time #
	##########################################################################
	if (not defined $position_hash{$chr_and_pos})
	{
		$position_hash{$chr_and_pos} = $ALT_base_overall;

		if (substr($chr_and_pos,-2) eq "_5")
		{
			print "Check 1a\n";
			print "chr_and_pos = $chr_and_pos\n";
			$answer=<STDIN>;
		}

		if (defined $genotypes_hash{$chr_and_pos})
		{
			print "chr_and_pos: $chr_and_pos\n";
			print "genotypes_hash is defined without position_hash being defined\n";
			print "This shouldn't happen\n\n";
			$answer=<STDIN>;
		}

		if (not defined $genotypes_hash{$chr_and_pos})
		{
			$genotypes_hash{$chr_and_pos} = $genotypes_string;
		}
	}
	else
	
	####################################################################
	# First VCF file...                                                #
	# If position_hash has been defined already at this position then  #
	# we need to store it at chr13_6322441_SNP_1, _2 etc               #
	#                                                                  #
	# SNP and IND at the same position will have a different key as    #
	# they will end in _SNP and _IND, but it's theretically possible   #
	# to have another position. The variable variant_at_same_pos_count #
	# counts any of these. There haven't been any so far               #
	####################################################################
	{
		$variant_at_same_pos_count = 0;
		$chr_and_pos_duplicates = $chr_and_pos;

		while(defined $position_hash{$chr_and_pos})
		{
			$variant_at_same_pos_count = $variant_at_same_pos_count + 1;
			$chr_and_pos = $chr_and_pos_original."_".$variant_at_same_pos_count;
			$chr_and_pos_duplicates = $chr_and_pos_original."_".$variant_at_same_pos_count;

			if ($variant_at_same_pos_count > $max_variants_at_same_pos){$max_variants_at_same_pos = $variant_at_same_pos_count}
		}

		print "Check 9: File 1.  chr_andPos = $chr_and_pos\n\n";
		$answer=<STDIN>;

		############################################################
		# Define position_hash and genotypes_hash at the same time #
		############################################################
		$position_hash{$chr_and_pos_duplicates} = $ALT_base_overall;

		if (substr($chr_and_pos_duplicates,-2) eq "_5")
		{
			print "Check 1b\n";
			print "chr_and_pos_duplicates = $chr_and_pos_duplicates\n";
			$answer=<STDIN>;
		}

		if (defined $genotypes_hash{$chr_and_pos_duplicates})
		{
			print "chr_and_pos_duplicates: $chr_and_pos_duplicates\n";
			print "genotypes_hash is defined without position_hash being defined\n";
			print "This shouldn't happen\n\n";
			$answer=<STDIN>;
		}

		$genotypes_hash{$chr_and_pos_duplicates} = $genotypes_string;

		#
		# Plan A
		# Rest chr_and_pos back to original #
		#
		$chr_and_pos = $chr_and_pos_original;

	}

	###########################################################
	# First VCF file...                                       #
	# Warn if position_hash is defined but not genotypes_hash #
	###########################################################
	if ((defined $position_hash{$chr_and_pos}) && (not defined $genotypes_hash{$chr_and_pos}))
	{
		print "FILE 1: If position_hash is defined then genotypes_hash should be defined\n\n";
		print "$position_hash{$chr_and_pos}\n";
		print "$genotypes_hash{$chr_and_pos}\n\n";
		$answer=<STDIN>;
	}

	####################################################
	# First VCF file...                                #
	# Update INFO field with details of variant caller #
	####################################################
	$INFO = $INFO.";CALLER=$variant_caller_array[$file_count]";
	$ID = $chr_and_pos;

	$single_line_new = "chr$chromosome\t$position\t$ID\t$REF_base\t$ALT_base_overall\t$QUAL\t$FILTER\t$INFO\t$FORMAT";

	for($sample_count = 1; $sample_count <= $no_of_data_columns; $sample_count++)
	{
		$single_line_new = "$single_line_new\t$genotype_block_array[$sample_count]";
    } # next array_count

	print POOLED "$single_line_new\n";

	$pooled_snp_count = $pooled_snp_count + 1;
	$found_in_pooled[$pooled_snp_count]="XXXXXX";

	$variants_in_file_1_count = $variants_in_file_1_count + 1;


	########################################################################
	# First VCF file...                                                    #
	# Store the locus_genotype string in the hash table pooled_hash        #
	# This is only for variants going into the pooled file                 #
	########################################################################
	#$locus_genotype_string = "$chr_and_pos"."-"."$genotypes_string";

	if ($variant_type eq "indel"){$locus_genotype_string = "$chr_and_pos"."_ind-"."$genotypes_string"}
	if ($variant_type eq "snp"){$locus_genotype_string = "$chr_and_pos"."_snp-"."$genotypes_string"}

	$pooled_hash{$locus_genotype_string} = "$file_count"."_"."$chr_and_pos"."_1";

}#read_line_if_first_vcf_file


################################################################################
# Subsequent VCF files...                                                      #
# If it is an existing chr and position then check that genotypes are the same #
# Check for basic chr_and_pos and also for any which have had _1, _2 added     #
# (up to the number of the files.)                                             #
################################################################################
sub read_line_if_subsequent_vcf_file
{
	if ($show eq "true") #1347
	{
		print "\n\nPosition: $position\tVariant caller: $variant_caller\n";
		print "Genotype: $genotypes_string\n";
		print "Check all duplicates for whether position_hash is defined:\n\n";

		for ($duplicates_count = 0; $duplicates_count < $no_of_vcf_files; $duplicates_count++)
		{
			if ($duplicates_count > 0){$chr_and_pos = $chr_and_pos_original."_".$duplicates_count}

			if (defined $position_hash{$chr_and_pos})
			{
				print "\t$duplicates_count \tposition_hash{$chr_and_pos} = $position_hash{$chr_and_pos}\n";
			}
			if (not defined $position_hash{$chr_and_pos})
			{
				print "\t$duplicates_count \tposition_hash{$chr_and_pos} = NOT DEFINED\n";
			}
			if (defined $genotypes_hash{$chr_and_pos})
			{
				print "\t\t$duplicates_count \tgenotype_hash{$chr_and_pos} = $genotypes_hash{$chr_and_pos}\n";
			}
			if (not defined $genotypes_hash{$chr_and_pos})
			{
				print "\t\t$duplicates_count \genotypes_hash{$chr_and_pos} = NOT DEFINED\n";
			}
		} # duplicates_count loop
		
		print "\n";

		$chr_and_pos = $chr_and_pos_original; # plan A
	}


	####################################################################################
	# Subsequent VCF files....                                                         #
	# Loop round all possibilities of the chr_and_pos_1, _2 etc to look for duplicates #
	#                                                                                  #
	# The loop is to LOOK for duplicates with a match                                  #
	####################################################################################
	$matching_genotype = "false";
	$matching_position = "false";

	for ($duplicates_count = 0; $duplicates_count < $no_of_vcf_files; $duplicates_count++)
	{
		if ($duplicates_count > 0){$chr_and_pos = $chr_and_pos_original."_".$duplicates_count}

		if ($position eq "39185825xx"){$show = "true"}

		if ($show eq "truex")
		{
			print "Check 8 - subsequent VCF files.  \n";
			print "\tduplicates_count: $duplicates_count\tchr_and_pos: $chr_and_pos\n";
			if (defined $position_hash{$chr_and_pos}){print "\tposition_hash{$chr_and_pos} defined\n";}
			else
			{print "\tposition_hash{$chr_and_pos} not defined\n";}
			print "\tmatching_position :$matching_position     matching_genotype: $matching_genotype\n\n";
			$answer=<STDIN>;
		} # SHOW


		######################################################################
		# If key is found here, then if matching_genotype is true then       #
		# they are the same and don't need to be re-added                    #
		######################################################################
		if (defined $position_hash{$chr_and_pos})
		{
			$matching_position = "true"; # start assuming this

			# debug
			if (not defined $genotypes_hash{$chr_and_pos}){print "genotypes_hash{$chr_and_pos} should be defined 1\n"; $answer=<STDIN>;}

			$already_defined_array[$file_count] = $already_defined_array[$file_count] + 1;

			if ($show eq "truex")
			{
				print "Check 6 - position_hash is already defined\n";
				print "\tduplicates_count = $duplicates_count\n";
				print "\tposition_hash{$chr_and_pos} = $position_hash{$chr_and_pos}\n";
				$answer=<STDIN>;
			}

			############################################################
			# If not defined then you can end the loop here, as higher #
			# values shouldn't be defined if lower values are not.     #
			############################################################
			if (not defined $genotypes_hash{$chr_and_pos})
			{
				#last; #add later
			}

		
			if ($show eq "truex")
			{
				if (defined $genotypes_hash{$chr_and_pos})
				{
					print "Check 7 -  genotypes_hash should also be already defined (and it is)\n";
					print "\tduplicates_count = $duplicates_count\n";
					print "\tgenotypes_hash{$chr_and_pos} = $genotypes_hash{$chr_and_pos}\n";
				}
				else
				{
					print "Check 7 - ERROR!!!  genotypes_hash should also be already defined\n";
					print "\tduplicates_count = $duplicates_count\n";
					print "\tgenotypes_hash{$chr_and_pos} = NOT DEFINED\n";
				}
				$answer=<STDIN>;
			}

			####################################################################
			# Subsequent VCF files...                                          #
			# If genotypes don't match then it needs to be added as a new line #
			# WHAT HAPPENS TO THE POSITION_HASH FOR THIS LINE?
			# 
			# Try calling it _1, _2 in output file for debugging.....??
			######################################################

			if ($genotypes_string eq $genotypes_hash{$chr_and_pos})
			{
				$matching_genotype = "true";
			}

			if ($show eq "truex")
			{
				print "Check - do genotypes match?\n";
				print "\tGenotypes string:           \t\t$genotypes_string\n";
				print "\tGenotypes_hash($chr_and_pos): \t$genotypes_hash{$chr_and_pos}\n";
				print "\tMatching genotype: $matching_genotype\n\n";

				if ($genotypes_string eq $genotypes_hash{$chr_and_pos})
				{
					print "Yes they do match\n\n";
				}
				else
				{print "NO they don't match\n\n";}
			}

			#########################################################################################
			# If genotypes do match then this is an exact duplicate so we don't need to do anything #
			# although I suppose we could record that it was found by more than one variant caller? #
			#########################################################################################

		} # if position_hash is defined already

		################################################################################
		# If we have found a position and genotype matching then no need to go further #
		################################################################################
		if (($matching_position eq "true") && ($matching_genotype eq "true")){ last; }

	} # duplicates_count loop to check for matching position and genotype


	if ($show eq "true")
	{
		print "\nCheck - After looping through all possible duplicates we find that:\n";
		print "\tmatching_position :$matching_position     matching_genotype: $matching_genotype\n\n";

		if (($matching_position eq "true") && ($matching_genotype eq "true"))
		{
			print "\tBoth are true so we don't need to write this to the pooled file\n\n";
		}
		if (($matching_position eq "true") && ($matching_genotype eq "false"))
		{
			print "\tWe have a new genotype at an old position, so we want to add this to the pooled file\n\n";
		}
		if (($matching_position eq "false") && ($matching_genotype eq "false"))
		{
			print "\tNeither are true so we need to write this new position/genotype to the pooled file\n\n";
		}
		if (($matching_position eq "true") && ($matching_genotype eq "false"))
		{
			print "\tERROR!!! This shouldn't happen\n\n";
		}
		$answer=<STDIN>;
	}


	#####################
	# Reset to original #
	#####################
	$chr_and_pos = $chr_and_pos_original;


	#################################################################
	# After looping through all the possible duplicates, if there   #
	# is no matching genotype then it has to be added as a new line #
	# to the output pooled VCF file                                 #
	#################################################################
	if (($matching_genotype eq "false") && ($matching_position eq "true"))
	{
		if ($position eq "39013405xx"){$show = "true"}else{$show = "false"} ## current check TEMP!!!!!

		if ($show eq "true")
		{
			print "Check 7a -  position matches but not genotype\n";
			$answer=<STDIN>;
		}


		#############################################################################
		# Now we know there is no matching variant                                  #
		# we need to add the new one with a different chr_and_pos position_hash key #
		#                                                                           #
		# So find the first free one (start at 1)                                   #
		#############################################################################
		$variant_at_same_pos_count = 0;

		while(defined $position_hash{$chr_and_pos})
		{
			if ($show eq "true")
			{
				print "Looking for free position_hash slot...\n";
				print "position_hash{$chr_and_pos}: $position_hash{$chr_and_pos}\n\n";
			}

			$variant_at_same_pos_count = $variant_at_same_pos_count + 1;
			$chr_and_pos = $chr_and_pos_original."_".$variant_at_same_pos_count;

			if ($variant_at_same_pos_count > $max_variants_at_same_pos){$max_variants_at_same_pos = $variant_at_same_pos_count}
		}

		if ($show eq "true")
		{
			if (not defined position_hash{$chr_and_pos})
			{
				print "This chr_and_pos $chr_and_pos is not used so we will use that. ($duplicates_count)\n";
			}
		}

		if ($show eq "true")
		{
			print "Check 7b -  next free duplicate found...\n";
			print "\tduplicates_count = $duplicates_count\n";
			print "\t\chr_and_pos = $chr_and_pos\n\n";

			$answer=<STDIN>;
		}

		##############################################################
		# So here we define the new (different) genotypes_hash       #
		# which needs a new position_hash otherwise it will be lost? #
		##############################################################
		$position_hash{$chr_and_pos} = $ALT_base_overall;

		if (defined $genotypes_hash{$chr_and_pos})
		{
			print "Check 3.\n";
			print "chr_and_pos: $chr_and_pos\n";
			print "genotypes_hash is defined without position_hash being defined\n";
			print "This shouldn't happen\n\n";
			$answer=<STDIN>;
		}

		# Store new genotype in genotypes_hash
		$genotypes_hash{$chr_and_pos} = $genotypes_string;


		######################################################
		# Update INFO field with details of variant caller   #
		# and also record a mismatch (with file 1 only note) #
		######################################################
		$INFO = $INFO.";CALLER=$variant_caller_array[$file_count];MISMATCH=YES";
		$ID = $chr_and_pos;

		$single_line_new = "chr$chromosome\t$position\t$ID\t$REF_base\t$ALT_base_overall\t$QUAL\t$FILTER\t$INFO\t$FORMAT";

		for($sample_count = 1; $sample_count <= $no_of_data_columns; $sample_count++)
		{
			$single_line_new = "$single_line_new\t$genotype_block_array[$sample_count]";
        }

        ###################################################################
        # Subsequent VCF files...                                         #
        # Add new line with mismatched genotype to the bottom of the file #
        ###################################################################
		print POOLED "$single_line_new\n";

		$pooled_snp_count = $pooled_snp_count + 1;

		$lines_added_array[$file_count] = $lines_added_array[$file_count] + 1;

		$defined_but_different_array[$file_count] = $defined_but_different_array[$file_count] + 1;

	
		########################################################################
		# Subsequent VCF files...                                              #
		# Store the locus_genotype string in the hash table pooled_hash        #
		# This is only for variants going into the pooled file                 #
		########################################################################
		#$locus_genotype_string = "$chr_and_pos"."-"."$genotypes_string";

		if ($variant_type eq "indel"){$locus_genotype_string = "$chr_and_pos"."_ind-"."$genotypes_string"}
		if ($variant_type eq "snp"){$locus_genotype_string = "$chr_and_pos"."_snp-"."$genotypes_string"}

		$pooled_hash{$locus_genotype_string} = "$file_count"."_"."$chr_and_pos";


		########################################
		# Set this array to any string for now #
		########################################
		$found_in_pooled[$pooled_snp_count]="XXXXXX";

	} # subsequent VCF files -- if (($matching_genotype eq "false") && ($matching_position eq "true"))



	#################################################################
	# If there is no matching position then it has to be added as a #
	# new line to the output pooled VCF file                        #
	#################################################################
	if ($matching_position eq "false")
	{
		if (defined $genotypes_hash{$chr_and_pos})
		{
			print "Check 3 <<<<<<< ERROR!\n";
			print "\tchr_and_pos: $chr_and_pos\n";
			print "\tMatching position is false, so there should be no defined genotypes_hash\n\n";
			$answer=<STDIN>;
		}


		if ($show eq "true")
		{
			if (not defined $genotypes_hash{$chr_and_pos})
			{
				print "Check 3b.\n";
				print "\tchr_and_pos: $chr_and_pos\n";
				print "\tMatching position is false (i.e position_hash is not defined), and there is no defined genotypes_hash\n";
				print "\tSo we need to add this line to the pooled file\n\n";

				$answer=<STDIN>;
			}
		} # show

		###################################################################
        # Subsequent VCF files...                                         #
        # Add new line with mismatched genotype to the bottom of the file #
        ###################################################################
        $INFO = $INFO.";CALLER=$variant_caller_array[$file_count]";
        $ID = $chr_and_pos;

        $single_line_new = "chr$chromosome\t$position\t$ID\t$REF_base\t$ALT_base_overall\t$QUAL\t$FILTER\t$INFO\t$FORMAT";

        for($sample_count = 1; $sample_count <= $no_of_data_columns; $sample_count++)
		{
			$single_line_new = "$single_line_new\t$genotype_block_array[$sample_count]";
        }

		print POOLED "$single_line_new\n";

		$pooled_snp_count = $pooled_snp_count + 1;

		$lines_added_array[$file_count] = $lines_added_array[$file_count] + 1;

		########################################################################
		# Subsequent VCF files...   (brand new position)                       #
		# Store the locus_genotype string in the hash table pooled_hash        #
		# This is only for variants going into the pooled file                 #
		########################################################################
		if ($variant_type eq "indel"){$locus_genotype_string = "$chr_and_pos"."_ind-"."$genotypes_string"}
		if ($variant_type eq "snp"){$locus_genotype_string = "$chr_and_pos"."_snp-"."$genotypes_string"}

		$pooled_hash{$locus_genotype_string} = "$file_count"."_"."$chr_and_pos"."_2"; # Not sure what this -2 is.

		##############################################################################
		# We also need to store this position in position_hash AND in genotypes_hash #
		##############################################################################
		$position_hash{$chr_and_pos} = $ALT_base_overall;
		$genotypes_hash{$chr_and_pos} = $genotypes_string;

		if ($show eq "true")
		{
			print "Check 3c.\n";
			print "\tSo we have now added this line to the pooled file\n\n";

			print "\t\tposition_hash($chr_and_pos):  $position_hash{$chr_and_pos}\n";
			print "\t\tgenotypes_hash($chr_and_pos): $genotypes_hash{$chr_and_pos}\n\n";

			$answer=<STDIN>;
	
		} # show

	} # $matching_position eq "false" (i.e. if NO variants at this position then you must write this line to the POOLED file)


	############################
	# Reset chr_and_pos plan A #
	############################
	$chr_and_pos = $chr_and_pos_original;

	######################################################
	# Subsequenct VCF files...                           #
	# If positon_hash is defined, we also need to check  #
	# chr_and_pos's that have had _1, _2 added to them   #
	######################################################

}#read_line_if_subsequent_vcf_file



############################################################
# Parsing of #CHROM data line to get a list of input files #
############################################################
sub get_samples_from_chrom_line
{
	@chrom_line_array = split(/\s+/,$single_line);
	$chrom_line_array_size = (scalar @chrom_line_array) - 1;
	$no_of_samples_chrom_line = $chrom_line_array_size - 8; # ignore first columns 0-9
	$sample_name_string[$file_count] ="";

	for ($array_count = 9; $array_count <= $chrom_line_array_size; $array_count++)
	{
		# Sample names go into sample_name_array
		$sample_name_array[$array_count - 8] = $chrom_line_array[$array_count];

		####################################################################
		# Also store sample names in a string for comparison between files #
		####################################################################
		$sample_name_string[$file_count] = $sample_name_string[$file_count].$chrom_line_array[$array_count]."\t";
	}
}#get_samples_from_chrom_line