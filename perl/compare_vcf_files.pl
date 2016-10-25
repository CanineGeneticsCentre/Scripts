#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	compare_vcf_files  						                            #     
#									                                    #
#	THIS PERL SCRIPT COMPARES VCF FILES                                 #
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

my $version							= "6";

my $debug							= "false";

my $answer							= "";
my $pooled_vcf_file					= "";
my $single_line						= "";
my $start_storing					= "false";
my $chromosome						= "";
my $position						= "";
my $chr_and_pos						= ""; #chromosome and position linked e.g. 12_6456378
my $REF_base						= "";
my $ALT_base_overall_1				= "";
my $ALT_base_overall_2				= "";
my $QUAL							= "";
my $FILTER							= "";
my $run_title						= "";
my $genotypes_1_string				= "";
my $genotypes_2_string				= "";
my $depth_2_string					= "";
my $depth							= "";
my $myString						= "";
my $GT_string						= "";
my $DP_string 						= "";
my $variant_type					= ""; # "snp" or "indel"

#files
my $vcf_file_1						= "";
my $vcf_file_2						= "";
my $command_log						= "";
my $discordant_variants_file		= "";
my $concordant_variants_file		= "";
my $only_in_1_variants_file			= "";
my $only_in_2_variants_file			= "";
my $in_both_variants_file			= "";
my $discordant_variants_details		= "";

#Boolean
my $all_genotypes_match				= ""; #"true" or "false"
my $sample_name_mismatch			= ""; #"true" or "false"  - checks if sample names are the same for the two VCF files

#Counts of variant types
my $defined_count_1					= 0; # no of counts in pos_1_hash
my $pos_1_hash_size					= 0; # no of keys in pos_1_hash
my $pos_2_hash_size					= 0; # no of keys in pos_2_hash
my $concordant_position_count		= 0; # Each position will have a number of genotpyes (same as no of samples)
my $discordant_position_count		= 0;
my $discordant_snps_count			= 0; # This refers to position not individual genotype
my $discordant_indels_count			= 0; # This refers to position not individual genotype
my $only_in_file_1_count			= 0;
my $only_in_file_2_count			= 0;
my $variants_in_file_1_count		= 0;
my $variants_in_file_2_count		= 0;
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
my $total_variants_in_file_1		= 0;
my $total_variants_in_file_1_check	= 0;
my $total_variants_in_file_2		= 0;
my $file_1_and_2_count				= 0;
my $file_1_only_calc				= 0;
my $file_count						= 0;
my $check_genotypes_count			= 0; # All genotypes counted

#Other
my $no_of_vcf_files					= 0;
my $line_count						= 0;
my $chrom_line_array_size			= 0;
my $array_count						= 0;
my $no_of_samples_chrom_line		= 0;
my $array_size_1					= 0;
my $array_size_8					= 0;
my $array_size_9					= 0;
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
my @sample_name_array_1				= ();
my @sample_name_array_2				= ();
my @myArray1						= ();
my @myArray8						= ();
my @myArray9						= ();
my @ALT_allele_array_1				= ();
my @ALT_allele_array_2				= ();
my @genotype_block_array			= ();
my @genotype_all_1_array			= ();
my @genotype_all_2_array			= ();
my @genotypes_1_array				= (); # split genotype string into individual genotypes whilst comparing
my @genotypes_2_array				= (); # split genotype string into individual genotypes whilst comparing
my @depth_2_array					= (); # split depth string into individual genotypes whilst comparing

my %pos_1_hash						= (); # hash table for positions in VCF file 1
my %pos_2_hash						= (); # hash table for positions in VCF file 2
my %genotypes_1_hash				= (); # hash table for genotypes in VCF file 1

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "              compare_vcf_files   \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "    - This program compares two VCF files (e.g. produced by different variant callers)\n";


print color 'reset';

&print_message("Here is a list of VCF files in this directory","message");
print "\n";
system("ls -1 *.vcf");
print "\n";

&print_message("First VCF file:   ","input");

$answer=<STDIN>;
chomp $answer;                           if ($answer eq ""){$answer = "AS_HC_UG.vcf"}
$vcf_file_1 = $answer;

&print_message("Second VCF file: ","input");

$answer=<STDIN>;
chomp $answer;                          if ($answer eq ""){$answer = "AS_HC_HC.vcf"}

$vcf_file_2 = $answer;

$run_title = &get_prefix("$vcf_file_1")."_vs_".&get_prefix("$vcf_file_2");

$command_log = "compare_vcf_files_"."$run_title"."_command_log.out";

$discordant_variants_file = "discordant_variants_"."$run_title".".out";
$concordant_variants_file = "concordant_variants_"."$run_title".".out";
$only_in_1_variants_file = "only_in_1_variants_"."$run_title".".out";
$only_in_2_variants_file = "only_in_2_variants_"."$run_title".".out";
$discordant_variants_details = "discordant_variants_details_"."$run_title".".out";



###################################
# Open the various output files   #
###################################
open (DVF, ">$discordant_variants_file") || die "Cannot open $discordant_variants_file";
open (CVF, ">$concordant_variants_file") || die "Cannot open $concordant_variants_file";
open (O1VF, ">$only_in_1_variants_file") || die "Cannot open $only_in_1_variants_file";
open (O2VF, ">$only_in_2_variants_file") || die "Cannot open $only_in_2_variants_file";
open (DVD, ">$discordant_variants_details") || die "Cannot open $discordant_variants_details";

open (DEBUG, ">debug.vcf") || die "Cannot open file";

#########################
# Open command log file #
#########################
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "compare_vcf_files  version $version\n\n";	
print COMMAND_LOG "VCF file 1: $vcf_file_1\n";
print COMMAND_LOG "VCF file 2: $vcf_file_2\n\n";		

print DVD "Chr\tPosition\tType\tVariant_ID\tSample\tDepth\t$vcf_file_1\t$vcf_file_2\n";


##########################
# Open the VCF file 1    #
##########################
open (VCF1, "$vcf_file_1") || die "Cannot open $vcf_file_1";

$line_count = 0;
$variants_in_file_1_count = 0;

while ($single_line = <VCF1> ) 
{
	chomp $single_line;
	&chomp_all ($single_line);

	$line_count = $line_count + 1;

	######################################################
	# VCF1:                                              #
	# Only do this next bit if 'start_storing is True    #   
	# i.e. if the line containing #CHROM has been found  #
	######################################################
	if ($start_storing eq "true")
	{
		#########################################################
        # Split whole line at TABs into an array myArray1 (1-9) #
        #########################################################
		@myArray1 = split (/\t/,$single_line);
		$array_size_1 = scalar @myArray1;
		$no_of_data_columns = $array_size_1 - 9;


		##################################################
        # VCF: Get first fixed columns of the VCF file   #
        # CHR, POS, REF, ALT, QUAL, FILTER               #
        ##################################################
		$chromosome = $myArray1[0];
        $position = $myArray1[1];
        $REF_base = $myArray1[3];
        $ALT_base_overall_1 = $myArray1[4];
        $QUAL = $myArray1[5];
        $FILTER = $myArray1[6];

        $variants_in_file_1_count = $variants_in_file_1_count + 1;


        ############################################
        # VCF1:  Various conversions on chromosome #
        ############################################
        # Convert chrX to chr39 for dog 
		if ($chromosome eq "chrX"){$chromosome = "chr39"}
         
		# Remove string 'chr'
        if (substr($chromosome,0,3) eq "chr"){$chromosome = substr($chromosome,3,99)}
		
		# If chromosome is an integer and less than 10 then add a leading zero
		if ($chromosome =~ /^\d+?$/){if ($chromosome < 10){$chromosome = "0".$chromosome}}


		#################################################
		# VCF1: Split ALT_base_overall_1 at commas      #
		# to get all the different alleles in the array #
		#################################################
		@ALT_allele_array_1 = split(",",$ALT_base_overall_1);
		$no_of_alt_alleles_1 = (scalar @ALT_allele_array_1);


		###########################################################
		# VCF1: Decide whether it is a SNP or an Indel VCF1       #
		###########################################################
		$variant_type = "snp";
		for ($allele_count = 0; $allele_count <$no_of_alt_alleles_1; $allele_count++)
		{
			if (length($REF_base) != length($ALT_allele_array_1[$allele_count]))
			{
				$variant_type = "indel";
			}
		}

		#########################
		# Count SNPs and Indels #
		#########################
		if ($variant_type eq "snp"){$snps_in_file_1_count = $snps_in_file_1_count + 1}
		if ($variant_type eq "indel"){$indels_in_file_1_count = $indels_in_file_1_count + 1}
		if (($variant_type ne "snp") && ($variant_type ne "indel")){$unknowns_in_file_1_count = $unknowns_in_file_1_count + 1}



        ###########################################################
		# VCF1: Store in hash table (faster than array)           #
		# If there are two variants for one position, store the   #
		# subsequent ones one as 12_12134562_1, 12_12134562_2 etc #
		###########################################################

		$chr_and_pos = $chromosome."_".$position;


		###########################################################
		# If it is an indel add "_IND" to the $chr_and_pos string #
		# and if it is a SNP add "_SNP"                           #
		###########################################################
		if ($variant_type eq "indel"){$chr_and_pos = $chr_and_pos."_IND"}
		if ($variant_type eq "snp"){$chr_and_pos = $chr_and_pos."_SNP"}


		##############################################################################
		# Deal with situation where there is a SNP and an Indel at the same position #
		# (using the "_SNP" and "_IND" suffixes may sort this anyway)                #
		##############################################################################
		$variant_at_same_pos_count = 0;

		if (not defined $pos_1_hash{$chr_and_pos})
		{
			$pos_1_hash{$chr_and_pos} = $ALT_base_overall_1;
		}
		else
		{
			while(defined $pos_1_hash{$chr_and_pos})
			{
				$variant_at_same_pos_count = $variant_at_same_pos_count + 1;
				$chr_and_pos = $chr_and_pos."_".$variant_at_same_pos_count;

				if ($variant_at_same_pos_count > $max_variants_at_same_pos){$max_variants_at_same_pos = $variant_at_same_pos_count}
			}
			$pos_1_hash{$chr_and_pos} = $ALT_base_overall_1;
		}

		
		
		################################################################################
		# VCF1: Parse myArray1(9)                                                      #
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
         if ($array_size_1 > 8)
		{  
		 
            ################################################
            # VCF1: Read the FORMAT in myArray1(8)         #
            #                                              #
            # This looks like this:  GT:AD:DP:GQ:PL        #
            # The genotype data is then read in this order #
            ################################################
            @myArray8 = split(":", $myArray1[8]);
            $array_size_8 = scalar @myArray8;

 			
 			##########################################################
            # Get the genotype block (looks like 1/1:0,1:1:3:34,3,0) #
            ##########################################################
            if ($array_size_8 > 0)
            {
                for($sample_count = 1; $sample_count <= $no_of_data_columns; $sample_count++)
				{
                    $genotype_block_array[$sample_count] = $myArray1[$sample_count + 8];
                }
            }


            ################################################################
            # VCF1: Now use order of fields in 'FORMAT' (myArray8)         #
            # to parse the genotype data fields of                         # 
            # each sample (order may not be the same for every VCF file)   #
            ################################################################
            $genotypes_1_string = "";
            $GT_string = "";

            for($sample_count = 1; $sample_count <= $no_of_data_columns; $sample_count++)
			{
				@myArray9 = split(":",$genotype_block_array[$sample_count]);
            	$array_size_9 = scalar @myArray9;

                $GT_string = "";

				$GT_string = $myArray9[0]; # Assume GT is first in block
				

				if ($genotypes_1_string eq "") {$genotypes_1_string = $GT_string}
				else
				{$genotypes_1_string = $genotypes_1_string.",".$GT_string}

            } # next array_count

            $genotypes_1_hash{$chr_and_pos} = $genotypes_1_string;

            #$answer=<STDIN>;

        } # if 9th element is present (genotypes) VCF1

	} # if start_storing is true


	#########################################################
    # VCF1: Has the start of the CHR data been reached yet? #
	# The line before the actual data starts with '#CHROM'  #
	# This part of the code checks that line, which has the #
	# (very important) headers for the sample columns.      #
    #########################################################
	
	if (index($single_line,"#CHROM") > -1)
	{
		$start_storing = "true";
		
		&print_message("Reading VCF file 1: $vcf_file_1","message");


		############################################################
		# Parsing of #CHROM data line to get a list of input files #
		############################################################
		@chrom_line_array = split(/\s+/,$single_line);
		$chrom_line_array_size = (scalar @chrom_line_array) - 1;
		$no_of_samples_chrom_line = $chrom_line_array_size - 8; # ignore first columns 0-9
		
		for ($array_count = 9; $array_count <= $chrom_line_array_size; $array_count++)
		{
			# Sample names go into sample_name_array_1
			$sample_name_array_1[$array_count - 8] = $chrom_line_array[$array_count];
			$file_count = $array_count - 8;
		}

	} # end of if CHROM line
} # reading VCF file 1

close VCF1;

$total_variants_in_file_1 = $variants_in_file_1_count;

&print_both("Finished reading VCF file 1: $vcf_file_1\tNumber of variants: $total_variants_in_file_1\tNumber of samples: $no_of_samples_chrom_line\n");


##########################
# Get size of hash array #
##########################
$pos_1_hash_size = scalar keys %pos_1_hash;


######################
# Reading VCF file 2 #
######################
open (VCF2, "$vcf_file_2") || die "Cannot open $vcf_file_2";
$line_count = 0;

$start_storing = "false";

while ($single_line = <VCF2> ) 
{
	chomp $single_line;
	&chomp_all ($single_line);

	$line_count = $line_count + 1;


	######################################################
	# VCF2:                                              #
	# Only do this next bit if 'start_storing is True    #   
	# i.e. if the line containing #CHROM has been found  #
	######################################################
	if ($start_storing eq "true")
	{
		#########################################################
        # Split whole line at TABs into an array myArray1 (1-9) #
        #########################################################
		@myArray1 = split (/\t/,$single_line);
		$array_size_1 = scalar @myArray1;
		$no_of_data_columns = $array_size_1 - 9;


		#############################################
        # Get first fixed columns of the VCF file   #
        # CHR, POS, REF, ALT, QUAL, FILTER          #
        #############################################
		$chromosome = $myArray1[0];
        $position = $myArray1[1];
        $REF_base = $myArray1[3];
        $ALT_base_overall_2 = $myArray1[4];
        $QUAL = $myArray1[5];
        $FILTER = $myArray1[6];

        $variants_in_file_2_count = $variants_in_file_2_count + 1;

		###########################################
        # VCF2: Various conversions on chromosome #
        ###########################################
        # Convert chrX to chr39 for dog 
		if ($chromosome eq "chrX"){$chromosome = "chr39"}
         
		# Remove string 'chr'
        if (substr($chromosome,0,3) eq "chr"){$chromosome = substr($chromosome,3,99)}
		
		# If chromosome is an integer and less than 10 then add a leading zero
		if ($chromosome =~ /^\d+?$/){if ($chromosome < 10){$chromosome = "0".$chromosome}}


		################################################
		# VCF2: Split ALT_base_overall_2 at commas     #
		# to get all the different alleles in an array #
		################################################
		@ALT_allele_array_2 = split(",",$ALT_base_overall_2);
		$no_of_alt_alleles_2 = (scalar @ALT_allele_array_2);


		###########################################################
		# VCF2: Decide whether it is a SNP or an Indel            #
		###########################################################
		$variant_type = "snp";
		for ($allele_count = 0; $allele_count <$no_of_alt_alleles_2; $allele_count++)
		{
			if (length($REF_base) != length($ALT_allele_array_2[$allele_count]))
			{
				$variant_type = "indel";
			}
		}

		###############################
		# VCF2: Count SNPs and Indels #
		###############################
		if ($variant_type eq "snp"){$snps_in_file_2_count = $snps_in_file_2_count + 1}
		if ($variant_type eq "indel"){$indels_in_file_2_count = $indels_in_file_2_count + 1}
		if (($variant_type ne "snp") && ($variant_type ne "indel")){$unknowns_in_file_2_count = $unknowns_in_file_2_count + 1}


        ###########################################################
		# VCF2: Store in hash table (faster than array)           #
		# If there are two variants for one position, store the   #
		# subsequent ones one as 12_12134562_1, 12_12134562_2 etc #
		##########################################################
		$chr_and_pos = $chromosome."_".$position;


		###########################################################
		# If it is an indel add "_IND" to the $chr_and_pos string #
		###########################################################
		if ($variant_type eq "indel"){$chr_and_pos = $chr_and_pos."_IND"}
		if ($variant_type eq "snp"){$chr_and_pos = $chr_and_pos."_SNP"}


		##############################################################################
		# VCF2:                                                                      #
		# Deal with situation where there is a SNP and an Indel at the same position #
		# (using the "_SNP" and "_IND" suffixes may sort this anyway)                #
		##############################################################################

		$variant_at_same_pos_count = 0;
		if (not defined $pos_2_hash{$chr_and_pos})
		{
			$pos_2_hash{$chr_and_pos} = $ALT_base_overall_2;
		}
		else
		{
			while(defined $pos_2_hash{$chr_and_pos})
			{
				$variant_at_same_pos_count = $variant_at_same_pos_count + 1;
				$chr_and_pos = $chr_and_pos."_".$variant_at_same_pos_count;

				if ($variant_at_same_pos_count > $max_variants_at_same_pos){$max_variants_at_same_pos = $variant_at_same_pos_count}
			}
			$pos_2_hash{$chr_and_pos} = $ALT_base_overall_2;
		}



		################################################################################
		# VCF2: Parse myArray1(9)                                                      #
        # These are the genotype columns (one for each sample)                         # 
        # Split the 9th item in the array at COLONS                                    #
        #                                                                              #
        # The 8th column should start with GT but it could be in any order so this     #
		# field is used to determine the order of the data in the 9th item             #
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


            ##########################################################
            # Get the genotype block (looks like 1/1:0,1:1:3:34,3,0) #
            ##########################################################
            if ($array_size_8 > 0)
            {
                for($sample_count = 1; $sample_count <= $no_of_data_columns; $sample_count++)
				{
                    $genotype_block_array[$sample_count] = $myArray1[$sample_count + 8];
                }
            }


            ################################################################
            # VCF2: Now use order of fields in 'FORMAT' (myArray8)         #
            # to parse the genotype block fields of                        # 
            # each sample (order may not be the same for every VCF file)   #
            ################################################################
            $genotypes_2_string = "";
            $depth_2_string = "";
            $GT_string = "";
            $DP_string = "";

            for($sample_count = 1; $sample_count <= $no_of_data_columns; $sample_count++)
			{
				@myArray9 = split(":",$genotype_block_array[$sample_count]);
            	$array_size_9 = scalar @myArray9;

                $GT_string = "";

				$GT_string = $myArray9[0]; # Assume GT is first in block
				

				if ($genotypes_2_string eq "") {$genotypes_2_string = $GT_string}
				else
				{$genotypes_2_string = $genotypes_2_string.",".$GT_string}


				################################
				# If genotype is or is not ./. #
				################################
				$DP_string = "0";
				if ($genotype_block_array[$sample_count] ne "./.")
				{
					if ($array_size_9 > 2){$DP_string = $myArray9[2];} # Assume DP is third in block.  Array size checks this if not there
				}
				if ($genotype_block_array[$sample_count] eq "./.")
				{
					$DP_string = "0"; # No depth string if genotype is ./.
				}

				################################################
				# DP_string can be just a dot, so replace by 0 #
				################################################
				if ($DP_string eq "."){$DP_string = "0"}


				################################
				# Now make up string of depths #
				################################
				if ($depth_2_string eq "") {$depth_2_string = $DP_string}
				else
				{$depth_2_string = $depth_2_string.",".$DP_string}

            } # next array_count

        } # if 9th element is present (genotypes)  VCF2


		######################################################
        # Check if position is in hash array from VCF file 1 #
		######################################################
		$chr_and_pos = $chromosome."_".$position;

		###########################################################
		# If it is an indel add "_IND" to the $chr_and_pos string #
		###########################################################
		if ($variant_type eq "indel"){$chr_and_pos = $chr_and_pos."_IND"}
		if ($variant_type eq "snp"){$chr_and_pos = $chr_and_pos."_SNP"}

		if (defined $pos_1_hash{$chr_and_pos})
		{
			$defined_count_1 = $defined_count_1 + 1;

			$ALT_base_overall_1 = $pos_1_hash{$chr_and_pos};

			#################################################
			# Store depth information (where is this from?) #
			#################################################
			if ($depth_2_string ne "")
			{
				@depth_2_array = split(",",$depth_2_string);
			}
			


			##############################################
			# VCF2: Check each genotype in this position #
			##############################################

			########################################################
			# Split the genotype strings into individual genotypes #
			########################################################
			@genotypes_1_array = split(",",$genotypes_1_hash{$chr_and_pos});
			@genotypes_2_array = split(",",$genotypes_2_string);

			$array_size_genotypes_1 = scalar @genotypes_1_array;


			###############################################
			# Loop through all genotypes at this position #
			###############################################
			$all_genotypes_match = "true";
			for ($array_count = 0; $array_count < $array_size_genotypes_1; $array_count++)
			{

				if ($depth_2_string ne ""){$depth = $depth_2_array[$array_count]}
				else
				{$depth= "x"}

				$check_genotypes_count = $check_genotypes_count + 1;

				################################
				# If the genotypes don't match #
				################################
				if ($genotypes_1_array[$array_count] ne $genotypes_2_array[$array_count])
				{
					############################################################################
					# Write individual discordant genotypes to discordant variant details file #
					############################################################################
					print DVD "$chromosome\t$position\t$variant_type\t$chr_and_pos\t$sample_name_array_1[$array_count+1]\t$depth\tg$genotypes_1_array[$array_count]\tg$genotypes_2_array[$array_count]\n";

					##############################
					# Record type of discordancy #
					##############################
					if (($genotypes_1_array[$array_count] eq "./.") || ($genotypes_2_array[$array_count] eq "./."))
					{
							$discordant_one_null_vcount = $discordant_one_null_vcount + 1;
					}
					if (($genotypes_1_array[$array_count] ne "./.") && ($genotypes_2_array[$array_count] ne "./."))
					{
							$discordant_neither_null_vcount = $discordant_neither_null_vcount + 1;
					}
					$all_genotypes_match = "false";
				} # loop through genotypes at this position


				################################
				# If the genotypes do match    #
				################################
				if ($genotypes_1_array[$array_count] eq $genotypes_2_array[$array_count])
				{
					##############################
					# Record type of concordancy #
					##############################
					if ($genotypes_1_array[$array_count] eq "./.")
					{
							$concordant_null_vcount = $concordant_null_vcount + 1;
					}
					if ($genotypes_1_array[$array_count] ne "./.")
					{
							$concordant_neither_null_vcount = $concordant_neither_null_vcount + 1;
					}
				}
			}

			if ($all_genotypes_match eq "true"){$concordant_position_count = $concordant_position_count + 1;}
			if ($all_genotypes_match eq "false")
			{
				$discordant_position_count = $discordant_position_count + 1;
				if ($variant_type eq "indel"){$discordant_indels_count = $discordant_indels_count + 1}
				if ($variant_type eq "snp"){$discordant_snps_count = $discordant_snps_count + 1}

				print DVF "$chr_and_pos\t$ALT_base_overall_1\t$ALT_base_overall_2\n";
				print DVF "$position\t$genotypes_1_hash{$chr_and_pos}\t$genotypes_2_string\n";
			}

		} # if in hash VCF 1

		if (not defined $pos_1_hash{$chr_and_pos})
		{
			$only_in_file_2_count = $only_in_file_2_count + 1;

			print O2VF "$chromosome\t$position\t$chr_and_pos\t$ALT_base_overall_2\n";
		}


	} # if start_storing is true


	#########################################################
    # Has the start of the CHR data been reached yet?       #
	# The line before the actual data starts with '#CHROM'  #
	# This part of the code checks that line, which has the #
	# (very important) headers for the sample columns.      #
    #########################################################
	
	if (index($single_line,"#CHROM") > -1)
	{
		$start_storing = "true";
		
		&print_message("Reading VCF file 2: $vcf_file_2","message");

		############################################################
		# Parsing of #CHROM data line to get a list of input files #
		############################################################
		@chrom_line_array = split(/\s+/,$single_line);
		$chrom_line_array_size = (scalar @chrom_line_array) - 1;
		$no_of_samples_chrom_line = $chrom_line_array_size - 8; # ignore first columns 0-9
		
		for ($array_count = 9; $array_count <= $chrom_line_array_size; $array_count++)
		{
			# Sample names go into sample_name_array_1
			$sample_name_array_2[$array_count - 8] = $chrom_line_array[$array_count];
			$file_count = $array_count - 8;
		}

	} # end of if CHROM line

} # reading VCF file 2


$total_variants_in_file_2 = $variants_in_file_2_count;
&print_both("Finished reading VCF file 2: $vcf_file_2\tNumber of variants: $total_variants_in_file_2\tNumber of samples: $no_of_samples_chrom_line\n");

##########################
# Get size of hash array #
##########################
$pos_2_hash_size = scalar keys %pos_2_hash;


#########################################################
# Show the list of sample names heading the columns in  #
# the VCF file to check they are the same for each file #
#########################################################
&print_message("Lists of sample names in the two files","message");
print COMMAND_LOG "Lists of sample names in the two files: \n\n";
$sample_name_mismatch = "false";
for ($sample_count = 1; $sample_count <=$no_of_samples_chrom_line; $sample_count++)
{
	&print_both("$sample_name_array_1[$sample_count]\t$sample_name_array_2[$sample_count]\n");
	if ($sample_name_array_1[$sample_count] ne $sample_name_array_2[$sample_count]){$sample_name_mismatch = "true"}
}


if ($sample_name_mismatch eq "true")
{
	&print_message("SAMPLE NAME MISMATCH","warning");

	&print_both("Sample names do not match.  Do you want to continue anyway? (y/n)  ");

	$answer=<STDIN>;
	if (substr(lc($answer),0,1) eq "n"){exit;}
}



######################################################################
# Now we want to run down file 1 (pos_1_hash) and list any variants  #
# which are only in file 1, i.e. not in pos_2_hash                   #
######################################################################
$variants_in_file_1_count = 0;
$only_in_file_1_count = 0;
$file_1_and_2_count = 0;

while ( ($chr_and_pos, $ALT_base_overall_1) = each(%pos_1_hash) ) 
{
    #print "Key: $position => $ALT_base_overall_1\n";
    $variants_in_file_1_count  = $variants_in_file_1_count  + 1;

	if (defined $pos_2_hash{$chr_and_pos})
    {
    	$file_1_and_2_count = $file_1_and_2_count + 1;
	}
    if (not defined $pos_2_hash{$chr_and_pos})
    {
    	$only_in_file_1_count = $only_in_file_1_count + 1;

    	#################
    	# Write to file #
    	#################
    	print O1VF "$chromosome\t$position\t$chr_and_pos\t$ALT_base_overall_1\n";
    }

    # Check for discordants here too? Compare to hash of 2 vs 1 discordants?

} # checking again through VCF 1


########################################
# Record this number as a second check #
########################################
$total_variants_in_file_1_check = $variants_in_file_1_count;


$file_1_only_calc = $variants_in_file_1_count - $concordant_position_count -$discordant_position_count;


############################################
# Need to add 1 to this as one variant = 0 #
############################################
$max_variants_at_same_pos = $max_variants_at_same_pos + 1;

close VCF2;
close DVF;
close CVF;
close O1VF;
close O2VF;
close DVD;
close DEBUG;

#########################
# Calculate some values #
#########################
$percent_discordant = ($discordant_neither_null_vcount / $concordant_neither_null_vcount) * 100;

$percent_discordant = sprintf("%.1f", $percent_discordant);


####################################################
# Warn if these calculations do not seem to add up #
####################################################
if ($only_in_file_1_count != $file_1_only_calc)
{
	&print_message("This looks a bit odd","warning");

	print "Number of variants found only in file 1: \t$only_in_file_1_count\n\n";
	print "This should be the same as variants in file 1 ($total_variants_in_file_1) minus concordant calls ($concordant_position_count) minus discordant calls ($discordant_position_count)\n\n";
	print "Maybe you should check the files a bit\n\n";
}

&print_message("Finished comparing VCF files","message");

print color 'bold yellow';
&print_both("FILE DETAILS:\n\n");
print color 'reset';

&print_both("  VCF file 1:                	               \t$vcf_file_1\n");
&print_both("  VCF file 2:                	               \t$vcf_file_2\n\n\n");

print color 'bold yellow';
&print_both("VARIANT POSITIONS:\n\n");
print color 'reset';

&print_both("  Total variant positions in file 1:            \t$total_variants_in_file_1 \t(SNPs: $snps_in_file_1_count \tIndels: $indels_in_file_1_count)\n");
&print_both("  Total variant positions in file 2:            \t$total_variants_in_file_2 \t(SNPs: $snps_in_file_2_count \tIndels: $indels_in_file_2_count)\n\n");


if ($debug eq "true"){print "   DEBUG: Size of hash array 1: $pos_1_hash_size\tSize of hash array 2: $pos_2_hash_size\n\n";}
if ($debug eq "true"){print "   DEBUG: Total variants in file 1 check: $total_variants_in_file_1_check\n\n";}

if (($unknowns_in_file_1_count > 0) || ($unknowns_in_file_2_count > 0))
{
	&print_both("  Unknown variants in file 1: $unknowns_in_file_1_count\n");
	&print_both("  Unknown variants in file 2: $unknowns_in_file_2_count\n\n");
}

&print_both("  Variant positions in both files:              \t$file_1_and_2_count\n");
&print_both("  Variant positions only in file 1:             \t$only_in_file_1_count\n");
&print_both("  Variant positions only in file 2:             \t$only_in_file_2_count\n\n");

&print_both("  Variant positions in both files, concordant:  \t$concordant_position_count\n");
&print_both("  Variant positions in both files, discordant:  \t$discordant_position_count \t(SNPs: $discordant_snps_count \tIndels: $discordant_indels_count)\n\n\n");


print color 'bold yellow';
&print_both("GENOTYPE DETAILS: (each variant position has $no_of_samples_chrom_line genotypes from $no_of_samples_chrom_line samples)\n\n");
print color 'reset';

&print_both("  Concordant genotypes (neither null):  \t$concordant_neither_null_vcount\n");
&print_both("  Concordant genotypes (both null):     \t$concordant_null_vcount\n\n");

&print_both("  Discordant genotypes (neither null):  \t$discordant_neither_null_vcount\n");
&print_both("  Discordant genotypes (one null):      \t$discordant_one_null_vcount\n\n");

&print_both("  Percentage discordant (neither null): \t$percent_discordant\n\n");


&print_message("OUTPUT FILES","message");

&print_both("  Concordant variants file:        \t $concordant_variants_file\n");
&print_both("  Discordant variants file:        \t $discordant_variants_file\n\n");

&print_both("  Only in 1 variants file:         \t $only_in_1_variants_file\n");
&print_both("  Only in 2 variants file:         \t $only_in_2_variants_file\n\n");

&print_both("  Details of discordant variants:  \t $discordant_variants_details\n\n\n");

&print_both( "Maximum number of variants at a single base position:  \t$max_variants_at_same_pos\n\n");

&print_both("Command log file: $command_log\n");

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