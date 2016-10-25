#!/usr/local/bin/perl

################################################################
#     plink2fastmixmod v9                                      #
#                                                              #
#     Converts PLINK .ped and .map files into the required     #
#     format for the Fast Moxed Model program                  #
#                                                              #
#     (Transforms from one row per sample to one row per SNP)  #
################################################################

#####################################################################
# you will need to put your Genotype data into a matrix with columns 
# representing markers and rows individuals. The SNP data should be 
# coded {0,1,2} (or if you prefer {-1,0,1}) with the intermediate 
# integer representing the heterozygote genotype. 

# You will also need to put your phenotype data into a vector 
# with entries representing individuals.
#####################################################################

##########################################################
# This version takes missing data (ie zero) into account #
##########################################################

# Mike Boursnell Feb 2009
# Animal Health Trust
# Newmarket
# UK

use strict;
use Getopt::Long;
use Term::ANSIColor;


my $version					= "10";

#####################################
# Define variables                  #
#####################################
my @item         		= ();
my @newline_affected      	= ();
my @newline_normal      	= ();
my @id_array			= ();
my @status_array		= ();
my @SNP_name_array		= ();
my @chromosome_array		= ();
my @position_array		= ();
my @row1 			= ();
my @row2 			= ();
my @row3 			= ();
my @row4 			= ();
my @standard_allele_A		= ();
my @standard_allele_B		= ();
my @mean_genotype_array		= ();

my @array_2d 			= (\@row1, \@row2);
my @array_2d_recoded		= (\@row3, \@row4);

my $prefix       		= "";
my $analysis_suffix		= "";
my $mapfile			= "";
my $pedfile			= "";
my $outfile			= "";
my $phenofile			= "";
my $mapfile_temp		= "";
my $pedfile_temp		= "";
my $genotype	 		= "";
my $genotype_string		= "";
my $genotype_bases		= "";
my $answer       		= "";
my $single_line			= "";
my $time_string			= "";
my $date_string			= "";
my $sample_name			= "";
my $pedigree			= "";
my $id				= "";
my $father			= "";
my $mother			= "";
my $gender			= "";
my $status			= "";
my $SNP_name			= "";
my $chromosome			= "";
my $position			= "";
my $mystring			= "";
my $chomped_all_string		= "";
my $allele_1			= "";
my $allele_2			= "";
my $allele_A			= "";
my $allele_B			= "";
my $genotype			= "";
my $ans				= "";
my $genotype_recoded		= "";
my $command			= "";
my $species_code		= "";
my $species			= "dog";
my $pheno_string		= "";
my $parameters_ok		= "false";

my $line_count    		= 1;
my $colcount     		= 0;
my $showcount    		= 0;
my $checkcount   		= 0;
my $snp_count			= 0;
my $id_count			= 0;
my $total_lines  		= 0;
my $total_items  		= 0;
my $total_items_on_first_row	= 0;
my $total_SNPs_on_first_row 	= 0;
my $total_SNPs			= 0;
my $year			= 0;
my $month			= 0;
my $day				= 0;
my $hour			= 0;
my $min				= 0;
my $sec				= 0;
my $correct_genotype_count 	= 0;
my $incorrect_genotype_count 	= 0;
my $total_genotypes_read	= 0;
my $help         		= 0;
my $snp_genotype_sum		= 0;
my $snp_valid_genotype_count	= 0;
my $snp_any_genotype_count	= 0;
my $snp_genotype_mean		= 0;
my $missing_data_count		= 0;
my $snp_missing_count		= 0;
my $snp_total_in_mapfile	= 0;
my $item_count			= 0;
my $means_used_count		= 0;
my $headers			= 1; # 1=add row headers, 0=don't add row headers
my $maf_value			= 0.05;
my $geno_value			= 0.1;
my $mind_value			= 0.1;




############################
# process -f flag (if any) #
############################

GetOptions("f=s"=>\$prefix);

print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      plink2fastmixmod       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";

print "  - This program converts PLINK .ped and .map files into  \n";
print "    the format required for the program Fast Mixed Model.\n\n";

print color 'bold yellow';

print "  NOTE: 	you MUST run Fast Mixed Model in R immediately  \n";
print "    		after running this program as it uses temporary files)\n\n";


print color 'reset';


#############################
# Get the input file name   #
#############################

if ($prefix eq "")
{
	print "Prefix for PLINK .map and .ped files:  ";
	$prefix = <STDIN>;
	chomp $prefix;
	print"\n";
}


##################################
# Make up map and ped file names #
##################################

$mapfile = $prefix.".map";
$pedfile = $prefix.".ped";



########################
# Check if files exist #
########################
if (! -f $pedfile)
{ 
	print "The PED file $pedfile does not exist\n\n";
	exit;
} 
if (! -e $mapfile)
{ 
	print "The MAP file $mapfile does not exist\n\n";
	exit;
} 


################################
# Get suffix for this analysis #
################################

print "Suffix for this analysis (e.g. a, b, c...):      ";
$analysis_suffix = <STDIN>;
chomp $analysis_suffix;
print "\n";


######################
# Make up file names #
######################

$outfile = $prefix."_".$analysis_suffix."_fmm.txt";
$phenofile = $prefix."_".$analysis_suffix."_pheno.txt";

$mapfile_temp = $prefix."_".$analysis_suffix.".map";
$pedfile_temp = $prefix."_".$analysis_suffix.".ped";



#############################
# Ask if it is dog or horse #
#############################

print color 'yellow';
print "####################################\n";
print "# Which species are you analysing? #\n";
print "####################################\n\n";
print color 'reset';
print "   <1> Dog\n";
print "   <2> Horse\n\n";

$species_code = <STDIN>;
chomp $species_code;

if ($species_code == 1){$species = "dog"}
if ($species_code == 2){$species = "horse"}



#############################
# Ask about parameters      #
#############################
 while ($parameters_ok eq "false")
{
	print "\n\n";
	print color 'yellow';
	print "####################\n";
	print "# PLINK parameters #\n";
	print "####################\n";
	print color 'reset';

	print "\nThese are the default parameters in use:    \n\n";
	print "    --maf $maf_value\n";
	print "    --geno $geno_value\n";
	print "    --mind $mind_value\n";

	print "    extra command line options (e.g. phenotype): $pheno_string\n";
	print "\nWould you like to change these? (y/n)    ";

	$ans = <STDIN>;
	chomp $ans;

	if ($ans eq "y" || $ans eq "Y")
	{
		print "\n";
		$parameters_ok = "false";
	
		print "    --maf value (press return to leave at $maf_value):   ";
		$ans = <STDIN>;
		chomp $ans;
		if ($ans ne ""){$maf_value = $ans;}

		print "    --geno value (press return to leave at $geno_value):   ";
		$ans = <STDIN>;
		chomp $ans;
		if ($ans ne ""){$geno_value = $ans;}

		print "    --mind value (press return to leave at $mind_value):   ";
		$ans = <STDIN>;
		chomp $ans;
		if ($ans ne ""){$mind_value = $ans;}


		print "    Extra command line options (press return to leave as $pheno_string):    ";
		$ans = <STDIN>;
		chomp $ans;
		if ($ans ne ""){$pheno_string = $ans;}

	}
	else
	{
		$parameters_ok = "true";
	}

} # end of ask about parameters


##############################################################
# Use PLINK to re-order SNPs in genome order                 #
# (This permanently re-orders the original map and ped files #
# which is a good thing. Doesn't effect PLINK analysis)      #
##############################################################

$command = "plink --noweb --".$species." --allow-no-sex --file $prefix --recode --out $prefix";
print("$command\n");
system("$command");



#########################################################################
# Use PLINK to create new map and ped files with the filtering in place #
#########################################################################

$command = "plink --noweb --".$species." --allow-no-sex --make-founders --nonfounders --file $prefix $pheno_string --maf $maf_value --geno $geno_value --mind $mind_value --recode --out $prefix"."_".$analysis_suffix;

print("$command\n");
system("$command");


############################
# Open files for output    #
############################

open (OUT, ">$outfile") || die "Cannot open $outfile";
open (PHE, ">$phenofile") || die "Cannot open $phenofile";


############################
# Open the files for INPUT #
############################

open (PED, "$pedfile_temp") || die "Cannot open $pedfile_temp";
open (MAP, "$mapfile_temp") || die "Cannot open $mapfile_temp";


################################
# Read data in MAP file to get # 
# the list of SNP names        #
################################

$snp_count = 0;
&print_message("Reading MAP file...$mapfile_temp");

while ($single_line = <MAP>)
{

	# Remove the end of line character #
	chomp $single_line;
	&chomp_all ($single_line);

	@item=split(/\t/,$single_line);

	$snp_count = $snp_count + 1;

	$chromosome = $item[0];
	$SNP_name = $item[1];
	$position = $item[3];

	$chromosome_array[$snp_count]=$chromosome;
	$SNP_name_array[$snp_count]=$SNP_name;
	$position_array[$snp_count]=$position;


}

close MAP;

$snp_total_in_mapfile = $snp_count;

print " Total SNPs in map file:  $snp_total_in_mapfile\n\n\n";

############################
# Read data in PED file    #
############################

&print_message("Reading PED file...$pedfile_temp");


###########################################
# Work down the lines in the ped file     #
###########################################
$line_count=1;

while ($single_line = <PED>) {

	# Remove the end of line character #
	chomp $single_line;
	&chomp_all ($single_line);


	#########################################################
        # Split line at TABS OR SPACES into the array 'item'    #
	#########################################################
        @item=split(/\s+/,$single_line);

	$pedigree = $item[0];
	$id = $item[1];
	$father = $item[2];
	$mother = $item[3];
	$gender = $item[4];
	$status = $item[5];


	#####################################################
	# Check how many items there are on the row.        #
	#####################################################
	$item_count = 0;
	while (length $item[$item_count]>0)
	{
		$item_count = $item_count + 1;
	}
	$total_items=$item_count;
	$total_SNPs = ($total_items - 6)/2;

	print "Line $line_count.  \tTotal SNPs on row = $total_SNPs\n";

	$command = "tput cuu1";
        system("$command");

	#############################################
	# Read the genotype columns  into array_2d  #
	#############################################
	for ($snp_count = 1;$snp_count <= $total_SNPs;$snp_count ++)
	{
		$genotype = $item[($snp_count * 2) + 4]." ".$item[($snp_count * 2) + 5];

		$array_2d[$line_count][$snp_count] = $genotype;
	}


	##################################
	# Assign each id to the id array #
	##################################
	#chomp $id;
	$id_array[$line_count] = $id;


	##########################################
	# Assign each status to the status array #
	##########################################
	#chomp $status;
	$status_array[$line_count] = $status;





	######################################################
	# Store the number of items on the first line        #
	# so we can check that all other lines have          #
	# the same number of items                           #
	#                                                    #
	######################################################

	if ($line_count==1) 
	{
		$item_count = 0;
		while (length $item[$item_count]>0)
		{
			$item_count = $item_count + 1;
		}

		$total_items_on_first_row = $item_count;

		$total_SNPs_on_first_row = ($total_items_on_first_row - 6) / 2;

		print " Number of SNPs on the first row is $total_SNPs_on_first_row\n\n";

	}



	######################################################
	# For each subsequent line, warn if the number of    #
	# items on the row doesn't match the number on the   #
	# first line.                                        #
	######################################################
	if ($item_count!=$total_items_on_first_row)
	{
		print "\n";
		print "##############\n";
		print "# WARNING!!! #\n";
		print "##############\n\n";
		print "Number of items on line $line_count is $item_count\n";
		print "which is not the same as the $total_items_on_first_row items on the first row.\n\n";
		print "Press return to continue\n";
		$answer= <STDIN>;
		print "\n";
	}

	######################################################
	# Move on to next line                               #
	######################################################
	$line_count = $line_count + 1;

} # end of while loop

$total_lines=$line_count -1;

print " Total number of data lines in the ped file: $total_lines\n\n\n";

&print_message("Recording standard alleles for each SNP");


#######################################################
# Find which alleles are present for each of the SNPs #
#######################################################

for ($snp_count = 1;$snp_count <= $total_SNPs;$snp_count ++)
{

	$allele_A="";
	$allele_B="";

	#########################################################################
	# Note: If data is missing (0 0) then allele_A will be assigned as zero #
	#########################################################################
	for ($line_count=1; $line_count <= $total_lines;$line_count ++)
	{
		$genotype = $array_2d[$line_count][$snp_count];

		$allele_1 = substr($genotype,0,1);
		$allele_2 = substr($genotype,2,1);

		#Consider allele_1
		if ($allele_1 ne "0")
		{
			#Add first allele
			if (($allele_A eq "") && ($allele_B eq "")) {$allele_A = $allele_1}

			#Add second allele (if different)
			if (($allele_A ne "") && ($allele_B eq "") && ($allele_1 ne $allele_A)) {$allele_B = $allele_1}
		}

		if ($allele_2 ne "0")
		{
			#Consider allele_2
			#Add first allele
			if (($allele_A eq "") && ($allele_B eq "")) {$allele_A = $allele_2}

			#Add second allele (if different)
			if (($allele_A ne "") && ($allele_B eq "") && ($allele_2 ne $allele_A)) {$allele_B = $allele_2}

		}

		#Exit loop if both filled
		if (($allele_A ne "") && ($allele_B ne "")) {last;}

	}

	################################################################
	# Now store the standard alleles in the array standard_alleles #
	################################################################

	$standard_allele_A[$snp_count] = $allele_A;
	$standard_allele_B[$snp_count] = $allele_B;

} # end of finding which are standard alleles




#####################################
# Write output file column headings #
#####################################


if ($headers == 1)
{
	for ($snp_count = 1;$snp_count <= $total_SNPs;$snp_count ++)
	{
		print OUT "$SNP_name_array[$snp_count]";
		if ($snp_count < $total_SNPs){print OUT " "}
	}
	print OUT "\n";
}

#######################################################
# Work out the recoded genotype for each SNP          #
# (and the mean genotype)                             #
# (This uses the data stored in the 2d array)         #
#######################################################

print "Working out the recoded and mean genotype for each SNP...\n\n";

for ($snp_count = 1;$snp_count <= $total_SNPs;$snp_count ++)
{

	#Set to zero for each SNP (used for calculating mean genotype)
	$snp_any_genotype_count = 0;
	$snp_valid_genotype_count = 0;
	$snp_genotype_sum = 0;
	$snp_missing_count = 0;

	#print "SNP $snp_count/$total_SNPs\n";

	if (int($snp_count/1000) == $snp_count/1000)
	{
		print "SNP $snp_count/$total_SNPs\n";
		$command = "tput cuu1";
        	system("$command");
	}

	

	# Work down the samples
	for ($line_count=1; $line_count <= $total_lines;$line_count ++)
	{
		$genotype = $array_2d[$line_count][$snp_count];
		$snp_any_genotype_count = $snp_any_genotype_count + 1;

		#Count if this genotype is missing
		if ($genotype eq "0 0"){$snp_missing_count = $snp_missing_count + 1}


		#set default recoded genotype as "X"
		$array_2d_recoded[$line_count][$snp_count] = "X";
		$genotype_recoded = "Y";

		$total_genotypes_read = $total_genotypes_read + 1;

		$allele_1 = substr($genotype,0,1);
		$allele_2 = substr($genotype,2,1);
		$allele_A = $standard_allele_A[$snp_count];
		$allele_B = $standard_allele_B[$snp_count];

		####################################################
		# Only recode and use for mean if not missing data #
		####################################################
		if ($allele_1 ne "0" && $allele_2 ne "0")
		{

			#Recoding
			if (($allele_1 eq $allele_A) && ($allele_2 eq $allele_A)) {$genotype_recoded = "0"}

			if (($allele_1 eq $allele_B) && ($allele_2 eq $allele_B)) {$genotype_recoded = "2"}

			if (($allele_1 eq $allele_A) && ($allele_2 eq $allele_B)) {$genotype_recoded = "1"}

			if (($allele_1 eq $allele_B) && ($allele_2 eq $allele_A)) {$genotype_recoded = "1"}

			###############################################
			# Record if genotype recoded successfully     #
			# Also add to the array for recoded genotypes #
			###############################################
			if ($genotype_recoded eq "0" || $genotype_recoded eq "1" || $genotype_recoded eq "2")
			{
				$correct_genotype_count = $correct_genotype_count + 1;
				$snp_genotype_sum = $snp_genotype_sum + $genotype_recoded;
				$snp_valid_genotype_count = $snp_valid_genotype_count + 1;

				$array_2d_recoded[$line_count][$snp_count] = $genotype_recoded;

			}

			#Record if genotype recoded unsuccessfully
			if ($genotype_recoded ne "0" && $genotype_recoded ne "1" && $genotype_recoded ne "2")
			{
				$incorrect_genotype_count = $incorrect_genotype_count + 1;
			}

			#print "Line: $line_count  SNP: $snp_count\n";
			#print "Alleles found : $allele_1 $allele_2   Standard alleles: $allele_A $allele_B\n";
			#print "Genotype recoded: $genotype_recoded\n\n";
			#$ans = <STDIN>;
		}


		if($array_2d_recoded[$line_count][$snp_count] eq "X")
		{
			$missing_data_count = $missing_data_count + 1;
		}

		
	}

	################################################
	# Calculate 'mean genotype' and store in array #
	################################################
	
	$snp_genotype_mean =0;
	if ($snp_valid_genotype_count > 0)
	{
		$snp_genotype_mean = $snp_genotype_sum / $snp_valid_genotype_count;
		$snp_genotype_mean = int ($snp_genotype_mean * 100)/100;
	}
	$mean_genotype_array[$snp_count] = $snp_genotype_mean;


}



########################################
# Now write this to the output file    #
# (adding mean genotype if none found) #
########################################

&print_message("Writing output file");

for ($line_count=1; $line_count <= $total_lines;$line_count ++)
{

	#Write ID to output file
	if ($headers == 1) {print OUT "$id_array[$line_count] "}

	for ($snp_count = 1;$snp_count <= $total_SNPs;$snp_count ++)
	{

		$genotype_recoded = $array_2d_recoded[$line_count][$snp_count];

		$allele_1 = substr($genotype,0,1);
		$allele_2 = substr($genotype,2,1);
		$allele_A = $standard_allele_A[$snp_count];
		$allele_B = $standard_allele_B[$snp_count];

		###################################################################
		# If genotype shows that data was missing, set to 'mean genotype' #
		###################################################################
		if ($genotype_recoded eq "X")
		{
			$genotype_recoded = $mean_genotype_array[$snp_count];
			$means_used_count = $means_used_count + 1;
		}

		#######################################
		# Add recoded genotype to output file #
		#######################################

		#print OUT "G: $genotype($allele_A$allele_B)$genotype_recoded, ";
		print OUT "$genotype_recoded";
		if ($snp_count < $total_SNPs){print OUT " "}
	}

	#Write EOL to output file
	print OUT "\n";

	##############################################
	# Write the status to the phenotype OUT file #
	##############################################

	print PHE "$status_array[$line_count]\n";

}


#######################################################
# Show how many genotypes have been read successfully #
#######################################################

print color 'bold cyan';
print "\n######################\n";
print "# Recoding completed #\n";
print "######################\n\n";
print "Total no of SNPs:		  \t$total_SNPs\n";
print "Total no of samples:		  \t$total_lines\n\n";

print "Total number of genotypes read:     \t$total_genotypes_read\n\n";
print "No of missing genotypes:            \t$missing_data_count\n\n";
print "Genotypes recoded successfully:     \t$correct_genotype_count\n";
print "Genotypes not recoded successfully: \t$incorrect_genotype_count\n";
print "Genotypes where mean was used:	   \t$means_used_count\n\n";

print "Output file (FMM genotypes):	   \t$outfile\n";
print "Output file (phenotypes):	   \t$phenofile\n";
print "Output file (map)		   \t$mapfile_temp\n\n\n";


print "############################################################\n";
print "#      To complete Fast Mixed Model analysis using R,      #\n";
print "#      read the file 'Fast Mixed Model Instructions'.      #\n";
print "#                                                          #\n";
print "#       Type source ('/home/genetics/scripts/FMM.R')       #\n";
print "############################################################\n\n";

print color 'reset';

###########################
#     Close all files     #
###########################
close PED;
close OUT;
close PHE;


exit;

##########################################################################################


##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

foreach (@_) {s/\n//g}  
foreach (@_) {s/\r//g}  

}

######################################
# Subroutine to print screen message #
######################################

sub print_message
{
	my $message 		= "";
	my $message_length 	= "";
	my $pos_count		= 0;
	
	$message = $_[0];
	$message_length = length($message);
	
	print "\n\n";
	print color 'yellow';
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++)
	{
		print "#";
	}
	print "\n#    $message    #\n";
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++)
	{
		print "#";
	}
	print "\n\n";
	
	print color 'reset';

}


