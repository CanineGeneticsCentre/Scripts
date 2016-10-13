#!/usr/local/bin/perl 

################################################################
#     check_duplicates                                         #
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

#####################################
# Define variables                  #
#####################################
my @item         		= ();
my @id_array			= ();
my @status_array		= ();
my @SNP_name_array		= ();
my @chromosome_array		= ();
my @position_array		= ();
my @row1 			= ();
my @row2 			= ();
my @row3 			= ();
my @row4 			= ();
my @percentage_array		= ();

my @array_2d 			= (\@row1, \@row2);


my $infile       		= "";
my $mapfile			= "";
my $pedfile			= "";
my $outfile			= "";
my $phenofile			= "";
my $genotype	 		= "";
my $answer       		= "";
my $single_line			= "";
my $pedigree			= "";
my $id				= "";
my $father			= "";
my $mother			= "";
my $gender			= "";
my $status			= "";
my $SNP_name			= "";
my $chromosome			= "";
my $position			= "";
my $genotype			= "";
my $ans				= "";
my $sample_name_1		= "";
my $sample_name_2		= "";
my $genotype_1			= "";
my $genotype_2			= "";
my $command			= "";
my $some_over_threshold 	= "False";


my $line_count    		= 1;
my $item_count    		= 0;
my $snp_count			= 0;
my $total_lines  		= 0;
my $total_items  		= 0;
my $total_items_on_first_row	= 0;
my $total_SNPs_on_first_row 	= 0;
my $total_SNPs			= 0;
my $snp_total_in_mapfile	= 0;
my $sample_1			= 0;
my $sample_2			= 0;
my $genotype_match_count 	= 0;
my $genotype_comparison_count	= 0;
my $percentage_match		= 0;
my $threshold			= 0;
my $no_of_comparisons 		= 0;
my $total_percentage_match	= 0;
my $average_percentage		= 0;
my $x				= 0;
my $y				= 0;
my $standard_deviation		= 0;
my $count			= 0;
my $sum_x			= 0;


############################
# process -f flag (if any) #
############################

GetOptions("f=s"=>\$infile);

print "\n";
print "################################################################\n";
print "#     check_duplicates v2                                      #\n";
print "#                                                              #\n";
print "#     Checks through PLINK .ped files to see if any of the     #\n";
print "#     samples have identical genotypes.                        #\n";
print "#                                                              #\n";
print "################################################################\n\n";



#############################
# Get the input file name   #
#############################

if ($infile eq "")
{
	print "Prefix for PLINK .map and .ped files:  ";
	$infile = <STDIN>;
	chomp $infile;
	print"\n";
}

$mapfile = $infile.".map";
$pedfile = $infile.".ped";
$outfile = $infile."_duplicates.txt";


############################
# Open files for output    #
############################

open (OUT, ">$outfile") || die "Cannot open $outfile";


############################
# Open the files for INPUT #
############################

open (PED, "$pedfile") || die "Cannot open $pedfile";
open (MAP, "$mapfile") || die "Cannot open $mapfile";

#################
# Get threshold #
#################
print "What is threshold percentage match above which you want to flag (default = 98) ";

$threshold = <STDIN>;
chomp $threshold;

if ($threshold eq "") 
{
	$threshold = 98;
	print "Threshold set to default (98%)\n";
}

######################
# Record in LOG file #
######################

print OUT "Checking for duplicates in genotyped samples\n\n";
print OUT "Ped file: \t$pedfile\n";
print OUT "Map file: \t$mapfile\n";
print OUT "Threshold: \t$threshold\n\n";


################################
# Read data in MAP file to get # 
# the list of SNP names        #
################################

$snp_count = 0;
print "Reading MAP file...\n\n";

while ($single_line = <MAP>)
{

	$single_line = &chomp_all ($single_line);

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


print "$snp_total_in_mapfile SNPs in the MAP file\n";
print "\n\n";


############################
# Read data in PED file    #
############################

print "Reading PED file...\n\n";

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

		$total_items_on_first_row=$item_count;
		$total_SNPs_on_first_row = ($total_items_on_first_row - 6) / 2;

		print "Number of SNPs on the first row is $total_SNPs_on_first_row\n\n";

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

print "\nTotal no of SNPs:		  \t$total_SNPs\n";
print "Total no of samples:		  \t$total_lines\n\n";



################################################
#        Check for duplicate genotypes         #
################################################

print "Checking for duplicates...\n\n";
print "Sample 1 \tSample 2 \tPercentage match\n";
print "======== \t======== \t================\n";

print OUT "Sample 1 \tSample 2 \tPercentage match\n";
print OUT "======== \t======== \t================\n";

$total_percentage_match = 0;
$no_of_comparisons = 0;

for ($sample_1=1;$sample_1 <=$total_lines;$sample_1++)
{

	print " Checking sample $id_array[$sample_1] \t($sample_1/$total_lines)\n";

	$command = "tput cuu1";
        system("$command");
	

	for ($sample_2=$sample_1 + 1;$sample_2<=$total_lines;$sample_2++)
	{

		$sample_name_1 = $id_array[$sample_1];
		$sample_name_2 = $id_array[$sample_2];

		$genotype_match_count 	= 0;
		$genotype_comparison_count = 0;

		for ($snp_count = 1;$snp_count<=$total_SNPs;$snp_count++)
		{

			$SNP_name = $SNP_name_array[$snp_count];

			$genotype_1=$array_2d[$sample_1][$snp_count];
			$genotype_2=$array_2d[$sample_2][$snp_count];

			$genotype_comparison_count = $genotype_comparison_count + 1;

			if ($genotype_1 eq $genotype_2)
			{
				$genotype_match_count = $genotype_match_count + 1;
			}
		}

		#print "Samples: $sample_name_1 vs $sample_name_2   Matches: $genotype_match_count/$genotype_comparison_count\n";

		$percentage_match = 0;

		if ($genotype_comparison_count > 0 )
		{
			$percentage_match = ($genotype_match_count / $genotype_comparison_count) * 100;
			$no_of_comparisons = $no_of_comparisons + 1;
			$percentage_array[$no_of_comparisons] = $percentage_match;
		}

		$total_percentage_match = $total_percentage_match + $percentage_match;

		if ($percentage_match >= $threshold)
		{
			
			print "                                                                                   \n";

			$command = "tput cuu1";
        		system("$command");

			print     "$sample_name_1    \t$sample_name_2    \t$percentage_match\n";
			print OUT "$sample_name_1    \t$sample_name_2    \t$percentage_match\n";
			$some_over_threshold = "True";		

		}


		

		
	}



}

################################
# Calculate average percentage #
################################

$average_percentage = $total_percentage_match/$no_of_comparisons;

$average_percentage = int($average_percentage * 100)/100;

################################
# Calculate standard deviation #
################################

$sum_x = 0;
for ($count = 1;$count <=$no_of_comparisons; $count ++)
{
	$x = ($percentage_array[$count] - $average_percentage)**2;
	$sum_x = $sum_x + $x;
}

#########################
# Divide sum_x by (n-1) #
# then take square root #
#########################
$y = $sum_x/($no_of_comparisons - 1);
$standard_deviation = sqrt ($y);
$standard_deviation = int ($standard_deviation * 100)/100;

#######################################################
# Show how many genotypes have been read successfully #
#######################################################

print "\n\nTotal no of SNPs:		  \t$total_SNPs\n";
print "Total no of samples:		  \t$total_lines\n\n";
print "Average percentage match:	  \t$average_percentage\n";
print "Standard deviaton:		  \t$standard_deviation\n\n";

print OUT "\n\nTotal no of SNPs:		  \t$total_SNPs\n";
print OUT "Total no of samples:		  \t$total_lines\n\n";
print OUT "Average percentage match:	  \t$average_percentage\n";
print OUT "Standard deviaton:		  \t$standard_deviation\n\n";


if ($some_over_threshold eq "True")
{
print "Output file:\t$outfile\n\n";
}

if ($some_over_threshold eq "False")
{
print "None of the genotype comparisons were over the threshold of $threshold%\n\n";
print OUT "None of the genotype comparisons were over the threshold of $threshold%\n\n";
}



###########################
#     Close all files     #
###########################
close PED;
close OUT;


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

