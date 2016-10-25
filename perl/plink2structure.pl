#!/usr/bin/perl

################################################################
#     plink2structure v8                                       #
#                                                              #
#     Converts files from PLINk format (map and ped files)     #
#     to a format which can be imported into STRUCTURE         #
#                                                              #
################################################################


# Mike BoursnellMar 2009

use strict;
use Getopt::Long;


#####################################
# Define variables                  #
#####################################
my @item_ped        		= ();
my @item_temp			= ();
my @piece			= ();
my @names     			= ();
my @newnames			= ();
my @id_array			= ();
my @snp_array			= ();
my @map_array			= ();
my @item_snp			= ();
my @item_aff			= ();
my @SNP_position_array		= ();
my @sample_name_array		= ();
my @status_array		= ();
my @snp_call_rate_array		= ();
my @snp_monomorphic_array	= ();
my @row1 			= ();
my @row2 			= ();
my @array_2d 			= (\@row1, \@row2);

my $pedfile       		= "";
my $outfile      		= "";
my $logfile			= "";
my $mapfile			= "";
my $mutfile			= "";
my $listfile			= "";
my $affectfile			= "";
my $single_line			= "";
my $new_single_line 		= "";
my $list_ID			= "";
my $sanger_ID			= "";
my $date_string			= "";
my $SNP_name			= "";
my $allele_1			= "";
my $allele_2			= "";
my $new_ID			= "";
my $ans				= "";
my $pedigree			= "";
my $pedigree_new		= "";
my $ID				= "";
my $father			= "";
my $mother			= "";
my $gender			= "";
my $status			= "";
my $log				= "";
my $snp				= "";
my $snp_name			= "";
my $sample_name			= "";
my $status			= "";
my $name 			= "";
my $new_status			= "";
my $genotype			= "";
my $new_status_string		= "";
my $prefix			= "";
my $mystring			= "";
my $first_genotype		= "";
my $check_char			= "";
my $check_ord			= "";
my $separator			= "";

my $help         		= 0;
my $linecount    		= 0;
my $colcount     		= 0;
my $count			= 0;
my $showcount    		= 0;
my $item_count    		= 0;
my $snp_count			= 0;
my $sample_count		= 0;
my $total_items   		= 0;
my $total_items_row_1		= 0;
my $total_items_temp		= 0;
my $renamed_count		= 0;
my $total_lines  		= 0;
my $total_samples		= 0;
my $total_SNPs			= 0;
my $total_in_list		= 0;
my $lines_counted		= 0;
my $found_on_list		= 0;
my $names_not_used_count	= 0;
my $copy_if_no_renaming		= 0;
my $mut_total			= 0;
my $snp_total			= 0;
my $aff_total			= 0;
my $mut_count			= 0;
my $store_mut_count		= 0;
my $position			= 0;
my $mutation			= 0;
my $SNP_position		= 0;
my $affected			= 0;
my $dot_pos			= 0;
my $snp_missing_data_count	= 0;
my $snp_missing_proportion	= 0;
my $snps_passed_count		= 0;
my $snps_failed_count		= 0;
my $snps_failed_mono_count	= 0;
my $minimum_call_rate		= 0.9;
my $remove_monomorphics		= 1;
my $snp_monomorphic_count	= 0;
my $pos				= 0;
my $found_tab			= 0;
my $mod				= 0;
my $no_gaps			= 0;
my $no_SNPs_output		= 0;

my $missing_data_value		= -9;


print"\n";
print"###################################################\n";
print"#                                                 #\n";
print"#     plink2structure v8                          #\n";
print"#                                                 #\n";
print"#     Converts ped file from PLINK to a format    #\n";
print"#     which can be used by STRUCTURE              #\n";
print"#                                                 #\n";
print"#     (also uses the PLINK map file)              #\n";
print"#                                                 #\n";
print"###################################################\n\n";



##########################################
# Get the file names and open the files  #
##########################################

print "Prefix of .map and .ped files:    ";
$prefix = <STDIN>;
chomp $prefix;

$pedfile = $prefix.".ped";
$mapfile = $prefix.".map";
$outfile = $prefix.".structure";

open (PED, "$pedfile") || die "Cannot open input file: $pedfile";
open (MAP, "$mapfile") || die "Cannot open input file: $mapfile";
open (OUT, ">$outfile")|| die "Cannot create output file: $!";

print "\nWhat is the minimum call rate for each SNP?\n";
print "(SNPs with a call rate lower than this are not included in the output file)\n\n";
print "Press return to leave at $minimum_call_rate:   ";

$ans = <STDIN>;
chomp $ans;
if ($ans ne "") {$minimum_call_rate = $ans}


###################################################
# Read in all SNPs in map file                    #
# and put them in snp_array                       #
# (This is so they can be used as column headers) #
###################################################
print "Reading MAP file...\n\n";
print "showing first 5 SNPs\n\n";

$count =0;
@map_array= <MAP>;
foreach $single_line (@map_array)
{
	$count = $count + 1;
	@item_snp=split(/\t/,$single_line);
	$snp_name= $item_snp[1];
	chomp $snp_name;

	if ($count <= 5)
	{
		print "$snp_name\n";
	}

	$snp_array[$count]=$snp_name;

}
close MAP;

$snp_total = scalar (@map_array);


print "...\n\n";
print "$snp_total SNPs in the MAP file\n";
print "\n\n";


##################################################
# Check whether PED file has tabs or just spaces #
##################################################

print "######################################################\n";
print "# Checking PED file to see if it uses tabs or spaces #\n";
print "######################################################\n\n";

$found_tab = 0;

while ($single_line = <PED>) {

	chomp $single_line;

	$linecount = $linecount + 1;

	for ($pos=0;$pos <=length ($single_line);$pos++)
	{
		$check_char= substr $single_line,$pos,1;
		$check_ord = ord($check_char);

		if ($check_ord == 9)
		{
			$found_tab = 1;
		}
	}


} # end of while loop

if ($found_tab == 1)
{
	print "     This looks like a ped file from Progeny with tabs separating the genotypes.\n";
	print "     I shall assume this is a Progeny PLINK file and proceed accordingly...\n";
	print "     >";
	$ans = <STDIN>;
}


if ($found_tab == 0)
{
	print "     This looks like a ped file from PLINK with spaces separating the genotypes.\n";
	print "     I shall assume this is a PLINK ped file and proceed accordingly...\n";
	print "     >";
	$ans = <STDIN>;
}

close PED;



#
# Re-open PED file #
#

open (PED, "$pedfile") || die "Cannot open input file: $pedfile";


###########################################
# Work down the lines in the PED file     #
###########################################
$linecount=0;
while ($single_line = <PED>) {

	# Remove the end of line character #
	chomp $single_line;
	&chomp_all ($single_line);

	$linecount = $linecount + 1;

	################################################
	# If there are no tabs, add them at this point #
	################################################

	if ($found_tab == 0)
	{
		#Slit line at spaces
		@item_temp=split(/ /,$single_line);
		$total_items_temp = scalar @item_temp;

		# Reconstruct the line using tabs between the genotypes
		$new_single_line = "";
		for ($item_count =0;$item_count <= $total_items_temp;$item_count++)
		{
			$separator = chr(9);
			if ($item_count == 0) {$separator = ""}
			if ($item_count > 0 && $item_count <=6) {$separator = chr(9)}

			#from pos=6 onwards the separator alternates between space and tab
			if ($item_count > 6)
			{
				$mod = ($item_count % 2);
				if ($mod == 0) {$separator = chr(9)}
				if ($mod == 1) {$separator = " "}
			}

			if ($item_count == $total_items_temp) {$separator = ""}

			$new_single_line = $new_single_line.$separator.$item_temp[$item_count];
			
		}

		# Replace the original line with the modified version with the tabs
		$single_line = $new_single_line;
	}

	#####################################################
        # Split line at tabs into the array 'item_ped'   #
	#####################################################
        @item_ped=split(/\t/,$single_line);

	$pedigree  = $item_ped[0];
	$ID = $item_ped[1];
	$father  = $item_ped[2];
	$mother  = $item_ped[3];
	$gender  = $item_ped[4];
	$status  = $item_ped[5];


	#################################
	# Store ID and status in arrays #
	#################################

	$id_array[$linecount] = $ID;
	$status_array[$linecount] = $status;

	#####################################################
	# Check how many items there are on the row.        #
	#####################################################
	$item_count = 0;

	while (length $item_ped[$item_count]>0)
	{
		$item_count = $item_count + 1;
	}

	$total_items = $item_count;
	$total_SNPs = $total_items - 6;

	print "Line: $linecount   \tItems: $item_count   \tSNPs: $total_SNPs\n";

	##################################
	# Store genotype in the 2D array #
	##################################

	for ($snp_count = 1; $snp_count <= $total_SNPs; $snp_count ++)
	{
		$genotype = $item_ped[$snp_count + 5];
		chomp $genotype;

		if ($genotype eq "0 0")
		{
			$genotype = $missing_data_value." ".$missing_data_value;
		}
		$array_2d[$linecount][$snp_count] = $genotype;
	
	}


	
	######################################################
	# Store the number of items on the first line        #
	# so we can check that all other lines have          #
	# the same number of items                           #
	######################################################
	if ($linecount==1) 
	{
		$total_items_row_1=$total_items;
		print "Number of items on the first row is $total_items_row_1.\n\n";
	}


	######################################################
	# Warn if the number of items on the row             #
	# doesn't match the number on the first line.        #
	######################################################
	if ($total_items!=$total_items_row_1)
	{
		print "##############\n";
		print "# WARNING!!! #\n";
		print "##############\n\n";
		print "Number of items on line $linecount is $item_count\n";
		print "which is not the same as the $total_items_row_1 items on the first row.\n\n";
		$ans= <STDIN>;
		print "\n";
	}


} # end of while loop

$total_samples = $linecount;
#$total_SNPs = $snp_count;

print "\n\n";

print "Number of samples in the PED file:  \t$total_samples\n";
print "Number of SNPs in the PED file:     \t$total_SNPs \n\n";


###################################################
# Work out frequency of missing data for each SNP #
###################################################

print "Calculating frequency of missing data for each SNP...\n";

$snps_passed_count = 0;

for ($snp_count = 1;$snp_count <= $total_SNPs; $snp_count ++)
{
	$snp_missing_data_count	= 0;
	$snp_missing_proportion = 0;

	$first_genotype = "";
	$snp_monomorphic_array[$snp_count] = 1;
	for ($sample_count = 1; $sample_count <=$total_samples;$sample_count ++)
	{
		$genotype = $array_2d[$sample_count][$snp_count];
		

		#Assign first real genotype
		if (($genotype ne "-9 -9") && ($first_genotype eq ""))
		{
			$first_genotype = $genotype;
		}

		# Check for non-monomorphic-ness
		

		if (($genotype ne $first_genotype) && ($genotype ne "-9 -9"))
		{
			$snp_monomorphic_array[$snp_count] = 0;
		}

		# Count missing genotypes
		if ($genotype eq "-9 -9")
		{
			$snp_missing_data_count = $snp_missing_data_count + 1;
		}
		
	}

	$snp_missing_proportion = $snp_missing_data_count / $total_samples;


	$snp_call_rate_array[$snp_count] = 1 - $snp_missing_proportion;

	if ($snp_monomorphic_array[$snp_count] == 1) {$snp_monomorphic_count = $snp_monomorphic_count + 1}

	if (($snp_call_rate_array[$snp_count] >= $minimum_call_rate) && ($snp_monomorphic_array[$snp_count] == 0))
	{
		$snps_passed_count = $snps_passed_count + 1;

	}

}


########################################################
# Now write the SNP names on the first row             #
########################################################

for ($snp_count = 1;$snp_count <= $snp_total;$snp_count ++)
{
	if (($snp_call_rate_array[$snp_count] >= $minimum_call_rate) && ($snp_monomorphic_array[$snp_count] == 0))
	{
		print OUT "$snp_array[$snp_count]";
		if ($snp_count < $total_SNPs) 
		{
			print OUT "\t";
		}

	}

	
}
print OUT "\n";


#############################################
# Print out the array-2d to the output file #
#############################################

for ($sample_count = 1; $sample_count <=$total_samples;$sample_count ++)
{
	print OUT "$id_array[$sample_count]\t";
	print OUT "$status_array[$sample_count]\t";

	$no_gaps = 2;
	$no_SNPs_output = 0;

	for ($snp_count = 1; $snp_count <=$total_SNPs;$snp_count ++)
	{
		if (($snp_call_rate_array[$snp_count] >= $minimum_call_rate) && ($snp_monomorphic_array[$snp_count] == 0))
		{
			print OUT "$array_2d[$sample_count][$snp_count]";

			$no_SNPs_output = $no_SNPs_output + 1;

			if ($snp_count < $total_SNPs)
			{
				print OUT "\t";
				$no_gaps = $no_gaps + 1;
			}
		}


		#Count SNPs which failed on call rate
		if (($snp_call_rate_array[$snp_count] < $minimum_call_rate) && ($sample_count ==1))
		{
			$snps_failed_count = $snps_failed_count + 1;
		}

		# Count SNPs which failed because they were monomorphic
		if (($snp_monomorphic_array[$snp_count] == 1) && ($sample_count ==1))
		{
			$snps_failed_mono_count = $snps_failed_mono_count + 1;
		}



	}
	print OUT "\n";
}



###########################
#     Close both files    #
###########################

close MAP;
close PED;
close OUT;

print "\n";
print "Input ped file:                    \t$pedfile\n";
print "Input map file:                    \t$mapfile\n";
print "Output ped file:                   \t$outfile\n\n";

print "Missing data character:            \t$missing_data_value\n\n";

print "Minimum call rate:                 \t$minimum_call_rate\n";
print "Number of monomorphic SNPs:        \t$snp_monomorphic_count\n";
print "Number of SNPs failed call rate:   \t$snps_failed_count\n\n";

print "Number of samples:                 \t$total_samples\n";
print "Number of SNPs remining:           \t$snps_passed_count\n\n";



exit;

##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

my $right_char		= "";
my $len			= 0;


$mystring = $_[0];


$len = length ($mystring);
$right_char = substr $mystring,1,$len;

#If CR or LF then remove them #
if (ord($right_char)==13 || ord($right_char)==10 || ord($right_char)==0)
{
	chop $mystring;
}

$len = length ($mystring);
$right_char = substr $mystring,1,$len;

#If CR or LF then remove them #
if (ord($right_char)==13 || ord($right_char)==10 || ord($right_char)==0)
{
	chop $mystring;
}

$mystring = $mystring;

}


##########################################################
# Subroutine to show what characters are in a string     #
##########################################################
sub show_chars
{

my $char		= "";
my $ord			= "";
my $len			= 0;
my $pos			= 0;


$mystring = $_[0];

for ($pos=1; $pos <= length ($mystring); $pos ++)
{
	$char= substr $mystring,$pos,1;
	$ord = ord($char);
	print ">$char<($ord), ";
}


}

