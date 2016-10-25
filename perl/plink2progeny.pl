#!/usr/bin/perl

################################################################
#     plink2progeny v2                                         #
#                                                              #
#     Converts files from PLINK format (map and ped files)     #
#     to a format which can be imported into Progeny           #
#                                                              #
################################################################


# Mike Boursnell Aug 2007

use strict;
use Getopt::Long;


#####################################
# Define variables                  #
#####################################
my @item_ped        		= ();
my @item_list			= ();
my @piece			= ();
my @names     			= ();
my @newnames			= ();
my @time_array			= ();
my @snp_array			= ();
my @map_array			= ();
my @aff_array			= ();
my @item_snp			= ();
my @item_aff			= ();
my @SNP_position_array		= ();
my @sample_name_array		= ();
my @status_array		= ();

my $infile			= "";
my $pedfile       		= "";
my $outfile      		= "";
my $logfile			= "";
my $mapfile			= "";
my $mutfile			= "";
my $listfile			= "";
my $affectfile			= "";
my $single_line			= "";
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

my $help         		= 0;
my $linecount    		= 0;
my $colcount     		= 0;
my $count			= 0;
my $showcount    		= 0;
my $itemcount    		= 0;
my $total_items   		= 0;
my $total_items_row_1		= 0;
my $renamed_count		= 0;
my $total_lines  		= 0;
my $total_samples		= 0;
my $total_items  		= 0;
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


print"\n";
print"######################################################\n";
print"#                                                    #\n";
print"#     plink2progeny v2                               #\n";
print"#                                                    #\n";
print"#     Converts ped and map files from PLINK to       #\n";
print"#     a format which can be imported into Progeny    #\n";
print"#                                                    #\n";
print"#     The format is one row per sample, and can      #\n";
print"#     be imported using the 'Custom Import' format   #\n";
print"#     of Progeny.                                    #\n";
print"#                                                    #\n";
print"######################################################\n\n";



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


open (PED, "$pedfile") || die "Cannot open input file: $pedfile";
open (MAP, "$mapfile") || die "Cannot open input file: $mapfile";


print "Name of output file: ";
$outfile = <STDIN>;
chomp $outfile;
print"\n";
open (OUT, ">$outfile")|| die "Cannot create output file: $!";


#######################################
# Read in all SNPs in map file        #
# and put them in snp_array           #
#######################################
print "Reading MAP file...\n\n";
$count =0;
@map_array= <MAP>;
foreach $single_line (@map_array)
{
	$count = $count + 1;
	@item_snp=split(/\t/,$single_line);
	$snp_name= $item_snp[1];
	chomp $snp_name;
	$snp_array[$count]=$snp_name;
	#print "$count:   snp_name: $snp_name       snp_array[$count]: $snp_array[$count]\n";
}
close MAP;

$snp_total = scalar (@map_array);
print "  $snp_total SNPs in the MAP file\n";
print "\n\n";

########################################################
# Now write the headers for first six columns          #
########################################################

print OUT "Ped\tID\tFather\tMother\tSex\tStatus\t";


########################################################
# Now write the headers to the columns in the ped file #
# (one for each SNP)                                   #
########################################################

for ($count=1;$count <=$snp_total;$count++)
{
	print OUT "$snp_array[$count]\t";
}
print OUT "\n";

##############################
# Read data in PED file    #
##############################

print "Reading PED file...\n";

###########################################
# Work down the lines in the PED file     #
###########################################
$linecount=0;
while ($single_line = <PED>) {

	# Remove the end of line character #
	chomp $single_line;
	$linecount = $linecount + 1;

	##########################################
	# Write this single_line to the OUT file #
	##########################################

	print OUT "$single_line\n";

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

	print "\nSample: $ID     ";


	#####################################################
	# Check how many items there are on the row.        #
	#####################################################
	$itemcount = 0;

	while (length $item_ped[$itemcount]>0)
	{
		$itemcount = $itemcount + 1;
	}

	$total_items = $itemcount;

	
	#print "Row $linecount has $total_items items\n";
	#system("tput cuu 1");

	
	######################################################
	# Store the number of items on the first line        #
	# so we can check that all other lines have          #
	# the same number of items                           #
	######################################################
	if ($linecount==1) 
	{
		$total_items_row_1=$itemcount;
		#print "Number of items on the first row is $total_items_row_1.\n\n";
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
		print "Number of items on line $linecount is $itemcount\n";
		print "which is not the same as the $total_items_row_1 items on the first row.\n\n";
		$ans= <STDIN>;
		print "\n";
	}



} # end of while loop

$total_lines=$linecount;


print "\n\n";

print "  $total_lines samples in the PED file\n\n";




###########################
#     Close both files    #
###########################

close MAP;
close PED;
close OUT;

print "\n";
print "Input ped file:     	$pedfile\n";
print "Input map file:     	$mapfile\n";
print "Output progeny file:     $outfile\n\n";

exit;

