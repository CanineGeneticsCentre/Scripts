#!/usr/local/bin/perl

################################################################
#     Check pedfile format v2                                            #
#                                                              #
#     This reads in a file and gives information               #
#     about the format                                         #
################################################################


# Mike Boursnell Aug 2007

use strict;
use Getopt::Long;

#####################################
# Define variables                  #
#####################################
my @item         		= ();
my @item_map        		= ();
my @newline      		= ();
my @linearray    		= ();
my @piece			= ();

my $prefix			= "";
my $pedfile			= "";
my $mapfile			= "";
my $outFile      		= "";
my $answer       		= "";
my $single_line			= "";
my $separator 			= "";
my $char			= "";
my $chromosome			= "";
my $chromosome_map		= "";
my $SNP_name			= "";
my $position			= "";
my $cM				= "";
my $pedigree			= "";
my $id				= "";
my $father			= "";
my $mother			= "";
my $gender			= "";
my $status			= "";
my $SNP_name			= "";

my $help         		= 0;
my $linecount    		= 0;
my $colcount     		= 0;
my $showcount    		= 0;
my $itemcount    		= 0;
my $checkcount   		= 0;
my $total_lines  		= 0;
my $total_items  		= 0;
my $charcount			= 0;
my $no_tabs			= 0;
my $no_spaces			= 0;
my $no_both			= 0;
my $snp_count			= 0;
my $line_count			= 0;
my $no_columns			= 0;
my $no_snps			= 0;
my $total_lines_map		= 0;


print"\n";
print"################################################################\n";
print"#     Check pedfile v3                                         #\n";
print"#                                                              #\n";
print"#     This checks a pair of PED and MAP files                  #\n";
print"#                                                              #\n";
print"################################################################\n\n";


#############################
# Get the file name         #
#############################

print "Prefix for PED and MAP files:  ";
$prefix = <STDIN>;
chomp $prefix;

$pedfile = $prefix.".ped";
$mapfile = $prefix.".map";



print"\n";

###########################
# Check the MAP file      #
###########################

open (MAP, "$mapfile") || die "Cannot open: $pedfile";

print "Checking MAP file...\n\n";

print "CHR\tSNP\t\tcM\t\tPOS\n\n";

while ($single_line = <MAP>)
{

	$line_count = $line_count + 1;

	chomp $single_line;

	@item_map=split(/\s+/,$single_line);
	$chromosome_map = $item_map[0];
	$SNP_name = $item_map[1];
	$cM = $item_map[2];
	$position = $item_map[3];

	if ($line_count < 11)
	{
		print "$chromosome_map\t$SNP_name\t\t$cM\t\t$position\n";
	}

	if (lc($cM) eq "nan")
	{
		print "\n\nFound Nan in cM column.\n\n";
		
	}
}

close MAP;

$total_lines_map = $line_count;

print "\nTotal number of lines in MAP file: $total_lines_map\n\n";



#################################
# Count total columns and lines #
# (delimiter tab or space)      #
#################################
$line_count = 0;

print "Checking PED file (columns split with tabs or spaces)...\n";
open (IN, "$pedfile");

print "\nLINE\tFID \tIID     \FATHER\tMOTHER\tGENDER\tSTATUS\tCOLUMNS\tSNPs\n";

while ($single_line = <IN>)
{
	$line_count = $line_count + 1;

	#########################################################
        # Split line at TABS OR SPACES into the array 'item'    #
	#########################################################
        @item=split(/\s+/,$single_line);

	$no_both = @item;
	$pedigree = $item[0];
	$id = $item[1];
	$father=$item[2];
	$mother = $item[3];
	$gender = $item[4];
	$status = $item[5];

	$no_snps = ($no_both - 6)/2;

	if ($line_count == 1) {$no_columns = $no_both}

	if ($line_count % 10 == 0) {print "\nLINE\tFID \tIID     \FATHER\tMOTHER\tGENDER\tSTATUS \tCOLUMNS\tSNPs\n"};

	print "$line_count\t$pedigree\t$id\t$father\t$mother\t$gender\t$status\t$no_both\t$no_snps\n";

	if ($no_both != $no_columns)
	{
		print "\n\nERROR.  Found $no_both columns but expected $no_columns (as in first line)\n\n";
	}


}

close IN;

if ($total_lines_map != $no_columns)
{
	print "\n\nERROR.  Found $no_both columns but expected $total_lines_map (number of SNPs in MAP file)\n\n";
}


print "\n\n";

print "Number of lines in PED file:    \t$total_lines\n";
print "Number of columns in PED file:  \t$no_columns\n\n";

print "Number of lines in MAP file:    \t$total_lines_map\n\n";

###########################
#     Close      file     #
###########################
close IN;

exit;

