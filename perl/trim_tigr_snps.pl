#!/usr/local/bin/perl 


# Mike Boursnell Jul 2010
# Animal Health Trust
# Newmarket
# UK

use strict;
use Getopt::Long;

#####################################
# Define variables                  #
#####################################
my @item         		= ();
my @item_temp			= ();
my @newline_affected      	= ();
my @newline_normal      	= ();
my @snp_array			= ();
my @id_array			= ();
my @status_array		= ();
my @SNP_name_array		= ();
my @chromosome_array		= ();
my @position_array		= ();
my @third_array			= ();
my @row1 			= ();
my @row2 			= ();
my @row3 			= ();
my @row4 			= ();
my @standard_allele_A		= ();
my @standard_allele_B		= ();
my @mean_genotype_array		= ();

my @array_2d 			= (\@row1, \@row2);
my @array_2d_recoded		= (\@row3, \@row4);

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
my $new_SNP_name		= "";
my $chromosome			= "";
my $position			= "";
my $third			= "";
my $genotype			= "";
my $ans				= "";
my $snp_name_1			= "";
my $snp_name_2			= "";


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
my $changed_count		= 0;
my $unchanged_count		= 0;


############################
# process -f flag (if any) #
############################

GetOptions("f=s"=>\$infile);

print "\n";
print "################################################################\n";
print "#     Trim TIGR SNPs v1                                        #\n";
print "#                                                              #\n";
print "#     This edits SNP names like TIGR1292837_rs3451122          #\n";
print "#     and removes the part of the name after the underscore:   #\n";
print "#                                                              #\n";
print "#     So:  TIGR1292837_rs3451122 --> TIGR1292837               #\n";
print "#                                                              #\n";
print "################################################################\n\n";



#############################
# Get the input file name   #
#############################

if ($infile eq "")
{
	print "Prefix for PLINK .map file:  ";
	$infile = <STDIN>;
	chomp $infile;
	print"\n";
}

$mapfile = $infile.".map";
$pedfile = $infile.".ped";
$outfile = $infile."_trimmed.map";


############################
# Open files for output    #
############################

open (OUT, ">$outfile") || die "Cannot open $outfile";


############################
# Open the files for INPUT #
############################

open (MAP, "$mapfile") || die "Cannot open $mapfile";



################################
# Read data in MAP file to get # 
# the list of SNP names        #
################################

$snp_count = 0;
print "Reading MAP file...\n\n";

while ($single_line = <MAP>)
{

	chomp $single_line;

	@item=split(/\t/,$single_line);

	$snp_count = $snp_count + 1;

	$chromosome = $item[0];
	$SNP_name = $item[1];
	$third = $item[2];
	$position = $item[3];

	$chromosome_array[$snp_count]=$chromosome;
	$SNP_name_array[$snp_count]=$SNP_name;
	$third_array[$snp_count]=$third;
	$position_array[$snp_count]=$position;

}

close MAP;

$snp_total_in_mapfile = $snp_count;

print "$snp_total_in_mapfile SNPs in the MAP file\n";
print "\n";

###########################################
# Now trim names and write to output file #
###########################################

print "Editing file...\n";

for ($snp_count=1; $snp_count<=$snp_total_in_mapfile ;$snp_count++)
{
	$SNP_name = $SNP_name_array[$snp_count];
	$new_SNP_name = $SNP_name;

	if (index($SNP_name,"_") > 0)
	{
		if (index($SNP_name,"TIGR") == 0)
		{
			$new_SNP_name = substr($SNP_name,0,index($SNP_name,"_"));
		}
	} 


	if ($SNP_name eq $new_SNP_name)
	{
		$unchanged_count = $unchanged_count + 1;
	}

	if ($SNP_name ne $new_SNP_name)
	{
		$changed_count = $changed_count + 1;
	}

	print OUT "$chromosome_array[$snp_count]\t";
	print OUT "$new_SNP_name\t";
	print OUT "$third_array[$snp_count]\t";
	print OUT "$position_array[$snp_count]\n";

}
close OUT;


print "\n";

print "No of SNPs changed:\t$changed_count\n";

print "No of SNPs unchanged:\t$unchanged_count\n";

print "Output file:\t\t$outfile\n\n\n";
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

