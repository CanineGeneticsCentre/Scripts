#!/usr/local/bin/perl 

################################################################
#     Run EMMAX v2 (not HD)                                    #
#                                                              #
#     Runs Efficient Mixed-Model Association eXpedited (EMMAX) #
#                                                              #
################################################################


# Mike Boursnell Sept 2010
# Animal Health Trust
# Newmarket
# UK
# mike.boursnell@aht.org.uk

use strict;
use Getopt::Long;

#####################################
# Define variables                  #
#####################################

my $prefix			= "";
my $command			= "";
my $outfile			= "";
my $phenofile			= "";
my $kinfile			= "";
my $kinfile_BN			= "";
my $kinfile_IBS			= "";
my $kintype			= "IBS";
my $pedfile			= "";
my $pedigree			= "";
my $id				= "";
my $status			= "";
my $single_line			= "";
my $parameters_ok		= "false";
my $ans				= "";
my $pheno_string		= "";
my $species			= "";
my $tfamfile			= "";
my $tpedfile 			= "";
my $mapfile			= "";
my $species_code		= "";
my $chromosome			= "";
my $SNP_name			= "";
my $position			= "";
my $psfile			= "";
my $final_outfile		= "";

my $geno_value			= 0.1;
my $maf_value			= 0.05;
my $mind_value			= 0.1;
my $line_count			= 0;
my $beta			= 0;
my $P_value			= 0;
my $minus_log_P			= 0;
my $snp_count_tped		= 0;
my $snp_count_ps		= 0;

my @item         		= ();
my @chr_array			= ();
my @snp_array			= ();
my @pos_array			= ();


###############################
# Process flags               #
###############################

GetOptions("file=s"=>\$prefix,"out=s"=>\$outfile);


print "\n";
print "###################################################################\n";
print "#     run_emmax v2   (not HD)                                     #\n";
print "#                                                                 #\n";
print "#     Runs Efficient Mixed-Model Association eXpedited (EMMAX)    #\n";
print "#                                                                 #\n";
print "#     (This version produces the output file in a single set of   #\n";
print "#      columns, so for the HD array will be too long to be        #\n";
print "#      used directly in Excel 2003)                               #\n";
print "###################################################################\n\n";


#############################
# Get the input file name   #
#############################

if ($prefix eq "")
{
	print "Prefix for pair of PLINK .map and .ped files:  ";
	$prefix = <STDIN>;
	chomp $prefix;
	print"\n";
}



######################
# Make up file names #
######################
$outfile = $prefix;
$phenofile = $prefix.".pheno";
$pedfile = $prefix.".ped";
$mapfile = $prefix.".map";
$tfamfile = $outfile.".tfam";
$tpedfile = $outfile.".tped";
$psfile = $outfile.".ps";
$final_outfile = $outfile.".emmax_out";

$kinfile_IBS = $prefix.".hIBS.kinf";
$kinfile_BN = $prefix.".hBN.kinf";

if ($kintype eq "IBS"){$kinfile = $kinfile_IBS}
if ($kintype eq "BN"){$kinfile = $kinfile_BN}



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



#############################
# Ask if it is dog or horse #
#############################

print "####################################\n";
print "# Which species are you analysing? #\n";
print "####################################\n\n";
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
	print "####################\n";
	print "# PLINK parameters #\n";
	print "####################\n";

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


#################
#TURN LOGGER ON #
#################

$| = 1;

open(STDOUT, "| tee $prefix".".emmax_log");



#########################################################
# Convert plink files into transposed ped and fam files #
#########################################################


$command = "plink --noweb --".$species." --allow-no-sex --make-founders --nonfounders --file $prefix $pheno_string --maf $maf_value --geno $geno_value --mind $mind_value --recode12 --output-missing-genotype 0 --transpose --out $prefix";

print("$command\n");
system("$command");

print "###########################################################\n";
print "# Converted plink files into transposed ped and fam files #\n";
print "###########################################################\n\n";



#################################################################
# Make pheno file from tfam file by using the first two columns #
# and the 6th column (phenotype)                                #
#################################################################

print "Creating phenotype file from pedfile...\n\n";

open (TFAM, "$tfamfile") || die "Cannot open $tfamfile";
open (PHENO, ">$phenofile") || die "Cannot open $phenofile";


while ($single_line = <TFAM>) {

	chomp $single_line;
	#&chomp_all ($single_line);

	#########################################################
        # Split line at TABS OR SPACES into the array 'item'    #
	#########################################################
        @item=split(/\s+/,$single_line);

	$pedigree = $item[0];
	$id = $item[1];
	$status = $item[5];

	print PHENO "$pedigree\t$id\t$status\n";

} # end of while loop

close TFAM;
close PHENO;



###############################
# Create kinship matrix (IBS) #
###############################

$command = "emmax-kin -v -h -s -d 10 $prefix";

print("$command\n");
system("$command");

print "################################\n";
print "# Created kinship matrix (IBS) #\n";
print "################################\n";



###############################
# Create kinship matrix (BN)  #
###############################

$command = "emmax-kin -v -h -d 10 $prefix";

print("$command\n");
system("$command");

print "################################\n";
print "# Created kinship matrix (BN)  #\n";
print "################################\n";




##############
# Run EMMAX  #
##############

$command = "emmax -v -d 10 -t $prefix -p $phenofile -k $kinfile -o $outfile";

print("$command\n");
system("$command");


###################################################
# Reformat output file for plotting in PLINK_plot #
###################################################

##########################################
# First get map positions from tped file #
# (use this rather than the map file as  #
# it is filtered and in order)           #
##########################################


open (TPED, "$tpedfile") || die "Cannot open $tpedfile";

$line_count = 0;

while ($single_line = <TPED>) {

	chomp $single_line;
	$line_count = $line_count + 1;

	#########################################################
        # Split line at TABS OR SPACES into the array 'item'    #
	#########################################################
        @item=split(/\s+/,$single_line);

	$chromosome = $item[0];
	$SNP_name = $item[1];
	$position = $item[3];

	$chr_array[$line_count] = $chromosome;
	$snp_array[$line_count] = $SNP_name;
	$pos_array[$line_count] = $position;

	#print PHENO "$pedigree\t$id\t$status\n";

} # end of while loop

$snp_count_tped = $line_count;

close TPED;


#####################################################
# Now open EMMAX output .ps file and merge          #
# these two to produce an output with the EMMAX     #
# data and the map positions and minus log P values #
#####################################################

print "Writing final output file...\n\n";

open (PS, "$psfile") || die "Cannot open $psfile";
open (FINAL, ">$final_outfile") || die "Cannot open $final_outfile";

########################
# Write column headers #
########################

print FINAL "INDEX\tCHR\tSNP\tPOS\tP_value\tminus_logP\n";


$line_count = 0;

while ($single_line = <PS>) {

	chomp $single_line;
	$line_count = $line_count + 1;

	#########################################################
        # Split line at TABS OR SPACES into the array 'item'    #
	#########################################################
        @item=split(/\s+/,$single_line);

	$SNP_name = $item[0];
	$beta = $item[1];
	$P_value = $item[2];


	#########################
	# Calculate minus log P #
	#########################

	$minus_log_P = log($P_value)/log(10) * -1;



	if ($SNP_name ne $snp_array[$line_count])
	{
		print "\n\nERROR!  SNP in TPED file should match SNP in PS file!\n\n";
		print FINAL "\n\nERROR!  SNP in TPED file should match SNP in PS file!\n\n";
		close PS;
		close FINAL;
		exit;
	}

	print FINAL "$line_count\t$chr_array[$line_count]\t$SNP_name\t$pos_array[$line_count]\t$P_value\t$minus_log_P\n";

} # end of while loop

$snp_count_ps = $line_count;

close FINAL;
close PS;

print "Number of SNPs in TPED file:\t$snp_count_tped\n";
print "Number of SNPs in PS file:  \t$snp_count_ps\n\n\n";



print "###########################################\n";
print "#        EMMAX analysis completed         #\n";
print "###########################################\n";

print "The output files from EMMAX are these:\n\n";

print "\t".$outfile.".reml\tREML output. The last line denotes the pseudo-heritability estimates.\n";
print "\t".$outfile.".ps  \tP-value file. Columns are SNP, beta, P-value.\n\n";

print "There is also a file which has chromosome and position columns, and a minus log P column.\n";
print "This can be plotted using PLINK_plot.\n\n";

print "\t$final_outfile\n\n";

print "I haven't finished the HD version yet but I'm working on it...\n\n";

print "\n\n\n";


#############################
# Turn off logging          #
#############################

close(STDOUT);

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

