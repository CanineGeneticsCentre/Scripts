#!/usr/bin/perl 

################################################################
#     Run EMMAX HD                                             #
#                                                              #
#     Runs Efficient Mixed-Model Association eXpedited (EMMAX) #
#                                                              #
################################################################


# Mike Boursnell Dec 2010
# Animal Health Trust
# Newmarket
# UK
# mike.boursnell@aht.org.uk

use strict;
use Getopt::Long;
use Statistics::Distributions;  # For CHI-SQUARED


#####################################
# Define variables                  #
#####################################
my $version				= "8";

my $command_line_mode	= "false";
my $prefix				= "";
my $command				= "";
my $outfile				= "";
my $phenofile			= "";
my $qq_file				= "";
my $kinfile				= "";
my $kinfile_BN			= "";
my $kinfile_IBS			= "";
my $kintype				= "IBS";
my $pedfile				= "";
my $PDF_file			= "";
my $psfile				= "";
my $gnufile				= "";
my $pedigree			= "";
my $id					= "";
my $status				= "";
my $single_line			= "";
my $parameters_ok		= "false";
my $ans					= "";
my $pheno				= "";
my $pheno_name			= ""; # if pheno file has columns
my $pheno_string		= "";
my $species				= "";
my $tfamfile			= "";
my $tpedfile 			= "";
my $mapfile				= "";
my $species_code		= "";
my $chromosome			= "";
my $SNP_name			= "";
my $position			= "";
my $psfile				= "";
my $final_outfile		= "";
my $final_outfile_hd	= "";
my $last_chromosome		= "";
my $is_y_chromosome 	= "";
my $lambda_string		= "";
my $title				= "";
my $species_string_plink = "";

my $geno_value			= 0.1;
my $maf_value			= 0.05;
my $mind_value			= 0.1;
my $line_count			= 0;
my $beta				= 0;
my $P_value				= 0;
my $minus_log_P			= 0;
my $snp_count_tped		= 0;
my $snp_count_ps		= 0;
my $rows_per_chromosome_count	= 0;
my $max_row_count		= 0;
my $array_count			= 0;
my $chromosome_count	= 0;
my $row_count			= 0;
my $total_no_chromosomes= 0;
my $chi_sq				= 0;
my $expected_chi_sq		= 0;
my $rank				= 0;
my $expected_P_value	= 0;
my $total_lines			= 0;
my $median				= 0;
my $lambda				= 0;
my $array_size			= 0;

my @item         				= ();
my @chr_array					= ();
my @snp_array					= ();
my @pos_array					= ();
my @rows_per_chromosome_array	= ();
my @new_line_array				= ();
my @minus_log_P_array			= ();
my @P_value_array				= ();
my @chi_sq_array				= ();
my @chi_sq_array_reverse		= ();
my @chi_sq_array_expected		= ();


###############################
# Process flags               #
###############################

GetOptions("file:s"=>\$prefix, "out:s"=>\$outfile, "species:s"=>\$species, "geno:s"=>\$geno_value, "maf:s"=>\$maf_value, "mind:s"=>\$mind_value, "pheno:s"=>\$pheno, "column:s"=>\$pheno_name);


print "\n";
print "###################################################################\n";
print "#     run_emmax HD                                                #\n";
print "#                                                                 #\n";
print "#     Runs Efficient Mixed-Model Association eXpedited (EMMAX)    #\n";
print "#                                                                 #\n";
print "###################################################################\n\n";

print "\n\nVersion $version\n\n";

#############################################################
# If you have specified the file name in the command line   #
# then assume this is being used in command line mode       #
#############################################################
if ($prefix ne "")
{
	$command_line_mode = "true";
	print "File:             \t$prefix\n";
	print "Output file:      \t$outfile\n";
	print "Species:          \t$species\n";
	print "GENO:             \t$geno_value\n";
	print "MAF:              \t$maf_value\n";
	print "MIND:             \t$mind_value\n";
	print "Pheno file:       \t$pheno\n";
	print "Pheno-name:       \t$pheno_name\n\n";
}



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

if ($outfile eq ""){$outfile = $prefix}

######################
# Make up file names #
######################

$phenofile = $prefix.".pheno";
$pedfile = $prefix.".ped";
$mapfile = $prefix.".map";
$tfamfile = $outfile.".tfam";
$tpedfile = $outfile.".tped";
$psfile = $outfile.".ps";
$final_outfile = $outfile.".emmax_out";
$final_outfile_hd = $outfile.".hd_emmax_out";

$kinfile_IBS = $prefix.".hIBS.kinf";
$kinfile_BN = $prefix.".hBN.kinf";

if ($kintype eq "IBS"){$kinfile = $kinfile_IBS}
if ($kintype eq "BN"){$kinfile = $kinfile_BN}

$qq_file = $prefix.".emmax_qq";

#############################
# Work out phenotype string #
#############################
if ($command_line_mode eq "true")
{
	if (($pheno ne "") && ($pheno_name ne "")){$pheno_string = "--pheno $pheno --pheno-name $pheno_name"}
	if (($pheno ne "") && ($pheno_name eq "")){$pheno_string = "--pheno $pheno"}
	if (($pheno eq "") && ($pheno_name eq "")){$pheno_string = ""}
}


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

if ($species eq "")
{
	print "####################################\n";
	print "# Which species are you analysing? #\n";
	print "####################################\n\n";
	print "   <1> Dog\n";
	print "   <2> Horse\n";
	print "   <3> Human\n\n";

	$species_code = <STDIN>;
	chomp $species_code;

	if ($species_code eq "1"){$species = "dog"}
	if ($species_code eq "2"){$species = "horse"}
	if ($species_code eq "3"){$species = "human"}
}

# PLINK species string
if ($species eq "dog"){$species_string_plink = "--dog"}
if ($species eq "horse"){$species_string_plink = "--horse"}
if ($species eq "human"){$species_string_plink = ""}
if ($species eq "cat"){$species_string_plink = ""}


#############################
# Ask about parameters      #
#############################
if ($command_line_mode eq "true"){$parameters_ok = "true"}

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

open(STDOUT, "| tee $prefix".".emmax_log.out");



#########################################################
# Convert plink files into transposed ped and fam files #
#########################################################

$command = "plink --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --file $prefix $pheno_string --maf $maf_value --geno $geno_value --mind $mind_value --recode12 --output-missing-genotype 0 --transpose --out $prefix";

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

print "###########################################\n";
print "# Created phenotype file $phenofile     \n";
print "###########################################\n\n";

close TFAM;
close PHENO;

###############################
# Create kinship matrix (IBS) #
###############################

$command = "/opt/emmax/emmax-kin -v -h -s -d 10 $prefix";

print("$command\n");
system("$command");

print "################################\n";
print "# Created kinship matrix (IBS) #\n";
print "################################\n\n";



###############################
# Create kinship matrix (BN)  #
###############################

$command = "/opt/emmax/emmax-kin -v -h -d 10 $prefix";

print("$command\n");
system("$command");

print "################################\n";
print "# Created kinship matrix (BN)  #\n";
print "################################\n\n";




##############
# Run EMMAX  #
##############

$command = "/opt/emmax/emmax -v -d 10 -t $prefix -p $phenofile -k $kinfile -o $outfile";

print("$command\n");
system("$command");


print "\n";
print "########################\n";
print "# EMMAX run completed. #\n";
print "########################\n\n";

print "Getting map positions from tped file...\n\n";


######################################################
# Reformat output file for plotting in PLINK_plot_HD #
######################################################


##########################################
# First get map positions from tped file #
# (use this rather than the map file as  #
# it is filtered and in order)           #
##########################################


open (TPED, "$tpedfile") || die "Cannot open $tpedfile";

$line_count = 0;
$rows_per_chromosome_count = 0;
$max_row_count = 0;

while ($single_line = <TPED>) {

	chomp $single_line;
	$line_count = $line_count + 1;
	$rows_per_chromosome_count = $rows_per_chromosome_count + 1;

	#########################################################
        # Split line at TABS OR SPACES into the array 'item'    #
	#########################################################
        @item=split(/\s+/,$single_line);

	$chromosome = $item[0];
	$SNP_name = $item[1];
	$position = $item[3];


	#################################
	# Deal with X and Y chromosomes #
	#################################
	if ($species eq "dog")
	{
		if ($chromosome eq "X"){$chromosome=39}
		if ($chromosome eq "Y"){$chromosome=40;$is_y_chromosome = "yes"}
	}
	if ($species eq "horse")
	{
		if ($chromosome eq "X"){$chromosome=32}
		if ($chromosome eq "Y"){$chromosome=33; $is_y_chromosome = "yes"}
	}
	if ($species eq "human")
	{
		if ($chromosome eq "X"){$chromosome=23}
		if ($chromosome eq "Y"){$chromosome=24; $is_y_chromosome = "yes"}
	}
	if ($species eq "cat")
	{
		if ($chromosome eq "X"){$chromosome=19}
		if ($chromosome eq "Y"){$chromosome=20; $is_y_chromosome = "yes"}
	}



	$chr_array[$line_count] = $chromosome;
	$snp_array[$line_count] = $SNP_name;
	$pos_array[$line_count] = $position;

	################################################
	# If this is the first row of a new chromosome #
	################################################
	if ($chromosome ne $last_chromosome)
	{
		$rows_per_chromosome_array[$last_chromosome]= $rows_per_chromosome_count;

		if ($rows_per_chromosome_count > $max_row_count)
		{
			$max_row_count = $rows_per_chromosome_count;
		}

		$rows_per_chromosome_count = 0;
	}


	$last_chromosome = $chromosome;

} # end of while loop



#########################################
# Store no of rows for final chromosome #
#########################################
$rows_per_chromosome_array[$chromosome]= $rows_per_chromosome_count + 1;


$snp_count_tped = $line_count;

close TPED;


################################
# Decide number of chromosomes #
################################
if ($species eq "dog")
{
	if ($is_y_chromosome eq "yes"){$total_no_chromosomes = 40} else {$total_no_chromosomes = 39}
}
if ($species eq "horse")
{
	if ($is_y_chromosome eq "yes"){$total_no_chromosomes = 33} else {$total_no_chromosomes = 32}
}
if ($species eq "human")
{
	if ($is_y_chromosome eq "yes"){$total_no_chromosomes = 24} else {$total_no_chromosomes = 23}
}
if ($species eq "cat")
{
	if ($is_y_chromosome eq "yes"){$total_no_chromosomes = 20} else {$total_no_chromosomes = 19}
}


##########################################################
# Now open EMMAX output .ps file and merge with the data #
# from the TPED file to produce an ARRAY with the EMMAX  #
# data and the map positions and minus log P values      #
##########################################################

print "Formatting data for output file...\n\n";

open (PS, "$psfile") || die "Cannot open $psfile";
open (FINALHD, ">$final_outfile_hd") || die "Cannot open $final_outfile_hd";
open (FINAL, ">$final_outfile") || die "Cannot open $final_outfile";


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
		print FINALHD "\n\nERROR!  SNP in TPED file should match SNP in PS file!\n\n";
		close PS;
		close FINALHD;
		exit;
	}

	###################################
	# Store all this output in arrays #
	###################################

	$P_value_array[$line_count] = $P_value;
	$minus_log_P_array[$line_count] = $minus_log_P;

	#print FINALHD "$line_count\t$chr_array[$line_count]\t$SNP_name\t$pos_array[$line_count]\t$P_value\t$minus_log_P\n";

} # end of while loop

$snp_count_ps = $line_count;


####################################################
# Convert to HD format (one column per chromosome) #
# Work through each chromosome, adding a           #
# bit to each newline as you work down             #
####################################################

#print FINALHD "$line_count\t$chr_array[$line_count]\t$SNP_name\t$pos_array[$line_count]\t$P_value\t$minus_log_P\n";

$array_count =0;
$line_count = 0;


for ($chromosome_count = 1;$chromosome_count<=$total_no_chromosomes;$chromosome_count++)
{


	######################################################
	# Add the correct number of rows for each chromosome #
	######################################################


	for ($row_count = 1; $row_count <= $rows_per_chromosome_array[$chromosome_count]; $row_count++)
		{
			$array_count = $array_count + 1;


			if ($chromosome_count == 1)
			{
				$new_line_array[$row_count] = $new_line_array[$row_count]."$array_count\t$chr_array[$array_count]";
			}
			if ($chromosome_count > 1)
			{ 
				$new_line_array[$row_count] = $new_line_array[$row_count]."\t$array_count\t$chr_array[$array_count]";
			}

			# Add other bits to new line
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$snp_array[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$pos_array[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$P_value_array[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$minus_log_P_array[$array_count]";


		}

	##################################################
	# Now fill up remaining newline arrays with tabs #
	##################################################
	for ($row_count = $rows_per_chromosome_array[$chromosome_count] + 1; $row_count <= $max_row_count; $row_count++)
	{
		if ($chromosome_count == 1) {$new_line_array[$row_count] = $new_line_array[$row_count]."\t\t\t\t\t"}
		if ($chromosome_count > 1) {$new_line_array[$row_count] = $new_line_array[$row_count]."\t\t\t\t\t\t"}
	}
	
}


######################################################################
# Write from final array to an output file in columns for Excel 2003 #
# (One column per chromosome)                                        #
######################################################################


for ($chromosome_count = 1;$chromosome_count<=$total_no_chromosomes - 1;$chromosome_count++)
	{
		print FINALHD "INDEX\tCHR$chromosome_count\tSNP\tPOS\tP_value\tminus_logP_EMMAX\t";

	}

	if ($is_y_chromosome eq "yes")
	{
		print FINALHD "INDEX\tCHRY\tSNP\tPOS\tP_value\tminus_logP\t";
	}

	print FINALHD "INDEX\tCHRL\tSNP\tPOS\tP_value\tminus_logP\n";

for ($row_count = 1;$row_count <= $max_row_count;$row_count++)
{
	print FINALHD "$new_line_array[$row_count]\n";
}


############################################################################
# Write from final array to an output file in single column for Excel 2007 #
# or GNUplot                                                               #
############################################################################

print FINAL "CHR\tSNP\tPOS\tP_value\tminus_logP_EMMAX\n";

for ($array_count = 1; $array_count <=$snp_count_tped; $array_count++)
{
	print FINAL "$chr_array[$array_count]\t$snp_array[$array_count]\t$pos_array[$array_count]\t$P_value_array[$array_count]\t$minus_log_P_array[$array_count]\n";
}



close FINAL;
close FINALHD;
close PS;


###########################
# Make output for QQ plot #
###########################

print "Calculating data for QQ plot...\n\n";

# First make chisq values from P values
for ($array_count = 0; $array_count <= $snp_count_ps; $array_count++)
{
	if (($P_value_array[$array_count] * 1) > 0)
	{
		$chi_sq = Statistics::Distributions::chisqrdistr(1,$P_value_array[$array_count]);
		$chi_sq_array[$array_count] = $chi_sq;
	}
}


##################
# Sort the array #
##################

print "Sorting file...\n\n";

@chi_sq_array_reverse = reverse sort {$a <=> $b} (@chi_sq_array);



#################################
# Print first 10 lines of array #
#################################

print "Showing the first 10 lines of the sorted array...\n\n";

for($line_count = 0; $line_count <=10; $line_count++)
{
	print "$line_count\t$chi_sq_array_reverse[$line_count]\n";
}

print "\n";



##################
#Stats reminders #
##################
#$expected_chi_sq = Statistics::Distributions::chisqrdistr(1,$expected_P_value);
#$expected_P_value = Statistics::Distributions::chisqrprob(1,$expected_chi_sq);



##############################################
# Now calculate Exp_p with this formula:     #
#                                            #
# Exp-p = (rank - 0.5)/total_no_of_lines     #
# and then getting the CHISQ of that P-value #
##############################################
$total_lines = $snp_count_ps;

for($array_count = 0; $array_count <= $total_lines; $array_count++)
{
	$rank = $array_count + 1;
	$expected_P_value = ($rank - 0.5) / $total_lines;

	if ($expected_P_value <= 1)
	{
		$expected_chi_sq = Statistics::Distributions::chisqrdistr(1,$expected_P_value);
		$chi_sq_array_expected[$array_count] = $expected_chi_sq;
	}

	
}


###############################################
# Calculate the median of the observed values #
###############################################
$array_size = @chi_sq_array_reverse;
$median = $chi_sq_array_reverse[int($array_size/2)];
$lambda = $median / 0.456;

$lambda_string = sprintf("%.3f",$lambda);

print "Median: \t$median\n";
print "Lambda: \t$lambda_string\n\n\n"; 
 



############################
# Open files for output    #
############################

open (OUT, ">$qq_file") || die "Cannot open $qq_file";

print OUT "EMMAX\tBLANK\tEXP_CHI\tOBS_CHI\n";

for ($line_count = 0; $line_count <= $total_lines; $line_count++)
{
	if (($chi_sq_array_reverse[$line_count] * 1) > 0)
	{
		print OUT "$chi_sq_array_reverse[$line_count]\t\t$chi_sq_array_expected[$line_count]\t$chi_sq_array_reverse[$line_count]\n";
	}
}

print OUT "\t\t\t\t\n";
print OUT "\t\t\t$lambda_string";

close OUT;


###############################################
# Make GNUplot file for pdf figure of QQ plot #
###############################################

$gnufile = "gnufile.xtxt";
$psfile = $prefix."_emmax_qq.ps";
$PDF_file = $prefix."_emmax_qq.pdf";
$title = "$prefix:  QQ plot of EMMAX -corrected chi-squared values";

open (GG, ">$gnufile") || die "can't open $gnufile\n" ;

##################################################
# Write the gnuplot commands to the command file #
##################################################
print GG "set terminal postscript color solid\n";
print GG "set title  \"$title\" \n" ; 
print GG "set key outside\n"; 
print GG "set xlabel  \"Expected\" \n" ; 
print GG "set ylabel  \"Observed\" \n" ; 
print GG "set nokey\n";

print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 3\n";
print GG "set style line 2 lt -1 lw 1 linecolor 2\n";

print GG "set label 1 'Inflation factor: $lambda_string' at graph 0.5,0.1\n";

print GG "plot '$qq_file' using 3:3 with lines ls 2 title 'NULL', '$qq_file' using 2:3 with points ls 1 title 'QQ'" ;


###############################
# Run the commands from linux #
###############################
print "Running GNUplot to create pdf file...\n\n";

system "gnuplot < $gnufile > $psfile" ;
#system "fixgreen  $psfile" ;
system "ps2pdf  $psfile " ;



print "#########################################################\n";
print "#               EMMAX analysis completed                #\n";
print "#########################################################\n\n";

print "Number of SNPs in TPED file:\t$snp_count_tped\n";
print "Number of SNPs in PS file:  \t$snp_count_ps\n";
print "Total number of chromosomes: \t$total_no_chromosomes\n";

if ($is_y_chromosome eq "yes")
{
	print "Y chromosome found:\t\tYes\n"
} 
else
{
	print "Y chromosome found:\t\tNo\n"
}
print "Inflation factor:\t\t$lambda_string\n\n";

#sprintf("%.3f", $number);

print "The output files from EMMAX are these:\n\n";

print ">>\t".$outfile.".reml\tREML output. The last line denotes the pseudo-heritability estimates.\n";
print ">>\t".$outfile.".ps  \tP-value file. Columns are SNP, beta, P-value.\n\n";

print "There is also an output file which has chromosome and position columns, and a minus log P column.\n";

print "This must be plotted using PLINK_plot_HD (the data has separate groups of columns for each chromosome)\n\n";

print ">>\t$final_outfile_hd\n\n";


print "There is also an output file which has chromosome and position columns, and a minus log P column.\n";

print "This must be plotted using PLINK_plot (the data has is not formatted in separate columns for each chromosome)\n\n";

print ">>\t$final_outfile\n\n";


print "There is also an output file for plotting a QQ-plot:\n\n";

print ">>\t$qq_file    (This can be used in the Excel 2003 file QQ_plot)\n\n";

print "There is a PDF file with a picture of the QQ-plot:\n\n";

print ">>\t$PDF_file\n\n";


print "(If you want a single column of data and you have less than 65000 SNPs , use run_emmax_not_HD)\n\n";






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

