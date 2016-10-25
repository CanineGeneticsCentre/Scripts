#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use File::Basename ;

#####################################
# Define variables                  #
#####################################

my $command			= "";
my $prefix			= "";
my $mapfile			= "";
my $mapfile_recoded		= "";
my $pedfile			= "";
my $pedfile_recoded		= "";
my $genfile			= "";
my $samplefile			= "";
my $ans				= "";
my $parameters_ok		= "false";
my $species_code		= "";
my $species			= "";
my $analysis_suffix		= "";
my $results_directory		= "";
my $results_file		= "";
my $filename			= "";
my $logfile			= "";
my $final_readme_file		= "";
my $directory_suffix		= "";

my $maf_value			= 0.05;
my $geno_value			= 0.05;
my $mind_value			= 0.1;
my $pheno_string		= "";



############################
# process -f flag (if any) #
############################

GetOptions("f=s"=>\$prefix);


print "\n";
print "#########################################\n";
print "#          RUNSNPTEST QTL v5            #\n";
print "#                                       #\n";
print "#   This program converts PLINK files   #\n";
print "#   to SNPTEST format, and then runs    #\n";
print "#   SNPTEST association software in     #\n";
print "#   QTL mode.                           #\n";
print "#                                       #\n";
print "#########################################\n";
print "\n";


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

#####################
# Check files exist #
#####################

unless (-e "$prefix.ped")
{
	print "\n###############\n";
	print "# FILE ERROR! #\n";
	print "###############\n\n";
	print " Can't find file $prefix.ped\n\n";
	exit;
}
#####################
# Check files exist #
#####################

unless (-e "$prefix.map")
{
	print "\n###############\n";
	print "# FILE ERROR! #\n";
	print "###############\n\n";
	print " Can't find file $prefix.map\n\n";
	exit;
}



################################
# Get suffix for this analysis #
################################

print "Suffix for this analysis (e.g. a, b, c... or maybe the particular phenotype chosen):      ";
$analysis_suffix = <STDIN>;
chomp $analysis_suffix;

######################
# Set all file names #
######################

#Input ped and map files
$pedfile = $prefix.".ped";
$mapfile = $prefix.".map";

if ($analysis_suffix)
{
	$analysis_suffix = "_".$analysis_suffix;
	$filename = $prefix.$analysis_suffix;

}

else
{
	$filename = $prefix."_recoded";
}

$pedfile_recoded = $filename.".ped";
$mapfile_recoded = $filename.".map";


################
# Output files #
################
#Output file prefix

$genfile = $filename.".gen"; 
$samplefile = $filename.".sample"; 

$logfile = $filename.".log";
$results_file= $filename."_snptest.results";
$final_readme_file = $filename."_snptest.readme";
$results_directory = $filename."_snptest_results";



#############################
# Ask if it is dog or horse #
#############################

print "Which species are you analysing:    \n\n";
print "   <1> Dog\n";
print "   <2> Horse        ";


$species_code = <STDIN>;
chomp $species_code;
print"\n\n";

if ($species_code == 1){$species = "dog"}
if ($species_code == 2){$species = "horse"}


#############################
# Ask about parameters      #
#############################
 while ($parameters_ok eq "false")
{
	print "These are the default parameters in use:    \n\n";
	print "    --maf $maf_value\n";
	print "    --geno $geno_value\n";
	print "    --mind $mind_value\n";
	print "    phenotype string: $pheno_string\n";
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

		print "    PLINK-compatible string to use phenotype file (press return to leave as $pheno_string):    ";
		$ans = <STDIN>;
		chomp $ans;
		if ($ans ne ""){$pheno_string = $ans;}


	}
	else
	{
		$parameters_ok = "true";
	}

}


#########################
# Open my results file  #
#########################

open(RESULTS,">$results_file") or die "Can't open file $results_file";

print RESULTS "#########################################\n";
print RESULTS "   SNPTEST results file for $filename\n";
print RESULTS "#########################################\n\n";
print RESULTS scalar localtime;
print RESULTS "\n\n";

print RESULTS "Parameters in place:\n\n";
print RESULTS "--maf $maf_value\n";
print RESULTS "--geno $geno_value\n";
print RESULTS "--mind $mind_value\n";

print RESULTS "\n";
print "\n";



######################################
# Create the results folder          #
######################################

while (-e $results_directory)
{
	$directory_suffix = $directory_suffix + 1;
	$results_directory = $filename."_snptest_results"."_".$directory_suffix;
}

$command = "mkdir $results_directory";

print("$command\n");
system("$command");

print RESULTS "Creating new directory for SNPTEST results $results_directory\n";
print RESULTS " $command\n\n";



print "\n";
print "####################################################################\n";
print "# Running PLINK to recode the data into filtered MAP and PED files #\n";
print "# for input into the SNPTEST package...                            #\n";
print "####################################################################\n";
print "\n";

################################################
# Run PLINK to recode the ped file with these  #
# new filters and phenotypes in place.         #
#                                              #
# The new recoded and filtered file is then    #
# used to be run through SNPTEST.              #
#                                              #
# This avoids having to use SNPTEST's own      #
# filters which are quite limited.             #
################################################


$command =  "plink --noweb --".$species." --allow-no-sex --make-founders --nonfounders";
$command .= " --maf $maf_value --geno $geno_value --mind $mind_value";
$command .= " --recode --tab";
$command .= " --file $prefix $pheno_string  --out $filename";

print("$command\n");
system("$command");

print RESULTS "Running PLINK to recode the data into filtered MAP and PED files.\n";
print RESULTS " $command\n\n";



print "\n\n";
print "###################################\n";
print "#    GENERAL CHECK ON THE DATA     #\n";
print "#                                  #\n";
print "# Please check at this point that  #\n";
print "# the values you used for geno and #\n";
print "# maf have not filtered off too    #\n";
print "# many SNPs or samples.            #\n";
print "###################################\n";
print "\n\nPress return to continue (q to quit):    ";
$ans = <STDIN>;
chomp $ans;

if (($ans eq "q") || ($ans eq "Q")) 
{
	exit;
}





#######################
#  Run GTOOL for QTL  #
#######################
print "\n\nRunning GTOOL for QTL...\n\n";

$command = "/opt/snptest/gtool -P";  
$command .= " --map $mapfile_recoded";
$command .= " --ped $pedfile_recoded";
$command .= " --discrete_phenotype 0";
$command .= " --og $filename.gen";
$command .= " --os $filename.sample";


print("$command\n");
system("$command");


print RESULTS "Running GTOOL to convert PLINK files into SNPTEST files\n";
print RESULTS " $command\n\n";


###################
#  Rename files   #
###################
print "\nRenaming files...\n\n";

print RESULTS "Renamed various files\n";

$command = "mv $filename.gen_1 $filename.gen_controls";
print("$command\n");
system("$command");
print RESULTS " $command\n";

$command = "mv $filename.gen_2 $filename.gen_cases";
print("$command\n");
system("$command");
print RESULTS " $command\n";

$command = "mv $filename.sample_1 $filename.sample_controls";
print("$command\n");
system("$command");
print RESULTS " $command\n";

$command = "mv $filename.sample_2 $filename.sample_cases";
print("$command\n");
system("$command");
print RESULTS " $command\n\n";



#############################################
#  Run SNPTEST - Basic association for QTLs #
#############################################
print "\nRunning SNPTEST - Basic association...\n\n";

$command = "/opt/snptest/snptest";  
$command .= " -controls $filename.gen $filename.sample";
$command .= " -o $filename"."_QTL_frequentist.out";
$command .= " -frequentist 1 2 3 4 5";
$command .= " -qt";

print("$command\n");
system("$command");

print RESULTS "Running SNPTEST to carry out a basic QTL association analysis.\n";
print RESULTS " $command\n\n";


################################################
#  Run SNPTEST - Bayesian association for QTLs #
################################################
print "\nRunning SNPTEST - Bayesian association...\n\n";

$command = "/opt/snptest/snptest";  
$command .= " -controls $filename.gen $filename.sample";
$command .= " -o $filename"."_QTL_bayesian.out";
$command .= " -bayesian 1 2 3 4 5";
$command .= " -qt";

print("$command\n");
system("$command");

print RESULTS "Running SNPTEST to carry out a Bayesian QTL association analysis.\n";
print RESULTS " $command\n\n";

print "\nQTL ANALYSIS COMPLETED.\n\n";


print "There are two output files, both of which can be plotted using the Excel sheet SNPTEST_plot:\n\n";

print "$filename"."_QTL_frequentist.out\n";
#print "$filename"."_QTL_bayesian.out\n\n";
print "(NO BAYESIAN FILE AT PRESENT)\n\n";

print "Look at the file $filename.log for details of how the filtering has affected the numbers of samples and SNPs\n\n";



#######################################
# Copy files  into the results folder #
#######################################
print "Copying all relevant files into $results_directory\n\n";
print RESULTS "Copying all relevant files into $results_directory\n\n";

&move_to_results_folder (".log");
&move_to_results_folder ("_QTL_bayesian.out");
&move_to_results_folder ("_QTL_frequentist.out");
&move_to_results_folder ("_snptest.results");


$command = "mv  gtool.log $results_directory/gtool.log";
print("$command\n");
system("$command");

###########################################
# Remove some unwanted intermediate files #
###########################################

print "\n\nRemoving some unwanted files\n\n";
print RESULTS "\n\nRemoving some unwanted files\n\n";

&delete_unwanted (".gen");
&delete_unwanted (".sample");
&delete_unwanted (".ped");
&delete_unwanted (".map");
&delete_unwanted (".hh");
&delete_unwanted (".irem");


#################################################################
# Write the final README file to tell user which files are what #
#################################################################

open(README,">$final_readme_file") or die "Can't open file $final_readme_file";

print README "###############################\n";
print README "# This lists the output files #\n";
print README "# and what to do with them    #\n";
print README "###############################\n\n";

print README "(Please also read any instructions in the Genotyping Analysis/Instructions/ folder)\n\n";

print README "Filename: $filename\n";
print README scalar localtime;
print README "\n\n";

print README "-----------------------------------------------------------------------------------------\n";

# Normal #
print README "The following file has the results of a basic association analysis.\n";
print README "(This can be plotted using the Excel sheet SNPTEST_plot)\n\n";

print README " > $filename"."_QTL_frequentist.out\n\n";

print README "-----------------------------------------------------------------------------------------\n";

# Bayesian #
print README "The following file has the results of a Bayesian association analysis.\n";
print README "(This can be plotted using the Excel sheet SNPTEST_plot)\n\n";

print README " > $filename"."_QTL_bayesian.out\n\n";

print README "NB. THE BAYESIAN OPTION DOES NOT WORK WITH QTLS AT PRESENT!\n\n";
print README "-----------------------------------------------------------------------------------------\n";

close README;

&move_to_results_folder ("_snptest.readme");


###################################################################################
#                  End message to say the program has finished                    #
###################################################################################

print"\n\n\n";
print "##############################################\n";
print "##############################################\n";
print "##       SNPTEST QTL analysis complete      ##\n";
print "##############################################\n";
print "##############################################\n";
print "\n\nCheck the following files:\n\n";
print "(These have been moved into a folder called $filename"."_snptest_results)\n\n";

print "    1.  $filename"."_snptest.readme   	\t-  for how to use the output files\n";
print "    2.  $filename".".log      		\t-  for a log of the PLINK recoding\n";
print "    3.  $filename"."_snptest.results  	\t-  for more information.\n\n";
print "    4.  $filename"."_QTL_frequentist.out \t-  for basic association results.\n\n";
print "    5.  $filename"."_QTL_bayesian.out  	\t-  for Bayesian association results.\n\n";






exit;

#===================================================================================================================


#############################################
#                                           #
# Subroutine to move file to results folder #
#                                           #
#############################################

sub move_to_results_folder
{
	my $suffix = "";

	$suffix = $_[0];
	$command = "mv  $filename".$suffix." $results_directory/$filename".$suffix;
	print("$command\n");
	system("$command");
	print RESULTS ("$command\n");

}

#############################################
#                                           #
# Subroutine to remove unwanted files       #
#                                           #
#############################################

sub delete_unwanted
{
	my $suffix = "";

	$suffix = $_[0];
	$command = "rm  $filename".$suffix;
	print("$command\n");
	system("$command");
	print RESULTS ("$command\n");

}

