#!/usr/bin/perl -w

use strict;
use Getopt::Std ;
use File::Basename ;

####################
# Define variables #
####################

my $filename			= "";
my $prefix			= "";
my $ans				= "";
my $species_code		= "";
my $species			= "";
my $command			= "";
my $logfile			= "";
my $logfile_all			= "";
my $results_file		= "";
my $final_readme_file		= "";
my $single_line			= "";
my $errors_found		= "";
my $results_directory		= "";
my $suffix			= "";
my $analysis_suffix		= "";
my $QTL				= "";
my $parameters_ok		= "false";
my $pheno_string		= "";

my $error_pos			= 0;
my $haplotype_count		= 0;

my $mperm_value			= 10000;
my $geno_value			= 0.1;
my $maf_value			= 0.05;
my $mind_value			= 0.1;
my $start_haplotype		= 1;
my $end_haplotype		= 8;
my $directory_suffix		= 0;
my $QTL_code			= 1;


print "\n";
print "###############################################\n";
print "#        Run eigensoft_all v7                 #\n";
print "#                                             #\n";
print "#     This program runs all the programs      #\n";
print "#     in the EIGENSOFT 3.0 package.           #\n";
print "#                                             #\n";
print "#  1. CONVERTF (to convert PLINK files)       #\n";
print "#                                             #\n";
print "#  2. SMARTPCA (PCA analysis)                 #\n";
print "#                                             #\n";
print "#  3. EIGENSTRAT                              #\n";
print "#                                             #\n";
print "#  A list of all the output files and what    #\n";
print "#  to do with them can be found  at the end   #\n";
print "#  of the run                                 #\n";
print "#                                             #\n";
print "###############################################\n";
print "\n";


############################
# Get the file name prefix #
############################

print "Prefix of ped and map file:      ";
$prefix = <STDIN>;
chomp $prefix;
print "\n";


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


#############################
# Ask if it is dog or horse #
#############################

print "Which species are you analysing:    \n\n";
print "   <1> Dog\n";
print "   <2> Horse\n\n";
print "   > ";

$species_code = <STDIN>;
chomp $species_code;
print"\n\n";

if ($species_code == 1){$species = "dog"}
if ($species_code == 2){$species = "horse"}



#############################
# Ask if it is a QTL        #
#############################

print "What type of trait are you analysing:    \n\n";
print "   <1> Binary\n";
print "   <2> Quantitative\n\n";
print "   > ";


$QTL_code = <STDIN>;
chomp $QTL_code;
print"\n\n";

if ($QTL_code == 1){$QTL = "no"}
if ($QTL_code == 2){$QTL = "yes"}


#####################################
# Get any phenotype string required #
#####################################

print "Specify external phenotype file (if any)\n\n";
print "  If you want to use an external phenotype file, type a PLINK-compatible string to use phenotype file.\n";
print "  (NOTE: this should only specify a single phenotype. i.e. don't use --all-pheno)\n\n";
print "  Press return for no input, i.e. to use the phenotype already in the .ped file\n\n";
print "   > ";
$pheno_string = <STDIN>;
chomp $pheno_string ;
print"\n\n";


################################
# Get suffix for this analysis #
################################

print "Enter a suffix for this analysis\n\n";

print "   (e.g. a, b, c... or maybe the particular phenotype chosen)\n\n";
print "   > ";
$analysis_suffix = <STDIN>;
chomp $analysis_suffix;

if ($analysis_suffix)
{
	$analysis_suffix= "_".$analysis_suffix;
}

print "\n";



######################
# Make up file names #
######################

$filename = $prefix.$analysis_suffix;

$logfile = $filename.".log";
$results_file= $filename."_eigen.results";
$final_readme_file = $filename."_eigen.readme";
$results_directory = $filename."_eigen_results";



print "\n";
print "####################################################################\n";
print "# Running PLINK to recode the data into new MAP and PED files      #\n";
print "# for input into the Eigensoft package...                          #\n";
print "####################################################################\n";
print "\n";

################################################
# Run PLINK to recode the ped file, with any   #
# phenotypes that might have been specified    #
################################################


$command =  "plink --noweb --".$species." --allow-no-sex --make-founders";
$command .= " --recode --tab";
$command .= " --file $prefix $pheno_string  --out $filename";

print("$command\n");
system("$command");

#########################################################
# Rename log file to plinklog                           #
#########################################################
$command = "mv $filename".".log $filename".".plinklog";
print("$command\n");
system("$command");


if ($pheno_string ne "")
{
	print "\n\n";
	print "###################################\n";
	print "#    GENERAL CHECK ON THE DATA     #\n";
	print "#                                  #\n";
	print "# Please check at this point that  #\n";
	print "# the phenotype file has been      #\n";
	print "# recognised correctly by PLINK    #\n";
	print "###################################\n";
	print "\n\nPress return to continue (q to quit):    ";
	$ans = <STDIN>;
 	chomp $ans;
}


###################################################
# Open file for results (ie list of programs run  #
###################################################

open(RESULTS,">$results_file") or die "Can't open file $results_file";

print RESULTS "#########################################\n";
print RESULTS "   Eigensoft results file for $filename\n";
print RESULTS "#########################################\n\n";
print RESULTS scalar localtime;
print RESULTS "\n\n";



######################################
# Create the results folder          #
######################################

while (-e $results_directory)
{
	$directory_suffix = $directory_suffix + 1;
	$results_directory = $filename."_eigen_results"."_".$directory_suffix;
}

$command = "mkdir $results_directory";

print("$command\n");
system("$command");

print RESULTS "Creating new directory for eigensoft results $results_directory\n";
print RESULTS " $command\n\n";



################################################################
# Eigenstrat analysis                                          #
################################################################

#########################################################
# Convertf - to convert PLINK files to Eigenstrat files #
#########################################################
$command = "perl /home/genetics/scripts/runconvertf_".$species.".pl -f $filename";
print("$command\n");
system("$command");

print RESULTS "Run convertf_"."$species to convert PLINK files to Eigenstrat files\n";
print RESULTS " $command\n\n";
	
print "\n\n";
print "#########################################################################\n";
print "# Completed running convertf to convert PLINK files to Eigenstrat files #\n";
print "#########################################################################\n";
print "\n\n\n";



#########################################################
# Smartpca - to carry out PCA analysis                  #
#########################################################

if ($QTL eq "no") 
{
	$command = "perl /home/genetics/scripts/runsmartpca_".$species.".pl -f $filename"
}


if ($QTL eq "yes") 
{
	$command = "perl /home/genetics/scripts/runsmartpcaQTL_".$species.".pl -f $filename";
}

print("$command\n");
system("$command");

print RESULTS "Run SMARTPCA to carry out Principal Components Analysis\n";
print RESULTS " $command\n\n";
	
print "\n\n\n";
print "##########################################################################\n";
print "# Completed running SMARTPCA to carry out Principal Components Analysis  #\n";
print "##########################################################################\n";
print "\n\n\n";



#########################################################
# Eigenstrat - to account for structure                 #
#########################################################

if ($QTL eq "no") 
{
	$command = "perl /home/genetics/scripts/runeigenstrat_".$species.".pl -f $filename"
}


if ($QTL eq "yes") 
{
	$command = "perl /home/genetics/scripts/runeigenstratQTL_".$species.".pl -f $filename";
}

print("$command\n");
system("$command");

print RESULTS "Ran EIGENSTRAT to produce CHISQs adjusted for stratification\n";
print RESULTS " $command\n\n";

#########################################################
# Rename log file to eigenlog                           #
#########################################################
$command = "mv $filename".".log $filename".".eigenlog";
print("$command\n");
system("$command");

print RESULTS "Renamed Eigenstrat .log file to .eigenlog\n";
print RESULTS " $command\n\n";
	
print "\n\n\n";
print "##############################################################################\n";
print "# Completed running EIGENSTRAT to produce CHISQs adjusted for stratification #\n";
print "##############################################################################\n";
print "\n\n\n";


#################################################################
# Write the final README file to tell user which files are what #
#################################################################


open(README,">$final_readme_file") or die "Can't open file $final_readme_file";

print README "###############################\n";
print README "# RUN EIGENSOFT_ALL OUTPUT    #\n";
print README "#                             #\n";
print README "# This lists the output files #\n";
print README "# and what to do with them    #\n";
print README "###############################\n\n";


#Eigenstrat #
print README "The following files have the results of running Eigensoft (smartpca and eigenstrat) on the data.\n";
print README "(These programs are NOT part of the PLINK package but use the same files)\n\n";

print README "This file shows a plot of the population structure, with the two main Principal Components on the axes:\n";
print README "(This should look very similar to the MDS plot from PLINK)\n\n";
print README " > ".$filename."_PCA_plot.pdf\n\n\n";


print README "This file has information from the PCA analysis, including statistical measures of the significance of the stratification.\n\n";

print README " > ".$filename.".smartpca_log\n\n\n";


print README "This file has CHISQ association data, uncorrected and also corrected for stratification.\n\n";
print README "(Use the Excel sheet Eigenstrat_plot to plot these two columns of data)\n\n";

print README " > ".$filename.".eigenstrat.out\n\n\n";


print README "-----------------------------------------------------------------------------------------\n";

print README "\nIf there any problems with this program, please see Mike.\n\n";


#######################################
# Copy files  into the results folder #
#######################################
print "Copying all relevant files into $results_directory\n\n";
print RESULTS "Copying all files into $results_directory\n\n";

&move_to_results_folder (".eigenlog");
&move_to_results_folder (".smartpca_log");
&move_to_results_folder (".plinklog");
&move_to_results_folder ("_PCA_plot.pdf");
&move_to_results_folder (".eigenstrat.out");
&move_to_results_folder (".eval");
&move_to_results_folder (".evec");
&move_to_results_folder (".ind");
&move_to_results_folder (".pheno");
&move_to_results_folder ("_eigen.readme");
&move_to_results_folder (".pca");

#################################
# Finally move the results file #
#################################
&move_to_results_folder ("_eigen.results");

###########################################
# Remove some unwanted intermediate files #
###########################################

print "\n\nRemoving some unwanted files\n\n";
print RESULTS "Removing some unwanted files\n\n";


#remove eigenstrat files
&delete_unwanted (".eigenstratgeno");
&delete_unwanted (".par");

&delete_unwanted ("_PCA_plot.ps");
&delete_unwanted ("_PCA_plot.xtxt");
&delete_unwanted (".snp");

&delete_unwanted (".nof");
&delete_unwanted (".hh");
&delete_unwanted (".nosex");
&delete_unwanted (".map");
&delete_unwanted (".ped");



###################################################################################
#                  End message to say the program has finished                    #
###################################################################################

print"\n\n\n";
print "##############################################\n";
print "##############################################\n";
print "##           All analyses complete          ##\n";
print "##############################################\n";
print "##############################################\n";
print "\n\nCheck the following files:\n\n";
print "(These have been moved into a folder called $results_file)\n\n";
print "    1.  $filename"."_eigen.readme   \t-  for how to use the output files\n";
print "    2.  $filename"."_eigen.logall   \t-  for a log of all the analyses\n";
print "    3.  $filename"."_eigen.results  \t-  for more information.\n\n";
print "    4.  $filename".".eigenstrat.out  \t-  for adjusted chis-squared values.\n\n";

if ($QTL eq "no")
{
	print "    Data was analysed as a binary trait.\n\n\n";
}

if ($QTL eq "yes")
{
	print "    Data was analysed as a quantitative trait.\n\n\n";
}


close README;

exit;
###################################################################################################################################################


############################################
#                                          #
# Subroutine to add logfile to logfile_all #
#                                          #
############################################

sub add_to_logfile
{


# open the logfile #
open(LOG,$logfile) or die "Can't open file $logfile";

$errors_found = "false";


while ($single_line = <LOG>) 
{

	chomp $single_line;

	# Check to see if the the line has an ERROR
	$error_pos= index($single_line,'ERROR');
	


	if ($error_pos > -1)
	{
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		print "+                                                     +\n";
		print "+        The Logfile contains an ERROR.               +\n";
		print "+                                                     +\n";
		print "+       Please check and fix the problem              +\n";
		print "+                                                     +\n";
		print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
		
		print LOG "ERROR found in PLINK logfile. Program terminated\n";
		$ans = <STDIN>;
		exit;
		
	}

	# Add line to overall logfile
	print LOGALL "$single_line\n";


} # end of while loop

if ($errors_found eq "true"){print RESULTS "No errors found in PLINK logfile\n\n"}

close LOG;

}

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


