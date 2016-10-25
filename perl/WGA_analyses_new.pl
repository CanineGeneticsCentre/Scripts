#!/usr/bin/perl

use strict;
use Getopt::Std ;
use File::Basename ;

####################
# Define variables #
####################

my $filename			= "";

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
my $pheno_file			= "none";
my $pheno_string		= "";
my $parameters_ok		= "false";
my $choice_ok			= "false";
my $analysis_scope		= "brief";

my $do_assoc			= "y";
my $do_assoc_mperm		= "y";
my $do_eigenstrat		= "y";
my $do_hap_assoc_mperm		= "n";
my $do_models			= "n";
my $do_mds			= "y";
my $do_qq_plot			= "y";
my $do_homozygosity		= "n";
my $do_plink_clustering		= "n";
my $do_which_samples		= "y";
my $do_missingness		= "y";
my $do_ibs_test			= "y";

my $error_pos			= 0;
my $haplotype_count		= 0;

my $mperm_value			= 10000;
my $geno_value			= 0.1;
my $maf_value			= 0.05;
my $mind_value			= 0.1;
my $start_haplotype		= 1;
my $end_haplotype		= 8;
my $directory_suffix		= 0;


print "\n";
print "###############################################\n";
print "#               WGA analyses v29              #\n";
print "#                                             #\n";
print "#     This program runs all the standard      #\n";
print "#        PLINK analyses on your data.         #\n";
print "#                                             #\n";
print "#  It also runs Eigenstrat, and also produces #\n";
print "#  the files for the Homozygosity Mapping     #\n";
print "#  spreadsheet.                               #\n";
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
$filename = <STDIN>;
chomp $filename;
print "\n";

########################
# Check if file exists #
########################
if (! -e "$filename.ped")
{ 
	print "File $filename.ped does not exist\n\n";
	exit;
} 
if (! -e "$filename.map")
{ 
	print "File $filename.map does not exist\n\n";
	exit;
} 



################################
# Get suffix for this analysis #
################################

print "Suffix for this analysis (e.g. a, b, c...):      ";
$analysis_suffix = <STDIN>;
chomp $analysis_suffix;
print "\n";

###########################################
# Copy map and ped file with new suffix #
###########################################

$command = "cp  $filename".".map $filename"."_".$analysis_suffix.".map";
system("$command");
$command = "cp  $filename".".ped $filename"."_".$analysis_suffix.".ped";
system("$command");

$filename = $filename."_".$analysis_suffix;

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


######################
# Make up file names #
######################

$logfile_all = $filename.".logall";
$logfile = $filename.".log";
$results_file= $filename.".results";
$final_readme_file = $filename.".readme";
$results_directory = $filename."_results";



##################################################
# Ask about Analysis scope (brief or full works) #
##################################################

print "\n\n";
print "####################################################\n";
print "# What scope of analysis do you want to carry out? #\n";
print "####################################################\n\n";
print "<1> Basic - This runs assoc, mperm, QQ-plot, mds, Eigenstrat\n";
print "<2> Full  - This runs the whole works (and can take a long time)\n";
print "<3> Pick and choose - This allows you to choose which analyses you do\n\n";

$ans = <STDIN>;
chomp $ans;
print"\n\n";

if ($ans == 1){$analysis_scope = "basic"}
if ($ans == 2){$analysis_scope = "full"}
if ($ans == 3){$analysis_scope = "choose"}

##########################
# Analysis scope = basic #
##########################
if ($analysis_scope eq "basic")
{
 $do_assoc		= "y";
 $do_assoc_mperm	= "y";
 $do_eigenstrat		= "y";
 $do_hap_assoc_mperm	= "n";
 $do_models		= "n";
 $do_mds		= "y";
 $do_qq_plot		= "y";
 $do_homozygosity	= "n";
 $do_plink_clustering	= "n";
 $do_which_samples	= "y";
 $do_missingness	= "y";
 $do_ibs_test		= "y";
}

##########################
# Analysis scope = full  #
##########################
if ($analysis_scope eq "full")
{
 $do_assoc		= "y";
 $do_assoc_mperm	= "y";
 $do_eigenstrat		= "y";
 $do_hap_assoc_mperm	= "y";
 $do_models		= "y";
 $do_mds		= "y";
 $do_qq_plot		= "y";
 $do_homozygosity	= "y";
 $do_plink_clustering	= "y";
 $do_which_samples	= "y";
 $do_missingness	= "y";
 $do_ibs_test		= "y";
}



####################################################
# If choose then ask which things you want to run  #
####################################################
if ($analysis_scope eq "choose"){

 while ($choice_ok eq "false")
{
	print "\n\n";
	print "##########################################\n";
	print "# This is the list of possible analyses. #\n";
	print "##########################################\n\n";

	print "Association:\n\n";

	print "    Association analysis:            			\t$do_assoc\n";
	print "    Association mperm analysis:      			\t$do_assoc_mperm\n";
	print "    Association mperm on all models (dom, rec etc):      \t$do_models\n";
	print "    Haplotype mperm association   			\t$do_hap_assoc_mperm\n\n\n";

	print "Stratification:\n\n";

	print "    IBS test in PLINK					\t$do_ibs_test\n";
	print "    Eigenstrat analysis for stratification:              \t$do_eigenstrat\n";
	print "    Clustering in PLINK for stratification:              \t$do_plink_clustering\n";
	print "    PLINK output for MDS plot:             		\t$do_mds\n";
	print "    PLINK output for QQ plot:             		\t$do_qq_plot\n\n\n";

	print "Other analyses:\n\n";

	print "    Make files for Excel Homozygosity Mapping:		\t$do_homozygosity\n";
	print "    PLINK missingness analysis:				\t$do_missingness\n";



	print "\nWould you like to change these? (y/n)    ";

	$ans = <STDIN>;
	chomp $ans;
	$ans = lc $ans;

	if ($ans eq "y")
	{
		$choice_ok = "false";
	
		print "\n\n";
		print "################################\n";
		print "# Association analysis options #\n";
		print "################################\n\n";

		print "Run association analysis?        	 \t(press return to leave at $do_assoc) ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_assoc = "y"};
		if ($ans eq "n"){$do_assoc = "n"};

		print "Run association mperm analysis?  	 \t(press return to leave at $do_assoc_mperm)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_assoc_mperm = "y"};
		if ($ans eq "n"){$do_assoc_mperm = "n"};

		print "Run mperm with all models (dom, rec etc)?  \t(press return to leave at $do_models)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_models = "y"};
		if ($ans eq "n"){$do_models = "n"};

		print "Run association analysis on haplotypes?    \t(press return to leave at $do_hap_assoc_mperm)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_hap_assoc_mperm = "y"};
		if ($ans eq "n"){$do_hap_assoc_mperm = "n"};


		print "\n\n";
		print "###################################\n";
		print "# Stratification analysis options #\n";
		print "###################################\n\n";

		print "Do IBS test in PLINK?             	\t(press return to leave at $do_ibs_test)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_ibs_test = "y"};
		if ($ans eq "n"){$do_ibs_test = "n"};

		print "Run Eigenstrat analysis?         	\t(press return to leave at $do_eigenstrat)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_eigenstrat = "y"};
		if ($ans eq "n"){$do_eigenstrat = "n"};


		print "Run PLINK with auto-clustering?  	\t(press return to leave at $do_plink_clustering)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_plink_clustering = "y"};
		if ($ans eq "n"){$do_plink_clustering = "n"};


		print "Run PLINK for MDS plot?   		\t(press return to leave at $do_mds)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_mds = "y"};
		if ($ans eq "n"){$do_mds = "n"};


		print "Run PLINK for QQ plot?    		\t(press return to leave at $do_qq_plot)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_qq_plot = "y"};
		if ($ans eq "n"){$do_qq_plot = "n"};


		print "\n\n";
		print "##########################\n";
		print "# Other analysis options #\n";
		print "##########################\n\n";

		print "Make files for 'Homozygosity Mapping'?	\t(press return to leave at $do_homozygosity)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_homozygosity = "y"};
		if ($ans eq "n"){$do_homozygosity = "n"};

		print "Do PLINK missingness analysis?		\t(press return to leave at $do_missingness)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_missingness = "y"};
		if ($ans eq "n"){$do_missingness= "n"};



	}
	else
	{
		$choice_ok = "true";
	}

}
} # end of choose


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
	print "    --mperm $mperm_value\n";

	if ($do_hap_assoc_mperm eq "y")
	{
		print "    haplotypes $start_haplotype-$end_haplotype\n";
	}
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

		print "    --mperm value (press return to leave at $mperm_value):   ";
		$ans = <STDIN>;
		chomp $ans;
		if ($ans ne ""){$mperm_value = $ans;}


		if ($do_hap_assoc_mperm eq "y")
		{
			print "    Smallest haplotype value (press return to leave at $start_haplotype):   ";
			$ans = <STDIN>;
			chomp $ans;
			if ($ans ne ""){$start_haplotype = $ans;}

			print "    Largest haplotype value (press return to leave at $end_haplotype):    ";
			$ans = <STDIN>;
			chomp $ans;
			if ($ans ne ""){$end_haplotype = $ans;}
		}

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

######################################################################
# Warn that if you have external phenotypes you can't do missingness #
######################################################################
if (($do_missingness eq "y") && ($pheno_string ne ""))
{
	print "If you have no phenotypes in your ped file, you can't do missingness on it.\n\n";

	print "Omit missingness analysis for now? (Y/N)    ";
	$ans = <STDIN>;
	if (($ans eq "y") || ($ans eq "Y"))
	{
		$do_missingness = "n";
	}
	if (($ans ne "y") && ($ans ne "Y"))
	{
		exit;
	}

}


########################
# Open overall logfile #
########################

open(LOGALL,">$logfile_all") or die "Can't open file $logfile_all";


########################
# Open simple logfile  #
########################

open(RESULTS,">$results_file") or die "Can't open file $results_file";
print RESULTS "#########################################\n";
print RESULTS "   Results file for $filename\n";
print RESULTS "#########################################\n\n";
print RESULTS scalar localtime;
print RESULTS "\n\n";

print RESULTS "Parameters in place:\n\n";
print RESULTS "--maf $maf_value\n";
print RESULTS "--geno $geno_value\n";
print RESULTS "--mind $mind_value\n";
print RESULTS "--mperm $mperm_value\n\n";
print RESULTS "haplotypes from $start_haplotype to $end_haplotype\n";
print RESULTS "\n";
print "\n";


######################################
# Create the results folder          #
######################################

while (-e $results_directory)
{
	$directory_suffix = $directory_suffix + 1;
	$results_directory = $filename."_results"."_".$directory_suffix;
}

$command = "mkdir $results_directory";

print("$command\n");
system("$command");

print RESULTS "Creating new directory for results $results_directory\n";
print RESULTS " $command\n\n";


#################################################################
#                                                               #
#                          Now run PLINK                        #
#                                                               #
#################################################################


###########################################################################
# General check that the data is working and that the filters are correct #
###########################################################################

print "\n\n";
print "############################################\n";
print "#    RUNNING GENERAL CHECK ON THE DATA     #\n";
print "#                                          #\n";
print "# Please check when this has run that the  #\n";
print "# values you used for geno, maf and mind   #\n";
print "# have not filtered off too many SNPs or   #\n";
print "# samples.                                 #\n";
print "############################################\n\n";


$command = "plink --noweb --".$species." --allow-no-sex --maf $maf_value --geno $geno_value --mind $mind_value --file $filename $pheno_string --make-founders --nonfounders --out $filename";
print("$command\n");
system("$command");


print "\n\n\n";
print "###################################\n";
print "#    GENERAL CHECK ON THE DATA     #\n";
print "#                                  #\n";
print "# Please check at this point that  #\n";
print "# the values you used for geno and #\n";
print "# maf have not filtered off too    #\n";
print "# many SNPs or samples.            #\n";
print "###################################\n";
print "\n\nPress return to continue.\n";
$ans = <STDIN>;


#########################################################
# Check on the missingness of the data                  #
#########################################################

if ($do_missingness eq "y")
{

 $command = "plink --noweb --".$species." --allow-no-sex --file $filename $pheno_string --nonfounders --make-founders --test-missing --out $filename";
 print("$command\n");
 system("$command");

 # Add logfile to logfile_all
 &add_to_logfile;

 print RESULTS "Run --test-missing on the data\n";
 print RESULTS " $command\n\n";

 &move_to_results_folder (".missing");

 print "\n\n\n";
 print "##########################################\n";
 print "# Checked on the missingness of the data #\n";
 print "##########################################\n";
 print "\n\n\n";

}


#########################################################
# First convert the .map and .ped files in binary files #
#########################################################


$command = "plink --noweb --".$species." --allow-no-sex --file $filename $pheno_string --make-founders --nonfounders --make-bed  --maf $maf_value --geno $geno_value --mind $mind_value --out $filename";
print("$command\n");
system("$command");

# Add logfile to logfile_all
&add_to_logfile;

print RESULTS "Convert map and ped files into binary file (also filter with MAF, GENO and MIND)\n";
print RESULTS " $command\n\n";

print "\n\n\n";
print "####################################\n";
print "# Converted files to binary files  #\n";
print "####################################\n";
print "\n\n\n";



################################################################################
# Run 'which_samples_were_used' to make a list of samples in the .ped file     #
################################################################################
$command = "perl /home/genetics/scripts/which_samples_were_used.pl -f $filename";
print("$command\n");
system("$command");

print RESULTS "Run which_samples_were_used to keep a record of which samples were in the .ped file\n";
print RESULTS " $command\n\n";

&move_to_results_folder (".samples_used.txt");

	
print "\n\n\n";
print "########################################################################################\n";
print "# Ran which_samples_were_used to keep a record of which samples were in the .ped file  #\n";
print "########################################################################################\n";
print "\n\n\n";



#########################################################
# Run standard association analysis                     #
#########################################################

if ($do_assoc eq "y")
{

$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --assoc --out $filename";
print("$command\n");
system("$command");

# Add logfile to logfile_all
&add_to_logfile;

print RESULTS "Run standard association analysis\n";
print RESULTS " $command\n\n";

&move_to_results_folder (".assoc");

print "\n\n\n";
print "#####################################\n";
print "# Ran standard association analysis #\n";
print "#####################################\n";
print "\n\n\n";

}

#########################################################
# Produce QQ-plot output                                #
#########################################################

if ($do_qq_plot eq "y")
{

 $command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --assoc --adjust --qq-plot --log10 --out $filename";
 print("$command\n");
 system("$command");

 # Add logfile to logfile_all
 &add_to_logfile;

 print RESULTS "Run PLINK to get data for QQ plot\n";
 print RESULTS " $command\n\n";

 &move_to_results_folder (".assoc");
 &move_to_results_folder (".assoc.adjusted");

 print "\n\n\n";
 print "##########################################\n";
 print "# Produced file with data for QQ plot #\n";
 print "##########################################\n";
 print "\n\n\n";

}

################################################################################
# Run plink2homozygosity to create files for Homozygosity Mapping spreadsheet  #
################################################################################

if ($do_homozygosity eq "y")
{

 $command = "perl /home/genetics/scripts/plink2homozygosity.pl -f $filename";
 print("$command\n");
 system("$command");

 print RESULTS "Run PLINK2HOMOZYGOSITY to create files for Homozygosity Mapping spreadsheet\n";
 print RESULTS " $command\n\n";
	
 &move_to_results_folder ("_affected.txt");
 &move_to_results_folder ("_normal.txt");

 print "\n\n\n";
 print "################################################################################\n";
 print "# Ran plink2homozygosity to create files for Homozygosity Mapping spreadsheet  #\n";
 print "################################################################################\n";
 print "\n\n\n";

}


#########################################################
# Produce MDS data to look at clustering                #
# (doesn't do it for any particular phenotype)          #
#########################################################

if ($do_mds eq "y")
{

 $command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --mds-plot 2 --out $filename";
 print("$command\n");
 system("$command");

 # Add logfile to logfile_all
 &add_to_logfile;

 print RESULTS "Run MDS (multi-dimensional-scaling) analysis corrected for multiple testing to look at clustering\n";
 print RESULTS " $command\n\n";

 &move_to_results_folder (".mds");

 print "\n\n\n";
 print "##########################################\n";
 print "# Ran MDS analysis to look at clustering #\n";
 print "##########################################\n";
 print "\n\n\n";

}


############################################################
# Run association corrected for multiple testing  (MPERM)  #
############################################################

if ($do_assoc_mperm eq "y")
{

 $command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --assoc --mperm $mperm_value --out $filename";
 print("$command\n");
 system("$command");

 # Add logfile to logfile_all
 &add_to_logfile;

 print RESULTS "Run association analysis corrected for multiple testing\n";
 print RESULTS " $command\n\n";

 &move_to_results_folder (".assoc.mperm");

 print "\n\n\n";
 print "################################################################################\n";
 print "# Ran standard association analysis corrected for multiple testing using mperm #\n";
 print "################################################################################\n";
 print "\n\n\n";

}


############################################################
# Only include all models if analysis scope is set to full #
############################################################

if ($do_models eq "y")
{

#########################################################
# Run association with all models (dom,rec etc)         #
#########################################################


$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --model --cell 0 --out $filename";
print("$command\n");
system("$command");

# Add logfile to logfile_all
&add_to_logfile;

print RESULTS "Run PLINK to look at different models (dom, rec etc) \n";
print RESULTS " $command\n\n";

&move_to_results_folder (".model");

print "\n\n\n";
print "#############################################\n";
print "# Looked at different models (dom, rec etc) #\n";
print "#############################################\n";
print "\n\n\n";



#########################################################
# Run association mperm with best model                 #
#########################################################


$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --model --cell 0 --mperm $mperm_value --out $filename";
print("$command\n");
system("$command");

# Add logfile to logfile_all
&add_to_logfile;

print RESULTS "Run PLINK to look at best model with mperm \n";
print RESULTS " $command\n\n";

&move_to_results_folder (".model");
&move_to_results_folder (".model.best.mperm");

print "\n\n\n";
print "#############################################\n";
print "# Looked at best model with mperm           #\n";
print "#############################################\n";
print "\n\n\n";





#########################################################
# Run mperm association with the dominant model         #
#########################################################


$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --model-dom --cell 0 --mperm $mperm_value --out $filename";
print("$command\n");
system("$command");

# Add logfile to logfile_all
&add_to_logfile;

print RESULTS "Run PLINK to look at dominant model with mperm \n";
print RESULTS " $command\n\n";

&move_to_results_folder (".model");
&move_to_results_folder (".model.dom.mperm");

print "\n\n\n";
print "#############################################\n";
print "# Looked at dominant model with mperm       #\n";
print "#############################################\n";
print "\n\n\n";




#########################################################
# Run mperm association with the recessive model        #
#########################################################


$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --model-rec --cell 0 --mperm $mperm_value --out $filename";
print("$command\n");
system("$command");

# Add logfile to logfile_all
&add_to_logfile;

print RESULTS "Run PLINK to look at recessive model with mperm \n";
print RESULTS " $command\n\n";

&move_to_results_folder (".model");
&move_to_results_folder (".model.rec.mperm");

print "\n\n\n";
print "#############################################\n";
print "# Looked at recessive model with mperm      #\n";
print "#############################################\n";
print "\n\n\n";




#########################################################
# Run mperm association with the genotypic model        #
#########################################################


$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --model-gen --cell 0 --mperm $mperm_value --out $filename";
print("$command\n");
system("$command");

# Add logfile to logfile_all
&add_to_logfile;

print RESULTS "Run PLINK to look at genotypic model with mperm \n";
print RESULTS " $command\n\n";

&move_to_results_folder (".model");
&move_to_results_folder (".model.gen.mperm");

print "\n\n\n";
print "#############################################\n";
print "# Looked at genotypic model with mperm      #\n";
print "#############################################\n";
print "\n\n\n";





#########################################################
# Run mperm association with the trend model            #
#########################################################


$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --model-trend --cell 0 --mperm $mperm_value --out $filename";
print("$command\n");
system("$command");

# Add logfile to logfile_all
&add_to_logfile;

print RESULTS "Run PLINK to look at Cochran Armitage trend model with mperm \n";
print RESULTS " $command\n\n";

&move_to_results_folder (".model");
&move_to_results_folder (".model.trend.mperm");

print "\n\n\n";
print "#############################################\n";
print "# Looked at trend model with mperm          #\n";
print "#############################################\n";
print "\n\n\n";


} # do_models = y




#############################################################
# Only do plink clustering if analysis scope is set to full #
#############################################################

if ($do_plink_clustering eq "y")
{

#########################################################
# Allow PLINK to cluster the samples                    #
#########################################################


$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --cluster --cc --ppc 0.01 --out $filename";

print("$command\n");
system("$command");

# Add logfile to logfile_all
&add_to_logfile;

print RESULTS "Run PLINK to produce clusters with at least one case and one control and ppc=0.01 (read PLINK documentation)...\n";
print RESULTS " $command\n\n";

print "\n\n\n";
print "##########################################\n";
print "# Allowed PLINK to cluster the samples   #\n";
print "##########################################\n";
print "\n\n\n";




#########################################################
# Run association only within these clusters            #
#########################################################


$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --mh --within $filename".".cluster2 --out $filename";

print("$command\n");
system("$command");

# Add logfile to logfile_all
&add_to_logfile;

print RESULTS "Run MH association within the cluster2 file...\n";
print RESULTS " $command\n\n";

&move_to_results_folder (".cmh");


print "\n\n\n";
print "#################################################\n";
print "# Looked for assication within these clusters   #\n";
print "#################################################\n";
print "\n\n\n";


###########################################################
# Run association only within these clusters using mperm  #
###########################################################


$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --mh --mperm $mperm_value --within $filename".".cluster2 --out $filename";

print("$command\n");
system("$command");

# Add logfile to logfile_all
&add_to_logfile;

print RESULTS "Run MH association within the cluster2 file using mperm...\n";
print RESULTS " $command\n\n";

&move_to_results_folder (".cmh");
&move_to_results_folder (".cmh.mperm");
&move_to_results_folder (".cluster2");

print "\n\n\n";
print "############################################################\n";
print "# Looked for assication within these clusters using Mperm  #\n";
print "############################################################\n";
print "\n\n\n";


} # of analysis = full or do_plink_clustering = "y"



############################################################
# Only do haplotypes if analysis scope is set to full      #
############################################################

if ($do_hap_assoc_mperm eq "y")
{



################################################################
# Haplotype analysis - from $start_haplotype to $end_haplotype #
################################################################

if ($start_haplotype > 0 && $start_haplotype < $end_haplotype)
{

    for ($haplotype_count = $start_haplotype;$haplotype_count <= $end_haplotype;$haplotype_count ++)
    {

	#############################################
	# Haplotype analysis  - generate haplotypes #
	#############################################
	
	$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename $pheno_string --hap-window $haplotype_count --hap-impute --out $filename"."_".$haplotype_count;
	print("$command\n");
	system("$command");
	
	# Add logfile to logfile_all
	&add_to_logfile;
	
	print RESULTS "Run Haplotyping analysis on window of $haplotype_count\n";
	print RESULTS " $command\n\n";
	
	print "\n\n\n";
	print "######################################################\n";
	print "# Imputed the haplotypes (window = $haplotype_count) #\n";
	print "######################################################\n";
	print "\n\n\n";



	#######################################################################
	# Haplotype analysis  - make binary files and delete non-binary files #
	#######################################################################
	
	$command = "plink --noweb --".$species." --allow-no-sex --file $filename"."_".$haplotype_count.".impute $pheno_string --make-bed --out $filename"."_".$haplotype_count.".impute";
	print("$command\n");
	system("$command");
	
	# Add logfile to logfile_all
	&add_to_logfile;
	
	print RESULTS "Made binary files for haplotyping on window of $haplotype_count\n";
	print RESULTS " $command\n\n";
	
	print "\n\n\n";
	print "#########################################################\n";
	print "# Made binary files for haplotyping  (window = $haplotype_count) #\n";
	print "#########################################################\n";
	print "\n\n\n";

	##################################
	# Delete the non-binary versions #
	##################################

	&delete_unwanted ("_".$haplotype_count.".impute.ped");
	
	
	#######################################################
	# Haplotype analysis  - association with permutation  #
	#######################################################
	
	
	$command = "plink --noweb --".$species." --allow-no-sex --bfile $filename"."_".$haplotype_count.".impute $pheno_string --assoc --mperm $mperm_value --out $filename"."_".$haplotype_count.".impute";
	print("$command\n");
	system("$command");
	
	# Add logfile to logfile_all
	&add_to_logfile;
	
	print RESULTS "Run Haplotyping association with mperm on window of $haplotype_count\n";
	print RESULTS " $command\n\n";
	

	print "\n\n\n";
	print "#########################################################\n";
	print "# Ran Haplotyping to use in Hap_Mperm_Plot (window = $haplotype_count) #\n";
	print "#########################################################\n";
	print "\n\n\n";

	
	########################################################
	# Convert impute.map and impute.assoc.mperm files to   #
	# a single file which can be plotted in Hap_Mperm_plot #
	########################################################

 	$command = "perl /home/genetics/scripts/hap_mperm_convert.pl -f $filename"."_".$haplotype_count." -s $species -o $filename"."_{".$haplotype_count."}";
 	print("$command\n");
 	system("$command");

	print RESULTS "Converted impute.map and impute.assoc.mperm files to single _hap_mperm.txt file\n";
	print RESULTS " $command\n\n";
	
	&move_to_results_folder ("_".$haplotype_count.".impute.map");
	&move_to_results_folder ("_".$haplotype_count.".impute.assoc.mperm");

	&move_to_results_folder ("_{".$haplotype_count."}");


	print "#########################################################\n";
	print "# Converted output files with hap_mperm_convert  (window = $haplotype_count) #\n";
	print "#########################################################\n";


	##########################
	# Deleted unwanted files #
	##########################
	&delete_unwanted ("_".$haplotype_count.".impute.bed");
	&delete_unwanted ("_".$haplotype_count.".impute.bim");
	&delete_unwanted ("_".$haplotype_count.".impute.fam");


    } # end of $haplotype_count loop

} # end of if start_haplotype > 0


} # end of analysis = full or hap_assoc_mperm = y



if ($do_ibs_test eq "y")
{
 #########################################################
 # IBS TEST in PLINK                                     #
 #########################################################

 $command = "plink --noweb --".$species." --allow-no-sex --file $filename $pheno_string --nonfounders --make-founders --ibs-test --out $filename";
 print("$command\n");
 system("$command");

 print RESULTS "Run ibs-test in PLINK\n";
 print RESULTS " $command\n\n";
	
 print "\n\n\n";
 print "###########################################################\n";
 print "# Ran ibs-test in PLINK                                   #\n";
 print "###########################################################\n";
 print "\n\n\n";

####################################
# Rename log file to ibs-test file #
####################################
 $command = "mv $filename".".log $filename".".ibs-test"; 
 print("$command\n");
 system("$command");

 &move_to_results_folder (".ibs-test");

}


################################################################
# Eigenstrat analysis                                          #
################################################################

if ($do_eigenstrat eq "y")
{

 #########################################################
 # Convertf - to convert PLINK files to Eigenstrat files #
 #########################################################
 $command = "perl /home/genetics/scripts/runconvertf_".$species.".pl -f $filename";
 print("$command\n");
 system("$command");

 print RESULTS "Run convertf to convert PLINK files to Eigenstrat files\n";
 print RESULTS " $command\n\n";
	
 print "\n\n\n";
 print "###########################################################\n";
 print "# Ran convertf to convert PLINK files to Eigenstrat files #\n";
 print "###########################################################\n";
 print "\n\n\n";




 #########################################################
 # Smartpca - to carry out PCA analysis                  #
 #########################################################
 $command = "perl /home/genetics/scripts/runsmartpca_".$species.".pl -f $filename";
 print("$command\n");
 system("$command");

 print RESULTS "Run SMARTPCA to carry out Principal Components Analysis\n";
 print RESULTS " $command\n\n";
	
 print "\n\n\n";
 print "############################################################\n";
 print "# Ran SMARTPCA to carry out Principal Components Analysis  #\n";
 print "############################################################\n";
 print "\n\n\n";



 #########################################################
 # Eigenstrat - to account for structure                 #
 #########################################################
 $command = "perl /home/genetics/scripts/runeigenstrat_".$species.".pl -f $filename";
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

 &move_to_results_folder (".eigenlog");
 &move_to_results_folder ("_plot.pdf");
 &move_to_results_folder (".smartpca_log");
 &move_to_results_folder (".eigenstrat.out");
 &move_to_results_folder (".chisq");

	
 print "\n\n\n";
 print "################################################################\n";
 print "# Ran EIGENSTRAT to produce CHISQs adjusted for stratification #\n";
 print "################################################################\n";
 print "\n\n\n";

} # end of if do eigenstrat = y


close LOGALL;

&move_to_results_folder (".logall");
&move_to_results_folder (".assoc");


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

# Missingness #
print README "The following file has a list of all the command lines used to run these analyses:\n";

print README " > ".$filename.".results\n";

print README "-----------------------------------------------------------------------------------------\n";

# Missingness #
if ($do_missingness eq "y")
{
print README "The following file can be used to look at any missing data\n";
print README "to make sure any 'missingness' is the same in the cases and controls\n\n";

print README " > ".$filename.".missing\n";


print README "Open this file in Excel and sort by the P-value column.\n";
print README "This tells you if there is a significance difference between cases and controls for any SNP.\n\n";
print README "-----------------------------------------------------------------------------------------\n";
}


# Basic association #
if ($do_assoc eq "y")
{
 print README "The following file has the simple association data.\n";
 print README "It can be plotted using the Excel sheet PLINK_plot.\n\n";

 print README " > ".$filename.".assoc\n\n";

 print README "-----------------------------------------------------------------------------------------\n";
}


# MPERM association #
if ($do_assoc_mperm eq "y")
{
 print README "The following file has the association data, taking into account multiple testing.\n";
 print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";

 print README " > ".$filename.".assoc.mperm\n\n";

 print README "-----------------------------------------------------------------------------------------\n";
}


# Clustering #
if ($do_mds eq "y")
{
print README "The following MDS file has the results of a clustering analysis (multi-dimensional scaling).\n";
print README "If you plot this in Excel you can visualise any clusters.\n";
print README "(Do a scatter plot, and plot a separate series for cases and controls).\n\n";

print README " > ".$filename.".mds\n\n";

print README "-----------------------------------------------------------------------------------------\n";
}




if ($do_models eq "y")
{

# Simple association with all models#
print README "The following file has the association data for five different models.\n";
print README "Normal (allelic), Recessive, Dominant, Trend and Genotypic.\n";
print README "It can also be plotted using the Excel sheet PLINK_plot.\n";
print README "(Note you must read itin this file with the special button on the Front Page)\n\n";

print README " > ".$filename.".model\n\n";

print README "-----------------------------------------------------------------------------------------\n";


# Mperm association with the best model#
print README "The following file has the association data for the best of five different models at each SNP.\n";
print README "Normal (allelic), Recessive, Dominant, Trend and Genotypic.\n";
print README "(This also takes into account multiple testing by using mperm)\n";
print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";

print README " > ".$filename.".model.best.mperm\n\n";

print README "-----------------------------------------------------------------------------------------\n";



# Mperm association with the dominant model#
print README "The following file has the association data for the dominant model.\n";
print README "(This also takes into account multiple testing by using mperm)\n";
print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";

print README " > ".$filename.".model.dom.mperm\n\n";

print README "-----------------------------------------------------------------------------------------\n";



# Mperm association with the recessive model#
print README "The following file has the association data for the recessive model.\n";
print README "(This also takes into account multiple testing by using mperm)\n";
print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";

print README " > ".$filename.".model.rec.mperm\n\n";

print README "-----------------------------------------------------------------------------------------\n";


# Mperm association with the genotypic model#
print README "The following file has the association data for the genotypic model.\n";
print README "(This also takes into account multiple testing by using mperm)\n";
print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";

print README " > ".$filename.".model.gen.mperm\n\n";

print README "-----------------------------------------------------------------------------------------\n";


# Mperm association with the trend model#
print README "The following file has the association data for the Cochran Armitage trend model.\n";
print README "(This also takes into account multiple testing by using mperm)\n";
print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";

print README " > ".$filename.".model.trend.mperm\n\n";

print README "-----------------------------------------------------------------------------------------\n";
}

# Stratification #
if ($do_ibs_test eq "y")
{
print README "The following file has the results of using the PLINK --ibs-test option.\n";
print README "This shows how statistically different the cases and control populations are.\n\n";

print README " > ".$filename.".ibs_test\n";

print README "-----------------------------------------------------------------------------------------\n";
}




# Clustering 2 #
if ($do_plink_clustering eq "y")
{
print README "The following files are the results of allowing PLINK to cluster the samples and\n";
print README "then carrying out a MH association only within clusters.\n";
print README "This file containes the cluster data...\n\n";

print README " > ".$filename.".cluster2\n\n";

print README "These files can be plotted in PLINK_plot.\n\n";

print README " > ".$filename.".cmh\n\n";

print README " > ".$filename.".cmh.mperm\n\n";

print README "-----------------------------------------------------------------------------------------\n";
}



# Output for QQ-plot #
if ($do_qq_plot eq "y")
{
print README "The following file can be used to produce a QQ-plot in Excel using QQ_plot.\n\n";

print README " > ".$filename.".assoc.adjusted\n";
print README "-----------------------------------------------------------------------------------------\n";
}


if ($do_which_samples eq "y")
{
# which_samples_were_used #
print README "The following file keeps a list of the samples used in the .ped file:\n\n";

print README " > ".$filename.".samples_used.txt\n";


print README "-----------------------------------------------------------------------------------------\n";
}


# Homozygosity Mapping files #
if ($do_homozygosity eq "y")
{
print README "The following files can be used in the Excel Homozygosity Mapping spreadsheet\n";
print README "Open them in two separate worksheets called Affected and Normal\n";
print README "(This will not use any phenotypes which are in a separate file from the .ped file)\n\n";

print README " > ".$filename."_affected.txt\n";
print README " > ".$filename."_normal.txt\n\n";

print README "-----------------------------------------------------------------------------------------\n";
}


if ($do_hap_assoc_mperm eq "y")
{

# Haplotyping #
print README "The following files have the results of the haplotype analysis with permutation.\n";
print README "(This has been done for window sizes (N) = $start_haplotype-$end_haplotype)\n";
print README "Use the Excel sheet Hap_Mperm_plot to visualise these results.\n\n";

print README " > ".$filename."_N.impute.map\n";
print README " > ".$filename."_N.impute.assoc.mperm\n\n";

if ($species eq "dog")
{
print README "This file has been converted to work directly in Hap_Mperm_plot with out having to read it in with the button.\n\n";
print README " > ".$filename."_N_hap_mperm.txt\n\n"
}

print README "-----------------------------------------------------------------------------------------\n";

}




#Eigenstrat #
if ($do_eigenstrat eq "y")
{
print README "The following files have the results of running Eigensoft (smartpca and eigenstrat) on the data.\n";
print README "(These programs are NOT part of the PLINK package but use the same files)\n";
print README "(NOTE: This will not use any phenotypes which are in a separate file from the .ped file)\n\n";

print README "This file shows a plot of the population structure, with the two main Principal Components on the axes:\n";
print README "(This should look very similar to the MDS plot from PLINK)\n\n";
print README " > ".$filename."_plot.pdf\n\n";

print README "This file has information from the Eigensoft analysis, including statistical measures of the significance of the stratification.\n\n";

print README " > ".$filename.".eigenlog\n\n";

print README "This file has CHISQ association data, uncorrected and also corrected for stratification.\n";
print README "(Use the Excel sheet Eigenstrat_plot to plot these two columns of data)\n\n";

print README " > ".$filename.".eigenstrat_out OR ".$filename.".chisq\n\n";

print README "-----------------------------------------------------------------------------------------\n";
}


print README "\nIf there any problems with this program, please see Mike Boursnell.\n\n";


###################################################
# Copy readme and results into the results folder #
###################################################

&move_to_results_folder (".readme");
&move_to_results_folder (".results");


###########################################
# Remove some unwanted intermediate files #
###########################################

print "\n\n";
print "################################\n\n";
print "# Removing some unwanted files #\n";
print "################################\n\n";
print RESULTS "Removing some unwanted files...\n\n";




if ($do_mds eq "y" ||$do_plink_clustering eq "y")
{
	&delete_unwanted (".cluster0");
	&delete_unwanted (".cluster3");
}

#remove plink files
&delete_unwanted (".log");
&delete_unwanted (".ind");
&delete_unwanted (".map");
&delete_unwanted (".ped");
&delete_unwanted (".hh");
&delete_unwanted (".irem");

#remove all PLINK files of these types
delete_all_of_this_type (".nosex");
delete_all_of_this_type (".nof");

#####################################################################################
# Delete the binary versions of the files (can easily regenerate these if required) #
#####################################################################################
  &delete_unwanted (".fam");
  &delete_unwanted (".bim");
  &delete_unwanted (".bed");

##########################
# Remove haplotype files #
##########################

if ($do_hap_assoc_mperm eq "y")
{
	for ($haplotype_count = $start_haplotype;$haplotype_count <= $end_haplotype;$haplotype_count ++)
	{
	#&delete_unwanted ("_".$haplotype_count.".impute.ped");
	#&delete_unwanted ("_".$haplotype_count.".impute.map");
	#&delete_unwanted ("_".$haplotype_count.".impute.fam");
	&delete_unwanted ("_".$haplotype_count.".impute.nosex");
	&delete_unwanted ("_".$haplotype_count.".impute.assoc");
	&delete_unwanted ("_".$haplotype_count.".impute.log");
	&delete_unwanted ("_".$haplotype_count.".log");
	&delete_unwanted ("_".$haplotype_count.".hh");
	}
}

###########################
# Remove eigenstrat files #
###########################
if ($do_eigenstrat eq "y")
{
	&delete_unwanted (".eigenstratgeno");
	&delete_unwanted (".par");
	&delete_unwanted (".pca");
	&delete_unwanted (".eval");
	&delete_unwanted (".evec");
	&delete_unwanted (".pheno");
	&delete_unwanted ("_plot.ps");
	&delete_unwanted ("_plot.xtxt");
	&delete_unwanted (".snp");
}

#remove plink2homozygosityfiles
if ($do_homozygosity eq "y") 
{
	&delete_unwanted (".p2h_log");
}

#remove gmon.out file
&delete_file ("gmon.out");

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
print "(These have been moved into a folder called $filename"."_results)\n\n";
print "    1.  $filename".".readme   \t-  for how to use the output files\n";
print "    2.  $filename".".logall   \t-  for a log of all the analyses\n";
print "    3.  $filename".".results  \t-  for more information.\n\n";

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

print LOGALL "------------------------------------------------------------------------------------------------\n";
print LOGALL "\n\n $command\n\n";

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

	if (-e "$filename".$suffix)
	{ 
 	
		$command = "rm  $filename".$suffix;
		print("$command\n");
		system("$command");
		print RESULTS ("$command\n");

	}
	else
	{
		print ">>> The file $filename".$suffix." doesn't exist <<<\n";
	}
}


#############################################
#                                           #
# Subroutine to remove all files of a type  #
#                                           #
#############################################

sub delete_all_of_this_type
{
	my $suffix = "";

	$suffix = $_[0];


	$command = "rm  *".$suffix;
	print("$command\n");
	system("$command");
	print RESULTS ("$command\n");


}


#############################################
#                                           #
# Subroutine to remove unwanted files       #
#                                           #
#############################################

sub delete_file
{
	my $thisfile = "";

	$thisfile = $_[0];

	if (-e "$thisfile")
	{ 
 	
		$command = "rm  $thisfile";;
		print("$command\n");
		system("$command");
		print RESULTS ("$command\n");

	}
	else
	{
		print ">>> The file $thisfile doesn't exist <<<\n";
	}

}


