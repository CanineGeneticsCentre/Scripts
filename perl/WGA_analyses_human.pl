#!/usr/bin/perl

#########################################################################
#									                                    #      
#	WGA_analyses           						                        #     
#									                                    #
#	Runs a series of programs to carry out a complete analysis of       #
#   Genome Wide Association Scan SNP genotyping data                    #
#									                                    #
#########################################################################

#############################
# Mike Boursnell            #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# mike.boursnell@aht.org.uk #
#############################
use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;

# Constants
my $version				= 35;

####################
# Define variables #
####################

my $pheno_file			= "none";
my $pheno_string		= "";
my $parameters_ok		= "false";
my $choice_ok			= "false";
my $analysis_scope		= "brief";
my $from				= 'WGA_analyses@samba64.aht.org.uk';

# Default choices for analysis
my $do_assoc			= "y";
my $do_assoc_mperm		= "y"; 
my $do_eigenstrat		= "n";
my $do_emmax			= "y";
my $do_hap_assoc		= "n";
my $do_hap_assoc_mperm	= "n";
my $do_models			= "n";
my $do_mds				= "y";
my $do_qq_plot			= "y";
my $do_homozygosity		= "n";
my $do_plink_clustering	= "n";
my $do_which_samples	= "y";
my $do_missingness		= "y";
my $do_ibs_test			= "y";

my $filename			= "";
my $ans					= "";
my $species_code		= "";
my $species				= "";
my $command				= "";
my $single_line			= "";
my $errors_found		= "";
my $results_directory	= "";
my $suffix				= "";
my $analysis_suffix		= "";
my $email_address		= "";


my $error_pos			= 0;
my $haplotype_count		= 0;

my $mperm_value			= 10000;
my $geno_value			= 0.1;
my $maf_value			= 0.05;
my $mind_value			= 0.1;
my $start_haplotype		= 1;
my $end_haplotype		= 8;
my $directory_suffix	= 0;

#File names
my $command_log			= "";
my $logfile				= "";
my $logfile_all			= "";
my $results_file		= "";
my $final_readme_file	= "";
my $prefix				= "";


print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "       WGA analyses       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This program can run all the standard PLINK analyses on your GWAS data\n";
print "  - It can also run analyses by other programs such as EMMAX, Eigenstrat etc.\n\n";

print color 'bold cyan';

print "    (The input files are PLINK format ped and map files)\n\n";

print color 'reset';


############################
# Get the file name prefix #
############################

print "Prefix of ped and map files:      ";
$prefix = <STDIN>;
chomp $prefix;
print "\n";




########################
# Check if file exists #
########################
if (! -e "$prefix.ped")
{ 
	print "File $prefix.ped does not exist\n\n";
	exit;
} 
if (! -e "$prefix.map")
{ 
	print "File $prefix.map does not exist\n\n";
	exit;
} 



################################
# Get suffix for this analysis #
################################

print "Suffix for this analysis (e.g. a, b, c...):      ";
$analysis_suffix = <STDIN>;
chomp $analysis_suffix;
print "\n";


######################################
# E-mail address for notifications   #
######################################

print "\nEnter your email address if you want to be notified when analyses are complete:    ";
print "(Otherwise press 'return')\n\n";
print ">  ";

$email_address = <STDIN>;
chomp $email_address;

# Some short cuts
if ($email_address eq "m"){$email_address = 'mike.boursnell@aht.org.uk';}
if ($email_address eq "o"){$email_address = 'oliver.forman@aht.org.uk';}


###########################################
# Copy map and ped file with new suffix   #
###########################################
$filename = $prefix."_".$analysis_suffix;

&run_unix_command("cp  $prefix".".map $filename".".map");
&run_unix_command("cp  $prefix".".ped $filename".".ped");


#############################
# Ask if it is dog or horse #
#############################
print color 'yellow';
print "\n\n####################################\n";
print "# Which species are you analysing? #\n";
print "####################################\n\n";
print color 'reset';
print "   <1> Dog\n";
print "   <2> Horse\n";
print "   <3> Human\n\n";


$species_code = <STDIN>;
chomp $species_code;

if ($species_code == 1){$species = "--dog"}
if ($species_code == 2){$species = "--horse"}
if ($species_code == 2){$species = ""}


######################
# Make up file names #
######################

$logfile_all = $filename."_logall.out";
$logfile = $filename.".log";
$results_file= $filename."_results.out";
$final_readme_file = $filename."_README.out";
$results_directory = $filename."_results";



##################################################
# Ask about Analysis scope (brief or full works) #
##################################################

print "\n";
print color 'yellow';
print "####################################################\n";
print "# What scope of analysis do you want to carry out? #\n";
print "####################################################\n\n";
print color 'reset';
print "   <1> Basic - This runs assoc, mperm, QQ-plot, mds, EMMAX, IBS, missingness\n";
print "   <2> Full  - This runs the whole works (and can take a long time)\n";
print "   <3> Pick and choose - This allows you to choose which analyses you do\n\n";

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
	 $do_assoc				= "y";
	 $do_assoc_mperm		= "y";
	 $do_eigenstrat			= "n";
	 $do_emmax				= "y";
	 $do_hap_assoc			= "n";
	 $do_hap_assoc_mperm	= "n";
	 $do_models				= "n";
	 $do_mds				= "y";
	 $do_qq_plot			= "y";
	 $do_homozygosity		= "n";
	 $do_plink_clustering	= "n";
	 $do_which_samples		= "y";
	 $do_missingness		= "y";
	 $do_ibs_test			= "y";
}

##########################
# Analysis scope = full  #
##########################
if ($analysis_scope eq "full")
{
	 $do_assoc				= "y";
	 $do_assoc_mperm		= "y";
	 $do_eigenstrat			= "y";
	 $do_emmax				= "y";
	 $do_hap_assoc			= "y";
	 $do_hap_assoc_mperm	= "y";
	 $do_models				= "y";
	 $do_mds				= "y";
	 $do_qq_plot			= "y";
	 $do_homozygosity		= "y";
	 $do_plink_clustering	= "y";
	 $do_which_samples		= "y";
	 $do_missingness		= "y";
	 $do_ibs_test			= "y";
}



####################################################
# If choose then ask which things you want to run  #
####################################################
if ($analysis_scope eq "choose"){

 while ($choice_ok eq "false")
{
	print "\n\n";
	print color 'yellow';
	print "##########################################\n";
	print "# This is the list of possible analyses. #\n";
	print "##########################################\n\n";
	print color 'reset';
	
	print color 'yellow';
	print "Association:\n\n";
	print color 'reset';
	
	print "    Association analysis:                                \t$do_assoc\n";
	print "    Association mperm analysis:                          \t$do_assoc_mperm\n";
	print "    Association mperm on all models (dom, rec etc):      \t$do_models\n";
	print "    Haplotype association:                               \t$do_hap_assoc\n";
	print "    Haplotype mperm association:                         \t$do_hap_assoc_mperm\n\n\n";

	print color 'yellow';
	print "Stratification:\n\n";
	print color 'reset';
	
	print "    PLINK output for MDS plot:                           \t$do_mds\n";
	print "    PLINK output for QQ plot:                            \t$do_qq_plot\n";
	print "    EMMAX analysis for stratification:                   \t$do_emmax\n";
	print "    IBS test in PLINK:                                   \t$do_ibs_test\n";
	print "    Eigenstrat analysis for stratification:              \t$do_eigenstrat\n";
	print "    Clustering in PLINK for stratification:              \t$do_plink_clustering\n\n\n";
	

	print color 'yellow';
	print "Other analyses:\n\n";
	print color 'reset';
	
	print "    Make files for Excel Homozygosity Mapping:           \t$do_homozygosity\n";
	print "    PLINK missingness analysis:                          \t$do_missingness\n";


	print "\nWould you like to change these? (y/n)    ";

	$ans = <STDIN>;
	chomp $ans;
	$ans = lc $ans;

	if ($ans eq "y")
	{
		$choice_ok = "false";
	
		print "\n\n";
		print color 'yellow';
		print "################################\n";
		print "# Association analysis options #\n";
		print "################################\n\n";
		print color 'reset';
		
		print "    Run association analysis?        	          \t(press return to leave at $do_assoc)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_assoc = "y"};
		if ($ans eq "n"){$do_assoc = "n"};

		print "    Run association mperm analysis?  	          \t(press return to leave at $do_assoc_mperm)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_assoc_mperm = "y"};
		if ($ans eq "n"){$do_assoc_mperm = "n"};

		print "    Run mperm with all models (dom, rec etc)?    \t(press return to leave at $do_models)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_models = "y"};
		if ($ans eq "n"){$do_models = "n"};

		print "    Run association analysis on haplotypes?      \t(press return to leave at $do_hap_assoc)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_hap_assoc = "y"};
		if ($ans eq "n"){$do_hap_assoc = "n"};
		
		print "    Run association MPERM analysis on haplotypes? \t(press return to leave at $do_hap_assoc_mperm)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_hap_assoc_mperm = "y"};
		if ($ans eq "n"){$do_hap_assoc_mperm = "n"};


		print "\n\n";
		print color 'yellow';
		print "###################################\n";
		print "# Stratification analysis options #\n";
		print "###################################\n\n";
		print color 'reset';

		
		print "    Run PLINK for MDS plot?                      \t(press return to leave at $do_mds)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_mds = "y"};
		if ($ans eq "n"){$do_mds = "n"};

		
		print "    Run PLINK for QQ plot?                       \t(press return to leave at $do_qq_plot)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_qq_plot = "y"};
		if ($ans eq "n"){$do_qq_plot = "n"};
		
		
		print "    Run EMMAX analysis?                           \t(press return to leave at $do_emmax)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_emmax = "y"};
		if ($ans eq "n"){$do_emmax = "n"};
		
		
		print "    Do IBS test in PLINK?             	          \t(press return to leave at $do_ibs_test)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_ibs_test = "y"};
		if ($ans eq "n"){$do_ibs_test = "n"};

		
		print "    Run Eigenstrat analysis?                       \t(press return to leave at $do_eigenstrat)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_eigenstrat = "y"};
		if ($ans eq "n"){$do_eigenstrat = "n"};


		print "    Run PLINK with auto-clustering?  	          \t(press return to leave at $do_plink_clustering)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_plink_clustering = "y"};
		if ($ans eq "n"){$do_plink_clustering = "n"};

		print "\n\n";
		print color 'yellow';
		print "##########################\n";
		print "# Other analysis options #\n";
		print "##########################\n\n";
		print color 'reset';

		print "    Make files for 'Homozygosity Mapping'?          \t(press return to leave at $do_homozygosity)  ";
		$ans = <STDIN>;
		chomp $ans;
		$ans = lc $ans;
		if ($ans eq "y"){$do_homozygosity = "y"};
		if ($ans eq "n"){$do_homozygosity = "n"};

		print "    Do PLINK missingness analysis?                  \t(press return to leave at $do_missingness)  ";
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
	print color 'yellow';
	print "####################\n";
	print "# PLINK parameters #\n";
	print "####################\n";
	print color 'reset';
	
	print "\nThese are the default parameters in use:    \n\n";
	print "    --maf $maf_value\n";
	print "    --geno $geno_value\n";
	print "    --mind $mind_value\n";
	print "    --mperm $mperm_value\n";

	if (($do_hap_assoc eq "y") || ($do_hap_assoc_mperm eq "y"))
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


		if (($do_hap_assoc eq "y") || ($do_hap_assoc_mperm eq "y"))
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

########################
# Open overall logfile #
########################

open(LOGALL,">$logfile_all") or die "Can't open file $logfile_all";



######################################
# Create the results folder          #
######################################

while (-e $results_directory)
{
	$directory_suffix = $directory_suffix + 1;
	$results_directory = $filename."_results"."_".$directory_suffix;
}

&run_unix_command("mkdir $results_directory","Create the results directory");


#########################
# open Command Log file #
#########################
$command_log = "$filename"."_wga_analyses_command_log.out";
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

print COMMAND_LOG "COMMAND LOG for WGA_analyses version $version\n\n";
print COMMAND_LOG "PLINK prefix:          \t$prefix\n";
print COMMAND_LOG "Filename with suffix:  \t$filename\n\n";

print COMMAND_LOG scalar localtime;
print COMMAND_LOG "\n\n";

print COMMAND_LOG "Parameters in place:\n\n";
print COMMAND_LOG "--maf $maf_value\n";
print COMMAND_LOG "--geno $geno_value\n";
print COMMAND_LOG "--mind $mind_value\n";
print COMMAND_LOG "--mperm $mperm_value\n\n";
print COMMAND_LOG "Haplotypes from $start_haplotype to $end_haplotype\n\n";
print COMMAND_LOG "Extra command line options for PLINK: \t$pheno_string\n\n";

print "\n";


#################################################################
#                                                               #
#                          Now run PLINK                        #
#                                                               #
#################################################################


###########################################################################
# General check that the data is working and that the filters are correct #
###########################################################################
print color 'bold green';
print "\n\n";
print "############################################\n";
print "#    RUNNING GENERAL CHECK ON THE DATA     #\n";
print "#                                          #\n";
print "# Please CHECK when this has run that the  #\n";
print "# values you used for geno, maf and mind   #\n";
print "# have not filtered off too many SNPs or   #\n";
print "# samples.                                 #\n";
print "############################################\n\n";
print color 'reset';

&run_unix_command("plink --noweb ".$species." --allow-no-sex --maf $maf_value --geno $geno_value --mind $mind_value --file $filename $pheno_string --make-founders --nonfounders --out $filename");

print color 'bold green';

print "\n";
print "################################################################\n";
print "#                 GENERAL CHECK ON THE DATA                    #\n";
print "#                                                              #\n";
print "# Please check at this point that the values you used for GENO #\n";
print "# and MAF have not filtered off too many SNPs or samples.      #\n";
print "################################################################\n";

print "\nPress return to continue.\n";
$ans = <STDIN>;

print color 'reset';

if ($pheno_string ne "")
{
	#####################################################
	# Recode the ped file to include any phenotype data #
	#####################################################
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --file $filename $pheno_string --make-founders --nonfounders --recode --out $filename");

	&print_message("Recoded file using command line option $pheno_string");
}

#########################################################
# Check on the missingness of the data                  #
#########################################################

if ($do_missingness eq "y")
{
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --file $filename $pheno_string --nonfounders --make-founders --test-missing --out $filename","Run --test-missing on the data");

	# Add logfile to logfile_all
	&add_to_logfile;

	&move_to_results_folder (".missing");

	&print_message("Checked on the missingness of the data");

}


#########################################################
# Convert the .map and .ped files in binary files       #
#########################################################

&run_unix_command("plink --noweb ".$species." --allow-no-sex --file $filename $pheno_string --make-founders --nonfounders --make-bed  --maf $maf_value --geno $geno_value --mind $mind_value --out $filename","Convert map and ped files into binary file (also filter with MAF, GENO and MIND)");


# Add logfile to logfile_all
&add_to_logfile;

&print_message("Converted files to binary files");


################################################################################
# Run 'which_samples_were_used' to make a list of samples in the .ped file     #
################################################################################
&run_unix_command("perl /home/genetics/scripts/which_samples_were_used.pl -f $filename",
"Run which_samples_were_used to keep a record of which samples were in the .ped file");

&move_to_results_folder (".samples_used.out");

&print_message("Ran which_samples_were_used to keep a record of which samples were in the .ped file");

#########################################################
# Run standard association analysis with no filtering   #
#########################################################

if ($do_assoc eq "y")
{
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --file $filename $pheno_string --assoc --maf 0 --out $filename"."_unfiltered",
	"Run standard association analysis - unfiltered");

	# Add logfile to logfile_all
	&add_to_logfile;

	&move_to_results_folder ("_unfiltered.assoc");

	&delete_unwanted ("_unfiltered.log");
	
	&print_message("Ran standard association analysis - unfiltered");
	
}


#########################################################
# Run standard association analysis                     #
#########################################################

if ($do_assoc eq "y")
{
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --assoc --out $filename",
	"Run standard association analysis");

	# Add logfile to logfile_all
	&add_to_logfile;

	&move_to_results_folder (".assoc");

	&print_message("Ran standard association analysis");
}


#########################################################
# Produce QQ-plot output                                #
#########################################################

if ($do_qq_plot eq "y")
{
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --assoc --adjust --qq-plot --log10 --out $filename","Run PLINK to get data for QQ plot");

	# Add logfile to logfile_all
	&add_to_logfile;

	&move_to_results_folder (".assoc");
	&move_to_results_folder (".assoc.adjusted");

	&print_message("Produced file with data for QQ plot");
}


################################################################################
# Run plink2homozygosity to create files for Homozygosity Mapping spreadsheet  #
################################################################################

if ($do_homozygosity eq "y")
{
	&run_unix_command("perl /home/genetics/scripts/plink2homozygosity.pl -f $filename","Run PLINK2HOMOZYGOSITY to create files for Homozygosity Mapping spreadsheet");

	&move_to_results_folder ("$filename"."_affected.txt");
	&move_to_results_folder ("$filename"."_normal.txt");

	&print_message("Ran plink2homozygosity to create files for Homozygosity Mapping spreadsheet");
}


#########################################################
# Produce MDS data to look at clustering                #
# (doesn't do it for any particular phenotype)          #
#########################################################

if ($do_mds eq "y")
{
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --mds-plot 2 --out $filename","Run MDS (multi-dimensional-scaling) analysis corrected for multiple testing to look at clustering");

	&add_to_logfile;

	&move_to_results_folder (".mds");

	&print_message("Ran MDS analysis to look at clustering");
}


############################################################
# Run association corrected for multiple testing  (MPERM)  #
############################################################

if ($do_assoc_mperm eq "y")
{
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --assoc --mperm $mperm_value --out $filename","Run association analysis corrected for multiple testing");

	&add_to_logfile;

	&move_to_results_folder (".assoc.mperm");

	&print_message("Ran standard association analysis corrected for multiple testing using mperm");
}


############################################################
# Only include all models if analysis scope is set to full #
############################################################

if ($do_models eq "y")
{
	#########################################################
	# Run association with all models (dom,rec etc)         #
	#########################################################
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --model --cell 0 --out $filename","Run PLINK to look at different models (dom, rec etc)");

	&add_to_logfile;
	&move_to_results_folder (".model");
	&print_message("Looked at different models (dom, rec etc)");

	
	#########################################################
	# Run association mperm with best model                 #
	#########################################################
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --model --cell 0 --mperm $mperm_value --out $filename",
	"Run PLINK to look at best model with mperm");

	&add_to_logfile;
	&move_to_results_folder (".model");
	&move_to_results_folder (".model.best.mperm");

	&print_message("Looked at best model with mperm");


	#########################################################
	# Run mperm association with the dominant model         #
	#########################################################
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --model-dom --cell 0 --mperm $mperm_value --out $filename",
	"Run PLINK to look at dominant model with mperm");

	&add_to_logfile;
	&move_to_results_folder (".model");
	&move_to_results_folder (".model.dom.mperm");

	&print_message("Looked at dominant model with mperm");


	#########################################################
	# Run mperm association with the recessive model        #
	#########################################################
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --model-rec --cell 0 --mperm $mperm_value --out $filename",
	"Run PLINK to look at recessive model with mperm");

	&add_to_logfile;
	&move_to_results_folder (".model");
	&move_to_results_folder (".model.rec.mperm");

	&print_message("Looked at recessive model with mperm");


	#########################################################
	# Run mperm association with the genotypic model        #
	#########################################################
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --model-gen --cell 0 --mperm $mperm_value --out $filename",
	"Run PLINK to look at genotypic model with mperm");

	&add_to_logfile;
	&move_to_results_folder (".model");
	&move_to_results_folder (".model.gen.mperm");

	&print_message("Looked at genotypic model with mperm");


	#########################################################
	# Run mperm association with the trend model            #
	#########################################################
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --model-trend --cell 0 --mperm $mperm_value --out $filename",
	"Run PLINK to look at trend model with mperm");
	
	&add_to_logfile;
	&move_to_results_folder (".model");
	&move_to_results_folder (".model.trend.mperm");

	&print_message("Looked at trend model with mperm");

} # do_models = y



#############################################################
# Only do plink clustering if analysis scope is set to full #
#############################################################

if ($do_plink_clustering eq "y")
{
	#########################################################
	# Allow PLINK to cluster the samples                    #
	#########################################################
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --cluster --cc --ppc 0.01 --out $filename",
	"Run PLINK to produce clusters with at least one case and one control and ppc=0.01 (read PLINK documentation)");

	&add_to_logfile;

	&print_message("Allowed PLINK to cluster the samples");


	#########################################################
	# Run association only within these clusters            #
	#########################################################
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --mh --within $filename".".cluster2 --out $filename",
	"Run MH association within the cluster2 file");

	&add_to_logfile;
	&move_to_results_folder (".cmh");

	&print_message("Looked for assication within these clusters");


	###########################################################
	# Run association only within these clusters using mperm  #
	###########################################################
	&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --mh --mperm $mperm_value --within $filename".".cluster2 --out $filename",
	"Run MH association within the cluster2 file using mperm");

	&add_to_logfile;
	&move_to_results_folder (".cmh");
	&move_to_results_folder (".cmh.mperm");
	&move_to_results_folder (".cluster2");

	&print_message("Looked for assication within these clusters using Mperm");

} # of analysis = full or do_plink_clustering = "y"


############################################################
# Haplotype association analysis                           #
############################################################

if ($do_hap_assoc eq "y")
{

	######################################################################
	# Haplotype ASSOC analysis - from $start_haplotype to $end_haplotype #
	######################################################################

	if ($start_haplotype > 0 && $start_haplotype < $end_haplotype)
	{

		for ($haplotype_count = $start_haplotype;$haplotype_count <= $end_haplotype;$haplotype_count ++)
		{

			##########################################################
			# Haplotype analysis  - generate haplotypes              #
			# For straight association analysis there is no need     #
			# to impute the haplotypes first, just run the analysis. #
			##########################################################
			&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --hap-window $haplotype_count --hap-assoc --out $filename"."_"."$haplotype_count",
			"Run Haplotyping Assoc analysis on window of $haplotype_count");

			#plink --noweb --dog --bfile NAME --hap-window 4 --hap-assoc --out NAME_w4

			&add_to_logfile;
			
			&print_message("Ran Haplotyping Assoc analysis  (window = $haplotype_count)");


			########################################################
			# Convert impute.map and impute.assoc.mperm files to   #
			# a single file which can be plotted in Hap_Mperm_plot #
			########################################################

			&run_unix_command("perl /home/genetics/scripts/hap_assoc_convert.pl -file $filename"."_"."$haplotype_count -hap_length $haplotype_count -species $species -map $filename.map",
			"Converted assoc.hap file to a _hap_assoc_{n}.txt file which can be plotted in 'Haplotype Association Plot'");
			
			
			&run_unix_command("perl /home/genetics/scripts/hap_mperm_convert.pl -file $filename -hap_length $haplotype_count -species $species -map $filename.map",
			"Converted impute.map and impute.assoc.mperm files to single _hap_mperm_{n}.txt file which can be plotted in 'Haplotype Association Plot'");

			
			&move_to_results_folder ("_".$haplotype_count.".assoc.hap");
			&move_to_results_folder ("_hap_assoc_{".$haplotype_count."}.txt");

			&print_message("Converted output files with hap_assoc_convert  (window = $haplotype_count)");


			##########################
			# Deleted unwanted files #
			##########################
			&delete_unwanted ("_".$haplotype_count.".log");
			&delete_unwanted ("_".$haplotype_count.".hh");

		} # end of $haplotype_count loop

	} # end of if start_haplotype > 0

} # if do_hap_assoc = y



############################################################
# Haplotype association analysis  MPERM                    #
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
			&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename $pheno_string --hap-window $haplotype_count --hap-impute --out $filename"."_"."$haplotype_count",
			"Run Haplotyping analysis on window of $haplotype_count");

			&add_to_logfile;
			
			&print_message("Imputed the haplotypes (window = $haplotype_count)");


			#######################################################################
			# Haplotype analysis  - make binary files and delete non-binary files #
			#######################################################################
			
			&run_unix_command("plink --noweb ".$species." --allow-no-sex --file $filename"."_".$haplotype_count.".impute $pheno_string --make-bed --out $filename"."_".$haplotype_count.".impute",
			"Made binary files for haplotyping on window of $haplotype_count");

			&add_to_logfile;
			
			&print_message("Made binary files for haplotyping  (window = $haplotype_count)");

			
			##################################
			# Delete the non-binary versions #
			##################################
			&delete_unwanted ("_".$haplotype_count.".impute.ped");
			
			
			#######################################################
			# Haplotype analysis  - association with permutation  #
			#######################################################
			&run_unix_command("plink --noweb ".$species." --allow-no-sex --bfile $filename"."_".$haplotype_count.".impute $pheno_string --assoc --mperm $mperm_value --out $filename"."_".$haplotype_count.".impute",
			"Ran Haplotype analysis to use in Hap_Mperm_Plot");
			
			&add_to_logfile;
			
			&print_message("Ran Haplotyping to use in Hap_Mperm_Plot (window = $haplotype_count)");


			########################################################
			# Convert impute.map and impute.assoc.mperm files to   #
			# a single file which can be plotted in Hap_Mperm_plot #
			########################################################

			&run_unix_command("perl /home/genetics/scripts/hap_mperm_convert.pl -file $filename -hap_length $haplotype_count -species $species -map $filename.map",
			"Converted impute.map and impute.assoc.mperm files to single _hap_mperm_{n}.txt file which can be plotted in 'Haplotype Association Plot'");

			
			&move_to_results_folder ("_".$haplotype_count.".impute.map");
			&move_to_results_folder ("_".$haplotype_count.".impute.assoc.mperm");
			&move_to_results_folder ("_hap_mperm_{".$haplotype_count."}.txt");

			&print_message("Converted output files with hap_mperm_convert  (window = $haplotype_count)");


			##########################
			# Deleted unwanted files #
			##########################
			&delete_unwanted ("_".$haplotype_count.".impute.bed");
			&delete_unwanted ("_".$haplotype_count.".impute.bim");
			&delete_unwanted ("_".$haplotype_count.".impute.fam");
			&delete_unwanted ("_".$haplotype_count.".impute.nosex");
			&delete_unwanted ("_".$haplotype_count.".impute.assoc");
			&delete_unwanted ("_".$haplotype_count.".impute.log");
			&delete_unwanted ("_".$haplotype_count.".log");
			&delete_unwanted ("_".$haplotype_count.".hh");

		} # end of $haplotype_count loop

	} # end of if start_haplotype > 0

} # if do_hap_assoc_mperm = y



if ($do_ibs_test eq "y")
{
	 #########################################################
	 # IBS TEST in PLINK                                     #
	 #########################################################
	 &run_unix_command("plink --noweb ".$species." --allow-no-sex --file $filename $pheno_string --nonfounders --make-founders --ibs-test --out $filename",
	 "Ran ibs-test in PLINK");

	 &print_message("Ran ibs-test in PLINK");

	####################################
	# Rename log file to ibs-test file #
	####################################
	 $command = "mv $filename".".log $filename".".ibs-test"; 
	 system("$command");

	 &move_to_results_folder (".ibs-test");

}


################################################################
# EMMAX analysis                                               #
################################################################

if ($do_emmax eq "y")
{
	#############################
	# Run PERL script run_emmax #
	#############################
	&run_unix_command("perl /home/genetics/test_scripts/run_emmax_07.pl --species $species --file $filename $pheno_string --maf $maf_value --geno $geno_value --mind $mind_value",
	"Ran EMMAX to produce CHISQs adjusted for stratification");

	&print_message("Ran EMMAX to produce CHISQs adjusted for stratification");
		
	# Rename ps file to emmax_pvalues.out
	$command = "mv $filename".".ps $filename"."_emmax_pvalues.out"; 
	system("$command");
	 
	&move_to_results_folder (".emmax_log.out");
	&move_to_results_folder (".emmax_out");
	&move_to_results_folder (".hd_emmax_out");
	&move_to_results_folder (".emmax_qq");
	&move_to_results_folder (".hBN.kinf");
	&move_to_results_folder (".hIBS.kinf");
	
	&move_to_results_folder ("$filename"."_emmax_qq.pdf");
	
	# Delete intermediate files
	&delete_file("$filename".".tped");
	&delete_file("$filename".".tfam");
	&delete_file("$filename".".pheno");
	&delete_file("$filename".".reml");
	&delete_file ("$filename"."_emmax_qq.ps");
	&delete_file ("gnufile.xtxt");
	&delete_file("$filename"."_emmax_pvalues.out");
	
	 
} # end of if do EMMAX = y


################################################################
# Eigenstrat analysis                                          #
################################################################

if ($do_eigenstrat eq "y")
{
	 #########################################################
	 # Convertf - to convert PLINK files to Eigenstrat files #
	 #########################################################
	 &run_unix_command("perl /home/genetics/scripts/runconvertf_".$species.".pl -f $filename",
	 "Run convertf to convert PLINK files to Eigenstrat files");

	 &print_message("Ran convertf to convert PLINK files to Eigenstrat files");

	 
	 #########################################################
	 # Smartpca - to carry out PCA analysis                  #
	 #########################################################
	 &run_unix_command("perl /home/genetics/scripts/runsmartpca_".$species.".pl -f $filename",
	 "Ran SMARTPCA to carry out Principal Components Analysis");
		
	 &print_message("Ran SMARTPCA to carry out Principal Components Analysis");

	 
	 #########################################################
	 # Eigenstrat - to account for structure                 #
	 #########################################################
	 &run_unix_command("perl /home/genetics/scripts/runeigenstrat_".$species.".pl -f $filename",
	 "Ran EIGENSTRAT to produce CHISQs adjusted for stratification");


	 # Rename log file to eigenlog
	 $command = "mv $filename".".log $filename".".eigenlog"; 
	 system("$command");

	 &move_to_results_folder (".eigenlog");
	 &move_to_results_folder ("$filename"."_PCA_plot.pdf");
	 &move_to_results_folder (".smartpca_log");
	 &move_to_results_folder (".eigenstrat.out");
	 &move_to_results_folder (".chisq");

	 &print_message("Ran EIGENSTRAT to produce CHISQs adjusted for stratification");
	 
	####################################
	# Remove unwanted Eigenstrat files #
	####################################
	&delete_unwanted (".eigenstratgeno");
	&delete_unwanted (".par");
	&delete_unwanted (".pca");
	&delete_unwanted (".eval");
	&delete_unwanted (".evec");
	&delete_unwanted (".pheno");
	&delete_unwanted ("_PCA_plot.ps");
	&delete_unwanted ("_PCA_plot.xtxt");
	&delete_unwanted (".snp");
	&delete_unwanted (".ind");
	
} # end of if do eigenstrat = y


close LOGALL;

&move_to_results_folder ("$logfile_all");
&move_to_results_folder (".assoc");


#################################################################
# Write the final README file to tell user which files are what #
#################################################################


open(README,">$final_readme_file") or die "Can't open file $final_readme_file";

print README "################################################\n";
print README "# List of the output files from WGA_analyses   #\n";
print README "# and some information on what to do with them #\n";
print README "################################################\n\n";

print README "(Please also read any instructions in the Genotyping Analysis/Instructions Folder)\n\n";

print README "Filename: \t$filename\n";
print README "Date:     \t";
print README scalar localtime;
print README "\n\n";

print README "-----------------------------------------------------------------------------------------\n";

# Missingness #
print README "The following file has a list of all the command lines used to run these analyses:\n\n";

print README " > $command_log\n\n";

print README "-----------------------------------------------------------------------------------------\n";

# Missingness #
if ($do_missingness eq "y")
{
	print README "The following file can be used to look at any missing data\n";
	print README "to make sure any 'missingness' is the same in the cases and controls\n\n";
	print README " > ".$filename.".missing\n\n";
	print README "Open this file in Excel and sort by the P-value column.\n";
	print README "This tells you if there is a significance difference between cases and controls for any SNP.\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}


# Basic association #
if ($do_assoc eq "y")
{
	 print README "The following file has the simple association data.\n";
	 print README "It can be plotted using the Excel sheet PLINK_plot\n\n";
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
	print README "(Note you must read it in this file with the special button on the Front Page)\n";
	print README "(This currently doesn't work. If you want to use it ask Mike to fix it.)\n\n";

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
	print README " > ".$filename.".ibs_test\n\n";
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
	print README " > ".$filename.".assoc.adjusted\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}


if ($do_which_samples eq "y")
{
	# which_samples_were_used #
	print README "The following file keeps a list of the samples used in the .ped file, and their disease statuses:\n\n";
	print README " > ".$filename.".samples_used.out\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}


# Homozygosity Mapping files #
if ($do_homozygosity eq "y")
{
	print README "The following files can be used in the Excel 'Homozygosity Mapping' spreadsheet\n";
	print README "and then in Excel sheets 'Haplotyping Display' and 'Analyse Critical Region'\n\n";
	print README "Open them in two separate worksheets called Affected and Normal\n";
	print README "(WARNING: This will not use any phenotypes which are in a separate file from the .ped file)\n\n";
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
		print README "The files with numbers in curly brackets have been converted to work directly in PLINK_plot (or Hap_Mperm_plot) with out having to read it in with the button.\n\n";
		print README " > e.g. ".$filename.".hap_mperm_{4}.txt\n\n"
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
	print README " > ".$filename."_PCA_plot.pdf\n\n";
	print README "This file has information from the Eigensoft analysis, including statistical measures of the significance of the stratification.\n\n";
	print README " > ".$filename.".smartpca_log\n\n";
	print README "This file has CHISQ association data, uncorrected and also corrected for stratification.\n";
	print README "(Use the Excel sheet Eigenstrat_plot to plot these two columns of data)\n\n";
	print README " > ".$filename.".eigenstrat_out OR ".$filename.".chisq\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}

#EMMAX #
if ($do_emmax eq "y")
{
	print README "The following files have the results of running EMMAX on the data.\n\n";
	print README "This is a mixed model approach to correct for population structure. \n";
	print README "(NOTE: We slightly prefer to use Fast Mixed Model, but the EMMAX results are very similar.)\n\n";
	
	print README "This file shows a QQ plot of the EMMAX-corrected data\n\n";
	print README " > ".$filename."_emmax_qq.pdf\n\n";
	
	print README "This file has the data to plot your own QQ plots\n\n";
	print README " > ".$filename.".emmax_qq\n\n";
	
	print README "This file has the adjusted EMMAX P-value data to produce an association plot after EMMAX correction\n\n";
	print README " > ".$filename.".emmax_out\n\n";
	
	print README "There is also a version of this file for direct plotting in PLINK_plot:\n\n";
	print README " > ".$filename.".hd_emmax_out\n\n";
	
	print README "These file have kinship data.  rPlease read the EMMAX documentation for full details.\n\n";
	print README " > ".$filename.".hBN.kinf and ".$filename.".hBN.kinf\n\n";
	
	print README "-----------------------------------------------------------------------------------------\n";
}



print README "\nIf there any problems with this program, please see Mike Boursnell.\n\n";
print README "e-mail: mike.boursnell AT aht.org.uk\n\n";

print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	

#######################################################
# Copy READEM and COMMAND_LOG into the results folder #
#######################################################

&print_message("Moving some files to the results folder");
&move_to_results_folder ("$final_readme_file");
&move_to_results_folder ("$command_log");


###########################################
# Remove some unwanted intermediate files #
###########################################

&print_message("Removing some unwanted files");
print COMMAND_LOG "Removing some unwanted files...\n\n";


if ($do_plink_clustering eq "y")
{
	&delete_unwanted (".cluster0");
	&delete_unwanted (".cluster1");
	&delete_unwanted (".cluster3");
}

#remove plink files
&delete_unwanted (".log");
&delete_unwanted (".ind");
&delete_unwanted (".hh");
&delete_unwanted (".irem");
#&delete_unwanted (".map");  #Useful to keep this in results folder
#&delete_unwanted (".ped");  #Useful to keep this in results folder

#move ped and map files to results folder
&move_to_results_folder (".ped");
&move_to_results_folder (".map");

#remove all PLINK files of these types
delete_all_of_this_type (".nosex");
delete_all_of_this_type (".nof");

#####################################################################################
# Delete the binary versions of the files (can easily regenerate these if required) #
#####################################################################################
 &delete_unwanted (".fam");
 &delete_unwanted (".bim");
 &delete_unwanted (".bed");


#remove plink2homozygosityfiles
if ($do_homozygosity eq "y") 
{
	&delete_unwanted (".p2h_log");
}

#remove gmon.out file
if (-e "gmon.out"){&delete_file ("gmon.out")}

if ($email_address ne "")
{
	open(MAIL, "|/usr/sbin/sendmail -t");

		## Mail Header
		print MAIL "To: $email_address\n";
		print MAIL "From: $from\n";
		print MAIL "Subject: WGA_analyses on $filename has finished\n\n";
		
		## Mail Body
		print MAIL "Your Whole Genome Association analysis on $filename is complete\n\n";
		print MAIL "PERL script:       \tWGA_analyses version $version\n\n";
		print MAIL "Results directory: \t$results_directory\n\n";
		print MAIL "Check the following files:\n\n";
		print MAIL "    1.  How to use the output files:        \t$final_readme_file\n";
		print MAIL "    2.  A list of all the unix commands:    \t$command_log\n";
		print MAIL "    3.  Program logs of the analyses:       \t$logfile_all\n\n";
	
	close(MAIL);
}

###################################################################################
#                  End message to say the program has finished                    #
###################################################################################
print color 'bold green';

print"\n\n\n";
print "###############################################################\n";
print "###############################################################\n";
print "##                 All GWAS analyses complete                ##\n";
print "###############################################################\n";
print "###############################################################\n\n\n";

print "Check the following files:\n\n";
print "(These have been moved into a folder called $filename"."_results)\n\n";
print "    1.  How to use the output files:        \t$final_readme_file\n";
print "    2.  A list of all the unix commands:    \t$command_log\n";
print "    3.  Program logs of the analyses:       \t$logfile_all\n\n";

print color 'reset';
close README;
close COMMAND_LOG;

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
		
	}

	# Add line to overall logfile
	print LOGALL "$single_line\n";


} # end of while loop

if ($errors_found eq "true"){print COMMAND_LOG "No errors found in PLINK logfile\n\n"}

close LOG;

}

#############################################
# Subroutine to move file to results folder #
#############################################

sub move_to_results_folder
{
	my $suffix = "";
	$suffix = $_[0];
	
	
	# if the argument is a suffix added to the filename (ie begins with a dot or underscore)
	if ((index($suffix,".") == 0) || (index($suffix,"_") == 0))
	{
		$command = "mv  $filename".$suffix." $results_directory/$filename".$suffix;
		system("$command");
		print COMMAND_LOG ("$command\n");
		print COMMAND_LOG ("Moved $filename".$suffix." to results directory\n");
	}
	
	# if the argument is a complete filename (ie doesn't begin with a dot)
	if ((index($suffix,".") != 0) && (index($suffix,"_") != 0))
	{
		$command = "mv  $suffix $results_directory/$suffix";
		system("$command");
		print COMMAND_LOG ("$command\n");
		print COMMAND_LOG ("Moved $suffix to results directory\n");
	}	
}

#############################################
# Subroutine to remove unwanted files       #
#############################################

sub delete_unwanted
{
	my $suffix = "";
	$suffix = $_[0];

	if (-e "$filename".$suffix)
	{ 
		$command = "rm  $filename".$suffix;
		system("$command");
		#print COMMAND_LOG ("$command\n"); # Leave these out as they clutter the command_log file
	}
	else
	{
		print "    >>> The file $filename".$suffix." doesn't exist <<<\n";
	}
}


#############################################
# Subroutine to remove all files of a type  #
#############################################

sub delete_all_of_this_type
{
	my $suffix = "";
	$suffix = $_[0];
	$command = "rm  *".$suffix;
	print("$command\n");
	system("$command");
	print COMMAND_LOG ("$command\n");
}


#############################################
# Subroutine to remove unwanted files       #
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
		#print COMMAND_LOG ("$command\n");
	}
	else
	{
		print "    >>> The file $thisfile doesn't exist <<<\n";
	}

}

#############################################
# Subroutine to execute unix command        #
#############################################

sub run_unix_command
{
	my $unix_command = "";
	my $information = "";
	$unix_command = $_[0];
	$information = $_[1]; # Optional?

	print "\n";
	print("$unix_command\n");
	system("$unix_command");
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	print COMMAND_LOG "*** $information ***\n";
	print COMMAND_LOG "$unix_command\n";

}

######################################
# Subroutine to print screen message #
######################################

sub print_message
{
	my $message 		= "";
	my $message_length 	= "";
	my $pos_count		= 0;
	
	$message = $_[0];
	$message_length = length($message);
	
	print "\n\n";
	print color 'yellow';
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++)
	{
		print "#";
	}
	print "\n#    $message    #\n";
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++)
	{
		print "#";
	}
	print "\n\n";
	
	print color 'reset';

}
