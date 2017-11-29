#!/usr/bin/perl

#########################################################################
#									                                    #      
#	GWAS_analyses           						                        #     
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
my $version				= "49";
my $testing_mode		= "off"; # Stops at each stage to allow you to check. off or on

####################
# Define variables #
####################

my $pheno_name			= "none";
my $plink_pheno_string	= "";
my $pheno_file_string	= "";
my $pheno_name_string	= "";
my $analysis_scope		= "brief";
my $from				= 'GWAS_analyses@samba64.aht.org.uk';


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
my $do_sex_check		= "y";

#Boolean
my $parameters_ok		= "false";
my $phenotypes_ok		= "false";
my $choice_ok			= "false";
my $use_plink2_if_poss	= "yes"; # yes or no
my $split_x				= "yes"; # if you want to split X and Y into X,Y and XY

my $answer					= "";
my $species_code		= "";
my $species				= "";
my $species_string_perl		= "";
my $species_string_plink	= "";
my $command				= "";
my $single_line			= "";
my $errors_found		= "";
my $results_folder	= "";
my $suffix				= "";
my $analysis_suffix		= "";
my $email_address		= "";
my $plink_extra_string	= "";

#File names
my $command_log					= "";
my $split_x_command_log			= "";
my $logfile						= "";
my $logfile_all					= "";
my $results_file				= "";
my $final_readme_file			= "";
my $prefix						= "";
my $pedfile						= "";
my $mapfile						= "";
my $prefix_with_suffix			= "";
my $prefix_with_suffix_unsplit	= "";
my $prefix_split				= "";
my $pheno_file					= "none";
my $last_pheno_file				= "none";
my $model_best_mperm_file		= ""; # from PLINK2

my $error_pos					= 0;
my $haplotype_count				= 0;

my $mperm_value					= 100000;
my $mperm_value_haplotyping		= 10000;
my $geno_value					= 0.1;
my $maf_value					= 0.05;
my $mind_value					= 0.1;
my $start_haplotype				= 1;
my $end_haplotype				= 8;
my $directory_suffix			= 0;
my $no_of_samples				= 0;
my $no_of_phenotypes			= 0;
my $phenotype_count		= 0;
my $pheno_column		= "";
my $pheno_file_input	= "";
my $snp_name			= "";
my $position			= "";
my $assembly			= "";



# Arrays
my @pheno_file_array	= ();
my @item				= ();


print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "       GWAS analyses       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This program can run all the standard PLINK analyses on your GWAS data\n";
print "  - It can also run analyses by other programs such as EMMAX etc.\n\n";

print color 'bold cyan';

print "    (The input files are PLINK format ped and map files)\n\n";

print color 'reset';


###########################################
# Get the file names and check the files  #
###########################################

if ($prefix eq "")
{
	until ((-e $pedfile) && (-e $mapfile))
	{
		print_message("Prefix for PLINK .map and .ped files","input");
		print ">  ";
		$prefix = <STDIN>;
		chomp $prefix;

		###########################################
		# If user has typed in full PED file name #
		###########################################
		if (index ($prefix,".ped") > -1){$prefix = &get_prefix($prefix);}

		if ($prefix eq "ls"){print "\n";system ("ls *.ped");print "\n"}
		if ($prefix ne "ls")
		{

			$pedfile = $prefix.".ped";
			$mapfile = $prefix.".map";

			if (! -e $pedfile){print "\n>>>>>>>>  PED file $pedfile not found.  Try again.  <<<<<<<<\n\n";}

			if (! -e $mapfile){print "\n>>>>>>>>  MAP file $mapfile not found.  Try again.  <<<<<<<<\n\n";}
		}
	}
}


###########################################
# Get suffix for this analysis            #
###########################################

&print_message("Choose a suffix for this analysis (e.g. a, b, c...):      ","input");

$analysis_suffix = "ls";
until ($analysis_suffix ne "ls")
{
	$analysis_suffix = <STDIN>;
	chomp $analysis_suffix;

	if ($analysis_suffix eq "ls")
	{
		system ("ls -d $prefix*results");
		print "\n";
	}
	print "\n";
}


###########################################
# E-mail address for notifications        #
###########################################

&print_message("Please enter your email address","input");

print "(Just typing the first letter of your first name is a short cut)\n\n";

$email_address = <STDIN>;
$email_address = lc $email_address;
chomp $email_address;


###########################################
# Some short cuts for e-mails             #
###########################################
if ($email_address eq "b"){$email_address = 'rebekkah.hitti@aht.org.uk';}
if ($email_address eq "c"){$email_address = 'chris.jenkins@aht.org.uk';}
if ($email_address eq "g"){$email_address = 'graham.newland@aht.org.uk';}
if ($email_address eq "j"){$email_address = 'james.oliver@aht.org.uk';}
if ($email_address eq "m"){$email_address = 'mike.boursnell@aht.org.uk';print 'mike.boursnell@aht.org.uk'; print"\n\n";}
if ($email_address eq "o"){$email_address = 'oliver.forman@aht.org.uk';}
if ($email_address eq "r"){$email_address = 'rebekkah.hitti@aht.org.uk';}
if ($email_address eq "s"){$email_address = 'sally.ricketts@aht.org.uk';}


###########################################
# Make new filename with suffix           #
###########################################
$prefix_with_suffix = $prefix."_".$analysis_suffix;


###########################################
# Ask what species                        #
###########################################

&print_message("Which species are you analysing?","input");

print "   <1> Dog\n";
print "   <2> Horse\n";
print "   <3> Human\n\n";


$species_code = <STDIN>;
chomp $species_code;

if (substr($species_code,0,1) eq "1" ){$species = "dog";  $species_string_plink = "--dog";$species_string_perl = "--species dog"}
if (substr($species_code,0,1) eq "2" ){$species = "horse";$species_string_plink = "--horse";$species_string_perl = "--species horse"}
if (substr($species_code,0,1) eq "3" ){$species = "human";$species_string_plink = "";$species_string_perl = "--species human"}
if (substr($species_code,0,1) eq "4" ){$species = "cat";  $species_string_plink = "";$species_string_perl = "--species cat"; print "\nSpecies is cat\n\n";}



#########################################################
# Check (if dog) whether MAP file is canfam2 or canfam3 #
#########################################################
if ((-e $mapfile) && ($species eq "dog"))
{
 	open(MAP,"$mapfile") or die "Can't open file $mapfile";

	 while ($single_line = <MAP>) 
	{
		chomp $single_line;
		@item=split(/\s+/,$single_line);
		$snp_name = $item[1];
		$position = $item[3];

		#print "$position\n";

		if ($snp_name eq "BICF2G630708027")
		{
			$assembly = "NOT canfam2 or canfam3";
			if ($position eq "3713883"){$assembly = "canfam2"; last}
			if ($position eq "714375"){$assembly = "canfam3"; last}
		}
	}
	close MAP;

	if ($assembly eq "canfam2")
	{
		&print_message( "WARNING!!   The MAP file $mapfile is based on canine assembly $assembly","message");
		print "If this is a mistake, quit the program now.\n\n";
		$answer=<STDIN>;
	}
	if ($assembly eq "NOT canfam2 or canfam3")
	{
		&print_message( "WARNING!!   The MAP file $mapfile doesn't seem to be based on canfam2 or canfam3","message");
		print "If this is a mistake, quit the program now.\n\n";
		$answer=<STDIN>;
	}
	
} # dog


##################################################
# Ask about using PLINK2                         #
##################################################

&print_message("Do you want to use the new PLINK2 where possible?","input");

print "  This speeds up analyses a lot, (for example MPERM)\n\n";
print "  It also allows better results in cases like check-sex etc.\n\n";

print "   <1> Yes - use PLINK2 where possible           [DEFAULT]\n";
print "   <2> No -  stick to original PLINK\n\n";

$answer = <STDIN>;
chomp $answer;
print"\n\n";
if (substr($answer,0,1) eq "2" )
{$use_plink2_if_poss = "no"} 
else 
{$use_plink2_if_poss = "yes";$mperm_value = 100000} # If you use PLINK2 the mperm analysis runs MUCH faster


##################################################
# Ask about using PLINK2  split-x                #
##################################################

&print_message("Do you want to split the X and Y chromosomes into the correct PLINK format of X, Y and XY?","input");

print "  (This means that the pseudoautosomal region is not excluded)\n\n";

print "   <1> Yes - use correct PLINK format for X and Y            [DEFAULT]\n";
print "   <2> No -  stick to original X and Y chromosomes [will get lots of heterozygous haploid SNPs excluded]\n";
print "   <3> The files are already split\n\n";

$answer = <STDIN>;
chomp $answer;
print"\n\n";

if ($answer eq ""){$answer = "1"} # default

if (substr($answer,0,1) eq "1" ){$split_x = "yes"}
if (substr($answer,0,1) eq "2" ){$split_x = "no"}
if (substr($answer,0,1) eq "3" ){$split_x = "no"}

if ($split_x eq "yes")
{
	$prefix_with_suffix_unsplit = $prefix_with_suffix;

	$prefix_with_suffix = $prefix_with_suffix."_split";
}

##################################################
# Make up some file names                        #
##################################################

$logfile_all = $prefix_with_suffix."_logall.out";
$logfile = $prefix_with_suffix.".log";
$results_file= $prefix_with_suffix."_results.out";
$final_readme_file = $prefix_with_suffix."_README.out";
$results_folder = $prefix_with_suffix."_results";


##################################################
# Ask about Analysis scope (brief or full works) #
##################################################

&print_message("What scope of analysis do you want to carry out?","input");

print "   <1> Basic - This runs assoc, mperm, QQ-plot, mds, EMMAX, IBS, missingness              [DEFAULT]\n";
print "   <2> Full  - This runs the whole works (and can take a long time)\n";
print "   <3> Pick and choose - This allows you to choose which analyses you do\n\n";

$answer = <STDIN>;
chomp $answer;
print"\n\n";

if ($answer eq ""){$answer = "1"} # default

if (substr($answer,0,1) eq "1" ){$analysis_scope = "basic"}
if (substr($answer,0,1) eq "2" ){$analysis_scope = "full"}
if (substr($answer,0,1) eq "3" ){$analysis_scope = "choose"}


##############################################
# Set default 'basic' values before choosing #
##############################################
$do_sex_check			= "y";
$do_assoc				= "y";
$do_assoc_mperm			= "y";
$do_emmax				= "y";
$do_mds					= "y";
$do_qq_plot				= "y";
$do_homozygosity		= "y";
$do_plink_clustering	= "n";
$do_which_samples		= "y";
$do_missingness			= "y";
$do_ibs_test			= "y";

$do_hap_assoc			= "n";
$do_hap_assoc_mperm		= "n";
$do_models				= "n";
$do_eigenstrat			= "n";
$plink_extra_string		= "";


##############################################
# Analysis scope = basic                     #
##############################################
if ($analysis_scope eq "basic")
{
	 $do_sex_check			= "y";
	 $do_assoc				= "y";
	 $do_assoc_mperm		= "y";
	 $do_eigenstrat			= "n";
	 $do_emmax				= "y";
	 $do_hap_assoc			= "n";
	 $do_hap_assoc_mperm	= "n";
	 $do_models				= "n";
	 $do_mds				= "y";
	 $do_qq_plot			= "y";
	 $do_homozygosity		= "y";
	 $do_plink_clustering	= "n";
	 $do_which_samples		= "y";
	 $do_missingness		= "y";
	 $do_ibs_test			= "y";
	 $plink_extra_string		= "";
}

##############################################
# Analysis scope = full                      #
##############################################
if ($analysis_scope eq "full")
{
	 $do_sex_check			= "y";
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
	 $do_plink_clustering	= "n"; # Don't do this any more as standard
	 $do_which_samples		= "y";
	 $do_missingness		= "y";
	 $do_ibs_test			= "y";
	 $plink_extra_string	= "";
}



######################################################
# If "choose" then ask which things you want to run  #
######################################################
if ($analysis_scope eq "choose")
{
	 while ($choice_ok eq "false")
	{
		print "\n\n";

		&print_message("This is the list of possible analyses.","message");

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
		print "    Eigenstrat analysis for stratification:              \t$do_eigenstrat\n\n\n";
		
		print color 'yellow';
		print "Other analyses:\n\n";
		print color 'reset';
		
		print "    Make files for Excel Homozygosity Mapping:           \t$do_homozygosity\n";
		print "    PLINK missingness analysis:                          \t$do_missingness\n\n\n";

		print color 'yellow';
		print "Extra command line options:\n\n";
		print color 'reset';

		print "    Extra command line options:                          \t$plink_extra_string\n";


		print "\nWould you like to change these? (y/n)    ";

		$answer = <STDIN>;
		chomp $answer;
		$answer = lc $answer;

		if ($answer eq "y")
		{
			$choice_ok = "false";
		
			print "\n\n";
			print color 'yellow';
			print "################################\n";
			print "# Association analysis options #\n";
			print "################################\n\n";
			print color 'reset';
			
			print "    Run association analysis?        	          \t(press return to leave at $do_assoc)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_assoc = "y"};
			if ($answer eq "n"){$do_assoc = "n"};

			print "    Run association mperm analysis?  	          \t(press return to leave at $do_assoc_mperm)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_assoc_mperm = "y"};
			if ($answer eq "n"){$do_assoc_mperm = "n"};

			print "    Run mperm with all models (dom, rec etc)?    \t(press return to leave at $do_models)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_models = "y"};
			if ($answer eq "n"){$do_models = "n"};

			print "    Run association analysis on haplotypes?      \t(press return to leave at $do_hap_assoc)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_hap_assoc = "y"};
			if ($answer eq "n"){$do_hap_assoc = "n"};
			
			print "    Run association MPERM analysis on haplotypes? \t(press return to leave at $do_hap_assoc_mperm)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_hap_assoc_mperm = "y"};
			if ($answer eq "n"){$do_hap_assoc_mperm = "n"};


			print "\n\n";
			print color 'yellow';
			print "###################################\n";
			print "# Stratification analysis options #\n";
			print "###################################\n\n";
			print color 'reset';

			
			print "    Run PLINK for MDS plot?                      \t(press return to leave at $do_mds)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_mds = "y"};
			if ($answer eq "n"){$do_mds = "n"};

			
			print "    Run PLINK for QQ plot?                       \t(press return to leave at $do_qq_plot)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_qq_plot = "y"};
			if ($answer eq "n"){$do_qq_plot = "n"};
			
			
			print "    Run EMMAX analysis?                           \t(press return to leave at $do_emmax)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_emmax = "y"};
			if ($answer eq "n"){$do_emmax = "n"};
			
			
			print "    Do IBS test in PLINK?             	          \t(press return to leave at $do_ibs_test)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_ibs_test = "y"};
			if ($answer eq "n"){$do_ibs_test = "n"};

			
			print "    Run Eigenstrat analysis?                       \t(press return to leave at $do_eigenstrat)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_eigenstrat = "y"};
			if ($answer eq "n"){$do_eigenstrat = "n"};


			print "\n\n";
			print color 'yellow';
			print "##########################\n";
			print "# Other analysis options #\n";
			print "##########################\n\n";
			print color 'reset';

			print "    Make files for 'Homozygosity Mapping'?          \t(press return to leave at $do_homozygosity)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_homozygosity = "y"};
			if ($answer eq "n"){$do_homozygosity = "n"};

			print "    Do PLINK missingness analysis?                  \t(press return to leave at $do_missingness)  ";
			$answer = <STDIN>;
			chomp $answer;
			$answer = lc $answer;
			if ($answer eq "y"){$do_missingness = "y"};
			if ($answer eq "n"){$do_missingness= "n"};

			print "\n\n";
			print color 'yellow';
			print "###################################\n";
			print "# Extra command line options      #\n";
			print "###################################\n\n";
			print color 'reset';

			print "    Add extra command line options?                  \t(press return to leave at $plink_extra_string)  ";

			print "    Type in extra command line options, e.g. --ci 0.95\n\n";
			print "     >  ";
			$answer = <STDIN>;
			chomp $answer;

			if ($answer ne ""){$plink_extra_string = $answer};
			if ($answer eq ""){$plink_extra_string= ""}

		}
		else
		{
			$choice_ok = "true";
		}

	} # while ($choice_ok eq "false")

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

	print "\nWould you like to change these? (y/n)    ";

	$answer = <STDIN>;
	chomp $answer;

	if ($answer eq "y" || $answer eq "Y")
	{
		print "\n";
		$parameters_ok = "false";
	
		print "    --maf value [lower value keeps more SNPs] (press return to leave at $maf_value):   ";
		$answer = <STDIN>;
		chomp $answer;
		if ($answer ne ""){$maf_value = $answer;}

		print "    --geno value [higher value keeps more SNPs] (press return to leave at $geno_value):   ";
		$answer = <STDIN>;
		chomp $answer;
		if ($answer ne ""){$geno_value = $answer;}

		print "    --mind value [higher value keeps more individuals] (press return to leave at $mind_value):   ";
		$answer = <STDIN>;
		chomp $answer;
		if ($answer ne ""){$mind_value = $answer;}

		print "    --mperm value (press return to leave at $mperm_value):   ";
		$answer = <STDIN>;
		chomp $answer;
		if ($answer ne ""){$mperm_value = $answer;}


		if (($do_hap_assoc eq "y") || ($do_hap_assoc_mperm eq "y"))
		{
			print "    Smallest haplotype value (press return to leave at $start_haplotype):   ";
			$answer = <STDIN>;
			chomp $answer;
			if ($answer ne ""){$start_haplotype = $answer;}

			print "    Largest haplotype value (press return to leave at $end_haplotype):    ";
			$answer = <STDIN>;
			chomp $answer;
			if ($answer ne ""){$end_haplotype = $answer;}

			if ($do_hap_assoc_mperm eq "y")
			{
				print "    --mperm value for haplotyping (press return to leave at $mperm_value_haplotyping):   ";
				$answer = <STDIN>;
				chomp $answer;
				if ($answer ne ""){$mperm_value_haplotyping = $answer;}

			}
		} # if (($do_hap_assoc eq "y") || ($do_hap_assoc_mperm eq "y"))

	}
	else
	{
		$parameters_ok = "true";
	}

} # end of ask about parameters


##############################
# Ask about phenotype file   #
##############################
 while ($phenotypes_ok eq "false")
 {

	print "\n\n";
	print color 'yellow';
	print "###########################\n";
	print "# External phenotype file #\n";
	print "###########################\n";
	print color 'reset';

	print "\nThese are the external phenotype file settings:    \n\n";
	print "    --pheno      (name of phenotype file)    \t$pheno_file\n";
	print "    --pheno-name (name of phenotype column)  \t$pheno_name\n";

	print "\nWould you like to change these? (y/n)    ";

	$answer = <STDIN>;
	chomp $answer;

	if (lc($answer) eq "y")
	{
		print "\n";
		$phenotypes_ok = "false";
		$pheno_file_input = "ls";

		if ($pheno_file ne "ls"){$last_pheno_file = $pheno_file;}
		$pheno_file = "none";

		until (-e $pheno_file)
		{
			print "    External phenotype file (--pheno) (press return to leave at $last_pheno_file):   ";
			$pheno_file_input = <STDIN>;
			chomp $pheno_file_input;

			if ($pheno_file_input ne ""){$pheno_file = $pheno_file_input}
			if ($pheno_file_input eq ""){$pheno_file = $last_pheno_file}

			if ($pheno_file eq "ls"){print "\n";system ("ls *.txt");print "\n"}
			if ($pheno_file ne "ls"){if (! -e $pheno_file){print "\n\n>>>>>>>>  File $pheno_file not found.  Try again.  <<<<<<<<\n\n";}}
		}


		#####################################
		# Look at columns in phenotype file #
		#####################################
		if (-e $pheno_file)
		{
			open(PHENO,"$pheno_file") or die "Can't open file $pheno_file";

			@pheno_file_array = <PHENO>;

			$no_of_samples = (scalar @pheno_file_array) - 1;

			$single_line = $pheno_file_array[0];
			@item=split(/\s+/,$single_line);

			$no_of_phenotypes = (scalar @item) - 2;

			print "Phenotype file $pheno_file:\n\n";

			print "No. of samples:     \t$no_of_samples\n";
			print "No. of phenotypes:  \t$no_of_phenotypes\n\n";

			&print_message("Which of these phenotypes do you want to use for this run [$prefix"."_"."$analysis_suffix]?","input");

			for ($phenotype_count = 1; $phenotype_count <= $no_of_phenotypes; $phenotype_count++)
			{
				print "  <$phenotype_count> \t$item[$phenotype_count + 1]\n";
			}

			print "\n> ";
			$pheno_column = <STDIN>;
			chomp $pheno_column;

			$pheno_name = $item[$pheno_column + 1];
		}
	}
	else
	{
		$phenotypes_ok = "true";
	}

	$pheno_file_string = "--pheno $pheno_file";
	$pheno_name_string = "--pheno-name $pheno_name";

	if (($pheno_file eq "") || ($pheno_file eq "none")){$pheno_file_string = ""; $pheno_name_string = "";}

	if ($pheno_file_string eq ""){$pheno_name_string = ""}

	$plink_pheno_string = $pheno_file_string." ".$pheno_name_string;

	print "plink_pheno_string: $plink_pheno_string\n\n";

	# If the plink_pheno_string is a space, make it nothing "" #
	if ($plink_pheno_string eq " "){$plink_pheno_string = ""}

 } # while ($phenotypes_ok eq "false")


#################################################################
# Open overall logfile                                          #
#################################################################
open(LOGALL,">$logfile_all") or die "Can't open file $logfile_all";



#################################################################
# Create the results folder                                     #
#################################################################
$directory_suffix = 0;
while (-e $results_folder)
{
	$directory_suffix = $directory_suffix + 1;
	$results_folder = $prefix_with_suffix."_results"."_".$directory_suffix;
}

&run_unix_command("mkdir $results_folder","Create the results directory");


#################################################################
# open Command Log file                                         #
#################################################################
$command_log = "$prefix_with_suffix"."_GWAS_analyses_command_log.out";
open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

print COMMAND_LOG "#########################################\n";
print COMMAND_LOG "COMMAND LOG for perl script GWAS_analyses\n";
print COMMAND_LOG "#########################################\n\n";
print COMMAND_LOG "  Version $version\n\n";

print COMMAND_LOG "PLINK prefix:          \t$prefix\n";
print COMMAND_LOG "Filename with suffix:  \t$prefix_with_suffix\n\n";

if ($species eq "dog") {print COMMAND_LOG "Assembly is $assembly\n\n";}

print COMMAND_LOG scalar localtime;
print COMMAND_LOG "\n\n";

print COMMAND_LOG "Parameters in place:\n\n";
print COMMAND_LOG "--maf $maf_value\n";
print COMMAND_LOG "--geno $geno_value\n";
print COMMAND_LOG "--mind $mind_value\n";
print COMMAND_LOG "--mperm $mperm_value\n";
print COMMAND_LOG "--mperm $mperm_value_haplotyping (if used for mperm haplotyping)\n\n";

if (($do_hap_assoc eq "y") || ($do_hap_assoc_mperm eq "y"))
{
	print COMMAND_LOG "Haplotypes from $start_haplotype to $end_haplotype\n\n";
}
print COMMAND_LOG "Phenotype string for PLINK:           \t$plink_pheno_string\n\n";
print COMMAND_LOG "Extra command line options for PLINK: \t$plink_extra_string\n\n";

print "\n";

system("cp $command_log /home/genetics/command_logs/$command_log");

#################################################################
#                                                               #
#                     Now run PLINK   ( start )                 #
#                                                               #
#################################################################


###########################################################################
# General check that the data is working and that the filters are correct #
###########################################################################

print color 'bold yellow';
print "\n\n";
print "############################################################\n";
print "#            RUNNING GENERAL CHECK ON THE DATA             #\n";
print "#                                                          #\n";
print "#         Please CHECK when this has run that the          #\n";
print "#         values you used for geno, maf and mind           #\n";
print "#         have not filtered off too many SNPs or           #\n";
print "#         samples.                                         #\n";
print "############################################################\n\n";
print color 'reset';

&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --maf $maf_value --geno $geno_value --mind $mind_value --file $prefix $plink_pheno_string --make-founders --nonfounders --out temp","First quick check on the data");

print color 'bold yellow';

print "\n";

if ($plink_pheno_string eq "")
{
	print "################################################################\n";
	print "#                 GENERAL CHECK ON THE DATA                    #\n";
	print "#                                                              #\n";
	print "# Please check at this point that the values you used for GENO #\n";
	print "# MIND and MAF have not filtered off too many SNPs or samples. #\n";
	print "################################################################\n";
}

if ($plink_pheno_string ne "")
{
	print "################################################################\n";
	print "#                 GENERAL CHECK ON THE DATA                    #\n";
	print "#                                                              #\n";
	print "# Please check at this point that the values you used for GENO #\n";
	print "# and MAF have not filtered off too many SNPs or samples.      #\n";
	print "#                                                              #\n";
	print "#              ------------------------------                  #\n";
	print "#                                                              #\n";
	print "# ALSO SINCE YOU HAVE USED AN EXTERNAL PHENOTYPE FILE, MAKE    #\n";
	print "# THE NUMBERS OF CASES AND CONTROLS LOOK RIGHT!                #\n";
	print "################################################################\n";
}

print "\nPress return to continue.\n";
$answer = <STDIN>;


print color 'reset';

if ($plink_pheno_string ne "")
{
	#####################################################
	# Recode the ped file to include any phenotype data #
	#####################################################
	print_message("Recoding file using command line option $plink_pheno_string","message");

	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --file $prefix $plink_pheno_string --make-founders --nonfounders --recode --out $prefix_with_suffix","Recode file using command line option $plink_pheno_string");

	&print_message("Recoded file using command line option $plink_pheno_string","message");
}

if ($plink_pheno_string eq "")
{
	################################################################################################################
	# No need to recode the ped file to include any phenotype data - but make a copy with the new name plus suffix #
	################################################################################################################

	if ($split_x eq "yes")
	{
		&print_message("Making a copy of PED and MAP files to the new name $prefix_with_suffix_unsplit ");
		&run_unix_command("cp $prefix".".map $prefix_with_suffix_unsplit".".map","Make copy of MAP file with run suffix");
		&run_unix_command("cp $prefix".".ped $prefix_with_suffix_unsplit".".ped","Make copy of PED file with run suffix");
	}
	if ($split_x eq "no")
	{
		&print_message("Making a copy of PED and MAP files to the new name $prefix_with_suffix ");
		&run_unix_command("cp $prefix".".map $prefix_with_suffix".".map","Make copy of MAP file with run suffix");
		&run_unix_command("cp $prefix".".ped $prefix_with_suffix".".ped","Make copy of PED file with run suffix");
	}

	&print_message("Made copy of PED and MAP file to the new name $prefix_with_suffix_unsplit ");
}

###########################################################################
# Split X and Y chromosomes into X, Y and XY                              #
###########################################################################

if ($split_x ne "no")
{
	&print_message("Splitting X and Y into X, Y and XY","message");

	$split_x_command_log = "$prefix_with_suffix_unsplit"."_split_x_command_log.out";

	&run_unix_command("perl /home/genetics/scripts/run_plink2_split_x.pl --file $prefix_with_suffix_unsplit","Run perl script run_plink2_split_x.pl");

	&move_to_results_folder("$split_x_command_log");
}


#########################################################
# Check on the missingness of the data                  #
#########################################################

if ($do_missingness eq "y")
{
	&print_message("Checking on the missingness of the data");

	if ($use_plink2_if_poss eq "no")
	{&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --file $prefix_with_suffix $plink_pheno_string --nonfounders --make-founders --test-missing --out $prefix_with_suffix","Run --test-missing on the data");}

	if ($use_plink2_if_poss eq "yes")
	{&run_unix_command("plink2 $species_string_plink --allow-no-sex --file $prefix_with_suffix $plink_pheno_string --nonfounders --make-founders --test-missing --out $prefix_with_suffix","Run --test-missing on the data using PLINK2");}

	# Add logfile to logfile_all
	&add_to_logfile;

	&move_to_results_folder (".missing");

	&print_message("Checked on the missingness of the data");
}


#########################################################
# Convert the .map and .ped files in binary files       #
#########################################################
&print_message("Converting files to binary files");

if ($use_plink2_if_poss eq "no")
{
	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --file $prefix_with_suffix $plink_pheno_string --make-founders --nonfounders --make-bed  --maf $maf_value --geno $geno_value --mind $mind_value --out $prefix_with_suffix","Convert map and ped files into binary file (also filter with MAF, GENO and MIND)");
}

if ($use_plink2_if_poss eq "yes")
{
	&run_unix_command("plink2 $species_string_plink --allow-no-sex --file $prefix_with_suffix $plink_pheno_string --make-founders --nonfounders --make-bed  --maf $maf_value --geno $geno_value --mind $mind_value --out $prefix_with_suffix","Convert map and ped files into binary file (also filter with MAF, GENO and MIND)");
}

&add_to_logfile;

&print_message("Converted files to binary files");


################################################################################
# Run 'which_samples_were_used' to make a list of samples in the .ped file     #
################################################################################
if ($do_which_samples eq "y")
{
	&print_message("Running the perl script 'which_samples_were_used' to keep a record of which samples were in the PED file");

	&run_unix_command("perl /home/genetics/scripts/which_samples_were_used.pl -f $prefix_with_suffix",
	"Run which_samples_were_used to keep a record of which samples were in the .ped file");

	&move_to_results_folder (".samples_used.out");

	&print_message("Ran the perl script 'which_samples_were_used' to keep a record of which samples were in the PED file");
}


################################################################################
# Run standard association analysis with no filtering                          #
################################################################################

if ($do_assoc eq "y")
{
	&print_message("Running standard association analysis - unfiltered");

	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --file $prefix_with_suffix $plink_pheno_string --assoc --maf 0 --out $prefix_with_suffix"."_unfiltered",
	"Run standard association analysis - unfiltered");

	# Add logfile to logfile_all
	&add_to_logfile;

	&move_to_results_folder ("_unfiltered.assoc");

	&move_to_results_folder (".model");

	&delete_unwanted ("_unfiltered.log");
	&delete_unwanted ("_unfiltered.hh");

	&print_message("Ran standard association analysis - unfiltered");
	
}


#########################################################
# Run standard association analysis                     #
#########################################################

if ($do_assoc eq "y")
{
	&print_message("Running standard association analysis");

	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --assoc --out $prefix_with_suffix",
	"Run standard association analysis");

	# Add logfile to logfile_all
	&add_to_logfile;

	&move_to_results_folder (".assoc");

	&print_message("Ran standard association analysis");
}

#########################################################
# Run sex check                                         #
#########################################################

if ($do_sex_check eq "y")
{
	&print_message("Running sex check");

	if ($use_plink2_if_poss eq "yes")
	{
		&run_unix_command("plink2 $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --check-sex ycount 0.5 0.8 --out $prefix_with_suffix",
	"Run sex check with PLINK2");
	}

	if ($use_plink2_if_poss eq "no")
	{
		&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --check-sex --out $prefix_with_suffix",
	"Run sex check");
	}
	

	# Add logfile to logfile_all
	&add_to_logfile;

	&move_to_results_folder (".sexcheck");

	&print_message("Ran sex check");
}



#########################################################
# Produce QQ-plot output                                #
#########################################################

if ($do_qq_plot eq "y")
{
	&print_message("Producing 'assoc.adjusted' file with data for QQ plot");

	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --assoc --adjust --qq-plot --log10 --out $prefix_with_suffix","Run PLINK to get data for QQ plot");

	# Add logfile to logfile_all
	&add_to_logfile;

	&move_to_results_folder (".assoc");
	&move_to_results_folder (".assoc.adjusted");

	&print_message("Produced 'assoc.adjusted' file with data for QQ plot");
}


#########################################################
# Produce MDS data to look at clustering                #
# (doesn't do it for any particular phenotype)          #
#########################################################

if ($do_mds eq "y")
{
	if ($use_plink2_if_poss eq "no")
	{
		&print_message("Running MDS analysis to look at clustering");

		&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --bfile $prefix_with_suffix $plink_pheno_string --mds-plot 2 --out $prefix_with_suffix","Run MDS (multi-dimensional-scaling) analysis corrected for multiple testing to look at clustering");

		&add_to_logfile;

		&move_to_results_folder (".mds");

		&print_message("Ran MDS analysis to look at clustering");
	}
	if ($use_plink2_if_poss eq "yes")
	{
		&print_message("Running MDS analysis with PLINK2 to look at clustering");

		&run_unix_command("plink2  $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --cluster --mds-plot 2 --out $prefix_with_suffix","Run MDS (multi-dimensional-scaling) analysis corrected for multiple testing to look at clustering");

		&add_to_logfile;

		&move_to_results_folder ($prefix_with_suffix.".mds");

		&print_message("Ran MDS analysis with PLINK2 to look at clustering");

		##########################
		# Deleted unwanted files #
		##########################
		&delete_unwanted ("_".$haplotype_count.".log");
		&delete_unwanted ("_".$haplotype_count.".hh");
	}

	&delete_unwanted (".cluster0");
	&delete_unwanted (".cluster1");
	&delete_unwanted (".cluster2");
	&delete_unwanted (".cluster3");
}


############################################################
# Run association corrected for multiple testing  (MPERM)  #
############################################################

if ($do_assoc_mperm eq "y")
{
	if ($use_plink2_if_poss eq "no")
	{
		print_message("Running association analysis corrected for multiple testing using MPERM");

		&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --bfile $prefix_with_suffix $plink_pheno_string --assoc --mperm $mperm_value --out $prefix_with_suffix","Run association analysis corrected for multiple testing");

		&add_to_logfile;

		&move_to_results_folder (".assoc.mperm");

		&print_message("Running standard association analysis corrected for multiple testing using MPERM with PLINK2");
	}
	if ($use_plink2_if_poss eq "yes")
	{
		print_message("Running association analysis corrected for multiple testing using MPERM");

		&run_unix_command("plink2 $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --model mperm=$mperm_value --out $prefix_with_suffix","Run association analysis corrected for multiple testing with PLINK2");

		&add_to_logfile;

		#&move_to_results_folder (".model.best.mperm");

		&move_to_results_folder ($prefix_with_suffix.".model.best.mperm");

		&print_message("Ran standard association analysis corrected for multiple testing using MPERM with PLINK2");
	}																
}


################################################################################
# Run plink2homozygosity to create files for Homozygosity Mapping spreadsheet  #
################################################################################

if ($do_homozygosity eq "y")
{
	&print_message("Running plink2homozygosity to create files for Homozygosity Mapping spreadsheet");

	&run_unix_command("perl /home/genetics/scripts/plink2homozygosity.pl -f $prefix_with_suffix","Run PLINK2HOMOZYGOSITY to create files for Homozygosity Mapping spreadsheet");

	&move_to_results_folder ("$prefix_with_suffix"."_affected.txt");
	&move_to_results_folder ("$prefix_with_suffix"."_normal.txt");

	&print_message("Ran plink2homozygosity to create files for Homozygosity Mapping spreadsheet");
}


############################################################
# Only include all models if analysis scope is set to full #
############################################################

if ($do_models eq "y")
{
	#########################################################
	# Run association with all models (dom,rec etc)         #
	#########################################################
	&print_message("Looking at different models (dominant, recessive etc)");

	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --bfile $prefix_with_suffix $plink_pheno_string --model --cell 0 --out $prefix_with_suffix","Run PLINK to look at different models (dom, rec etc)");

	&add_to_logfile;
	&move_to_results_folder (".model");

	&print_message("Looked at different models (dominant, recessive etc)");

	
	#########################################################
	# Run association mperm with best model                 #
	#########################################################
	&print_message("Looking at best model with mperm");

	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --model --cell 0 --mperm $mperm_value --out $prefix_with_suffix",
	"Run PLINK to look at best model with mperm");

	&add_to_logfile;
	&move_to_results_folder (".model");
	&move_to_results_folder (".model.best.mperm");

	&print_message("Looked at best model with mperm");


	#########################################################
	# Run mperm association with the dominant model         #
	#########################################################
	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --bfile $prefix_with_suffix $plink_pheno_string --model-dom --cell 0 --mperm $mperm_value --out $prefix_with_suffix",
	"Run PLINK to look at dominant model with mperm");

	&add_to_logfile;
	&move_to_results_folder (".model");
	&move_to_results_folder (".model.dom.mperm");

	&print_message("Looked at dominant model with mperm");


	#########################################################
	# Run mperm association with the recessive model        #
	#########################################################
	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --model-rec --cell 0 --mperm $mperm_value --out $prefix_with_suffix",
	"Run PLINK to look at recessive model with mperm");

	&add_to_logfile;
	&move_to_results_folder (".model");
	&move_to_results_folder (".model.rec.mperm");

	&print_message("Looked at recessive model with mperm");


	#########################################################
	# Run mperm association with the genotypic model        #
	#########################################################
	&print_message("Looking at genotypic model with mperm","message");

	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --bfile $prefix_with_suffix $plink_pheno_string --model-gen --cell 0 --mperm $mperm_value --out $prefix_with_suffix",
	"Run PLINK to look at genotypic model with mperm");

	&add_to_logfile;
	&move_to_results_folder (".model");
	&move_to_results_folder (".model.gen.mperm");

	&print_message("Looked at genotypic model with mperm","message");


	#########################################################
	# Run mperm association with the trend model            #
	#########################################################
	&print_message("Looking at trend model with mperm","message");

	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --model-trend --cell 0 --mperm $mperm_value --out $prefix_with_suffix",
	"Run PLINK to look at trend model with mperm");
	
	&add_to_logfile;
	&move_to_results_folder (".model");
	&move_to_results_folder (".model.trend.mperm");

	&print_message("Looked at trend model with mperm","message");

} # do_models = y



#########################################################
# IBS TEST in PLINK                                     #
#########################################################
if ($do_ibs_test eq "y")
{
	 
	 if ($use_plink2_if_poss eq "no")
	 {
	 	&print_message("Running ibs-test in PLINK","message");

	 	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --file $prefix_with_suffix $plink_pheno_string --nonfounders --make-founders --ibs-test --out $prefix_with_suffix",
	 "Ran ibs-test in PLINK");

	 	&print_message("Ran ibs-test in PLINK");
	}

	if ($use_plink2_if_poss eq "yes")
	 {
	 	&print_message("Running ibs-test in PLINK2","message");

	 	&run_unix_command("plink2 $species_string_plink --allow-no-sex --make-founders --nonfounders --file $prefix_with_suffix $plink_pheno_string --nonfounders --make-founders --ibs-test --out $prefix_with_suffix",
	 "Ran ibs-test in PLINK");

	 	&print_message("Ran ibs-test in PLINK2","message");
	}

	####################################
	# Rename log file to ibs-test file #
	####################################
	 $command = "mv $prefix_with_suffix".".log $prefix_with_suffix".".ibs-test"; 
	 system("$command");

	 &move_to_results_folder (".ibs-test");

} # ibs-test


################################################################
# EMMAX analysis                                               #
################################################################

if ($do_emmax eq "y")
{
	#############################
	# Run PERL script run_emmax #
	#############################
	print_message("Running EMMAX to produce CHISQs adjusted for stratification","message1");

	&run_unix_command("perl /home/genetics/scripts/run_emmax.pl $species_string_perl --file $prefix_with_suffix $plink_pheno_string --maf $maf_value --geno $geno_value --mind $mind_value",
	"Ran EMMAX to produce CHISQs adjusted for stratification");

	&print_message("Ran EMMAX to produce CHISQs adjusted for stratification","message");
		
	# Rename ps file to emmax_pvalues.out
	$command = "mv $prefix_with_suffix".".ps $prefix_with_suffix"."_emmax_pvalues.out"; 
	system("$command");
	 
	&move_to_results_folder (".emmax_log.out");
	&move_to_results_folder (".emmax_out");
	&move_to_results_folder (".hd_emmax_out");
	&move_to_results_folder (".emmax_qq");
	&move_to_results_folder (".hBN.kinf");
	&move_to_results_folder (".hIBS.kinf");
	
	&move_to_results_folder ("$prefix_with_suffix"."_emmax_qq.pdf");
	
	# Delete intermediate files
	&delete_file("$prefix_with_suffix".".tped");
	&delete_file("$prefix_with_suffix".".tfam");
	&delete_file("$prefix_with_suffix".".pheno");
	&delete_file("$prefix_with_suffix".".reml");
	&delete_file ("$prefix_with_suffix"."_emmax_qq.ps");
	&delete_file ("gnufile.xtxt");
	#&delete_file("$prefix_with_suffix"."_emmax_pvalues.out");
	
	 
} # end of if do EMMAX = y


#############################################################
# Only do plink clustering if analysis scope is set to full #
#############################################################

if ($do_plink_clustering eq "y")
{
	#########################################################
	# Allow PLINK to cluster the samples                    #
	#########################################################
	&print_message("Allowing PLINK to cluster the samples","message");

	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --cluster --cc --ppc 0.01 --out $prefix_with_suffix",
	"Run PLINK to produce clusters with at least one case and one control and ppc=0.01 (read PLINK documentation)");

	&add_to_logfile;

	&print_message("Allowed PLINK to cluster the samples","message");


	#########################################################
	# Run association only within these clusters            #
	#########################################################
	&print_message("Looking for association within these clusters","message");

	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --mh --within $prefix_with_suffix".".cluster2 --out $prefix_with_suffix",
	"Run MH association within the cluster2 file");

	&add_to_logfile;
	&move_to_results_folder (".cmh");

	&print_message("Looked for association within these clusters","message");


	###########################################################
	# Run association only within these clusters using mperm  #
	###########################################################
	&print_message("Looking for association within these clusters using Mperm","message");

	&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --mh --mperm $mperm_value --within $prefix_with_suffix".".cluster2 --out $prefix_with_suffix",
	"Run MH association within the cluster2 file using mperm");

	&add_to_logfile;
	&move_to_results_folder (".cmh");
	&move_to_results_folder (".cmh.mperm");
	&move_to_results_folder (".cluster2");

	&delete_unwanted (".cluster0");
	&delete_unwanted (".cluster1");
	&delete_unwanted (".cluster2");
	&delete_unwanted (".cluster3");

	&print_message("Looked for association within these clusters using Mperm","message");

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
			&print_message("Running Haplotyping Assoc analysis  (window = $haplotype_count)","message");

			&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --hap-window $haplotype_count --hap-assoc --out $prefix_with_suffix"."_"."$haplotype_count",
			"Run Haplotyping Assoc analysis on window of $haplotype_count");

			#plink $plink_extra_string --noweb --dog --bfile NAME --hap-window 4 --hap-assoc --out NAME_w4

			&add_to_logfile;
			
			&print_message("Ran Haplotyping Assoc analysis  (window = $haplotype_count)","message");


			########################################################
			# Convert impute.map and impute.assoc.mperm files to   #
			# a single file which can be plotted in Hap_Mperm_plot #
			########################################################
			&print_message("Converting output files with hap_assoc_convert  (window = $haplotype_count)","message");

			if ($species ne ""){$species_string_perl eq "-species $species"} else {$species_string_perl eq ""}

			&run_unix_command("perl /home/genetics/scripts/hap_assoc_convert.pl -file $prefix_with_suffix"."_"."$haplotype_count -hap_length $haplotype_count $species_string_perl -map $prefix_with_suffix.map",
			"Converted assoc.hap file to a _hap_assoc_{n}.txt file which can be plotted in 'Haplotype Association Plot'");
			
			
			&run_unix_command("perl /home/genetics/scripts/hap_mperm_convert.pl -file $prefix_with_suffix -hap_length $haplotype_count $species_string_perl -map $prefix_with_suffix.map",
			"Converted impute.map and impute.assoc.mperm files to single _hap_mperm_{n}.txt file which can be plotted in 'Haplotype Association Plot'");

			
			&move_to_results_folder ("_".$haplotype_count.".assoc.hap");
			&move_to_results_folder ("_".$haplotype_count."_hap_assoc_{".$haplotype_count."}.txt");

			&print_message("Converted output files with hap_assoc_convert  (window = $haplotype_count)","message");


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
			&print_message("Imputing the haplotypes (window = $haplotype_count)","message");

			&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix $plink_pheno_string --hap-window $haplotype_count --hap-impute --out $prefix_with_suffix"."_"."$haplotype_count",
			"Run Haplotyping MPERM analysis on window of $haplotype_count");

			&add_to_logfile;
			
			&print_message("Imputed the haplotypes (window = $haplotype_count)","message");


			#######################################################################
			# Haplotype analysis  - make binary files and delete non-binary files #
			#######################################################################
			&print_message("Making binary files for MPERM haplotyping  (window = $haplotype_count)","message");


			&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --file $prefix_with_suffix"."_".$haplotype_count.".impute $plink_pheno_string --make-bed --out $prefix_with_suffix"."_".$haplotype_count.".impute",
			"Made binary files for MPERM haplotyping on window of $haplotype_count");

			&add_to_logfile;
			
			&print_message("Made binary files for MPERM haplotyping  (window = $haplotype_count)","message");

			
			##################################
			# Delete the non-binary versions #
			##################################
			&delete_unwanted ("_".$haplotype_count.".impute.ped");
			
			
			#######################################################
			# Haplotype analysis  - association with permutation  #
			#######################################################
			&print_message("RUNNING Haplotyping MPERM analysis to use in Hap_Mperm_Plot (window = $haplotype_count)","message");

			&run_unix_command("plink $plink_extra_string --noweb $species_string_plink --allow-no-sex --make-founders --nonfounders --bfile $prefix_with_suffix"."_".$haplotype_count.".impute $plink_pheno_string --assoc --mperm $mperm_value_haplotyping --out $prefix_with_suffix"."_".$haplotype_count.".impute",
			"Ran Haplotype MPERM analysis to use in Hap_Mperm_Plot");
			
			&add_to_logfile;
			
			&print_message("Ran Haplotyping MPERM analysis to use in Hap_Mperm_Plot (window = $haplotype_count)","message");


			########################################################
			# Convert impute.map and impute.assoc.mperm files to   #
			# a single file which can be plotted in Hap_Mperm_plot #
			########################################################

			if ($species ne ""){$species_string_perl eq "-species $species"} else {$species_string_perl eq ""}

			print_message("Converting output files with hap_mperm_convert  (window = $haplotype_count)","message");

			&run_unix_command("perl /home/genetics/scripts/hap_mperm_convert.pl -file $prefix_with_suffix -hap_length $haplotype_count $species_string_perl -map $prefix_with_suffix.map",
			"Converted impute.map and impute.assoc.mperm files to single _hap_mperm_{n}.txt file which can be plotted in 'Haplotype Association Plot'");

			
			&move_to_results_folder ("_".$haplotype_count.".impute.map");
			&move_to_results_folder ("_".$haplotype_count.".impute.assoc.mperm");
			&move_to_results_folder ("_hap_mperm_{".$haplotype_count."}.txt");

			&print_message("Converted output files with hap_mperm_convert  (window = $haplotype_count)","message");


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



################################################################
# Eigenstrat analysis                                          #
################################################################

if ($do_eigenstrat eq "y")
{
	 #########################################################
	 # Convertf - to convert PLINK files to Eigenstrat files #
	 #########################################################
	 &print_message("Running CONVERTF to convert PLINK files to Eigenstrat files","message");

	 &run_unix_command("perl /home/genetics/scripts/runconvertf_".$species.".pl -f $prefix_with_suffix",
	 "Run convertf to convert PLINK files to Eigenstrat files");

	 &print_message("Ran CONVERTF to convert PLINK files to Eigenstrat files","message");

	 
	 #########################################################
	 # Smartpca - to carry out PCA analysis                  #
	 #########################################################
	 &print_message("Running SMARTPCA to carry out Principal Components Analysis","message");

	 &run_unix_command("perl /home/genetics/scripts/runsmartpca_".$species.".pl -f $prefix_with_suffix",
	 "Ran SMARTPCA to carry out Principal Components Analysis");
		
	 &print_message("Ran SMARTPCA to carry out Principal Components Analysis","message");

	 
	 #########################################################
	 # Eigenstrat - to account for structure                 #
	 #########################################################
	 print_message("Running EIGENSTRAT to produce CHISQs adjusted for stratification","message");

	 &run_unix_command("perl /home/genetics/scripts/runeigenstrat_".$species.".pl -f $prefix_with_suffix",
	 "Ran EIGENSTRAT to produce CHISQs adjusted for stratification");


	 # Rename log file to eigenlog
	 $command = "mv $prefix_with_suffix".".log $prefix_with_suffix".".eigenlog"; 
	 system("$command");

	 &move_to_results_folder (".eigenlog");
	 &move_to_results_folder ("$prefix_with_suffix"."_PCA_plot.pdf");
	 &move_to_results_folder (".smartpca_log");
	 &move_to_results_folder (".eigenstrat.out");
	 &move_to_results_folder (".chisq");

	 &print_message("Ran EIGENSTRAT to produce CHISQs adjusted for stratification","message");
	 
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


&print_message_command_log("Move log files etc to the Results folder","message");

&move_to_results_folder ("$logfile_all");



#################################################################
# Write the final README file to tell user which files are what #
#################################################################


open(README,">$final_readme_file") or die "Can't open file $final_readme_file";

print README "################################################\n";
print README "# List of the output files from GWAS_analyses   #\n";
print README "# and some information on what to do with them #\n";
print README "################################################\n\n";

print README "(Please also read any instructions in the Genotyping Analysis/Instructions Folder)\n\n";

print README "Filename: \t$prefix_with_suffix\n";
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
	print README " > ".$prefix_with_suffix.".missing\n\n";
	print README "Open this file in Excel and sort by the P-value column.\n";
	print README "This tells you if there is a significance difference between cases and controls for any SNP.\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}


# Sex check #
if ($do_sex_check eq "y")
{
	 print README "The following file has the sex check data.\n";
	 print README "Look at this to see if there any samples that look as if they might be the wrong sex\n\n";
	 print README "If PLINK2 has been used, look especially at the column counting the number of Y SNPs detected\n\n";
	 print README " > ".$prefix_with_suffix.".sexcheck\n\n";
	 print README "-----------------------------------------------------------------------------------------\n";
}


# Basic association #
if ($do_assoc eq "y")
{
	 print README "The following file has the simple association data.\n";
	 print README "It can be plotted using the Excel sheet PLINK_plot\n\n";
	 print README " > ".$prefix_with_suffix.".assoc\n\n";
	 print README "-----------------------------------------------------------------------------------------\n";
}


# MPERM association #
if ($do_assoc_mperm eq "y")
{
	 print README "The following file has the association data, taking into account multiple testing.\n";
	 print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";
	 print README " > ".$prefix_with_suffix.".assoc.mperm\n\n";
	 print README "-----------------------------------------------------------------------------------------\n";
}


# Clustering #
if ($do_mds eq "y")
{
	print README "The following MDS file has the results of a clustering analysis (multi-dimensional scaling).\n";
	print README "If you plot this in Excel you can visualise any clusters.\n";
	print README "(Do a scatter plot, and plot a separate series for cases and controls).\n\n";
	print README " > ".$prefix_with_suffix.".mds\n\n";
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

	print README " > ".$prefix_with_suffix.".model\n\n";

	print README "-----------------------------------------------------------------------------------------\n";


	# Mperm association with the best model#
	print README "The following file has the association data for the best of five different models at each SNP.\n";
	print README "Normal (allelic), Recessive, Dominant, Trend and Genotypic.\n";
	print README "(This also takes into account multiple testing by using mperm)\n";
	print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";
	print README " > ".$prefix_with_suffix.".model.best.mperm\n\n";
	print README "-----------------------------------------------------------------------------------------\n";


	# Mperm association with the dominant model#
	print README "The following file has the association data for the dominant model.\n";
	print README "(This also takes into account multiple testing by using mperm)\n";
	print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";
	print README " > ".$prefix_with_suffix.".model.dom.mperm\n\n";
	print README "-----------------------------------------------------------------------------------------\n";

	# Mperm association with the recessive model#
	print README "The following file has the association data for the recessive model.\n";
	print README "(This also takes into account multiple testing by using mperm)\n";
	print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";
	print README " > ".$prefix_with_suffix.".model.rec.mperm\n\n";
	print README "-----------------------------------------------------------------------------------------\n";

	# Mperm association with the genotypic model#
	print README "The following file has the association data for the genotypic model.\n";
	print README "(This also takes into account multiple testing by using mperm)\n";
	print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";
	print README " > ".$prefix_with_suffix.".model.gen.mperm\n\n";
	print README "-----------------------------------------------------------------------------------------\n";

	# Mperm association with the trend model#
	print README "The following file has the association data for the Cochran Armitage trend model.\n";
	print README "(This also takes into account multiple testing by using mperm)\n";
	print README "It can also be plotted using the Excel sheet PLINK_plot.\n\n";
	print README " > ".$prefix_with_suffix.".model.trend.mperm\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}

# Stratification #
if ($do_ibs_test eq "y")
{
	print README "The following file has the results of using the PLINK --ibs-test option.\n";
	print README "This shows how statistically different the cases and control populations are.\n\n";
	print README " > ".$prefix_with_suffix.".ibs_test\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}

# Clustering 2 #
if ($do_plink_clustering eq "y")
{
	print README "The following files are the results of allowing PLINK to cluster the samples and\n";
	print README "then carrying out a MH association only within clusters.\n";
	print README "This file containes the cluster data...\n\n";
	print README " > ".$prefix_with_suffix.".cluster2\n\n";
	print README "These files can be plotted in PLINK_plot.\n\n";
	print README " > ".$prefix_with_suffix.".cmh\n\n";
	print README " > ".$prefix_with_suffix.".cmh.mperm\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}



# Output for QQ-plot #
if ($do_qq_plot eq "y")
{
	print README "The following file can be used to produce a QQ-plot in Excel using QQ_plot.\n\n";
	print README " > ".$prefix_with_suffix.".assoc.adjusted\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}


if ($do_which_samples eq "y")
{
	# which_samples_were_used #
	print README "The following file keeps a list of the samples used in the .ped file, and their disease statuses:\n\n";
	print README " > ".$prefix_with_suffix.".samples_used.out\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}


# Homozygosity Mapping files #
if ($do_homozygosity eq "y")
{
	print README "The following files can be used in the Excel 'Homozygosity Mapping' spreadsheet\n";
	print README "and then in Excel sheets 'Haplotyping Display' and 'Analyse Critical Region'\n\n";
	print README "Open them in two separate worksheets called Affected and Normal\n";
	print README "(WARNING: This will not use any phenotypes which are in a separate file from the .ped file)\n\n";
	print README " > ".$prefix_with_suffix."_affected.txt\n";
	print README " > ".$prefix_with_suffix."_normal.txt\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}


if ($do_hap_assoc_mperm eq "y")
{

	# Haplotyping #
	print README "The following files have the results of the haplotype analysis with permutation.\n";
	print README "(This has been done for window sizes (N) = $start_haplotype-$end_haplotype)\n";
	print README "Use the Excel sheet Hap_Mperm_plot to visualise these results.\n\n";
	print README " > ".$prefix_with_suffix."_N.impute.map\n";
	print README " > ".$prefix_with_suffix."_N.impute.assoc.mperm\n\n";
	
	if ($species eq "dog")
	{
		print README "The files with numbers in curly brackets have been converted to work directly in PLINK_plot (or Hap_Mperm_plot) with out having to read it in with the button.\n\n";
		print README " > e.g. ".$prefix_with_suffix.".hap_mperm_{4}.txt\n\n"
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
	print README " > ".$prefix_with_suffix."_PCA_plot.pdf\n\n";
	print README "This file has information from the Eigensoft analysis, including statistical measures of the significance of the stratification.\n\n";
	print README " > ".$prefix_with_suffix.".smartpca_log\n\n";
	print README "This file has CHISQ association data, uncorrected and also corrected for stratification.\n";
	print README "(Use the Excel sheet Eigenstrat_plot to plot these two columns of data)\n\n";
	print README " > ".$prefix_with_suffix.".eigenstrat_out OR ".$prefix_with_suffix.".chisq\n\n";
	print README "-----------------------------------------------------------------------------------------\n";
}

#EMMAX #
if ($do_emmax eq "y")
{
	print README "The following files have the results of running EMMAX on the data.\n\n";
	print README "This is a mixed model approach to correct for population structure. \n";
	print README "(NOTE: We slightly prefer to use Fast Mixed Model, but the EMMAX results are very similar.)\n\n";
	
	print README "This file shows a QQ plot of the EMMAX-corrected data\n\n";
	print README " > ".$prefix_with_suffix."_emmax_qq.pdf\n\n";
	
	print README "This file has the data to plot your own QQ plots\n\n";
	print README " > ".$prefix_with_suffix.".emmax_qq\n\n";
	
	print README "This file has the adjusted EMMAX P-value data to produce an association plot after EMMAX correction\n\n";
	print README " > ".$prefix_with_suffix.".emmax_out\n\n";
	
	print README "There is also a version of this file for direct plotting in PLINK_plot:\n\n";
	print README " > ".$prefix_with_suffix.".hd_emmax_out\n\n";
	
	print README "These file have kinship data.  rPlease read the EMMAX documentation for full details.\n\n";
	print README " > ".$prefix_with_suffix.".hBN.kinf and ".$prefix_with_suffix.".hBN.kinf\n\n";
	
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



###########################################
# Remove some unwanted intermediate files #
###########################################

&print_message("Removing some unwanted files");

print COMMAND_LOG "\nRemoving some unwanted files...\n\n";

#remove plink files
&delete_unwanted (".log");
&delete_unwanted (".ind");
&delete_unwanted (".hh");
&delete_unwanted (".irem");

&delete_unwanted ("temp.hh");
&delete_unwanted ("temp.irem");
&delete_unwanted ("temp.log");

#move ped and map files to results folder
&print_message_command_log("Moving copies of the PED and MAP files to the results folder");
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
		print MAIL "Subject: GWAS_analyses on $prefix_with_suffix has finished\n\n";
		
		## Mail Body
		print MAIL "Your Genome Wide Association analysis on $prefix_with_suffix is complete\n\n";
		print MAIL "PERL script:       \tGWAS_analyses version $version\n\n";
		print MAIL "Results directory: \t$results_folder\n\n";
		print MAIL "Check the following files:\n\n";
		print MAIL "    1.  How to use the output files:        \t$final_readme_file\n";
		print MAIL "    2.  A list of all the unix commands:    \t$command_log\n";
		print MAIL "    3.  Program logs of the analyses:       \t$logfile_all\n\n";

		if($assembly eq "canfam2")
		{
			print  MAIL "=========================================\n";
			print  MAIL "Remember the output files are all canfam2\n";
			print  MAIL "=========================================\n\n";

			print  MAIL "To use the canfam3 assembly SNP positions you can either:\n\n";
			print  MAIL " - convert output assoc files using 'convert_canfam2_assoc_files'\n\n";
			print  MAIL " - convert MAP file using 'convert_canfam2_map_files' and re-analyse\n\n";
		}
	
	close(MAIL);
}

close README;
close COMMAND_LOG;


######################################################
# Copy command log to command log folder in genetics #
######################################################
print "cp $command_log /home/genetics/command_logs/$command_log\n";
system("cp $command_log /home/genetics/command_logs/$command_log");

print "\n\n\n";

###################################################################################
#                  End message to say the program has finished                    #
###################################################################################
print color 'bold white';

print"\n\n\n";
print "###############################################################\n";
print "###############################################################\n";
print "##                 All GWAS analyses complete                ##\n";
print "###############################################################\n";
print "###############################################################\n\n\n";

print "Check the following files:\n\n";
print "(These have been moved into a folder called $prefix_with_suffix"."_results)\n\n";
print "    1.  How to use the output files:        \t$final_readme_file\n";
print "    2.  A list of all the unix commands:    \t$command_log\n";
print "    3.  Program logs of the analyses:       \t$logfile_all\n\n\n";


if($assembly eq "canfam2")
{
	print  "=========================================\n";
	print  "Remember the output files are all canfam2\n";
	print  "=========================================\n\n";

	print  "To use the canfam3 assembly SNP positions you can either:\n\n";
	print  " - convert output assoc files using 'convert_canfam2_assoc_files'\n\n";
	print  " - convert MAP file using 'convert_canfam2_map_files' and re-analyse\n\n";
}


print color 'reset';



######################################################
# Move command log to results folder                 #
######################################################

&move_to_results_folder ("$command_log");

exit;
###################################################################################################################################################


##########################################################
# Subroutine to add single PLINK logfiles to logfile_all #
##########################################################

sub add_to_logfile
{
	# create log file if necessary
	if (! -e $logfile)
	{
		open(LOG,">$logfile") or die "Can't create file $logfile";
		close LOG;
	}

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
			
			print LOG "ERROR found in PLINK logfile. Program will continue to run\n";
			
		}

		# Add line to overall logfile
		print LOGALL "$single_line\n";

	} # end of while loop

	if ($errors_found eq "true"){print COMMAND_LOG "No errors found in PLINK logfile\n\n"}

	close LOG;
}

#############################################
# Subroutine to remove unwanted files       #
#############################################

sub delete_unwanted
{
	my $_suffix 				= "";
	my $_file_to_be_deleted		= "";
	$_suffix 					= $_[0];

	##########################################################################################
	# if the argument is a suffix added to the filename (ie begins with a dot or underscore) #
	##########################################################################################

	if ((index($_suffix,".") == 0) || (index($_suffix,"_") == 0))
	{
		$_file_to_be_deleted = $prefix_with_suffix.$_suffix;
	}

	##########################################################################################
	# if the argument is a complete filename (ie doesn't begin with a dot or underscore)     #
	##########################################################################################

	if ((index($_suffix,".") != 0) && (index($_suffix,"_") != 0))
	{
		$_file_to_be_deleted = $_suffix;
	}

	if (-e "$_file_to_be_deleted")
	{ 
		$command = "rm  $_file_to_be_deleted";
		system("$command");
	}
	else
	{
		print "    >>> The file $_file_to_be_deleted doesn't exist <<<\n";
	}
}

#############################################
# Subroutine to move file to results folder #
#############################################

sub move_to_results_folder
{
	my $_suffix = "";
	my $_file_to_be_moved = "";

	$_suffix = $_[0];
	
	##########################################################################################
	# if the argument is a suffix added to the filename (ie begins with a dot or underscore) #
	##########################################################################################

	if ((index($_suffix,".") == 0) || (index($_suffix,"_") == 0))
	{
		$_file_to_be_moved = $prefix_with_suffix.$_suffix;
	}
	
	##########################################################################################
	# if the argument is a complete filename (ie doesn't begin with a dot or underscore)     #
	##########################################################################################

	if ((index($_suffix,".") != 0) && (index($_suffix,"_") != 0))
	{
		$_file_to_be_moved = $_suffix;
	}	

	if (! -e "$_file_to_be_moved")
	{
		print "\n$_file_to_be_moved could not be found to be moved to $results_folder\n";
		print COMMAND_LOG "\n$_file_to_be_moved could not be found to be moved to $results_folder\n";
	}
	
	if (-e "$_file_to_be_moved")
	{
		$command = "mv  $_file_to_be_moved $results_folder/$_file_to_be_moved";
		print "\n$_file_to_be_moved was moved to $results_folder\n";
		system("$command");
		print COMMAND_LOG "\n  Move to results folder:  $_file_to_be_moved was moved to $results_folder\n";
	}
}





#############################################
# Subroutine to remove unwanted files       #
#############################################

sub delete_unwanted_plink2
{
	my $_suffix = "";
	$_suffix = $_[0];

	if (-e "$prefix_with_suffix".$_suffix)
	{ 
		$command = "rm  $prefix_with_suffix".$_suffix;
		system("$command");
		#print COMMAND_LOG ("$command\n"); # Leave these out as they clutter the command_log file
	}
	else
	{
		print "    >>> The file $prefix_with_suffix".$_suffix." doesn't exist <<<\n";
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

	
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	
	&print_message_command_log("$information","message");

	print COMMAND_LOG "$unix_command\n";

	print "\n";
	print("$unix_command\n");
	system("$unix_command");

	if ($testing_mode eq "on")
	{
		print "\n";
		print "#########################################################\n";
		print "      Testing mode is 'on'. Press RETURN to continue     \n";
		print "#########################################################\n";
		$answer=<STDIN>;
	}

}

######################################
# Subroutine to print screen message #
######################################
sub print_message
{
	my $_message_length 	= "";
	my $_pos_count		= 0;
	my $_char			= "";
	
	my $_message = $_[0];
	my $_style = $_[1];
	
	$_message_length = length($_message);
	
	if ($_style eq ""){$_char = "#"}
	if ($_style eq "input"){$_char = "~"}
	if ($_style eq "message"){$_char = "#"}
	if ($_style eq "warning"){$_char = "!"}
	
	print "\n\n";
	print color ' bold yellow';
	if ($_style eq "warning"){print color ' bold red'}
	if ($_style eq "input"){print color ' bold white'}
	if ($_style eq "message1"){print color ' bold green'}
	
	for ($_pos_count = 1;$_pos_count <=($_message_length + 10);$_pos_count++){print $_char}
	
	print "\n$_char    $_message    $_char\n";
	
	for ($_pos_count = 1;$_pos_count <=($_message_length + 10);$_pos_count++){print $_char}
	
	print "\n\n";
	print color 'reset';

}#

######################################
# Subroutine to print screen message #
######################################

sub print_message_command_log
{
	my $_message 		= "";
	my $_message_length 	= "";
	my $_pos_count		= 0;
	
	$_message = $_[0];
	$_message_length = length($_message);
	
	for ($_pos_count = 1;$_pos_count <=($_message_length + 10);$_pos_count++)
	{
		print COMMAND_LOG "#";
	}
	print COMMAND_LOG "\n#    $_message    #\n";
	for ($_pos_count = 1;$_pos_count <=($_message_length + 10);$_pos_count++)
	{
		print COMMAND_LOG "#";
	}
	print COMMAND_LOG "\n\n";
	

}

###########################################################
# Subroutine to get filename before file type suffix      #
# (i.e. prefix) (e.g. to get "test" from "test.fasta")    #
###########################################################

sub get_prefix
{
	my $_filename = "";

	$_filename = $_[0];
	
	if (index($_filename,".") > 0)
	{
		$_filename = substr($_filename, 0, index($_filename,"."));
	}
	if (index($_filename,".") == -1)
	{
		$_filename = $_filename;
	}
}

sub delay
{
	my $_count = 0;

	for ($_count = 1; $_count <100000; $_count++)
	{

	}
}

##########################################################
# Subroutine to pause until user hits 'return'           #
##########################################################
sub pause
{
	my $_answer = "";
	print "\n Press RETURN to continue\n";
	$_answer=<STDIN>;
}
