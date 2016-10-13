 #!/usr/bin/perl -w

############################################################################
#									                                       #      
#	RUN VARIANT AFFECT PREDICTOR    		                               #     
#									                                       #
#	THIS PERL SCRIPT WILL RUN Varient Effector Predictor on all VCF files  #
#									                                       #
############################################################################

#############################
# Mike Boursnell July 2012  #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use Getopt::Long;
use File::Basename ;
use Term::ANSIColor;

# Constants
my $version					= "15";
my $vep_path				= "/usr/bin/vep";
my $fork_string				= " --fork 4";

# VEP Database names (species name is used so no update required)
my $canine_database			= "canis_familiaris";
my $equine_database			= "equus_caballus";


#Filenames
my $vcf_file				= "";
my $command_log				= "";
my $output_file				= "";

#Numbers

# Strings
my $prefix					= "";
my $command					= "";
my $answer					= "";
my $ref						= "";
my $ref_seq_name			= "";
my $chromosome				= "";
my $start_time				= "";
my $end_time				= "";
my $run_time				= "";

#Boolean
my $use_defined_region					= "yes";
my $VCF_specified_in_command_line		= "no";
my $species_specified_in_command_line	= "no"; 
my $region_specified_in_command_line	= "no"; 

#Arrays


###############################
# Process flags               #
###############################

GetOptions("file:s"=>\$vcf_file,"ref:s"=>\$ref,"chr:s"=>\$chromosome);


print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "         run_variant_effect_predictor       \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "    - This program runs the variant_effect_predictor (VEP) using the ensembl database\n\n";
print "      This takes a VCF file of SNP or Indel variants and works out what the effect of the mutation will be.\n";
print "      (For example is the mutation non_synonymous, in a splice region etc etc)\n\n";

print color 'reset';

if ($vcf_file ne "")  {print "VCF file specified in command line: \t$vcf_file\n"; $VCF_specified_in_command_line = "yes"}
if ($ref ne "")       {print "Species specified in command line:  \t$ref\n"; $species_specified_in_command_line = "yes"}
if ($chromosome ne ""){print "Region specified in command line:   \t$chromosome\n"; $region_specified_in_command_line = "yes"}


#####################
# Check for symlink #
#####################

if (-d "$ENV{HOME}/.vep")
{
      print "\n.vep directory found\n";
}

if (! -d "$ENV{HOME}/.vep")
{
    print ".vep directory not found\n\n";

    print "Creating symlink...\n\n";

	system ('ln -s /opt/vepcache ~/.vep');
	system ("cd ~/.vep");
	system ("ls");

	exit;
}


if ($vcf_file eq "")
{
	until (-e $vcf_file)
	{
		&print_message("Name of the input VCF file","input");
		$vcf_file = <STDIN>;
		chomp $vcf_file;
		
		if ($vcf_file eq "ls")
		{
			print "\n";
			system ("ls *.vcf");
			print "\n";
		}
		if ($vcf_file ne "ls")
		{
			if (! -e "$vcf_file"){print "\nFile doesn't exist. Try again...    \n";}
		}
		
	}
	
	#########################################
	# Make default name for output VCF file #
	#########################################	
	$prefix = &get_prefix ($vcf_file);
	$output_file=$prefix."_VEP.vcf";
	
	&print_message("Name of the output VCF file","input");
	print "(default = $output_file)\n\n";

	$answer = <STDIN>;
	chomp $answer;
	if ($answer ne "") {$output_file = $answer}


	#########################################
	# Add .vcf suffix if required           #
	#########################################
	if (index($output_file,".vcf") == -1) {$output_file = $output_file.".vcf";}
		
	if (-e $output_file)
	{
		print "\nAn output file called $output_file already exists.\n\n";
		print "Do you want to continue and overwrite this file? (y/n)     ";
		$answer = <STDIN>;
		chomp $answer;
		if (lc($answer) eq "n")
		{
			print "Type in new name for the output file:  ";
			$answer=<STDIN>;
			chomp $answer;
			if ($answer ne ""){$output_file = $answer}else{exit;}
		}
	}
} # if vcf_file not specified in command line



$command_log = $prefix."_VEP_command_log.out";

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "COMMAND LOG for run_variant_effect_predictor version $version\n\n";

if ($VCF_specified_in_command_line eq "yes")    {print COMMAND_LOG "VCF file specified in command line: \t$vcf_file\n\n";}
if ($species_specified_in_command_line eq "yes"){print COMMAND_LOG "Species specified in command line:  \t$ref\n\n";}
if ($region_specified_in_command_line eq "yes") {print COMMAND_LOG "Region specified in command line:   \t$chromosome\n\n";}


if ($vcf_file ne "")  {print "VCF file specified in command line: \t$vcf_file\n"; $VCF_specified_in_command_line = "yes"}
if ($ref ne "")       {print "Species specified in command line:  \t$ref\n"; $species_specified_in_command_line = "yes"}
if ($chromosome ne ""){print "Region specified in command line:   \t$chromosome\n"; $region_specified_in_command_line = "yes"}


if ($ref eq "")
{
	&print_message("Which reference sequence do you want to use?","input");
	print "   <1> Dog   (canFam3)\n";
	print "   <2> Horse (EquCab2)\n";

	$answer = <STDIN>;chomp $answer;

	if ($answer eq ""){print "\n\nYou didn't choose a reference sequence.\n\n";exit;}

	if (substr($answer,0,1) eq "1" ){$ref = $canine_database; $ref_seq_name = "canfam3";}
	if (substr($answer,0,1) eq "2" ){$ref = $equine_database; $ref_seq_name = "equcab2";}

}


if ($chromosome eq "")
{
	&print_message("Would you like to focus on alignments to specific chromosome?","input");

	print "   Enter 1 for YES\n";
	print "   Enter 2 for NO\n\n";

	$answer = <STDIN>;
	chomp $answer;

	if (substr($answer,0,1) eq "1"){$use_defined_region  = "yes"}
	if (substr($answer,0,1) eq "2"){$use_defined_region = "no"}

	if ($use_defined_region eq "yes")
	{
		&print_message("What is the number of the chromosome?","input");
		print "> ";

		$chromosome = <STDIN>;
		chomp $chromosome;
		
	}
} # if chromosome is not specified in command line

if ($chromosome eq "none"){$use_defined_region = "no"} # Used by fastq2vcf to specify no region


################################
# add chr to chromosome number #
################################
if (index($chromosome,"chr") == -1) {$chromosome = "chr"."$chromosome";}




print COMMAND_LOG "Reference sequence:                 \t$ref\n";

if ($use_defined_region eq "yes")
{
	print COMMAND_LOG "Chromosome:                         \t$chromosome\n";
}

print COMMAND_LOG "Input VCF file:                     \t$vcf_file\n";
print COMMAND_LOG "Output VCF file (annotated by VEP): \t$output_file\n\n";

&print_message("Making annotation files with VariantEffectPredictor","message");

$start_time = time();

##############################################################
# If the VCF file exists, run vep (variant_effect_predictor) #
##############################################################
if (-e "$vcf_file")
{
	if ($use_defined_region eq "yes")
	{
		$command = "$vep_path  --species $ref  --format vcf $fork_string --vcf  --chr $chromosome --input $vcf_file --output_file $output_file --verbose --force_overwrite --offline --everything";
		print COMMAND_LOG "$command\n";
	}
	if ($use_defined_region eq "no")
	{
		$command = "$vep_path --species $ref  --format vcf $fork_string --vcf --input $vcf_file --output_file $output_file --verbose --force_overwrite --offline --everything";
		print COMMAND_LOG "$command\n";
	}

	print("\n$command\n");
	system("$command");
}

if (! -e "$vcf_file")
{
	print "\n\nFILE CAN'T BE FOUND: $vcf_file\n\n"; 
}

$end_time = time();
$run_time = $end_time - $start_time;

printf COMMAND_LOG "\nRun time: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];

close COMMAND_LOG;

&print_message("Variant Effect Predictor finished","message");

print  "Input VCF file:                     \t$vcf_file\n";
print  "Output VCF file (annotated by VEP): \t$output_file\n\n";

printf "Run time: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
print "\n\n";


exit;

####################################################################
#                                                                  #
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test" from "test.fasta")                           #
#                                                                  #
####################################################################

sub get_prefix
{
	my $filename = "";

	$filename = $_[0];
	if (index($filename,".") > 0)
	{
		$filename = substr($filename, 0, index($filename,"."));
	}
	if (index($filename,".") == -1)
	{
		$filename = $filename;
	}
}


######################################
# Subroutine to print screen message #
######################################

sub print_message
{
	my $message_length 	= "";
	my $pos_count		= 0;
	my $char			= "";
	
	my $message = $_[0];
	my $style = $_[1];
	
	$message_length = length($message);
	
	if ($style eq ""){$char = "#"}
	if ($style eq "input"){$char = "~"}
	if ($style eq "message"){$char = "#"}
	if ($style eq "warning"){$char = "!"}
	
	print "\n\n";
	print color ' bold yellow';
	if ($style eq "warning"){print color ' bold red'}
	if ($style eq "input"){print color ' bold white'}
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n$char    $message    $char\n";
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n\n";
	print color 'reset';

}

