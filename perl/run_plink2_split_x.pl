 #!/usr/bin/perl -w

############################################################################
#									                                       #      
#	RUN PLINK2 split-x                      		                               #     
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
use Getopt::Long ;
use File::Basename ;
use Term::ANSIColor;

my $version						= "1";

# Constants
my $edge_of_pseudo_region		= 6717628;
my $end_of_chromosome_x			= 126883977;


# File names
my $ped_file					= "";
my $map_file					= "";
my $ped_file_new				= "";
my $map_file_new				= "";
my $bed_file_new				= "";
my $bim_file_new				= "";
my $fam_file_new				= "";
my $command_log					= "";
my $plink_log_file				= "";

#Strings
my $prefix						= "";
my $prefix_new					= "";
my $assembly					= ""; # canfam2 or canfam3

#Numbers
my $hh_before_split_x 			= 0;
my $hh_after_split_x 			= 0;

##########################################################################################################
# THIS IS THE COMMAND USED                                                                               #
#                                                                                                        #
# plink2 --dog --file BE_SRM_merged_4 --out BE_SRM_merged_4_split --make-bed --split-x 6717628 126883977 #
##########################################################################################################

###############################
# Process flag                #
###############################

GetOptions("file:s"=>\$prefix);

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "              run_plink2_split-x                 \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "    - This program runs the PLINK2 option split-x\n\n";
print "      This takes your exisitIng PED and MAP files and creates a new pair in which\n";
print "      the X and Y chromosomes are handled correctly\n\n";
print "      This makes the following chromosomes:\n\n";
print "      X:   SNPs from non-pseudo-autosomal region of X\n\n";
print "      Y:   SNPs from non-pseudo-autosomal region of Y\n\n";
print "      XY:  SNPs from the pseudo-autosomal region which runs up to $edge_of_pseudo_region\n\n";

print color 'reset';

#######################################
# If file name is on the command line #
#######################################
if ($prefix ne "")
{
	$ped_file = $prefix.".ped";
	$map_file = $prefix.".map";
	if (! -e "$ped_file"){print "\nFile $ped_file doesn't exist. Try again...    \n\n"; exit;}
	if (! -e "$map_file"){print "\nFile $map_file doesn't exist. Try again...    \n\n"; exit;}
}

if ($prefix eq "")
{
	until ((-e $ped_file) && (-e $map_file))
	{
		&print_message("Prefix of PED and MAP files","input");
		$prefix = <STDIN>;
		chomp $prefix;

		$prefix=&get_prefix($prefix);

		$ped_file = $prefix.".ped";
		$map_file = $prefix.".map";
		
		if ($prefix eq "ls")
		{
			print "\n";
			system ("ls *.ped");
			print "\n";
		}
		if ($prefix ne "ls")
		{
			if (! -e "$ped_file"){print "\nFile $ped_file doesn't exist. Try again...    \n\n"; exit;}
			if (! -e "$map_file"){print "\nFile $map_file doesn't exist. Try again...    \n\n"; exit;}
		}
	}
} # if prefix not specified in command line


#######################
# New file names      #
#######################
$prefix_new = $prefix."_split";

$ped_file_new = $prefix_new.".ped";
$map_file_new = $prefix_new.".map";

$bed_file_new = $prefix_new.".bed";
$bim_file_new = $prefix_new.".bim";
$fam_file_new = $prefix_new.".fam";

#######################
# Make log file names #
#######################
$command_log = $prefix."_split_x_command_log.out";
$plink_log_file = $prefix_new.".log";

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";

####################
# Check if canfam2 #
####################

$assembly = &check_if_mapfile_is_canfam2($map_file);

if ($assembly ne "canfam2")
{
	&print_message("This program currently only works for canfam2","warning");

	print "The MAP file $map_file seems to be from canfam3\n\n";

	print "Ask Mike to make a canfam3 version\n\n";
	exit;
}


print COMMAND_LOG "COMMAND LOG for run_plink2_split_x Version $version\n\n";

print COMMAND_LOG "Input PED file:        \t$ped_file\n";
print COMMAND_LOG "Input MAP file:        \t$map_file\n\n";

print COMMAND_LOG "Output PED file:       \t$ped_file_new\n";
print COMMAND_LOG "Output MAP file:       \t$map_file_new\n\n";

&print_message("Splitting X and Y chromosomes for $prefix","message");

&run_unix_command("plink2 --dog --allow-no-sex --nonfounders --make-founders --file $prefix --out $prefix_new --make-bed --split-x $edge_of_pseudo_region $end_of_chromosome_x");	

&add_plink_log_to_command_log(1);

&print_message("Converting binary files back to PED and MAP","message");

&run_unix_command("plink2 --dog --allow-no-sex --nonfounders --make-founders --bfile $prefix_new --recode --out $prefix_new");	

&add_plink_log_to_command_log(2);

&print_message("Splitting finished","message");

#######################
# Delete binary files #
#######################

&delete_file("$bed_file_new");
&delete_file("$bim_file_new");
&delete_file("$fam_file_new");


print  "Input PED file:        \t$ped_file\n";
print  "Input MAP file:        \t$map_file\n\n";

print  "Output PED file:       \t$ped_file_new\n";
print  "Output MAP file:       \t$map_file_new\n\n";


&print_message("Number of heterozygous haploid genotypes present before split-x","message");
&print_both("$hh_before_split_x");

&print_message("Number of heterozygous haploid genotypes present after split-x ","message");
&print_both("$hh_after_split_x\n\n");

&print_both("These are probably bad SNPS.  You can check this in Genome Studio\n\n");

if ((! -e $ped_file_new) || (! -e $ped_file_new))
{
	&print_message("It looks as if the PED and MAP files for $prefix_new were not created","warning");
}

close COMMAND_LOG;

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
	if ($style eq "message"){print COMMAND_LOG "\n\n";}

	print color ' bold yellow';
	if ($style eq "warning"){print color ' bold red'}
	if ($style eq "input"){print color ' bold white'}
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char;if ($style eq "message"){print COMMAND_LOG $char}}
	
	print "\n$char    $message    $char\n";
	if ($style eq "message"){print COMMAND_LOG "\n$char    $message    $char\n";}
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char;if ($style eq "message"){print COMMAND_LOG $char}}
	
	print "\n\n";
	if ($style eq "message"){print COMMAND_LOG "\n\n";}

	print color 'reset';

}



sub run_unix_command
{
	my $unix_command = "";
	$unix_command = $_[0];
		
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------\n";
	print COMMAND_LOG "$unix_command\n";
	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------\n";

	print("$unix_command\n\n");
	system("$unix_command");
}

############################################
#                                          #
# Subroutine to add logfile to logfile_all #
#                                          #
############################################

sub add_plink_log_to_command_log
{
	my $_pass = $_[0];

	my $_single_line 	= "";
	my $_error_pos		= 0;
	my $_het_hap_pos	= 0;
	my $_errors_found	= "";

	# open the logfile #
	open(PLINK_LOG,$plink_log_file) or die "Can't open PLINK logfile $plink_log_file";

	$_errors_found = "false";

	print COMMAND_LOG "\n";

	while ($_single_line = <PLINK_LOG>) 
	{
		chomp $_single_line;

		# Check to see if the the line has an ERROR
		$_error_pos= index($_single_line,'ERROR');

		if ($_error_pos > -1)
		{
			print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
			print "+                                                     +\n";
			print "+        The Logfile contains an ERROR.               +\n";
			print "+                                                     +\n";
			print "+       Please check and fix the problem              +\n";
			print "+                                                     +\n";
			print "+++++++++++++++++++++++++++++++++++++++++++++++++++++++\n";
			
			print PLINK_LOG "ERROR found in PLINK logfile. Program will continue to run\n";

			$_errors_found = "true";
		}

		$_het_hap_pos = index($_single_line,"het. haploid genotypes present");

		if ($_het_hap_pos > 0)
		{
			if ($_pass == 1)
			{
				&print_message("Number of heterozygous haploid genotypes present before split-x","message");
				print "$_single_line\n\n";

				$hh_before_split_x = $_single_line;
			}
			if ($_pass == 2)
			{
				&print_message("Number of heterozygous haploid genotypes present after split-x","message");
				print "$_single_line\n\n";

				$hh_after_split_x = $_single_line;
			}
		}

		# Add line to overall logfile
		print COMMAND_LOG "$_single_line\n";


	} # end of while loop

	close PLINK_LOG;

	if ($_errors_found eq "true"){print COMMAND_LOG "No errors found in PLINK logfile\n\n"}

} # add_plink_log_to_command_log

####################################################
# Subroutine to print to screen and to COMMAND_LOG #
####################################################
sub print_both
{
	my $_message = $_[0];

	print "$_message";
	print COMMAND_LOG "$_message";
}


##############################################
# Subroutine to delete file                  #
##############################################
sub delete_file
{
	my $file_to_be_deleted 	= "";
	my $_command			= "";

	$file_to_be_deleted = $_[0];
	
	if (! -e "$file_to_be_deleted")
	{
		print "\n$file_to_be_deleted could not be found to be deleted\n";
		print COMMAND_LOG "\n$file_to_be_deleted could not be found to be deleted\n";
	}
	
	if (-e "$file_to_be_deleted")
	{
		$_command = "rm  $file_to_be_deleted";
		system("$_command");
		print COMMAND_LOG "\n$file_to_be_deleted was deleted\n";
	}
}

#########################################################
# Check (if dog) whether MAP file is canfam2 or canfam3 #
#########################################################

sub check_if_mapfile_is_canfam2
{
	my $_mapfile 		= $_[0];
	my $_single_line	= "";
	my $_snp_name 		= "";
	my $_position 		= "";
	my $_assembly		= "";
	my $_answer			= "";

	my @_item			= ();

	if (-e $_mapfile)
	{
	 	open(MAP,"$_mapfile") or die "Can't open file $_mapfile";

		 while ($_single_line = <MAP>) 
		{
			chomp $_single_line;
			@_item=split(/\s+/,$_single_line);
			$_snp_name = $_item[1];
			$_position = $_item[3];

			#print "$_position\n";

			if ($_snp_name eq "BICF2G630708027")
			{
				$_assembly = "NOT canfam2 or canfam3";
				if ($_position eq "3713883"){$_assembly = "canfam2"; last}
				if ($_position eq "714375"){$_assembly = "canfam3"; last}
			}
		}
		close MAP;

		$_assembly = $_assembly;
	}
	else
	{
		&print_message("Can't find MAP file $map_file","warning");
		exit;
	}
} # sub check_if_mapfile_is_canfam2