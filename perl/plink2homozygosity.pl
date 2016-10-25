#!/usr/bin/perl

################################################################
#     plink2homozygosity v15                                   #
#                                                              #
#     Converts PLINK .ped and .map files into the required     #
#     format for the Homozygosity Mapping spreadsheet          #
#                                                              #
#     (Transforms from one row per sample to one row per SNP)  #
################################################################


# Mike Boursnell Jan 2012

use strict;
use Getopt::Long;
use Term::ANSIColor;

my $version					= "15";

#####################################
# Define variables                  #
#####################################
my @item_ped        		= ();
my @item_temp				= ();
my @item					= ();
my @piece					= ();
my @names     				= ();
my @newnames				= ();
my @id_array				= ();
my @snp_array				= ();
my @map_array				= ();
my @item_snp				= ();
my @item_aff				= ();
my @SNP_position_array		= ();
my @sample_name_array		= ();
my @status_array			= ();
my @status_array_sorted		= ();
my @snp_call_rate_array		= ();
my @snp_monomorphic_array	= ();
my @chromosome_array		= ();
my @SNP_name_array			= ();
my @position_array			= ();
my @row1 					= ();
my @row2 					= ();
my @array_2d 				= (\@row1, \@row2);

# Filenames
my $pedfile       			= "";
my $outfile      			= "";
my $logfile					= "";
my $mapfile					= "";
my $mutfile					= "";
my $listfile				= "";
my $affectfile				= "";
my $outfile_affected		= "";
my $outfile_normal			= "";
my $outfile_neither			= "";

#General variables
my $answer					= "";
my $single_line				= "";
my $new_single_line 		= "";
my $list_ID					= "";
my $sanger_ID				= "";
my $date_string				= "";
my $SNP_name				= "";
my $allele_1				= "";
my $allele_2				= "";
my $new_ID					= "";
my $ans						= "";
my $pedigree				= "";
my $pedigree_new			= "";
my $ID						= "";
my $father					= "";
my $mother					= "";
my $gender					= "";
my $status					= "";
my $log						= "";
my $snp						= "";
my $sample_name				= "";
my $status					= "";
my $name 					= "";
my $new_status				= "";
my $genotype				= "";
my $new_status_string		= "";
my $prefix					= "";
my $mystring				= "";
my $first_genotype			= "";
my $check_char				= "";
my $check_ord				= "";
my $separator				= "";
my $command					= "";

# Boolean
my $remove_chromosome_zero	= ""; # Boolean - yes or no
my $no_stops				= "false"; # if it being run as part of wga_analyses then don't stop

my $chromosome			= "";
my $genotype_string		= "";
my $genotype_bases		= "";
my $QTL					= "False";

my $help         		= 0;
my $linecount    		= 0;
my $colcount     		= 0;
my $count			= 0;
my $showcount    		= 0;
my $item_count    		= 0;
my $snp_count			= 0;
my $sample_count		= 0;
my $total_items   		= 0;
my $total_items_row_1		= 0;
my $total_items_temp		= 0;
my $renamed_count		= 0;
my $total_lines  		= 0;
my $total_samples		= 0;
my $total_SNPs			= 0;
my $total_in_list		= 0;
my $lines_counted		= 0;
my $found_on_list		= 0;
my $names_not_used_count	= 0;
my $copy_if_no_renaming		= 0;
my $mut_total			= 0;
my $snp_total_maps		= 0;
my $aff_total			= 0;
my $mut_count			= 0;
my $store_mut_count		= 0;
my $position			= 0;
my $mutation			= 0;
my $SNP_position		= 0;
my $affected			= 0;
my $dot_pos			= 0;
my $snp_missing_data_count	= 0;
my $snp_missing_proportion	= 0;
my $snps_passed_count		= 0;
my $snps_failed_count		= 0;
my $snps_failed_mono_count	= 0;
my $minimum_call_rate		= 0.9;
my $pos				= 0;
my $no_gaps			= 0;
my $no_SNPs_output		= 0;
my $affected_count		= 0;
my $normal_count		= 0;
my $neither_count		= 0;
my $total_affected		= 0;
my $total_normal		= 0;
my $total_neither		= 0;
my $median			= 0;
my $array_size			= 0;

my $missing_data_value		= -9;
my $missing_data_value_out	= "";


print color 'reset';


print color 'bold magenta';

print "\n\n";
print "##############################\n";
print color 'bold white';
print "      plink2homozygosity   \n";
print color 'bold magenta';
print "##############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This PERL script converts PLINK .ped and .map files into the two files\n";
print "    required for the Homozygosity Mapping spreadsheet.\n\n";

print "    (Transforms from one row per sample to one row per SNP\n";
print "    and splits into an affected file and a normal file) .\n\n";
print color 'reset';



############################
# process -f flag (if any) #
############################

GetOptions("file=s"=>\$prefix);

#############################################################################################################
# If the file name is specified in the command line assume the script is being run from within wga_analyses #
#############################################################################################################
if ($prefix ne ""){$no_stops = "true"}else{$no_stops = "false"}


##########################################
# Get the file names and open the files  #
##########################################

if ($prefix eq "")
{
	until ((-e $pedfile) && (-e $mapfile))
	{
		print_message("Prefix for PLINK .map and .ped files","input");
		print ">  ";
		$prefix = <STDIN>;
		chomp $prefix;

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








$pedfile = $prefix.".ped";
$mapfile = $prefix.".map";
$logfile = $prefix."_p2h.log";

$outfile_affected = $prefix."_affected.txt";
$outfile_normal = $prefix."_normal.txt";
$outfile_neither = $prefix."_nophenotype.txt";

if ($no_stops eq "false")
{
	&print_message("What do you want to do if the chromosome is labelled as '0'","input");

	print "   <1>  Add it to the output file\n\n";
	print "   <2>  Don't add it to the output file\n\n";

	$answer=<STDIN>;
	chomp $answer;

	if (substr($answer,0,1) eq "1" ){$remove_chromosome_zero = "no"}
	if (substr($answer,0,1) eq "2" ){$remove_chromosome_zero = "yes"}
}
else
{$remove_chromosome_zero = "yes"}

open (PED, "$pedfile") || die "Cannot open input file: $pedfile";
open (MAP, "$mapfile") || die "Cannot open input file: $mapfile";




###################################################
# Read in all SNPs in map file                    #
# and put them in snp_array                       #
# (This is so they can be used as column headers) #
###################################################
print "Reading MAP file...\n\n";
print "Showing first 5 SNPs\n\n";

$count =0;
@map_array= <MAP>;
foreach $single_line (@map_array)
{
	&chomp_all ($single_line);
	@item=split(/\t/,$single_line);

	

	$chromosome = $item[0];
	$SNP_name = $item[1];
	$position = $item[3];

	if ($remove_chromosome_zero eq "no")
	{
		$linecount = $linecount + 1;

		$chromosome_array[$linecount]=$chromosome;
		$SNP_name_array[$linecount]=$SNP_name;
		$position_array[$linecount]=$position;
	}

	if ($remove_chromosome_zero eq "yes")
	{
		if ($chromosome ne "0")

		{
			$linecount = $linecount + 1;
			$chromosome_array[$linecount]=$chromosome;
			$SNP_name_array[$linecount]=$SNP_name;
			$position_array[$linecount]=$position;
		}
	}


	if ($linecount <= 5)
	{
		print "$SNP_name\n";
	}

}
close MAP;

$snp_total_maps = scalar (@map_array);


print "...\n\n";
print "$snp_total_maps SNPs in the MAP file\n\n";


####################
# Re-open PED file #
####################

open (PED, "$pedfile") || die "Cannot open input file: $pedfile";


###########################################
# Work down the lines in the PED file     #
# to store it all in a 2D array           #
###########################################
$linecount=0;
while ($single_line = <PED>) {

	# Remove the end of line character #
	chomp $single_line;
	
#&chomp_all ($single_line);

	$linecount = $linecount + 1;


	#########################################################
        # Split line at TABS OR SPACES into the array 'item'    #
	#########################################################
        @item_ped=split(/\s+/,$single_line);


	######################################################
	# Store the number of items on the first line        #
	# so we can check that all other lines have          #
	# the same number of items                           #
	######################################################
	if ($linecount==1) 
	{
		$total_items_row_1=@item_ped;
		$total_SNPs = ($total_items_row_1 -6)/2;

		print "Number of items on the first row of the .ped file is $total_items_row_1\tTotal SNPs: $total_SNPs\n\n";
	}


	$pedigree  = $item_ped[0];
	$ID = $item_ped[1];
	$father  = $item_ped[2];
	$mother  = $item_ped[3];
	$gender  = $item_ped[4];
	$status  = $item_ped[5];


	#################################
	# Store ID and status in arrays #
	#################################

	$id_array[$linecount] = $ID;
	$status_array[$linecount] = $status;


	#####################################################
	# Check how many items there are on the row.        #
	#####################################################
	$total_items = @item_ped;


	##################################
	# Store alleles in the 2D array   #
	###################################

	for ($snp_count = 1; $snp_count <= $total_SNPs; $snp_count++)
	{
		$allele_1 = $item_ped[($snp_count * 2) + 4];
		$allele_2 = $item_ped[($snp_count * 2) + 5];

		chomp $allele_1;
		chomp $allele_2;

		$genotype="$allele_1\t$allele_2";

		if ($allele_1 eq '1'){$allele_1 = "A"}
		if ($allele_1 eq '2'){$allele_1 = "C"}
		if ($allele_1 eq '3'){$allele_1 = "G"}
		if ($allele_1 eq '4'){$allele_1 = "T"}
		if ($allele_1 eq '0'){$allele_1 = ""}

		if ($allele_2 eq '1'){$allele_2 = "A"}
		if ($allele_2 eq '2'){$allele_2 = "C"}
		if ($allele_2 eq '3'){$allele_2 = "G"}
		if ($allele_2 eq '4'){$allele_2 = "T"}
		if ($allele_2 eq '0'){$allele_2 = ""}

		$genotype_bases="$allele_1\t$allele_2";


		if ($genotype_string eq "0 0")
		{
			$genotype_bases = "$missing_data_value_out\t$missing_data_value_out";
		}

		$array_2d[$linecount][$snp_count] = $genotype_bases;

	}


	print "Sample: $linecount\n";

	######################################################
	# Warn if the number of items on the row             #
	# doesn't match the number on the first line.        #
	######################################################
	if ($total_items!=$total_items_row_1)
	{
		print "##############\n";
		print "# WARNING!!! #\n";
		print "##############\n\n";
		print "Number of items on line $linecount is $item_count\n";
		print "which is not the same as the $total_items_row_1 items on the first row.\n\n";
		$ans= <STDIN>;
		print "\n";
	}


} # end of while loop

$total_samples = $linecount;

print "\n";

&print_message("Number of samples and numbers of SNPs","message");


print "Number of samples in the PED file:  \t$total_samples\n";
print "Number of SNPs in the PED file:     \t$total_SNPs \n\n";


##############################################################
# Check in the status array to see if the phenotype is a QTL #
##############################################################
print "Checking whether phenotype is a QTL:       ";

for ($sample_count = 1;$sample_count <= $total_samples;$sample_count++)
{
	$status=$status_array[$sample_count];

	if ($status ne "0" && $status ne "1" && $status ne "2" and $status ne $missing_data_value)
	{
		$QTL = "True";
	}
}
if ($QTL eq "True")
{
	print "treating trait as quantitative\n";
}

if ($QTL eq "False")
{
	print "treating trait as binary\n";
}

#######################################################
# Check the number of affected and unaffected samples #
#######################################################
for ($sample_count = 1;$sample_count <= $total_samples;$sample_count++)
{
	if ($status_array[$sample_count]==2)
	{
		$affected_count = $affected_count + 1;

	}
	if ($status_array[$sample_count]==1)
	{
		$normal_count = $normal_count + 1;
	}
	if (($status_array[$sample_count] != 1)  && (($status_array[$sample_count] != 2)))
	{
		$neither_count = $neither_count + 1;
	}

}
$total_affected = $affected_count;
$total_normal = $normal_count;
$total_neither = $neither_count;

&print_message("Counts of Affected, Normal or Neither","message");

print "\nTotal affected:      \t$total_affected\n";
print "Total normal:         \t$total_normal\n\n";

if($total_neither > 0)
{
	print "Total neither:        \t$total_neither\n\n";
}
print "Total samples:        \t$total_samples\n\n";

if (($total_affected == 0) && ($total_normal > 0))
{
	print "##########################################\n";
	print "#    No affected samples in ped file!    #\n";
	print "##########################################\n\n";
	
	print "There are no affected phenotypes in the .ped file.\n\n";
	print "No affected.out file will be produced\n\n";
	print "Press 'return' to continue  >";
	
	$ans = <STDIN>;
	print "\n";
}

if ($total_affected > 0 && $total_normal == 0)
{
	print "############################################\n";
	print "#    No unaffected samples in ped file!    #\n";
	print "############################################\n\n";
	
	print "There are no unaffected phenotypes in the .ped file.\n\n";
	print "No unaffected.out file will be produced\n\n";
	print "Press 'return' to continue  >";

	$ans = <STDIN>;
	print "\n";
}

if ($total_affected == 0 && $total_normal == 0)
{
	print "####################################\n";
	print "#    No phenotypes in ped file!    #\n";
	print "####################################\n\n";
	
	print "There are no affected or unaffected phenotypes in the .ped file.\n\n";
	print "It is possible that all phenotypes for these samples are held in an external phenotype file.\n\n";
	print "You can use the PLINK function --recode with --pheno to make a new ped file which includes phenotypes.\n\n";
	print "An output file will be created with all the un-phenotyped samples ($outfile_neither)\n\n";
	print "Press 'return' to continue  >";
	
	$ans = <STDIN>;
	print "\n";
}

################################################
# Print out the array-2d to the output file(s) #
# (new format is one row per SNP)              #
################################################




###################
# Trait is binary #
###################
if ($QTL eq "False")
{

	###########################
	# Open any files required #
	###########################

	if($total_affected > 0)
	{
		open (OUTA, ">$outfile_affected")|| die "Cannot create output file: $!";
	}

	if($total_normal > 0)
	{
		open (OUTN, ">$outfile_normal")|| die "Cannot create output file: $!";
	}

	if($total_neither > 0)
	{
		open (OUTX, ">$outfile_neither")|| die "Cannot create output file: $!";
	}

	if($total_affected > 0)
	{
		print "Writing output file for affected samples...\n";
	}
	if($total_normal> 0)
	{
		print "Writing output file for unaffected samples...\n";
	}
	if($total_neither > 0)
	{
		print "Writing output file for samples with no phenotype...\n";
	}

#############################################
# Write first line of Sample names          #
#                                           #
# if status=2 then write to affected file   #
# if status=1 then write to normal file     #
#############################################

# Write headers for the first three columns for both files #
if($total_affected > 0){print OUTA "CHR\tSNP_ID\tPOSITION\t";}
if($total_normal > 0){print OUTN "CHR\tSNP_ID\tPOSITION\t";}
if($total_neither > 0){print OUTX "CHR\tSNP_ID\tPOSITION\t";}

# Write out the headers for the sample columns
for ($sample_count = 1;$sample_count <= $total_samples;$sample_count++)
{
	if ($status_array[$sample_count]==2)
	{
		print OUTA "$id_array[$sample_count]-A1\t$id_array[$sample_count]-A2\t";
	}
	if ($status_array[$sample_count]==1)
	{
		print OUTN "$id_array[$sample_count]-A1\t$id_array[$sample_count]-A2\t";
	}
	if (($status_array[$sample_count] != 1)  && (($status_array[$sample_count] != 2)))
	{
		print OUTX "$id_array[$sample_count]-A1\t$id_array[$sample_count]-A2\t";
	}
}

if($total_affected > 0){print OUTA "\n";}
if($total_normal > 0){print OUTN "\n";}
if($total_neither > 0){print OUTX "\n";}

$affected_count = 0;
$normal_count = 0;
$neither_count = 0;

print "\n";

for ($snp_count = 1;$snp_count <= $total_SNPs; $snp_count ++)
{
	if($total_affected > 0)
	{
		print OUTA "$chromosome_array[$snp_count]\t";
		print OUTA "$SNP_name_array[$snp_count]\t";
		print OUTA "$position_array[$snp_count]\t";
	}

	if($total_normal > 0)
	{
		print OUTN "$chromosome_array[$snp_count]\t";
		print OUTN "$SNP_name_array[$snp_count]\t";
		print OUTN "$position_array[$snp_count]\t";
	}
	
	if($total_neither > 0)
	{
		print OUTX "$chromosome_array[$snp_count]\t";
		print OUTX "$SNP_name_array[$snp_count]\t";
		print OUTX "$position_array[$snp_count]\t";
	}

	$affected_count = 0;
	$normal_count = 0;
	$neither_count = 0;
	
	if ($snp_count % 10000 ==0)
	{
		print "SNP: $snp_count\n";
	}

	for ($sample_count = 1; $sample_count <=$total_samples;$sample_count ++)
	{

		if ($status_array[$sample_count]==2)
		{
			print OUTA "$array_2d[$sample_count][$snp_count]";
			$affected_count = $affected_count + 1;
			if ($affected_count < $total_affected)
			{
				print OUTA "\t";
			}

		}

		if ($status_array[$sample_count]==1)
		{
			print OUTN "$array_2d[$sample_count][$snp_count]";
			$normal_count = $normal_count + 1;

			if ($normal_count < $total_normal)
			{
				print OUTN "\t";
			}

		}
		
		if (($status_array[$sample_count] != 1)  && (($status_array[$sample_count] != 2)))
		{
			print OUTX "$array_2d[$sample_count][$snp_count]";
			$neither_count = $neither_count + 1;

			if ($neither_count < $total_neither)
			{
				print OUTX "\t";
			}

		}

	}

if($total_affected > 0){print OUTA "\n";}
if($total_normal > 0){print OUTN "\n";}
if($total_neither > 0){print OUTX "\n";}

}

} #If QTL is False



#########################
# Trait is quantitative #
#########################
if ($QTL eq "True")
{

##########################################
# Open any files required                #
# (for QTL only affected and unaffected) #
##########################################

open (OUTA, ">$outfile_affected")|| die "Cannot create output file: $!";
open (OUTN, ">$outfile_normal")|| die "Cannot create output file: $!";

print "\nWriting two output files for quantitative data...\n\n";

###################################################
# Calculate median value of the trait to act as a #
# split point between the two files               #
###################################################

# Make a sorted COPY of the status array

@status_array_sorted = sort(@status_array);


if ($total_samples % 2 == 0) 
{ 
	$median =  ($status_array_sorted[$total_samples/2] + $status_array_sorted[$total_samples/2 + 1]) / 2;
} 
else 
{ 
	$median =  $status_array_sorted[int($total_samples/2)]; 
} 

print "Median value is $median\n\n";
print "Samples with trait value greater than $median go into $outfile_affected\n";
print "Samples with trait value less than $median go into $outfile_normal\n\n";


####################################################
# Write first line of Sample names                 #
#                                                  #
# if status > median then write to affected file   #
# if status < median then write to normal file     #
####################################################

# Write headers for the first three columns for both files #
print OUTA "CHR\tSNP_ID\tPOSITION\t";
print OUTN "CHR\tSNP_ID\tPOSITION\t";

# Write out the headers for the sample columns
for ($sample_count = 1;$sample_count <= $total_samples;$sample_count++)
{
	if ($status_array[$sample_count] >= $median)
	{
		print OUTA "$id_array[$sample_count]-A1\t$id_array[$sample_count]-A2\t";
		$affected_count = $affected_count + 1;

	}
	if ($status_array[$sample_count] < $median)
	{
		print OUTN "$id_array[$sample_count]-A1\t$id_array[$sample_count]-A2\t";
		$normal_count = $normal_count + 1;
	}

}

print OUTA "\n";
print OUTN "\n";

$total_affected = $affected_count;
$total_normal = $normal_count;

print "Total affected (trait >= $median):      \t$total_affected\n";
print "Total normal (trait < $median):         \t$total_normal\n\n";

$affected_count = 0;
$normal_count = 0;

for ($snp_count = 1;$snp_count <= $total_SNPs; $snp_count ++)
{
	print OUTA "$chromosome_array[$snp_count]\t";
	print OUTN "$chromosome_array[$snp_count]\t";

	print OUTA "$SNP_name_array[$snp_count]\t";
	print OUTN "$SNP_name_array[$snp_count]\t";

	print OUTA "$position_array[$snp_count]\t";
	print OUTN "$position_array[$snp_count]\t";

	$affected_count = 0;
	$normal_count = 0;

	for ($sample_count = 1; $sample_count <=$total_samples;$sample_count ++)
	{

		if ($status_array[$sample_count] >= $median)
		{
			print OUTA "$array_2d[$sample_count][$snp_count]";
			$affected_count = $affected_count + 1;
			if ($affected_count < $total_affected)
			{
				print OUTA "\t";
			}

		}

		if ($status_array[$sample_count] < $median)
		{
			print OUTN "$array_2d[$sample_count][$snp_count]";
			$normal_count = $normal_count + 1;

			if ($normal_count < $total_normal)
			{
				print OUTN "\t";
			}

		}

	}

	print OUTA "\n";
	print OUTN "\n";

}

} #If QTL is True


###########################
#     Close all files     #
###########################
close PED;
close OUT;
close OUTA;
close OUTN;
close OUTX;
close OUTQ;

print "\nConverting files to DOS format...\n";

######################################
# Use todos to convert to DOS format #
######################################

if ($total_affected > 0)
{
	$command = "todos $outfile_affected";
	print("$command\n");
	system("$command");
}

if ($total_normal > 0)
{
	$command = "todos $outfile_normal";
	print("$command\n");
	system("$command");
}

if ($total_neither > 0)
{
	$command = "todos $outfile_neither";
	print("$command\n");
	system("$command");
}

print "\n";

&print_message("Output files","message");

if ($total_affected > 0)
{
	print "   Output file (affecteds):     \t$outfile_affected\n";
}

if ($total_normal > 0)
{
	print "   Output file (normals):       \t$outfile_normal\n";
}

if ($total_neither > 0)
{
	print "   Output file (no phenotype):  \t$outfile_neither\n";
}

    print "   Log file:                    \t$logfile\n\n";

print "   (Output files are in DOS format)\n\n";


exit;


##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

foreach (@_) {s/\n//g}  
foreach (@_) {s/\r//g}  

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
