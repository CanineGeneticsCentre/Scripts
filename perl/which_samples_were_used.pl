#!/usr/bin/perl

################################################################
#     which samples were used?                                 #
#                                                              #
#     Reads PLINK .ped and .map files and makes a list of      #
#     which samples were used                                  #
#                                                              #
################################################################


# Mike Boursnell Feb 2009
# Animal Health Trust
# Newmarket
# UK

use strict;
use Getopt::Long;

#####################################
# Define variables                  #
#####################################
my @item         		= ();
my @newline_affected    = ();
my @newline_normal     	= ();
my @SNP_array			= ();
my @id_array			= ();
my @status_array		= ();
my @SNP_name_array		= ();
my @chromosome_array	= ();
my @position_array		= ();

my $infile       		= "";
my $logfile				= "";
my $mapfile				= "";
my $pedfile				= "";
my $outfile_list		= "";
my $outfile_normal		= "";
my $genotype	 		= "";
my $genotype_string		= "";
my $genotype_bases		= "";
my $answer       		= "";
my $single_line			= "";
my $time_string			= "";
my $date_string			= "";
my $sample_name			= "";
my $pedigree			= "";
my $id					= "";
my $father				= "";
my $mother				= "";
my $gender				= "";
my $status				= "";
my $SNP_name			= "";
my $chromosome			= "";
my $position			= "";
my $mystring			= "";
my $chomped_all_string	= "";
my $allele_1			= "";
my $allele_2			= "";

my $linecount    		= 1;
my $colcount     		= 0;
my $showcount    		= 0;
my $itemcount    		= 0;
my $checkcount   		= 0;
my $snp_count			= 0;
my $id_count			= 0;
my $total_lines  		= 0;
my $total_items  		= 0;
my $total_items_on_first_row	= 0;
my $total_SNPs_on_first_row 	= 0;
my $total_SNPs			= 0;
my $year				= 0;
my $month				= 0;
my $day					= 0;
my $hour				= 0;
my $min					= 0;
my $sec					= 0;
my $correct_genotype_count 	= 0;
my $incorrect_genotype_count 	= 0;
my $total_genotypes_read	= 0;
my $help         		= 0;
my $affected_count		= 0;
my $normal_count		= 0;

############################
# process -f flag (if any) #
############################

GetOptions("f=s"=>\$infile);



print"\n";
print"################################################################\n";
print"#     which_samples_were_used v2                               #\n";
print"#                                                              #\n";
print"#     Reads PLINK .ped file and makes a list of                #\n";
print"#     which samples were used.                                 #\n";
print"#                                                              #\n";
print"################################################################\n\n";


#############################
# Get the input file name   #
#############################

if ($infile eq "")
{
	print "Prefix for PLINK .ped file:  ";
	$infile = <STDIN>;
	chomp $infile;
	print"\n";
}

$pedfile = $infile.".ped";
$outfile_list = $infile.".samples_used.out";


############################
# Open files for output    #
############################

open (OUT, ">$outfile_list") || die "Cannot open $outfile_list";


###############################
# Open the ped file for INPUT #
###############################

open (PED, "$pedfile") || die "Cannot open $pedfile";


$linecount=1;

##############################
# Read data in INPUT file    #
##############################

print "Reading PED file...\n\n";

###########################################
# Work down the lines in the file         #
###########################################
while ($single_line = <PED>) {

	# Remove the end of line character #
	chomp $single_line;


	######################################################
        # Split line at tabs or spaces into the array 'item' #
	######################################################
	@item=split(/\s+/,$single_line);


	#####################################################
	# Check how many items there are on the row.        #
	#####################################################
	$itemcount = 0;
	while (length $item[$itemcount]>0)
	{
		$itemcount = $itemcount + 1;
	}
	$total_items=$itemcount;
	$total_SNPs = $total_items - 6;


	##############################
	# Read the first six columns #
	##############################
	$pedigree = $item[0];
	$id = $item[1];
	$father = $item[2];
	$mother = $item[3];
	$gender = $item[4];
	$status = $item[5];


	##################################
	# Assign each id to the id array #
	##################################
	chomp $id;
	$id_array[$linecount] = $id;


	##########################################
	# Assign each status to the status array #
	#########################################
	$status_array[$linecount] = $status;


	######################################################
	# Store the number of items on the first line        #
	# so we can check that all other lines have          #
	# the same number of items                           #
	#                                                    #
	######################################################

	if ($linecount==1) 
	{
		$itemcount = 0;
		while (length $item[$itemcount]>0)
		{
			$itemcount = $itemcount + 1;
		}

		$total_items_on_first_row=$itemcount;
		$total_SNPs_on_first_row = $total_items_on_first_row - 6;

	}


	######################################################
	# For each subsequent line, warn if the number of    #
	# items on the row doesn't match the number on the   #
	# first line.                                        #
	######################################################
	if ($itemcount!=$total_items_on_first_row)
	{
		print "\n";
		print "##############\n";
		print "# WARNING!!! #\n";
		print "##############\n\n";
		print "Number of items on line $linecount is $itemcount\n";
		print "which is not the same as the $total_items_on_first_row items on the first row.\n\n";
		print "Press return to continue\n";
		$answer= <STDIN>;
		print "\n";
	}


	######################################################
	# Move on to next line                               #
	######################################################
	$linecount = $linecount + 1;

} # end of while loop

$total_lines=$linecount -1;

######################
# Show some details  #
######################




#############################################
# Write file of Sample names                #
#############################################

print OUT "Samples used in $pedfile\n\n";

#########################
# First print affecteds #
#########################
for ($id_count = 1;$id_count <=$total_lines;$id_count++)
{
	if ($status_array[$id_count]==2)
	{
		print OUT "$id_array[$id_count]\t$status_array[$id_count]\n";
		$affected_count = $affected_count + 1;
	}

}

#########################
# Then print normals    #
#########################
for ($id_count = 1;$id_count <=$total_lines;$id_count++)
{
	if ($status_array[$id_count]==1)
	{
		print OUT "$id_array[$id_count]\t$status_array[$id_count]\n";
		$normal_count = $normal_count + 1;
	}

}

print OUT "\n";

###########################
#     Close all files     #
###########################
close PED;
close OUT;


print "Total no of samples:  \t$total_lines\n";
print "No of affecteds:      \t$affected_count\n";
print "No of normals:        \t$normal_count\n\n";


print "Output file containing a list of all samples used:\t$outfile_list\n\n";


