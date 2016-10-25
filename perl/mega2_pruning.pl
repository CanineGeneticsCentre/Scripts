#!/usr/bin/perl

################################################################
#     mega2 pruning v1                                         #
#                                                              #
#     Allows user to prune out a list of SNPs                  #
#     from the mega2 input map file                            #
#                                                              #
################################################################


# Mike Boursnell Oct 2007

use strict;
use Getopt::Long;

#####################################
# Define variables                  #
#####################################

#arrays

my @item         		= ();
my @item_list			= ();
my @snp_array			= ();
my @newline      		= ();
my @linearray    		= ();
my @id_array			= ();
my @id_found			= ();
my @remove_line			= ();

#strings

my $infile       		= "";
my $outfile      		= "";
my $listfile			= "";
my $answer       		= "";
my $single_line			= "";
my $individual			= "";
my $id				= "";
my $id_line			= "";
my $pedigree			= "";
my $prefix       		= "";
my $pedfilename  		= "";
my $pedfilenew			= "";
my $mapfilename  		= "";
my $mapfilenew	 		= "";
my $datafilename  		= "";
my $datafilenew  		= "";
my $chromosome			= "";
my $chromosome_padded		= "";
my $left			= "";
my $SNP_name			= "";
my $snp_to_check		= "";

#numbers

my $help         		= 0;
my $linecount    		= 0;
my $listcount			= 0;
my $total_in_list		= 0;
my $itemcount    		= 0;
my $id_count			= 0;
my $id_total			= 0;
my $total_lines  		= 0;
my $total_SNPs			= 0;
my $total_items  		= 0;
my $pedcount			= 1;	# the number of the pedigree
my $totperson			= 0; 	# total number of individual in the pedigree
my $i				= 0;	# simple counter
my $position			= 0;    # map position from mega2 map file
my $included_count		= 0;	# count of lines included in the new output file
my $removed_count		= 0;	# count of markers removed

print"\n";
print"################################################################\n";
print"#     Mega2 pruning                                            #\n";
print"#                                                              #\n";
print"#     Allows user to 'prune' out a list of SNPs                #\n";
print"#     from the map file for mega2                              #\n";
print"#                                                              #\n";
print"################################################################\n\n";


##############################
# Get the map file name      #
##############################
print "Name of map file: ";
$mapfilename = <STDIN>;
chomp $mapfilename;
print"\n";

##############################
# Get the new map file name  #
##############################
print "Name for pruned map file: ";
$mapfilenew = <STDIN>;
chomp $mapfilenew;
print"\n";

##############################
# Check mapfile file exists  #
##############################
if (!-e $mapfilename) {
print "$mapfilename does not exist\n";
exit;
}


####################################################
# Get the file with a list of SNPs to remove       #
####################################################
print "Name of file containing the list of SNPs to be removed (e.g. plink.prune.out): ";
$listfile = <STDIN>;
chomp $listfile;
print"\n";


##############################
# Check SNP list file exists #
##############################
if (!-e $listfile) {
print "SNP list file $listfile does not exist\n";
exit;
}

###################################################
# Read list of SNPs from list file                #
###################################################

open (LIST, "$listfile");

print "List of SNPs to be removed:\n\n";

while ($single_line = <LIST> ) 
	{
	 chomp $single_line;

	#####################################################
        # Split line at tabs into the array 'item'          #
	#####################################################

        @item_list =split(/\t/,$single_line);

	print "$item_list[0]\n";

	$snp_array[$listcount]=$item_list[0];

	$listcount=$listcount + 1;
	}

close(LIST);

$total_in_list = $listcount;

print "\nNo of SNPs in list is $total_in_list\n\n";





###########################
# Open the files for INPUT #
###########################

open (MAP, "$mapfilename");
open (MAPOUT, ">$mapfilenew");


###########################################
# Work down the lines in the 'line' array #
###########################################
while ($single_line = <MAP>) {

	# Remove the end of line character #
	chomp $single_line;

	#####################################################
        # Split line at tabs into the array 'item'          #
	#####################################################
        @item=split(/\t/,$single_line);

	#####################################################
        # Store the whole line in the array linearray       #
	#####################################################
        $linearray[$linecount]=$single_line;

	#####################################################
        # Store various columns in the required arrays      #
	#####################################################

	$chromosome=$item[0];
        $position=$item[1];
        $SNP_name=$item[2];

	#####################################################
	# Noe check down all snp_array in the list to see   #
	# if this one is in the list to be removed          #
	#####################################################

	for ($listcount=0;$listcount<$total_in_list;$listcount++)
	{

		$snp_to_check=$snp_array[$listcount];

		if ($SNP_name eq $snp_to_check)
		{
			print "$SNP_name on line $linecount marked for removal\n";

			$id_found[$listcount]=1;

			#############################################
			# mark this line for removal #
			#############################################

			$remove_line[$linecount]=1;
			$removed_count = $removed_count + 1;

		}

	}


	#####################################################
	# Move on to next line                              #
	#####################################################
	$linecount = $linecount + 1;

} # end of while loop

$total_lines = $linecount;
$total_SNPs = $total_lines - 1;

############################################################
#               Rewrite data to OUTPUT file                #
############################################################


for ($linecount=0;$linecount<$total_lines;$linecount++)
{

 	if ($remove_line[$linecount] == 0){
		print MAPOUT "$linearray[$linecount]\n";
		$included_count = $included_count + 1;	
	}
}

# subtract first line
$included_count = $included_count - 1;


###########################
#     Close both files    #
###########################
close MAP;
close MAPOUT;

print "\n\n";
print "Markers included: \t$included_count\n";
print "Markers removed: \t$removed_count\n";
print "Original total: \t$total_SNPs\n\n";
exit;

