#!/usr/bin/perl

################################################################
#     Split Prefile v1                                         #
#                                                              #
#     This reads in a linkgae format prefile                   #
#     and splits it into separate pedigrees                    #
#                                                              #
################################################################


# Mike Boursnell Aug 2007

use strict;
use Getopt::Long;

#####################################
# Define variables                  #
#####################################

#arrays

my @item         		= ();
my @newline      		= ();
my @linearray    		= ();
my @newlinearray		= ();
my @pedigree			= ();
my @new_pedigree		= ();
my @id				= ();
my @father			= ();
my @mother			= ();
my @sex				= ();
my @affection			= ();


#strings

my $infile       		= "";
my $outfile      		= "";
my $answer       		= "";
my $single_line			= "";


#numbers

my $help         		= 0;
my $linecount    		= 0;
my $colcount     		= 0;
my $showcount    		= 0;
my $itemcount    		= 0;
my $checkcount   		= 0;
my $total_lines  		= 0;
my $total_items  		= 0;
my $pedcount			= 1;	# the number of the pedigree
my $totpedigrees		= 0;	# total number of pedigrees in all
my $no_assigned			= 0;	# number of lines which have been successfully assinged a new pedigree number
my $last_no_assigned		= 0;
my $totperson			= 0; 	# total number of individual in the pedigree
my $not_completed		= 1;	# set to 1 meaning True
my $i				= 0;	# simple counter
my $j				= 0;	# simple counter
my $pass			= 0; 	# may not be used?

print"\n";
print"################################################################\n";
print"#     Split Prefile v1                                         #\n";
print"#                                                              #\n";
print"#     This reads in a linkage format prefile                   #\n";
print"#     and splits it into separate pedigrees                    #\n";
print"#                                                              #\n";
print"################################################################\n\n";


#############################
# Get the file name         #
#############################
print "Name of input file containing unsplit prefile: ";
$infile = <STDIN>;
chomp $infile;
print"\n";


###########################
# Open the file for INPUT #
###########################

open (IN, "$infile");
open (OUT, ">prefile_split.txt");


##############################
# Read data in INPUT file    #
# into the line array        #
##############################

print "Reading input file...\n\n";

###########################################
# Work down the lines in the 'line' array #
###########################################
while ($single_line = <IN>) {

	# Remove the end of line character #
	chomp $single_line;

	#####################################################
        # Split line at spaces into the array 'item'        #
	#####################################################
        @item=split(/\ +/,$single_line);


	#####################################################
        # Store the whole line in the array linearray       #
	#####################################################
        $linearray[$linecount]=$single_line;

	#####################################################
        # Store various columns in the required arrays      #
	#####################################################
        $pedigree[$linecount]=$item[0];
        $id[$linecount]=$item[1];
        $father[$linecount]=$item[2];
        $mother[$linecount]=$item[3];
	$sex[$linecount]=$item[4];
	$affection[$linecount]=$item[5];

	#print "Line: $linecount   Pedigree: $pedigree[$linecount]  ID: $id[$linecount]  Father: $father[$linecount]  Mother: $mother[$linecount]  Sex: $sex[$linecount]  Affection: $affection[$linecount]   ";


	#####################################################
	# Check how many items there are on the row.        #
	#####################################################
	$itemcount = 0;
	while (length $item[$itemcount]>0)
	{
		$itemcount = $itemcount + 1;
	}

	#print "  $itemcount items\n";

	######################################################
	# Store the number of items on the first line        #
	# so we can check that all other lines have          #
	# the same number of items                           #
	######################################################
	if ($linecount==0) 
	{
		$total_items=$itemcount;
		print "Number of items on the first row is $total_items.\n";
	}


	######################################################
	# Warn if the number of items on the row             #
	# doesn't match the number on the first line.        #
	######################################################
	if ($itemcount!=$total_items)
	{
		print "##############\n";
		print "# WARNING!!! #\n";
		print "##############\n\n";
		print "Number of items on line $linecount is $itemcount\n";
		print "which is not the same as the $total_items items on the first row.\n\n";
		$answer= <STDIN>;
		print "\n";
	}

	######################################################
	# Make up the new line minus the pedigree column     #
	######################################################

	for ($i=1;$i<$total_items;$i++){
		$newlinearray[$linecount]=$newlinearray[$linecount]." ".$item[$i];
		}

	######################################################
	# Move on to next line                               #
	######################################################
	$linecount = $linecount + 1;

} # end of while loop

$total_lines=$linecount;
$totperson = $linecount;

print"\nTotal number of individuals in the file: $totperson\n";


############################################################
#      Now split into separate pedigrees by altering       #
#               the pedigree column only                   #
############################################################


# Set all pedigree numbers to zero #

#print "Clearing pedigree numbers...\n";
for ($i=0;$i<$totperson;$i++){
     $pedigree[$i]= 0;
}

#set the first pedigree to 1
$pedigree[0]=1;


print "\n\n";

################################################################
# Now get the loop going, starting at pedigree 1               #
#                                                              #
# When the no_assigned is the same for two consecutive passes  #
# then we go to the next pedigree                              #
################################################################


while ($not_completed == 1){

	$last_no_assigned = $no_assigned;

	# now work your way down the list of individuals
	#print "Before i loop.\n";
	for ($i=0;$i<$totperson;$i++){

		# if pedigree equals pedcount then..
		if ($pedigree[$i] == $pedcount){
		##############################################
            	# Assign the same pedigree to all the other  #
            	# individuals mentioned on this line         #
            	##############################################


			# If father ID is non-zero then check down list
			if ($father[$i]!=0){
				for ($j=0;$j<$totperson;$j++){
			
					if ($id[$j] == $father[$i]){

						# add new pedigree ID
						$pedigree[$j]=$pedcount;
						last;
					}

				}

			}


			# If mother ID is non-zero then check down list
			if ($mother[$i]!=0){
				for ($j=0;$j<$totperson;$j++){
			
					if ($id[$j] == $mother[$i]){
						# add new pedigree ID
						$pedigree[$j]=$pedcount;
						last;
					}

				}

			}


			# Assign same pedigree to all lines mentioning this individual as father or mother
			for ($j=0;$j<$totperson;$j++){
			
				if ($id[$i] eq $father[$j]){
					$pedigree[$j]=$pedcount;
				}
				if ($id[$i] eq $mother[$j]){
					$pedigree[$j]=$pedcount;
				}

			} # j loop

		} # if pedigree = pedcount

	} # i loop



	###############################################
	#  Now count number assigned                  #
	###############################################

	$last_no_assigned = $no_assigned;
	$no_assigned = 0;

	for ($j=0;$j<$totperson;$j++){
		if ($pedigree[$j] > 0) {
			$no_assigned = $no_assigned + 1;
		}
	}

	#print "\nNo assigned: $no_assigned    (Last no assigned: $last_no_assigned)\n\n";

	###############################################
	#  if no_assigned has stayed the same for two #
	#  passes then it is not going to increase    #
	#  any more so go to the next pedcount.       #
	###############################################

	if ($last_no_assigned == $no_assigned) {

		$pedcount = $pedcount + 1;


		#####################################################################
		# Run down the list to find next person with an unassigned pedigree #
		#####################################################################

		$j = 1;
		$not_completed = 0; # i.e. FALSE

		#print "Finding next unassigned pedigree...\n";

		while (($not_completed == 0) && ($j <= $totperson)) {

			if ($pedigree[$j]eq 0){

				$not_completed = 1;
				$pedigree[$j] = $pedcount;
				#print "Next blank pedigree filled at j=$j with pedcount=$pedcount\n";
			}

			$j = $j + 1;

		} # while loop

	} # if

	if ($not_completed == 0){

		$totpedigrees = $pedcount - 1;
		if ($totpedigrees > 1){
			print "Split into $totpedigrees pedigrees completed\n";
		}

	} # if not_completed = FALSE



} # while not_completed loop

############################################################
#               Rewrite data to OUTPUT file                #
############################################################


for ($linecount=0;$linecount<$totperson;$linecount++)
{
	print OUT "$pedigree[$linecount]$newlinearray[$linecount]\n";
}



###########################
#     Close both files    #
###########################
close IN;
close OUT;


print "\nThe split pedigree is in the file:  prefile_split.txt\n\n";
exit;

