#!/usr/bin/perl

######################################################
#      Progeny2mega2                                 #
#                                                    #
#                                                    #
######################################################

# Mike Boursnell Oct 2007

use strict;
use File::Copy;
use Getopt::Long;

# Define variables
my $prefix       	= "";
my $pedfilename  	= "";
my $pedfilenew		= "";
my $mapfilename  	= "";
my $mapfilenew	 	= "";
my $datafilename  	= "";
my $datafilenew  	= "";
my $chromosome		= "";
my $chromosome_padded	= "";
my $left		= "";
my $help         	= 0;


print"\n";

print"#############################################################\n";
print"#                                                           #\n";
print"#      progeny2mega2                                        #\n";
print"#                                                           #\n";
print"#      This renames the files exported from Progeny Lab     #\n";
print"#      so that they have the default file names             #\n";
print"#      expected by mega2 on the UNIX server.                #\n";
print"#                                                           #\n";
print"#############################################################\n\n";

#############################
# Get command item options  #
#############################
GetOptions("i=s" => \$prefix,
"h" => \$help);
if (($help == 1) or ($prefix eq "")) {

############################
# Get the file name prefix #
############################
print "What is the name of the linkage export from Progeny Lab.\n";
print "(e.g. the part of the file name before the underscore in 'linkage export 2_pedin.2')\n\n";
print "> ";
$prefix = <STDIN>;
chomp $prefix;
print"\n";

############################
# Get the chromosome       #
############################
print "Which chromosome: \n\n";
print "> ";
$chromosome = <STDIN>;
chomp $chromosome;
print"\n";

}

########################################################
# Create chromosome suffix by adding zero if necessary #
########################################################

$left = substr($chromosome,0,1);

if (($chromosome < 10) && ($left ne "0")){
	$chromosome_padded = "0".$chromosome;
}

else {
	$chromosome_padded = $chromosome;
}

################################################
# Generate existing file names from the prefix #
################################################

$pedfilename = $prefix."_pedin.".$chromosome;
$pedfilenew = "pedin.".$chromosome_padded;

$mapfilename = $prefix."_map.".$chromosome;
$mapfilenew = "map.".$chromosome_padded;

$datafilename = $prefix."_datain.".$chromosome;
$datafilenew = "datain.".$chromosome_padded;

###################################
# Check pedfile (at least) exists #
# and only proceed if it does     #
###################################

if (-e $pedfilename){


	#######################################
	# Copy files to default mega2 names   #
	#######################################

	print "Files renamed as follows:\n\n";


	copy($pedfilename,$pedfilenew);
	print "$pedfilename --> $pedfilenew\n";

	copy($mapfilename,$mapfilenew);
	print "$mapfilename --> $mapfilenew\n";

	copy($datafilename,$datafilenew);
	print "$datafilename --> $datafilenew\n";

#############################################################
# Now run fromdos on the files in case they need converting #
# (you need to run the Unix fromdos using 'system')         #
#############################################################

system "fromdos $pedfilenew";
system "fromdos $mapfilenew";
system "fromdos $datafilenew";


} # if pedfile exists

else {
	print "\nERROR: $pedfilename does not exist.\n\n";
}
exit;
