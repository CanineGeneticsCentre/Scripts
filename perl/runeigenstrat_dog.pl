#!/usr/bin/perl -w
#

#########################################
# Run smart eigenstrat modified by Mike #
#########################################

use Getopt::Std;


print "\n";
print "###########################\n";
print "#  Runeigenstrat_dog v2   #\n";
print "#  (non-QTL mode)         #\n";
print "#                         #\n";
print "#  (runs smarteigenstrat) #\n";
print "###########################\n\n";

getopts('f:',\%opts);

$f = $opts{"f"};

unless ($f)
{
	print "You need to specify the -f flag\n\n";
	exit;
}


######################
# Make up file names #
######################

$genofilename = $f.".eigenstratgeno";
$indfilename  = $f.".ind";
$snpfilename  = $f.".snp";
$pcafilename  = $f.".pca";
$outfilename  = $f.".eigenstrat.out";
$logfilename  = $f.".log";
$parfile = $f.".par";

#########################
# Some other parameters #
#########################
$qtmode = "NO";   # whether this is a QTL or not
$k = 10;          # Number of Principal Components to correct for

#################################
# Tell user the list of files   #
#################################
print ("Input files:\n");
print ("  Genotype file:   \t\t$genofilename \n");
print ("  SNP file:        \t\t$snpfilename  \n");
print ("  Individual file: \t\t$indfilename  \n");
print ("  PCA file:        \t\t$pcafilename  \n\n");

print ("Output files:\n");
print ("  Log file:        \t\t$logfilename  \n");
print ("  Output file:     \t\t$outfilename  \n\n");

print ("Other parameters:\n");
print ("  qtmode:             \t\t$qtmode\n");
print ("  k:             \t\t$k\n\n");


############################################
# Write parameter file for smarteigenstrat #
############################################

open(PAR, ">$parfile") || die("Error:  unable to open $parfile\n");
print PAR "genotypename:  $genofilename\n";
print PAR "snpname:       $snpfilename\n";
print PAR "indivname:     $indfilename\n";
print PAR "pcaname:       $pcafilename\n";
print PAR "outputname:    $outfilename\n";
print PAR "numpc:         $k\n";
print PAR "qtmode:        $qtmode\n";
close(PAR);

#################################################
# run smarteigenstrat using this parameter file #
#################################################
print "Running the smarteigenstrat program...\n";
$cmd = "/opt/eigensoft/dog/smarteigenstrat -p $parfile >$logfilename";
print "$cmd\n";
system($cmd);

print "\n\n";
print "#################################\n";
print "# Eigenstrat analysis completed #\n";
print "#################################\n";

print "\n  --> Use $outfilename in the Excel sheet Eigenstrat_plot to look at the corrected and uncorrectd CHISQ values\n\n";









