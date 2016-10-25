#!/usr/bin/perl

use Getopt::Std ;
use File::Basename ;

print "\n";
print "#######################################\n";
print "#     RUNCONVERTF perl script v3      #\n";
print "#     (HORSE version)                 #\n";
print "#######################################\n";
print "\n";

### process f flag

my @flaglist = ("f");
$x = @ARGV;
for($n=0; $n<$x; $n++)
{
  foreach $flag (@flaglist) 
  {
    if($ARGV[$n] eq "-$flag") { $specified{$flag} = 1; }
  }
}

#This section has been altered since we don't want to use all these flags
foreach $flag ("f")
{
 unless($specified{$flag}) { die("You need to specify the -f flag\n\n"); }
}

getopts('f:',\%opts);

#If you specify a filename with the -f flag then 
#all the other input or output files have the same prefix

$f = $opts{"f"};

################
# Input files  #
################

#Input ped file
$pedfile = $f.".ped";

#Input map file
$mapfile = $f.".map";


################
# Output files #
################
#Output file prefix
$parfile = $f.".par"; 
$genofile = $f.".eigenstratgeno";
$snpfile = $f.".snp";
$indfile = $f.".ind";


##############################
# Open output parameter file #
##############################
open (PAR, ">$parfile") || die "Cannot open: $parfile";

print PAR "genotypename:    $pedfile\n";
print PAR "snpname:         $mapfile\n";
print PAR "indivname:       $pedfile\n";
print PAR "outputformat:    EIGENSTRAT\n";
print PAR "genotypeoutname: $genofile\n";
print PAR "snpoutname:      $snpfile\n";
print PAR "indivoutname:    $indfile\n";
print PAR "familynames:     NO";

close PAR;

#################################
# Tell user the list of files   #
#################################
print ("Filenames:\n");
print ("  Ped file:    \t$pedfile\n");
print ("  Map file:    \t$pedfile\n");


###################
#  Run convertf   #
###################
print "\nRunning CONVERTF_HORSE...\n\n";

$command = "/opt/eigensoft/horse/convertf";
$command .= " -p $parfile";

print("$command\n");
system("$command");

print "\n";
print "#####################################################\n";
print "# Created the following files for use with smartpca #\n";
print "#####################################################\n";
print "\n";

print "Genotype file:  \t$genofile\n";
print "Ind file:       \t$indfile\n";
print "SNP file:       \t$snpfile\n";
print "\n";

