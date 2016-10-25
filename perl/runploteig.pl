#!/usr/bin/perl

use Getopt::Std ;
use File::Basename ;

print "\n";
print "#######################################\n";
print "#     RUNPLOTEIG perl script v2.0     #\n";
print "#######################################\n";
print "\n";

### process flags
# -w poplist is compute eigenvectors using populations in poplist, then project
# -y poplistplot is use populations in poplistplot for plot
# -z badsnpfile is use badsnpname: badsnpfile in call to smartpca
my @flaglist = ("f","i","a","b","k","o","p","e","l","m","t","s","w","y","z");
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

getopts('f:i:a:b:k:o:p:e:l:m:t:s:w:y:z:',\%opts);

#If you specify a filename with the -f flag then 
#all the other input or output files have the same prefix

$f = $opts{"f"};

################
# Input files  #
################
#Input genotype file
$i = $f.".eigenstratgeno"; if($specified{"i"}) { $i = $opts{"i"}; }

#Input SNP file
$a = $f.".snp"; if($specified{"a"}) { $a = $opts{"a"}; }

#Input indiv file
$b = $f.".ind"; if($specified{"b"}) { $b = $opts{"b"}; }

#Input value of K
$k = 10; if($specified{"k"}) { $k = $opts{"k"}; }

#############################
# Get the vectors to plot   #
#############################

print "Number of X eigenvector to plot: ";
$xplot = <STDIN>;
chomp $xplot  ;
print"\n";

print "Number of Y eigenvector to plot: ";
$yplot = <STDIN>;
chomp $yplot  ;
print"\n";

################
# Output files #
################
#Output file prefix
$o = $f; if($specified{"o"}) { $o = $opts{"o"}; }

#Output plot file prefix
$p = $f.$xplot.$yplot."_plot";; if($specified{"p"}) { $p = $opts{"p"}; }
 

# Name out file with eigenvalues
$evec = "$o.evec";




############################################
#  make string of populations for ploteig  #
############################################
$popstring = "";
open(EVEC,$evec) || die("OOPS couldn't open file $evec for reading");
while($line = <EVEC>)
{
  chomp($line);
  my @array = split(/[\t ]+/,$line);
  $x = @array;
  if($array[1] =~ /eigvals/) { next; } # eigvals header line
  $pop = $array[$x-1];
  if($popfound{$pop}) { next; }
  $popstring = $popstring . "$pop:";
  $popfound{$pop} = 1;
}
close(EVEC);
chop($popstring); # remove last ":"

print "\nString of populations: $popstring\n\n";

if($specified{"y"}) 
{
  ### make string of populations for ploteig based on -y flag input
  $popstring = "";
  open(Y,$y) || die("COF");
  while($line = <Y>)
  {
    chomp($line);
    $popstring .= "$line:";
  }
  chop($popstring);
}

### cax ploteig


$command = "/opt/eigensoft/ploteig";           # opt is new location
$command .= " -i $evec";
$command .= " -c $xplot:$yplot ";
$command .= " -p $popstring ";
$command .= " -x ";
$command .= " -o $p.xtxt "; # must end in .xtxt
print "\nRunning PLOTEIG (produces a pdf to show plot of values)...\n";
print("$command\n");
system("$command");


print "\nCOMPLETED.\n\n";
print "  --> Look at $p.pdf for a plot of the main two principal components\n\n";
