#!/usr/bin/perl -w

use Getopt::Std ;
use File::Basename ;

my $q			= "YES";

print "\n";
print "#########################################\n";
print "#     RUNSMARTPCA_QTL perl script v5    #\n";
print "#     (DOG version)                     #\n";
print "#########################################\n";
print "\n";

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
#Input genotype file
$i = $f.".eigenstratgeno"; 

#Input SNP file
$a = $f.".snp";

#Input indiv file
$b = $f.".ind"; 

#Input value of K
$k = 10; 

################
# Output files #
################
#Output file prefix
$o = $f; 

#Output plot file prefix
$p = $f."_plot";
 
#Out file of all eigenvalues
$e = $f.".eval";
 
#Output log file
$l = $f.".smartpca_log"; 

#Optional flags
$m = 5; 
$t = 10; 
$s = 6.0; 


# Parameter file for running smartpca
$parfile = "$o.par"; 

# Name out file with eigenvalues
$evec = "$o.evec";

# Name out file with PCA values
$pca = "$o.pca";

#################################
# Tell user the list of options #
#################################
print ("Input files:\n");
print ("  Genotype file:   \t\t$i\n");
print ("  SNP file:        \t\t$a\n");
print ("  Individual file: \t\t$b\n\n");

print ("Output files:\n");
print ("  Eigenvectors output file: \t$evec\n");
print ("  Eigenvalues output file:  \t$e\n");

print ("\nOther parameters:\n");
print ("  altnormstyle:    \t\tNO\n");
print ("  numoutevec:      \t\t$k\n");
print ("  numoutlieriter:  \t\t$m\n");
print ("  numoutlierevec:  \t\t$t\n");
print ("  outliersigmathresh: \t\t$s\n");

print ("  qtmode:           \t\t$q\n");
print "\n\n";

########################################
#  Create parameter file for smartpca  #
########################################

open(PARFILE,">$parfile") || die("OOPS couldn't open file $parfile for writing");

print PARFILE ("genotypename: $i\n");
print PARFILE ("snpname: $a\n");
print PARFILE ("indivname: $b\n");

print PARFILE ("evecoutname: $evec\n");
print PARFILE ("evaloutname: $e\n");

print PARFILE ("altnormstyle: NO\n");
print PARFILE ("numoutevec: $k\n");
print PARFILE ("numoutlieriter: $m\n");
print PARFILE ("numoutlierevec: $t\n");
print PARFILE ("outliersigmathresh: $s\n");

print PARFILE ("qtmode: $q\n");

close(PARFILE);


##################
#  Run smartpca  #
##################

$command = "/opt/eigensoft/smartpca_dog";          # Use version 2.0 eigensoft for QTLs
$command .= " -p $parfile >$l";

print "Running SMARTPCA_DOG for QTL... (principal components analysis)\n";
print("$command\n");
system("$command");

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

print "Popstring:  $popstring\n\n";

###############
# Run ploteig #
###############
$command = "/opt/eigensoft/ploteig";           # MUST put bin directory in path
$command .= " -i $evec";
$command .= " -c 1:2 ";
$command .= " -p $popstring ";
$command .= " -x ";
$command .= " -o $p.xtxt "; # must end in .xtxt

print "\nRunning PLOTEIG (produces a pdf to show plot of values)...\n";
print("$command\n");
system("$command");


#######################################################################
#  translate .evec file to .pca file expected by eigenstrat program   #
#  Note: .evec file does not contain entries for outliers             #
#       .pca  file does contain entries (set to all 0.0) for outliers #
#######################################################################

$command = "perl /home/genetics/scripts/evec2pca.perl $k $evec $b $pca";
print "\nRunning EVEC2PCA...(converts evec file to pca file)\n";
print("$command\n");
system("$command");


#######################################################################
#  translate .ind file to .pheno file expected by eigenstrat program  #
#  (NOTE: this uses the QTL version which I have written)             #
#######################################################################

$command = "perl /home/genetics/scripts/ind2phenoQTL.pl $b $f.pheno";
print "\nRunning IND2PHENOQTL...(converts indiv file to pheno file for EIGENSTRAT)\n";
print("$command\n");
system("$command");
print "\n";


#######################################################################
#  If labels are Case and Control only, compute correlations between  #
#  Case/Control status and each eigenvector.  Append to logfile.      #
#######################################################################
if(($popstring eq "Case:Control") || ($popstring eq "Control:Case"))
{
  open(LOG,">>$l") || die("OOPS couldn't open file $l for appending");
  print LOG ("\n");
  for($x=0; $x<$k; $x++) # compute correlation for evec $x
  {
    open(EVEC,$evec) || die("OOPS couldn't open file $evec for reading");
    $sum1=0; $sumx=0; $sumxx=0; $sumy=0; $sumyy=0; $sumxy=0;
    $line = <EVEC>; chomp($line); # eigvals line
    while($line = <EVEC>)
    {
      chomp($line);
      my @array = split(/[\t ]+/,$line);
      $this = $array[2+$x];
      $sumy += $this;
      $sumyy += $this*$this;
      $sum1 += 1;
      if($line =~ /Case/) # Case is 1, Control is 0
      {
        $sumx += 1;
        $sumxx += 1;
        $sumxy += $this;
      }
    }
    close(EVEC);
    $meanx = $sumx/$sum1;
    $meany = $sumy/$sum1;
    if($sum1 == 0) { next; }
    $sdevx = sqrt($sumxx/$sum1 - $meanx*$meanx);
    $sdevy = sqrt($sumyy/$sum1 - $meany*$meany);
    if($sdevx * $sdevy == 0) { next; }
    $corr = ($sumxy/$sum1) / ($sdevx*$sdevy);
    $x1 = $x+1;
    printf LOG ("Correlation between eigenvector $x1 (of $k) and Case/Control status is %.03f\n",$corr);
  }
  close(LOG);
}



print "\nCOMPLETED.\n\n";
print "  --> Look at $p.pdf for a plot of the main two principal components\n\n";
print "  --> Run EIGENSTRAT on the output files.\n\n";
