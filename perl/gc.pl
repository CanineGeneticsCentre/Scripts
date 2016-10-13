#!/usr/bin/perl

$P = $ARGV[0];
$out = $ARGV[1];

print "\n\n";
print "##########################\n";
print "#          GC v1         #\n";
print "#                        #\n";
print "# Calculates lambda from #\n";
print "# eigenstrat .chisq file #\n";
print "#                        #\n";
print "##########################\n";
print "\n\n";

########################################
# get data from Eigenstrat .chisq file #
# (into two arrays: chisq1, chisq2)    #
########################################

###########################################
# Expects two arguments, input and output #
###########################################



if ($P eq "")
{
	print "Eigenstrat chisq file:      ";
	$P = <STDIN>;
	chomp $P;

	

	print "\n";
}

if ($out eq "")
{
	print "Output file:      ";
	$out = <STDIN>;
	chomp $out;
	print "\n";
}

$m=0;
open(P,"$P") || die("Can't open file $P");

while($line = <P>) { if($line =~ /Chisq/) { last; } } # header lines
while($line = <P>)
{
  chomp($line);
  if($line =~ /NA/) 
  {
    $chisq1[$m] = -100;
    $chisq2[$m] = -100;
    $m++;  
    next;
  }
  my @array = split(/[\t ]+/,$line);
  $chisq1[$m] = $array[0]; # Chisq
  $chisq2[$m] = $array[1]; # EIGENSTRAT
  $m++; 
  $mvalid++;
}
close(P);
$nSNP = $m;
$nSNPvalid = $mvalid;

print "\nThere are $nSNP SNPs in the file $P\n\n";

print "Computing lambda...\n\n";

###################################
# Compute $lambda1 (Chisq column) #
###################################

$CHISQTHRESH = 0.456;
$step = 0.25;
$oktoreducestep = 0;

for($iter=0; $iter<20; $iter++)
{
  $mm = 0;
  for($m=0; $m<$nSNP; $m++) 
  { 
    if($chisq1[$m] > $CHISQTHRESH) { $mm++; } 
  }

  $frac = $mm/$nSNPvalid; # frac of SNPs exceeding CHISQTHRESH
  if($frac > 0.5) { $CHISQTHRESH += $step; }
  else { $CHISQTHRESH -= $step; $oktoreducestep = 1; }
  if($oktoreducestep) { $step *= 0.5; }
}

$lambda1 = $CHISQTHRESH/0.456; # 0.456 is median if no inflation
if($lambda1 < 1) { $lambda1 = 1; } # not allowed to be less than 1



########################################
# Compute $lambda2 (Eigenstrat column) #
########################################
$CHISQTHRESH = 0.456;
$step = 0.25;
$oktoreducestep = 1;
for($iter=0; $iter<20; $iter++)
{
  $mm = 0;
  for($m=0; $m<$nSNP; $m++) 
  { 
    if($chisq2[$m] > $CHISQTHRESH) { $mm++; } 
  }
  $frac = $mm/$nSNPvalid; # frac of SNPs exceeding CHISQTHRESH
  if($frac > 0.5) { $CHISQTHRESH += $step; }
  else { $CHISQTHRESH -= $step; $oktoreducestep = 1; }
  if($oktoreducestep) { $step *= 0.5; }
}
$lambda2 = $CHISQTHRESH/0.456; # 0.456 is median if no inflation
if($lambda2 < 1) { $lambda2 = 1; } # not allowed to be less than 1


# output
open(OUT,">$out") || die("COF2");
print OUT ("Chisq EIGENSTRAT\n");
printf OUT ("lambda=%.03f lambda=%.03f\n",$lambda1,$lambda2);
for($m=0; $m<$nSNP; $m++)
{
  if($chisq1[$m] < 0) { print OUT ("NA NA\n"); next; }
  printf OUT ("%.04f %.04f\n",$chisq1[$m]/$lambda1,$chisq2[$m]/$lambda2);
}

print  ("Chisq        EIGENSTRAT\n");
printf ("lambda=%.03f lambda=%.03f\n",$lambda1,$lambda2);

print "Input file:    \t$P\n";
print "Output file:   \t$out\n\n";
