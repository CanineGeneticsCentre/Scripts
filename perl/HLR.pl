use strict;

use POSIX "log10";
use Math::BigInt;
use Math::BigFloat;
use Math::Round ('round');
use Getopt::Long;

my $out_count; # counter to help visualise the output process

use Time::HiRes ('time');

my ($hapassocf,$famf,$test);
GetOptions("f=s",\$hapassocf,"t",\$test,"fam=s",\$famf);

my $mb = 1000*1000;

my @ohs = ("chr","start","end","snpcnt","size","f_a","f_u","chisq","p","common","ucsc","haplotype","snps");

print "####################################\n";
print "#  Haplotype Likelihood Ratio v2   #\n";
print "#                                  #\n";
print "#  Elinor Karlsson                 #\n";
print "####################################\n\n";

#####################
# Check files exist #
#####################

unless (-e $hapassocf) 
{
	print "Can't find hap assoc file $hapassocf\n";
} 

unless (-e $famf) 
{
	print "Can't find fam file $famf\n";
} 




#################################################################################
# In maximum likelihood search, check alphas between 0 and 1 in 0.05 increments #
#################################################################################

my @alphas;
my $i = 0;
for (my $i = 0; $i <= 1; $i += 0.05){
  $i = round($i*100)/100;
  push @alphas, $i;
}


#############################################################
# Get numbers of cases and controls from the plink fam file #
#############################################################

my %cnts;
open IN,$famf or die;
while (<IN>){
  chomp;
  my $st = (split /\s+/)[5];
  $cnts{all}++;
  if ($st == 1){
    $cnts{u}++;
  }
  elsif ($st == 2){
    $cnts{a}++;
  }
}
close IN;

print "\nParsed dog counts from $famf ($cnts{a} cases and $cnts{u} controls)\n\n";


################################################################
# This array stores factorials in advance to speed up run time #
################################################################
my @factorials;

foreach my $i (0 .. $cnts{a}*2+10){
  my $fact = Math::BigInt->bfac($i);
  push @factorials, $fact;
}
print "Got factorials through ".($cnts{a}*2)."\n";


# my $time = time();
my $out = $hapassocf.".LOD.txt";
if ($test){
  $out =~ s/\.txt$/.test.txt/;
}

### this hash stores information on the frequencies of each haplotype (allele) in each haplotype region in cases and controls
my %freqs;

### this hash stores info on the size and position of each haplotype region tested
my %haplotype_info;

##############################################################################################
# parse plink SNP assoc file for allele frequencies in cases and controls (1 SNP haplotypes) #
# (checks to see if the file name ends in '.assoc')                                          #
##############################################################################################

if ($hapassocf =~ m/\.assoc$/){
  my $lncnt = 0;
  my @hs;
  print "Parsing $hapassocf (plink SNP assoc file = 1 SNP haplotypes)\n";
  open IN, $hapassocf or die "Dying in ".__FILE__." at line ".__LINE__.":can't open $hapassocf\n";

  while(<IN>){
    chomp;
    print "Line $lncnt\n";
    $lncnt++;
    s/^\s+//;
    if (@hs == 0){
      $_ = lc $_;
      @hs = split /\s+/;
    }
    else {
      my @d = split /\s+/;
      my %h;
      foreach my $i (0 .. @d-1){
	$h{@hs[$i]} = @d[$i];
#	print "($i) @hs[$i] = @d[$i]\n";
      }

      $freqs{$h{snp}}{$h{a1}}{a} = $h{f_a};
      $freqs{$h{snp}}{$h{a1}}{u} = $h{f_u};
      $freqs{$h{snp}}{$h{a2}}{a} = 1-$h{f_a};
      $freqs{$h{snp}}{$h{a2}}{u} = 1-$h{f_u};

      my ($chr,$pos) = split /\./, $h{snp} or die;

      $chr =~ s/^chr//;
      $haplotype_info{$h{snp}}{chr} = $chr;
      $haplotype_info{$h{snp}}{start} = $pos;
      $haplotype_info{$h{snp}}{end} = $pos;
    }

    if ($lncnt > 10 && $test){
      last;
    }
  }
  close IN;
  print "Parsed $hapassocf (1 SNP haplotypes)\n\n";
}


################################################################################
# parse plink haplotype association file for frequencies in cases and controls #
# (if the filename doesn't end in '.assoc' e.g. ends in '.hap')                #
################################################################################

else {
  my $lncnt = 0;
  my @hs;
  print "Parsing $hapassocf\n";
  open IN, $hapassocf or die "Dying in ".__FILE__." at line ".__LINE__.":can't open $hapassocf\n";
  while(<IN>){
    chomp;
    s/^\s+//;
    if (@hs == 0){
      $_ = lc $_;
      @hs = split /\s+/;
    }
    else {
      my @d = split /\s+/;
      my %h;
      my $temp;
      foreach my $i (0 .. @d-1){
	$h{@hs[$i]} = @d[$i];
	$temp .= "($i) @hs[$i] = @d[$i]\n";
      }
#      die "Dying in ".__FILE__." at line ".__LINE__.":\n$temp\n";
      if (m/OMNIBUS/i){ 
	$lncnt++;
        print "Line $lncnt\n";
      }
      else {

	$freqs{$h{snps}}{$h{haplotype}}{a} = $h{f_a};
	$freqs{$h{snps}}{$h{haplotype}}{u} = $h{f_u};
	my @snps = split /\|/, $h{snps};
	if (! exists $haplotype_info{$h{snps}}){
	  foreach my $snp (@snps){
	    if (! $haplotype_info{$h{snps}}{chr}){
	      $haplotype_info{$h{snps}}{chr} = $snp;
	      $haplotype_info{$h{snps}}{chr} =~ s/\.[0-9]+$//;
	      $haplotype_info{$h{snps}}{chr} =~ s/^chr//i;
	    }
	    my $p = $snp;
	    $p =~ s/^.+\.([0-9]+)$/$1/;
	    if (!$haplotype_info{$h{snps}}{start} || $haplotype_info{$h{snps}}{start} > $p){
	      $haplotype_info{$h{snps}}{start} = $p;
	    }
	    if (!$haplotype_info{$h{snps}}{end} || $haplotype_info{$h{snps}}{end} < $p){
	      $haplotype_info{$h{snps}}{end} = $p;
	    }
	  }
	}
      }
    }
    if ($test && $lncnt > 10){
      last;
    }
  }
  close IN;
  print "Parsed $hapassocf\n\n";
}


############################################################################
# These hashes just make the output prettier by sorting by genome position #
############################################################################

print "\nSorting by genome position\n";

my %order_by_chrpos;

foreach my $name (keys %haplotype_info){
  my $chr = $haplotype_info{$name}{chr};
  $chr =~ s/X/39/;
  my $pos = $haplotype_info{$name}{start};
  push @{$order_by_chrpos{$chr}{$pos}}, $name;
}
my %order_by_idx;
my $idx = 1;

foreach my $chr (sort {$a <=> $b} keys %order_by_chrpos){
  foreach my $pos (sort {$a <=> $b} keys %{$order_by_chrpos{$chr}}){
    foreach my $n (sort @{$order_by_chrpos{$chr}{$pos}}){
      $order_by_idx{$n} = $idx;
      $idx++;
    }
  }
}

############################
# Write to the output file #
############################

$out_count = 0;

print "Writing LOD scores to $out\n\n";
open OUT, ">$out" or die;
print OUT "Type\tName\tHap\tAlpha\tLOD\tHaps\tExp\tObs\n";
foreach my $str (sort {$order_by_idx{$a} <=> $order_by_idx{$b}} keys %freqs){
  my @obs;
  my @obsfs;
  my @exp;
  my @haps;

  my $obscnt = 0;
  my $expsum = 0;

  

  
  foreach my $hap (sort {$freqs{$str}{$b}{a} <=> $freqs{$str}{$a}{a}} keys %{$freqs{$str}}){
    ### observed frequencies are frequecies in cases
    push @obs, round(2*$cnts{a}*$freqs{$str}{$hap}{a});
    $obscnt += round(2*$cnts{a}*$freqs{$str}{$hap}{a});
    push @obsfs, $freqs{$str}{$hap}{a};

    ### expected frequencies are frequecies in controls
    push @exp, $freqs{$str}{$hap}{u};
    $expsum += $freqs{$str}{$hap}{u};
    push @haps, $hap;
  }

  my $expstr = join ",", @exp;
  my $obsstr = join ",",@obsfs;

  ##############################################################################################
  # Correct for frequencies of zero using pseudocounts (probabilities of zero are not allowed) #
  ##############################################################################################

  @exp = correct(\@exp,$cnts{u}+$cnts{a});

  #############################################################################
  # Calculate the denominator of the likelihood ratio - 
  # the probability of the observed frequences given the expected frequencies #
  #############################################################################

  my $denprob = getMultProb(\@obs,\@exp);
  
  my $bestlod = log10(1);
  my $besthap = undef;
  my $bestalpha = 1;
  my $bstr;
  
  if ($denprob == 0){   ### if denominator is 0, indicates something odd happened
    die "Dying in ".__FILE__." at line ".__LINE__.": Calculated denominator probabilty of zero\n\tObserved = (".(join ",", @obs).")\n\tExpected = (".(join ",", @exp).")\n";

  }
  

  ##############################################################################
  # Calculate the nominator using maximum likelihood values for the associated #
  # haplotype and alpha (test each pair of haplotype and alpha)                #
  ##############################################################################

  else {
    foreach my $hapi (0 .. @exp-1){ 
      foreach my $alpha (@alphas){
	my @new; ### New is the expected frequencies based on haplotype and alpha
	my $sum = 0;
	my $old = 0;
	foreach my $hapj (0 .. @exp-1){
	  my $nfreq = @exp[$hapj] * (1-$alpha);
	  if ($hapi == $hapj){
	    $nfreq += 1*$alpha;
	  }
	  push @new, $nfreq;
	  $sum += $nfreq;
	  $old += @exp[$hapj];
	}
	@new = correct(\@new,$cnts{u}+$cnts{a});
	my $numprob = getMultProb(\@obs,\@new);
	
	if ($numprob == 0){   ### if numerator is 0, indicates something odd happened
	 die "Dying in ".__FILE__." at line ".__LINE__.": Calculated numerator probabilty of zero\n\tObserved = (".(join ",", @obs).")\n\tExpected = (".(join ",", @new).")\n";
	  
	}
	else {
	  $numprob->bdiv($denprob);
	  $numprob->blog(10);
	  my $lod = $numprob->bstr();
	  if ($lod > $bestlod){		  
	    $bstr = "($hapi) Alpha = $alpha | LOD = $lod | Cases = $cnts{a} | Controls = $cnts{u}\n"; # | $numprob | $denprob\n";
	    $bestlod = $lod;
	    $besthap = @haps[$hapi];
	    $bestalpha = $alpha;
	  }
	}
      }
    }
  }

  ########################################################
  # Added by Mike to see progress of writing output file #
  ########################################################
  $out_count++;
  print "Line out: $out_count\n";

  my @vals = ("risk",$str,$besthap,$bestalpha,$bestlod,(join ",", @haps),$expstr,$obsstr);
  print OUT "".(join "\t", @vals)."\n";
  
}
close OUT;
print "\nMade $out\n";
#print "Time elapsed = ".(time()-$time)." seconds\n";


#######################################################################
# Calculate probability of obs (array of counts) given expected 
# (array of frequencies/probabilities) using multinomial distribution #
#######################################################################

sub getMultProb {
  my ($obs,$exp) = @_;
  my $t1 = time();  
  my $numer = Math::BigInt->new(1);
  my $den = 0;
  
  foreach my $obs (@$obs){
    $den += $obs;
  }

  foreach my $obs (@$obs){
    if (@factorials < $obs+1){
      die "No factorial for $obs\n";
    }
    my $fact = $factorials[$obs]; #Math::BigInt->bfac($obs);
    $numer->bmul($fact);
  }

  if (@factorials < $den-1){
    die "No factorial stored for $den\n";
  }

  $den = $factorials[$den];
  my $Mn = Math::BigFloat->new($den/$numer); 
  my $prob = Math::BigFloat->new(1);

  foreach my $i (0 .. @$obs-1){
    my $power = Math::BigFloat->new($$exp[$i]);
    $power->bpow($$obs[$i]);
    $prob->bmul($power);
  }  

  ### Round; preserve just 6 digits (improves speed)
  $prob->bround(6);
  $prob->bmul($Mn);
  
  ### Round; preserve just 4 digits (improves speed)
  $prob->bround(4);
  return $prob;
}


########################################################################
# This corrects for probabilities of zero using a pseudocount (Wilkes) #
########################################################################

sub correct {
  my ($arr, $cnt) = @_;
  my $chrs = $cnt*2;
  my $ok = 1;
  
  foreach my $x (@$arr){
    if ($x < 1/$chrs){
      $ok = 0;
    }
  }
  my @new;
  if ($ok){
    @new = @$arr;
  }
  else {
    my $nchrs = 0;
    foreach my $var (@$arr){
      push @new, ($var*$chrs)+1;
      $nchrs += ($var*$chrs)+1;
    }
    my $tot = 0;
    foreach my $n (@new){
      $n /= $nchrs;
      $n = sprintf("%.8f", $n);
      $tot += $n;

    }
  }

  return @new;
}


