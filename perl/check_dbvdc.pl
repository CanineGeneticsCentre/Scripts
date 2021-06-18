#!/usr/bin/env perl

use strict;
use warnings;

#use List::Util qw[min max];

my %snps = ();

open SAMPLES, "samples.list";
my @samples = <SAMPLES>;
my $samples_count = scalar(@samples);
close SAMPLES;

open(RESULTS, "snps.out");
while (<RESULTS>){
	chomp $_;
	my ($chr, $pos) = split("\t", $_);
  $snps{$chr.':'.$pos} = $_;
}
close RESULTS;

open(IN, $ARGV[0]);
open(OUT, ">>dbvdc.".$samples_count.".out");
while (<IN>){
	chomp $_;
	my ($chr, $pos) = split("\t", $_);
  if(defined $snps{$chr.':'.$pos}){
    my $row = $snps{$chr.':'.$pos};
    my @cols = split("\t", $row);
    my $ref_base = $cols[2];
    my $alt_base = $cols[3];

    print OUT join("\t", splice(@cols, 0, 4));

    #foreach my $gt (splice(@cols, 5)){
    foreach my $gt (@cols){
      $gt =~ s/0/$ref_base/g;
      $gt =~ s/1/$alt_base/g;
      $gt =~ s/\./X/g;
      $gt =~ s/\//\t/g;
      print OUT "\t".$gt;
    }
    print OUT "\n";
  }
  else{
    # No result for the CHR:POS
    print OUT join("\t", $chr, $pos). "\t-" x (2*$samples_count + 2)."\n";
  }

}
close IN;
close OUT;




# open(IN, $ARGV[0]);
# while (<IN>){
# 	chomp $_;
# 	my ($chr, $pos) = split("\t", $_);
# 	my @rows = `tabix $TOSSO $chr:$pos-$pos`;
# 	print STDERR "tabix $TOSSO $chr:$pos-$pos\n";
# 	if (scalar @rows == 0){
# 		#print STDERR "Investigate tabix +/- 100bp to look at indels...\n\ttabix $TOSSO $chr:".($pos-100)."-".($pos+100)."\n";
# 		print join("\t", $chr, $pos). "\t-" x (2*(scalar @samples) + 2)."\n";
# #		my $dsPos = $pos-100;
# #		my $upPos = $pos+100;
# #		my @extra_rows = `tabix $TOSSO $chr:$dsPos-$upPos`;
# #		my %seen = ();
# #		foreach my $r (@extra_rows){
# #			my @cols = split("\t", $extra_rows[0]);
# #			next if (exists $seed{$cols[2]});
# #			$seen{$cols[2]}++;
# #			my $bp_length = max(length($cols[2]),length($cols[3]));
# #			print join ("\t", $pos, $cols[1], $bp_length)."\n" if ($cols[1] + $bp_length >= $pos && $cols[1] < $pos);
# #		}
# #		exit(0);
# 	}
# 	else{
# 		#next;
# 		my $r = $rows[0];
# 		#foreach my $r (@rows){
# 			my @cols = split("\t", $r);
# 			my $ref_base = $cols[$REF];
# 			my $alt_base = $cols[$ALT];
# 			print join("\t", $chr, $pos, $cols[$REF], $cols[$ALT]);
# 			foreach my $sample (splice(@cols, $first_sample , ($last_sample-$first_sample))){
# 				my @info = split(':', $sample);
# 				my $gt = $info[0];
# 				$gt =~ s/0/$cols[$REF]/g;
# 				$gt =~ s/1/$cols[$ALT]/g;
# 				$gt =~ s/\./X/g;
# 				$gt =~ s/\//\t/g;
# 				print "\t".$gt;
# 			}
# 			print "\n";
# 		#}
# 	}
# }
# close IN;
