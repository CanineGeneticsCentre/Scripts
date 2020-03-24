#!/usr/bin env perl

$RUN_NAME=$ARGV[0];
$SAMPLE=$ARGV[1];

$OUTPUT_PATH=$ENV{'RESULTS'}."/$SAMPLE/$RUN_NAME";

open(TABLE, $OUTPUT_PATH."/".$SAMPLE.".table") or die "Unable to find file ".$OUTPUT_PATH."/".$SAMPLE.".table\n";
open (PED, ">".$OUTPUT_PATH."/".$SAMPLE.".ped");
open (MAP, ">".$OUTPUT_PATH."/".$SAMPLE.".map");

print "PED file = ".$OUTPUT_PATH."/".$SAMPLE.".ped\n";
print "MAP file = ".$OUTPUT_PATH."/".$SAMPLE.".map\n";

print PED join("\t", "FAM", $SAMPLE, 0,0,0,0);

my $first = 1;
while (<TABLE>){
	#CHROM   POS     hd_chip.ID      REF     ALT     FILTER  GS_32059.GT
	#1       14548   BICF2G630707759 A       <NON_REF>       homREF  A/A
	if ($first){ $first = 0; next; }
	chomp $_;
	my @tmp = split("\t", $_);
	my @genotype = split('/', $tmp[-1]);
	print PED "\t".$genotype[0]." ".$genotype[1];
	print MAP join("\t", $tmp[0], $tmp[2], 0, $tmp[1])."\n";
}
print PED "\n";

close TABLE;
close PED;
close MAP;