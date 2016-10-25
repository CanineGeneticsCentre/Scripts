#!/usr/bin/perl -w


use strict;
use Getopt::Std ;
use File::Basename ;

my $line_count			= 0;
my $hit_count			= 0;
my $lines_to_keep		= 0;
my $keep_count			= 0;
my $position 			= 0;
my $position_start		= 0;

my $name				= "";
my $chr_01				= "";
my $chr_02				= "";
my $sam_file			= "";
my $output_file			= "";
my $single_line			= "";
my $answer				= "";
my $chromosome			= "";

my @item         		= ();

print "\n\nMAKE SMALLER SAM\n\n";

print "\nInput SAM file:      ";
$sam_file = <STDIN>;
chomp $sam_file;

print "\nOuput SAM file:      ";
$output_file = <STDIN>;
chomp $output_file;

print "Which chromosome: ";
$chromosome = <STDIN>;
chomp $chromosome;

if (index($chromosome,"chr") == -1) {$chromosome = "chr"."$chromosome";}

chomp $chromosome;

print "How many lines with this chromosome ($chromosome) do you want to copy to the new file: ";
$lines_to_keep = <STDIN>;
chomp $lines_to_keep;

print "Starting at which base position: ";
$position_start = <STDIN>;
chomp $position_start;

############################
# Open the files for INPUT #
############################

open (IN, "$sam_file") || die "Cannot open $sam_file";

############################
# Open files for output    #
############################

open (OUT, ">$output_file") || die "Cannot open $output_file";


while ($single_line = <IN>)
{

	chomp $single_line;

	@item=split(/\t/,$single_line);

	$name = $item[0];
	$chr_01 = $item[2];
	$position = $item[3];
	$chr_02 = $item[6];
	
	$line_count = $line_count + 1;
	
	if (($line_count % 50000) == 0 )
	{
		print "Line: $line_count\tPosition: $position\n";
	}
	
	
	if ($line_count < 58)
	{
		print OUT "$single_line\n";
	}

	if ($chr_01 eq $chromosome)
	{

		if ($chr_02 eq "=")
		{
			$hit_count = $hit_count + 1;
			
			print "Position: $position\tPosition start: $position_start\n";
			
			if (($position > $position_start) && ($keep_count <= $lines_to_keep))
			{
				print "    Transferred: $keep_count\n";
				print OUT "$single_line\n";
				$keep_count = $keep_count + 1;
			}
			
			
			if ($hit_count > $lines_to_keep)
			{
				close IN;
				close OUT;

				print "\n\nFINISHED\n\n";

				print "$keep_count lines were transferred\n";
				print "Output file is $output_file\n";
				exit;
			}
		
		} #
	}
}

close IN;
close OUT;

print "\n\nFINISHED\n\n";

print "There were $hit_count hits\n\n";
print "$keep_count lines were transferred\n";
print "Output file is $output_file\n";

		