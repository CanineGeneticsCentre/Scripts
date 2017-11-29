#!/usr/bin/env perl

##################################################################################
#									                                             #
#	convert_ens_ucsc				                                             #
#									                                             #
#	This PERL script converts chrs of files from UCSC/Ensembl format 			 #
#   to Ensembl/UCSC format	 													 #
#									                                             #
##################################################################################

##############################
# Ellen Schofield Dec 2016   #
# Animal Health Trust        #
# Newmarket                  #
# UK                         #
# ellen.schofield@aht.org.uk #
##############################

use strict;
use warnings;
use File::Basename;
use Getopt::Long;
use PerlIO::gzip;


my ($conversion_file, $in, $out);
my $length = 24;
my ($bed, $vcf, $gff, $gtf);
my %conversion = ();

GetOptions (
	"length=i"	=> \$length,		# numeric
	"file=s"	=> \$conversion_file,	# string
	"in=s"		=> \$in,				# string
	"out=s"		=> \$out,				# string
	"bed"		=> \$bed,			# flag
	"vcf" 		=> \$vcf,			# flag
	"gff"		=> \$gff,			# flag
	"gtf"		=> \$gtf,			# flag
);

unless (defined $conversion_file && -e $conversion_file){
	print STDERR "You need to define a conversion file - either ucsc2ensembl or ensmebl2ucsc\n";
	print STDERR "\n$0 --file <cf3_ucsc2ensembl.txt / cf3_ensembl2ucsc.txt> --in <input file> [--bed --vcf --gff -gtf --out <output file>]\n\n";
	exit(0);
}

if (defined $out){
	close STDOUT;
	#open (STDOUT, ">$out") or die "Could not write to $out: $!";
	if ($in =~ /gz$/){ unless($out =~ /gz$/){ $out .= '.gz';} open STDOUT, ">:gzip", $out or die "Could not write to $out: $!"; }
	else{ open (STDOUT, ">$out") or die "Could not write to $out: $!"; }
}
else{
	if ($in =~ /gz$/){
		warn "Unable to write output to gzipped format like input format. You must supply an output file for this option\n";
		exit(0);
	}
}

open (FH, $conversion_file) or die "Unable to open file $conversion_file\n";
while(<FH>){
	chomp $_;
	my ($key, $value) = split("\t", $_);
	$conversion{$key} = $value;
}
close FH;

if (defined $bed){ convert_bed($in);}
elsif (defined $gff){ convert_gff($in);}
elsif (defined $gtf){ convert_gtf($in);}
elsif (defined $vcf){ convert_vcf($in);}
else{
	print STDERR "You need to define an input/output filetype - bed vcf gff gtf...\n";
	print STDERR "\n$0 --file <cf3_ucsc2ensembl.txt / cf3_ensembl2ucsc.txt> --in <input file> [--bed --vcf --gff -gtf --out <output file>]\n\n";
	exit(0);
}

close STDOUT;


sub _convert_bed_gtf_gff{
	my $file = $_[0];
	my $i = 0;
	open(IN, $in) or die "Unable to open file $in\n";
	while (<IN>){
		$i++;
		if ($_ =~ /^#/){ print $_; next;}
		my @a = split("\t", $_);
		if (!exists $conversion{$a[0]}){
			warn "ERROR (line $i) - $a[0] not defined in conversion table... line skipped\n";
			next;
		}
		$a[0] = $conversion{$a[0]};
		print join("\t", @a);
	}
	close IN;
}
	
sub convert_vcf{
	print STDERR "Converting VCF file - ".$_[0]."\n";
	my $file = $_[0];
	my $i = 0;
	if ($in =~ /gz$/){
		open(IN, "gunzip -c $file | ") or die "Unable to open file $file\n";
	}
	else { open(IN, $file) or die "Unable to open file $file\n"; }
	
	while(<IN>){
		$i++;
		if ($_ =~ /^##contig=/){
			my ($chr) = $_ =~ /##contig=<ID=(.*),length/;
			if (!exists $conversion{$chr}){
				warn "ERROR (line $i) - $chr not defined in conversion table... line skipped\n";
				next;
			}
			$_ =~ s/$chr/$conversion{$chr}/;
			print $_;
			next;
		}
		elsif ($_ =~ /^#/){  print $_; next;}
		
		my @a = split("\t", $_);
		if (!exists $conversion{$a[0]}){
			warn "ERROR (line $i) - $a[0] not defined in conversion table... line skipped\n";
			next;
		}
		$a[0] = $conversion{$a[0]};
		print join("\t", @a);
	}
	close IN;
}
	
sub convert_bed{
	print STDERR "Converting BED file - ".$_[0]."\n";
	_convert_bed_gtf_gff($_[0]);
}
	
sub convert_gff{
	print STDERR "Converting GFF file - ".$_[0]."\n";
	_convert_bed_gtf_gff($_[0]);
}
	
sub convert_gtf{
	print STDERR "Converting GTF file - ".$_[0]."\n";
	_convert_bed_gtf_gff($_[0]);
}