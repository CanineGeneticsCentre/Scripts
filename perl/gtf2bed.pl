#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;

open IN, $ARGV[0] or die "Unable to open $ARGV[0]\n";

my %trans_hash = ();
my %gene_names = ();
my @exon_starts = ();
my @exon_lengths = ();
my $previous_trans = '';

my %biotype_colors = (
	lincRNA => '255,128,14',
	miRNA => '0,0,0',
	misc_RNA => '0,0,0',
	Mt_rRNA => '0,0,0',
	Mt_tRNA => '0,0,0',
	processed_pseudogene => '128,128,128',
	protein_coding => '31,119,180',
	pseudogene => '128,128,128',
	rRNA => '0,0,0',
	snoRNA => '0,0,0',
	snRNA => '0,0,0',
);


while(<IN>){
	next if ($_ =~ /^#/);
	chomp $_;
	my %atts = ();
	my ($chr, undef, $term, $start, $stop, undef, $strand, undef, $atts) = split("\t", $_);
	foreach my $a (split(/;\s?/, $atts)){
		my ($k, $v) = $a =~ /(.*) "(.*)"/;
		$atts{$k} = $v;
	}
	
	if ($term eq 'gene'){
		my $gene_id = $atts{gene_id};
		my $gene_name = defined $atts{gene_name} ? $atts{gene_name} : $gene_id;
		$gene_names{$gene_id} = $gene_name;
		next;
	}
		
	
	my $transcript_id = $atts{transcript_id};
	
	if(exists $trans_hash{$transcript_id}){
		if ($term eq 'exon'){
			if ($strand eq '+'){
				push (@exon_starts, ($start - $trans_hash{$transcript_id}{start}));
				push (@exon_lengths, ($stop - $start + 1));
			}
			elsif ($strand eq '-'){
				unshift (@exon_starts, ($start - $trans_hash{$transcript_id}{start}));
				unshift (@exon_lengths, ($stop - $start + 1));				
			}
		}
		if ($term eq 'start_codon'){
			if ($strand eq '+'){ $trans_hash{$transcript_id}{tc_start} = ($start-1); }
			elsif ($strand eq '-'){ $trans_hash{$transcript_id}{tc_stop} = $stop; }
		}
		if ($term eq 'stop_codon'){ 
			if ($strand eq '+'){ $trans_hash{$transcript_id}{tc_stop} = $stop; }
			elsif ($strand eq '-'){ $trans_hash{$transcript_id}{tc_start} = ($start-1); }
		}
		$previous_trans = $transcript_id;
	}
	else{
		if ($previous_trans ne ""){
			$trans_hash{$previous_trans}{exon_starts} = join(',', @exon_starts);
			$trans_hash{$previous_trans}{exon_lengths} = join(',', @exon_lengths);
			$trans_hash{$previous_trans}{exon_count} = scalar @exon_lengths;
		}
		my $gene_name = $gene_names{$atts{gene_id}};
		my $biotype = defined $atts{gene_biotype} ? $atts{gene_biotype} : "";
		$trans_hash{$transcript_id} = {
			id			=> $transcript_id,
			name		=> $gene_name,
			biotype		=> $biotype,
			chr			=> $chr,
			start		=> $start,
			stop		=> $stop,
			strand		=> $strand,
			tc_start	=> $stop,
			tc_stop		=> $stop
		};
		@exon_starts = ();
		@exon_lengths = ();
	}
}
$trans_hash{$previous_trans}{exon_starts} = join(',', @exon_starts);
$trans_hash{$previous_trans}{exon_lengths} = join(',', @exon_lengths);
$trans_hash{$previous_trans}{exon_count} = scalar @exon_lengths;
#print Dumper %trans_hash;
#print "\n\n";

print "#Ensembl Genes\n";
foreach my $tc_id (keys %trans_hash){
	print join ("\t", $trans_hash{$tc_id}{chr}, ($trans_hash{$tc_id}{start}-1), $trans_hash{$tc_id}{stop}, $trans_hash{$tc_id}{name}, 0, $trans_hash{$tc_id}{strand}, $trans_hash{$tc_id}{tc_start}, $trans_hash{$tc_id}{tc_stop}, $biotype_colors{$trans_hash{$tc_id}{biotype}}, $trans_hash{$tc_id}{exon_count}, $trans_hash{$tc_id}{exon_lengths}, $trans_hash{$tc_id}{exon_starts})."\n";
}


