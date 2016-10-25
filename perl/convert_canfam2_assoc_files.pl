#!/usr/bin/perl -w

use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;

#Constants
my $version						= "5";

# Variables
my $canfam2_file				= "";
my $canfam3_file				= "";
my $canfam3_markers_file		= "/home/genetics/SNPs/canfam3/canfam3_canineHD_SNPs.txt";
my $single_line					= "";
my $chromosome 					= "";
my $last_chromosome				= "0";
my $position_cf2 				= "";
my $position_cf3 				= "";
my $last_position_cf3			= "";
my $centimorgans 				= "";
my $snp_id_cf2					= "";
my $snp_id_cf3					= "";
my $answer						= "";
my $position_cf3_markers		= "";
my $rs_name						= "";
my $col_4						= "";
my $col_5						= "";
my $col_6						= "";
my $col_7						= "";
my $chi_sq						= "";
my $file_type					= "assoc"; # This could be mperm as well.
my $prefix						= "";

my $line_count					= 0;
my $array_size					= 0;
my $undefined_count 			= 0;
my $defined_count 				= 0;
my $last_offset					= 0;
my $position_cf3_guess			= 0;
my $guess_error					= 0;
my $p_value						= 0;
my $OR							= 0;
my $emp1						= 0;
my $emp2						= 0;
my $give_last_position_count	= 0;
my $give_1000_count				= 0;
my $total_SNPs_in_output		= 0;
my $total_SNPs_in_input			= 0;

my @item						= ();
my %hash						= ();


print color 'reset';
print color 'bold magenta';
print "\n\n";
print "#######################################\n";
print color 'bold white';
print "      convert_canfam2_assoc_files      \n";
print color 'bold magenta';
print "#######################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This program converts PLINK assoc and assoc.mperm files\n";
print "    from canFam2 positions to canFam3 positions\n\n";

print "    NB. Not all canFam2 SNPs map onto canFam3, but any that do not are assigned\n";
print "    a position based on the adjacent canFam3 SNPs, so are not missed out.\n\n";

print "    The data file is $canfam3_markers_file\n\n";

print color 'reset';
print "\n\n";

print "canFam2 PLINK file to convert: ";
$canfam2_file = <STDIN>;
chomp $canfam2_file;

###########################
# What type of file is it #
###########################
if (index($canfam2_file,".assoc") > -1){$file_type = "assoc";$prefix = &get_prefix($canfam2_file);$canfam3_file = $prefix."_cf3.assoc"}
if (index($canfam2_file,".assoc.mperm") > -1){$file_type = "mperm";$prefix = &get_prefix($canfam3_file);$canfam3_file = $prefix."_cf3.mpermassoc"}

if (! -e $canfam2_file)
{
	print "\n\n";
	print "###################\n";
	print "# FILE NOT FOUND  #\n";
	print "###################\n\n";
	print "Can't find input file $canfam2_file\n\n";
	exit;
}	


###############################
# Open canfam3 markers file   #
# and create hash table       #
###############################
print "\nReading in canfam3 markers file /home/genetics/SNPs/canfam3/canfam3_canineHD_SNPs.txt...\n\n";

open (MARKERS, "$canfam3_markers_file") || die "Cannot open markers file $canfam3_markers_file";
while ($single_line = <MARKERS>) {

# 	Part of the markers file...

#	MARKER		CHROM	DIST	RS
#	TIGRP2P259	1		249580	rs8993730

	chomp $single_line;
	@item=split(/\s+/,$single_line);

	$snp_id_cf3 = $item[0];
	$chromosome = $item[1];
	$position_cf3_markers = $item[2];
	
	$array_size = scalar @item;
	
	if ($array_size == 3){$rs_name = ""}
	if ($array_size == 4)
	{
		$rs_name = $item[3];
		$snp_id_cf3 = $snp_id_cf3."_".$rs_name;
	}
	
	$hash{$snp_id_cf3} = $position_cf3_markers; # This allows you to look up the SNP by name and find its position in the CF3 markers file
}
close MARKERS;


###############################
# Open input and output files #
###############################
print "Converting $canfam2_file to $canfam3_file...\n\n";

open (IN, "$canfam2_file") || die "Cannot open input file $canfam2_file";
open (OUT, ">$canfam3_file") || die "Cannot open output file $canfam3_file";
open (UNDEFINED, ">undefined.out") || die "Cannot open output file undefined.out";

print UNDEFINED "CHR\tSNP\tPOS\tP\n";

 while ($single_line = <IN>) {
	$line_count = $line_count + 1;
	chomp $single_line;
	@item=split(/\s+/,$single_line);
	
	$array_size = scalar @item;
	
	if ($file_type eq "assoc")
	{
		$chromosome = $item[1];
		$snp_id_cf2 = $item[2];
		$position_cf2 = $item[3];
		$col_4 = $item[4];
		$col_5 = $item[5];
		$col_6 = $item[6];
		$col_7 = $item[7];
		$chi_sq = $item[8];
		$p_value = $item[9];
		$OR = $item[10];
		
	}
	
	if ($file_type eq "mperm")
	{
		$chromosome = $item[1];
		$snp_id_cf2 = $item[2];
		$emp1 = $item[3];
		$emp2 = $item[4];
		
		#print "$single_line\n";
		#print "CHR:$chromosome\tSNP:$snp_id\tEMP1:$emp1\tEMP2:$emp2\n";
		#$answer=<STDIN>;
	}
	
	#####################################
	# Check on first line (header line) #
	#####################################
	if ($line_count ==1)
	{
		if ($file_type eq "assoc")
		{
			print OUT "$chromosome   $snp_id_cf2   $position_cf2   $col_4   $col_5   $col_6   $col_7   $chi_sq   $p_value   $OR   INFO\n";
			
			if ($p_value ne "P")
			{
				print "##############\n";
				print "# FILE ERROR #\n";
				print "##############\n\n";
				print "P value column seems to be in the wrong place\n\n";
				print "The header is $p_value rather then 'P'\n\n";
				print "See Mike\n\n";
				close IN;
				close OUT;
				close UNDEFINED;
				exit;
			}
		} # assoc file
		if ($file_type eq "mperm")
		{
			print OUT "$chromosome   $snp_id_cf2   $emp1   $emp2   INFO\n";
		
			if ($emp1 ne "EMP1")
			{
				print "##############\n";
				print "# FILE ERROR #\n";
				print "##############\n\n";
				print "EMP1 column seems to be in the wrong place\n\n";
				print "The header is $emp1 rather then 'EMP1'\n\n";
				print "See Mike\n\n";
				close IN;
				close OUT;
				close UNDEFINED;
				exit;
			}
		} # assoc file
		
	} # first line only
	

	##############################################
	# Look up new canFam3 position in hash table #
	##############################################


	#################################################################
	# If the canfam2 SNP is found (by SNP name) in the canfam3 list #
	#################################################################
	if (defined $hash{$snp_id_cf2})
	{
		$position_cf3 = $hash{$snp_id_cf2};
		$defined_count = $defined_count + 1;
		$last_position_cf3 = $position_cf3 + 1;
		
		if ($file_type eq "assoc")
		{
			print OUT "$chromosome   $snp_id_cf2   $position_cf3   $col_4   $col_5   $col_6   $col_7   $chi_sq   $p_value   $OR   DEFINED\n";
		}
		if ($file_type eq "mperm")
		{
			print OUT "$chromosome   $snp_id_cf2   $emp1   $emp2   DEFINED\n";
		}
		
		# Store last chromosome so we can can see when it changes
		$last_chromosome = $chromosome;
	}

	#######################################################
	# If the canfam2 SNP is NOT found in the canfam3 list #
	#######################################################
	else
	{
		if ($line_count > 1)
		{
			if ($chromosome eq $last_chromosome)
			{
				$position_cf3 = $last_position_cf3;  # actually last position plus 1
				if ($file_type eq "assoc")
				{
					print OUT "$chromosome   $snp_id_cf2   $position_cf3   $col_4   $col_5   $col_6   $col_7   $chi_sq   $p_value   $OR   LAST\n";
				} # assoc file
				
				if ($file_type eq "mperm")
				{
					print OUT "$chromosome   $snp_id_cf2   $emp1   $emp2   LAST\n";
				} # mperm file

				$give_last_position_count = $give_last_position_count + 1;

			} # If ($chromosome eq $last_chromosome)
			
			if ($chromosome ne $last_chromosome)
			{
				$position_cf3 = "1000";
				if ($file_type eq "assoc")
				{
					print OUT "$chromosome   $snp_id_cf2   $position_cf3   $col_4   $col_5   $col_6   $col_7   $chi_sq   $p_value   $OR   THOU\n";
				} # assoc file
				
				if ($file_type eq "mperm")
				{
					print OUT "$chromosome   $snp_id_cf2   $emp1   $emp2   THOU\n";
				} # mperm file

				$give_1000_count = $give_1000_count + 1;
			} # If ($chromosome eq $last_chromosome)
			
			
			
			$undefined_count = $undefined_count + 1;
			print UNDEFINED "$chromosome\t$snp_id_cf2\t$position_cf3\t$p_value\n";
		}
	}
	
	if (($line_count % 10000) == 0)
	{
		if ($file_type eq "assoc")
		{
			print "Line: $line_count\tcanFam2 position: $position_cf2\tcanFam3 position: $position_cf3\n"
		}
		if ($file_type eq "mperm")
		{
			print "Line: $line_count\t$snp_id_cf2\tcanFam3 position: $position_cf3\n"
		}
	}
}

close IN;
close OUT;
close UNDEFINED;


print "\n";
print "#############################\n";
print "#     FINISHED CONVERSION   #\n";
print "#############################\n\n";

print "Input canfam2 file:  \t$canfam2_file\n\n";
print "Output canfam3 file: \t$canfam3_file\n\n";

$total_SNPs_in_output = $defined_count + $give_last_position_count + $give_1000_count;
$total_SNPs_in_input = $line_count - 1;

print "No. of SNPs in canfam2 file:                                               \t$total_SNPs_in_input\n\n";

print "No. of SNPs with direct match in canfam3:                                  \t$defined_count\n";
print "No. of SNPs with no direct match in canfam3, but used last position + 1:   \t$give_last_position_count\n";
print "No. of SNPs with no direct match in canfam3, but used position as 1000:    \t$give_1000_count\n\n";

print "No. of SNPs converted one way or another                                   \t$total_SNPs_in_output\n\n";
 
 
 
 exit;

 ####################################################################
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test" from "test.fasta")                           #
####################################################################

sub get_prefix
{
	my $filename = "";

	$filename = $_[0];
	if (index($filename,".") > 0)
	{
		$filename = substr($filename, 0, index($filename,"."));
	}
	if (index($filename,".") == -1)
	{
		$filename = $filename;
	}
}
