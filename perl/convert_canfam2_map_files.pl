#!/usr/bin/perl -w

use strict;
use Term::ANSIColor;

#Constants
my $version						= "4";

# Variables
my $canfam2_file				= "";
my $canfam3_file				= "";
my $conversion_file				= ""; # TEMP see below
my $single_line					= "";
my $chromosome 					= "";
my $position_cf2 				= "";
my $position_cf3 				= "";
my $centimorgans 				= "";
my $snp_id_cf2 					= "";
my $snp_id_cf3 					= "";
my $answer						= "";
my $position_markers			= "";
my $rs_name						= "";
my $array_type					= ""; # SNP20 or HD array
my $conversion_file				= "";
my $prefix						= "";
my $renamed_input_file			= "";

my $line_count					= 0;
my $array_size					= 0;
my $undefined_count 			= 0;
my $defined_count 				= 0;
my $last_offset					= 0;
my $position_cf3_guess			= 0;
my $guess_error					= 0;
my $no_of_map_lines				= 0;

my @item						= ();
my @map_file_array				= ();

my %hash_position				= ();
my %hash_snp_id					= ();

my $hash_position_ref 			= ();

print color 'bold magenta';

print "\n\n";
print "#################################\n";
print color 'bold white';
print "    convert canfam2 map files       \n";
print color 'bold magenta';
print "#################################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";

print "  - This program converts plink canine MAP files \n\n";

print "    It make three types of conversion:\n\n";

print "    1. Converts canfam2 SNP20 arrays to canfam3 \n";
print "    2. Converts canfam2 HD arrays to canfam3 \n";
print "    3. Converts canfam3 HD arrays with the new manifest to a version with more SNP position assignments\n\n";

print color 'reset';



until (-e $canfam2_file)
{
	&print_message("canFam2 PLINK map file to convert: ","input");

	$canfam2_file = <STDIN>;
	chomp $canfam2_file;

	if ($canfam2_file eq "ls"){system "ls *.map"}

	if ($canfam2_file ne "ls")
	{
		if (! -e $canfam2_file)
		{
			print "\n\n  >>>>>>>>>   Can't find input file $canfam2_file   <<<<<<<<<<<\n\n";
		}
	}
}

#####################
# Count the markers #
#####################
open (IN, "$canfam2_file") || die "Cannot open input file $canfam2_file";
@map_file_array = <IN>;
close IN;
$no_of_map_lines = scalar @map_file_array;


########################
# Name the output file #
########################
$prefix = &get_prefix($canfam2_file);
$canfam3_file = $prefix."_cf3.map";
$renamed_input_file = $prefix."_cf2.map";


&print_message("Map file has $no_of_map_lines lines.  Which type of MAP file is this?","input");

print "  <1>  SNP20 array    (22,362 SNPs)\n";
print "  <2>  canineHD array canfam2 (173,662 SNPs)\n";
print "  <3>  canineHD array canfam3 (173,662 SNPs) - new canfam3 manifest\n\n";

$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1" ){$array_type = "SNP20"; $conversion_file = "/home/genetics/SNPs/SNP20_canfam3_conversion.txt";$renamed_input_file = $prefix."_cf2.map";}
if (substr($answer,0,1) eq "2" ){$array_type = "HD"; $conversion_file = "/home/genetics/SNPs/canfam2_canfam3_conversion.txt";$renamed_input_file = $prefix."_cf2.map";}
if (substr($answer,0,1) eq "3" ){$array_type = "HD_CF3"; $conversion_file = "/home/genetics/SNPs/canfam3_new_manifest_conversion.txt";$renamed_input_file = $prefix."_cf3_original.map";}


#####################################
# Copy original map file as CF2.map #
#####################################
system ("cp $canfam2_file $renamed_input_file");


#############################################
# Open canfam2 to canfam3 conversion file   #
# and create hash table for pos and snp_id  #
#############################################
print "\nReading in conversion_file...\n\n";

open (MARKERS, "$conversion_file") || die "Cannot open markers file $conversion_file";
while ($single_line = <MARKERS>) 
{

	chomp $single_line;
	@item=split(/\s+/,$single_line);

	# 4 columns #
	$chromosome = $item[0];
	$snp_id_cf2 = $item[1];
	$position_cf2 = $item[2];
	$snp_id_cf3 = $item[3];
	$position_cf3 = $item[4];

	$array_size = scalar @item;
	

	if ($array_size == 6)
	{
		$rs_name = $item[4];
		$snp_id_cf2 = $snp_id_cf2."_".$rs_name;
	}
	
	$hash_position{$snp_id_cf2} = $position_cf3;
	$hash_snp_id{$snp_id_cf2} = $snp_id_cf3;

}
close MARKERS;


###############################
# Open input and output files #
###############################
if ($array_type eq "SNP20"){print "Converting SNP20 MAP file $canfam2_file to $canfam3_file...\n\n";}
if ($array_type eq "HD"){print "Converting canineHD MAP file $canfam2_file to $canfam3_file...\n\n";}
if ($array_type eq "HD_CF3"){print "Converting canineHD (new canfam3 manifest) MAP file $canfam2_file to $canfam3_file...\n\n";}


open (IN, "$canfam2_file") || die "Cannot open input file $canfam2_file";
open (OUT, ">$canfam3_file") || die "Cannot open output file $canfam3_file";
open (UNDEFINED, ">undefined.out") || die "Cannot open output file undefined.out";

 while ($single_line = <IN>) {

	chomp $single_line;
	@item=split(/\s+/,$single_line);
	$chromosome = $item[0];
	$snp_id_cf2 = $item[1];
	$centimorgans = $item[2];
	$position_cf2 = $item[3];
	$line_count = $line_count + 1;
	


	##############################################
	# Look up new canFam3 position in hash table #
	##############################################
	if ((defined $hash_position{$snp_id_cf2}) && (defined $hash_snp_id{$snp_id_cf2}))
	{
		$position_cf3 = $hash_position{$snp_id_cf2};
		$snp_id_cf3 = $hash_snp_id{$snp_id_cf2};
		$defined_count = $defined_count + 1;
		
		print OUT "$chromosome\t$snp_id_cf3\t$centimorgans\t$position_cf3\n";
	}
	else
	{
		$position_cf3 = "ND";
		$undefined_count = $undefined_count + 1;
		print UNDEFINED "$chromosome\t$snp_id_cf2\t$centimorgans\t$position_cf3\n";
	}
	
	if (($line_count % 1000) == 0)
	{
		print "Line: $line_count\tcanFam2 position: $position_cf2";
		print "\tcanFam3 position: $position_cf3\n"
	}
}

close IN;
close OUT;
close UNDEFINED;


&print_message("FINISHED CONVERSION","message");

print "No SNPs converted:   \t$defined_count/$line_count\n";
print "No SNPs unconverted: \t$undefined_count/$line_count\n\n";

print color ' bold yellow';

print "Input canfam2 file:   \t$canfam2_file\n";
print "Output canfam3 file:  \t$canfam3_file\n\n";

print color 'reset';

print "(A copy of the original canfam2 file called $renamed_input_file has also been made)\n\n";


 
 
 exit;

 ######################################
# Subroutine to print screen message #
######################################
sub print_message
{
	my $message_length 	= "";
	my $pos_count		= 0;
	my $char			= "";
	
	my $message = $_[0];
	my $style = $_[1];
	
	$message_length = length($message);
	
	if ($style eq ""){$char = "#"}
	if ($style eq "input"){$char = "~"}
	if ($style eq "message"){$char = "#"}
	if ($style eq "warning"){$char = "!"}
	
	print "\n\n";
	print color ' bold yellow';
	if ($style eq "warning"){print color ' bold red'}
	if ($style eq "input"){print color ' bold white'}
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n$char    $message    $char\n";
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n\n";
	print color 'reset';

}#
 
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
