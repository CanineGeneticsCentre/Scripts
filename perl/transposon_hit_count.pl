# !/usr/bin/perl -w

###############################################################################    
#	transposon hit count      						                          #     
#									                                          #
#	Works out when reads in a SAM hit inside genes or not.                    #
#   For use as a replacment for dexseq_count.py                               #
#                                                                             #
#   For Amy Charbonneau                                                       # 
###############################################################################


#############################
# Mike Boursnell Aug 2013   #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use warnings;
use Term::ANSIColor;

# Version of program
my $version						= "5";

# File names
my $sam_file					= "";
my $gff_file					= "";
my $output_file					= ""; # Output file for list of genes and hits
my $track_file					= ""; # GFF3 file created for use as a track in IGV
my $track_file_genes			= ""; # GFF3 file created for use as a track in IGV to show genes only
my $track_file_transposons		= ""; # GFF3 file created for use as a track in IGV to show transposons only

# Colours for IGV tracks
my $transposon_hit_colour		= "#bc0ebe"; # Purple
my $transposon_miss_colour		= "#0c1d9d"; # Dark blue
my $essential_gene_colour		= "#a91526"; # Red
my $non_essential_gene_colour	= "#15a91f"; # Green

# Boolean
my $hit_gene					= "false"; # Does transposon hit gene?

my $single_line					= "";
my $next_line					= "";
my $prefix						= "";
my $answer						= "";
my $qname 						= "";
my $flag 						= "";
my $sam_strand					= "";
my $chromosome   				= "";
my $position     				= "";
my $mapq     					= "";
my $cigar   					= "";
my $rnext	 					= "";
my $pnext 						= "";
my $tlen   						= "";
my $seq     					= "";
my $qual     					= "";
my $transposon_id				= "";
my $flag_number					= 0;
my $total_no_hits				= 0;
my $transposon_count			= 0;
my $percentage_done				= 0;

#GFF file
my $attribute_string			= "";
my $strand						= "";
my $feature						= "";
my $attribute					= "";
my $locus_tag					= "";
my $start						= "";
my $end							= "";
my $attribute_array_size		= 0;
my $equals_pos					= 0;
my $gene_count					= 0;
my $total_no_genes				= 0;


# SAM file
my $flag_hex					= "";
my $flag_set					= "";
my $read_is_mapped				= ""; # false or true, depends on FLAG bit 2
my $read_length					= 0;
my $no_of_lines_sam_file		= 0;
my $no_of_lines_gff_file		= 0;
my $no_of_mapped_lines			= 0; # No of lines up to first non-mapped read
my $array_size					= 0;
my $line_count					= 0;
my $array_count					= 0;

#TRACK file
my $end_position				= 0;

#Arrays
my @start_array					= (); # stores start position of gene from GFF file
my @end_array					= (); # stores end position of gene from GFF file
my @locus_tag_array				= (); # stores locus_tag (i.e. name of gene)
my @gene_strand_array			= (); # stores strand of gene from GFF file
my @hit_count_array				= ();
my @file_array					= ();
my @item						= ();
my @attribute_array				= (); # from GFF file

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "#######################################\n";
print color 'bold white';
print "          transposon_hit_count         \n";
print color 'bold magenta';
print "#######################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This program uses a GFF file and a SAM file as input.\n\n";
print "    It counts the number of read starts in the SAM file  \n";
print "    that fall within the genes specified in the GFF file.\n\n";

print "    In this version the genes used in the GFF file have 'gene' in the 'type' column,\n";
print "    and the gene name used is in 'locus_tag' of the attributes.\n\n";

print color 'reset';
print "\n";

#####################
# Get GFF file name #
#####################
until (-e $gff_file)
{
	print "Name of the GFF file:                    ";
	$gff_file = <STDIN>;
	chomp $gff_file;
	
	if ($gff_file eq "ls")
	{
		print "\n";
		system ("ls *.gff");
		print "\n";
	}
	if ($gff_file ne "ls")
	{
		if (! -e "$gff_file"){print "\nFile doesn't exist. Try again...    \n";}
	}	
}

#####################
# Get SAM file name #
#####################
until (-e $sam_file)
{
	print "\nName of the SAM file:                    ";
	
	$sam_file = <STDIN>;
	chomp $sam_file;
	
	if ($sam_file eq "ls")
	{
		print "\n";
		system ("ls *.sam");
		print "\n";
	}
	if ($sam_file ne "ls")
	{
		if (! -e "$sam_file"){print "\nFile doesn't exist. Try again...    \n";}
	}	
}


############################
# Get name for output file #
############################
print "\nChoose a name for the output file:       ";
$output_file = <STDIN>;
chomp $output_file;

if (index($output_file,".") == -1)
{
	$output_file = $output_file.".out";
}

########################################
# Make name for GFF file for IGV track #
########################################
$prefix = &get_prefix ($output_file);
$track_file = "$prefix".".gff3";
$track_file_genes = "$prefix"."_genes.gff3";
$track_file_transposons = "$prefix"."_transposons.gff3";


#############################################
# Read GFF file into an array (file_array)  #
#############################################

open (IN_GFF, "$gff_file") || die "Cannot open $gff_file\n\n";
	@file_array = <IN_GFF>;
	$no_of_lines_gff_file = scalar @file_array;
close IN_GFF;


###########################################
# Looping through file_array for GFF file #
###########################################

print "\n\n#########################\n";
print "#  Processing GFF file  #\n";
print "#########################\n\n";

for ($line_count = 0;$line_count <= $no_of_lines_gff_file - 1;$line_count++)
{
	$single_line = $file_array[$line_count];
	&chomp_all($single_line);
	
	@item=split(/\t/,$single_line);
	$array_size = scalar @item;
	
	# If line doesn't start with hash (ie past header lines)
	if (index($single_line,"#") != 0)
	{
		$feature = $item[2];
		$start = $item[3];
		$end = $item[4];
		$strand = $item[6];
		$attribute_string = $item[8];
		
		if ($feature eq "gene")
		{
			# Split the attribute string at semicolons
			@attribute_array = split (";",$attribute_string);
			
			####################################################
			# Loop through all elements of the attribute array #
			####################################################
			for $attribute (@attribute_array)
			{
				if (index($attribute,"locus_tag") > -1)
				{
					$gene_count = $gene_count + 1;
					$hit_count_array[$gene_count] = 0;
					
					$equals_pos = index($attribute,"=");
					$locus_tag = substr($attribute,$equals_pos + 1,99);
					
					$start_array[$gene_count] = $start;
					$end_array[$gene_count] = $end;
					$locus_tag_array[$gene_count] = $locus_tag;
					
					# Store strand of gene for IGV track in gene_strand_array
					$gene_strand_array[$gene_count] = $strand;
				}
			}
		} # if gene
	} # if line does not start with a hash #
	
} # $line_count loop


$total_no_genes = $gene_count;

print "GFF file $gff_file analysed\n\n";
print "$gene_count genes found\n\n";

print " > Press return to continue\n\n";
$answer = <STDIN>;


#############################################
# Read SAM file into an array (file_array)  #
#############################################

print "#########################\n";
print "#  Processing SAM file  #\n";
print "#########################\n\n";

open (IN_SAM, "$sam_file") || die "Cannot open $sam_file\n\n";
	@file_array = <IN_SAM>;
	$no_of_lines_sam_file = scalar @file_array;
close IN_SAM;

print "No. of lines in SAM file $sam_file: $no_of_lines_sam_file\n\n";


###################################################
# Open Track file for GFF3 file to be used in IGV #
###################################################
open (TRACK_BOTH, ">$track_file") || die "Cannot open $track_file\n\n";
open (TRACK_GENES, ">$track_file_genes") || die "Cannot open $track_file_genes\n\n";
open (TRACK_TRANSPOSONS, ">$track_file_transposons") || die "Cannot open $track_file_transposons\n\n";


#######################################
# Write headers for GFF track file    #
# (not critical what this is exactly) #
#######################################
print TRACK_BOTH "##gff-version 3\n";
#print TRACK_BOTH "#!gff-spec-version 1.20\n";
#print TRACK_BOTH "#!processor NCBI annotwriter\n";
#print TRACK_BOTH "##sequence-region NC_012471.1 1 2253793\n";
#print TRACK_BOTH "##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=553482\n";
print TRACK_BOTH "#track name=GENES and TRANSPOSONS\n";

print TRACK_GENES "##gff-version 3\n";
#print TRACK_GENES "#!gff-spec-version 1.20\n";
#print TRACK_GENES "#!processor NCBI annotwriter\n";
#print TRACK_GENES "##sequence-region NC_012471.1 1 2253793\n";
#print TRACK_GENES "##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=553482\n";
print TRACK_GENES "#track name=GENES\n";

print TRACK_TRANSPOSONS "##gff-version 3\n";
#print TRACK_TRANSPOSONS "#!gff-spec-version 1.20\n";
#print TRACK_TRANSPOSONS "#!processor NCBI annotwriter transposons\n";
#print TRACK_TRANSPOSONS "##sequence-region NC_012471.1 1 2253793\n";
#print TRACK_TRANSPOSONS "##species http://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=553482\n";
print TRACK_TRANSPOSONS "#track name=TRANSPOSONS\n";


###########################################
# Looping through file_array for SAM file #
###########################################

# Check how far mapped reads go down the file
for ($line_count = 0;$line_count<=$no_of_lines_sam_file - 1; $line_count++)
{
	$single_line = $file_array[$line_count];
	chomp $single_line;
	
	@item=split(/\t/,$single_line);
	$array_size = scalar @item;
	
	if ($array_size > 11)
	{
		$position = $item[3];
		if ($position <1)
		{
			$no_of_mapped_lines = $line_count + 1;
			last;
		}
	}
}

print "No. of mapped reads in SAM file: $no_of_mapped_lines\n\n";

for ($line_count = 0;$line_count<=$no_of_lines_sam_file - 1; $line_count++)
{
	$single_line = $file_array[$line_count];
	chomp $single_line;
	
	@item=split(/\t/,$single_line);
	$array_size = scalar @item;
	

	########################
	# Deal with full lines #
	########################
	if ($array_size > 11)
	{
		#Read first 11 columns of SAM file (well, only ones we want)
		#$qname 			= $item[0];
		$flag 			= $item[1];
		#$chromosome   	= $item[2];
		$position     	= $item[3];
		$mapq     		= $item[4];
		$cigar   		= $item[5];
		#$rnext	 		= $item[6];
		#$pnext 		= $item[7];
		#$tlen   		= $item[8];
		$seq     		= $item[9];
		#$qual     		= $item[10];

		
		################################################################
		# Get read length so can correct for strand                    #
		# (if strand is reversed then start position is the other end) #
		################################################################
		$read_length = length($seq);
		
		
		#####################################
		# Read FLAG to determine strand     #
		# 0x0010 - should be read bitwise ? #
		#####################################
		$strand = "+"; # default
		$flag_number = $flag; # bitwise comparisons only work with numbers not strings
		
		#######################################################
		# Use bitwise comparison to look at bit 4 (2^^4 = 16) #
		# Bit 4 if to determine if the strand is reversed     #
		#######################################################
		if ($flag_number & 16){$strand = "-"} else {$strand = "+"}
		
		
		#######################################################
		# Use bitwise comparison to look at bit 4 (2^^4 = 16) #
		# Bit 4 if to determine if the strand is mapped       #
		#######################################################
		#if ($flag_number & 4){$read_is_mapped = "true"} else {$read_is_mapped = "false"}
		
		
		#############################################################
		# If reversed then add read length to get other end of read #
		#############################################################
		if ($strand eq "-")
		{
			$position = $position + $read_length - 1;  # IS this right?
		}
		
		if (($line_count % 10000) == 0)
		{
			$percentage_done = ($line_count / $no_of_mapped_lines) * 100;
			$percentage_done = sprintf('%.1f', $percentage_done);
			
			print "Line: $line_count \tPosition: $position \tPercent done: $percentage_done %\n";
		}
		

		############################################
		# Is position inside a gene from GFF file? #
		# Check through list of GFF genes          #
		############################################

		$hit_gene = "false";
		
		if ($position > 0 )
		{
			for ($gene_count = 1; $gene_count <=$total_no_genes; $gene_count++)
			{
				# If the start of the gene is further along than 'position' then exit loop with last
				if ($position < $start_array[$gene_count])
				{
					last;
				}
				
				if (($position >= $start_array[$gene_count]) && ($position <= $end_array[$gene_count]))
				{
					$hit_gene = "true";
					$hit_count_array[$gene_count] = $hit_count_array[$gene_count] + 1;
					
					$total_no_hits = $total_no_hits + 1;
				}
				
				
			}
			
			
			########################################################
			# Work out start and end of transposon marker in track #
			# (make the transposons lok 2 bases long in IGV)       #
			########################################################
			if ($strand eq "+")
			{
				$end_position = $position + 2;
			}
			if ($strand eq "-")
			{
				$end_position = $position;
				$position = $end_position - 2;
			}
		
			######################################
			# Add to TRACK file for IGV          #
			# (different colours if hit or miss) #
			######################################
			
			$transposon_count = $transposon_count + 1;
			
			$transposon_id = "T"."$transposon_count";
			
			if ($hit_gene eq "true"){
				print TRACK_BOTH "chr1\tRefseq\tHIT\t$position\t$end_position\t.\t$strand\t.\tID=$transposon_id;locus_tag=TT;color=$transposon_hit_colour\n";
			}
			if ($hit_gene eq "false"){
				print TRACK_BOTH "chr1\tRefseq\tHIT\t$position\t$end_position\t.\t$strand\t.\tID=$transposon_id;locus_tag=TT;color=$transposon_miss_colour\n";
			}
			if ($hit_gene eq "true"){
				print TRACK_TRANSPOSONS "chr1\tRefseq\tHIT\t$position\t$end_position\t.\t$strand\t.\tID=$transposon_id;color=$transposon_hit_colour\n";
			}
			if ($hit_gene eq "false"){
				print TRACK_TRANSPOSONS "chr1\tRefseq\tHIT\t$position\t$end_position\t.\t$strand\t.\tID=$transposon_id;color=$transposon_miss_colour\n";
			}

		}
		
	} # if array_size > 10
}


##############################################################
# Write the list of genes with hit counts to the output file #
# and the gene positions to the IGV track file.              #
##############################################################
open (OUT, ">$output_file") || die "Cannot open $output_file\n\n";

print "\n\nWriting hit counts to $output_file...\n\n";

print OUT "GENE\tHITS\n";

for ($gene_count = 1; $gene_count <= $total_no_genes; $gene_count++)
{
	print OUT "$gene_count\t$locus_tag_array[$gene_count]\t$hit_count_array[$gene_count]\n";
	
	$strand = $gene_strand_array[$gene_count];
	
	# If gene has been hit then it is non-essential
	if ($hit_count_array[$gene_count] > 0)
	{
		print TRACK_BOTH "chr1\tRefseq\tCDS\t$start_array[$gene_count]\t$end_array[$gene_count]\t.";
		print TRACK_BOTH "\t$strand\t.\tID=$locus_tag_array[$gene_count];color=$non_essential_gene_colour\n";
		
		print TRACK_GENES "chr1\tRefseq\tCDS\t$start_array[$gene_count]\t$end_array[$gene_count]\t.";
		print TRACK_GENES "\t$strand\t.\tID=$locus_tag_array[$gene_count];color=$non_essential_gene_colour\n";
		
	}
	
	# If gene has not been hit then it is marked by colour as "essential"	
	if ($hit_count_array[$gene_count] == 0)
	{
		print TRACK_BOTH "chr1\tRefseq\tCDS\t$start_array[$gene_count]\t$end_array[$gene_count]\t.";
		print TRACK_BOTH "\t$strand\t.\tID=$locus_tag_array[$gene_count];color=$essential_gene_colour\n";
		
		print TRACK_GENES "chr1\tRefseq\tCDS\t$start_array[$gene_count]\t$end_array[$gene_count]\t.";
		print TRACK_GENES "\t$strand\t.\tID=$locus_tag_array[$gene_count];color=$essential_gene_colour\n";
		
	}
}

print OUT "Total number of transposon hits:  \t$transposon_count\n";
print OUT "Total number of hits in genes :   \t$total_no_hits\n\n";

close OUT;
close TRACK;
close TRACK_BOTH;
close TRACK_TRANSPOSONS;


#####################################
# Tell user final results files etc #
#####################################
print "\n\n";
print "############################################\n";
print "#     Finished counting transposon hits    #\n";
print "############################################\n\n";

print "Input GFF file:                       \t$gff_file\n";
print "Input SAM file:                       \t$sam_file\n";
print "Output transposon hits file:          \t$output_file\n\n";

print "Output GFF track file:                \t$track_file\n";
print "Output GFF gene track file:           \t$track_file_genes\n";
print "Output GFF transposon track file:     \t$track_file_transposons\n\n";

print "Total number of transposon hits:      \t$transposon_count\n";
print "Total number of hits in genes :       \t$total_no_hits\n\n";

exit;



##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

foreach (@_) {s/\n//g}  
foreach (@_) {s/\r//g}  

}


####################################################################
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test" from "test.out")                             #
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
