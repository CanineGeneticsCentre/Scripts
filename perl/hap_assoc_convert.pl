#!/usr/local/bin/perl 

################################################################
#     Hap Assoc Convert                                        #
#                                                              #
#     Uses the output files from doing Haplotype Association   #
#     in PLINK and converts them into a plottable format       #
#                                                              #
################################################################


# Mike Boursnell April 2013
# Animal Health Trust
# Newmarket
# UK

use strict;
use Getopt::Long;
my $version					= 11;
my $show					= "false";

#####################################
# Define variables                  #
#####################################

#File names
my $prefix       			= "";
my $mapfile					= "";
my $impute_mapfile			= "";
my $mpermfile				= "";
my $pedfile					= "";
my $outfile					= "";
my $phenofile				= "";
my $assoc_hap_file			= "";

#Arrays
my @item         			= ();
my @item_snp				= ();
my @newline_affected      	= ();
my @newline_normal      	= ();
my @snp_array				= ();
my @id_array				= ();
my @status_array			= ();
my @P_value_array_window	= ();
my @haplotype_array_window	= ();
my @chromosome_array		= ();
my @position_array			= ();
my @win_string_array		= ();
my @haplotype_array			= ();


my @P_value_array			= ();
my @haplotype_array			= ();
my @new_line_array			= ();
my @win_per_chr_array		= ();
my @minus_log_P_array		= ();

#Hash table for chromosome and positions
# (based on SNP name)
my %hash_chr				= ();
my %hash_pos				= ();

my $genotype	 		= "";
my $answer       		= "";
my $single_line			= "";
my $pedigree			= "";
my $id					= "";
my $father				= "";
my $mother				= "";
my $gender				= "";
my $status				= "";
my $SNP_name			= "";
my $chromosome			= "";
my $last_chromosome		= "";
my $position			= "";
my $last_position		= "";
my $ans					= "";
my $win_string			= "";
my $haplotype			= "";
my $logfile				= "";
my $species				= "";
my $is_y_chromosome		= "no";

my $hap_length					= 0; # Haplotype length
my $line_count    				= 1;
my $item_count    				= 0;
my $snp_count					= 0;
my $total_lines  				= 0;
my $total_items  				= 0;
my $total_items_on_first_row	= 0;
my $total_SNPs_on_first_row 	= 0;
my $total_SNPs					= 0;
my $snp_total_in_mapfile		= 0;
my $sample_1					= 0;
my $sample_2					= 0;
my $map_array_count				= 0;
my $first_underscore_pos		= 0;
my $second_underscore_pos		= 0;
my $win_pos						= 0;
my $stat						= 0;
my $win_count_assoc_hap			= 0;
my $array_count					= 0;
my $window						= 0;
my $last_window					= 1;
my $lowest_P				    = 0;
my $count						= 0;
my $haplotype_count				= 0;
my $lowest_count				= 0;
my $minus_log_P					= 0;
my $chromosome_count			= 0;
my $row_count					= 0;
my $max_row_count_chr1			= 0;
my $max_row_count				= 0;
my $window_count				= 0;
my $check_window_count			= 0;
my $total_no_chromosomes		= 0;
my $win_per_chr_count			= 0;

#New for assoc.hap
my $P_value						= 0;
my $haplotype_string			= "";
my $SNP_string					= "";
my $first_SNP					= "";
my $last_first_SNP				= "";

###############################
# process -file flag (if any) #
# and species flag (if any)   #
###############################

GetOptions("file=s"=>\$prefix,"hap_length=s"=>\$hap_length,"species=s"=>\$species,"map=s"=>\$mapfile,"out=s"=>\$outfile);


print "\n";
print "###################################################################\n";
print "#     Hap Assoc Convert                                           #\n";
print "#                                                                 #\n";
print "#     Uses the output files from doing Haplotype Association      #\n";
print "#     in PLINK and converts them into a plottable format.         #\n";
print "#                                                                 #\n";
print "###################################################################\n\n";

print "Version:\t$version\n\n";


if ($prefix ne "")
{
	print "OPTIONS FROM COMMAND LINE\n\n";
	print "prefix:  \t$prefix\n";
	print "hap_length:  \t$hap_length\n";
	print "species:  \t$species\n";
	print "mapfile:  \t$mapfile\n";
	print "outfile:  \t$outfile\n\n";
}


#############################
# Get the input file name   #
#############################
if ($prefix ne "")
{
	$assoc_hap_file = $prefix.".assoc.hap";
}

if ($prefix eq "")
{
	until (-e $assoc_hap_file)
	{
		print "\nPrefix for assoc.hap file:      ";
		$prefix = <STDIN>;
		chomp $prefix;
		$assoc_hap_file = $prefix.".assoc.hap";
		
		if ($prefix eq "ls"){print "\n";system ("ls *.assoc.hap")}
		
		if ($prefix ne "ls")
		{
			if (! -e $assoc_hap_file){print "\n\n>>>>>>>>  File $assoc_hap_file not found.  Try again.  <<<<<<<<\n\n";}
		}
	}
}	
		
if ($species eq "")
{
print "\nWhich species?\n\n";
	print "<1> Dog\n";
	print "<2> Horse\n";
	print "<3> Human\n";
	
	$answer = <STDIN>;
	chomp $answer;
	
	if ($answer eq "1"){$species = "dog"}
	if ($answer eq "2"){$species = "horse"}
	if ($answer eq "3"){$species = "human"}

	print"\n";
}

if ($hap_length == 0)
{
	print "\nHaplotype length:  ";
	$hap_length = <STDIN>;
	chomp $hap_length;
	print"\n";
}

if ($mapfile eq "") 
{
	until (-e $mapfile)
	{
		print "PLINK map file:      ";
		$mapfile = <STDIN>;
		chomp $mapfile;
		if (index($mapfile,".map") == -1)
		{
			$mapfile=$mapfile.".map";
		}
		
		if ($mapfile eq "ls.map"){print "\n";system ("ls *.map")}
		
		if ($mapfile ne "ls.map")
		{
			if (! -e $mapfile){print "\n\n>>>>>>>>  File $mapfile not found.  Try again.  <<<<<<<<\n\n";}
		}
	}
}


if ($outfile eq "")
{
	$outfile = $prefix."_hap_assoc_{".$hap_length."}.txt";
}

if (! -e $mapfile){print "ERROR: Can't find PLINK MAP file $mapfile\n\n";exit;}

print "Prefix:           \t$prefix\n";
print "Assoc hap file:   \t$assoc_hap_file\n";
print "Haplotype length: \t$hap_length\n";
print "Map file:         \t$mapfile\n";
print "Output file:      \t$outfile\n\n";


############################
# Open files for output    #
############################
open (OUT, ">$outfile") || die "Cannot open $outfile";


############################
# Open the files for INPUT #
############################
open (MAP, "$mapfile") || die "Cannot open $mapfile";
open (ASSOC_HAP, "$assoc_hap_file") || die "Cannot open $assoc_hap_file";



#######################
# Read PLINK map file #
#######################

print "Reading PLINK MAP file...\n\n";

while ($single_line = <MAP>)
{
	chomp $single_line;
	# split line at tabs or spaces
	@item=split(/\s+/,$single_line);
	
	$chromosome = $item[0];
	$SNP_name = $item[1];
	$position = $item[3];
	
	$hash_pos{$SNP_name} = $position;
	$hash_chr{$SNP_name} = $chromosome;
	
} # end of reading PLINK map file

close MAP;

##########################
# Read data in MAP file  # 
##########################

$snp_count = 0;
print "Reading assoc.hap file...\n\n";

$last_window = 0;
$win_count_assoc_hap = 0;

$line_count=0;
$max_row_count =0;
$max_row_count_chr1 =0;

while ($single_line = <ASSOC_HAP>)
{
	chomp $single_line;

	$line_count = $line_count + 1;

	# split line at tabs or spaces
	@item=split(/\s+/,$single_line);

	$win_string = $item[1];
	$haplotype_string = $item[2];
	$P_value = $item[7];
	$SNP_string = $item[8];
	
	# Split SNP string at "|"
	@item_snp = split(/\|/,$SNP_string);
	$first_SNP = $item_snp[0];
	
	##################################################
	# Ignore OMNIBUS windows by setting P-value high #
	##################################################
	if ($haplotype_string eq "OMNIBUS") {$P_value = 1}
	
	#############################################
	# Process win_string to get a window number #
	#############################################
	$window = substr ($win_string,3);
	
	###################################################
	# Get chromosome and position from first SNP name #
	# by using the hash lookup from the .map file     #
	###################################################
	if (defined $hash_chr{$first_SNP})
	{
		$chromosome  = $hash_chr{$first_SNP};
	}
	if (defined $hash_pos{$first_SNP})
	{
		$position  = $hash_pos{$first_SNP};
	}
	
	#################################
	# Deal with X and Y chromosomes #
	#################################
	if ($species eq "dog")
	{
		if ($chromosome eq "X"){$chromosome=39}
		if ($chromosome eq "Y"){$chromosome=40;$is_y_chromosome = "yes"}
	}
	if ($species eq "horse")
	{
		if ($chromosome eq "X"){$chromosome=32}
		if ($chromosome eq "Y"){$chromosome=33; $is_y_chromosome = "yes"}
	}
	if ($species eq "human")
	{
		if ($chromosome eq "X"){$chromosome=23}
		if ($chromosome eq "Y"){$chromosome=24; $is_y_chromosome = "yes"}
	}
	if ($species eq "cat")
	{
		if ($chromosome eq "X"){$chromosome=19}
		if ($chromosome eq "Y"){$chromosome=20; $is_y_chromosome = "yes"}
	}

	if ($show eq "true")
	{
		print "\nLine: $line_count\tWIN: $window\tCHR: $chromosome\tPOS: $position\tSNP: $first_SNP\tHAP: $haplotype_string\tP: $P_value\n\n";
	}

	if ($window % 5000 == 0){print "Window: $window\n"}
	
	#######################################################
    # If this is a row with a new window number then      #
    # write it to a new array  			                  #
    #######################################################
	if ($window > $last_window)
	{
		if ($show eq "true")
		{print "\n\t#           NEW WINDOW SIZE             #\n"}
		
		###########################################################
		# Count number of actual windows for this chromosome      #
		# (Remember they are not always consecutive.  Some window #
		# numbers can be missed out)                              #
		###########################################################
		$win_per_chr_count = $win_per_chr_count + 1;
		$win_per_chr_array[$chromosome] = $win_per_chr_array[$chromosome] + 1;
		
		######################################################
		# Store the highest number of rows seen for any      #
		# of the chromosomes.  This is needed for the blank  #
		# spaces in the output file                          #
		######################################################
		if ($win_per_chr_array[$chromosome] > $max_row_count)
		{
			$max_row_count = $win_per_chr_array[$chromosome];
		}
		
		##########################
		# Find lowest P value    #
		##########################
		$lowest_P = 1;
		$lowest_count = 1;
		for ($count=1; $count <=$haplotype_count; $count++)
		{
			if ($P_value_array_window[$count] < $lowest_P)
			{
				$lowest_P = $P_value_array_window[$count];
				$lowest_count = $count;
			}
		}
		
		#############################
		# Calculate minus log P     #
		#############################
        if ($P_value_array_window[$lowest_count] > 0)
		{		
			$minus_log_P = (log($P_value_array_window[$lowest_count]) / log(10)) * -1;
		}

		if ($show eq "true")
		{
			print "Lowest P-value: $lowest_P\n";
			print "Lowest -logP: $minus_log_P\n\n";
		}
		
		##################################
		# Convert haplotype 1234 to ACGT #
		##################################
		#$haplotype_array[$lowest_count] =~ tr/1234/ACGT/;

		
		#############################################
		# Store lowest value line in the array      #
		# to be used later to write the output file #
		#############################################

		$win_count_assoc_hap = $win_count_assoc_hap + 1;

		$win_string_array[$win_count_assoc_hap]=$last_window;
		$chromosome_array[$win_count_assoc_hap]=$last_chromosome;
		$position_array[$win_count_assoc_hap]=$last_position;
		$snp_array[$win_count_assoc_hap]=$last_first_SNP;
		$haplotype_array[$win_count_assoc_hap]=$haplotype_array_window[$lowest_count];
		$P_value_array[$win_count_assoc_hap]=$lowest_P;
		$minus_log_P_array[$win_count_assoc_hap] = $minus_log_P;

		if ($chromosome == "1")
		{
			$max_row_count_chr1 = $max_row_count_chr1 + 1;
		}
		
	
		if ($show eq "true")
		{
			print "\tWIN: $last_window\tCHR: $last_chromosome\tPOS: $last_position\tSNP: $last_first_SNP\tHAP: $haplotype_array_window[$lowest_count]\tP: $lowest_P\tlogP: $minus_log_P\n\n";
			
			print "\tLowest vales for saving to array:\n\n";
			
			print "\t\tHaplotype:      \t$haplotype_array[$win_count_assoc_hap]\n";
			print "\t\tLowest P:        \t$P_value_array[$win_count_assoc_hap]\n";
			print "\t\tLowest log P:    \t$minus_log_P_array[$win_count_assoc_hap]\n";
			print "\t==================================================================\n\n";
			
			$answer=<STDIN>;
			print "\n";
		}
		
		########################
		# Reset arrays to zero #
		########################

		for ($count = 0; $count <= 50; $count++)
		{
			$P_value_array_window[$count] = 0;
			$haplotype_array_window[$count] = 0;
		}
		$haplotype_count = 0;
		
	
	} # if ($window > $last_window)

	
	#########################################################
	# If window is same or different....                    #
	# Start making a list of the haplotypes for each window #
	# (Also store P-value for each haplotype)               #
	#########################################################

	$haplotype_count = $haplotype_count + 1;

	$haplotype_array_window[$haplotype_count] = $haplotype_string;
	$P_value_array_window[$haplotype_count] = $P_value;

	if ($show eq "true")
	{
		print "Haplotypes: ";
		for ($count = 1; $count <=$haplotype_count;$count++)
		{
			print "$count: $haplotype_array_window[$count]  ";
		}
		print "\n";
		print "P-values: ";
		for ($count = 1; $count <=$haplotype_count;$count++)
		{
			print "$count: $P_value_array_window[$count]  ";
		}
		print "\n";
	}
	
	####################################
	# Store the last mentions of these #
	####################################
	$last_window = $window;
	$last_chromosome = $chromosome;
	$last_position = $position;
	$last_first_SNP = $first_SNP;
	
} # reading assoc.hap  file


close ASSOC_HAP;


################################
# Decide number of chromosomes #
################################
if ($species eq "dog")
{
	if ($is_y_chromosome eq "yes"){$total_no_chromosomes = 40} else {$total_no_chromosomes = 39}
}
if ($species eq "horse")
{
	if ($is_y_chromosome eq "yes"){$total_no_chromosomes = 33} else {$total_no_chromosomes = 32}
}
if ($species eq "human")
{
	if ($is_y_chromosome eq "yes"){$total_no_chromosomes = 24} else {$total_no_chromosomes = 23}
}
if ($species eq "cat")
{
	if ($is_y_chromosome eq "yes"){$total_no_chromosomes = 20} else {$total_no_chromosomes = 19}
}


print "Highest chromosome: $total_no_chromosomes\n\n";

print "Total number of windows in map file: $win_count_assoc_hap\n\n";



#####################
# After end of file #
#####################
print "After end of file process last row...\n\n";

##########################
# Find lowest P value    #
##########################
$lowest_P = 1;
$lowest_count = 1;
for ($count=1; $count <=$haplotype_count; $count++)
{
	if ($P_value_array_window[$count] < $lowest_P)
	{
		$lowest_P = $P_value_array_window[$count];
		$lowest_count = $count;
	}
}

#############################
# Calculate minus log P     #
#############################
if ($P_value_array_window[$lowest_count] > 0)
{		
	$minus_log_P = (log($P_value_array_window[$lowest_count]) / log(10)) * -1;
}


##################################
# Convert haplotype 1234 to ACGT #
##################################
#$haplotype_array[$lowest_count] =~ tr/1234/ACGT/;


########################################
# Store lowest value line in the array #
########################################

$win_count_assoc_hap = $win_count_assoc_hap + 1;

$chromosome_array[$win_count_assoc_hap]=$last_chromosome;
$win_string_array[$win_count_assoc_hap]=$last_window;
$haplotype_array[$win_count_assoc_hap]=$haplotype_array[$lowest_count];
$P_value_array[$win_count_assoc_hap]=$lowest_P;
$minus_log_P_array[$win_count_assoc_hap] = $minus_log_P;


########################
# Reset arrays to zero #
########################

for ($count = 0; $count <= 50; $count++)
{
	$P_value_array_window[$count] = 0;
	$haplotype_array_window[$count] = 0;
}
$haplotype_count = 0;



#################################
# Show part of ASSOC.HAP file   #
#################################
print "Showing part of array from ASSOC.HAP file...\n\n";

print "WIN\tCHR\tHAP\tP\tlogP\tPOS\tSNP\n";
for ($array_count=1;$array_count <= 10;$array_count++)
{
	print "$win_string_array[$array_count]\t";
	print "$chromosome_array[$array_count]\t";
	print "$haplotype_array[$array_count]\t";
	print "$P_value_array[$array_count]\t";
	print "$minus_log_P_array[$array_count]\t";
	print "$position_array[$array_count]\t";
	print "$snp_array[$array_count]\n";
}


print "\n$win_count_assoc_hap windows in the ASSOC.HAP file\n\n";


##########################################
# Work through each chromosome, adding a #
# bit to each newline as you work down   #
##########################################

$array_count =0;
$line_count = 0;

print "Species: $species\tNo. of chromosomes:\t$total_no_chromosomes\tMaximum row count:\t$max_row_count\n\n";

for ($chromosome_count = 1;$chromosome_count<=$total_no_chromosomes;$chromosome_count++)
	{

	######################################################
	# Add the correct number of rows for each chromosome #
	######################################################
	
	for ($row_count = 1; $row_count <= $win_per_chr_array[$chromosome_count]; $row_count++)
		{
			$array_count = $array_count + 1;

			
			
			if ($chromosome_count == 1)
			{
				$new_line_array[$row_count] = $new_line_array[$row_count]."$win_string_array[$array_count]";
				#$new_line_array[$row_count] = $new_line_array[$row_count]."$chromosome_array[$array_count]";
			}
			if ($chromosome_count > 1)
			{ 
				$new_line_array[$row_count] = $new_line_array[$row_count]."\t$win_string_array[$array_count]";
				#$new_line_array[$row_count] = $new_line_array[$row_count]."\t$chromosome_array[$array_count]";
			}
			
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$chromosome_array[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$snp_array[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$position_array[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$haplotype_array[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$minus_log_P_array[$array_count]";
			
		}

	##################################################
	# Now fill up remaining newline arrays with tabs #
	##################################################
	for ($row_count = $win_per_chr_array[$chromosome_count] + 1; $row_count <= $max_row_count; $row_count++)
	{
		if ($chromosome_count == 1) {$new_line_array[$row_count] = $new_line_array[$row_count]."\t\t\t\t\t"}
		if ($chromosome_count > 1) {$new_line_array[$row_count] = $new_line_array[$row_count]."\t\t\t\t\t\t"}
	}
	
}

######################################################################
# Write from final array to an output file in columns for Excel 2003 #
# (One column per chromosome)                                        #
######################################################################

print "Writing to output file $outfile...\n\n\n";

######################################
# Headers across the top of the file #
######################################

for ($chromosome_count = 1;$chromosome_count<=$total_no_chromosomes - 1;$chromosome_count++)
{
	print OUT "WIN\tCHR$chromosome_count\tSNP\tPOS\tHAP\tlogP\t";

}

if ($is_y_chromosome eq "yes")
{
	print OUT "WIN\tCHRY\tSNP\tPOS\tHAP\tlogP\t";
}

print OUT "WIN\tCHRL\tSNP\tPOS\tHAP\tlogP\n";


######################################
# Now the rest of the file           #
######################################

for ($row_count = 1;$row_count <= $max_row_count;$row_count++)
{
	print OUT "$new_line_array[$row_count]\n";
}


print "######################################\n";
print "# Finished ASSOC HAP file conversion #\n";
print "######################################\n\n";

print "Output file: \t$outfile\n\n";

print "(This file can be used directly in the Excel sheet Hap_Assoc_plot.)\n\n";

###########################
#     Close all files     #
###########################


close OUT;

exit;

##########################################################################################


##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

foreach (@_) {s/\n//g}  
foreach (@_) {s/\r//g}  

}

