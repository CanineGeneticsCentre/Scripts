#!/usr/local/bin/perl 

################################################################
#     Hap Mperm Convert                                        #
#                                                              #
#     Uses the output files from doing Haplotype Association   #
#     in PLINK and converts them into a plottable format       #
#                                                              #
################################################################


# Mike Boursnell Jun 2010
# Animal Health Trust
# Newmarket
# UK

use strict;
use Getopt::Long;
my $version					= "12";

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

#Arrays
my @item         			= ();
my @item_temp				= ();
my @newline_affected      	= ();
my @newline_normal      	= ();
my @snp_array				= ();
my @id_array				= ();
my @status_array			= ();
my @SNP_name_array			= ();

my @chromosome_array_mperm	= ();
my @position_array_mperm	= ();
my @win_string_array_mperm	= ();
my @haplotype_array_mperm	= ();

my @chromosome_array_map	= ();
my @position_array_map		= ();
my @win_string_array_map	= ();
my @haplotype_array_map		= ();

my @EMP1_array_mperm		= ();
my @EMP2_array_mperm		= ();
my @EMP1_array				= ();
my @EMP2_array				= ();
my @haplotype_array			= ();
my @minus_log_EMP2_array	= ();
my @new_line_array			= ();
my @win_per_chr_array		= ();

#Hash table for SNP names and positions
my %hash					= ();

my $show				= "false";

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
my $ans					= "";
my $snp_name_1			= "";
my $snp_name_2			= "";
my $win_string			= "";
my $haplotype			= "";
my $logfile				= "";
my $species				= "dog";
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
my $EMP1						= 0;
my $EMP2						= 0;
my $win_count_mapfile			= 0;
my $win_count_mpermfile			= 0;
my $array_count					= 0;
my $window						= 0;
my $last_window					= 1;
my $lowest_EMP1					= 0;
my $count						= 0;
my $haplotype_count				= 0;
my $lowest_count				= 0;
my $minus_log_EMP2				= 0;
my $chromosome_count			= 0;
my $row_count					= 0;
my $max_row_count				= 0;
my $window_count				= 0;
my $check_window_count			= 0;
my $win_per_chr_count			= 0;
my $total_no_chromosomes		= 0;



###############################
# process -file flag (if any) #
# and species flag (if any)   #
###############################

GetOptions("file=s"=>\$prefix,"hap_length=s"=>\$hap_length,"species=s"=>\$species,"map=s"=>\$mapfile,"out=s"=>\$outfile);


print "\n";
print "###################################################################\n";
print "#     Hap Mperm Convert                                           #\n";
print "#                                                                 #\n";
print "#     Uses the output files from doing Haplotype Association      #\n";
print "#     in PLINK and converts them into a plottable format.         #\n";
print "#                                                                 #\n";
print "###################################################################\n\n";

print "Version:\t$version\n\n";


#############################
# Get the input file name   #
#############################

if ($prefix eq "")
{
	print "Prefix for pair of PLINK impute.map and .impute.assoc.mperm files:  ";
	$prefix = <STDIN>;
	chomp $prefix;
	print"\n";
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
	if ($answer eq "4"){$species = "cat"}

	print"\n";
}


if ($hap_length eq "")
{
	print "Haplotype length:  ";
	$hap_length = <STDIN>;
	chomp $hap_length;
	print"\n";
}

if ($mapfile eq "") 
{
	print "PLINK map file:  ";
	$mapfile = <STDIN>;
	chomp $mapfile;
	print"\n";
}

$impute_mapfile = $prefix."_$hap_length.impute.map";
$mpermfile = $prefix."_$hap_length.impute.assoc.mperm";

if ($outfile eq "")
{
	$outfile = $prefix."_hap_mperm_{".$hap_length."}.txt";
}

if (! -e $mapfile){print "ERROR: Can't find PLINK MAP file $mapfile\n\n";exit;}

print "Prefix:     \t$prefix\n";
print "Species:    \t$species\n";
print "Map file:   \t$impute_mapfile\n";
print "Mperm file: \t$mpermfile\n";
print "Output file:\t$outfile\n\n";


############################
# Open files for output    #
############################
open (OUT, ">$outfile") || die "Cannot open $outfile";


############################
# Open the files for INPUT #
############################
open (MAP, "$mapfile") || die "Cannot open $mapfile";
open (IMPUTE_MAP, "$impute_mapfile") || die "Cannot open $impute_mapfile";
open (MPERM, "$mpermfile") || die "Cannot open $mpermfile";


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
	
	$hash{$position} = $SNP_name;
	
} # end of reading PLINK map file

##########################
# Read data in MAP file  # 
##########################

$snp_count = 0;
print "Reading impute MAP file...\n\n";

$last_window = 0;
$win_count_mapfile = 0;

$line_count=0;
$max_row_count =0;

while ($single_line = <IMPUTE_MAP>)
{
	chomp $single_line;

	$line_count = $line_count + 1;

	# split line at tabs or spaces
	@item=split(/\s+/,$single_line);

	
	$chromosome = $item[0];
	$win_string = $item[1];
	$position = $item[3];

	if (defined $hash{$position})
	{
		$SNP_name  = $hash{$position};
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


	######################
	# Process win_string #
	######################

	$win_pos = index $win_string,"WIN";
	$first_underscore_pos = index $win_string,"_";

	$window = substr ($win_string,$win_pos + 3, $first_underscore_pos - ($win_pos + 3));


	#######################################################
    # If this is a row with a new window number then      #
    # write it to a new array  			                  #
    #######################################################
	if ($window > $last_window)
	{
		$win_count_mapfile = $win_count_mapfile  + 1;

		$chromosome_array_map[$win_count_mapfile]=$chromosome;
		$win_string_array_map[$win_count_mapfile]=$window;
		$position_array_map[$win_count_mapfile]=$position;
		$snp_array[$win_count_mapfile]=$SNP_name;

		$line_count = 0;

		if ($chromosome == "1")
		{
			$max_row_count = $max_row_count + 1;
		}
		
	}

	$last_window = $window; 
	
} # reading imput map file


close IMPUTE_MAP;


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

print "Total number of windows in map file: $win_count_mapfile\n\n";


#########################
# Show part of map file #
#########################
print "Showing part of impute.map file...\n\n";

print "WIN\tCHR\tPOS   \tSNP\n";
for ($map_array_count=1;$map_array_count < 10;$map_array_count++)
{
	print "$win_string_array_map[$map_array_count]\t";
	print "$chromosome_array_map[$map_array_count]\t";
	print "$position_array_map[$map_array_count]\t";
	print "$snp_array[$map_array_count]\n";
}



###################################
# Read data in ASSOC. MPERM file  # 
###################################

$snp_count = 0;
print "\n\nReading ASSOC.MPERM file...\n\n";

$win_count_mpermfile = -1;
$last_window = 1;
$haplotype_count = 0;
$ans = "y";

while ($single_line = <MPERM>)
{
	chomp $single_line;

	################################
	# Split line at tabs or spaces #
	################################
	@item=split(/\s+/,$single_line);

	$chromosome = $item[1];
	$win_string = $item[2];
	$EMP1 = $item[3];
	$EMP2 = $item[4];


	##########################################
	# Process win_string in assoc.mperm file #
	##########################################

	$win_pos = index $win_string,"WIN";
	$first_underscore_pos = index $win_string,"_";
	$second_underscore_pos = index $win_string,"_",$first_underscore_pos + 1;

	$window = substr ($win_string,$win_pos + 3, $first_underscore_pos - ($win_pos + 3));
	$haplotype = substr ($win_string, $first_underscore_pos + 1, $second_underscore_pos - $first_underscore_pos - 1);

	
	#######################################################
    # If this is a row with a new window number then      #
	# store all the values associated with the            #
	# previous group of haplotypes to the results sheet   #
	#######################################################

	if ($window > $last_window)
	{
		
		###########################################################
		# Count number of actual windows for this chromosome      #
		# (Remember they are not always consecutive.  SOme window #
		# numbers can be missed out)                              #
		###########################################################
		$win_per_chr_count=$win_per_chr_count + 1;
		$win_per_chr_array[$chromosome] = $win_per_chr_array[$chromosome] + 1;

		##########################
		# Find lowest EMP1 value #
		##########################
		$lowest_EMP1 = 1;
		$lowest_count = 1;
		for ($count=1; $count <=$haplotype_count; $count++)
		{
			if ($EMP1_array[$count] < $lowest_EMP1)
			{
				if ($show eq "true"){print "LOWEST\n"}
				$lowest_EMP1 = $EMP1_array[$count];
				$lowest_count = $count;
			}
		}
		
		
		
		
		#############################
		# Calculate minus log EMP2  #
		#############################
        if ($EMP2_array[$lowest_count] > 0)
		{		
			$minus_log_EMP2 = (log($EMP2_array[$lowest_count]) / log(10)) * -1;
		}



		##################################
		# Convert haplotype 1234 to ACGT #
		##################################
		$haplotype_array[$lowest_count] =~ tr/1234/ACGT/;



		########################################
		# Store lowest value line in the array #
		########################################

		$win_count_mpermfile = $win_count_mpermfile + 1;

		$chromosome_array_mperm[$win_count_mpermfile]=$last_chromosome;
		$win_string_array_mperm[$win_count_mpermfile]=$last_window;
		$haplotype_array_mperm[$win_count_mpermfile]=$haplotype_array[$lowest_count];
		$EMP1_array_mperm[$win_count_mpermfile]=$lowest_EMP1;
		$EMP2_array_mperm[$win_count_mpermfile]=$EMP2_array[$lowest_count];
		$minus_log_EMP2_array[$win_count_mpermfile] = $minus_log_EMP2;

		if ($show eq "true")
		{
			print "win_count_mpermfile: $win_count_mpermfile\n";
			print "lowest_count:  $lowest_count\n";
			print "win_string_array_mperm: $chromosome_array_mperm[$win_count_mpermfile]\n";
			print "chromosome_array_mperm: $chromosome_array_mperm[$win_count_mpermfile]\n";
			print "haplotype_array_mperm: $haplotype_array_mperm[$win_count_mpermfile]\n";
			print "EMP1_array_mperm: $EMP1_array_mperm[$win_count_mpermfile]\n";
			print "EMP2_array_mperm: $EMP2_array_mperm[$win_count_mpermfile]\n";
			print "minus_log_EMP2_array: $minus_log_EMP2_array[$win_count_mpermfile]\n";
			
			$ans = <STDIN>;
			
		}

		########################
		# Reset arrays to zero #
		########################

		for ($count = 0; $count <= 50; $count++)
		{
			$EMP1_array[$count] = 0;
			$EMP2_array[$count] = 0;
			$haplotype_array[$count] = 0;
		}
		$haplotype_count = 0;

	} 

	#########################################################
	# Start making a list of the haplotypes for each window #
	# (Also store EMP1 and EMP2 for each haplotype)         #
	#########################################################

	$haplotype_count = $haplotype_count + 1;

	$haplotype_array[$haplotype_count] = $haplotype;
	$EMP1_array[$haplotype_count] = $EMP1;
	$EMP2_array[$haplotype_count] = $EMP2;

	$last_window = $window;
	$last_chromosome = $chromosome;

}

#####################
# After end of file #
#####################
#print "After end of file process last row...\n\n";

		##########################
		# Find lowest EMP1 value #
		##########################
		$lowest_EMP1 = 1;
		for ($count=1; $count <=$haplotype_count; $count++)
		{
			if ($EMP1_array[$count] < $lowest_EMP1)
			{
				$lowest_EMP1 = $EMP1_array[$count];
				$lowest_count = $count;
			}
		}

		#############################
		# Calculate minus log EMP2  #
		#############################
        	if ($EMP2_array[$lowest_count] > 0)
		{		
			$minus_log_EMP2 = (log($EMP2_array[$lowest_count]) / log(10)) * -1;
		}


		##################################
		# Convert haplotype 1234 to ACGT #
		##################################
		$haplotype_array[$lowest_count] =~ tr/1234/ACGT/;


		########################################
		# Store lowest value line in the array #
		########################################

		$win_count_mpermfile = $win_count_mpermfile + 1;

		$chromosome_array_mperm[$win_count_mpermfile]=$last_chromosome;
		$win_string_array_mperm[$win_count_mpermfile]=$last_window;
		$haplotype_array_mperm[$win_count_mpermfile]=$haplotype_array[$lowest_count];
		$EMP1_array_mperm[$win_count_mpermfile]=$lowest_EMP1;
		$EMP2_array_mperm[$win_count_mpermfile]=$EMP2_array[$lowest_count];
		$minus_log_EMP2_array[$win_count_mpermfile] = $minus_log_EMP2;


		########################
		# Reset arrays to zero #
		########################

		for ($count = 0; $count <= 50; $count++)
		{
			$EMP1_array[$count] = 0;
			$EMP2_array[$count] = 0;
			$haplotype_array[$count] = 0;
		}
		$haplotype_count = 0;

close MPERM;

#################################
# Show part of ASSOC.MPERM file #
#################################
print "Showing part of ASSOC.MPERM file...\n\n";

print "WIN\tCHR\tHAP\tEMP1\tEMP2\tlogEMP2\tPOS\tSNP\n";
for ($array_count=1;$array_count < 10;$array_count++)
{
	print "$win_string_array_mperm[$array_count]\t";
	print "$chromosome_array_mperm[$array_count]\t";
	print "$haplotype_array_mperm[$array_count]\t";
	print "$EMP1_array_mperm[$array_count]\t";
	print "$EMP2_array_mperm[$array_count]\t";
	print "$minus_log_EMP2_array[$array_count]\t";
	print "$position_array_map[$array_count]\t";
	print "$snp_array[$array_count]\n";
}


print "\n\n$win_count_mapfile windows in the MAP file\n";
print "$win_count_mpermfile windows in the ASSOC.MPERM file\n\n";


##########################################
# Work through each chromosome, adding a #
# bit to each newline as you work down   #
##########################################

$array_count =0;
$line_count = 0;


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
				$new_line_array[$row_count] = $new_line_array[$row_count]."$win_string_array_mperm[$array_count]";
				#$new_line_array[$row_count] = $new_line_array[$row_count]."$chromosome_array_mperm[$array_count]";
			}
			if ($chromosome_count > 1)
			{ 
				$new_line_array[$row_count] = $new_line_array[$row_count]."\t$win_string_array_mperm[$array_count]";
				#$new_line_array[$row_count] = $new_line_array[$row_count]."\t$chromosome_array_mperm[$array_count]";
			}
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$chromosome_array_mperm[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$snp_array[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$position_array_map[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$haplotype_array_mperm[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$minus_log_EMP2_array[$array_count]";
			


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
		print OUT "WIN\tCHR$chromosome_count\tSNP\tPOS\tHAP\tlogEMP2\t";

	}

	if ($is_y_chromosome eq "yes")
	{
		print OUT "WIN\tCHRY\tSNP\tPOS\tHAP\tlog EMP2\t";
	}

	print OUT "WIN\tCHRL\tSNP\tPOS\tHAP\tlog EMP2\n";


######################################
# Now the rest of the file           #
######################################

for ($row_count = 1;$row_count <= $max_row_count;$row_count++)
{
	print OUT "$new_line_array[$row_count]\n";
}



print "######################################\n";
print "# Finished HAP MPERM file conversion #\n";
print "######################################\n\n";

print "Output file: $outfile\n\n";

print "(This file can be used directly in the Excel sheet Hap_Mperm_plot.)\n\n";

###########################
#     Close all files     #
###########################
close MAP;
close MPERM;
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

