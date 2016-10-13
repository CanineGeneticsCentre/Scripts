#!/usr/local/bin/perl 

################################################################
#     SNPTEST Convert v3                                       #
#                                                              #
#     Uses the output files fromSNPTEST                        #
#     in PLINK and converts them into a plottable format       #
#                                                              #
################################################################


#############################
# Mike Boursnell Aug 2010   #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# mike.boursnell@aht.org.uk #
#############################

use strict;
use Getopt::Long;

#####################################
# Define variables                  #
#####################################
my @item         		= ();

my @SNP_name_array		= ();
my @chromosome_array		= ();
my @position_array		= ();

my @snptest_add_array		= ();
my @snptest_add_array		= ();
my @snptest_dom_array		= ();
my @snptest_rec_array		= ();
my @snptest_gen_array		= ();
my @snptest_het_array		= ();
my @new_line_array		= ();
my @new_line_2_array		= ();
my @lines_per_chr_array		= ();


my $infile       		= "";
my $outfile			= "";
my $outfile_2			= "";
my $answer       		= "";
my $single_line			= "";
my $pedigree			= "";
my $id				= "";
my $father			= "";
my $mother			= "";
my $gender			= "";
my $status			= "";
my $SNP_name			= "";
my $chromosome			= "";
my $last_chromosome		= "";
my $position			= "";
my $ans				= "";
my $snp_name_1			= "";
my $snp_name_2			= "";
my $win_string			= "";
my $haplotype			= "";
my $logfile			= "";
my $species			= "dog";
my $is_y_chromosome		= "no";
my $file_type			= "";
my $root_infile			= "";


my $line_count    		= 1;
my $item_count    		= 0;
my $snp_count			= 0;
my $total_lines  		= 0;
my $total_items  		= 0;
my $total_items_on_first_row	= 0;
my $total_SNPs_on_first_row 	= 0;
my $total_SNPs			= 0;
my $snp_total_in_mapfile	= 0;
my $sample_1			= 0;
my $sample_2			= 0;
my $map_array_count		= 0;
my $first_underscore_pos	= 0;
my $second_underscore_pos	= 0;
my $win_pos			= 0;
my $stat			= 0;
my $EMP1			= 0;
my $EMP2			= 0;
my $win_count_mapfile		= 0;
my $win_count_mpermfile		= 0;
my $array_count			= 0;
my $window			= 0;
my $last_window			= 1;
my $lowest_EMP1			= 0;
my $count			= 0;
my $haplotype_count		= 0;
my $lowest_count		= 0;
my $minus_log_EMP2		= 0;
my $chromosome_count		= 0;
my $row_count			= 0;
my $max_row_count		= 0;
my $window_count		= 0;
my $check_window_count		= 0;
my $max_chromosome_count	= 0;
my $snptest_add		= 0;
my $snptest_dom		= 0;
my $snptest_rec		= 0;
my $snptest_gen		= 0;
my $snptest_het		= 0;


###############################
# process -file flag (if any) #
# and species flag (if any)   #
###############################

GetOptions("file=s"=>\$infile,"species=s"=>\$species);


print "\n";
print "###########################################################################\n";
print "#     SNPTEST Convert v4                                                  #\n";
print "#                                                                         #\n";
print "#     Converts SNPTEST output files frequentist.out and bayesian.out      #\n";
print "#     so they can be plotted in Excel 2003.                               #\n";
print "#                                                                         #\n";
print "###########################################################################\n\n";


#############################
# Get the input file name   #
#############################

if ($infile eq "")
{
	print "Name of SNPTEST file:  ";
	$infile = <STDIN>;
	chomp $infile;
	print"\n";
}

###########################################
# Get root of filename so the output file #
# can be named with xl before the dot.    #
###########################################

$root_infile = $infile;

if (index ($infile,".") > 0)
{
	$root_infile = substr($infile, 0, rindex($infile,"."));
}


$outfile = $root_infile.".xlout";
$outfile_2 = $root_infile."2.xlout";


############################
# Open files for output    #
############################

open (OUT, ">$outfile") || die "Cannot open $outfile";
open (OUT2, ">$outfile_2") || die "Cannot open $outfile_2";

############################
# Open the file for INPUT  #
############################

open (IN, "$infile") || die "Cannot open $infile";



############################
# Read data in input file  # 
############################

$snp_count = 0;
print "Reading input file...\n";

$last_window = 0;
$win_count_mapfile = 0;

$line_count=0;
$max_row_count =0;

$array_count = 0;


while ($single_line = <IN>)
{

	chomp $single_line;

	$line_count = $line_count + 1;

	##############################################
	# Work out file type bayesian or frequentist #
	##############################################

	if ($line_count == 1)
	{
		if (index($single_line,"frequen") > 0)
		{
			$file_type = "frequentist";
		}
		if (index($single_line,"bayes") > 0)
		{
			$file_type = "bayesian";
		}
		if ($file_type eq "")
		{
			print "\n\nFile is not Bayesian or Frequentist SNPTEST file.\n\n\n";
			exit;

		}
		print "\nFile type is $file_type\n\n";
	}
	
	################################
	# split line at tabs or spaces #
	################################
	@item=split(/\s+/,$single_line);


	#################################################
	# Assign elements of the array to our variables #
	#################################################

	$array_count = $array_count + 1;

	if ($array_count % 10000 ==0)
	{
		print "$array_count\n";
	}

	$chromosome = $item[0];
	$SNP_name = $item[1];
	$position = $item[2];

	$snptest_add = $item[41];
	$snptest_dom = $item[42];
	$snptest_rec = $item[43];
	$snptest_gen = $item[44];
	$snptest_het = $item[45];

	if ($file_type eq "frequentist")
	{

		#####################
		# Take minus log 10 #
		#####################
        	if ($snptest_add > 0)
		{		
			$snptest_add = (log($snptest_add) / log(10)) * -1;
		}
		else
		{
			$snptest_add = 0;
		}

        	if ($snptest_dom > 0)
		{
			$snptest_dom = (log($snptest_dom) / log(10)) * -1;
		}		
		else
		{
			$snptest_dom = 0;
		}

        	if ($snptest_rec > 0)
		{		
			$snptest_rec = (log($snptest_rec) / log(10)) * -1;
		}		
		else
		{
			$snptest_rec = 0;
		}

        	if ($snptest_gen > 0)
		{		
			$snptest_gen = (log($snptest_gen) / log(10)) * -1;
		}		
		else
		{
			$snptest_gen = 0;
		}

        	if ($snptest_het > 0)
		{		
			$snptest_het = (log($snptest_het) / log(10)) * -1;
		}		
		else
		{
			$snptest_het = 0;
		}

	}


	###################
	# Store in arrays #
	###################

	$chromosome_array[$array_count] = $chromosome;
	$SNP_name_array[$array_count] = $SNP_name;
	$position_array[$array_count] = $position;


	$snptest_add_array[$array_count] = $snptest_add;
	$snptest_dom_array[$array_count] = $snptest_dom;
	$snptest_rec_array[$array_count] = $snptest_rec;
	$snptest_gen_array[$array_count] = $snptest_gen;
	$snptest_het_array[$array_count] = $snptest_het;




	#################################################
	# Store the number of lines for each chromosome #
	#################################################

	$lines_per_chr_array[$chromosome] = $lines_per_chr_array[$chromosome] + 1;

	if ($chromosome == "1")
	{
		$max_row_count = $max_row_count + 1;
	}


	if ($species == "dog")
	{
		if ($chromosome eq "X"){$chromosome=39}
		if ($chromosome eq "Y")
			{
				$chromosome=40; 
				$is_y_chromosome = "yes";
			}
	}
	if ($species == "horse")
	{
		if ($chromosome eq "X"){$chromosome=32}
		if ($chromosome eq "Y"){$chromosome=330; $is_y_chromosome = "yes"}
	}



	
        
}

close IN;

print "\nNumber of lines: $array_count\n";

if ($species eq "dog")
{
	if ($is_y_chromosome eq "yes"){$max_chromosome_count = 40} else {$max_chromosome_count = 39}
}
if ($species eq "horse")
{
	if ($is_y_chromosome eq "yes"){$max_chromosome_count = 33} else {$max_chromosome_count = 32}
}




##########################################
# Work through each chromosome, adding a #
# bit to each newline as you work down   #
##########################################


$array_count = 1;
$line_count = 0;

print "\nWriting output files... \n\n";


for ($chromosome_count = 1;$chromosome_count<=$max_chromosome_count;$chromosome_count++)
{


	######################################################
	# Add the correct number of rows for each chromosome #
	######################################################
	for ($row_count = 1; $row_count <= $lines_per_chr_array[$chromosome_count]; $row_count++)
		{
			$array_count = $array_count + 1;


			#########################################################
			# newline_array is the first file with add, dom and rec #
			##########################################################
			if ($chromosome_count == 1)
			{
				$new_line_array[$row_count] = $new_line_array[$row_count].($array_count - 1);
				$new_line_2_array[$row_count] = $new_line_2_array[$row_count].($array_count - 1);
			}
			if ($chromosome_count > 1)
			{ 
				$new_line_array[$row_count] = $new_line_array[$row_count]."\t".($array_count - 1);
				$new_line_2_array[$row_count] = $new_line_2_array[$row_count]."\t".($array_count - 1);
			}

			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$SNP_name_array[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$position_array[$array_count]";

			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$snptest_add_array[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$snptest_dom_array[$array_count]";
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t$snptest_rec_array[$array_count]";


			##########################################################
			# newline_2_array is the second file with gen and het    #
			##########################################################

			$new_line_2_array[$row_count] = $new_line_2_array[$row_count]."\t$SNP_name_array[$array_count]";
			$new_line_2_array[$row_count] = $new_line_2_array[$row_count]."\t$position_array[$array_count]";

			$new_line_2_array[$row_count] = $new_line_2_array[$row_count]."\t$snptest_gen_array[$array_count]";
			$new_line_2_array[$row_count] = $new_line_2_array[$row_count]."\t$snptest_het_array[$array_count]";
			$new_line_2_array[$row_count] = $new_line_2_array[$row_count]."\t";



		}

	##################################################
	# Now fill up remaining newline arrays with tabs #
	##################################################
	for ($row_count = $lines_per_chr_array[$chromosome_count] + 1; $row_count <= $max_row_count; $row_count++)
	{
		if ($chromosome_count == 1) 
		{
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t\t\t\t\t";
			$new_line_2_array[$row_count] = $new_line_2_array[$row_count]."\t\t\t\t\t";
		}

		if ($chromosome_count > 1) 
		{
			$new_line_array[$row_count] = $new_line_array[$row_count]."\t\t\t\t\t\t";
			$new_line_2_array[$row_count] = $new_line_2_array[$row_count]."\t\t\t\t\t\t";
		}
	}
	
}


######################################################################
# Write from final array to an output file in columns for Excel 2003 #
# (One column per chromosome)                                        #
######################################################################


for ($chromosome_count = 1;$chromosome_count<=$max_chromosome_count - 1;$chromosome_count++)
	{
		if ($file_type eq "bayesian")
		{
			print OUT "CHR$chromosome_count\tSNP\tPOS\tBAY_ADD\tBAY_DOM\tBAY_REC\t";
			print OUT2 "CHR$chromosome_count\tSNP\tPOS\tBAY_GEN\tBAY_HET\t\t";
		}
		if ($file_type eq "frequentist")
		{
			print OUT "CHR$chromosome_count\tSNP\tPOS\tFRQ_ADD\tFRQ_DOM\tFRQ_REC\t";
			print OUT2 "CHR$chromosome_count\tSNP\tPOS\tFRQ_GEN\tFRQ_HET\t\t";
		}

	}

	if ($is_y_chromosome eq "yes")
	{
		if ($file_type eq "bayesian")
		{
			print OUT "CHRY\tSNP\tPOS\tBAY_ADD\tBAY_DOM\tBAY_REC\t";
			print OUT2 "CHRY\tSNP\tPOS\tBAY_GEN\tBAY_HET\t\t";
		}
		if ($file_type eq "frequentist")
		{
			print OUT "CHRY\tSNP\tPOS\tFRQ_ADD\tFRQ_DOM\tFRQ_REC\t";
			print OUT2 "CHRY\tSNP\tPOS\tFRQ_GEN\tFRQ_HET\t\t";
		}
	}

	if ($file_type eq "bayesian")
	{
		print OUT "CHRL\tSNP\tPOS\tBAY_ADD\tBAY_DOM\tBAY_REC\n";
		print OUT2 "CHRL\tSNP\tPOS\tBAY_GEN\tBAY_HET\t\n";
	}
	if ($file_type eq "frequentist")
	{
		print OUT "CHRL\tSNP\tPOS\tFRQ_ADD\tFRQ_DOM\tFRQ_REC\n";
		print OUT2 "CHRL\tSNP\tPOS\tFRQ_GEN\tFRQ_HET\t\n";
	}



for ($row_count = 1;$row_count <= $max_row_count;$row_count++)
{
	print OUT "$new_line_array[$row_count]\n";
	print OUT2 "$new_line_2_array[$row_count]\n";
}


print "############################\n";
print "# Finished file conversion #\n";
print "############################\n\n";

print "Output file with ADD, DOM and REC: \t$outfile\n";
print "Output file with GEN and HET:      \t$outfile_2\n\n";

print "(These files can be used directly in the Excel 2003 sheet SNPTEST_plot HD)\n\n";

###########################
#     Close all files     #
###########################
close IN;
close OUT;
close OUT2;

