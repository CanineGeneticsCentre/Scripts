#!/usr/bin/perl

################################
# which_sample_creates_peak v6 #
################################

#use Term::ReadKey;
use File::Basename ;
use strict;
use Getopt::Long;


#####################################
# Define variables                  #
#####################################

my $version					= "2";

my $pedfile 				= "";
my $logfile					= "";
my $mapfile					= "";
my $outfile					= "temp";
my $resultsfile				= "";
my $snpfile					= "";
my $command					= "";
my $single_line				= "";
my $pedigree				= "";
my $sample					= "";
my $remove_file				= "";
my $snp_name				= "";
my $ans						= "";
my $status					= "";
my $prefix					= "";
my $arg						= "";
my $flag					= "";
my $genotype				= "";
my $chromosome				= "";
my $position				= "";
my $too_many_samples 		= "False";


my $chisq			= 0;
my $chisq_for_this_snp		= 0;
my $max_chisq			= 0;
my $highest_chisq		= 0;
my $highest_chisq_relative	= 0;
my $total_chisq			= 0;
my $global_highest_chisq	= 0;
my $best_sample_relative	= 0;
my $best_sample_absolute	= 0;
my $best_sample			= 0;
my $initial_chisq		= 0;
my $p_value			= 0;
my $minus_log_p_value		= 0;
my $snp_found			= 0;

my $linecount			= 0;
my $total_samples		= 0;
my $sample_count		= 0;
my $no_to_remove		= 0;
my $count			= 0;
my $no_of_passes		= 6;
my $chisq_reduction		= 0;
my $flag_count			= 0;
my $total_snps			= 0;
my $snp_count			= 0;
my $i				= 0;
my $j				= 0;
my $no_selected_SNPs		= 0;

my @item         		= ();
my @sample_array		= ();
my @aht_no_array		= ();
my @best_remove			= ();
my @best_remove_relative	= ();
my @chisq_array			= ();
my @status_array		= ();
my @map_array			= ();
my @item_snp			= ();
my @snp_array			= ();
my @snp_array_1			= ();
my @chr_array_1			= ();
my @pos_array_1			= ();
my @genotype_array		= ();
my @selected_snp_array		= ();

my @row1 			= ();
my @row2 			= ();
my @row3 			= ();
my @row4 			= ();
my @array_2d 			= (\@row1, \@row2); # array to hold chisq reduction for SNP, sample

my @initial_chisq_array		= ();


my %specified			= ();
my %opts			= ();


print "\n";
print "#######################################################\n";
print "#     Which sample creates peak                       #\n";
print "#                                                     #\n";
print "#     This program removes samples one by             #\n";
print "#     one from the PLINK association analysis,        #\n";
print "#     and records the effect this has on the          #\n";
print "#     chisq value of each SNP.                        #\n";
print "#                                                     #\n";  
print "#     The results are saved in a tab-delimited file   #\n";
print "#     with a row per SNP and a column per sample.     #\n";
print "#                                                     #\n";
print "#     The data in the table is the reduction in       #\n";
print "#     chisq, so represents the effect each sample     #\n";
print "#     has on that particular SNP.                     #\n";
print "#                                                     #\n";
print "#######################################################\n";
print "\n";

print "Version $version\n\n";

#########################################
# process prefix and snp flags (if any) #
#########################################

GetOptions ("prefix=s" => \$prefix);


#######################################
# Get the prefix for the file names   #
#######################################

if ($prefix eq "")
{
	print "Prefix for PLINK .map and .ped files:  ";
	$prefix = <STDIN>;
	chomp $prefix;
	print"\n";
}


$too_many_samples = "False";

################################################
# If too many samples, get some information on #
# which selected SNPs are required.            #
# This is a list in a file                     #
################################################
if ($ans eq "y" || $ans eq "Y")
{
	$too_many_samples = "True";

	print "File containing a list of SNPs to check (less then 250):  ";
	$snpfile = <STDIN>;
	chomp $snpfile;
	print"\n";

	open (SNP, $snpfile)|| die "Cannot open input file: $snpfile";

	$no_selected_SNPs = 0;

	while ($single_line = <SNP>) 
	{
		$no_selected_SNPs = $no_selected_SNPs + 1;
		chomp $single_line;
		&chomp_all ($single_line);

		$selected_snp_array[$no_selected_SNPs]=$single_line;

	}
	close SNP;

	print "\nList of selected SNPs\n\n";
	for ($count=1;$count<=$no_selected_SNPs;$count++)
	{
		print "$count\t$selected_snp_array[$count]\n";
	}

	print "\n$no_selected_SNPs SNPs found in file\n\n";
	print "Press return to continue   ";
	$ans = <STDIN>;
}

##

#############
# Filenames #
#############

$pedfile = $prefix.".ped";
$mapfile = $prefix.".map";



$outfile="wscp_temp";

######################################
# Open the log file and results file #
######################################
$logfile = "wscp_".$prefix.".log";
$resultsfile = "wscp_results_".$prefix.".txt";
$remove_file = "wscp_remove.txt";
		
open (LOG, ">$logfile")|| die "Cannot open log file: $logfile";
print LOG "Which sample creates peak LOG FILE\n";
print LOG scalar localtime;
print LOG "\n";

open (RESULTS, ">$resultsfile")|| die "Cannot open results file: $resultsfile";




#######################################################
# Read in all SNPs in map file                        #
# and put them in snp_array                           #
# (This is so we can find which column the SNP is in) #
#######################################################


open (MAP, $mapfile)|| die "Cannot open input file: $mapfile";

print "\nReading MAP file...\n\n";
print "showing first 5 SNPs\n\n";

$count =0;
while ($single_line = <MAP>) 
{
	chomp $single_line;
	$count = $count + 1;
	@item_snp=split(/\t/,$single_line);
	$snp_name= $item_snp[1];
	chomp $snp_name;
	&chomp_all ($snp_name);

	if ($count <= 5)
	{
		print "$snp_name\n";
	}

	$snp_array[$count]=$snp_name;

}
close MAP;


$total_snps = $count;

print "$total_snps SNPs read from map file.\n\n";

print "Press return to continue  ";
$ans = <STDIN>;
print "\n";

#############################################
# Get the list of samples from the ped file #
# Also get the genotype for each sample     #
#############################################


open (PED, $pedfile)|| die "Cannot open input file: $pedfile";

while ($single_line = <PED>) 
{

	$linecount = $linecount + 1;
	# Remove the end of line character #
	chomp $single_line;

	#############################################
        # Split line at tabs into the array 'item'  #
	#############################################
        @item=split(/\t/,$single_line);

	$pedigree = $item[0];
	$sample = $item[1];
	$status = $item[5];

	$sample_array[$linecount] = $pedigree."\t".$sample;
	$aht_no_array[$linecount] = $sample;
	$status_array[$linecount] = $status;

	print "$linecount:   $sample_array[$linecount]\n";



} # end of while loop for reading PED file



$sample_array[0] = "9999\t9999";

$total_samples = $linecount;


	for ($sample_count = 0;$sample_count <=$total_samples;$sample_count++)
	{


			########################################
			# Create the temporary remove.txt file #
			########################################
	
			open (REMOVE, ">$remove_file")|| die "Cannot open remove file: $remove_file";
			print LOG "Samples to be removed: \n";
	
			# Add one line for each sample to be removed #
			# First the previous best_removes
			for ($count=1;$count < $no_to_remove;$count++)
			{
				print REMOVE "$sample_array[$best_remove[$count]]\n";
				print LOG "$count:   $sample_array[$best_remove[$count]]\n";
			}
	
			# Then add the sample which is being checked on this run
			print REMOVE "$sample_array[$sample_count]";
			print LOG "$no_to_remove:   $sample_array[$sample_count]\n\n";
	
	
			#################################################################
			#                          Now run PLINK                        #
			#################################################################
	
	
			$command = "plink --noweb --dog --allow-no-sex --file $prefix --assoc --remove $remove_file --out $outfile";
			print("$command\n");
			system("$command");
	
	

			#################################################################
			# Find chisq for all the SNPs and store in the 2d array         #
			#################################################################
			print "Finding out CHISQ value for all SNPs in temp.assoc\n";

			open (ASSOC, $outfile.".assoc") || die "Cannot open ASSOC file: $outfile.assoc";


			#########################################
			# Work down the lines in the assoc file #
			# recording chisq for each SNP          #
			#########################################
			$snp_count = 0;

			while ($single_line = <ASSOC>)
			{

				# Remove the end of line character #
				chomp $single_line;


				#########################################################
        			# Split line at TABS OR SPACES into the array 'item'    #
				#########################################################
        			@item=split(/\s+/,$single_line);

				$chromosome = $item[1];
				$snp_name = $item[2];
				$position = $item[3];
				$chisq_for_this_snp = $item[8];


				###############################################################
				# Store initial chisq for each SNP in the initial_chisq array #
				###############################################################
			
				if ($sample_count == 0) 
				{
					$initial_chisq_array[$snp_count] = $chisq_for_this_snp;
				}


				###############################################################
				# Store name of SNP in the snp_array_1                        #
				# and chromosome in chr_array_1                               #
				###############################################################
			
				if ($sample_count == 0) 
				{
					$initial_chisq_array[$snp_count] = $chisq_for_this_snp;
					$snp_array_1[$snp_count] = $snp_name;
					$chr_array_1[$snp_count] = $chromosome;
					$pos_array_1[$snp_count] = $position;
				}




				$snp_array[$snp_count] = $snp_name;

				#####################################
				# Calculate amount chisq is reduced #
				# and Store chisq in 2D array       #
				#####################################

				$chisq_reduction = $initial_chisq_array[$snp_count] - $chisq_for_this_snp;
				$array_2d[$snp_count][$sample_count] = $chisq_reduction;

				#if ($sample_count > 0)
				#{
				#if ($snp_count < 5)
				#{
				#	print "Sample count: $sample_count    SNP count: $snp_count\n";
				#	print "Initial_chisq: $initial_chisq_array[$snp_count]    Chisq: $chisq_for_this_snp   Reduction: $chisq_reduction\n";
				#	print "Array_2d: $array_2d[$snp_count][$sample_count]\n";
				#	$ans = <STDIN>;
				#}
				#}

				$snp_count = $snp_count + 1;

			} # end of while loop

			close ASSOC;

			print  "Sample: $sample_array[$sample_count]\tChisq: $chisq_for_this_snp\tStatus: $status_array[$sample_count]\n";

	
			print "\n";
			print "**************************************************************************************\n";
			print "Sample count:   \t$sample_count/$total_samples\n";
			print "Sample removed: \t$sample_array[$sample_count]\n";
			print "**************************************************************************************\n";


		

	} 
	# end of PLINK loop


#############################
# Write results in RESULTS  #
#############################

print "Total number of SNPs:     \t$total_snps\n";
print "Total number of samples:  \t$total_samples\n\n";


######################
# Write out all data #
######################

######################################################################
# If there are less than 256 samples, then we want to write the file #
# as SNPs down the side, and a samples across the top.               #
######################################################################
if ($too_many_samples eq "False")
{

	print "Too many samples is False\n";
	############################################
	# Write out sample names as column headers #
	############################################
	print RESULTS "CHR\tSNP\tPOS\t";

	for ($sample_count=1;$sample_count <= $total_samples;$sample_count++)
	{
		print RESULTS "$aht_no_array[$sample_count]\t";
	}

	print RESULTS "\n";


	for ($snp_count=1;$snp_count <= $total_snps; $snp_count++)
	{

		print RESULTS "$chr_array_1[$snp_count]\t";
		print RESULTS "$snp_array_1[$snp_count]\t";
		print RESULTS "$pos_array_1[$snp_count]\t";

		for ($sample_count=1;$sample_count <= $total_samples;$sample_count++)
		{

			print RESULTS "$array_2d[$snp_count][$sample_count]\t";

		}

		print RESULTS "\n";

	}
}


######################################################################
# If there are more than 256 samples, then we want to write the file #
# as samples down the side, and a selection of SNPs across the top   #
######################################################################
if ($too_many_samples eq "True")
{
		print "Too many samples is True\n";
	print RESULTS "Sample\t";

	#############################################
	# Write selected SNPs in the column headers #
	#############################################
	for ($snp_count=1;$snp_count <= $total_snps; $snp_count++)
	{
			# Check if SNP is one of the selected #
			for ($count=1;$count<=$no_selected_SNPs;$count++)
			{
				if ($snp_array[$snp_count] eq $selected_snp_array[$count])
				{
					print RESULTS "$snp_array[$snp_count]\t";
				}
			}

	}

	print RESULTS "\n";


	for ($sample_count=1;$sample_count <= $total_samples;$sample_count++)
	{

		#####################################
		# write sample name in first column #
		#####################################
		print RESULTS "$aht_no_array[$sample_count]\t";


		# Write selected SNPs in the columns #
		for ($snp_count=1;$snp_count <= $total_snps; $snp_count++)
		{

			# Check if SNP is one of the selected #
			for ($count=1;$count<=$no_selected_SNPs;$count++)
			{

				if ($snp_array[$snp_count] eq $selected_snp_array[$count])
				{
					print RESULTS "$array_2d[$snp_count][$sample_count]\t";
				}
			}


		}

		print RESULTS "\n";

	}
}




print RESULTS "\nFinished at ";
print RESULTS scalar localtime;
print "\n\nThe results are in $resultsfile\n\n";

close LOG;
close RESULTS;
close ASSOC;


exit;

##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

foreach (@_) {s/\n//g}  
foreach (@_) {s/\r//g}  

}


