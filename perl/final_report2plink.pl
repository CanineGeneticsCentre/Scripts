#!/usr/local/bin/perl

################################################################
#     Bead 2 PLINK v1                                          #
                                                               #
#                                                              #
#     Converts files from TALL Genome Studio format            #
#     to the format for PLINK                                  #
#                                                              #
################################################################


# Mike Boursnell June 2009

use strict;
use Getopt::Long;
use Term::ANSIColor;

my $version				= "4";

my $base_format			= "ACGT";


#####################################
# Define variables                  #
#####################################
my @item        		= ();
my @item_list			= ();
my @piece				= ();
my @names     			= ();
my @newnames			= ();
my @time_array			= ();
my @newline				= ();
my @firstline			= ();
my @standard_allele_A	= ();
my @standard_allele_B	= ();
my @major_allele		= ();
my @ID_array			= ();
my @SNP_array			= ();
my @row1 				= ();
my @row2 				= ();
my @row3 				= ();
my @row4 				= ();
my @array_2d 			= (\@row1, \@row2);
my @array_2d_recoded	= (\@row3, \@row4);

my %allele_A_hash		= ();
my %allele_B_hash		= ();

#Boolean
my $show				= "false"; # for debugging
my $add_column_headers	= "false"; # for debugging


#File names
my $infile       		= "";
my $listfile       		= "";
my $outfile      		= "";
my $pedfile				= "";
my $prefix				= "";
my $logfile				= "";

#Variables
my $single_line			= "";
my $list_ID				= "";
my $sample_ID			= "";
my $date_string			= "";
my $SNP_name			= "";
my $allele				= "";
my $allele_1			= "";
my $allele_2			= "";
my $allele_A			= "";
my $allele_B			= "";
my $allele_1_recoded 	= "";
my $allele_2_recoded 	= "";
my $allele_major		= "";
my $allele_temp			= "";
my $allele_string		= "";
my $genotype			= "";
my $genotype_recoded	= "";
my $genotype_string		= "";
my $new_ID				= "";
my $ans					= "";
my $first_SNP			= "";
my $mystring			= "";

my $allele_count				= 0;
my $header_line_array_size		= 0;
my $count						= 0;
my $column						= 0;
my $chosen_column				= 0;
my $help         				= 0;
my $line_count    				= 0;
my $listcount					= 0;
my $colcount     				= 0;
my $showcount    				= 0;
my $itemcount    				= 0;
my $checkcount   				= 0;
my $snp_count					= 1;
my $snp_out						= 0;
my $renamed_count				= 0;
my $total_lines  				= 0;
my $total_samples				= 0;
my $total_items  				= 0;
my $total_SNPs					= 0;
my $total_in_list				= 0;
my $year						= 0;
my $month						= 0;
my $day							= 0;
my $hour						= 0;
my $min							= 0;
my $sec							= 0;
my $found_on_list				= 0;
my $names_not_used_count		= 0;
my $copy_if_no_renaming			= 0;
my $first_time					= 1;
my $sample_count				= 0;
my $total_data_items 			= 0;
my $header_count				= 0;
my $genotype_count 				= 0;
my $A_count						= 0;
my $B_count						= 0;
my $valid_genotype_count 		= 0;
my $incorrect_genotype_count 	= 0;
my $snp_valid_genotype_count	= 0;
my $swap_count					= 0;
my $no_genotype_count 			= 0;



print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      final_report2plink             \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This PERL script converts GenomeStudio genotype data files\n\n";
print "    from the FinalReport format to PLINK format\n\n";
print "    It doesn't currently make a map file.\n\n";
print color 'reset';


#############################
# Get the file names        #
#############################
until (-e $infile)
{
	print "Name of input FinalReport file:    ";
	$infile = <STDIN>;
	chomp $infile;
	if ($infile eq ""){$infile = "FinalReport.txt"}  # TEMP!!!
	if ($infile eq "ls"){print "\n";system ("ls *inal*");print "\n"}
	if ($infile ne "ls"){if (! -e $infile){print "\n\n>>>>>>>>  File $infile not found.  Try again.  <<<<<<<<\n\n";}}
}

open (IN, "$infile") || die "Cannot open input file: $!";

print "Prefix for PLINK format output file: ";
$outfile = <STDIN>;
chomp $outfile;

$prefix = &get_prefix ($outfile);

$pedfile = $outfile.".ped";

print"\n";

#####################################
# Open the chosen file for output   #
#####################################
open (OUT, ">$pedfile")|| die "Cannot create output file: $pedfile";


###################################################
# Read data in INPUT file                         #
# This reads the input file so that the genotypes #
# are stored in a 2-dimensional array called      #
# $array_2d[$sample_count][$snp_count]            #
###################################################
print "Reading input file $infile...\n\n";

&print_message("THESE ARE THE 8 HEADER LINES","message");

###########################################
# Work down the lines in the input file   #
###########################################
while ($single_line = <IN>) {

	# Remove the end of line character #
	chomp $single_line;

	if ($line_count < 9){print "$single_line\n";} # Show first set of header lines

	###################################################################################################
	# This is the line with the column headings, so use this to ask the user which columns they want  #
	###################################################################################################
	if ($line_count == 9)
	{

		&print_message("THIS IS THE COLUMN HEADINGS LINE","message");

		@item=split(/\t/,$single_line);

		$header_line_array_size = scalar @item;

		print "Column 1\t\t$item[0]\t\tThis should be the 'SNP Name'\n";
		print "Column 2\t\t$item[1]\t\tThis should be the 'Sample ID'\n";

		for ($count = 2; $count < $header_line_array_size; $count ++)
		{
			$column = $count + 1;
			print "Column $column:\t\t$item[$count]\n";
		}
		
		print "\n";
		&print_message("Which allele columns do you want (type in column number of first of the pair): ","input");
		$ans = <STDIN>;

		if ($ans > 2)
		{
			$chosen_column = $ans;
			print "\nColumns chosen:  $item[$chosen_column-1]\t$item[$chosen_column]\n\n";

			print "If these are the correct allele columns press RETURN";
			$ans = <STDIN>;

			print "\n\nReading in data from $infile...\n\n";
		}

	}
	if ($line_count > 9)
	{
		

		################################################################
	    # Split line at tabes or spaces into the array 'item'          #
		################################################################
        @item=split(/\t/,$single_line);

		$SNP_name  = $item[0];
		$sample_ID = $item[1];
		$allele_1  = $item[$chosen_column - 1];
		$allele_2  = $item[$chosen_column];
		$genotype  = "$allele_1"."$allele_2";

		#$genotype = &chomp_all ($genotype);
		$genotype =~ s/\s+$//; # get rid of \r at end of line

		##################
		# If no genotype #
		##################
		if (($genotype eq "--") || ($genotype eq "  ") || ($genotype eq "00") || ($genotype eq "NN"))
		{
			$genotype = "XX";
			$no_genotype_count = $no_genotype_count + 1;
		}
		else
		{
			$valid_genotype_count = $valid_genotype_count + 1;
		}
		

		######################################################
		# Note name of first SNP (first SNP line is line 10) #
		######################################################
		if ($line_count == 10)
		{
			$first_SNP = $SNP_name;
			$sample_count = $sample_count + 1;
			$ID_array[$sample_count]= $sample_ID;
			$SNP_array[$snp_count] = $SNP_name;
		}


		###############################################################
		# Check each SNP name to see when first one comes round again #
		###############################################################

		if (($SNP_name eq $first_SNP) && ($line_count > 10))
		{
			print "Line: $line_count \t\tSample: $sample_ID\n";

			$sample_count = $sample_count + 1;
			$ID_array[$sample_count]= $sample_ID;
			$SNP_array[$snp_count] = $SNP_name;

			$first_time=0;
			$snp_count = 1;
		}


		##################################
		# Store genotype in the 2D array #
		##################################
		$genotype =~ s/\s+$//; # get rid of \r at end of line


		##################################
		# PLINK uses zero as 'no allele' #
		##################################
		if ($genotype eq "XX"){$genotype = "00"}

		$array_2d[$sample_count][$snp_count] = $genotype;
		$SNP_array[$snp_count] = $SNP_name;
		
		# increase SNP count by one

		$snp_count = $snp_count + 1;

	} # end of if line_count > 9

	$line_count = $line_count + 1;

} # end of while loop

$total_lines=$line_count;
$total_samples = $sample_count;

$snp_count = $snp_count - 1;

$total_SNPs = $snp_count;
print "\nNumber of SNPs: $total_SNPs\n";


##################################################################################################
# Now that all the file has been read we can loop through the SNPs and get all the data we need  #
##################################################################################################

print "\nLooping through all the SNPs...\n\n";

for ($snp_count = 1;$snp_count <= $total_SNPs;$snp_count++)
{
	if ($snp_count % 1000 == 0){print "SNP: $snp_count\n"}

	###########################################################
	# First find which 2 alleles are present at this position #
	###########################################################
	$allele_A="S";   $allele_B="T";

	for ($sample_count=1; $sample_count <= $total_samples;$sample_count++)
	{
		$genotype = $array_2d[$sample_count][$snp_count];
		$allele = "";

		if ($genotype ne "XX")
		{
			for ($allele_count=1;$allele_count<=2;$allele_count++)
			{
				$allele = substr($genotype,$allele_count -1,1);

				if (not defined $allele_A_hash{$snp_count}){$allele_A_hash{$snp_count} = $allele;}
				else
				{
					if (not defined $allele_B_hash{$snp_count})
					{
						if ($allele ne $allele_A_hash{$snp_count}){$allele_B_hash{$snp_count} = $allele;}
					}
				}
			} # allele_count
		}# if genotype not XX

	} # sample_count loop 1 for finding allele_A and allele_B

	####################################################
	# Check that two standard alleles aren't the same  #
	####################################################
	if ((defined $allele_A_hash{$snp_count}) && (defined $allele_A_hash{$snp_count}))
	{
		if ($allele_A_hash{$snp_count} eq $allele_B_hash{$snp_count})
		{
			print "Sample: $sample_count SNP: $snp_count allele_A_hash shouldn't be the same as allele_B_hash\n\n";
			$ans=<STDIN>;
		}
	}

	###############################################################################
	# Now store the standard alleles in allele_A and allele_B (for this SNP only) #
	###############################################################################
	$allele_A = $allele_A_hash{$snp_count};
	$allele_B = $allele_B_hash{$snp_count};


	###################################################################
	# Find out which allele is the major allele (sample_count loop 2) #
	###################################################################
	$A_count = 0;
	$B_count = 0;

	for ($sample_count=1; $sample_count <= $total_samples;$sample_count++)
	{
		$genotype = $array_2d[$sample_count][$snp_count];  # This could still be "XX"

		$allele_1 = substr($genotype,0,1);
		$allele_2 = substr($genotype,1,1);

		if ($allele_1 eq $allele_A) {$A_count = $A_count + 1}
		if ($allele_1 eq $allele_B) {$B_count = $B_count + 1}
		if ($allele_2 eq $allele_A) {$A_count = $A_count + 1}
		if ($allele_2 eq $allele_B) {$B_count = $B_count + 1}

	} # sample_count loop 2 for finding major allele

	##########################################################
	# Now store the MAJOR alleles in the array major_alleles #
	##########################################################
	if ($A_count >= $B_count)
	{$allele_major = $allele_A}
	else
	{$allele_major = $allele_B}


	##################################################
	# Make sure allele_B represents the major allele #
	# If not, swap allele_A and allele_B             #
	# (Do we want this in?)                          #
	##################################################
	if ($allele_B ne $allele_major)
	{
		$allele_temp = $allele_A;
		$allele_A = $allele_B;
		$allele_B = $allele_temp;
		$swap_count = $swap_count + 1;
	}


	#################################################################################
	# Now we have found alleleA and allele_B, and counted which is the major allele #
	# we need to put that swapped genotype into a new 2d array                      #
	#################################################################################

	for ($sample_count=1; $sample_count <= $total_samples;$sample_count++)
	{
		$genotype = $array_2d[$sample_count][$snp_count];
		$allele_1 = substr($genotype,0,1);
		$allele_2 = substr($genotype,1,1);

		$array_2d_recoded[$sample_count][$snp_count] = $genotype; # add later!!!!!
	}
} # snp_count loop



############################
# Write from array to file #
############################

print "\n\nWriting output file $pedfile ...\n\n";


if ($add_column_headers eq "true")
{
	#######################################################
	# Write the first six column headings of the PED file #
	#######################################################

	print OUT "FID\tIID\tFather\tMother\tSex\tStatus";


	#######################################
	# Write the column headings (SNPs)    #
	#######################################

	for ($snp_count=1;$snp_count<=$total_SNPs;$snp_count++)
	{
		print OUT "\t$SNP_array[$snp_count]";
		print OUT "\t$SNP_array[$snp_count]";
	}
	print OUT "\n";

} # add headers if true

#######################################
# Write the rest of the sample rows   #
#######################################
for ($sample_count=1;$sample_count<=$total_samples;$sample_count++)
{
	print OUT "$ID_array[$sample_count]\t$ID_array[$sample_count]\t0\t0\t0\t0";

	for ($snp_count=1;$snp_count<=$total_SNPs;$snp_count++)
	{
		$genotype = $array_2d[$sample_count][$snp_count];
		$allele_1 = substr($genotype,0,1);
		$allele_2 = substr($genotype,1,1);


		print OUT "\t$allele_1\t$allele_2";

	}

	print OUT "\n";
}

###########################
#     Close all files     #
###########################
close IN;
close OUT;

&print_message("Finished converting $infile to $pedfile","message");

print "Total number of SNPs found:                          \t$total_SNPs\n";
print "Total number of samples counted:                     \t$total_samples\n\n";
print "Valid genotype count:                                \t$valid_genotype_count\n";
print "Failed genotype count:                               \t$no_genotype_count\n";
print "Number of genotypes swapped so major allele was A:   \t$swap_count\n";

print"\n";
print "Input FinalReport file:                              \t$infile\n";
print "Output PED file:                                     \t$pedfile\n\n";


exit;





##########################################################
# Gets prefix of file. e.g. getes 'name' from 'name.txt' #
##########################################################

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
	if ($style eq "help"){$char = "+"}
	
	print "\n\n";
	print color ' bold yellow';
	if ($style eq "warning"){print color ' bold red'}
	if ($style eq "input"){print color ' bold white'}
	if ($style eq "help"){print color ' bold green'}
	if ($style eq "message"){print color ' bold yellow'}
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n$char    $message    $char\n";
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n\n";
	print color 'reset';

}




##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

foreach (@_) {s/\n//g}  
foreach (@_) {s/\r//g}  

}