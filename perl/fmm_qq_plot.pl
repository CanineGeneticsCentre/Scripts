#!/usr/local/bin/perl 

################################################################
#     FMM QQ plot                                              #
#                                                              #
#    Makes a QQ plot from FMM output data                      #
#                                                              #
################################################################


# Mike Boursnell Dec 2010
# Animal Health Trust
# Newmarket
# UK
# mike.boursnell@aht.org.uk

use strict;
use Getopt::Long;
use Statistics::Distributions;  # Fot CHI-SQUARED
use Term::ANSIColor;

my $version				= "2";
#####################################
# Define variables                  #
#####################################

my $prefix			= "";
my $command			= "";
my $outfile			= "";
my $phenofile			= "";
my $fmm_out_file	= "";
my $qq_file			= "";
my $kinfile			= "";
my $kinfile_BN			= "";
my $kinfile_IBS			= "";
my $kintype			= "IBS";
my $pedfile			= "";
my $PDF_file			= "";
my $psfile			= "";
my $gnufile			= "";
my $pedigree			= "";
my $id				= "";
my $status			= "";
my $single_line			= "";
my $parameters_ok		= "false";
my $ans				= "";
my $pheno_string		= "";
my $species			= "";
my $tfamfile			= "";
my $tpedfile 			= "";
my $mapfile			= "";
my $species_code		= "";
my $chromosome			= "";
my $SNP_name			= "";
my $position			= "";
my $psfile			= "";
my $final_outfile		= "";
my $final_outfile_hd		= "";
my $last_chromosome		= "";
my $is_y_chromosome 		= "";
my $lambda_string		= "";
my $title			= "";

my $geno_value			= 0.1;
my $maf_value			= 0.05;
my $mind_value			= 0.1;
my $line_count			= 0;
my $beta			= 0;
my $P_value			= 0;
my $minus_log_P			= 0;
my $snp_count_tped		= 0;
my $snp_count_ps		= 0;
my $rows_per_chromosome_count	= 0;
my $max_row_count		= 0;
my $array_count			= 0;
my $chromosome_count		= 0;
my $row_count			= 0;
my $max_chromosome_count	= 0;
my $chi_sq			= 0;
my $expected_chi_sq		= 0;
my $rank			= 0;
my $expected_P_value		= 0;
my $total_lines			= 0;
my $median			= 0;
my $lambda			= 0;
my $array_size			= 0;

my @item         		= ();
my @chr_array			= ();
my @snp_array			= ();
my @pos_array			= ();
my @rows_per_chromosome_array	= ();
my @new_line_array		= ();
my @minus_log_P_array		= ();
my @P_value_array		= ();
my @chi_sq_array		= ();
my @chi_sq_array_reverse	= ();
my @chi_sq_array_expected	= ();


###############################
# Process flags               #
###############################

GetOptions("file=s"=>\$fmm_out_file,"out=s"=>\$outfile);


print "\n";
print "###################################################################\n";
print "#     FMM QQ plot v1                                              #\n";
print "#                                                                 #\n";
print "#     Makes a QQ plot from the FMM output file                    #\n";
print "#                                                                 #\n";
print "###################################################################\n\n";


#############################
# Get the input file name   #
#############################

		
if ($fmm_out_file eq "")
{
	until (-e $fmm_out_file)
	{
		print "Name of FMM output file:  ";
		$fmm_out_file = <STDIN>;
		chomp $fmm_out_file;
		
		if ($fmm_out_file eq "ls"){print "\n";system ("ls *fmm.out")}
		
		if ($fmm_out_file ne "ls")
		{
			if (! -e $fmm_out_file){print "\n\n>>>>>>>>  File $fmm_out_file not found.  Try again.  <<<<<<<<\n\n";}
		}
			
		print"\n";
	}
}



######################
# Make up file names #
######################
$qq_file = $fmm_out_file."_fmm_qq";


open (FMM, "$fmm_out_file") || die "Cannot open FMM output file $fmm_out_file";


####################################################
# Open the FMM out put file and get the CHI values #
####################################################
$line_count = 0;

while ($single_line = <FMM>) 
{

	chomp $single_line;
	$line_count = $line_count + 1;
	@item=split(/\t/,$single_line);
	$chi_sq=$item[4];
	$chi_sq_array[$line_count] = $chi_sq;

}

$total_lines = $line_count;

print "$line_count lines read.\n\n";

###########################
# Make output for QQ plot #
###########################

&print_message("Calculating data for QQ plot");


##################
# Sort the array #
##################

&print_message("Sorting file");

@chi_sq_array_reverse = reverse sort {$a <=> $b} (@chi_sq_array);


#################################
# Print first 10 lines of array #
#################################

&print_message("Showing the first 10 lines of the sorted array");

for($line_count = 0; $line_count <=10; $line_count++)
{
	print "$line_count\t$chi_sq_array_reverse[$line_count]\n";
}

print "\n";



##################
#Stats reminders #
##################
#$expected_chi_sq = Statistics::Distributions::chisqrdistr(1,$expected_P_value);
#$expected_P_value = Statistics::Distributions::chisqrprob(1,$expected_chi_sq);



##############################################
# Now calculate Exp_p with this formula:     #
#                                            #
# Exp-p = (rank - 0.5)/total_no_of_lines     #
# and then getting the CHISQ of that P-value #
##############################################

for($array_count = 0; $array_count <= $total_lines; $array_count++)
{
	$rank = $array_count + 1;
	$expected_P_value = ($rank - 0.5) / $total_lines;

	if ($expected_P_value <= 1)
	{
		$expected_chi_sq = Statistics::Distributions::chisqrdistr(1,$expected_P_value);
		$chi_sq_array_expected[$array_count] = $expected_chi_sq;
	}

	
}


###############################################
# Calculate the median of the observed values #
###############################################
$array_size = @chi_sq_array_reverse;
$median = $chi_sq_array_reverse[int($array_size/2)];
$lambda = $median / 0.456;

$lambda_string = sprintf("%.3f",$lambda);

print "Median: \t$median\n";
print "Lambda: \t$lambda_string\n\n\n"; 
 


############################
# Open files for output    #
############################

open (OUT, ">$qq_file") || die "Cannot open $qq_file";

print OUT "FMM\tBLANK\tEXP_CHI\tOBS_CHI\n";

for ($line_count = 0; $line_count <= $total_lines; $line_count++)
{
	if (($chi_sq_array_reverse[$line_count] * 1) > 0)
	{
		print OUT "$chi_sq_array_reverse[$line_count]\t\t$chi_sq_array_expected[$line_count]\t$chi_sq_array_reverse[$line_count]\n";
	}
}

print OUT "\t\t\t\t\n";
print OUT "\t\t\t$lambda_string";

close OUT;


###############################################
# Make GNUplot file for pdf figure of QQ plot #
###############################################

$gnufile = "gnufile.xtxt";
$psfile = $fmm_out_file."_fmm_qq.ps";
$PDF_file = $fmm_out_file."_fmm_qq.pdf";
$title = "$fmm_out_file:  QQ plot of FMM -corrected chi-squared values";

print "Default title:  $fmm_out_file:  QQ plot of FMM-corrected chi-squared values\n\n";
print "\nType in the title for this graph: ";
print "(Press return to accept the default title)\n\n";
print "> ";

$ans=<STDIN>;
chomp $ans;

if ($ans ne ""){$title = $ans}



open (GG, ">$gnufile") || die "can't open $gnufile\n" ;

##################################################
# Write the gnuplot commands to the command file #
##################################################
print GG "set terminal postscript color solid\n";
print GG "set title  \"$title\" \n" ; 
print GG "set key outside\n"; 
print GG "set xlabel  \"Expected\" \n" ; 
print GG "set ylabel  \"Observed\" \n" ; 
print GG "set nokey\n";

print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 3\n";
print GG "set style line 2 lt -1 lw 1 linecolor 2\n";

print GG "set label 1 'Inflation factor: $lambda_string' at graph 0.5,0.1\n";

print GG "plot '$qq_file' using 3:3 with lines ls 2 title 'NULL', '$qq_file' using 2:3 with points ls 1 title 'QQ'" ;


###############################
# Run the commands from linux #
###############################
&print_message("Running GNUplot to create pdf file");

system "gnuplot < $gnufile > $psfile" ;
#system "fixgreen  $psfile" ;
system "ps2pdf  $psfile " ;


print color 'yellow';
print "######################################################\n";
print "#               FMM QQ plot completed                #\n";
print "######################################################\n\n";


print "Inflation factor:    \t$lambda_string\n\n";

print "PDF file of QQ-plot: \t$PDF_file\n\n";

print color 'reset';

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

######################################
# Subroutine to print screen message #
######################################

sub print_message
{
	my $message 		= "";
	my $message_length 	= "";
	my $pos_count		= 0;
	
	$message = $_[0];
	$message_length = length($message);
	
	print "\n\n";
	print color 'yellow';
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++)
	{
		print "#";
	}
	print "\n#    $message    #\n";
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++)
	{
		print "#";
	}
	print "\n\n";
	
	print color 'reset';

}


