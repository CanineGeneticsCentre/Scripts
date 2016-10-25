#!/usr/local/bin/perl

################################################################
#     pdf_summary                                              #
#                                                              #
#     Makes a summary pdf of all the data in the folder        #
#                                                              #
################################################################

# Mike Boursnell Feb 2009
# Animal Health Trust
# Newmarket
# UK

use strict;
use Getopt::Long;
use Statistics::Distributions;  # For CHI-SQUARED
use Term::ANSIColor;
use Cwd; # Current working directory

my $version					= "13";

#####################################
# Define variables                  #
#####################################

#Constants
my $x_axis_font_size		= "8";
my $label_offset_factor		= "10";
my $data_point_offset_y		= "0"; # Adds this to all data points (if positive moves them up the graph)
my $size_of_dot				= 300;

#File names
my $prefix					= "";
my $prefix_default			= "";
my $prefix_temp				= "";
my $gnufile					= "";
my $infile					= "";
my $psfile					= "";
my $PDF_file				= "";
my $pedfile					= "";
my $mds_tempfile			= "";
my $pedfile_up_one_directory = "";
my $samples_used_file		= "";
my $qq_file 				= "";
my $file_to_find			= "";

# data files used
my $fmm_out_file			= "";
my $fmm_temp_file			= ""; # Used to make fmm_out file look like emmax_out file

# Variables
my $bottom_axis_threshold	= 0.05;
my $chr_count				= 0;
my $total_no_chr			= 0;
my $total_no_lines			= 0;
my $line_count				= 0;
my $check_count				= 0;
my $x_label_position		= 0;
my $label_offset_x			= 0;
my $label_offset_y			= 0;
my $total_no_snps			= 0;
my $p_value					= 0;
my $max_p_value				= 0;
my $max_minus_log_p			= 0;
my $total_no_controls		= 0;
my $total_no_cases			= 0;
my $array_size				= 0;
my $median					= 0;
my $lambda					= 0;
my $total_unadj				= 0;
my $array_mean				= 0;
my $array_count				= 0;
my $chi_sq					= 0;
my $minus_log_p				= 0;
my $top_unadj				= 0;
my $P_value					= 0;
my $sum_chi_sq				= 0;
my $xpos					= 0;
my $ypos					= 0;
my $yinc					= 0;
my $file_type_count			= 0;
my $total_no_file_types		= 0;
my $expected_P_value		= 0;
my $total_lines_fmm_out				= 0; # can merge with total_no_lines??
my $rank					= 0;
my $expected_chi_sq			= 0;
my $file_count				= 0;

my $answer					= "";
my $snp_name				= "";
my $position				= "";
my $title					= "";
my $single_line				= "";
my $chromosome				= "";
my $last_chromosome			= "";
my $ans						= "";
my $pedigree				= "";
my $id						= "";
my $C1						= "";
my $C2						= "";
my $affection				= "";
my $qq						= "";
my $lambda_string			= "";
my $unadj					= "";
my $log_values				= "";
my $date					= "";
my $plot_choice				= "";
my $current_directory		= "";

# Boolean
my $some_files_found 		= "false";
my $use_pedfile				= "no";
my $use_samples_file		= "yes";
my $show					= "true";
my $command_line_mode		= "false"; # if it is being run from within GWAS_analyses

my @item					= ();
my @chr_start_array			= ();
my @id_array				= ();
my @pedigree_array			= ();
my @affection_array			= ();
my @qq_array				= ();
my @chi_sq_array			= ();
my @chi_sq_array_reverse 	= ();
my @chi_sq_array_expected	= ();
my @p_value_array			= ();
my @minus_log_p_array		= ();
my @file_type_array			= ();
my @choice_array			= ();


############################
# process -f flag (if any) #
############################

GetOptions("f=s"=>\$prefix);
if ($prefix ne ""){$command_line_mode	= "true"} else {$command_line_mode	= "false"}

print color 'bold magenta';



print "\n\n";
print "############################\n";
print color 'bold white';
print "      PDF summary       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";

print "  - This program creates a PDF summary document with plots from all the \n";

print "    GWAS analysis files in the current directory..\n\n";

print "    Or it can plot individual PDFs.\n\n";

print " It looks for the folllowing files:\n\n";

print "    example.assoc             \tPLINK association files\n";
print "    example_unfiltered.assoc  \tPLINK association files with no MAF, GENO, MIND filtering\n";
print "    example.assoc.mperm       \tPLINK association mperm files\n";
print "    example.assoc.adjusted    \tPLINK QQ data files\n";
print "    example.mds               \tPLINK MDS files\n\n";

print "    example.emmax_out         \tEMMAX association files\n";
print "    example.emmax_qq          \tEMMAX QQ plot data files\n\n";

print "    example_fmm.out           \tFMM association files\n\n";

print color 'reset';


########################################
# Find any ped files in this directory #
########################################
$current_directory=getcwd;
opendir (DIR, $current_directory);

while ($file_to_find = readdir(DIR)) 
{
	if ((substr($file_to_find, -4) eq ".ped") && !($file_to_find =~ m/^\./))
	{
		$file_count = $file_count + 1;
		$prefix_temp = $file_to_find;
	}
}

if ($file_count == 1)
{
	$prefix_default = &get_prefix($prefix_temp);
}

#############################
# Get the input file name   #
#############################

if ($prefix eq "")
{
	until ((-e $pedfile) || (-e $pedfile_up_one_directory))
	{
		if ($prefix_default ne ""){print "Prefix for PLINK .map and .ped files (default = $prefix_default):  ";}
		if ($prefix_default eq ""){print "Prefix for PLINK .map and .ped files:  ";}
		$prefix = <STDIN>;
		chomp $prefix;

		##########################################################
		# If there is only one file you can use return to get it #
		##########################################################
		if ($prefix eq ""){$prefix = $prefix_default}

		# Remove .ped if accidentally added #
		if (index($prefix,".ped") > -1)
		{
			$prefix = substr($prefix,0,length($prefix) - 4);
			print "PED REMOVED TO GIVE $prefix\n\n";
		}

		if ($prefix eq "ls"){print "\n";system ("ls *.ped");system ("ls ../*.ped")}

		if ($prefix ne "ls")
		{
			$pedfile = $prefix.".ped";
			$pedfile_up_one_directory = "../$pedfile";

			if ((! -e $pedfile) && (! -e $pedfile_up_one_directory)){print "\n\n>>>>>>>>  File $pedfile not found in this directory or the one above.  Try again.  <<<<<<<<\n\n";}
		}
		print"\n";
	}
}

chomp $prefix;

###########################
# Make up some file names #
###########################
$samples_used_file = $prefix.".samples_used.out";
$fmm_out_file = $prefix."_fmm.out";
$fmm_temp_file = $prefix."_fmm_temp.out";
$gnufile = "gnufile.txt";
$psfile = $prefix."_summary.ps";
$PDF_file = $prefix."_summary.pdf";
$qq_file = $fmm_out_file."_fmm_qq";


###########################################
# Warn user if PDF summary already exists #
###########################################

if ($command_line_mode eq "false")
{
	if (-e $PDF_file)
	{
		&print_message("An output file called $PDF_file already exists");

		print "Continuing will overwrite the file.  Type 'q' to quit if you want to rename or move this file.\n\n";
		print ">  ";

		$answer=<STDIN>;
		chomp $answer;
		if (lc $answer eq "q"){exit;}
	}
} # $command_line_mode = false


################################################
# MDS file will need PED file so look for this #
################################################
if (-e $prefix.".mds")
{ 
	
	######################################
	# Look for ped file one directory up #
	######################################
	$pedfile = $prefix.".ped";
	$pedfile_up_one_directory = "../$pedfile";

	if ((! -e $pedfile) && (! -e $samples_used_file))
	{
		&print_message("Affection files for MDS plot not found in this directory");
		$use_samples_file = "no";
		$use_pedfile = "no";
		
		print "\tYou need either a PED file or a samples_used file\n\n";
		print "\tLooking in the directory above this one...\n\n";
		
		if (-e $pedfile_up_one_directory)
		{
			print "\t\t$pedfile_up_one_directory PED file found one directory up.\n\n";
			
			print "\t\tDo you want to use this PED file? (y/n)  ";
			$ans =<STDIN>;
			chomp $ans;
			if (lc $ans eq "y"){$pedfile=$pedfile_up_one_directory;$use_pedfile="yes";$use_samples_file="no"}
			
		}
	}

}


####################################
# Make a list of file types found  #
####################################
$file_type_count = 0;

if (-e $prefix.".assoc")
{$file_type_count = $file_type_count + 1; $choice_array[$file_type_count] = "plink_assoc"; $file_type_array[$file_type_count] = "assoc for association Manhattan plot"}

if (-e $prefix."_unfiltered.assoc")
{$file_type_count = $file_type_count + 1; $choice_array[$file_type_count] = "plink_unfiltered_assoc"; $file_type_array[$file_type_count] = "assoc unfiltered for association Manhattan plot"}

if (-e $prefix.".assoc.mperm")
{$file_type_count = $file_type_count + 1; $choice_array[$file_type_count] = "plink_mperm_assoc";$file_type_array[$file_type_count] = "assoc.mperm for MPERM association plot"}

if (-e $prefix.".assoc.adjusted")
{$file_type_count = $file_type_count + 1; $choice_array[$file_type_count] = "plink_qq_plot";$file_type_array[$file_type_count] = "assoc.adjusted for QQ-plot"}

if (-e $prefix.".emmax_out")
{$file_type_count = $file_type_count + 1; $choice_array[$file_type_count] = "emmax_assoc";$file_type_array[$file_type_count] = "emmax_out for EMMAX association Manhattan plot"}

if (-e $prefix.".emmax_qq")
{$file_type_count = $file_type_count + 1; $choice_array[$file_type_count] = "emmax_qq_plot";$file_type_array[$file_type_count] = "emmax_qq for EMMAX QQ-plot"}

if (-e $prefix."_fmm.out")
{$file_type_count = $file_type_count + 1; $choice_array[$file_type_count] = "fmm_assoc";$file_type_array[$file_type_count] = "_fmm.out for Fast Mixed Model Manhattan plot"}


if (-e $prefix.".mds")
{
	if ((-e $pedfile)  && (-e $samples_used_file))
	{
		$file_type_count = $file_type_count + 1; $choice_array[$file_type_count] = "plink_mds_plot" ;$file_type_array[$file_type_count] = "mds for MDS plot";
	}
}

$total_no_file_types = $file_type_count;

print "\nThese files were found. Which of these do you want to plot?\n\n";

for ($file_type_count = 1; $file_type_count <= $total_no_file_types; $file_type_count++)
{
	print "<$file_type_count>  $file_type_array[$file_type_count]\n";
}

print "\n<a>  Plot all in a single summary PDF\n\n";

$answer=<STDIN>;chomp $answer;

if ((lc($answer) eq "a") || ($answer eq "")){$plot_choice = "all"}
else
{$plot_choice = $choice_array[$answer]}


######################################################
# Open gnufile to receive a list of gnuplot commands #
######################################################
open (GG, ">$gnufile") || die "can't open $gnufile\n" ;

print GG "set terminal postscript color solid\n";

##########################
# Draw the title page    #
##########################

&title_page;


###################################################
# Run all the subroutines, one for each file type #
###################################################

###########################
# PLINK association       #
###########################
if (($plot_choice eq "all") || ($plot_choice eq "plink_assoc"))
{
	if (-e $prefix.".assoc")
	{ 
		&plink_assoc;
		$some_files_found = "true";
	}
}


################################
# PLINK unfiltered association #
################################
if (($plot_choice eq "all") || ($plot_choice eq "plink_unfiltered_assoc"))
{
	if (-e $prefix."_unfiltered.assoc")
	{ 
		&plink_unfiltered_assoc;
		$some_files_found = "true";
	}
}



###########################
# PLINK mperm association #
###########################
if (($plot_choice eq "all") || ($plot_choice eq "plink_mperm_assoc"))
{
	if (-e $prefix.".assoc")
	{ 
		&plink_mperm_assoc;
		$some_files_found = "true";
	}
}

###########################
# PLINK MDS plot          #
###########################
if (($plot_choice eq "all") || ($plot_choice eq "plink_mds_plot"))
{
	if (-e $prefix.".mds")
	{ 
		&plink_mds_plot;
		$some_files_found = "true";
	}
}


###########################
# PLINK QQ plot           #
###########################
if (($plot_choice eq "all") || ($plot_choice eq "plink_qq_plot"))
{
	if (-e $prefix.".assoc.adjusted")
	{ 
		&plink_qq_plot;
		$some_files_found = "true";
	}
}


###########################
# EMMAX QQ plot           #
###########################
if (($plot_choice eq "all") || ($plot_choice eq "emmax_qq_plot"))
{
	if (-e $prefix.".emmax_qq")
	{ 
		&emmax_qq_plot;
		$some_files_found = "true";
	}
}


###########################
# EMMAX assoc plot        #
###########################
if (($plot_choice eq "all") || ($plot_choice eq "emmax_assoc"))
{
	if (-e $prefix.".emmax_out")
	{ 
		&emmax_assoc;
		$some_files_found = "true";
	}
}


###########################
# FMM assoc plot          #
###########################
if (($plot_choice eq "all") || ($plot_choice eq "fmm_assoc"))
{

	if (-e "$fmm_out_file")
	{ 
		&fmm_assoc;
		$some_files_found = "true";
	}
}


##############################
# Warn if no files found     #
##############################

if ($some_files_found eq "false")
{
	&print_message("No files with the prefix $prefix found!");
	exit;
}

##########################
# Draw the last page     #
##########################

&last_page;



###############################
# Run the commands from linux #
###############################

&print_message("Creating PDF file");

print "\tCreating postscript file...\n\n";

system "gnuplot < $gnufile > $psfile" ;
#system "fixgreen  $psfile" ;

print "\tConverting postscript to PDF...\n";

system "ps2pdf  $psfile " ;

close GG;

&print_message("FINISHED PROGRAM");

print "PDF summary file: $PDF_file\n\n";

exit;

########################################################################################################
########################################################################################################


sub plink_assoc
{
################################################################################
# PLINK assoc file subroutine                                                  #
################################################################################

&print_message("Creating PLINK association plot");

$infile = $prefix.".assoc";

$title = "$prefix:  PLINK association data";
$total_no_chr = 3;

###############################################
# plink_assoc: Open file to get some key data #
###############################################
$chr_count = 0;
$line_count = 0;
$max_p_value = 0;

$last_chromosome = "1";
open (IN, "$infile") || die "Cannot open $infile";
while ($single_line = <IN>)
{
	$line_count = $line_count + 1;
	chomp $single_line;

	@item=split(/\s+/,$single_line);

	$chromosome = $item[1];
	$p_value = $item[9];

	if ($chromosome ne $last_chromosome)
	{
		$chr_count = $chr_count + 1;
		$chr_start_array[$chr_count-1] = $line_count;
	}

	if ($p_value > $max_p_value){$max_p_value = $p_value}

	$last_chromosome = $chromosome;
}

close (IN);

$total_no_chr = $chr_count -1;
$total_no_snps = $line_count;

print "    Total no chrs: $total_no_chr\tTotal no SNPs: $total_no_snps\n\n";


###############################################################
# plink assoc: Write the gnuplot commands to the command file #
###############################################################

print GG "set title  \"$title\" \n"; 
print GG "set key outside\n"; 
print GG "set xlabel  \"Chromosome\" offset 0,-1.2\n"; 
print GG "set ylabel  \"-logP\" \n" ; 
print GG "set nokey\n";
print GG "unset xtics\n";
print GG "set ytics\n";


print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 1\n";
print GG "set style line 2 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 3\n";
print GG "set style line 3 lt -1 lw 1 linecolor 2\n";

###########################
# assoc     x-axis labels #
###########################

$label_offset_y = $max_p_value / $label_offset_factor;

for ($chr_count = 1;$chr_count <=$total_no_chr; $chr_count++)
{
	if ($chr_count < 10) {$label_offset_x = 35}
	if ($chr_count >9) {$label_offset_x = 70}

	$x_label_position = ($chr_start_array[$chr_count] + $chr_start_array[$chr_count+1])/2 - $label_offset_x;

	if ($chr_count == $total_no_chr) 
	{
		$x_label_position = ($chr_start_array[$chr_count] + $total_no_snps)/2 - $label_offset_x;
	}

	if ($chr_count % 2 == 0) 
	{
		print GG "set label $chr_count center '$chr_count' font 'arial,$x_axis_font_size' at $x_label_position,-$label_offset_y\n";
	}
	if ($chr_count % 2 == 1) 
	{
		print GG "set label $chr_count center '$chr_count' font 'arial,$x_axis_font_size' at $x_label_position,-($label_offset_y * 2)\n";
	}
}



print GG "plot ";
for ($chr_count = 1; $chr_count <= $total_no_chr; $chr_count++)
{
   print GG "'$infile' using 0:((\$1==$chr_count) ?(log10(\$9)*-1) : 1/0)";

	if ($chr_count % 2 == 0) {print GG " ls 1"}
	if ($chr_count % 2 == 1) {print GG " ls 2"}
	
	if ($chr_count < $total_no_chr)
	{
		print GG ", ";
	}
}
print GG "\n";


} # end of plink_assoc subroutine

###############################################################################################################


sub plink_unfiltered_assoc
{
################################################################################
# PLINK unfiltered assoc file subroutine                                       #
################################################################################

&print_message("Creating PLINK unfiltered association plot");

$infile = $prefix."_unfiltered.assoc";

$title = "$prefix:  PLINK unfiltered association data";
$total_no_chr = 3;

##########################################################
# plink_unfiltered_assoc: Open file to get some key data #
##########################################################
$chr_count = 0;
$line_count = 0;
$max_p_value = 0;

$last_chromosome = "1";
open (IN, "$infile") || die "Cannot open $infile";
while ($single_line = <IN>)
{
	$line_count = $line_count + 1;
	chomp $single_line;

	@item=split(/\s+/,$single_line);

	$chromosome = $item[1];
	$p_value = $item[9];

	if ($chromosome ne $last_chromosome)
	{
		$chr_count = $chr_count + 1;
		$chr_start_array[$chr_count-1] = $line_count;
	}

	if ($p_value > $max_p_value){$max_p_value = $p_value}

	$last_chromosome = $chromosome;
}

close (IN);

$total_no_chr = $chr_count -1;
$total_no_snps = $line_count;

print "    Total no chrs: $total_no_chr\tTotal no SNPs: $total_no_snps\n\n";


###############################################################
# plink assoc: Write the gnuplot commands to the command file #
###############################################################

print GG "set title  \"$title\" \n"; 
print GG "set key outside\n"; 
print GG "set xlabel  \"Chromosome\" offset 0,-1.2\n"; 
print GG "set ylabel  \"-logP\" \n" ; 
print GG "set nokey\n";
print GG "unset xtics\n";
print GG "set ytics\n";


print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 1\n";
print GG "set style line 2 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 3\n";
print GG "set style line 3 lt -1 lw 1 linecolor 2\n";

###########################
# assoc     x-axis labels #
###########################

$label_offset_y = $max_p_value / $label_offset_factor;

for ($chr_count = 1;$chr_count <=$total_no_chr; $chr_count++)
{
	if ($chr_count < 10) {$label_offset_x = 35}
	if ($chr_count >9) {$label_offset_x = 70}

	$x_label_position = ($chr_start_array[$chr_count] + $chr_start_array[$chr_count+1])/2 - $label_offset_x;

	if ($chr_count == $total_no_chr) 
	{
		$x_label_position = ($chr_start_array[$chr_count] + $total_no_snps)/2 - $label_offset_x;
	}

	if ($chr_count % 2 == 0) 
	{
		print GG "set label $chr_count center '$chr_count' font 'arial,$x_axis_font_size' at $x_label_position,-$label_offset_y\n";
	}
	if ($chr_count % 2 == 1) 
	{
		print GG "set label $chr_count center '$chr_count' font 'arial,$x_axis_font_size' at $x_label_position,-($label_offset_y * 2)\n";
	}
}


print GG "plot ";

for ($chr_count = 1; $chr_count <= $total_no_chr; $chr_count++)
{
   print GG "'$infile' using 0:((\$1==$chr_count) ?(log10(\$9)*-1) : 1/0)";

	if ($chr_count % 2 == 0) {print GG " ls 1"}
	if ($chr_count % 2 == 1) {print GG " ls 2"}
	
	if ($chr_count < $total_no_chr)
	{
		print GG ", ";
	}
}
print GG "\n";


} # end of plink_unfiltered_assoc subroutine

###############################################################################################################




sub plink_mperm_assoc
{
####################################################################################
# PLINK mperm assoc file subroutine                                                #
####################################################################################

&print_message("Creating Association Mperm plot");

$infile = $prefix.".assoc.mperm";




$title = "$prefix:  PLINK association mperm data";
$total_no_chr = 3;


###############################
# Open PLINK mperm.assoc file #
###############################
$chr_count = 0;
$line_count = 0;
$max_p_value = 0;

$last_chromosome = "1";
open (IN, "$infile") || die "Cannot open $infile";
while ($single_line = <IN>)
{
	$line_count = $line_count + 1;
	chomp $single_line;

	@item=split(/\s+/,$single_line);

	$chromosome = $item[1];
	$p_value = $item[4];

	if ($chromosome ne $last_chromosome)
	{
		$chr_count = $chr_count + 1;
		$chr_start_array[$chr_count-1] = $line_count;
	}

	if ($p_value > $max_p_value){$max_p_value = $p_value}

	$last_chromosome = $chromosome;
}

close (IN);

$total_no_chr = $chr_count - 1;
$total_no_snps = $line_count;

print "    Total no chrs: $total_no_chr\tTotal no SNPs: $total_no_snps\n\n";


###############################################################
# mperm.assoc: Write the gnuplot commands to the command file #
###############################################################

print GG "set title  \"$title\" \n" ; 
print GG "set key outside\n"; 
print GG "set xlabel  \"Chromosome\" offset 0,-1.2\n" ; 
print GG "set ylabel  \"-logP\" \n" ; 
print GG "set nokey\n";
print GG "unset xtics\n";
print GG "set ytics\n";


print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 1\n";
print GG "set style line 2 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 3\n";
print GG "set style line 3 lt -1 lw 1 linecolor 2\n";

################################
# mperm.assoc    x-axis labels #
################################

$label_offset_y = $max_p_value / $label_offset_factor;

for ($chr_count = 1;$chr_count <=$total_no_chr; $chr_count++)
{
	if ($chr_count < 10) {$label_offset_x = 35}
	if ($chr_count >9) {$label_offset_x = 70}

	$x_label_position = ($chr_start_array[$chr_count] + $chr_start_array[$chr_count+1])/2 - $label_offset_x;

	if ($chr_count == $total_no_chr) 
	{
		$x_label_position = ($chr_start_array[$chr_count] + $total_no_snps)/2 - $label_offset_x;
	}

	if ($chr_count % 2 == 0) 
		{
			print GG "set label $chr_count center '$chr_count' font 'arial,$x_axis_font_size' at $x_label_position,-$label_offset_y\n";
		}
		if ($chr_count % 2 == 1) 
		{
			print GG "set label $chr_count center '$chr_count' font 'arial,$x_axis_font_size' at $x_label_position,-($label_offset_y * 2)\n";
		}
}



print GG "plot ";
for ($chr_count = 1; $chr_count <= $total_no_chr; $chr_count++)
{
   print GG "'$infile' using 0:((\$1==$chr_count) ?(log10(\$4)*-1) : 1/0)";

	if ($chr_count % 2 == 0) {print GG " ls 1"}
	if ($chr_count % 2 == 1) {print GG " ls 2"}
	
	if ($chr_count < $total_no_chr)
	{
		print GG ", ";
	}
}
print GG "\n";


} # end of mperm.assoc subroutine

##############################################################################################################

sub plink_mds_plot
{
	####################################################################################
	# PLINK MDS file subroutine                                                        #
	####################################################################################

	&print_message("Creating MDS plot");

	$infile = $prefix.".mds";

	$title = "$prefix:  MDS plot";

	if ($use_samples_file eq "yes"){print "\tAffection file: $samples_used_file\n\n"}
	if ($use_pedfile eq "yes"){print "\tAffection file: $pedfile\n\n"}
	
	###########################
	# PLINK MDS plot          #
	###########################
	if ((!-e $pedfile)  && (! -e $samples_used_file))
	{ 
		print "  Neither ped file or samples used file were found so MDS plot cannot be completed.\n";
		return;
	}

	$total_no_lines = 0;
		
	if ($use_samples_file eq "yes")
	{
		
		open (SAMPLES, "$samples_used_file") || die "Cannot open $samples_used_file";
		$line_count = 0;
		
		print "    Reading sample used file  $samples_used_file to get affection data...\n\n";
		
		while ($single_line = <SAMPLES>)
		{
			$line_count = $line_count + 1;
			chomp $single_line;
			
			@item=split(/\s+/,$single_line);
			$id = $item[0];
			$affection  = $item[1];
			
			$id_array[$line_count] = $id;
			$affection_array[$line_count] = $affection;
			
			if ($affection eq "1"){$total_no_controls++}
			if ($affection eq "2"){$total_no_cases++}
		}
		$total_no_lines = $line_count;
		
		close SAMPLES;
	}
	
	
	#####################################################
	# If there is a pedfile but not a samples used file #
	#####################################################
	if ($use_pedfile eq "yes")
	{
		print "    Reading pedfile $pedfile to get affection data...\n\n";

		open (PED, "$pedfile") || die "Cannot open $pedfile";
		while ($single_line = <PED>)
		{
			$line_count = $line_count + 1;
			chomp $single_line;

			@item=split(/\s+/,$single_line);

			$pedigree = $item[0];
			$id = $item[1];
			$affection  = $item[5];

			$id_array[$line_count] = $id;
			$pedigree_array[$line_count] = $pedigree;
			$affection_array[$line_count] = $affection;

			if ($affection eq "1"){$total_no_controls++}
			if ($affection eq "2"){$total_no_cases++}
		}
		$total_no_lines = $line_count;
		
		close PED;
		
	} # end of reading ped file
	
	

	############################################################
	# If one of the affection files has been read in, go ahead #
	# (ie if the file is read and lines are counted)           #
	############################################################
	if ($total_no_lines > 0)
	{
		###########################################
		# Open temporary output file for MDS plot #
		###########################################
		$mds_tempfile = "mds_temp.txt";
		open (MDS_OUT, ">$mds_tempfile") || die "Cannot open $mds_tempfile";

		#######################################
		# Open mds file to get C1 and C2 data #
		#######################################
		$chr_count = 0;
		$line_count = 0;
		$max_p_value = 0;

		$last_chromosome = "1";
		open (IN, "$infile") || die "Cannot open $infile";
		while ($single_line = <IN>)
		{
			$line_count = $line_count + 1;
			chomp $single_line;

			@item=split(/\s+/,$single_line);

			$pedigree = $item[1];
			$id = $item[2];
			$C1 = $item[4];
			$C2 = $item[5];

			############################
			# Get affection of this id #
			############################
			if ($use_pedfile eq "yes")
			{
				for  ($check_count = 1; $check_count <=$total_no_lines; $check_count++)
				{
					if ($pedigree.$id eq $pedigree_array[$check_count].$id_array[$check_count])
					{
						$affection = $affection_array[$check_count];
					}
				}
			}
			if ($use_samples_file eq "yes")
			{
				for  ($check_count = 1; $check_count <=$total_no_lines; $check_count++)
				{
					if ($id eq $id_array[$check_count])
					{
						$affection = $affection_array[$check_count];
					}
				}
			}

			if ($line_count > 1){print MDS_OUT "$id\t$affection\t$C1\t$C2\n"}
		}


		close (IN);
		close (MDS_OUT);


		$total_no_chr = $chr_count - 1;
		$total_no_snps = $line_count;

		print "    Total no samples: $total_no_lines\tCases: $total_no_cases\tControls: $total_no_controls\n\n";

		############################################################
		# MDS plot: Write the gnuplot commands to the command file #
		############################################################

		print GG "set title  \"$title\" \n" ; 
		print GG "set xlabel  \"C1\"\n" ; 
		print GG "set ylabel  \"C2\"\n" ; 
		print GG "set key\n";
		print GG "unset xtics\n";
		print GG "unset ytics\n";

		print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.8 linecolor 2\n";
		print GG "set style line 2 lt 20 lw 1 pointtype 7 ps 0.8 linecolor 1\n";
		print GG "set style line 3 lt -1 lw 1 linecolor 2\n";

		print GG "plot ";

		print GG "'$mds_tempfile' using 3:((\$2=='1') ?\$4: 1/0) ls 1 title 'Controls',";

		print GG "'$mds_tempfile' using 3:((\$2=='2') ?\$4: 1/0) ls 2 title 'Cases'";
			
		print GG "\n";

	} # end of if mds file exists

} # end MDS file subroutine

##########################################################################################################################

sub plink_qq_plot
{
	####################################################################################
	# PLINK QQ plot subroutine                                                         #
	####################################################################################

	&print_message("Creating QQ plot");

	$infile = $prefix.".assoc.adjusted";

	$title = "$prefix:  QQ plot";


	###########################################################
	# Open assoc.adjusted file to get the UNADJ and QQ values #
	###########################################################

	$line_count = -1; # so that first line is line zero (not used later)

	print "    Reading file $infile...\n\n";

	open (IN, "$infile") || die "Cannot open $infile";
	while ($single_line = <IN>)
	{
		$line_count = $line_count + 1;
		chomp $single_line;

		@item=split(/\s+/,$single_line);

		$unadj = $item[3];

		if ($line_count == 1){$top_unadj = $unadj}

		if ($line_count > 0) {$qq_array[$line_count] = $unadj}
	}



	$total_no_lines = $line_count;


	######################################################
	# Calculate the inflation factor                     #
	# (depends on whether adjusted file is log10 or not) #
	######################################################

	if ($top_unadj > 1) {$log_values = "true"}
	if ($top_unadj < 1) {$log_values = "false"}

	print "    Top unadjusted score: $top_unadj\tLog values: $log_values\n\n";


	#########################################
	# First make chisq values from P values #
	#########################################
	for ($array_count = 1; $array_count <= $total_no_lines; $array_count++)
	{
		#if ($array_count % 1000 == 0 ){print "Array count: $array_count\n"}

		if ($log_values eq "true")
		{
			$P_value = 10**($qq_array[$array_count] * -1);
		
			#print "Array_count: $array_count\tqq_array: $qq_array[$array_count]\tP_value: $P_value\n";

			if ($P_value <= 1)
			{
				$chi_sq = Statistics::Distributions::chisqrdistr(1,$P_value);
				$chi_sq_array[$array_count] = $chi_sq;
			}
		}
		if ($log_values eq "false")
		{
			$P_value = $qq_array[$array_count];
		
			if ($P_value <= 1)
			{
				$chi_sq = Statistics::Distributions::chisqrdistr(1,$P_value);
				$chi_sq_array[$array_count] = $chi_sq;
			}
		}

		$sum_chi_sq = $sum_chi_sq + $chi_sq;

	}

	###############################################
	# Calculate the median of the observed values #
	###############################################

	$array_mean = $sum_chi_sq / $line_count;
	$array_mean = sprintf("%.3f",$array_mean);

	@chi_sq_array = reverse sort {$a <=> $b} (@chi_sq_array);

	$array_size = @chi_sq_array;
	$median = $chi_sq_array[int($array_size/2)];
	$lambda = $median / 0.456;

	$lambda_string = sprintf("%.3f",$lambda);

	print "\tMean chi_squared: \t$array_mean\n";
	print "\tMedian:           \t$median\n";
	print "\tLambda:           \t$lambda_string\n\n"; 


	###########################################################
	# QQ plot: Write the gnuplot commands to the command file #
	###########################################################
	print GG "unset label 1\n";
	print GG "unset label 2\n";
	print GG "set title  \"$title\" \n" ; 
	print GG "set xlabel  \"Expected\"\n" ; 
	print GG "set ylabel  \"Observed\"\n" ; 
	print GG "unset key\n";

	print GG "set label 1 'Inflation factor: $lambda_string' at graph 0.5,0.1 font 'Arial,12'\n";

	#print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.8 linecolor 2\n";
	#print GG "set style line 2 lt 20 lw 1 pointtype 7 ps 0.8 linecolor 1\n";

	print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 3\n";
	print GG "set style line 2 lt -1 lw 1 linecolor 2\n";


	print GG "plot ";

	print GG "'$infile' using 3:3 with lines ls 2 title 'NULL',";

	print GG "'$infile' using 5:3 with points ls 1 title 'QQ'";
		
	print GG "\n";


} # end of subroutine plink_qq_plot

#####################################################################################################################################
sub emmax_qq_plot
{
	####################################################################################
	# EMMAX QQ plot subroutine                                                         #
	####################################################################################

	&print_message("Creating EMMAX QQ plot");

	$infile = $prefix.".emmax_qq";

	$title = "$prefix:  EMMAX QQ plot";

	$title = "$prefix:  QQ plot of EMMAX -corrected chi-squared values";


	#################################################
	# Open .emmax_qq file to get the OBS_CHI values #
	#################################################

	$line_count = -1; # so that first line is line zero (not used later)

	print "    Reading file $infile...\n\n";

	open (IN, "$infile") || die "Cannot open $infile";
	while ($single_line = <IN>)
	{
		$line_count = $line_count + 1;
		chomp $single_line;

		@item=split(/\s+/,$single_line);


		$unadj = $item[0];  # This is the Observed chi-squareds in the first column

		if ($line_count == 1){$top_unadj = $unadj}

		if ($line_count > 0) {$chi_sq_array[$line_count] = $unadj}
	}


	###############################################
	# Calculate the median of the observed values #
	###############################################

	$array_mean = $sum_chi_sq / $line_count;
	$array_mean = sprintf("%.3f",$array_mean);

	@chi_sq_array = reverse sort {$a <=> $b} (@chi_sq_array);

	$array_size = @chi_sq_array;
	$median = $chi_sq_array[int($array_size/2)];
	$lambda = $median / 0.456;

	$lambda_string = sprintf("%.3f",$lambda);

	print "\tEMMAX Mean chi_squared: $array_mean\n\n";

	print "\tEMMAX Median: \t$median\n";
	print "\tEMMAX Lambda: \t$lambda_string\n\n\n"; 


	##################################################
	# Write the gnuplot commands to the command file #
	##################################################
	print GG "unset label 1\n";
	print GG "unset label 2\n";
	print GG "set terminal postscript color solid\n";
	print GG "set title  \"$title\" \n" ; 
	print GG "set key outside\n"; 
	print GG "set xlabel  \"Expected\" \n" ; 
	print GG "set ylabel  \"Observed\" \n" ; 
	print GG "set nokey\n";

	print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 3\n";
	print GG "set style line 2 lt -1 lw 1 linecolor 2\n";

	print GG "set label 1 'Inflation factor: $lambda_string' at graph 0.5,0.1 font 'Arial,12'\n";

	print GG "plot '$infile' using 3:3 with lines ls 2 title 'NULL', '$infile' using 2:3 with points ls 1 title 'QQ'\n";


} # end of emmax qq plot


sub emmax_assoc
	{
	################################################################################
	# EMMAX assoc file subroutine  (emmax_out)                                     #
	################################################################################

	&print_message("Creating EMMAX association plot");

	$infile = $prefix.".emmax_out";

	$title = "$prefix:  EMMAX-corrected association data";
	$total_no_chr = 3;

	###############################################
	# emmax_assoc: Open file to get some key data #
	###############################################
	$chr_count = 0;
	$line_count = 0;
	$max_p_value = 0;

	$last_chromosome = "1";
	open (IN, "$infile") || die "Cannot open $infile";
	while ($single_line = <IN>)
	{
		$line_count = $line_count + 1;
		chomp $single_line;

		@item=split(/\s+/,$single_line);

		$chromosome = $item[0];
		$p_value = $item[3];


		if ($chromosome ne $last_chromosome)
		{
			$chr_count = $chr_count + 1;
			$chr_start_array[$chr_count-1] = $line_count;
		}

		if ($p_value > $max_p_value){$max_p_value = $p_value}

		$last_chromosome = $chromosome;
	}

	close (IN);

	$total_no_chr = $chr_count -1;
	$total_no_snps = $line_count;

	print "    Total no chrs: $total_no_chr\tTotal no SNPs: $total_no_snps\n\n";


	##############################################################
	# Emmax_assoc Write the gnuplot commands to the command file #
	##############################################################

	print GG "set title  \"$title\" \n" ; 
	print GG "set key outside\n"; 
	print GG "set xlabel  \"Chromosome\" offset 0,-1.2\n" ; 
	print GG "set ylabel  \"-logP\" \n" ; 
	print GG "set nokey\n";
	print GG "unset xtics\n";
	print GG "set ytics\n";


	print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 1\n";
	print GG "set style line 2 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 3\n";
	print GG "set style line 3 lt -1 lw 1 linecolor 2\n";

	#################################
	# Emmax assoc     x-axis labels #
	#################################

	$label_offset_y = $max_p_value / $label_offset_factor;

	for ($chr_count = 1;$chr_count <=$total_no_chr; $chr_count++)
	{

		$x_label_position = ($chr_start_array[$chr_count] + $chr_start_array[$chr_count+1])/2;

		if ($chr_count == $total_no_chr) 
		{
			$x_label_position = ($chr_start_array[$chr_count] + $total_no_snps)/2;
		}

		if ($chr_count % 2 == 0) 
		{
			print GG "set label $chr_count center '$chr_count' font 'arial,$x_axis_font_size' at $x_label_position,-$label_offset_y\n";
		}
		if ($chr_count % 2 == 1) 
		{
			print GG "set label $chr_count center '$chr_count' font 'arial,$x_axis_font_size' at $x_label_position,-($label_offset_y * 2)\n";
		}
	}


	print GG "plot ";
	for ($chr_count = 1; $chr_count <= $total_no_chr; $chr_count++)
	{
	   #print GG "'$infile' using 0:((\$1==$chr_count) ?(log10(\$4)*-1) : 1/0)"; # used to use the p_value column. Now uses minus log P
	   print GG "'$infile' using 0:((\$1==$chr_count) ?(\$5) : 1/0)";

		if ($chr_count % 2 == 0) {print GG " ls 1"}
		if ($chr_count % 2 == 1) {print GG " ls 2"}
		
		if ($chr_count < $total_no_chr)
		{
			print GG ", ";
		}
	}
	print GG "\n";

} # end of emmax_assoc subroutine




sub fmm_assoc
	{
	################################################################################
	# EMMAX assoc file subroutine  (emmax_out)                                     #
	################################################################################

	&print_message("Creating FMM association plot");

	$infile = $fmm_out_file;

	#######################################################################
	# First reformat the FMM.out file so it looks like the emmax_out file #
	#######################################################################
	
	open (IN, "$infile") || die "Cannot open $infile";
	open (FMM_TEMP ,">$fmm_temp_file") || die "can't open $fmm_temp_file\n" ;

	# headers for columns
	print FMM_TEMP "CHR\tSNP\tPOS\tP_value\tminus_logP_FMM\n";
	$line_count = 0;

	#########################
	# Open the FMM out file #
	#########################
	while ($single_line = <IN>)
	{
		$line_count = $line_count + 1;
		chomp $single_line;
		@item=split(/\s+/,$single_line);

		#INDEX	CHR	SNP	POS	FMM_CHI
		$chromosome = $item[1];
		$snp_name = $item[2];
		$position = $item[3];
		$chi_sq = $item[4];

		$chi_sq_array[$line_count] = $chi_sq;

		# Convert to p-value and minus_log_P
		$p_value = Statistics::Distributions::chisqrprob(1,$chi_sq);
		$minus_log_p = log($p_value)/log(10) * -1;

		#CHR	SNP	POS	P_value	minus_logP_EMMAX
		if ($line_count > 1){print FMM_TEMP "$chromosome\t$snp_name\t$position\t$p_value\t$minus_log_p\n";}
	}

	$total_lines_fmm_out = $line_count;

	close IN;
	close FMM_TEMP;

	###########################################
	# Sort the chi_sq array (for the QQ plot) #
	###########################################
	@chi_sq_array_reverse = reverse sort {$a <=> $b} (@chi_sq_array);


	#################################
	# Print first 10 lines of array #
	#################################

	&print_message("Showing the first 10 lines of the reverse sorted chi-sq array");

	for($line_count = 0; $line_count <=10; $line_count++)
	{
		print "$line_count\t$chi_sq_array_reverse[$line_count]\n";
	}

	print "\n";

	##############################################
	# FMM QQ plot                                #
	# Now calculate Exp_p with this formula:     #
	#                                            #
	# Exp-p = (rank - 0.5)/total_no_of_lines     #
	# and then getting the CHISQ of that P-value #
	##############################################

	for($array_count = 0; $array_count <= $total_lines_fmm_out; $array_count++)
	{
		$rank = $array_count + 1;
		$expected_P_value = ($rank - 0.5) / $total_lines_fmm_out;

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

	###########################################
	# Open file for output for FMM QQ plot    #
	###########################################

	open (OUT, ">$qq_file") || die "Cannot open $qq_file";

	print OUT "FMM\tBLANK\tEXP_CHI\tOBS_CHI\n";

	for ($line_count = 0; $line_count <= $total_lines_fmm_out; $line_count++)
	{
		if (($chi_sq_array_reverse[$line_count] * 1) > 0)
		{
			print OUT "$chi_sq_array_reverse[$line_count]\t\t$chi_sq_array_expected[$line_count]\t$chi_sq_array_reverse[$line_count]\n";
		}
	}

	print OUT "\t\t\t\t\n";
	print OUT "\t\t\t$lambda_string";

	close OUT;



	################################################################################
	# Now fmm_temp_out looks like an emmax_out file so plotting should be the same #
	################################################################################

	$title = "$prefix:  Fast Mixed Model-corrected association data";
	$total_no_chr = 3;

	###############################################
	# fmm_assoc: Open file to get some key data #
	###############################################
	$chr_count = 0;
	$line_count = 0;
	$max_p_value = 0;
	$max_minus_log_p = 0;

	$last_chromosome = "1";
	open (IN, "$fmm_temp_file") || die "Cannot open $fmm_temp_file";
	while ($single_line = <IN>)
	{
		$line_count = $line_count + 1;
		chomp $single_line;

		@item=split(/\s+/,$single_line);

		$chromosome = $item[0];
		$p_value = $item[3];
		$minus_log_p = $item[4];

		if ($chromosome ne $last_chromosome)
		{
			$chr_count = $chr_count + 1;
			$chr_start_array[$chr_count-1] = $line_count;
		}

		if ($minus_log_p > $max_minus_log_p){$max_minus_log_p = $minus_log_p}

		$last_chromosome = $chromosome;
	}

	close (IN);


	$total_no_chr = $chr_count -1;
	$total_no_snps = $line_count;

	print "    Total no chrs: $total_no_chr\tTotal no SNPs: $total_no_snps\n\n";


	##############################################################
	# fmm_assoc Write the gnuplot commands to the command file   # carrot
	##############################################################

	print GG "set title  \"$title\" \n" ; 
	print GG "set key outside\n"; 
	print GG "set xlabel  \"Chromosome\" offset 0,-1.2\n" ; 
	print GG "set ylabel  \"-logP\" \n" ; 
	print GG "set nokey\n";
	print GG "unset xtics\n";
	print GG "set ytics\n";

	print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 1\n";
	print GG "set style line 2 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 3\n";
	print GG "set style line 3 lt -1 lw 1 linecolor 2\n";


	#################################
	# FMM assoc     x-axis labels   #
	#################################

	$label_offset_y = 1 / $label_offset_factor;
	$bottom_axis_threshold = $max_minus_log_p/$size_of_dot;


	print "Max -logP:             \t$max_minus_log_p\n";
	print "Bottom axis threshold: \t$bottom_axis_threshold\n\n";

	print "Label offset factor:   \t$label_offset_factor\n";
	print "Label offset y:        \t$label_offset_y\t\t(max_p_value/label_offset_factor)\n\n";


	for ($chr_count = 1;$chr_count <=$total_no_chr; $chr_count++)
	{
		#if ($chr_count < 10) {$label_offset_x = 500}
		#if ($chr_count > 9) {$label_offset_x = 1000}

		$x_label_position = ($chr_start_array[$chr_count] + $chr_start_array[$chr_count+1])/2;

		if ($chr_count == $total_no_chr) 
		{
			$x_label_position = ($chr_start_array[$chr_count] + $total_no_snps)/2;
		}

		if ($chr_count % 2 == 0) 
		{
			print GG "set label $chr_count center '$chr_count' font 'arial,$x_axis_font_size' at $x_label_position,-$label_offset_y\n";
		}
		if ($chr_count % 2 == 1) 
		{
			print GG "set label $chr_count center '$chr_count' font 'arial,$x_axis_font_size' at $x_label_position,-($label_offset_y * 2)\n";
		}
	}

	#set style rectangle {front|back} {fillcolor <colorspec>} {fs <fillstyle>} {lw|linewidth <lw>}

	print GG "set style rectangle front fillcolor rgb 'blue' fs solid noborder\n";
	print GG "set object 1 rect from graph 0,-0.005 to graph 1,-0.01 front\n";

	#print GG "set arrow from graph 0,0 to graph 0,1 nohead lc rgb \'red\' lw 1\n";  ##<<< could use a white line


	print GG "plot ";
	for ($chr_count = 1; $chr_count <= $total_no_chr; $chr_count++)
	{
	  # print GG "'$fmm_temp_file' using 0:((\$1==$chr_count) ?($data_point_offset_y + log10(\$4)*-1) : 1/0)"; # used to use the p_value column. Now uses minus log P
	  
	  # emmax_version:   print GG "'$infile' using 0:((\$1==$chr_count) ?(\$5) : 1/0)";


	  print GG "'$fmm_temp_file' using 0:((\$1==$chr_count) ? (\$5 + ((\$5 < $bottom_axis_threshold)) * 0.1) : 1/0)";

	  ###############################################################################################################
	  # Note: the expression ((x == y) ? A : B) means if x = y then do A else do B
	  #####################

		if ($chr_count % 2 == 0) {print GG " ls 1"}
		if ($chr_count % 2 == 1) {print GG " ls 2"}
		
		if ($chr_count < $total_no_chr)
		{
			print GG ", ";
		}
	}
	print GG "\n";


	######################
	# Now do FMM qq_plot #
	######################

	&print_message("Creating FMM QQ plot");

	$title = "$prefix:  QQ plot of FMM -corrected chi-squared values";

	print "\tFMM Mean chi_squared: \t$array_mean\n";
	print "\tFMM Median:           \t$median\n";
	print "\tFMM Lambda:           \t$lambda_string\n\n\n"; 


	################################################################
	# Write the gnuplot commands to the command file - FMM QQ plot #
	################################################################
	print GG "set terminal postscript color solid\n";
	print GG "set title  \"$title\" \n" ; 
	print GG "set key outside\n"; 
	print GG "set xlabel  \"Expected\" \n" ; 
	print GG "set ylabel  \"Observed\" \n" ; 
	print GG "set nokey\n";

	print GG "set style line 1 lt 20 lw 1 pointtype 7 ps 0.5 linecolor 3\n";
	print GG "set style line 2 lt -1 lw 1 linecolor 2\n";

	print GG "set label 1 'Inflation factor: $lambda_string' at graph 0.5,0.1\n";

	print GG "plot '$qq_file' using 3:3 with lines ls 2 title 'NULL', '$qq_file' using 2:3 with points ls 1 title 'QQ'\n";



} # end of fmm_assoc subroutine


#####################################################################################################################################
sub title_page
{
print GG "set style line 1 lt rgb '#FFFFFF'\n";

#######################
# Various title lines #
#######################
$date = localtime;

print GG "set label 1 'Summary of plots for:  \t$prefix' at graph 0.1,0.6 front font 'Times New Roman,26'\n";
print GG "set label 2 'Date:                  \t$date' at graph 0.1,0.5 front font 'Arial,12'\n";

print GG "unset key\n";
print GG "unset xtics\n";
print GG "unset ytics\n";
#print GG "unset border\n";

print GG "plot x ls 1\n";

}

#####################################################################################################################################
sub last_page
{
	print GG "set style line 1 lt rgb '#FFFFFF'\n";

	#######################
	# Various title lines #
	#######################
	$date = localtime;
	$xpos = 0.45;
	$ypos = 0.85;
	$yinc = 0.05;

	print GG "set label 1 'Summary of plots for:  \t$prefix' at graph $xpos,$ypos front font 'Times New Roman,18'\n";
	$ypos = $ypos - $yinc;
	print GG "set label 2 'Files found and plotted...' at graph $xpos,$ypos front font 'Arial,12'\n";
	$ypos = $ypos - $yinc;


	####################
	# Association file #
	####################
	
	$ypos = $ypos - $yinc;
	if (-e $prefix.".assoc")
	{ 
		print GG "set label 3 'Association: \t$prefix.assoc \tFile plotted' at graph $xpos,$ypos front font 'Arial,12'\n";
	}
	else{
		print GG "set label 3 'Association: \t$prefix.assoc \tFile not found' at graph $xpos,$ypos front font 'Arial,12'\n";
	}

	####################
	# MPERM file       #
	####################
	$ypos = $ypos - $yinc;
	if (-e $prefix.".assoc.mperm")
	{ 
		print GG "set label 4 'Association mperm: \t$prefix.assoc.mperm \tFile plotted' at graph $xpos,$ypos front font 'Arial,12'\n";
	}
	else{
		print GG "set label 4 'Association mperm: \t$prefix.assoc.mperm \tFile not found' at graph $xpos,$ypos front font 'Arial,12'\n";
	}

	####################
	# MDS file         #
	####################
	$ypos = $ypos - $yinc;
	if (-e $prefix.".mds")
	{ 
		if ($use_pedfile eq "yes")
		{
			print GG "set label 5 'MDS plot: \t$prefix.mds \tFile plotted using affection data from $pedfile' at graph $xpos,$ypos front font 'Arial,12'\n";
		}
		if ($use_samples_file eq "yes")
		{
			print GG "set label 5 'MDS plot: \t$prefix.mds \tFile plotted using affection data from $samples_used_file' at graph $xpos,$ypos front font 'Arial,12'\n";
		}
		if (($use_pedfile eq "no") && ($use_samples_file eq "no"))
		{
			print GG "set label 5 'MDS plot: \t$prefix.mds \tFile not plotted as affection files not found' at graph $xpos,$ypos front font 'Arial,12'\n";
		}
	}
	else
	{
		print GG "set label 5 'MDS plot: \t$prefix.mds \tFile not found' at graph $xpos,$ypos front font 'Arial,12'\n";
	}

	####################
	# QQ plot file     #
	####################
	$ypos = $ypos - $yinc;
	if (-e $prefix.".assoc.adjusted"){ 
		print GG "set label 6 'QQ plot: \t$prefix.assoc.adjusted \tFile plotted' at graph $xpos,$ypos front font 'Arial,12'\n";
	}
	else{
		print GG "set label 6 'QQ plot: \t$prefix.assoc.adjusted \tFile not found' at graph $xpos,$ypos front font 'Arial,12'\n";
	}


	######################
	# EMMAX QQ plot file #
	######################
	$ypos = $ypos - $yinc;
	if (-e $prefix.".emmax_qq"){ 
		print GG "set label 7 'EMMAX QQ plot: \t$prefix.emmax_qq \tFile plotted' at graph $xpos,$ypos front font 'Arial,12'\n";
	}
	else{
		print GG "set label 7 'EMMAX QQ plot: \t$prefix.emmax_qq \tFile not found' at graph $xpos,$ypos front font 'Arial,12'\n";
	}


	#############################
	# EMMAX manhattan plot file #
	#############################
	$ypos = $ypos - $yinc;
	if (-e $prefix.".emmax_out"){ 
		print GG "set label 8 'EMMAX manhattan plot: \t$prefix.emmax_out \tFile plotted' at graph $xpos,$ypos front font 'Arial,12'\n";
	}
	else{
		print GG "set label 8 'EMMAX manhattan plot: \t$prefix.emmax_out \tFile not found' at graph $xpos,$ypos front font 'Arial,12'\n";
	}


	#############################
	# FMM manhattan plot file   #
	#############################
	$ypos = $ypos - $yinc;
	if (-e $fmm_out_file){ 
		print GG "set label 9 'Fast Mixed Model manhattan plot: \t$fmm_out_file \tFile plotted' at graph $xpos,$ypos front font 'Arial,12'\n";
	}
	else{
		print GG "set label 9 'Fast Mixed Model manhattan plot: \t$fmm_out_file \tFile not found' at graph $xpos,$ypos front font 'Arial,12'\n";
	}



	##################################
	# Details such as number of SNPs #
	##################################
	$ypos = $ypos - $yinc;
	$ypos = $ypos - $yinc;
	print GG "set label 10 'Total number of SNPs: \t$total_no_snps' at graph $xpos,$ypos front font 'Arial,12'\n";
	$ypos = $ypos - $yinc;
	print GG "set label 11 'Total number of chromosomes: \t$total_no_chr' at graph $xpos,$ypos front font 'Arial,12'\n";







	print GG "unset key\n";
	print GG "unset xtics\n";
	print GG "unset ytics\n";
	print GG "unset xlabel\n";
	print GG "unset ylabel\n";
	print GG "unset title\n";
	#print GG "unset border\n";

	print GG "plot x ls 1\n";

} # last page summary of files used

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

