#!/usr/bin/perl -w


#########################################################################
#									#      
#	CNV ANALYSIS v3							#     
#									#
#	THIS PERL SCRIPT WILL COMPARE TWO BAM FILES TO LOOK FOR CNVs 	#
#									#
#########################################################################

#############################
# Oliver Forman JUN 2010    #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# oliver.forman@aht.org.uk  #
#############################


use strict;
use Getopt::Std ;
use File::Basename ;

####################
# Define variables #
####################

my $proceed	= "";
my $word	= "";
my $command	= "";
my $to		= "";
my $name	= "";
my $bam1	= "";
my $bam2	= "";
my $focus	= "";
my $focus_type	= "";
my $gregion	= "";
my $F		= "";
my @F		= "";	


########################
# Define non variables #
########################


my $title='Perl Mail demo';
my $from= 'oliver.forman@aht.org.uk';
my $reads4 = 'out.fastq';

my $inter = '$F[2]\t$F[3]';


#START TIMERS

BEGIN { our $start_run = time(); }


#################
#TURN LOGGER ON #
#################

$| = 1;

open(STDOUT, "| tee log.rtf");



    use Term::ANSIColor;
    print color 'magenta';

$command =  "clear";
system("$command");

print"  ___ _  ___   __    _             _         _    \n";
print" / __| '| ' ' / /   /_'  _ _  __ _| |_  _ __(_)___\n";
print"| (__| .` |' V /   / _ '| ' '/ _` | | || (_-< (_-<\n";
print" '___|_|'_| '_/   /_/ '_'_||_'__,_|_|'_, /__/_/__/\n";
print"                                     |__/		\n";
 


 
print "Read depth calculator and CNV detection pipeline\n";

print "\n";

print color 'reset';


###############
# Name Analysis
###############


print "\nPlease enter a name for this analysis (with no spaces):    ";
$name = <STDIN>;
chomp $name;


#####################
#Check if name exist
#####################

if (-e "results_$name/log.rtf")
{ 
	print "\nA results folder with the name results_$name already exists!\n\n";
	exit;
}


##############
#EMAIL ADDRESS
##############


print "\nPlease enter your email address:    ";
$to = <STDIN>;
chomp $to;


#############
#ENTER BAM 1
#############

print "\nPlease enter the name of bam file 1 (usually a case):    ";
$bam1 = <STDIN>;
chomp $bam1;

#############
#ENTER BAM 1
#############

print "\nPlease enter the name of bam file 2 (usually a control):    ";
$bam2 = <STDIN>;
chomp $bam2;


###########################
#FOCUS TO A GENOME REGION?
###########################


print "\nWould you like to focus on a particular genomic region?\n\n";

print " Enter 1 for yes\n";
print " Enter 2 for no\n\n";


$focus_type = <STDIN>;
chomp $focus_type;

if ($focus_type == 1){$focus = "focus_yes"}
if ($focus_type == 2){$focus = "focus_no"}


if ($focus eq "focus_yes")

{

print "\nPlease define your genomic region of interest (eg chr7:12000000-13000000)\n\n";
$gregion = <STDIN>;
chomp $gregion;

}

########
#SUMMARY
########


$command =  "clear";
system("$command");

print color 'bold white';

print"-------";
print"SUMMARY";
print"-------\n\n";	

print color 'reset';

print"ANALYSIS NAME:\t $name \n";
print"EMAIL:\t\t $to \n";
print"BAM FILE 1:\t $bam1 \n";
print"BAM FILE 2:\t $bam2 \n";
if ($focus eq "focus_yes")

{
print "Analysis region: $gregion \n";
}

if ($focus eq "focus_no")

{
print "All reads will be analysed \n";
}

print "\nPlease press enter to proceed (or 'Q' to quit):      ";
$proceed = <STDIN>;
chomp $proceed; 
 
if (lc $proceed eq "q")
{
 exit;
} 


#-------------------------------------------------------------------------------------------------

####################
#Make new Bam files#
####################

if ($focus eq "focus_yes")
{

print color 'bold white';

print "\nSelecting reads in user defined region for bam1...\n\n";

print color 'reset';

$command =  "/opt/samtools/samtools view $bam1 $gregion -b -o focus1.bam";

print("$command\n\n");
system("$command");



print color 'bold white';

print "\nSelecting reads in user defined region for bam2...\n\n";

print color 'reset';

$command =  "/opt/samtools/samtools view $bam2 $gregion -b -o focus2.bam";

print("$command\n\n");
system("$command");

######################
#Create best hits file
######################

print color 'bold white';

print "\nCreating best hits file for bam1...\n\n";

print color 'reset';

$command =  "/opt/samtools/samtools view -F 4 focus1.bam | perl -lane 'print \"$inter\"' > case.hits";
print("$command\n");
system("$command");



print color 'bold white';

print "\nCreating best hits file for bam2...\n\n";

print color 'reset';

$command =  "/opt/samtools/samtools view -F 4 focus2.bam | perl -lane 'print \"$inter\"' > control.hits";
print("$command\n");
system("$command");

}


#-------------------------------------------------------------------------------------------------


if ($focus eq "focus_no") 
{

######################
#Create best hits file
######################

print color 'bold white';

print "\nCreating best hits file for bam1...\n\n";

print color 'reset';

$command =  "/opt/samtools/samtools view -F 4 $bam1 | perl -lane 'print \"$inter\"' > case.hits";
print("$command\n");
system("$command");



print color 'bold white';

print "\nCreating best hits file for bam2...\n\n";

print color 'reset';

$command =  "/opt/samtools/samtools view -F 4 $bam2 | perl -lane 'print \"$inter\"' > control.hits";
print("$command\n");
system("$command");

}

#---------------------------------------------------------------------------------------------------

############
#Run cnv-seq
############



print color 'bold white';

print "\nRunning cnv-seq for window of 1k,5k,10k,20k...\n\n";

print color 'reset';

$command =  "perl /opt/cnv-seq/cnv-seq.pl --test case.hits  --ref control.hits --genome-size 3000000 --window-size 1000";

print("$command\n");
system("$command");

$command =  "perl /opt/cnv-seq/cnv-seq.pl --test case.hits  --ref control.hits --genome-size 3000000 --window-size 5000";

print("$command\n");
system("$command");


$command =  "perl /opt/cnv-seq/cnv-seq.pl --test case.hits  --ref control.hits --genome-size 3000000 --window-size 10000";

print("$command\n");
system("$command");


$command =  "perl /opt/cnv-seq/cnv-seq.pl --test case.hits  --ref control.hits --genome-size 3000000 --window-size 20000";

print("$command\n");
system("$command");


$command =  "perl /opt/cnv-seq/cnv-seq.pl --test case.hits  --ref control.hits --genome-size 3000000 --window-size 40000";

print("$command\n");
system("$command");


$command = "mkdir results_$name";
print("$command\n");
system("$command");


&move ("case.hits-vs-control.hits.window-1000.minw-4.cnv");
&move ("case.hits-vs-control.hits.window-5000.minw-4.cnv");
&move ("case.hits-vs-control.hits.window-10000.minw-4.cnv");
&move ("case.hits-vs-control.hits.window-20000.minw-4.cnv");
&move ("case.hits-vs-control.hits.window-40000.minw-4.cnv");
&move ("case.hits-vs-control.hits.window-1000.minw-4.count");
&move ("case.hits-vs-control.hits.window-5000.minw-4.count");
&move ("case.hits-vs-control.hits.window-10000.minw-4.count");
&move ("case.hits-vs-control.hits.window-20000.minw-4.count");
&move ("case.hits-vs-control.hits.window-40000.minw-4.count");
&move ("control.hits");
&move ("case.hits");
&move ("log.rtf");

my $end_run = time();
my $run_time = $end_run - our $start_run;


open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: CNV analysis COMPLETE ($name)\n\n";
## Mail Body
print MAIL "CNV analysis ($name) is complete\n\n";
print MAIL "Run time : $run_time seconds\n";
printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
close(MAIL);


print "\nANALYSIS COMPLETE! YOU HAVE BEEN NOTIFIED BY EMAIL\n\n";


print "Run time : $run_time seconds\n";
printf "%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];



#############################
# Turn off logging          #
#############################

close(STDOUT);

#############################################
#                                           #
# Subroutine to move file to results folder #
#                                           #
#############################################

sub move
{
	my $suffix = "";	

	$suffix = $_[0];
	$command = "mv  $suffix results_$name/$suffix";
	print("$command\n");
	system("$command");
}
