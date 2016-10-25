#!/usr/bin/perl -w


#########################################################################
#									#      
#	NGS ANALYSIS v 15.3						#     
#									#
#	THIS PERL SCRIPT WILL ANALYSE 454 and ILLUMINA NGS DATA		#
#									#
#########################################################################

#############################
# Oliver Forman Sep 2010    #
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

my $word	= "";
my $command	= "";
my $ref		= "";
my $reads	= "";
my $mem		= "";
my $proceed	= "";
my $data_type	= "";
my $data	= "";
my $reads2	= "";
my $to		= "";
my $name	= "";
my $platform	= "";
my $platform_type = "";
my $reads3 = 	"";
my $suffix	= "";
my $dup_type	= "";
my $var_type	= "";
my $var		= "";
my $maq		= "";
my $maq_type		= "";
my $ref_size	= "";
my $ref_ans	= "";
my $sv_type	= "";
my $sv		= "";
my $chromo	= "";
my $insert	= "";
my $region	= "";
my $target	= "";
my $target_ans	= "";
my $file_type	= "";
my $file	= "";
my $bam_file	= "";

########################
# Define non variables #
########################


my $title='Perl Mail demo';
my $from= 'oliver.forman@aht.org.uk';
my $subject='NGS ANALYSIS';
my $reads4 = 'out.fastq';
my $dup		= "yes";



#START TIMERS

BEGIN { our $start_run = time(); }

BEGIN { our $start_run2 = time(); }

BEGIN { our $start_run3 = time(); }


#################
#TURN LOGGER ON #
#################

$| = 1;

open(STDOUT, "| tee log.rtf");



    use Term::ANSIColor;
    print color 'bold cyan';


$command =  "clear";
system("$command");


print"\n\n";
print"*************************************************************\n";
print"*  _   _  ____ ____                      _           _      *\n";
print"* | ` | |/ ___/ ___|    __ _ _ __   __ _| |_   _ ___(_)___  *\n";
print"* |  `| | |  _'___ `   / _` | '_ ` / _` | | | | / __| / __| *\n";
print"* | |`  | |_| |___) | | (_| | | | | (_| | | |_| `__ ` `__ ` *\n";
print"* |_| `_|`____|____/   `__,_|_| |_|`__,_|_|`__, |___/_|___/ *\n";
print"*                                          |___/            *\n";
print"*************************************************************\n";
print"\n";

print color 'reset';

    use Term::ANSIColor;
    print color 'magenta';


print "                 ##########################\n";
print color 'bold white';
print "                 NGS data analysis pipeline\n";
print color 'reset';
print color 'magenta';
print "                 ##########################\n";
print "\n";

print color 'reset';


#############################
# Name Analysis and email   #
#############################

print "\n~~~~~~~~~~~~~~~";
print "\nYour details...";
print "\n~~~~~~~~~~~~~~~\n";

print "\nPlease enter a name for this analysis (with no spaces):    ";
$name = <STDIN>;
chomp $name;


########################
# Check if name exist #
########################

if (-e "results_$name/log.rtf")
{ 
	print "\nA results folder with the name results_$name already exists!\n\n";
	exit;
}



print "\nPlease enter your email address:    ";
$to = <STDIN>;
chomp $to;



##############################
#Ask if the data 454 or Illumina?
##############################

print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " Please select your platform type?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   Enter 1 for Illumina\n";
print "   Enter 2 for 454\n\n";


$platform_type = <STDIN>;
chomp $platform_type;

if ($platform_type == 1){$platform = "Illumina"}
if ($platform_type == 2){$platform = "454"}

if ($platform eq "Illumina")
{

	print "\nYou have selected Illumina\n"

}

if ($platform eq "454")
{

	print "\nYou have selected 454\n"

}



##############################
##############################
# START OF THE ILLUMINA SCRIPT
##############################
##############################


if ($platform eq "Illumina")

{


#
#ASK IF STARTING FILES ARE BAM OR FASTQ
#

print "\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " Are you starting with raw reads (.fastq) or aligned reads (.bam)?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   Enter 1 for RAW\n";
print "   Enter 2 for ALIGNED\n\n";

$file_type = <STDIN>;
chomp $file_type;

if ($file_type == 1){$file = "fastq"}
if ($file_type == 2){$file = "bam"}






if ($file eq "bam")
{


print "\nPlease input the name of your aligned reads file (aligned.bam):      ";
$bam_file = <STDIN>;
chomp $bam_file;


}

















































##############################
#Ask if the data is SE or PE data?
##############################


print "\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " Do you have a Paired-end or Single-end dataset?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   Enter 1 for SE\n";
print "   Enter 2 for PE\n\n";


$data_type = <STDIN>;
chomp $data_type;

if ($data_type == 1){$data = "SE"}
if ($data_type == 2){$data = "PE"}


if ($data eq "SE")
{

	print "\nYou have selected Single-end\n"

}

if ($data eq "PE")
{

	print "\nYou have selected Paired-end\n"

}



##################################
# Define data files              #
##################################

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
print "\nDefine data files and memory setting";
print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";


if ($file eq "fastq")
{

#############################
#Ask if ref is whole genome #
#############################

print "\nIs your reference genome >2 Gb?\n\n";


print "   Enter 1 for YES\n";
print "   Enter 2 for NO\n\n";


$ref_size = <STDIN>;
chomp $ref_size;

if ($ref_size == 1){$ref_ans = "whole"}
if ($ref_size == 2){$ref_ans = "partial"}

}

print "\nPlease input the name of your reference fasta file (eg ref):      ";
$ref = <STDIN>;
chomp $ref;






#####################################
# Add fasta if user hasn't added it #
#####################################

if (index($ref,".fasta") == -1 )
{
	$ref = $ref.".fasta";
}

if ($file eq "fastq")
{

if ($data eq "SE")
{
print "\nPlease input the name of your fastq sequence file (eg reads):      ";
$reads = <STDIN>;
chomp $reads;
}

#####################################
# Add fastq if user hasn't added it #
#####################################
if (index($reads,".fastq") == -1 )
{
	$reads = $reads.".fastq";
}


if ($data eq "PE")
{

print "\nPlease input the name of your 1st fastq sequence file (eg reads1):      ";
$reads = <STDIN>;
chomp $reads;

#####################################
# Add fastq if user hasn't added it #
#####################################
if (index($reads,".fastq") == -1 )
{
	$reads = $reads.".fastq";
}


print "\nPlease input the name of your 2nd fastq sequence file (eg reads2):      ";
$reads2 = <STDIN>;
chomp $reads2;

#####################################
# Add fastq if user hasn't added it #
#####################################

if (index($reads2,".fastq") == -1 )
{
	$reads2 = $reads2.".fastq";
}

}

}


print "\nPlease enter the memory setting for this analysis (default is -Xmx2g):      ";
$mem = <STDIN>;
chomp $mem;
 
if ($mem eq "")
{
 $mem = "-Xmx2g"
}




















if ($file eq "fastq")
{


########################
# Check if file exists #
########################

if (! -e "$ref")
{ 
	print "\nFile $ref does not exist\n\n";
	exit;
} 
if (! -e "$reads")
{ 
	print "\nFile $reads does not exist\n\n";
	exit;
} 

if  ($data eq "PE")
{
if (! -e "$reads2")
{ 
	print "\nFile $reads2 does not exist\n\n";
	exit;
} 
}


}


###########################################################################
# Would you like to focus the analysis on a specific region of the genome #
###########################################################################


print "\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Would you like to focus on alignments to specific region of the genome?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

print "(Eg sequence capture libraries or managing whole genome sequencing)\n\n";


print "   Enter 1 for YES\n";
print "   Enter 2 for NO\n\n";

$target = <STDIN>;
chomp $target;

if ($target == 1){$target_ans = "yes"}
if ($target == 2){$target_ans = "no"}


if ($target_ans eq "yes")

{
print "\nPlease define your region of interest (eg chr5:21000000-23000000):      ";
$region = <STDIN>;
chomp $region;

}



######################
# Annotate variants  #
######################

print "\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " Would you like annotate your SNP and INDEL calls?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

print "Variants will be annotated using information available in Ensembl\n\n"; 

print "   Enter 1 for YES\n";
print "   Enter 2 for NO\n\n";

$var_type = <STDIN>;
chomp $var_type;

if ($var_type == 1){$var = "yes"}
if ($var_type == 2){$var = "no"}


if ($var eq "yes")
{

	print "\nVariants will be annotated using the ensembl database\n";

}

if ($var eq "no")
{

	print "\nVariants will not be annotated\n";

}

if ($file eq "fastq")
{


################
#MAQ analysis
################


print "\n";
print "~~~~~~~~~~~~~~~~~~~~~~\n";
print " Include MAQ analysis?\n";
print "~~~~~~~~~~~~~~~~~~~~~~\n\n";


print color 'bold red';

print "PLEASE NOTE-\tMAQ ONLY WORKs ON READS SHORTER THAN 64 bp\n\n";

print color 'reset';

print "   Enter 1 for YES\n";
print "   Enter 2 for NO\n\n";


$maq_type = <STDIN>;
chomp $maq_type;

if ($maq_type == 1){$maq = "yes"}
if ($maq_type == 2){$maq = "no"}


if ($maq eq "yes")
{

	print "\nMAQ ANALYSIS WILL BE INCLUDED\n";

}

if ($maq eq "no")
{

	print "\nMAQ ANALYSIS WILL NOT BE INCLUDED\n";

}
}
}

#################
#PINDEL analysis#
#################

if ($data eq "PE")
{

print "\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " Perform structural variant analysis?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";

print "   Enter 1 for YES\n";
print "   Enter 2 for NO\n\n";


$sv_type = <STDIN>;
chomp $sv_type;

if ($sv_type == 1){$sv = "yes"}
if ($sv_type == 2){$sv = "no"}


if ($sv eq "yes")
{

	print "\nSTUCTURAL VARIANT ANALYSIS WILL BE INCLUDED\n\n";

print "\nPlease input the expected average insert size for your library (eg 200):      ";
$insert = <STDIN>;
chomp $insert;
print "\n\n";

}

if ($sv eq "no")
{

	print "\nSTUCTURAL VARIANT ANALYSIS WILL NOT BE INCLUDED\n";

}



if ($sv eq "yes")
{

if ($ref_ans eq "partial")
{

	print "You have indicated that the reference is a partial genome sequence\n";

print "\nPlease enter the chromosome number you are interested in:    ";
$chromo = <STDIN>;
chomp $chromo;


}

}





$command =  "clear";
system("$command");

print color "bold white";

print "================================================\n";
print "SUMMARY - PLEASE CHECK THESE DETAILS CAREFULLY!!\n";
print "================================================\n\n";


print color "bold green";

print "YOUR DETAILS :\n\n";

print color "bold white";

print "Name of this analysis:\t$name\n\n"; 


print "Your email address:\t$to\n\n";

print color "reset";


print "Please press enter to proceed (or 'Q' to quit):      ";
$proceed = <STDIN>;
chomp $proceed; 
 
if (lc $proceed eq "q")
{
 exit;
} 


$command =  "clear";
system("$command");

print color "bold white";

print "================================================\n";
print "SUMMARY - PLEASE CHECK THESE DETAILS CAREFULLY!!\n";
print "================================================\n\n";


print color "bold green";

print "SETTINGS :\n\n";

print color "bold white";


if  ($data eq "SE")
{
print "Single-end analysis\n\n";
}

if  ($data eq "PE")
{
print "Paired-end analysis\n\n";
}



if ($dup eq "yes")
{

	print "Duplicates will be removed\n\n";

}

if ($dup eq "no")
{

	print "Duplicates will not be removed\n\n";

}






if ($maq eq "yes")
{

	print "Include MAQ analysis\n\n";

}

if ($maq eq "no")
{

	print "MAQ analysis not included\n\n";

}

if ($var eq "yes")
{

	print "Variants will be annotated using the ensembl database\n\n";

}

if ($var eq "no")
{

	print "\nVariants will not be annotated\n\n";

}


if ($sv eq "yes")
{

	print "Structural variant analysis included - mean insert size $insert bp ";

}



if ($sv eq "yes")
{
if ($ref_ans eq "partial")
{

print "(Chromosome $chromo)\n\n"
}
}

if ($sv eq "yes")
{
if ($ref_ans eq "whole")
{

print "\n\n"
}
}


if ($sv eq "no")
{

	print "Structural variant analysis not included\n\n";

}






print "Memory setting: $mem\n\n";

print color "reset";

print "Please press enter to proceed (or 'Q' to quit):      ";
$proceed = <STDIN>;
chomp $proceed; 
 
if (lc $proceed eq "q")
{
 exit;
} 


$command =  "clear";
system("$command");

print color "bold white";

print "================================================\n";
print "SUMMARY - PLEASE CHECK THESE DETAILS CAREFULLY!!\n";
print "================================================\n\n";


print color "bold green";

print "DATA FILES :\n\n";

print color "bold white";

if ($target_ans eq "whole")
{

	print "whole genome reference\n\n";

}

if ($ref_ans eq "partial")
{

	print "partial genome reference\n\n";

}


print "Reference file:\t\t$ref\n\n";


if  ($data eq "SE")
{
print "Sequence reads file:\t$reads\n\n";
}
if  ($data eq "PE")
{
print "Sequence reads file 1:\t$reads\n\n";
print "Sequence reads file 2:\t$reads2\n\n";
}

 
print color "reset";


print "Please press enter to proceed (or 'Q' to quit):      ";

print color "bold red";
print "\nNOTE - RUN TIME MAY BE SEVERAL HOURS."; 

$proceed = <STDIN>;
chomp $proceed; 
 
if (lc $proceed eq "q")
{
 exit;
} 
      
}


# END OF THE ILLUMINA SCRIPT

print color "reset";

###############################
#                             #
#   START OF THE 454 SCRIPT   #
#                             #
###############################

if ($platform eq "454")

{

print "\n\n454 analysis is currently only for SE datasets\n\n";



##################################
# Define data files              #
##################################

print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~";
print "\nDefine data files and memory setting";
print "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

print "\nPlease input the name of your reference fasta file (eg ref):      ";
$ref = <STDIN>;
chomp $ref;


#####################################
# Add fasta if user hasn't added it #
#####################################

if (index($ref,".fasta") == -1 )
{
	$ref = $ref.".fasta";
}


print "\nPlease input the name of your sff sequence file (eg reads):      ";
$reads3 = <STDIN>;
chomp $reads3;

#####################################
#  Add sff if user hasn't added it  #
#####################################

if (index($reads3,".sff") == -1 )
{
	$reads3 = $reads3.".sff";
}




print "\nPlease enter the memory setting for this analysis (default is -Xmx2g):      ";
$mem = <STDIN>;
chomp $mem;
 
if ($mem eq "")
{
 $mem = "-Xmx2g"
}




########################
# Check if file exists #
########################

if (! -e "$ref")
{ 
	print "\nFile $ref does not exist\n\n";
	exit;
} 
if (! -e "$reads3")
{ 
	print "\nFile $reads3 does not exist\n\n";
	exit;
} 



print "\n\n\n\n\n================================================\n";
print           "SUMMARY - PLEASE CHECK THESE DETAILS CAREFULLY!!\n";
print           "================================================\n\n";



print "The platform is 454\n\n";

print "Name of this analysis:\t$name\n\n"; 

print "Reference file:\t\t$ref\n\n";

print "Sequence reads file:\t$reads3\n\n";

print "Memory setting:\t\t$mem\n\n"; 

print "Your email address:\t$to\n\n"; 



print "*************************************\n";
print "NOTE Run time may be several hours...\n";
print "*************************************\n\n";


print "Please press enter to proceed (or 'Q' to quit):      ";
$proceed = <STDIN>;
chomp $proceed; 
 
if (lc $proceed eq "q")
{
 exit;
} 


print "\n\nConvert data from sff to fastq format\n\n";



$command =  "/opt/sff2fastq/sff2fastq $reads3 -o out.fastq";

print("\n$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print "\nDONE!\n";

}




if ($maq eq "yes")

{

print "\nRunning MAQ analysis...\n\n";


$command =  "maq.pl easyrun -d MAQ_$name $ref $reads $reads2";

print("\n$command\n");
system("$command");

}




if ($file eq "fastq")
{




##################################
# 1 Index the reference sequence #
##################################


if ($ref_ans eq "partial")

{

$command =  "bwa64 index -a is $ref";

print("\n$command\n");
system("$command");

}

if ($ref_ans eq "whole")

{

$command =  "bwa64 index -a bwtsw $ref";

print("\n$command\n");
system("$command");

}




#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);



print "\n============================================";
print "\n1/17 The reference sequence has been indexed";
print "\n============================================\n\n";



######################################################
# 2 Finds suffix array (SA) coordinates of good hits for each read #
######################################################

#
#SINGLE END DATA !
#

if  ($data eq "SE")
{

$command =  "bwa64 aln $ref $reads > aln_sa.sai";

print("$command\n");
system("$command");


print	"\n====================================================================";
print "\n2/17 Suffix array coordinates have been found for good sequence hits";
print "\n====================================================================\n\n";

}

#
#PAIRED END DATA !
#

if  ($data eq "PE")
{

$command =  "bwa64 aln $ref $reads > aln_sa1.sai";

print("$command\n");
system("$command");


print "\n====================================================================================";
print "\n2a/17 Suffix array coordinates have been found for good sequence hits for seq file 1";
print "\n====================================================================================\n\n";

$command =  "bwa64 aln $ref $reads2 > aln_sa2.sai";

print("$command\n");
system("$command");

print "\n====================================================================================";
print "\n2b/17 Suffix array coordinates have been found for good sequence hits for seq file 2";
print "\n====================================================================================\n\n";

}




#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);



##############################################################
# 3 Convert SA coordinates into chromosomal coordinates 
##############################################################

if  ($data eq "SE")

{

$command =  "bwa64 samse $ref aln_sa.sai $reads > aligned.sam";

print("$command\n");
system("$command");

print "\n==============================================================================";
print "\n3/17 Suffix array coordinates have been converted into chromosomal coordinates";
print "\n==============================================================================\n\n";

}

if  ($data eq "PE")

{

$command =  "bwa64 sampe $ref aln_sa1.sai aln_sa2.sai $reads $reads2 > aligned.sam";

print("$command\n");
system("$command");

print "\n==============================================================================";
print "\n3/17 Suffix array coordinates have been converted into chromosomal coordinates";
print "\n==============================================================================\n\n";

}

if  ($platform eq "454")

{

$command =  "bwa64 bwasw $ref $reads4 > aligned.sam";

print("$command\n");
system("$command");

print "\n========================================================================";
print "\n3/17 454 Suffix array coords have been converted into chromosomal coords";
print "\n========================================================================\n\n";

}








#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


}


###########################################
# 4 Create the fasta sequence dictionary file          
###########################################

$command =  "/usr/bin/java16 -jar $mem /opt/picard/CreateSequenceDictionary.jar R=$ref O=$ref.dict";

print("$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print	"\n=================================================================";
print "\n4/17 A dictionary file has been created for the reference sequence";
print "\n=================================================================\n\n";




########################################################################################
# 5 Create the fasta index file       
########################################################################################


$command =  "/opt/samtools/samtools faidx $ref";

print("$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print	"\n=============================================================";
print "\n5/17 An index file has been created for the reference sequence";
print "\n=============================================================\n\n";



if ($file eq "fastq")
{

######################################
# 6 Sort the Sam file generated by BWA   
######################################

$command =  "/usr/bin/java16 $mem -jar /opt/picard/SortSam.jar I=aligned.sam O=aligned_sorted.sam SO=coordinate VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print	"\n=====================================================";
print "\n6/17 The aligned reads in SAM format have been sorted";
print "\n=====================================================\n\n";

########################################################################################
# 7 Convert the SAM file from BWA to a BAM file       
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/picard/SamFormatConverter.jar I=aligned_sorted.sam O=aligned_sorted.bam VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");

#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);

print	"\n=================================================";
print "\n7/17 The SAM file has been coverted to BAM format";
print "\n=================================================\n\n\n";

}


if ($file eq "bam")

{


#renaming bam file


$command = "cp $bam_file aligned_sorted.bam";
system("$command");


}



















###############################################################
# 7 b Duplicate removal
###############################################################


if  ($data eq "PE")
{

if ($dup eq "yes")
{

open(STDERR, "| tee Duplicate_info.rtf");

$command =  "/opt/samtools/samtools rmdup aligned_sorted.bam duplicates_removed.bam";

print("$command\n");
system("$command");


print	"\n========================";
print 	"\n7b/17 DUPLICATES REMOVED";
print 	"\n========================\n\n\n";

print	"\nre-naming original file...\n\n";

$command =  "mv aligned_sorted.bam aligned_plus_dup.bam";

print("$command\n");
system("$command");

print	"\nre-naming new file...\n\n";


$command =  "mv duplicates_removed.bam aligned_sorted.bam ";

print("$command\n\n");
system("$command");


close(STDERR);


}
}
	
if  ($data eq "SE")
{

if ($dup eq "yes")
{

open(STDERR, "| tee Duplicate_info.rtf");

$command =  "/opt/samtools/samtools rmdup -s aligned_sorted.bam duplicates_removed.bam";

print("$command\n");
system("$command");


print	"\n========================";
print 	"\n7b/17 DUPLICATES REMOVED";
print 	"\n========================\n\n\n";

print	"\nre-naming original file...\n\n";

$command =  "mv aligned_sorted.bam aligned_plus_dup.bam";

print("$command\n");
system("$command");

print	"\nre-naming new file...\n\n";


$command =  "mv duplicates_removed.bam aligned_sorted.bam ";

print("$command\n\n");
system("$command");


close(STDERR);


}
}





########################################################################################
# 8 Create an index for the Bam file    
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/picard/BuildBamIndex.jar I=aligned_sorted.bam VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");





if ($target_ans eq "yes")

{

$command =  "/opt/samtools/samtools view aligned_sorted.bam $region -b -o aligned_sorted_2.bam";

print("$command\n\n");
system("$command");


$command =  "mv aligned_sorted.bam all_reads_aligned_sorted.bam ";

print("$command\n\n");
system("$command");

$command =  "mv aligned_sorted_2.bam aligned_sorted.bam ";

print("$command\n\n");
system("$command");

$command =  "/usr/bin/java16 $mem -jar /opt/picard/BuildBamIndex.jar I=aligned_sorted.bam VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");

}






#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);

print	"\n====================================================";
print "\n8/17 An index file has been created for the BAM file";
print "\n====================================================\n\n\n";


if ($file eq "fastq")
{

########################################################################################
# 9 Count covariates    
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/gatk_v1/GenomeAnalysisTK.jar -R $ref -I aligned_sorted.bam -T CountCovariates -cov ReadGroupCovariate -cov QualityScoreCovariate -cov CycleCovariate -cov DinucCovariate -recalFile reads.csv --default_platform $platform";

print("$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);



print	"\n=================================";
print "\n9/17 Covariates have been counted";
print "\n=================================\n\n\n";



########################################################################################
# 10 TableRecalibration     
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/newgatk/GenomeAnalysisTK.jar -R $ref -I aligned_sorted.bam -T TableRecalibration -outputBam recal.bam -recalFile reads.csv --default_platform $platform";

print("$command\n");
system("$command");

#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print	"\n===========================================================";
print "\n10/17 Quality scores in the BAM file have been recalculated";
print "\n===========================================================\n\n";

########################################################################################
# 11 Index the new BAM file    
########################################################################################

$command =  "/opt/samtools/samtools index recal.bam";

print("$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print "\n=======================================";
print "\n11/17 The new BAM file has been indexed";
print "\n=======================================\n\n";

########################################################################################
# 12 Create intervals    
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/newgatk/GenomeAnalysisTK.jar -T RealignerTargetCreator -I recal.bam -R $ref -o forRealigner.intervals";

print("$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print "\n==========================================================================";
print "\n12/17 Small intervals possibly in need of realignment have been identified";
print "\n==========================================================================\n\n";

########################################################################################
# 13 Realigning    
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/newgatk/GenomeAnalysisTK.jar -I recal.bam -R $ref -T IndelRealigner -targetIntervals forRealigner.intervals --output cleaned.bam";

print("$command\n");
system("$command");

#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print "\n===================================";
print "\n13/17 The reads have been realigned";
print "\n===================================\n\n";

########################################################################################
# 14 Sort the Cleaned Bam file    
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/picard/SortSam.jar I=cleaned.bam O=cleaned_sorted.bam SO=coordinate VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");

#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print "\n====================================================";
print "\n14/17 The realigned cleaned Bam file has been sorted ";
print "\n====================================================\n\n";

########################################################################################
# 15 Create an index for the cleaned and sorted bam file   
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/picard/BuildBamIndex.jar I=cleaned_sorted.bam VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print "\n====================================================";
print "\n15/17 The cleaned sorted Bam file has been indexed ";
print "\n====================================================\n\n";

}

if ($file eq "bam")


{

$command =  "mv aligned_sorted.bam cleaned_sorted.bam ";

print("$command\n\n");
system("$command");

$command =  "mv aligned_sorted.bai cleaned_sorted.bai ";

print("$command\n\n");
system("$command");


}



if ($platform eq "Illumina")
{

########################################################################################
# 16 Indel genotyper   
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/newgatk/GenomeAnalysisTK.jar -T IndelGenotyperV2 -R $ref -I cleaned_sorted.bam -O indels.raw.bed -o detailed.output.bed --verbose";

print("$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print "\n================================";
print "\n16/17 Indel calls have been made";
print "\n================================\n\n";

}





if ($platform eq "454")
{

########################################################################################
# 16 Indel genotyper   
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/newgatk/GenomeAnalysisTK.jar -T IndelGenotyperV2 -R $ref -I cleaned_sorted.bam -O indels.raw.bed -o detailed.output.bed --verbose --window_size 800";

print("$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);


print "\n================================";
print "\n16/17 Indel calls have been made";
print "\n================================\n\n";

}






if ($platform eq "Illumina")
{

########################################################################################
# 17 Making SNP calls   
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/newgatk/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper -I cleaned_sorted.bam -stand_emit_conf 10.0 -varout snps.raw.vcf -stand_call_conf 50.0 --platform SOLEXA";

print("$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);



print "\n==============================";
print "\n17/17 SNP Calls have been made";
print "\n==============================\n";

}





if ($platform eq "454")
{

########################################################################################
# 17 Making SNP calls   
########################################################################################

$command =  "/usr/bin/java16 $mem -jar /opt/newgatk/GenomeAnalysisTK.jar -R $ref -T UnifiedGenotyper -I cleaned_sorted.bam -varout snps.raw.vcf -stand_emit_conf 10.0 -stand_call_conf 50.0 --platform ROCHE454";

print("$command\n");
system("$command");


#~~~~~~~~~~~~~~~#
#ERROR CHECK v2!#
#~~~~~~~~~~~~~~~#

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print "\nAN ERROR HAS BEEN FOUND IN THE LOG!\n   PERL SCRIPT WILL TERMINATE!\n\n" if /\bERROR\b/i;
	
if (/\bERROR\b/i)
{
open (MAIL, "|/usr/sbin/sendmail -t");
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: ERROR! $subject ($name)\n\n";
print MAIL "AN ERROR HAS OCCURED IN NGS ANALYSIS! ($name)\n\nPlease check the log.rtf for details.\n";

close(MAIL);
}

	exit if /\bERROR\b/i;
}
close(LOG);



print "\n==============================";
print "\n17/17 SNP Calls have been made";
print "\n==============================\n";

}










########################################
# Structural variant analysis - Pindel #
########################################


if ($sv eq "yes")

{

if ($ref_ans eq "partial")
{

print "\n\nStarting script to find structural variants...";

print "\n\nSorting Bam file by query name...";

$command =  "/usr/bin/java16 $mem -jar /opt/picard/SortSam.jar I=cleaned_sorted.bam O=queryname_sorted.bam SO=queryname VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");


print "\n\nBam sorted, now making input files for pindel\n\n";


$command =  "perl /opt/pindel/bam2pindel.pl -i queryname_sorted.bam -o sv -s $name -pi $insert";


print("$command\n");
system("$command");

print "\n\nRunning Pindel\n\n";

$command =  "/opt/pindel/pindel022 -f $ref -p sv_chr$chromo.txt -c chr$chromo -o pindel";


print("$command\n");
system("$command");

}



if ($ref_ans eq "whole")
{

print "\n\nStarting script to find structural variants...";

print "\n\nSorting Bam file by query name...";

$command =  "/usr/bin/java16 $mem -jar /opt/picard/SortSam.jar I=cleaned_sorted.bam O=queryname_sorted.bam SO=queryname VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");


print "\n\nBam sorted, now making input files for pindel\n\n";


$command =  "perl /opt/pindel/bam2pindel.pl -i queryname_sorted.bam -o sv -s $name -pi $insert";

print("$command\n");
system("$command");


print "\n\nMerging all Pindel files\n\n";


$command =  "cat sv_chr1.txt sv_chr2.txt sv_chr3.txt sv_chr4.txt sv_chr5.txt sv_chr6.txt sv_chr7.txt sv_chr8.txt sv_chr9.txt sv_chr10.txt sv_chr11.txt sv_chr12.txt sv_chr13.txt sv_chr14.txt sv_chr15.txt sv_chr16.txt sv_chr17.txt sv_chr18.txt sv_chr19.txt sv_chr20.txt sv_chr21.txt sv_chr22.txt sv_chr23.txt sv_chr24.txt sv_chr25.txt sv_chr26.txt sv_chr27.txt sv_chr28.txt sv_chr29.txt sv_chr30.txt sv_chr31.txt sv_chr32.txt sv_chr33.txt sv_chr34.txt sv_chr35.txt sv_chr36.txt sv_chr37.txt sv_chr38.txt sv_chrX.txt sv_chrM.txt sv_chrUn.txt>all.txt";

print("$command\n");
system("$command");




print "\n\nRunning Pindel\n\n";

$command =  "/opt/pindel/pindel022 -f $ref -p all.txt -c ALL -o pindel";


print("$command\n");
system("$command");

}

}


###############################
# RUN PINDEL TO VCF CONVERTER #
###############################


if ($sv eq "yes")

{


$command =  "/opt/pindel/pindel2vcf -p pindel_BP  -r $ref -R canfam2 -d 2006";
print("$command\n");
system("$command");

$command =  "/opt/pindel/pindel2vcf -p pindel_D  -r $ref -R canfam2 -d 2006";
print("$command\n");
system("$command");

$command =  "/opt/pindel/pindel2vcf -p pindel_INV  -r $ref -R canfam2 -d 2006";
print("$command\n");
system("$command");

$command =  "/opt/pindel/pindel2vcf -p pindel_LI  -r $ref -R canfam2 -d 2006";
print("$command\n");
system("$command");

$command =  "/opt/pindel/pindel2vcf -p pindel_SI  -r $ref -R canfam2 -d 2006";
print("$command\n");
system("$command");

$command =  "/opt/pindel/pindel2vcf -p pindel_TD  -r $ref -R canfam2 -d 2006";
print("$command\n");
system("$command");

print "\n\nVCF FILES MADE FOR PINDEL RESULTS\n\n";

}


####################################################################
# Checking depth of coverage, average insert size and GC bias
####################################################################


print "\n Checking depth of coverage....\n\n";

$command = "/usr/bin/java16 $mem -jar /opt/newgatk/GenomeAnalysisTK.jar -R $ref -I cleaned_sorted.bam -o depth -T DepthOfCoverage -omitBaseOutput";

print("$command\n");
system("$command");


####################################################

print "\n Calculating GC Bias....\n\n";

$command = "/usr/bin/java16 $mem -jar /opt/picard/CollectGcBiasMetrics.jar REFERENCE_SEQUENCE=$ref INPUT=cleaned_sorted.bam O=out1.junk CHART_OUTPUT=GC_bias.pdf VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");



####################################################

if ($dup eq "yes")
{

print "\n Collecting alignment summary metrics....\n\n";

$command = "/usr/bin/java16 $mem -jar /opt/picard/CollectAlignmentSummaryMetrics.jar I=aligned_plus_dup.bam O=Alignment_Summary.xls R=$ref VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");

}

if ($dup eq "no")
{

print "\n Collecting alignment summary metrics....\n\n";

$command = "/usr/bin/java16 $mem -jar /opt/picard/CollectAlignmentSummaryMetrics.jar I=aligned_sorted.bam O=Alignment_Summary.xls R=$ref VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");

}


###################################################


print "\n Plotting insert size histogram (for PE data only)....\n\n";


if ($dup eq "no")
{

$command = "/usr/bin/java16 $mem -jar /opt/picard/CollectInsertSizeMetrics.jar INPUT=aligned_sorted.bam O=out2.junk HISTOGRAM_FILE=INSERT_SIZE.pdf VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");

}


if ($dup eq "yes")
{

$command = "/usr/bin/java16 $mem -jar /opt/picard/CollectInsertSizeMetrics.jar INPUT=aligned_plus_dup.bam O=out2.junk HISTOGRAM_FILE=INSERT_SIZE.pdf MINIMUM_PCT=0.05 VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");

}


####################################################


######################################
# Create the results folder          #
######################################


print "\nCreating new directory for results....\n\n";


$command = "mkdir results_$name";

print("$command\n");
system("$command");

if ($sv eq "yes")

{

$command = "mkdir results_$name/Raw_SV_data";

print("$command\n");
system("$command");

}


print "\nSorting results into new directory....\n\n";


if ($dup eq "yes")
{


$command = "rm aligned_sorted.bam";
print("$command\n");
system("$command");


$command = "mv aligned_plus_dup.bam aligned_sorted.bam";
print("$command\n");
system("$command");

#rebuild bam index

$command =  "/usr/bin/java16 $mem -jar /opt/picard/BuildBamIndex.jar I=aligned_sorted.bam VALIDATION_STRINGENCY=LENIENT";

print("$command\n");
system("$command");


}


&move_to_results_folder ("cleaned_sorted.bam");
&move_to_results_folder ("cleaned_sorted.bai");
&move_to_results_folder ("snps.raw.vcf");
&move_to_results_folder ("aligned_sorted.bai");
&move_to_results_folder ("aligned_sorted.bam");
&move_to_results_folder ("indels.raw.bed");
&move_to_results_folder ("depth");
&move_to_results_folder ("depth.sample_summary");
&move_to_results_folder ("GC_bias.pdf");
&move_to_results_folder ("INSERT_SIZE.pdf");
&move_to_results_folder ("Duplicate_info.rtf");
&move_to_results_folder ("Alignment_Summary.xls");


$command = "mv results_$name/aligned_sorted.bam results_$name/raw_align_$name.bam";
system("$command");

$command = "mv results_$name/aligned_sorted.bai results_$name/raw_align_$name.bai";
system("$command");

$command = "mv results_$name/cleaned_sorted.bam results_$name/best_align_$name.bam";
system("$command");

$command = "mv results_$name/cleaned_sorted.bai results_$name/best_align_$name.bai";
system("$command");

$command = "mv results_$name/indels.raw.bed results_$name/INDELS_$name.bed";
system("$command");

$command = "mv results_$name/snps.raw.vcf results_$name/SNPS_$name.vcf";
system("$command");

$command = "mv results_$name/depth results_$name/depth_bp.xls";
system("$command");

$command = "mv results_$name/depth.sample_summary results_$name/depth_summary.xls";
system("$command");



if ($sv eq "yes")

{

&delete_file ("queryname_sorted.bam");

&move_to_results_folder ("pindel_BP");
&move_to_results_folder ("pindel_D");
&move_to_results_folder ("pindel_INV");
&move_to_results_folder ("pindel_LI");
&move_to_results_folder ("pindel_SI");
&move_to_results_folder ("pindel_TD");


$command = "mv results_$name/pindel_BP results_$name/Raw_SV_data/SV_Break_points.txt";
system("$command");

$command = "mv results_$name/pindel_D results_$name/Raw_SV_data/SV_Deletions.txt";
system("$command");

$command = "mv results_$name/pindel_INV results_$name/Raw_SV_data/SV_Inversions.txt";
system("$command");

$command = "mv results_$name/pindel_LI results_$name/Raw_SV_data/SV_Long_insertions.txt";
system("$command");

$command = "mv results_$name/pindel_SI results_$name/Raw_SV_data/SV_Non_template_seq_in_deleletion.txt";
system("$command");

$command = "mv results_$name/pindel_TD results_$name/Raw_SV_data/SV_Tandom_Dup.txt";
system("$command");


&move_to_results_folder ("pindel_BP.vcf");
&move_to_results_folder ("pindel_D.vcf");
&move_to_results_folder ("pindel_INV.vcf");
&move_to_results_folder ("pindel_LI.vcf");
&move_to_results_folder ("pindel_SI.vcf");
&move_to_results_folder ("pindel_TD.vcf");


$command = "mv results_$name/pindel_BP.vcf results_$name/SV_Break_points.vcf";
system("$command");

$command = "mv results_$name/pindel_D.vcf results_$name/SV_Deletions.vcf";
system("$command");

$command = "mv results_$name/pindel_INV.vcf results_$name/SV_Inversions.vcf";
system("$command");

$command = "mv results_$name/pindel_LI.vcf results_$name/SV_Long_insertions.vcf";
system("$command");

$command = "mv results_$name/pindel_SI.vcf results_$name/SV_Non_template_seq_in_deleletion.vcf";
system("$command");

$command = "mv results_$name/pindel_TD.vcf results_$name/SV_Tandom_Dup.vcf";
system("$command");




if ($ref_ans eq "partial")

{
&delete_file ("sv_chr$chromo.txt");
}

if ($ref_ans eq "whole")
{
&delete_file ("all.txt");
&delete_file ("sv_chr1.txt");
&delete_file ("sv_chr2.txt");
&delete_file ("sv_chr3.txt");
&delete_file ("sv_chr4.txt");
&delete_file ("sv_chr5.txt");
&delete_file ("sv_chr6.txt");
&delete_file ("sv_chr7.txt");
&delete_file ("sv_chr8.txt");
&delete_file ("sv_chr9.txt");
&delete_file ("sv_chr10.txt");
&delete_file ("sv_chr11.txt");
&delete_file ("sv_chr12.txt");
&delete_file ("sv_chr13.txt");
&delete_file ("sv_chr14.txt");
&delete_file ("sv_chr15.txt");
&delete_file ("sv_chr16.txt");
&delete_file ("sv_chr17.txt");
&delete_file ("sv_chr18.txt");
&delete_file ("sv_chr19.txt");
&delete_file ("sv_chr20.txt");
&delete_file ("sv_chr21.txt");
&delete_file ("sv_chr22.txt");
&delete_file ("sv_chr23.txt");
&delete_file ("sv_chr24.txt");
&delete_file ("sv_chr25.txt");
&delete_file ("sv_chr26.txt");
&delete_file ("sv_chr27.txt");
&delete_file ("sv_chr28.txt");
&delete_file ("sv_chr29.txt");
&delete_file ("sv_chr30.txt");
&delete_file ("sv_chr31.txt");
&delete_file ("sv_chr32.txt");
&delete_file ("sv_chr33.txt");
&delete_file ("sv_chr34.txt");
&delete_file ("sv_chr35.txt");
&delete_file ("sv_chr36.txt");
&delete_file ("sv_chr37.txt");
&delete_file ("sv_chr38.txt");
&delete_file ("sv_chrX.txt");
&delete_file ("sv_chrM.txt");
&delete_file ("sv_chrUn.txt");

}

}

if ($target_ans eq "yes")

{
&delete_file ("all_reads_aligned_sorted.bam");
}

##########
# README #
##########

open (MYFILE, '>>README.rtf'); 

print	MYFILE	"========================\n";
print	MYFILE	"Summary of results files\n";
print	MYFILE	"========================\n\n";

print	MYFILE	"PDF FILES\n\n";

print	MYFILE	"GC_ Bias.pfd	Histogram of GC content of aligned reads\n";
print	MYFILE	"INSERT_SIZE.pdf	Insert size histogram\n\n";

print	MYFILE	"ALIGNMENT FILES\n\n";

print	MYFILE	"raw_align.bam	Raw alignments to the reference\n";
print	MYFILE	"best_align.bam	Best alignments to the reference after processing by GATK\n\n";

print 	MYFILE	"SNP AND INDEL CALLS\n\n";

print	MYFILE	"INDELS.bed		List of InDels\n";
print	MYFILE	"SNPS.vcf		List of SNPs\n";
print	MYFILE	"annotated_indels.xls	List of InDels annotated using the ensembl database\n";
print	MYFILE	"annotated_snp.xls	List of SNPs annotated using the ensembl database\n\n";


print 	MYFILE	"STRUCTURAL VARIANT FILES\n\n";

print 	MYFILE	"SV_Break_points.txt\n";
print 	MYFILE	"SV_Deletions.txt\n";
print 	MYFILE	"SV_Inversions.txt\n";
print 	MYFILE	"SV_Long_insertions.txt\n";
print 	MYFILE	"SV_Non_template_seq_in_deleletion.txt\n";
print 	MYFILE	"SV_Tandom_Dup.txt\n\n";

print	MYFILE	"INFORMATION FILES\n\n";	

print	MYFILE	"Depth_summary	\tSummary of reads depth across the target region\n";
print	MYFILE	"Log.rtf	\t\tRun log for the NGS pipeline\n";
print	MYFILE	"Duplicate_info.rtf	fraction of PCR duplicates reads in the dataset\n\n";

print	MYFILE	"NOTES\n\n"; 
print	MYFILE	".bam files must have and associated .bai index file for loading into IGV\n";
print	MYFILE	"Alignment_Summary.xls can be used to calculate success of target enrichment\n\n";

close (MYFILE);

&move_to_results_folder ("README.rtf");

#########################################################################################


&delete_file ("$ref.dict");
&delete_file ("aligned.sam");
&delete_file ("aln_sa.sai");
&delete_file ("$ref.amb");
&delete_file ("$ref.ann");
&delete_file ("$ref.bwt");
&delete_file ("$ref.dict");
&delete_file ("$ref.fai");
&delete_file ("$ref.pac");
&delete_file ("$ref.rbwt");
&delete_file ("$ref.rpac");
&delete_file ("$ref.rsa");
&delete_file ("aligned_sorted.sam");
&delete_file ("$ref.sa");
&delete_file ("reads.csv");
&delete_file ("recal.bam");
&delete_file ("recal.bam.bai");
&delete_file ("forRealigner.intervals");
&delete_file ("cleaned.bam");
&delete_file ("detailed.output.bed");
&delete_file ("aln_sa2.sai");
&delete_file ("aln_sa1.sai");
&delete_file ("depth.sample_cumulative_coverage_counts");
&delete_file ("depth.sample_cumulative_coverage_proportions");
&delete_file ("depth.sample_statistics");
&delete_file ("*.dict");
&delete_file ("*.junk");








###########################
#                         # 
#   ANNOTATING VARIANTS   #
#                         #
###########################


if ($var eq "yes")
{

print "\n\nANNOTATING VARIANTS...\n\n";

$command = ". ~/.bashrc";

print("$command\n\n");
system("$command");



$command = "perl /opt/ensembl/variant_effect_predictor.pl --species canis_familiaris --input results_$name/SNPS_$name.vcf --format VCF --output_file annotated_snps_$name.xls -v";

print("$command\n\n");
system("$command");

}


if ($var eq "yes")
{

$command = ". ~/.bashrc";

print("$command\n\n");
system("$command");



$command = "perl /opt/ensembl/variant_effect_predictor.pl --species canis_familiaris --input results_$name/INDELS_$name.bed --format pileup --output_file annotated_indels_$name.xls -v";

print("$command\n\n");
system("$command");

print "\n\nVariant annotation complete\n\n";

}


if ($var eq "yes")
{

&move_to_results_folder ("annotated_snps_$name.xls");
&move_to_results_folder ("annotated_indels_$name.xls");

}

#RUN TIMER END

my $end_run = time();
my $run_time = $end_run - our $start_run;

open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $to\n";
print MAIL "From: $from\n";
print MAIL "Subject: $subject COMPLETE ($name)\n\n";
## Mail Body
print MAIL "Your next generation sequence analysis ($name) is complete\n\n";
print MAIL "For pipeline details see log.rtf\n\n";
print MAIL "For file information see README.rtf\n\n";
print MAIL "Run time : $run_time seconds\n";
printf MAIL"%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];
close(MAIL);


print "\nANALYSIS COMPLETE! YOU HAVE BEEN NOTIFIED BY EMAIL\n\n";


print "Run time : $run_time seconds\n";
printf "%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];


print "XXXXXXXXXXXXXXXXXXXXXXX\n";
print "  ERRORS AND WARNINGS  \n";
print "XXXXXXXXXXXXXXXXXXXXXXX\n\n";

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print if /\bwarning\b/i;
}
close(LOG);

print "\n";

open(LOG,"log.rtf") or die "Unable to open logfile:$!\n";
while(<LOG>){
	print if /\berror\b/i;
}
close(LOG);



print "\nPlease check log.rtf for full details\n\n";



#############################
# Turn off logging          #
#############################

close(STDOUT);

&move_log ("log.rtf");



#############################################
#                                           #
# Subroutine to move file to results folder #
#                                           #
#############################################

sub move_to_results_folder
{

	my $suffix = "";	

	$suffix = $_[0];
	$command = "mv  $suffix results_$name/$suffix";
	print("$command\n");
	system("$command");

}

#############################################
#                                           #
# Subroutine to delete files              r #
#                                           #
#############################################

sub delete_file
{

	my $suffix = "";	

	$suffix = $_[0];
	$command = "rm  $suffix";
	print("$command\n");
	system("$command");

}

#############################################
#                                           #
# Sub to move log			    #
#                                           #
#############################################

sub move_log
{

	my $suffix = "";	

	$suffix = $_[0];
	$command = "mv  $suffix results_$name/$suffix";
	system("$command");

}

