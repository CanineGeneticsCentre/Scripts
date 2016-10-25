#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	BAM MERGE       						                            #     
#									                                    #
#	THIS PERL SCRIPT WILL RUN PICARD MergeSamFiles   	                #
#									                                    #
#########################################################################

#############################
# Mike Boursnell July 2012  #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# oliver.forman@aht.org.uk  #
#############################

use strict;
use Getopt::Std ;
use File::Basename ;

my $list_count				= 0;
my $no_of_files				= 0;

my $list_file				= "";
my $command					= "";
my $bam_file				= "";
my $merge_string		    = "";
my $answer					= "";
my $tempdir					= "";
my $temp_dir_string			= "";
my $make_index_files		= "no";
my $output_file				= "merged.bam";

my @bam_file_array			= ();


 


 
 
print "\n\n";
print "##############################################\n";
print "#  PERL script to run Picard MergeSamFiles   #\n";
print "##############################################\n\n";

################################
# Make the temp java directory #
################################
$tempdir = "$ENV{HOME}/javatempdir";

if (-e $tempdir)
{
	print "\nJava temporary directory at $tempdir already exists, so no need to create it\n\n";
}
if (! -e $tempdir)
{
	unless(mkdir $tempdir) 
	{
		die "Unable to create directory $tempdir";
	}
}


until (-e "$list_file")
	{
		print "\nPlease input the name of your file with a list of file names of the .bam files:      ";
		$list_file = <STDIN>;
		chomp $list_file;

		if ($list_file eq "ls")
		{
			print "\nList of .txt files...\n\n";
			system ("ls *.txt");
			print "\n";
		}
		
		if ($list_file ne "ls")
		{
			if (! -e "$list_file"){print "\nFile doesn't exist. Try again...    \n";}
		}
	}

print "\nDo you want to make the index files from the BAM files?\n\n";

print "<1> YES\n";
print "<2> NO\n\n";

$answer=<STDIN>;
chomp $answer;

if ($answer eq "1") {$make_index_files = "yes"}
	
	
#############################################
# Make sure the list file is in Unix format #
#############################################

$command = "dos2unix $list_file";
print("\n$command\n");
system("$command");



####################################################
# Open the list file to get the list of file names #
####################################################
open (LIST, "$list_file") || die "Cannot open $list_file";
$list_count=1;

while ($bam_file = <LIST> ) 
{
	chomp $bam_file;
	
	#########################################################
	# Add file type suffix .bam if user hasn't added it     #
	#########################################################

	if (index($bam_file,".bam") == -1 ){$bam_file = $bam_file.".bam"}

	$bam_file_array[$list_count]=$bam_file;
	$list_count=$list_count + 1;
}

close LIST;

$no_of_files=$list_count - 1;

###################
# List file names #
###################
print "\n\nThere are $no_of_files BAM files in this file of file names.\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print "File $list_count	\t$bam_file_array[$list_count]\n";
}

########################
# Make index bai files #
########################

if ($make_index_files eq "yes")
{
	print "\n\nTCreating index files from BAM files...\n\n";

	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		print "File $list_count	\t$bam_file_array[$list_count]\n";
		$bam_file = $bam_file_array[$list_count];
		
		$command = "java -Xmx4g -Djava.io.tmpdir=$tempdir -jar /opt/picard/BuildBamIndex.jar INPUT=$bam_file VALIDATION_STRINGENCY=LENIENT";
		print "$command\n";
		system ("$command");
	}

}






################################################
# Make a string of the names for MergeSamFiles #
################################################
print "\n\nThere are $no_of_files VCF files in this file of file names.\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	$merge_string = $merge_string." INPUT=$bam_file_array[$list_count]";
}



$command = "java -Xmx4g -Djava.io.tmpdir=$tempdir -jar /opt/picard/MergeSamFiles.jar $merge_string OUTPUT=$output_file VALIDATION_STRINGENCY=LENIENT";
print "$command\n";
system ("$command");

print "\n\nFINISHED MERGING BAM FILES\n\n";
print "Output file: $output_file\n\n";
	