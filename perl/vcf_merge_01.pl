#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	VCF MERGE       						                            #     
#									                                    #
#	THIS PERL SCRIPT WILL RUN VCFTOOLS vcf-merge     	                #
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
my $vcf_file				= "";
my $vcf_merge_string		= "";

my @vcf_file_array			= ();

print "\n\n";
print "#############################################\n";
print "#  PERL script to run VCFTOOLS vcf-merge    #\n";
print "##############################################\n\n";

print "The input is a file with a list of the VCF file names.\n\n";
print "(Files must be compressed with bgzip and indexed with tabix)\n\n";

print "Name of the file:      ";
$list_file = <STDIN>;
chomp $list_file;

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

while ($vcf_file = <LIST> ) 
{
	chomp $vcf_file;
	
	#########################################################
	# Add file type suffix .vcf if user hasn't added it     #
	#########################################################

	if (index($vcf_file,".vcf") == -1 ){$vcf_file = $vcf_file.".vcf"}

	####################################
	# if not already gzipped do it now #
	####################################
	if (index($vcf_file,".gz") == -1 )
	{
		$command = "bgzip $vcf_file";
		print "$command\n";
		system ("$command");
		
		$command = "tabix $vcf_file";
		print "$command\n";
		system ("$command");
	}
	
	#########################################################
	# Add file type suffix .gz if user hasn't added it      #
	#########################################################

	if (index($vcf_file,".gz") == -1 ){$vcf_file = $vcf_file.".gz"}
	
	
	$vcf_file_array[$list_count]=$vcf_file;
	$list_count=$list_count + 1;
}

close LIST;

$no_of_files=$list_count - 1;

###################
# List file names #
###################
print "\n\nThere are $no_of_files VCF files in this file of file names.\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	print "File $list_count	\t$vcf_file_array[$list_count]\n";
}

############################################
# Make a string of the names for vcf-merge #
############################################
print "\n\nThere are $no_of_files VCF files in this file of file names.\n\n";

for ($list_count=1;$list_count<=$no_of_files;$list_count++)
{
	$vcf_merge_string = $vcf_merge_string.$vcf_file_array[$list_count]." ";
}

$command = "perl /opt/vcftools/perl/vcf-merge $vcf_merge_string > vcf_merge_out.vcf";
print "$command\n";
system ("$command");

print "\n\nFINISHED\n\n";
	