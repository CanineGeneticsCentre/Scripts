#!/usr/bin/perl
use warnings;
use strict;

#####################################################
#              check_bam_file                       #
#                                                   #
# This runs various checks on the BAM file          #
#                                                   #
#  QualityScoreDistribution                         #
#                                                   #
#####################################################

#############################
# Mike Boursnell Dec 2012   #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# mike.boursnell@aht.org.uk #
#############################

my $version					= "1";


my $bam_file				= "";
my $output_file				= "";
my $picard_qsd_pdf			= "";
my $picard_qsd_out			= "";
my $picard_casm_out			= "";
my $picard_cism_out			= "";
my $picard_cism_histogram	= "";
my $picard_mm				= "";
my $prefix					= "";

my $single_line				= "";
my $last_char				= "";
my $penult_char				= "";
my $answer					= "";
my $char					= "";
my $start_lines				= "";


print "\n\n";
print "####################\n";
print "# check_bam_file   #\n";
print "####################\n\n";

print "Version: $version\n\n";

print "This program runs various checks on your BAM file.\n\n";

print "  - picard/QualityScoreDistribution\n";

# Get file name #
until (-e $bam_file)
{
	print "\n\nName of input BAM file (type 'ls' to see a list of BAM files):  ";
	$bam_file = <STDIN>;
	chomp $bam_file;
	
	if ($bam_file ne "ls")
	{
		if (index($bam_file,".bam") == -1){$bam_file = $bam_file.".bam";}
	}
	if ($bam_file eq "ls")
	{
		print"\n";
		system ("ls *.bam");
	}
}

##########################################
# Make other file names for output files #
##########################################
$prefix = &get_prefix ($bam_file);

$picard_mm = $prefix;
$picard_qsd_pdf = $prefix."_QualityScoreDistribution.pdf";
$picard_qsd_out = $prefix."_QualityScoreDistribution.out";
$picard_cism_histogram = $prefix."_CollectInsertSizeMetrics_histogram.out";
$picard_casm_out = $prefix."_CollectAlignmentSummaryMetrics.out";
$picard_cism_out = $prefix."CollectInsertSizeMetrics.out";


print "\n\n";
###################################
# picard/QualityScoreDistribution #
###################################
print "#########################################\n";
print "# picard/QualityScoreDistribution       #\n";
print "#########################################\n\n";

&run_unix_command("java -Xmx4g -jar /opt/picard/QualityScoreDistribution.jar  I=$bam_file CHART_OUTPUT=$picard_qsd_pdf OUTPUT=$picard_qsd_out VALIDATION_STRINGENCY=SILENT");



#########################################
# picard/CollectAlignmentSummaryMetrics #
#########################################
print "#########################################\n";
print "# picard/CollectAlignmentSummaryMetrics #\n";
print "#########################################\n\n";
&run_unix_command("java -Xmx4g -jar /opt/picard/CollectAlignmentSummaryMetrics.jar  I=$bam_file  METRIC_ACCUMULATION_LEVEL=SAMPLE OUTPUT=$picard_casm_out VALIDATION_STRINGENCY=SILENT");




#########################################
# picard/CollectInsertSizeMetrics       #
#########################################
print "#########################################\n";
print "# picard/CollectInsertSizeMetrics       #\n";
print "#########################################\n\n";
&run_unix_command("java -Xmx4g -jar /opt/picard/CollectInsertSizeMetrics.jar  I=$bam_file  HISTOGRAM_FILE=$picard_cism_histogram METRIC_ACCUMULATION_LEVEL=SAMPLE OUTPUT=$picard_cism_out VALIDATION_STRINGENCY=SILENT");




#########################################
# picard/CollectMultipleMetrics         #
#########################################
print "#########################################\n";
print "# picard/CollectMultipleMetrics         #\n";
print "#########################################\n\n";
&run_unix_command("java -Xmx4g -jar /opt/picard/CollectMultipleMetrics.jar  I=$bam_file PROGRAM=CollectAlignmentSummaryMetrics PROGRAM=CollectInsertSizeMetrics PROGRAM=QualityScoreDistribution PROGRAM=MeanQualityByCycle  OUTPUT=$picard_mm VALIDATION_STRINGENCY=SILENT");

print "\n\n\n";
print "#######################################\n";
print "#    Finished analysis of BAM file    #\n";
print "#######################################\n\n";

print "BAM file: \t$bam_file\n\n";

print "Analyses carried out:\n\n";

print "\tpicard/picard/QualityScoreDistribution\n\n";

print "\tpicard/CollectAlignmentSummaryMetrics\n\n";

print "\tpicard/CollectInsertSizeMetrics\n\n";

print "\tpicard/CollectMultipleMetrics\n\n\n";


exit;



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

#############################################
# Subroutine to execute unix command        #
#############################################

sub run_unix_command
{
	my $unix_command = "";
	$unix_command = $_[0];

	print "\n";
	print("$unix_command\n");
	system("$unix_command");
	
	print "\n\n=============================================================================\n\n\n";

}

	