#!/usr/bin/perl -w
use strict;
use Term::ANSIColor;


##############################################################
# Converts BAM files from old Illumina quality scores        #
# to Sanger quality scores                                   #
# using GenomeAnalysisTK PrintReads with -fixMisencodedQuals #
##############################################################
my $version					= "2";

my $gatk_directory			= "gatk";
my $GATK_validation_stringency	= "LENIENT";

#File names
my $list_file				= "";
my $bam_file         		= "";
my $bam_file_new			= "";
my $chart_pdf				= ""; 

my $temp_dir_string			= " -Djava.io.tmpdir=javatempdir";
my $tempdir					= "javatempdir";
my $ref						= "";
my $ref_seq_name			= "";
my $species					= "";
my $dummy_dbsnp_file		= "";
my $mem						= "-Xmx4g";

#String variables
my $single_line				= "";
my $answer					= "";
my $create_quality_score_pdf = ""; # yes or no to using picard/QualityScoreDistribution to check output files

my $read_file_method		= "";
my $command					= "";
my $second_column_found		= "";
my $new_sample_name			= "";
my $prefix					= "";
my $testing_mode			= "no";
my $process_file 			= ""; # yes or no


#Number variables
my $line_count				= 0;
my $list_count				= 0;
my $loop_count				= 0;
my $pair_count				= 0;
my $no_of_files				= 0;
my $array_count				= 0;
my $array_size				= 0;

my @item					= ();
my @bam_file_array			= ();
my @bam_file_new_array		= ();
my @sample_name_array		= ();


print color 'reset';

use Term::ANSIColor;
print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      convert_bam_quals       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";

print "  - This program converts BAM file quality scores from Illumina to Sanger.\n\n";

print "  - (i.e. from old Illumina type scores to the new Illumina type scores)\n\n";

print "  - It does this using the Broad GenomeAnalysisTK PrintReads program\n\n";

print color 'reset';


####################################################
# Ask if user wants to create quality score graphs #
####################################################
print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " Do you want to create Quality Score graphs for the output files?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   <1>\tYES\n";
print "   <2>\tNO\n\n";

$answer = <STDIN>;
chomp $answer;
if ($answer eq ""){$answer = "1"} # default

$answer = substr($answer,0,1);

if (substr($answer,0,1) eq "1" ){$create_quality_score_pdf = "yes"}
if (substr($answer,0,1) eq "2" ){$create_quality_score_pdf = "no"}


################################
# Make the temp java directory #
################################
$tempdir = "$ENV{HOME}/javatempdir";
$temp_dir_string = " -Djava.io.tmpdir=$tempdir";
if (! -e $tempdir)
{
	unless(mkdir $tempdir){die "Unable to create temporary Java directory $tempdir";}
	$temp_dir_string = " -Djava.io.tmpdir=$tempdir";	
}


print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " How do you want to read the input files?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   <1>\tMULTIPLE FILES (using a file of file names)\n";
print "   <2>\tSINGLE FILE\n\n";

$answer = <STDIN>;
chomp $answer;
if ($answer eq ""){$answer = "1"} # default

$answer = substr($answer,0,1);

if (substr($answer,0,1) eq "1" ){$read_file_method = "multiple"}
if (substr($answer,0,1) eq "2" ){$read_file_method = "single"}



###################################
# Input file names                #
###################################

if ($read_file_method eq "single")
{
	until (-e $bam_file)
	{
		print "\nInput BAM file:               \t";
		
		$bam_file = <STDIN>;
		chomp $bam_file;
		
		if ($bam_file eq "ls"){print "\n";system ("ls *.bam")}
		
		if ($bam_file ne "ls")
		{
			if (! -e $bam_file){print "\n\n>>>>>>>>  File $bam_file not found.  Try again.  <<<<<<<<\n\n";}
		}
	}
	

		

	##############################################################
	# Assign to array 1 of file names but just use first element #
	# (array 2 is not used because no_of_pairs = 1)              #
	##############################################################
	$bam_file_array[1]=$bam_file;
	$no_of_files=1;
	
	
	##################################
	# Make default new BAM file name #
	##################################
	$prefix = &get_prefix ($bam_file);
	$bam_file_new = "$prefix"."_newquals.bam";
	$chart_pdf = "$prefix"."_quality_scores.pdf";
			
	print   "Output file (new BAM file):   (default=$bam_file_new)\t";
	
	$answer = <STDIN>;
	chomp $answer;
	
	if ($answer ne "")
	{
		$bam_file_new = $answer;
	}
	
	$bam_file_new_array[1]=$bam_file_new;

} # End of read = single


if ($read_file_method eq "multiple")
{
	until (-e "$list_file")
	{
		print "\nPlease input the name of your file with a list of file names of the BAM files:    (type 'ls' to get a list of .txt files)  ";
		$list_file = <STDIN>;
		chomp $list_file;

		if ($list_file eq "ls")
		{
			print "\n";
			system ("ls *.txt");
			print "\n";
		}
		if ($list_file ne "ls")
		{
			if (! -e "$list_file"){print "\nFile doesn't exist. Try again...    \n";}
		}
	}
	
	#############################################
	# Make sure the list file is in Unix format #
	#############################################

	$command = "dos2unix $list_file";
	system("$command");
	
	####################################################
	# Open the list file to get the list of file names #
	####################################################
	open (LIST, "$list_file") || die "Cannot open $list_file";
	$list_count=1;
	
	while ($single_line = <LIST> ) 
	{
		chomp $single_line;
		@item=split(/\t/,$single_line);
			
		$array_size = scalar @item;
		

		
		if ($array_size == 1)
		{
			$bam_file = $item[0];
	
			$bam_file_array[$list_count]=$bam_file;
			
			#Output file names
			$prefix = &get_prefix ($bam_file);
			$bam_file_new_array[$list_count] = "$prefix"."_newquals.bam";
			
			$list_count = $list_count + 1;
		}
		
		# If there is a third column to use as the sample name
		if ($array_size == 2)
		{
			$second_column_found = "true";
			
			$bam_file = $item[0];
			$new_sample_name = $item[1];
			
			$bam_file_array[$list_count]=$bam_file;
			$sample_name_array[$list_count] = $new_sample_name;
			
			#Output file names
			$prefix = &get_prefix ($bam_file);
			$bam_file_new_array[$list_count] = "$prefix"."_newquals.bam";
			
			$list_count = $list_count + 1;
		}
			
	} # while

	close LIST;

	$no_of_files=$list_count - 1;
	
	
	###################
	# List file names #
	###################
	print "\nThere are $no_of_files BAM files in this file of file names.\n\n";
	
	if ($second_column_found eq "true")
	{
		print "There is also a second column in the input file which will be used for the output file names.\n\n";
	}
	
	for ($list_count=1;$list_count<=$no_of_files;$list_count++)
	{
		if ($second_column_found eq "true")
		{
			print "File $list_count	\t$bam_file_array[$list_count]  \t$sample_name_array[$list_count]\n";
			print "                 \tOutput file: $bam_file_new_array[$list_count]\n";
		}
		if ($second_column_found eq "false")
		{
			print "File $list_count	\t$bam_file_array[$list_count]\n";
			print "                 \tOutput file: $bam_file_new_array[$list_count]\n";
		}
	}

} # End of if ($read_file_method eq "multiple")

##################################
# Define data files              #
##################################

print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print " Which reference sequence do you want to use?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   Enter 1 for CanFam3\n";
print "   Enter 2 for CanFam2\n";
print "   Enter 3 for EquCab2\n";
print "   Enter 4 for Human\n\n";

print "   Enter 5 for Strep. equi\n";
print "   Enter 6 for Strep. zoo\n\n";

print "   Enter 9 for other\n\n";


$answer = <STDIN>;
chomp $answer;

if (substr($answer,0,1) eq "1" ){$ref = "/home/genetics/canfam3/canfam3.fasta"; $ref_seq_name = "canfam3"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3/canfam3_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "2" ){$ref = "/home/genetics/canfam2/canfam2.fasta"; $ref_seq_name = "canfam2"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam2/canfam2_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "3" ){$ref = "/home/genetics/equcab2/equcab2.fasta"; $ref_seq_name = "equcab2"; $species = "equus_caballus";$dummy_dbsnp_file = "/home/genetics/equcab2/equcab2_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "4" ){$ref = "/home/genetics/human/human.fasta"; $ref_seq_name = "human"; $species = "homo_sapiens";$dummy_dbsnp_file = "/home/genetics/human/human_dummy_DBSNP.vcf"}

if (substr($answer,0,1) eq "5" ){$ref = "/home/genetics/strep_equi/strep_equi.fasta"; $ref_seq_name = "s_equi"; $species = "streptococcus_equi";$dummy_dbsnp_file = "/home/genetics/strep_equi/strep_equi_dummy_DBSNP.vcf"}
if (substr($answer,0,1) eq "6" ){$ref = "/home/genetics/strep_zoo/strep_zoo.fasta"; $ref_seq_name = "s_zoo"; $species = "streptococcus_zoo";$dummy_dbsnp_file = "/home/genetics/strep_zoo/strep_zoo_dummy_DBSNP.vcf"}


if (substr($answer,0,1) eq "u" ){$ref = "/home/genetics/canfam3nu/canfam3nu.fasta"; $ref_seq_name = "canfam3nu"; $species = "canis_familiaris";$dummy_dbsnp_file = "/home/genetics/canfam3nu/canfam3nu_dummy_DBSNP.vcf"}




print "\nConverting quality scores from OLD to NEW...\n\n";


for ($loop_count=1;$loop_count <=$no_of_files;$loop_count++)
{
		
		if ($read_file_method eq "single")
		{
			$bam_file = $bam_file_array[$loop_count];
			$bam_file_new = $bam_file_new_array[$loop_count];
			
			if ($bam_file ne ""){$process_file = "yes";}
		}
		
		if ($read_file_method eq "multiple")
		{
			$bam_file = $bam_file_array[$loop_count];
			$bam_file_new = $bam_file_new_array[$loop_count];
		
			$prefix = &get_prefix ($bam_file);
			
			$chart_pdf = "$prefix"."_quality_scores.pdf";
			
			print "\n\n";
			print "Loop $loop_count  File: $pair_count\n\n";
			print "BAM file:     \t$bam_file\n";
			print "New BAM file: \t$bam_file_new\n\n";
			
			if ($bam_file ne ""){$process_file = "yes";}

		} # read file method is multiple
		
	
	######################################################################################
	# Use GenomeAnalysisTK/PrintReads to convert the base quality scores                 #
	######################################################################################

	print	"\n\n";
	print   "\n====================================================================================================";
	print 	"\nFile $loop_count/$no_of_files  Using GenomeAnalysisTK/PrintReads to convert the base quality scores ";
	print 	"\n====================================================================================================\n\n";


	&run_unix_command ("java  $mem $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -T PrintReads  -I $bam_file -o $bam_file_new -R $ref  -fixMisencodedQuals -S $GATK_validation_stringency ");

	&record_output_file_size ("$bam_file_new");

	print	"\n====================================================================================================";
	print 	"\nFile $loop_count/$no_of_files   Used GenomeAnalysisTK/PrintReads to convert the base quality scores ";
	print 	"\n====================================================================================================\n\n";

	print "\n* * * \n";
	
	
	if ($create_quality_score_pdf eq "yes")
	{
		######################################################################################
		# Use picard/QualityScoreDistribution to check the base quality scores               #
		######################################################################################

		print	"\n\n";
		print   "\n=====================================================================================================";
		print 	"\nFile $loop_count/$no_of_files  Using picard/QualityScoreDistribution to check the base quality scores ";
		print 	"\n=====================================================================================================\n\n";


		&run_unix_command ("java  $mem $temp_dir_string -jar /opt/picard/QualityScoreDistribution.jar I=$bam_file_new O=junk.out CHART=$chart_pdf VALIDATION_STRINGENCY=SILENT");

		&record_output_file_size ("$chart_pdf");

		print	"\n=====================================================================================================";
		print 	"\nFile $loop_count/$no_of_files   Used picard/QualityScoreDistribution to check the base quality scores ";
		print 	"\n=====================================================================================================\n\n";

		print "\n* * * \n";
	}
} # end of main loop

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

}

##############################################
#                                            #
# Subroutine to record size of output file   #
#                                            #
##############################################
sub record_output_file_size
{
	my $outputfile = "";
	my $filesize = "";	

	$outputfile = $_[0];
	
	if (-e $outputfile)
	{
		$filesize = -s "$outputfile";
		print "\n  Output file: $outputfile\t\tSize: $filesize\n\n";
	}
	else
	{
		print "\n  Output file: $outputfile\t\tSize: Not found\n\n";
	}
	
}
