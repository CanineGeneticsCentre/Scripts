#!/usr/bin/perl
use warnings;
use strict;

#####################################################
#              estimate_fastq_insert_size           #
#                                                   #
# This runs various checks on the BAM file          #
#                                                   #
#  QualityScoreDistribution                         #
#                                                   #
#####################################################

#############################
# Mike Boursnell June 2013  #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# mike.boursnell@aht.org.uk #
#############################

my $version					= "1";

my $bowtie_path				= "/opt/bowtie2/bowtie2"; # Note "bowtie" means "bowtie2"
my $bowtie_index			= "/home/genetics/canfam3/bowtie_files/canfam3"; #NOTE CANFAM3 ONLY NOW
my $bwa_path				="/opt/bwa/bwa";
my $mem						= "-Xmx4g";
my $ref						= "/home/genetics/canfam3/canfam3.fasta";

my $fastq_file_1			= "";
my $fastq_file_2			= "";
my $fastq_file				= ""; # name made out of $fastq_file_1 and $fastq_file_2
my $fastq_file_1_small		= "";
my $fastq_file_2_small		= "";
my $aligned_sam				= "";
my $aligned_sorted_bam		= "";
my $aln_sa1_sai				= "";
my $aln_sa2_sai				= "";
my $char_1					= "";
my $char_2					= "";
my $guess					= "";
my $next_line_is_data		= "false";


my $number_pos_count		= 0;
my $number_pos				= 0;
my $pos						= 0;
my $array_size				= 0;
my $insert_size				= 0;
my $standard_deviation		= 0;
my $mean_insert_size		= 0;
my $insert_sum				= 0; # From SAM file
my $insert_count			= 0; # From SAM file
my $insert_mean				= 0; # From SAM file
my $line_count				= 0;
my $read_length				= 0;
my $mate_inner_dist			= 0;
my $mate_std_dev			= 0;

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

my @item					= ();

print "\n\n";
print "################################\n";
print "# estimate_fastq_insert_size   #\n";
print "################################\n\n";

print "Version: $version\n\n";

print "This program takes a sample of the lines in your FASTQ file,\n\n";
print "and runs the following analysis:\n\n";

print "  - picard/CollectInsertSizeMetrics\n\n";

print "This enables you to find the mean insert size and its SD (standard deviation)\n\n";



# Get file name #
until (-e $fastq_file_1)
{
	print "\n\nName of first FASTQ file (type 'ls' to see a list of FASTQ files):  ";
	$fastq_file_1 = <STDIN>;
	chomp $fastq_file_1;
	
	if ($fastq_file_1 eq "ls")
	{
		print"\n";
		system ("ls *.fastq");
	}
}

# Guess second file name #
until (-e $guess)
{
	for ($pos = 0; $pos <= length ($fastq_file_1); $pos++)
	{
		$char_1 = substr ($fastq_file_1,$pos,1);
		
		if ($char_1 eq "1")
		{
			$guess = substr($fastq_file_1,0,$pos)."2".substr($fastq_file_1,$pos+1,99);
		}
	}

}

# Get file name #
until (-e $fastq_file_2)
{
	print "\n\nName of second FASTQ file:         \t\t\t(default:  $guess)  ";
	
	$fastq_file_2 = <STDIN>;
	chomp $fastq_file_2;
	
	if ($fastq_file_2 eq "")
	{
		$fastq_file_2 = $guess;
	}
}



##########################################
# Make names for  files                  #
##########################################
$prefix = &get_prefix($fastq_file_1);
$fastq_file_1_small = $prefix."_small_temp.fastq";


$prefix = &get_prefix($fastq_file_2);
$fastq_file_2_small = $prefix."_small_temp.fastq";

# Make single name out of reads1 and reads2 #

for ($pos=0; $pos <= length($fastq_file_1); $pos ++)
{
	$char_1 = substr ($fastq_file_1,$pos,1);
	$char_2 = substr ($fastq_file_2,$pos,1);
	
	if (($char_1 eq "1") && ($char_2 eq "2"))
	{
		$number_pos = $pos;
		$number_pos_count = $number_pos_count + 1;
	}
}


$fastq_file = substr($fastq_file_1,0,$number_pos-1);
$fastq_file = $fastq_file.substr($fastq_file_1,$number_pos + 1,99);

$prefix = &get_prefix($fastq_file);
$aligned_sam = $prefix."_small_temp.sam";
$aligned_sorted_bam = $prefix."_small_temp.bam";

$aln_sa1_sai = $prefix."small_temp_sa1.sai";
$aln_sa2_sai = $prefix."small_temp_sa2.sai";
$picard_cism_histogram = $prefix."_insert_size.pdf";

$picard_cism_out = $prefix."_data.out";

print "\n\n";
print "Original file names:  \t$fastq_file_1\t\t$fastq_file_2\n";
print "Small file names:     \t$fastq_file_1_small\t\t$fastq_file_2_small\n\n";


#########################
# Make output file name #
#########################
$output_file = "$fastq_file"."_insert_size.out";

open (OUT, ">$output_file")|| die "Cannot create output file: $output_file";

print OUT "Output from PERL script estimate_insert_size\n\n";

print OUT print "FASTQ file names:  \t$fastq_file_1\t\t$fastq_file_2\n\n";

if ($prefix eq "carrots")
{
##########################################
# Make  small FASTQ files                #
##########################################

print "###############################\n";
print "# Making small FASTQ files... #\n";
print "###############################\n\n";

&run_unix_command ("head -n 40000 $fastq_file_1 > $fastq_file_1_small");
&run_unix_command ("head -n 40000 $fastq_file_2 > $fastq_file_2_small");


##########################################
# Make SAM files using bowtie2           #
##########################################

#&run_unix_command("$bowtie_path -x $bowtie_index -1 $fastq_file_1_small -2 $fastq_file_2_small -S $aligned_sam","Make SAM file with bowtie2");


##########################################
# Make BAM files using BWA               #
##########################################

&run_unix_command("$bwa_path aln $ref $fastq_file_1_small > $aln_sa1_sai","BWA 1a");
&run_unix_command("$bwa_path aln $ref $fastq_file_2_small > $aln_sa2_sai","BWA 1b");

&run_unix_command("$bwa_path sampe $ref $aln_sa1_sai $aln_sa2_sai $fastq_file_1_small $fastq_file_2_small > $aligned_sam","BWA 2");

&record_output_file_size("$aligned_sam");

		
##########################################
# Make BAM files using picard/SortSam    #
##########################################
#&run_unix_command("java $mem -jar /opt/picard/SortSam.jar I=$aligned_sam O=$aligned_sorted_bam SO=coordinate VALIDATION_STRINGENCY=SILENT","Make BAM file");

&run_unix_command("java $mem -jar /opt/picard/AddOrReplaceReadGroups.jar I=$aligned_sam O=$aligned_sorted_bam rgID='Title' LB='Lib' PL='ILLUMINA' PU=$prefix SM=$prefix SORT_ORDER=coordinate CREATE_INDEX=true VALIDATION_STRINGENCY=LENIENT","Make Bam");	
	
&record_output_file_size("$aligned_sorted_bam");


#########################################
# picard/CollectInsertSizeMetrics       #
#########################################

print "\n\n";

print "#########################################\n";
print "# picard/CollectInsertSizeMetrics       #\n";
print "#########################################\n\n";

&run_unix_command("java -Xmx4g -jar /opt/picard/CollectInsertSizeMetrics.jar  I=$aligned_sorted_bam  HISTOGRAM_FILE=$picard_cism_histogram METRIC_ACCUMULATION_LEVEL=ALL_READS OUTPUT=$picard_cism_out VALIDATION_STRINGENCY=SILENT");

&record_output_file_size("$picard_cism_out");


} # temporary



####################################
# Read SAM file to get insert size #
####################################

if (! -e $aligned_sam){print "Cannot open $aligned_sam\n\n"};

if (-e $aligned_sam)
{
	print "Analysing SAM file ...\n\n";
	
	open (SAM, "$aligned_sam");
	while ($single_line = <SAM> ) 
	{
		chomp $single_line;
		
		@item=split(/\s/,$single_line);
				
		$array_size = scalar @item;
		
		if ($array_size > 9)
		{
			$insert_size = $item[8];
			
			if (($insert_size > 10) && ($insert_size < 500))
			{
				$insert_sum = $insert_sum + $insert_size;
				$insert_count = $insert_count + 1;
			}
			
		} # array size 1
	}
	
	$insert_mean = $insert_sum / $insert_count;


}



###################################################
# Read Picard Output file to get insert size data #
###################################################
if (! -e $picard_cism_out){print "Cannot open $picard_cism_out\n\n"};


if (-e $picard_cism_out)
{
	print "Analysing Picard output file ...\n\n";
	
	open (DATA, "$picard_cism_out");
	while ($single_line = <DATA> ) 
	{
		chomp $single_line;
		
		@item=split(/\s/,$single_line);
		$array_size = scalar @item;
		
		if ($array_size > 4)
		{
			print "$array_size\t$single_line\n";
			if ($next_line_is_data eq "true")
			{
				$mean_insert_size = $item[4];
				$standard_deviation = $item[5];
			}
			if ($item[0] eq "MEDIAN_INSERT_SIZE")
			{
				$next_line_is_data = "true";
			}	
		}

	}

}

###################################################
# Read FASTQ file to get read length              #
###################################################
if (! -e $fastq_file_1_small){print "Cannot open $fastq_file_1_small\n\n"};

$line_count = 0;

if (-e $fastq_file_1_small)
{
	print "Analysing FASTQ file ...\n\n";
	
	open (FASTQ, "$fastq_file_1_small");
	while ($single_line = <FASTQ> ) 
	{
		if ($line_count < 3) 
		{
			chomp $single_line;
			
			$line_count = $line_count + 1;
			
			if ($line_count == 2)
			{
				$read_length = length $single_line;
				
			}
			
		}

	}
	
	close FASTQ;

}

############################################
# Calculate Mate Inner Distance for TopHat #
############################################

$mate_inner_dist = int ($mean_insert_size - ($read_length * 2));
$mate_std_dev = int ($standard_deviation);



#######################################
# Delete temporary intermediate files #
#######################################

print "Do you want to delete the temporary intermediate file? (y/n)  ";

$answer = <STDIN>;
chomp $answer;

if (lc $answer eq "y")
{
	&run_unix_command("rm $aligned_sam","Remove temporary SAM file");
	&run_unix_command("rm $fastq_file_1_small","Remove temporary samll FASTQ file 1");
	&run_unix_command("rm $fastq_file_1_small","Remove temporary small FASTQ file 2");
}


print "\n\n\n";
print "###################################################\n";
print "#    Finished insert size analysis of FASTQ files #\n";
print "###################################################\n\n";

print "FASTQ file 1:             \t$fastq_file_1\n";
print "FASTQ file 2:             \t$fastq_file_2\n\n\n";

print "Data from SAM file           \t$aligned_sam\n\n";

print "\tSAM insert sum:         \t$insert_sum\n";
print "\tSAM insert count:       \t$insert_count\n";
print "\tSAM insert mean:        \t$insert_mean\n\n\n";
	
print "Data from picard output file:\t$picard_cism_out\n\n";

print "\tInsert size PDF file:   \t$picard_cism_histogram\n\n";
print "\tInsert size data file:  \t$picard_cism_out\n\n";
print "\tMean insert size:       \t$mean_insert_size\n";
print "\tStandard deviation:     \t$standard_deviation\n\n";

print "Data from FASTQ file      \t$fastq_file_1_small\n\n\n";

print "\tRead length:            \t$read_length\n\n\n";

print "For TopHat:\n\n";

print "\t--mate-inner-dist:	     \t$mate_inner_dist\n";
print "\t--mate-std-dev:         \t$mate_std_dev\n\n\n";


print OUT "###################################################\n";
print OUT "#    Finished insert size analysis of FASTQ files #\n";
print OUT "###################################################\n\n";

print OUT "FASTQ file 1:             \t$fastq_file_1\n";
print OUT "FASTQ file 2:             \t$fastq_file_2\n\n\n";

print OUT "Data from SAM file        \t$aligned_sam\n\n";

print OUT "\tSAM insert sum:         \t$insert_sum\n";
print OUT "\tSAM insert count:       \t$insert_count\n";
print OUT "\tSAM insert mean:        \t$insert_mean\n\n\n";
	
print OUT "Data from picard output file $picard_cism_out\n\n";

print OUT "\tInsert size PDF file:   \t$picard_cism_histogram\n\n";
print OUT "\tInsert size data file:  \t$picard_cism_out\n\n";
print OUT "\tMean insert size:       \t$mean_insert_size\n";
print OUT "\tStandard deviation:     \t$standard_deviation\n\n\n";

print OUT "Data from FASTQ file      \t$fastq_file_1_small\n\n";

print OUT "\tRead length:            \t$read_length\n\n";

print OUT "For TopHat:\n\n";

print OUT "\t--mate-inner-dist:	     \t$mate_inner_dist\n";
print OUT "\t--mate-std-dev:         \t$mate_std_dev\n";

close OUT;

exit;


####################################################################
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test" from "test.fasta")                           #
# if name is test.dedup.bam then the long prefix is test.dedup     #
####################################################################
sub get_long_prefix
{
	my $filename = "";
	my $new_filename = "";
	my $dot_pos	= 0;
	
	$filename = $_[0];
	
	$dot_pos = rindex($filename,".");
	
	if ($dot_pos > 0)
	{
		$new_filename = substr($filename, 0, $dot_pos);
	}
	if ($dot_pos == -1)
	{
		$new_filename = $filename;
	}
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
	
	print "\n=============================================================================\n\n";

}

	
	
##############################################
# Subroutine to record size of output file   #
##############################################

sub record_output_file_size
{
	my $outputfile = "";
	my $filesize = "";	

	$outputfile = $_[0];
	
	if (-e $outputfile)
	{
		$filesize = -s "$outputfile";
		#print COMMAND_LOG "\n  Output file: $outputfile\t\tSize: $filesize\n";
		print "\n  Output file: $outputfile\t\tSize: $filesize\n\n"
	}
	else
	{
		$filesize=0;
		#print COMMAND_LOG "\n  Output file: $outputfile\t\tSize: Not found\n";
		print "\n  Output file: $outputfile\t\tSize: Not found\n\n"
	}


}