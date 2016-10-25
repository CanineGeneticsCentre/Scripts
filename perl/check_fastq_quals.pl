#!/usr/bin/perl -w
use strict;
use Term::ANSIColor;

#####################
# Check fastq quals #
#####################
my $version					= "4";

my $fastq_file         		= "";
my $single_line				= "";
my $answer					= "";
my $id_line					= "";
my $seq_line				= "";
my $plus_line				= "";
my $qual_line				= "";
my $char					= "";
my $sanger_scoring 			= "false";
my $looks_like_illumina 	= "";
my $looks_like_sanger 		= "";
my $some_less_than_64		= "false";
my $quick_check				= ""; # yes or no

my $line_count				= 0;
my $array_count				= 0;
my $mod						= 0;
my $level					= 0;
my $block_count				= 0;
my $pos						= 0;
my $qual_more_than_40_count	= 0;
my $qual_less_than_64_count	= 0;
my $qual_41_count			= 0;
my $level_size				= 0;
my $threshold				= 0;
my $ascii_score				= 0;
my $old_score				= 0;
my $new_score				= 0;
my $no_of_levels			= 40; # levels in histogram of scores.
my $max_qual_score			= 0;
my $min_qual_score			= 1000;
my $max_ascii_score_array	= 0;
my $max_ascii_score			= 0;
my $min_ascii_score			= 1000;
my $min_new_score			= 1000;
my $max_new_score			= 0;
my $min_old_score			= 1000;
my $max_old_score			= 0;
my $sanger_more_than_40_count	= 0;
my $sanger_41_count			= 0;

my @ascii_score_array		= ();
my @new_score_array			= ();
my @old_score_array			= ();

print "\n\n";
print color 'bold magenta';

print "################################\n";
print color 'bold white';
print "      check_fastq_quals       \n";
print color 'bold magenta';
print "################################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This PERL script checks whether FASTQ file quality scores are Old or New\n\n";
print "    Some old Illumina FASTQ files have a different way of coding the quality scores.\n";
print "    The GATK programs will not work if some files have the Old scores.\n\n";
print "    (All recent data should be in the New system).\n\n";
print color 'reset';



until (-e $fastq_file)
{
	print "Fastq file: ";
	$fastq_file = <STDIN>;
	chomp $fastq_file;

	
		
	if ($fastq_file eq "ls"){print "\n";system ("ls *.fastq");print "\n";}
			
	if ($fastq_file ne "ls")
	{
		if (index($fastq_file,".") == -1 ){$fastq_file = $fastq_file.".fastq"}
		
		if (! -e $fastq_file){print "\n\n>>>>>>>>  File $fastq_file not found.  Try again.  <<<<<<<<\n\n";}
	}
}


#print "\n\n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
#print "Do you want to do a quick check (i.e. only the first 100,000 lines)?   \n";
#print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";


#print "   <1>  YES\n";
#print "   <2>  NO\n\n";

#$answer = <STDIN>;
#chomp $answer;


#if (substr($answer,0,1) eq "1"){$quick_check = "yes"}
#if (substr($answer,0,1) eq "2"){$quick_check = "no"}

# Quick check is always enough
$quick_check = "yes";


print "\nAnalysing quality scores...\n\n";

open (IN, "$fastq_file") || die "Cannot open $fastq_file";

######################
# Set up score array #
######################
for ($array_count = 0;$array_count <=126;$array_count++)
{
	$ascii_score_array[$array_count] = 0;
	$old_score_array[$array_count] = 0;
	$new_score_array[$array_count] = 0;
}

############################################
# Read through FASTQ file a line at a time #
############################################
while ($single_line = <IN>)
{
	$line_count = $line_count + 1;
	$mod = $line_count % 4;
	
	
	
	chomp $single_line;
	
	if ($mod == 0)
	{
		$qual_line = $single_line;
	}
	
	####################################
	# Last line of block of four lines #
	####################################
	if ($mod == 0)
	{
		$block_count = $block_count	 + 1;
		
		if (($block_count % 10000) == 0) 
		{
			print "$block_count\n";
			
			if ($quick_check eq "yes")
			{
				if ($block_count > 50000)
				{
					last;
				}
			}
		} # show progress


		##########################################
		# Run along the string of quality scores #
		##########################################
		for ($pos = 0; $pos < length $qual_line; $pos++)
		{
			$char = substr ($qual_line, $pos, 1);

			if (length ($char) > 0) 
			{
				#####################################################
				#  Work out what old or new quality scores would be #
				#####################################################
				$ascii_score = ord($char);
				
				$old_score = ord($char)- 64; # Before Illumina 1.8
				$new_score = ord($char)- 33; # Sanger


				###############################################
				# Any ASCII less than 64 (means it is Sanger) #
				###############################################
				if ($ascii_score < 64){$some_less_than_64 = "true";}
								
			
				##################################################
				# Record New scores over 40 (the GATK threshold) #
				##################################################
				if ($new_score > 40) 
				{
					$sanger_more_than_40_count = $sanger_more_than_40_count + 1;
				}
				if ($new_score == 41)
				{
					$sanger_41_count = $sanger_41_count + 1;
				}
				

				########################################
				# Store scores in arrays for Histogram #
				########################################
				if ($ascii_score >= 0 ) {$ascii_score_array[$ascii_score] = $ascii_score_array[$ascii_score] + 1}
				if ($new_score >= 0 ) {$new_score_array[$new_score] = $new_score_array[$new_score] + 1}
				if ($old_score >= 0 ) {$old_score_array[$old_score] = $old_score_array[$old_score] + 1}

			} # if char <> ""
			

			# Get max and min quals #
			if ($ascii_score > $max_ascii_score)
			{
				$max_ascii_score = $ascii_score;
				
				$max_old_score = $max_ascii_score- 64;
				$max_new_score = $max_ascii_score- 33;
				
			}
			if ($ascii_score < $min_ascii_score)
			{
				$min_ascii_score = $ascii_score;
				
				$min_old_score = $min_ascii_score- 64;
				$min_new_score = $min_ascii_score- 33;
			}
		}
	}
}

#####################################
# List distribution of ASCII scores #
#####################################

print "\n\n\n";
print "Distribution of ASCII scores.\n\n";
print "=============================\n\n";

print "ASCII\tNumber\n";

for ($array_count = 30; $array_count <=126;$array_count++)
{
	if ($ascii_score_array[$array_count] > 0)
	{
		print "$array_count\t$ascii_score_array[$array_count]\n";
		if ($ascii_score_array[$array_count] > $max_ascii_score_array)
		{
			$max_ascii_score_array = $ascii_score_array[$array_count];
		}
	}
}
	
#############################
# Histogram of ASCII values #
#############################

$level_size = $max_ascii_score_array / $no_of_levels;
print "\n\n\n\n\n\n\n";

for ($level=0;$level<$no_of_levels;$level++)
{
	$threshold= ($no_of_levels - $level) * $level_size;

	for ($array_count = 30; $array_count <=126;$array_count++)
	{		
		if ($ascii_score_array[$array_count] >= $threshold){print "|";}else{print " ";}	
	}
	print "\n";
}

print "\n30........40........50........60........70........80........90........100.......110.......120.......130\n";

print "\n\n";

print "Histogram of ASCII values. \tASCII values $min_ascii_score-$max_ascii_score\n\n";

print " >> For OLD quality scores subtract 64, for NEW quality scores subtract 33 (from the ASCII value)\n\n";


print "    If quality scores in this file are NEW they would go from:       \t$min_new_score-$max_new_score\n";
print "    If quality scores in this file are OLD they would go from:       \t$min_old_score-$max_old_score\n\n";
	

#Decide which system
print "\nFactors making it look like NEW scores:\n\n";

if ($max_new_score <= 40)
{
	print "\tThere are no NEW scores over 40\n";
	$looks_like_sanger = "true";
}

if ($some_less_than_64 eq "true")
{
	print "\tThere are some ASCII values less than 64.\n";
	$looks_like_sanger = "true";
}

if ($looks_like_sanger eq "")
{
	print "\tNONE\n";
}

print "\n\nFactors making it look like OLD scores\n\n";


if (($max_new_score > 40) && ($max_old_score <=40))
{
	print "\tThe max NEW score is over 40 and the max OLD score is not over 40\n";
	$looks_like_illumina = "true";
}


if ($some_less_than_64 eq "false")
{
	print "\tThere are no ASCII values less than 64.\n";
	$looks_like_illumina = "true";
}

if ($looks_like_illumina eq "")
{
	print "\tNONE\n";
}


if ($sanger_more_than_40_count > 0)
{
	print "\n\nN.B. If you use this file as NEW scores there are $sanger_more_than_40_count scores more than 40.  This will cause a problem in GATK\n";
	$looks_like_illumina = "true";
}

print "\nFastq file checked:   $fastq_file\n\n";


