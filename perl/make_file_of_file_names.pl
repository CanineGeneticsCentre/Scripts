#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	Make file of file names    									        #
#									                                    #
#########################################################################

#############################
# Mike Boursnell Mar 2013   #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use warnings;

use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;
use Cwd;
my $version							= "6";

my $file_count						= 0;
my $check_count						= 0;
my $no_of_files						= 0;
my $found_count						= 0;
my $pos								= 0;
my $chrom_line_array_size			= 0;
my $no_of_samples_chrom_line		= 0;
my $array_count						= 0;

my $single_line						= "";
my $vcf_file						= "";
my $file_to_find 					= "";
my $current_directory				= "";
my $folder							= "";
my $bam_file						= "";
my $fastq_file						= "";
my $fastq_file_1					= "";
my $fastq_file_2					= "";
my $back_to_one						= "";
my $answer							= "";
my $files_to_use					= "";
my $command							= "";
my $char_1							= "";
my $char_2							= "";
my $guess							= "";
my $list_file_unsorted				= "fof_unsorted.txt";
my $list_file						= "fof.txt";
my $list_file_paired_unsorted		= "fof_paired_unsorted.txt";
my $list_file_paired				= "fof_paired.txt";

my @file_name_array					= ();
my @fastq_file_1_array				= ();
my @fastq_file_2_array				= ();
my @chrom_line_array				= ();
my @sample_name_array				= (); # array of samples freom VCF file #CHROM line

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "#######################################\n";
print color 'bold white';
print "              make_fof                \n";
print "      (make file of file names)                \n";
print color 'bold magenta';
print "#######################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - This program helps make a file of file names from BAM or FASTQ files\n\n";

print color 'reset';

print "\n\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";
print "Which type of file do you want to make into a file of file names?\n";
print "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n";
print "   <1>  BAM        \tUse all the BAM files in the folder\n";
print "   <2>  FASTQ      \tUse all the FASTQ files in the folder\n";
print "   <3>  FASTQ.GZ   \tUse all the FASTQ gzipped files in the folder\n";
print "   <4>  gVCF       \tUse all the gVCF files in the folder\n";
print "   <5>  gsVCF      \tUse all the gsVCF files in the folder\n";
print "   <6>  gVCF.gz    \tUse all the gVCF bgzipped files in the folder\n\n";
print "   <7>  VCF/VCF.gz \tUse the file names in the columns of a VCF file\n\n";


$answer=<STDIN>;
chomp $answer;
if (substr($answer,0,1) eq "1" ){$files_to_use = "bam"}
if (substr($answer,0,1) eq "2" ){$files_to_use = "fastq"}
if (substr($answer,0,1) eq "3" ){$files_to_use = "fastq_gz"}
if (substr($answer,0,1) eq "4" ){$files_to_use = "gVCF"}
if (substr($answer,0,1) eq "5" ){$files_to_use = "gsVCF"}
if (substr($answer,0,1) eq "6" ){$files_to_use = "gVCF_gz"}
if (substr($answer,0,1) eq "7" ){$files_to_use = "vcf"}



$current_directory=getcwd;
opendir (DIR, $current_directory);

open (OUT, ">$list_file_unsorted") || die "Cannot open $list_file_unsorted";

if ($files_to_use eq "vcf")
{
	until (-e $vcf_file)
	{
		print "   VCF file:  ";
		$vcf_file = <STDIN>;
		chomp $vcf_file;
		if ($vcf_file eq "ls"){print "\n";system ("ls *.vcf");print "\n"}
		if ($vcf_file ne "ls"){if (! -e $vcf_file){print "\n\n>>>>>>>>  File $vcf_file not found.  Try again.  <<<<<<<<\n\n";}}
	}	
	
	##################################################
	# Use bcftools to query VCF file for sample list #
	##################################################
	
	@sample_name_array = `bcftools query -l $vcf_file`;
	my $line = 1;
	
	print "\n\nThese are the sample names from the VCF file:\n\n";
	
	foreach my $sample (@sample_name_array){
		chomp $sample;
		print "$line: \t$sample\n";

		print OUT "$sample";

		################################
		# Add disease status as "Omit" #
		################################

		print OUT "\tOmit\n";
		$line++;
	}
} # if files_to_use = vcf

else{
	
	while (my $file_to_find = readdir(DIR)) 
	{
			if ($files_to_use eq "bam")
			{
				if ((substr($file_to_find, -4) eq ".bam") && !($file_to_find =~ m/^\./))
				{
					$file_count = $file_count + 1;
					$file_name_array[$file_count] = $file_to_find;
					
					print "Bam file $file_count: \t$file_to_find\n";
					print OUT "$file_to_find\n";
				}
			} # bam
			
			if ($files_to_use eq "fastq")
			{
				if ((substr($file_to_find, -6) eq ".fastq") && !($file_to_find =~ m/^\./))
				{
					$file_count = $file_count + 1;
					$file_name_array[$file_count] = $file_to_find;
					
					print "Fastq file $file_count: \t$file_to_find\n";
					print OUT "$file_to_find\n";
					
				}
			} # fastq
	
			if ($files_to_use eq "fastq_gz")
			{
				if ((substr($file_to_find, -9) eq ".fastq.gz") && !($file_to_find =~ m/^\./))
				{
					$file_count = $file_count + 1;
					$file_name_array[$file_count] = $file_to_find;
					
					print "Fastq.gz file $file_count: \t$file_to_find\n";
					print OUT "$file_to_find\n";
					
				}
			} # fastq_gz
	
			if ($files_to_use eq "gVCF")
			{
				if ((substr($file_to_find, -5) eq ".gVCF") && !($file_to_find =~ m/^\./))
				{
					$file_count = $file_count + 1;
					$file_name_array[$file_count] = $file_to_find;
					
					print "gVCF file $file_count: \t$file_to_find\n";
					print OUT "$file_to_find\n";
					
				}
			} # gVCF
	
			if ($files_to_use eq "gsVCF")
			{
				if ((substr($file_to_find, -6) eq ".gsVCF") && !($file_to_find =~ m/^\./))
				{
					$file_count = $file_count + 1;
					$file_name_array[$file_count] = $file_to_find;
					
					print "gsVCF file $file_count: \t$file_to_find\n";
					print OUT "$file_to_find\n";
					
				}
			} # gVCF
	
			if ($files_to_use eq "gVCF_gz")
			{
				if ((substr($file_to_find, -8) eq ".gVCF.gz") && !($file_to_find =~ m/^\./))
				{
					$file_count = $file_count + 1;
					$file_name_array[$file_count] = $file_to_find;
					
					print "gVCF.gz file $file_count: \t$file_to_find\n";
					print OUT "$file_to_find\n";
					
				}
			} # gVCF
	
	}
	
	$no_of_files = $file_count;
	
	
	###################################################
	# Process list of FASTQ files to arrange in pairs #
	###################################################
	
	if ($files_to_use eq "fastq")
	{
		print "\n\n";
		
		open (PAIRED, ">$list_file_paired_unsorted") || die "Cannot open $list_file_paired_unsorted";
		
		for ($file_count =1; $file_count <= $no_of_files; $file_count++)
		{
			$fastq_file_1 = $file_name_array[$file_count];
			
			for ($check_count =1; $check_count <= $no_of_files; $check_count++)
			{
				$fastq_file_2 = $file_name_array[$check_count];
				
				# Guess second file name #
				for ($pos = 0; $pos <= length ($fastq_file_1); $pos++)
				{
					$char_1 = substr ($fastq_file_1,$pos,1);
					
					if ($char_1 eq "1")
					{
						# Does same position in second file have a "2"?
						
						$char_2 = substr ($fastq_file_2,$pos,1);
						
						if ($char_2 eq "2")
						{
							$back_to_one = substr($fastq_file_2,0,$pos)."1".substr($fastq_file_2,$pos+1,99);
							
							if ($fastq_file_1 eq $back_to_one)
							{
								$found_count = $found_count + 1;
								
								$fastq_file_1_array[$found_count] = $fastq_file_1;
								$fastq_file_2_array[$found_count] = $fastq_file_2;
								
								print "Found pair: $fastq_file_1    \t $fastq_file_2\n";
							}
						} # if
					}
				} # guess 2nd file name
				
			} # check loop
		}
		
		
		
		# Add to OUT file
		for ($file_count =1; $file_count <= $found_count; $file_count++)
		{
			print PAIRED "$fastq_file_1_array[$file_count]\t$fastq_file_2_array[$file_count]\n";
		}
		
		close PAIRED;
		
		$command = "sort $list_file_paired_unsorted > $list_file_paired";
		system ("$command");
	
		system ("rm $list_file_paired_unsorted");
	}
	
	
	if ($files_to_use eq "fastq_gz")
	{
		print "\n\n";
		
		open (PAIRED, ">$list_file_paired_unsorted") || die "Cannot open $list_file_paired_unsorted";
		
		for ($file_count =1; $file_count <= $no_of_files; $file_count++)
		{
			$fastq_file_1 = $file_name_array[$file_count];
			
			for ($check_count =1; $check_count <= $no_of_files; $check_count++)
			{
				$fastq_file_2 = $file_name_array[$check_count];
				
				# Guess second file name #
				for ($pos = 0; $pos <= length ($fastq_file_1); $pos++)
				{
					$char_1 = substr ($fastq_file_1,$pos,1);
					
					if ($char_1 eq "1")
					{
						# Does same position in second file have a "2"?
						
						$char_2 = substr ($fastq_file_2,$pos,1);
						
						if ($char_2 eq "2")
						{
							$back_to_one = substr($fastq_file_2,0,$pos)."1".substr($fastq_file_2,$pos+1,99);
							
							if ($fastq_file_1 eq $back_to_one)
							{
								$found_count = $found_count + 1;
								
								$fastq_file_1_array[$found_count] = $fastq_file_1;
								$fastq_file_2_array[$found_count] = $fastq_file_2;
								
								print "Found pair: $fastq_file_1    \t $fastq_file_2\n";
							}
						} # if
					}
				} # guess 2nd file name
				
			} # check loop
		}
		
		
		
		# Add to OUT file
		for ($file_count =1; $file_count <= $found_count; $file_count++)
		{
			print PAIRED "$fastq_file_1_array[$file_count]\t$fastq_file_2_array[$file_count]\n";
		}
		
		close PAIRED;
		
		$command = "sort $list_file_paired_unsorted > $list_file_paired";
		system ("$command");
	
		system ("rm $list_file_paired_unsorted");
	} # fastq_gzipped
	
	closedir (DIR);
}

close OUT;
	
$command = "sort $list_file_unsorted > $list_file";
system ("$command");

system ("rm $list_file_unsorted");

print "\n\n############\n";
print "# FINISHED #\n";
print "############\n\n";

print "File of file names created:  fof.txt\n\n";

if (($files_to_use eq "fastq") || ($files_to_use eq "fastq_gz"))
{
	print "File of file names with paired FASTQ files:   \tfof_paired.txt\n\n\n";
}

if ($files_to_use eq "vcf")
{
	print "     -  The file has a second column for disease statuses.\n";
	print "     -  These can be 'Affected', 'Carrier', 'Normal' or 'Omit'.\n\n";

	print "     -  They are all currently set as 'Omit' so you will have to add your own statuses.\n\n";
}


exit;


##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

foreach (@_) {s/\n//g}  
foreach (@_) {s/\r//g}  

}

