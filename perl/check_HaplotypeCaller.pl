#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	check_HaplotypeCaller						                        #     
#									                                    #
#									                                    #
#########################################################################

#############################
# Mike Boursnell April 2016 #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;
use Cwd;

my $version							= "w3";
my $testing_mode					= "off";
#my $flank_size						= 1000;
my $flank_size						= 500;

my $memory_in_Gb				= 0; # memory setting for Java in gigabytes

my $no_of_threads_gatk_ug_nct	= "1";  # UnifiedGenotyper -nct (or HaplotypeCaller) (UG has to be 1)
my $no_of_threads_gatk_ug_nt	= "2";  # UnifiedGenotyper -nt (or HaplotypeCaller)

my $no_of_threads_gatk_pr_nct	= "2";  # PrintReads -nct                <== This one is used for PrintReads
my $no_of_threads_gatk_br_nct	= "2";  # BaseRecalibrator -nct          <== This one is used for BaseRecalibrator
my $no_of_threads_gatk_rtc_nt	= "2";
my $workstation 				= "unknown";
my $e_mail_from 				= ""; # Who e-mails come from

# Workstation version - (if 'w' in version name)
if (index($version,"w") > -1)
{
	$memory_in_Gb				= 45; # memory setting for Java in gigabytes
	$no_of_threads_gatk_ug_nct	= "8";  # UnifiedGenotyper -nct (or HaplotypeCaller) (UG has to be 1) # 8 seems to work now. Mike 30/3/2016
	$no_of_threads_gatk_ug_nt	= "16";  # UnifiedGenotyper -nt (but NOT HaplotypeCaller)
	$no_of_threads_gatk_pr_nct	= "8";  # PrintReads -nct                <== This one is used for PrintReads
	$no_of_threads_gatk_br_nct	= "8";  # BaseRecalibrator -nct          <== This one is used for BaseRecalibrator
	$no_of_threads_gatk_rtc_nt	= "16"; # RealignerTargetCreator         <== This one is used for RealignerTargetCreator
	$workstation 				= "true";
	$e_mail_from 				= 'NGS_analysis@gen-x1404-ws01.aht.org.uk'; # Who e-mails come from
}
else
# Samba64 version (if no 'w' in version name)
{
	$memory_in_Gb				= 4; # memory setting for Java in gigabytes
	$no_of_threads_gatk_ug_nct	= "1";  # UnifiedGenotyper -nct (or HaplotypeCaller) (UG has to be 1)
	$no_of_threads_gatk_ug_nt	= "2";  # UnifiedGenotyper -nt (but NOT HaplotypeCaller)
	$no_of_threads_gatk_pr_nct	= "2";  # PrintReads -nct                <== This one is used for PrintReads
	$no_of_threads_gatk_br_nct	= "2";  # BaseRecalibrator -nct          <== This one is used for BaseRecalibrator
	$no_of_threads_gatk_rtc_nt	= "2"; # RealignerTargetCreator         <== This one is used for RealignerTargetCreator
	$workstation 				= "false";
	$e_mail_from 				= 'NGS_analysis@samba64.org.uk'; # Who e-mails come from
}

#Constants
my $ref								= "/geneticsdata/alces-flight-image/canfam3/ensembl/canfam3.fasta";

#File names
my $bam_file						= "";
my $bamout_file						= "";
my $output_vcf						= "";
my $IGV_track_file					= "";
my $list_file						= "";

#Strings
my $quality_scores					= ""; # Can be "New" or "Old". Assigned by user during program
my $GATK_fix_quals_string			= "";
my $GATK_region_string				= "";
my $temp_dir_string					= "";
my $tempdir							= "";
my $staticdata						= "";
my $input							= "";
my $memory_string					= "";
my $region							= ""; # chr15:2500
my $position					= ""; # 2500 from above string
my $chromosome						= ""; # 15 from above string
my $position_left					= 0;
my $position_right					= 0;
my $answer							= "";
my $complete_processing_string		= "";
my $prefix							= "";
my $input_method					= ""; # single or multiple
my $command							= "";
my $single_line						= "";

#Numbers
my $array_size						= 0;
my $position_count					= 0;
my $no_of_positions					= 0;

#Arrays
my @item							= ();
my @chromosome_array				= ();
my @position_array					= ();

#Various Parameters (e.g. for the Unified Genotyper)
my $gatk_directory					= "gatk";
my $GATK_validation_stringency		= "STRICT";
my $stand_emit_conf					= 30;
my $stand_call_conf					= 30;
my $max_alt_alleles					= 6; # This is how many alleles the INDEL caller in UnifiedGenotyper can allow

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "              check_HaplotypeCaller        \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';

print "    - If you have a problem where IGV disagrees with HaplotypeCaller\n";
print "      this script will produce a 'BAMOUT' file of how HaplotypeCaller\n";
print "      sees the region which may help explain the discrepancy.\n\n";

print color 'reset';

################################
# Make the temp java directory #
# Check if staticdata exists   #
################################
$staticdata = "$ENV{HOME}/staticdata";
if (! stat($staticdata))
{
            $tempdir = "$ENV{HOME}/javatempdir"; # Individual's HOME space (if staticdata isn't online)
            $temp_dir_string = " -Djava.io.tmpdir=$tempdir";
            if (! -e $tempdir)
            {
                        unless(mkdir $tempdir){die "Unable to create temporary Java directory $tempdir";}
                        $temp_dir_string = " -Djava.io.tmpdir=$tempdir";            
            }
} else {
            $tempdir = "$ENV{HOME}/staticdata/javatempdir"; # Moved from individual's HOME space into staticdata 02/05/14
            $temp_dir_string = " -Djava.io.tmpdir=$tempdir";
            if (! -e $tempdir)
            {
                        unless(mkdir $tempdir){die "Unable to create temporary Java directory $tempdir";}
                        $temp_dir_string = " -Djava.io.tmpdir=$tempdir";            
            }
} 


###################################################
# AGet BAM file name                              #
###################################################
$bam_file = &get_file("What is the name of your BAM file:",".bam");
$prefix = &get_prefix($bam_file);


###################################################
# Ask how you want to read the positions to check #
###################################################
&print_message("How do you want to input the position(s) to check?","input");

print "  <1>  Single position in the form chr15:235356000\n";
print "  <2>  File containing a list of positions               [DEFAULT]\n\n";

$answer = <STDIN>;
chomp $answer;
if ($answer eq ""){$answer = "2"} # DEFAULT

if ($answer eq "1"){$input_method = "single"} else {$input_method = "multiple"}


if ($input_method eq "single")
{
	###########################################
	# Get position of problem in the BAM file #
	###########################################
	&print_message("At what position is the problem you want to look at?","input");
	print "  Enter the position in the form:    chr15:52905628\n\n";

	$region = <STDIN>;
	chomp $region;

	if (index($region,"chr") > -1)
	{
		$chromosome = substr($region,3,index($region,":")-3);
	} else {
		$chromosome = substr($region,0,index($region,":"));
	}
	if (index($region,":") > -1)
	{
		$position = substr($region,index($region,":")+1,99);
	}

	$no_of_positions = 1;
	$chromosome_array[1] = $chromosome;
	$position_array[1] = $position;

}

if ($input_method eq "multiple")
{
	###########################################
	# Get position of problem in the BAM file #
	###########################################
	&print_message("What is the name of the file with the positions to check?","input");
	print "  Lines in the file can be in this form:     chr15:52905628\n";
	print "  or two tab separated columns of chromosome and position (e.g. copied from Excel)\n\n";

	until (-e $list_file)
	{
		print "\n> ";$list_file = <STDIN>;
		chomp $list_file;

		if (($list_file eq "ls") || ($list_file eq "")) {print "\n";system ("ls *.txt")}

		if (($list_file ne "ls")  && ($list_file ne ""))
		{
			if (!-e $list_file){print "\n  >>>>>>>>>>>>>  ERROR.  File $list_file cannot be found <<<<<<<<<<<\n\n";}
		}
	}

	$command = "dos2unix $list_file";
	system("$command");

	open (LIST, "$list_file") || die "Cannot open $list_file";

	while ($single_line = <LIST>) 
	{
		chomp $single_line;

		@item = split(/\s+/,$single_line);
		$array_size = scalar @item;

		print "$single_line\n";

		if ($array_size == 1)
		{
			if (index($single_line,"chr") > -1)
			{
				$chromosome = substr($single_line,3,index($single_line,":")-3);
				if ($chromosome eq "39"){$chromosome = "X"}
			} else {
				$chromosome = substr($region,0,index($region,":"));
			}
			if (index($region,":") > -1)
			{
				$position = substr($single_line,index($single_line,":")+1,99);
			}
		}

		if ($array_size == 2)
		{
			$chromosome = $item[0];
			$position = $item[1];
			if ($chromosome eq "39"){$chromosome = "X"}
		}

		if ($array_size > 2)
		{
			&print_message("There should only be 1 or 2 columns in this file","warning");
			exit;
		}

		###################################
		# Save into array for chr and pos #
		###################################
		$position_count = $position_count + 1;
		$chromosome_array[$position_count] = $chromosome;
		$position_array[$position_count] = $position;

	} # while

	$no_of_positions = $position_count;
}


########################################################
# Force complete processing     (now fixed as default) #
########################################################
#&print_message("Do you want to force HaplotypeCaller to complete processing of this site?","input");

#print "(You'll have to read the docs to find out what this does. Default is 'Yes')\n\n";

#print "  <1>  Yes       [DEFAULT]\n";
#print "  <2>  No\n\n";

#$answer = <STDIN>;
#chomp $answer;

#if ($answer eq ""){$answer = "1"}
#if ($answer eq "1"){$complete_processing_string = "-forceActive -disableOptimizations "} else {$complete_processing_string = ""}

$complete_processing_string = "-forceActive -disableOptimizations ";



if (index($region,"chr") > -1)
{
	$chromosome = substr($region,3,index($region,":")-3);
} else {
	$chromosome = substr($region,0,index($region,":"));
}
if (index($region,":") > -1)
{
	$position = substr($region,index($region,":")+1,99);
}


###########################################
# Loop through all positions              #
###########################################
$GATK_region_string ="";
print "No.  \tRegion round position\n";
print "==   \t=====================\n";

for ($position_count = 1; $position_count <= $no_of_positions; $position_count++)
{
	$chromosome = $chromosome_array[$position_count];
	$position = $position_array[$position_count];

	$position_left = $position - $flank_size;
	$position_right = $position + $flank_size;

	print "$position_count\t$chromosome:$position_left-$position_right\n";

	$GATK_region_string = $GATK_region_string." -L ".$chromosome.":".$position_left."-".$position_right;
}


print "\nIf this looks OK, press 'return' to continue\n";

$answer=<STDIN>;


#########################################################
# Make up output file names                             #
#########################################################
if ($complete_processing_string eq "")
{
	$bamout_file = $prefix."_bamout_chr".$chromosome."_".$position.".bam";
}
else
{
	$bamout_file = $prefix."_bamout_chr".$chromosome."_".$position."_force.bam";
}

$IGV_track_file = $prefix."_IGV_track_chr".$chromosome."_".$position.".gff";
$output_vcf = $prefix."_output_vcf_".$position.".vcf";



#########################################################
# Make track file for IGV (to show where the region is) #
#########################################################
open (TRACK_FILE, ">$IGV_track_file")|| die "Cannot create output file: $IGV_track_file";

print TRACK_FILE "##gff-version 3\n";
print TRACK_FILE "#track name=BAMOUT\n";

#print TRACK_TRANSPOSONS "chr1\tRefseq\tHIT\t$position\t$end_position\t.\t$strand\t.\tID=$transposon_id;color=$transposon_miss_colour\n";

print TRACK_FILE "chr".$chromosome."\tRefseq\tHIT\t$position_left\t$position_right\t.\t+\t.\tID=REGION_TO_CHECK\n";

close TRACK_FILE;

&run_unix_command_single("java $memory_string $temp_dir_string -jar /opt/$gatk_directory/GenomeAnalysisTK.jar -R $ref -T HaplotypeCaller $complete_processing_string $GATK_fix_quals_string $GATK_region_string  -I $bam_file -bamout $bamout_file -o $output_vcf -S $GATK_validation_stringency");

&print_message("CHECK HAPLOTYPECALLER FINISHED","message");

print "Look at the original BAM file, and the BAMOUT file in IGV\n";
print "along with the GFF TRACK file created\n\n";

print "Original BAM file: \t$bam_file\n";
print "BAMOUT file:       \t$bamout_file\n";
print "IGV track file:    \t$IGV_track_file\n\n";
print "Genome position:   \tchr".$chromosome.":".$position."\n\n";
print "This shows how HaplotypeCaller actually sees this region.\n\n";
exit;


################################################################################################
#################################################################
# Subroutine to get a file name                                 #
# Argument 1: Text that user sees, asking for the file          #
# Argument 2: suffix of files to search for with ls e.g. ".bed" #
#################################################################
sub get_file
{
	my $input_message 	= $_[0];
	my $_file_type	 	= $_[1];
	my $_file_name		= "";
	my $_search_string	= "";

	if ($_file_type eq ""){$_file_type = "*"}

	until ((-e $_file_name) || (lc $_file_name eq "q"))
	{
		&print_message("$input_message","input");
		print "> ";

		$_file_name = <STDIN>;
		chomp $_file_name;

		# User types 'ls'
		if (($_file_name eq "ls") || ($_file_name eq "")) {print "\n";system ("ls *"."$_file_type")}

		# Starts with 'ls' followed by search string
		if (($_file_name ne "ls") && (index ($_file_name,"ls") == 0))
		{
			$_search_string = substr($_file_name,3,99);
			print "\n";
			system ("ls *$_search_string*");
			$_file_name = "ls";
			print "\n";
		}

		if (($_file_name ne "ls")  && (lc $_file_name ne "q") && ($_file_name ne ""))
		{

			if (!-e $_file_name){print "\n  >>>>>>>>>>>>>  ERROR.  File $_file_name cannot be found <<<<<<<<<<<\n\n";}

			if (-e $_file_name)
			{
				$_file_name = $_file_name;
			}
		} # not ls

	} # until loop

	print "\n";
	if ($_file_type eq ".txt"){system("dos2unix $_file_name")}

	$_file_name = $_file_name;
} # get_file

######################################
# Subroutine to print screen message #
######################################
sub print_message
{
	my $_message_length 	= "";
	my $_pos_count		= 0;
	my $_char			= "";
	
	my $_message = $_[0];
	my $_style = $_[1];
	
	$_message_length = length($_message);
	
	if ($_style eq ""){$_char = "#"}
	if ($_style eq "input"){$_char = "~"}
	if ($_style eq "message"){$_char = "#"}
	if ($_style eq "warning"){$_char = "!"}
	
	print "\n\n";
	print color ' bold yellow';
	if ($_style eq "warning"){print color ' bold red'}
	if ($_style eq "input"){print color ' bold white'}
	
	for ($_pos_count = 1;$_pos_count <=($_message_length + 10);$_pos_count++){print $_char}
	
	print "\n$_char    $_message    $_char\n";
	
	for ($_pos_count = 1;$_pos_count <=($_message_length + 10);$_pos_count++){print $_char}
	
	print "\n\n";
	print color 'reset';

}#

#############################################
# Subroutine to execute unix command        #
# Argument 1: command to execute            #
#############################################
sub run_unix_command_single
{
	my $unix_command = "";
	$unix_command = $_[0];
		
	print("$unix_command\n\n");
	system("$unix_command");
}



###########################################################
# Subroutine to get filename before file type suffix      #
# (i.e. prefix) (e.g. to get "test" from "test.fasta")    #
# New version using rindex rather than index.  This means #
# it can deal with files like filename.something.vcf      #
###########################################################

sub get_prefix
{
	my $_filename 	= $_[0];

	if (rindex($_filename,".") > 0){
		$_filename = substr($_filename, 0, rindex($_filename,"."));
	}
	elsif (rindex($_filename,".") == -1){
		$_filename = $_filename;
	}

	if (rindex($_filename,"/") > 0){
		$_filename = substr($_filename, 0, rindex($_filename,"/"));
	}
	elsif (rindex($_filename,"/") == -1){
		$_filename = $_filename;
	}

	die($_filename);
	$_filename = $_filename;

} # get_prefix
