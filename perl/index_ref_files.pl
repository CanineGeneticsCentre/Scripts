#!/usr/bin/perl -w

#########################################################################
#									                                    #      
#	INDEX REF FILES           						                    #     
#									                                    #
#									                                    #
#########################################################################

#############################
# Oliver Forman Sep 2010    #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
# oliver.forman@aht.org.uk  #
#############################

# May 2012
# Modifed my Mike Boursnell to work for just indexing the ref file
# and to work on the new server

use strict;
use File::Basename ;

# VERSION OF SOFTWARE #
my $version						= "8";

# Constants
my $bwa_path						= "/opt/bwa/bwa";

####################
# Define variables #
####################

#Constants
my $testing_mode				= "off";
my $make_dbsnp_file				= "yes";

my $list_count					= 0; #Counter for the list of input files
my $no_of_files					= 0; #No of files in the list of input files (or no of lines if Paired Ends)
my $loop_count					= 0; #Loops round once for each of the files in the list of files
my $read_length					= 0;
my $run_time					= 0;
my $results_folder_count		= 0; # Counts existing results folder to see if run has already partially run.
my $next_folder					= 0;
my $start_folder				= 1;
my $snp_count					= 0; # counter to name SNPs in dummy DBSNP file

my $second_dict_file			= "";
my $snp_id						= "";
my $dummy_dbsnp_file			= "";
my $command						= "";
my $ref							= "";
my $ref_prefix					= "";
my $reads						= "";
my $mem							= "";
my $proceed						= "";
my $data_type					= "";
my $data						= "";
my $reads2						= "";
my $email_address							= "";
my $run_title					= "Index Ref File";
my $platform					= "Illumina";
my $platform_type 				= "";
my $reads3 						= "";
my $suffix						= "";
my $annotate_variants			= "";
my $maq							= "";
my $ref_size					= "";
my $ref_genome					= "";
my $structural_variant_analysis	= "";
my $chromosome					= "";
my $insert						= "";
my $region						= "";
my $use_defined_region			= "";
my $input_file					= "";
my $bam_file					= "";
my $filename					= "";
my $list_file					= "";
my $read_file_method			= "";
my $answer						= "";
my $single_line					= "";
my $results_folder				= "";
my $file_to_be_moved			= "";
my $input_string				= "";
my $ref_seq_name				= "";  # Name of referecne sequence for pindel analysis
my $log_file					= "log.rtf";
my $exit_on_problem				= "";
my $sample_name					= "";  # Name of each sample. Used to give the .vcf file an individual sample name
my $use_preindexed_reference	= "";  # can be 'yes' or 'no'  Flags whether you want to use pred-indexed ref seq such as CanFam2 or EquCab1

my $command_log					= "";
my $bwt_algorithm				= "";

my @item						= ();



########################
# Define non variables #
########################

my $title='Perl Mail';
my $from= 'NGS_analysis@samba64.aht.org.uk';
my $subject='INDEXING REF FILES';

			
			
################
# START TIMERS #
################
BEGIN { our $start_run = time(); }
BEGIN { our $start_run2 = time(); }
BEGIN { our $start_run3 = time(); }

use Term::ANSIColor;
print color 'bold cyan';

$command =  "clear";
system("$command");


print color 'reset';

use Term::ANSIColor;
print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      INDEX_REF_FILES       \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";

print "  - This program indexes reference files for bwa\n\n";

print "  - Don't do this unless you KNOW the old ones are faulty\n\n";

print color 'reset';


######################################
# E-mail address for notifications   #
######################################

&print_message("Please enter your email address.","input");
$email_address = <STDIN>;
chomp $email_address;

if ($email_address eq "m"){$email_address = 'mike.boursnell@aht.org.uk';}
if ($email_address eq "o"){$email_address = 'oliver.forman@aht.org.uk';}
if ($email_address eq "k"){$email_address = 'karen.steward@aht.org.uk';}


#######################################
# Get name of Reference sequence file #
#######################################

until (-e "$ref")
{
	&print_message("Please input the name of your reference fasta file (full path)","input");
	$ref = <STDIN>;
	chomp $ref;

	if (! -e $ref){print "\n\n>>>>>>    $ref not found!.  Try again   <<<<<<\n\n";}
}

#########################################################
# Check if REF file exists (if a REF file is specified) #
########################################################

if (! -e "$ref")
{ 
	print "\nFile $ref does not exist\n\n";
	exit;
}


&print_message("Is your sequence larger than 2GB?","input");

print "  <1>  Yes - larger then 2GB\n";
print "  <2>  No - smaller than 2GB\n\n";

$answer = <STDIN>; chomp $answer;

if (substr($answer,0,1) eq "1" ){$bwt_algorithm = "bwtsw"}
if (substr($answer,0,1) eq "2" ){$bwt_algorithm = "is"}



if (-e "$ref.dict"){print "$ref.dict exists\n"};
if (-e "$ref.bwt"){print "$ref.bwt exists\n"};
if (-e "$ref.pac"){print "$ref.pac exists\n"};
if (-e "$ref.sa"){print "$ref.sa exists\n"};
if (-e "$ref.amb"){print "$ref.amb exists\n"};
if (-e "$ref.ann"){print "$ref.ann exists\n"};
if (-e "$ref.fai"){print "$ref.fai exists\n"};

if ((-e "$ref.dict") && (-e "$ref.bwt") && (-e "$ref.pac") && (-e "$ref.sa") && (-e "$ref.amb") && (-e "$ref.ann") && (-e "$ref.fai"))
{ 
	&print_message("WARNING!!","warning");
	
	print "\nIt looks as if the bwa-0.75 indexing files for $ref already exist\n\n";
	print "Do you really want to re-do the indexing?  (y/n)  ";
	$answer = <STDIN>;
	chomp $answer;

	if (lc $answer eq "n"){exit}
} 



#####################################################################
# Check if ref.dict already exists - picard will not run if it does #
#####################################################################
if (-e "$ref.dict")
{
	&print_message("File $ref.dict exists and will have to be deleted before running a new indexing.","message");

	print "  Please delete it first and then run the script again.";
	exit;
}


###################################################################
# Get bit before name of ref sequence to delete "test.dict" later #
###################################################################
$ref_prefix = &get_prefix ($ref);

#################################
# Make name of dummy dbsnp file #
#################################

$dummy_dbsnp_file = "$ref_prefix"."_dummy_DBSNP.vcf";


#########################
# open Command Log file #
#########################
$command_log = "index_ref_files_command_log.out";

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "COMMAND LOG for index_ref_files $version\n\n";
print COMMAND_LOG "Analysis name: $run_title\n\n";


##########################################
# 1 Index the reference sequence         #
##########################################

############################################
#  This command makes the following files: #
#  .ann, .amb, .pac, rpac                  #
#  .bwt, .rbwt, .sa, .rsa                  #
############################################

&run_unix_command("$bwa_path index -a $bwt_algorithm $ref","BWA INDEX");

&print_message("The reference sequence has been indexed using bwa $bwt_algorithm","message");




##########################################
# 2 Create Sequence Dictionary           #
##########################################

&run_unix_command("java $mem -jar /opt/picard/CreateSequenceDictionary.jar R=$ref O=$ref.dict","CreateSequenceDictionary");

&print_message("The dictionary has been created:  $ref.dict","message");

print "\n\n\n";

########################################################################################
# 3 Create the fasta index file       
########################################################################################

&run_unix_command("/opt/samtools/samtools faidx $ref","Create FASTA index .fai");

&print_message("The fasta index file has been created:  $ref.fai","message");

print "\n\n\n";

########################################################################################
# 4 Make a dummy DB SNP file       
########################################################################################
if ($make_dbsnp_file eq "yes")
{
	print "\nCreating dummy DBSNP file $dummy_dbsnp_file...\n\n";
	
	open (REF, "$ref") || die "Cannot open $ref";
	open (DBSNP, ">$dummy_dbsnp_file") || die "Cannot open $dummy_dbsnp_file";
	
	print DBSNP "##fileformat=VCFv4.0\n";
	print DBSNP "##source=index_ref_files\n";
	print DBSNP "##reference=$ref\n";
	print DBSNP "#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n";

	while ($single_line = <REF> ) 
	{
		chomp $single_line;

		# Look for >
		if (index($single_line,">") == 0)
		{
			@item=split(/\s+/,$single_line);
			
			$chromosome = $item[0];
			$chromosome = substr($chromosome,1);
			
			$snp_count = $snp_count + 1;
			if ($snp_count < 10){$snp_id = "SNP_0"."$snp_count";}
			if ($snp_count >= 10){$snp_id = "SNP_"."$snp_count";}
			
			print DBSNP "$chromosome\t3\t$snp_id\tT\tA\t.\t.\t.\n";
		}
		
	}
	
close REF;
close DBSNP;
	
}

&print_message("A dummy DBSNP file has been created:  $dummy_dbsnp_file","message");

########################################################################################
# 5 Copy the dict file example.fasta.dict --> example.dict (it seems both are needed)     
########################################################################################

$second_dict_file = "$ref_prefix".".dict";

print "\n\n\n";

&run_unix_command("cp $ref.dict $second_dict_file","Make copy of dict file");

&print_message("A second copy of the $ref.dict file has been made called $second_dict_file","message");


#####################
# Send mail to user #
#####################
open(MAIL, "|/usr/sbin/sendmail -t");

## Mail Header
print MAIL "To: $email_address\n";
print MAIL "From: $from\n";
print MAIL "Subject: INDEX REFERENCING ($run_title)\n\n";
## Mail Body
print MAIL "Your reference sequence $ref has been indexed using $bwa_path with algoritm $bwt_algorithm\n";
close(MAIL);

&print_message("INDEXING COMPLETE! YOU HAVE BEEN NOTIFIED BY EMAIL","message");

close (COMMAND_LOG);
exit;


#############################################
#                                           #
# Subroutine to delete files              r #
#                                           #
#############################################

sub delete_file
{

	my $file_to_be_deleted = "";	

	$file_to_be_deleted = $_[0];
	$command = "rm  $file_to_be_deleted";
	print("$command\n");
	system("$command");
	print COMMAND_LOG "$command\n";


}


#############################################
#                                           #
# Subroutine to execute unix command        #
#                                           #
#############################################

sub run_unix_command
{
	my $unix_command = "";
	my $step = "";	
	$unix_command = $_[0];
	$step = $_[1];
		
	print("$unix_command\n\n");

	print COMMAND_LOG "\n------------------------------------------------------------------------------------------------------------------------------------\n";
	print COMMAND_LOG "$step\n";
	print COMMAND_LOG "$unix_command\n";
	
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
	}
	else
	{
		$filesize="Not found";
	}
	print COMMAND_LOG "\n  Output file: $outputfile\t\tSize: $filesize\n";

}
	
		
#############################################
#                                           #
# Subroutine to move log     			    #
#                                           #
#############################################

sub move_log
{

	my $file_to_be_moved = "";	

	$file_to_be_moved = $_[0];
	$command = "mv  $file_to_be_moved $results_folder/$file_to_be_moved";
	system("$command");

}

####################################################################
#                                                                  #
# Subroutine to get filename before file type suffix (i.e. prefix) #
# (e.g. to get "test" from "test.fasta")                           #
#                                                                  #
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



######################################
# Subroutine to print screen message #
######################################

sub print_message
{
	my $message_length 	= "";
	my $pos_count		= 0;
	my $char			= "";
	
	my $message = $_[0];
	my $style = $_[1];
	
	$message_length = length($message);
	
	if ($style eq ""){$char = "#"}
	if ($style eq "input"){$char = "~"}
	if ($style eq "message"){$char = "#"}
	if ($style eq "warning"){$char = "!"}
	
	print "\n\n";
	print color ' bold yellow';
	if ($style eq "warning"){print color ' bold red'}
	if ($style eq "input"){print color ' bold white'}
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n$char    $message    $char\n";
	
	for ($pos_count = 1;$pos_count <=($message_length + 10);$pos_count++){print $char}
	
	print "\n\n";
	print color 'reset';

}

###############
# Help on bwa #
###############
sub bwa_help
{
	&print_message("Help on bwa","message");

	print "You really need to read about bwa on the web.\n\n";

	print "Go to http://bio-bwa.sourceforge.net\n\n";

	print "Basically the new version of bwa has extra features that are better than the older versions\n";
	print "but it needs a whole new set of indexing files.  The indexing files have the SAME NAMES\n";
	print "as the old ones but are DIFFERENT, so need to go in a separate folder.\n\n";
}