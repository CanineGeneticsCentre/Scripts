 #!/usr/bin/perl -w

############################################################################
#									                                       #      
#	RUN snpEff                      		                               #     
#									                                       #
#	THIS PERL SCRIPT WILL RUN Varient Effector Predictor on all VCF files  #
#									                                       #
############################################################################

#############################
# Mike Boursnell July 2012  #
# Animal Health Trust       #
# Newmarket                 #
# UK                        #
#############################

use strict;
use Getopt::Long ;
use File::Basename ;
use Term::ANSIColor;


# snpEff Database versions (must be changed when required!! update_required)
my $canine_database			= "CanFam3.1.79";
my $equine_database			= "EquCab2.79";
my $strep_equi_database		= "Streptococcus_equi_4047_uid59259";

#Constants
my $version						= "11";
my $java_memory					= "-Xmx4g";
my $config_file					= "/opt/snpeff/snpEff.config"; 
my $upstream_distance			= "200";  # distance from gene where snpEff says it is upstream or downstream

# File names
my $output_file					= "";
my $new_output_file				= "";
my $vcf_file					= "";
my $vcf_file_chr1_replaced		= ""; # file with 'chr1' in strep equi files replaced with 'NC_012471'
my $command_log					= "";

#Strings
my $prefix						= "";
my $command						= "";
my $answer						= "";
my $ref							= ""; # exact snpEff ref sequence name like "CanFam3.1.71"
my $ref_seq_name				= ""; # name like canfa2, canfam3
my $start_time					= "";
my $end_time					= "";
my $run_time					= "";

###############################
# Process flags               #
###############################

GetOptions("file:s"=>\$vcf_file,"ref:s"=>\$ref_seq_name);

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "                 run_snpEff                 \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "    - This program runs snpEff\n\n";
print "      This takes a VCF file of SNP or Indel variants and works out what the effect of the mutation will be.\n";
print "      (For example is the mutation non_synonymous, in a splice region etc etc)\n\n";

print color 'reset';

if ($vcf_file eq "")
{
	until (-e $vcf_file)
	{
		&print_message("Name of the input VCF file","input");
		$vcf_file = <STDIN>;
		chomp $vcf_file;
		
		if ($vcf_file eq "ls")
		{
			print "\n";
			system ("ls *.vcf");
			print "\n";
		}
		if ($vcf_file ne "ls")
		{
			if (! -e "$vcf_file"){print "\nFile doesn't exist. Try again...    \n";}
		}
	}

	if (-e $output_file)
	{
		&print_message("An output file called $output_file already exists","warning");
		print "Do you want to continue and overwrite this file? (y/n)     ";
		$answer = <STDIN>;
		chomp $answer;
		if (lc($answer) eq "n")
		{
			print "Existing output file: \t$output_file\n";
			print "New output file:       \t";

			$new_output_file = <STDIN>;

			$output_file = $new_output_file;
		}
	} # if output file exists

} # if vcf_file not specified in command line


$prefix = &get_prefix ($vcf_file);
	
$output_file= "$prefix"."_snpEff.vcf";
$vcf_file_chr1_replaced = "$prefix"."_chr1_replaced.vcf";

########################
# Get reference genome #
########################


##################################################################
# If reference sequence is already specified in the command line #
##################################################################
if ($ref_seq_name eq "canfam3"){$ref = $canine_database;}
if ($ref_seq_name eq "equcab2"){$ref = $equine_database;}
if ($ref_seq_name eq "strep_equi"){$ref = $strep_equi_database;}

if ($ref_seq_name eq "")
{
	&print_message("Which reference sequence do you want to use?","input");

	print "   <1>  Dog\n";
	print "   <2>  Horse\n\n";

	print "   <3>  Strep. equi\n";

	$answer = <STDIN>;
	chomp $answer;

	if (substr($answer,0,1) eq "1" ){$ref = $canine_database; $ref_seq_name = "canfam3";} 
	if (substr($answer,0,1) eq "2" ){$ref = $equine_database; $ref_seq_name = "equcab2";}
	if (substr($answer,0,1) eq "3" ){$ref = $strep_equi_database; $ref_seq_name = "strep_equi";}

} # if ref is not specified in the command line


$command_log = $prefix."_snpEff_command_log.out";

open (COMMAND_LOG, ">$command_log")|| die "Cannot create output file: $command_log";
print COMMAND_LOG "COMMAND LOG for run_snpEff version $version\n\n";

print COMMAND_LOG "Reference sequence:\t$ref \t$ref_seq_name\n";
print COMMAND_LOG "Input file:        \t$vcf_file\n\n";


###############################################
# If Strep equi convert 'chr1' to 'NC_012471' #
###############################################
if ($ref_seq_name eq "strep_equi")
{
	&print_message("Replacing 'chr1' with 'NC_012471' in strep equi VCF file","message");
	$command = "sed 's/chr1/NC_012471/' $vcf_file > $vcf_file_chr1_replaced";
	
	print COMMAND_LOG "$command\n";
	
	print("\n$command\n");
	system("$command");

	print "\nConverted $vcf_file to $vcf_file_chr1_replaced\n\n";
}

&print_message("NOTE: The version of the genome used is $ref","message");

&print_message("Making annotation files with snpEff","message");

$start_time = time();

if (-e "$vcf_file")
{
	#########################################
	# Run snpEff with config file specified #
	#########################################

	if ($ref_seq_name eq "strep_equi")
	{
		$command = "java $java_memory -jar /opt/snpeff/snpEff.jar -ud $upstream_distance  -config $config_file -verbose $ref $vcf_file_chr1_replaced > $output_file";
	}
	else
	{
		$command = "java $java_memory -jar /opt/snpeff/snpEff.jar -config $config_file -verbose $ref $vcf_file > $output_file";
	}

	print COMMAND_LOG "$command\n";
	
	print("\n$command\n");
	system("$command");
}

if (! -e "$vcf_file")
{
	print "\n\nFILE CAN'T BE FOUND: $vcf_file\n\n"; 
}


$end_time = time();
$run_time = $end_time - $start_time;

print "\n\n";
print "######################\n";
print "#  snpEff FINISHED   #\n";
print "######################\n\n";

print "Input file:     \t$vcf_file\n";
print "Output file:    \t$output_file\n\n";
print "Summary file:   \tsnpEff_summary.html\n\n";

print COMMAND_LOG "\n";
print COMMAND_LOG "Input file:     \t$vcf_file\n";
print COMMAND_LOG "Output file:    \t$output_file\n\n";
print COMMAND_LOG "Summary file:   \tsnpEff_summary.html\n\n";

printf COMMAND_LOG "\nRun time: \t%d days, %d hours, %d minutes and %d seconds\n",(gmtime $run_time)[7,2,1,0];

close COMMAND_LOG;

exit;

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

