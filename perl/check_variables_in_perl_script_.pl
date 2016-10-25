#!/usr/bin/perl -w

use strict;
use Term::ANSIColor;

# Constants
my $version 				= "1";
my $file_path				= "/home/genetics/test_scripts/";

#Variables
my $perl_file				= "";
my $perl_file_edited		= "";
my $single_line				= "";
my $single_word				= "";
my $answer					= "";
my $prefix					= "";
my $not_used_again			= ""; # true or false

my $array_size				= 0;
my $word_count				= 0;
my $variable_count			= 0;
my $line_count 				= 0;
my $no_of_lines				= 0;
my $no_of_variables			= 0;
my $no_of_variables_unused	= 0;

my @word_array				= (); # Split line into words
my @perl_file_array			= (); # whol eperl file read into this array
my @variable_array			= ();
my @variable_used_array		= (); # stores "yes" or "no" with each variable
my @variable_line_array		= (); # stores line that variable was declared on

print color 'reset';
print color 'bold magenta';
print "\n\n";
print "############################################\n";
print color 'bold white';
print "      check_variables_in_perl_script   \n";
print color 'bold magenta';
print "############################################\n\n";

print color 'yellow';
print "Version $version\n\n";
print "    - This program looks through a PERL file and marks all variables that are declared but never used again\n\n";
print "      An edited version of the PERL file is written with these varaibles commented out with '#####'\n\n";

print color 'reset';


############################
# Get the file name prefix #
############################

until (-e $perl_file)
{
	print "PERL file to check:    ";
	$perl_file = <STDIN>;
	chomp $perl_file;

	if ($perl_file ne "ls")
	{
		$perl_file = $file_path.$perl_file;
		if (index($perl_file,".pl") == -1){$perl_file = $perl_file.".pl"}
	}

	if ($perl_file eq "ls"){print "\n";system ("ls $file_path*.pl");print "\n"}
	if ($perl_file ne "ls"){if (! -e $perl_file){print "\n\n>>>>>>>>  File $perl_file not found.  Try again.  <<<<<<<<\n\n";}}
}


####################################
# Create name for edited PERL file #
####################################
$prefix = &get_prefix($perl_file);
$perl_file_edited = $prefix."_edited.pl";

##################
# Open PERL file #
##################
$single_word = "X";

open (PERL, $perl_file)|| die "Cannot open input file: $perl_file";
@perl_file_array = <PERL>;
close PERL;
$no_of_lines = scalar (@perl_file_array);

open (PERL_EDITED, ">$perl_file_edited")|| die "Cannot open input file: $perl_file_edited";

###############################
# Look for declared variables #
###############################
for ($line_count = 1; $line_count <=$no_of_lines; $line_count++)
{
	&chomp_all ($perl_file_array[$line_count - 1]);
	$single_line = $perl_file_array[$line_count - 1];

	if ($single_line ne "")
	{
		@word_array=split(/\s+/,$single_line);

		$array_size = scalar (@word_array);

		if ($array_size > 0 )
		{
			if($word_array[0] eq "my")
			{
				if (index($word_array[1],"\$") == 0)
				{
					$variable_count = $variable_count + 1;
					$variable_array[$variable_count] = $word_array[1];
					$variable_line_array[$variable_count] = $line_count;
				}
			}
		} # if array_size > 0
	}
}

$no_of_variables = $variable_count;


#####################################
# Check if variables are used again #
#####################################

print "\nVariables not used after initial declaration\n\n";

for ($variable_count = 1; $variable_count <=$no_of_variables; $variable_count++)
{
	$variable_used_array[$variable_count] = "no";

	for ($line_count = 1; $line_count <=$no_of_lines; $line_count++)
	{
		$single_line = $perl_file_array[$line_count - 1];
		if (index($single_line,"my") != 0)
		{
			if (index($single_line,$variable_array[$variable_count]) > -1)
			{
				$variable_used_array[$variable_count] = "yes";
				$variable_line_array[$variable_count] = -9;
			}
		} # if no my at start

	} # lines loop

	if ($variable_used_array[$variable_count] ne "yes")
	{
		print "\t$variable_array[$variable_count]\n";
		$no_of_variables_unused = $no_of_variables_unused + 1;
	}
} # variables loop

print "\n\n";

###############################
# Write edited file           #
###############################

for ($line_count = 1; $line_count <=$no_of_lines; $line_count++)
{
	$single_line = $perl_file_array[$line_count - 1];

	############################################
	# Check if any variables were on this line #
	############################################
	$not_used_again = "false";
	for ($variable_count = 1; $variable_count <=$no_of_variables; $variable_count++)
	{
		if ($variable_line_array[$variable_count] == $line_count)
		{
			$not_used_again = "true";
		}
	}

	if ($not_used_again eq "true")
	{
		print PERL_EDITED "###DELETE### $single_line\n";
	}
	else
	{
		print PERL_EDITED "$single_line\n";
	}

}

close PERL_EDITED;


print "Number of variables found:    \t$no_of_variables\n";
print "Number of variables unused:   \t$no_of_variables_unused\n\n";

print "Edited PERL script:           \t$perl_file_edited\n\n";

exit;

##########################################################
# Subroutine to remove both LF and CR from end of string #
##########################################################
sub chomp_all
{

foreach (@_) {s/\n//g}  
foreach (@_) {s/\r//g}  

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