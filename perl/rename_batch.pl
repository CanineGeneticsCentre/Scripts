#!/usr/bin/perl -w

################
# rename batch #
################


my $list_file			= "";
my $single_line			= "";
my $from_name			= "";
my $to_name				= "";
my $answer				= "";
my $command				= "";
my $version				= "3";

my @item				= ();

print "\n\n";
print "################\n";
print "# rename batch #\n";
print "################\n\n";

print "Version $version\n\n";

print "Input list of files in two tab-separated columns, first is FROM name and second is TO name\n\n";

$list_file = <STDIN>;
chomp $list_file;


if (-e $list_file)
{
	system ("dos2unix $list_file");
}

####################################################
# Open the list file to get the list of file names #
####################################################
open (LIST, "$list_file") || die "Cannot open $list_file";
$list_count=1;
$no_of_files=0;

print "\n\nList of files to rename:\n\n";

while ($single_line = <LIST> ) 
{
	chomp $single_line;

	@item=split(/\t/,$single_line);
	$from_name = $item[0];
	$to_name = $item[1];
	
	print "$from_name \t-->\t$to_name\n";
}

close LIST;

print "\nProceed to rename all these files? (Q to quit)   ";
$answer = <STDIN>;
chomp $answer;
if (lc $answer eq "q"){exit;}

print "\n\nRunning program\n\n";

#########################################
# Open the list file to do the renaming #
#########################################
open (LIST, "$list_file") || die "Cannot open $list_file";
$list_count=1;
$no_of_files=0;
while ($single_line = <LIST> ) 
{
	chomp $single_line;

	@item=split(/\t/,$single_line);
	$from_name = $item[0];
	$to_name = $item[1];
	
	print "From: $from_name    To: $to_name\n";
	
	################################
	# Check if file already exists #
	################################
	if (! -e $from_name)
	{
		print "\n>>>>>   File $from_name was not found    <<<<<\n\n";
	}
	if (-e $from_name)
	{
		if (-e $to_name)
		{
			print "\n\n";
			print "##########################\n";
			print "#        WARNING!!!!     #\n";
			print "##########################\n\n";
			
			print "File $to_name already exists\n\n";
			print "Are you sure you want to overwrite it? (y/n/q)   ";
			$answer=<STDIN>;
			chomp $answer;
			
			if (lc ($answer) eq "q") {print "\n\nExiting program because of q\n\n";exit;}
		
			if (lc ($answer) eq "y") 
			{
				$command = "mv $from_name $to_name";
				print "$command\n";
				system ("$command");
			}
		
		} # if to_name exists
		else
		{
			$command = "mv $from_name $to_name";
			print "$command\n";
			system ("$command");
		}
	
		
	} # if from_name exists
	
}

close LIST;

$no_of_files=$list_count - 1;

print "\n\nFinished program\n\n";