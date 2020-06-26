#!/usr/bin/perl -w

use strict;
use Getopt::Std ;
use File::Basename ;
use Term::ANSIColor;

################
# OMIM spider #
################

my $version						= "21";
my $API_key						= "y8QqoiJVQk6SO53fyeUjCg";  # update_required. The API key needs to be renewed every so often.

my $flank_size					= 60; # This can be changed. It sets the extra context" text that is printed to the output file
my $retnet_file					= "/home/genetics/information/retnet_gene_list.txt"; 
my $epilepsy_file				= "/home/genetics/information/epilepsy_gene_list.txt"; 
my $OMIA_file					= "/home/genetics/information/OMIA_gene_list.txt";

#Debugging
my $gene_to_debug				= "CDS2xxx";
my $UID_to_debug				= "113730xxx";
my $keyword_to_debug			= "xxxx";

#Files
my $gene_list_file				= "";
my $gene_list_file_cleaned		= ""; # Genes file cleaned up for varElect
my $keyword_list_file			= "";
my $gene_output_file			= "gene_output_OMIM.txt";
my $UID_output_file				= "UID_output_OMIM.txt";
my $output_file					= ""; # Results file
my $gene_results_output_file	= ""; # Results shown as a list of genes with the number of hits per gene
my $keywords_varElect_file		= ""; # Keywords file in varElect format


#Strings
my $gene_search_short_name		= "";
my $gene_search_long_name		= "";
my $gene_hit_short_name			= "";
my $gene_hit_long_name			= "";
my $gene_name_uid_line			= "";
my $gene_string					= "";
my $keyword						= "";
my $gene_line					= "";
my $UID_line					= "";
my $from_name					= "";
my $to_name						= "";
my $answer						= "";
my $command						= "";
my $UID_string					= "";
my $UID_symbol					= ""; # This is the asterisk or hash # at the start of the UID line
my $UID_type					= ""; # e.g. gene, disease etc
my $key							= "";
my $value						= "";
my $prefix						= "";
my $prefix_gene_list			= "";
my $prefix_keyword_list		= "";
my $last_output_string			= "x";
my $output_string				= "";
my $extracted_string			= "";
my $gene_line_forward			= ""; # Searching forward ingene output file
my $keyword_hits_string			= ""; # String of keyword hit counts for output gene file
my $retnet_gene					= ""; # retnet gene name
my $epilepsy_gene				= ""; # epilepsy gene name
my $OMIA_gene					= "";
my $OMIA_data_line				= "";
my $OMIA_phenotype				= "";

#Counters
my $gene_count					= 0;
my $retnet_gene_count			= 0;
my $retnet_gene_found_count		= 0;
my $epilepsy_gene_count			= 0;
my $epilepsy_gene_found_count	= 0;
my $OMIA_gene_count				= 0;
my $OMIA_gene_found_count		= 0;
my $keyword_count				= 0;

#Numbers
my $line_uid_pos				= 0;
my $UID_count					= 0;
my $right_angle_bracket_pos		= 0;
my $left_angle_bracket_pos		= 0;
my $space_pos					= 0;
my $semi_colon_pos				= 0;
my $a_tag_pos					= 0;
my $total_no_genes				= 0;
my $total_no_keywords			= 0;
my $item_count					= 0;
my $no_of_items					= 0;
my $no_of_genes					= 0;			
my $gene_line_count				= 0;
my $UID_line_count				= 0;
my $file_size					= 0;
my $key_word_found_count		= 0;
my $hit_point					= 0;
my $left_point					= 0;
my $right_point					= 0;
my $gene_output_array_size		= 0;
my $UID_output_array_size		= 0;
my $check_count					= 0;
my $found_count_for_this_gene	= 0;

#Boolean
my $keyword_found_in_UID_file 	= ""; # true or false
my $gene_list_comments_found	= "false";
my $keyword_list_comments_found	= "false";

#Arrays
my @gene_array							= ();
my @keyword_array						= ();
my @item_array							= ();
my @gene_output_array					= ();
my @UID_output_array					= ();
my @gene_search_short_name_results_array= ();
my @gene_search_long_name_results_array	= ();
my @gene_hit_short_name_results_array	= ();
my @gene_hit_long_name_results_array	= ();
my @found_count_for_this_gene_array		= ();
my @keyword_hits_array					= (); # Store the number of hits for each keyword (for current gene only)
my @keyword_hits_string_array			= ();


#Hashes
my %gene_hash					= ();
my %retnet_gene_hash					= ();
my %epilepsy_gene_hash					= ();
my %OMIA_gene_hash						= ();

print color 'reset';
print color 'bold magenta';

print "\n\n";
print "############################\n";
print color 'bold white';
print "      OMIM spider      \n";
print color 'bold magenta';
print "############################\n";
print "\n\n";

print color 'yellow';
print "Version $version\n\n";
print "  - Checks OMIM for a list of genes\n\n";
print "  - OMIM API key (valid from March 2017) is $API_key\n\n";


print color 'reset';

$gene_list_file = &get_file("What is the name of the file with a list of genes to check?",".txt");
$keyword_list_file = &get_file("What is the name of the file with a list of keywords to check?",".txt");
open(GENES,"$gene_list_file") or die "Can't open file $gene_list_file";


##################################################
# Extract clean list of genes from the gene list #
##################################################
 while ($gene_string = <GENES>) 
{
	chomp $gene_string;

	###########################################
	# Lines starting with a hash are comments #
	###########################################
	if (index ($gene_string,"#") == -1)
	{
		if (index($gene_string,"+") > -1)
		{
			@item_array = split (/\+/,$gene_string);

			$no_of_items = scalar @item_array;

			for ($item_count = 0; $item_count < $no_of_items; $item_count++)
			{
				$gene_count = $gene_count + 1;
				$gene_array[$gene_count] = $item_array[$item_count];
			}
		}
		else
		{
			$gene_count = $gene_count + 1;
			$gene_array[$gene_count] = $gene_string;
		}

	} # If doesn't start with a hash

	if (index ($gene_string,"#") > -1)
	{
		print "$gene_string\n";
		$gene_list_comments_found = "true";
	} # If this is a comment
}

close GENES;



####################################################################
# Make up a cleaned list of genes: No duplicates and no plusses    #
# If it is from vcf2excel via NGS SNP Viewer it may have these     #
####################################################################
$prefix_gene_list= &get_prefix($gene_list_file);
$gene_list_file_cleaned = $prefix_gene_list."_varElect.txt";

open(GENES_CLEANED,">$gene_list_file_cleaned") or die "Can't open file $gene_list_file_cleaned";

$no_of_genes = $gene_count;

print "\nList of genes to search\n";
print "=======================\n";

for ($gene_count = 1; $gene_count <= $no_of_genes; $gene_count++)
{
	$found_count_for_this_gene_array[$gene_count] = 0;
	$gene_hit_short_name_results_array[$gene_count] = "hit_gene_short_name_".$gene_count;
	$gene_hit_long_name_results_array[$gene_count] = "hit_gene_long_name_".$gene_count;


	if (not defined $gene_hash{$gene_array[$gene_count]})
	{
		$gene_hash{$gene_array[$gene_count]} = "defined";
		print "$gene_count\t$gene_array[$gene_count]\n";
	}
}


print "\n";
$gene_count = 0;
while (($gene_search_short_name, $value) = each %gene_hash) 
{
	$gene_count = $gene_count + 1;
	print GENES_CLEANED "$gene_search_short_name\n";
}

close GENES_CLEANED;

$no_of_genes = $gene_count;
print "Total number of genes: $no_of_genes\n";


open(KEYWORDS,"$keyword_list_file") or die "Can't open file $keyword_list_file";
$keyword_count =0;

$prefix = &get_prefix("$keyword_list_file");
$keywords_varElect_file = $prefix."_varElect.txt";

open(KEYWORDS_VARELECT,">$keywords_varElect_file") or die "Can't open file $keywords_varElect_file";


#############################
# Read list of keywords     #
#############################
print "\nList of keywords\n";
print "================\n";


while ($keyword = <KEYWORDS>) 
{
	chomp $keyword;

	###########################################
	# Lines starting with a hash are comments #
	###########################################
	if (index ($keyword,"#") == -1)
	{
		$keyword_count = $keyword_count + 1;
		$keyword_array[$keyword_count] = $keyword;

		$keyword =~ s/^\s+|\s+$//g;

		if (($keyword ne "") && ($keyword_count > 1))
		{
			print KEYWORDS_VARELECT ",\n";
		}
		print  "$keyword";
		print KEYWORDS_VARELECT "\"$keyword\"";
	}
	if (index ($keyword,"#") > -1)
	{
		print "$keyword\n";
		$keyword_list_comments_found = "true";
	} # If this is a comment

} # Keywords file loop


close KEYWORDS;

$total_no_keywords = $keyword_count;
close KEYWORDS_VARELECT;


###########################################
# Make name for output file               #
###########################################
$prefix_gene_list = &get_prefix($gene_list_file);
$prefix_keyword_list = &get_prefix($keyword_list_file);


$output_file = "OMIM_search_".$prefix_gene_list."_".$prefix_keyword_list."_Keywords.out";
$gene_results_output_file = "OMIM_search_".$prefix_gene_list."_".$prefix_keyword_list."_Genes.out";


print "\nKeywords output file: \t$output_file\n";
print "Genes output file:    \t$gene_results_output_file\n";

print "\nPress 'return' to continue with OMIM searches\n\n";
$answer=<STDIN>;


&read_retnet_gene_list;
print "$retnet_gene_count genes read from RETNET file\n\n";

&read_epilepsy_gene_list;
print "$epilepsy_gene_count genes read from EPILEPSY file\n\n";

&read_OMIA_gene_list;
print "$OMIA_gene_count genes read from OMIA file\n\n";





###########################################
# Make an output file                     #
###########################################
open(OUTPUT,">$output_file") or die "Can't open file $output_file";

print OUTPUT "Gene (search gene)\tGene (hit gene)\tLong name (hit gene)\tUID type\tUID\tKeyword\tText around hit\n";


$gene_count = 0;

while (($gene_search_short_name, $value) = each %gene_hash) 
{
	$gene_count = $gene_count + 1;
	$found_count_for_this_gene = 0; # <== set this to zero at the start of each gene

	############################
	# Clear keyword hits array #
	############################
	for ($keyword_count = 0; $keyword_count <=100; $keyword_count++)
	{
		$keyword_hits_array[$keyword_count] = 0;
	}
	

	print "===========================================================================\n";
	print "Gene $gene_count/$no_of_genes: $gene_search_short_name\n\n";

	#####################################################################
	# Access gene and write page to a temporary output file GENE_OUTPUT #
	#####################################################################

	########################################
	# Use wget with API key to create file #
	########################################
	&run_unix_command("wget \"https://api.omim.org/api/entry/search?search=$gene_search_short_name&start=0&limit=20&apiKey=$API_key\"  -q -O $gene_output_file", "no_print");


	#wget "http://api.omim.org/api/entry/search?search=duchenne&start=0&limit=20&apiKey=ED9CB85D8738D707FFCDEF22CD3D3DD0E3455A5A"

	##DEBUGGING ##
	if ($gene_search_short_name eq $gene_to_debug)
	{
		print "Stopping at gene: $gene_search_short_name\n";
		&pause;
	}


	#######################################################
	# Open gene output file and load into an array, then  #
	# search through it for 'mimNumber' and other things  #
	#######################################################

	open(GENE_OUTPUT,"$gene_output_file") or die "Can't open file $gene_output_file";
	$file_size = -s $gene_output_file;

	@gene_output_array = <GENE_OUTPUT> ; 
	close GENE_OUTPUT ;
	$gene_output_array_size = scalar @gene_output_array;


	#######################################################
	# Loop through the array from the gene output file    #
	#######################################################
	if ($gene_search_short_name eq $gene_to_debug)
	{
		print"\n\n================ Looking through Gene Output file for $gene_search_short_name =========================\n";
	}

	for ($gene_line_count=0;$gene_line_count < $gene_output_array_size; $gene_line_count++)
	{
		$gene_line = $gene_output_array[$gene_line_count];
		chomp $gene_line;

		# Look for <mimNumber>606966</mimNumber>

		#############################################
		# If a UID is found on the line (mimNumber) #
		#############################################
		if (index ($gene_line,"<mimNumber>") > -1)
		{
			
			$gene_name_uid_line = "";

			$UID_count = $UID_count + 1;

			$line_uid_pos = index($gene_line,"<mimNumber>");

			$right_angle_bracket_pos = index($gene_line,">",$line_uid_pos);

			$left_angle_bracket_pos = index($gene_line,"<",$right_angle_bracket_pos);

			$UID_string = substr($gene_line,$right_angle_bracket_pos + 1,$left_angle_bracket_pos - $right_angle_bracket_pos - 1);


			############################################################
			# Search forward a couple of lines for the preferred Title #
			############################################################
			$gene_hit_short_name = "SHORT_HIT";
			$gene_hit_long_name = "LONG_HIT";

			for ($check_count=1;$check_count<=3; $check_count++)
			{
				if (($gene_line_count + $check_count) < $gene_output_array_size)
				{
					$gene_line_forward =  $gene_output_array[$gene_line_count + $check_count];

					chomp $gene_line_forward;

					# <preferredTitle>CLAUDIN 12; CLDN12</preferredTitle>
					# <preferredTitle>HOLOCARBOXYLASE SYNTHETASE DEFICIENCY</preferredTitle> <== Note no gene short name
					# 0.........x.........x.........x.........x.........x   34-26 
					if (index($gene_line_forward,"<preferredTitle>") > -1)
					{
						$right_angle_bracket_pos = index($gene_line_forward,">");
						$left_angle_bracket_pos = index($gene_line_forward,"<",2);
						$semi_colon_pos = index($gene_line_forward,";",2);

						##########################################################
						# If there is a semi-colon (and presumably a short name) #
						##########################################################
						if ($semi_colon_pos > 0)
						{
							$gene_hit_long_name = substr($gene_line_forward,$right_angle_bracket_pos + 1,$semi_colon_pos - $right_angle_bracket_pos - 1);
							$gene_hit_short_name = substr($gene_line_forward,$semi_colon_pos + 2,$left_angle_bracket_pos - $semi_colon_pos -2);
						}
						else
						{
							$gene_hit_long_name = substr($gene_line_forward,$right_angle_bracket_pos + 1,$left_angle_bracket_pos - $right_angle_bracket_pos - 1);
							$gene_hit_short_name = "";
						}
						
						if ($gene_search_short_name eq $gene_hit_short_name)
						{
							$gene_search_long_name = $gene_hit_long_name;
						}


						if ($gene_search_short_name eq $gene_to_debug)
						{
							print "\n^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n";
							print "  Gene               \t$gene_search_short_name\n";
							print "  UID_string         \t$UID_string\n";
							print "  gene_line_forward: \t$gene_line_forward\n";
							print "  gene_long_name:    \t>$gene_hit_long_name<\n";
							print "  gene_short_name:   \t>$gene_hit_short_name<\n";

							&pause;
						}

						

					} # if word <preferredTitle> is present in line

				}	

			} # check loop (looking ahead)


			


			########################################
			# What do symbols before the UID mean? #
			########################################

			# An asterisk (*) before an entry number indicates a gene.
			# A number symbol (#) before an entry number indicates that it is a descriptive entry, usually of a phenotype, and does not represent a unique locus. The reason for the use of the number symbol is given in the first paragraph of the entry. Discussion of any gene(s) related to the phenotype resides in another entry(ies) as described in the first paragraph.			
			# A plus sign (+) before an entry number indicates that the entry contains the description of a gene of known sequence and a phenotype.
			# A percent sign (%) before an entry number indicates that the entry describes a confirmed mendelian phenotype or phenotypic locus for which the underlying molecular basis is not known.
			# No symbol before an entry number generally indicates a description of a phenotype for which the mendelian basis, although suspected, has not been clearly established or that the separateness of this phenotype from that in another entry is unclear.
			# A caret (^) before an entry number means the entry no longer exists because it was removed from the database or moved to another entry as indicated.
			# See also the description of symbols used in the disorder column of the OMIM Gene Map and Morbid Map.


			print "  UID: $UID_string\t$gene_hit_long_name\n";

			##########################################
			# Remove white spaces from string if any #
			##########################################
			$UID_string =~ s/\s+//g;
			

			######################################################
			# Use wget with API key to create UID output  file   #
			# This is what you get when you click on an UID line #
			######################################################
			&run_unix_command("wget \"http://api.omim.org/api/entry?mimNumber=$UID_string&include=all&apiKey=$API_key\"  -q -O $UID_output_file","no_print");


			######################################################
			#     Open UID output file and search for keywords   #
			######################################################


			$UID_line_count = 0;

			$file_size = -s $UID_output_file;

			if ($UID_string eq $UID_to_debug)
			{
				print "Pause at UID $UID_to_debug\n";
				&pause;
			}

			open(UID_OUTPUT,"$UID_output_file") or die "Can't open file $UID_output_file";

			
			if ($gene_search_short_name eq $gene_to_debug)
			{
				print"\n\n   ================ Looking through UID Output file $UID_string ($gene_hit_short_name) for $gene_search_short_name =========================\n\n";
			}

			#######################################################
			# Open UID output file and load into an array, then   #
			# search through it for keywords                      #
			# Search whole file with each keyword                 #
			#######################################################
			@UID_output_array = <UID_OUTPUT> ; 
			close UID_OUTPUT ;
			$UID_output_array_size = scalar @UID_output_array;

			##############################
			# Set Boolean found as false #
			##############################
			$keyword_found_in_UID_file = "false";

			for ($keyword_count = 1; $keyword_count <= $total_no_keywords; $keyword_count++)
			{
				$keyword = $keyword_array[$keyword_count];
				

				######################################
				# Loop through keywords for this UID #
				######################################
				if ($gene_search_short_name eq $gene_to_debug)
				{
					print"   Keyword ($keyword_count/$total_no_keywords) $keyword for Gene: $gene_search_short_name and UID: $UID_string ($gene_hit_short_name)\n";
				}


				#########################################################
				# Search through whole UID file array with this keyword #
				#########################################################
				for ($UID_line_count=0;$UID_line_count < $UID_output_array_size; $UID_line_count++)
			{
					$UID_line = $UID_output_array[$UID_line_count];
				chomp $UID_line;


				######################################
				# Get prefix e.g. asterisk, hash etc #
				######################################
				if (index($UID_line,"<prefix>") > -1)
				{
					$right_angle_bracket_pos = index($UID_line,">");
					$UID_symbol = substr($UID_line,$right_angle_bracket_pos +1,1);

					$UID_type = "phenotype no known mendelian basis";
					if ($UID_symbol eq "*"){$UID_type = "gene"}
					if ($UID_symbol eq "#"){$UID_type = "phenotype no locus"}
					if ($UID_symbol eq "+"){$UID_type = "gene and phenotype"}
					if ($UID_symbol eq "%"){$UID_type = "mendelian phenotype"}
					if ($UID_symbol eq "^"){$UID_type = "removed"}

					} # if <prefix>


					################################################################
					# If the keyword is found in the line from the UID output file #
					################################################################
					if (index(lc $UID_line,lc $keyword) > -1)
					{
						
						$keyword_found_in_UID_file = "true";
						$key_word_found_count = $key_word_found_count + 1;
						$found_count_for_this_gene = $found_count_for_this_gene + 1;  # <== keyword count for just this gene
						
						########################################
						# Store the number of hits per keyword #
						########################################
						$keyword_hits_array[$keyword_count] = $keyword_hits_array[$keyword_count] + 1;


						########################################
						# Get some surrounding text #
						########################################
						$hit_point = index(lc $UID_line,lc $keyword);
						$left_point = $hit_point - $flank_size;
						$right_point = $hit_point + $flank_size;

						if ($left_point < 1) {$left_point = 1};
						if ($right_point > length($UID_line)) {$right_point = length($UID_line)};

                        $extracted_string = substr($UID_line, $left_point, $right_point - $left_point + 1);


						#################################################################
						# Make up output string to check against previous output values #
						# If it's the same, don't bother to print to OUTPUT again       #
						#################################################################
						$output_string = $gene_search_short_name.$gene_name_uid_line.$UID_string.$keyword;

						if ($gene_search_short_name eq $gene_to_debug)
						{
							print "\n\t~~~~ Found a hit in the UID output file array ~~~~ Keyword:$keyword\tHit gene:$gene_hit_short_name\tGene from list:$gene_search_short_name  \n\n";
							print "\tUID_line:        \t$UID_line\n\n";
							print "\tKeyword:         \t$keyword\n";
							print "\tThis MIGHT BE written to output file\n\n";
							print "\tGENE:            \t$gene_search_short_name\n";
							print "\tUID_type:        \t$UID_type\n";
							print "\tgene_hit_long_name:  \t$gene_hit_long_name\n";
							print "\tgene_hit_short_name: \t$gene_hit_short_name\n\n";

							$answer=<STDIN>;
						}


						###########################################
						# Phenotype no locus - no gene name short #
						###########################################
						if ($UID_type eq "phenotype no locus")
						{
							$gene_hit_short_name = "No locus";
						}


						##################################################################
						# Write to output file if not the same as the last output string #
						##################################################################
						if ($output_string ne $last_output_string)
						{
							if ($gene_search_short_name ne $gene_hit_short_name)
							{
								if ($UID_type eq "gene"){ $UID_type = "related gene"; }
							}
							if ($gene_search_short_name eq $gene_hit_short_name)
							{
								if ($UID_type eq "gene"){ $UID_type = "same gene"; }
							}
							print color 'bold magenta';
							print "    >>>>>>>> $gene_search_short_name\tKeyword: $keyword\t$gene_hit_long_name\tGene hit: $gene_hit_short_name\tTotal keywords found: $key_word_found_count\n";
							print color 'reset';

							print OUTPUT "$gene_search_short_name\t$gene_hit_short_name\t$gene_hit_long_name\t$UID_type\t$UID_string\t$keyword\t$extracted_string\n";

						}## if not the same as the last


						$last_output_string = $output_string;

					} # If keyword is found

				} # UID line array loop


				####################################################
				# Now we have searched through this whole UID file #
				# What do we need to store from this?              #
				####################################################

				if ($keyword_found_in_UID_file eq "false")
				{
							if ($gene_search_short_name eq $gene_to_debug)
							{
						print "    Keyword $keyword was not found in the UID file from $UID_string\n\n";
					}
				}


				if ($keyword_found_in_UID_file eq "true")
				{
					if (($gene_search_short_name eq $gene_to_debug) || ($gene_hit_short_name eq $gene_to_debug))
					{
						print "    Keyword $keyword WAS found in the UID file from $UID_string\n\n";
						print "       > So store this name for the hit gene: $gene_hit_short_name\n\n";
							}

					#####################################################################
					# We want to note if the gene name is the same as the hit gene name #
					#####################################################################
					#print "    Keyword $keyword WAS found in the UID file from $UID_string\n\n";
					#print "       > So store this name for the hit gene: Short: $gene_hit_short_name   Long: $gene_hit_long_name  Search: $gene_search_short_name\n\n";


					$gene_hit_short_name_results_array[$gene_count] = $gene_hit_short_name;
					$gene_hit_long_name_results_array[$gene_count] = $gene_hit_long_name;

					#&pause;

				} # if (keyword found

			} # keyword loop

		} #if (index ($gene_line,"link_uid") > 0)



		#################################################
		# Make a string for the numbers of keyword hits #
		# This is for the Gene summary output file      #
		#################################################
		$keyword_hits_string = "";
		for ($keyword_count = 1; $keyword_count <= $total_no_keywords; $keyword_count++)
		{
			if ($keyword_hits_array[$keyword_count] > 0)
			{
				$keyword_hits_string = $keyword_hits_string.$keyword_hits_array[$keyword_count]."\t";
			}
			else
			{
				$keyword_hits_string = $keyword_hits_string."\t";
			}
		}


		#########################
		# Summary for this gene #
		#########################
		$found_count_for_this_gene_array[$gene_count] = $found_count_for_this_gene;
		$gene_search_short_name_results_array[$gene_count] = $gene_search_short_name;
		$gene_search_long_name_results_array[$gene_count] = $gene_search_long_name;

		#print "Gene search short name: $gene_search_short_name\n";
		#print "Gene search long name:  $gene_search_long_name\n";

		#&pause;

		################################################################
		# This stores the name of the gene hit                         #
		# However this is the last gene looked at so shouldn't be here #
		################################################################
		#$gene_hit_short_name_results_array[$gene_count] = $gene_hit_short_name;
		#$gene_hit_long_name_results_array[$gene_count] = $gene_hit_long_name;

		$keyword_hits_string_array[$gene_count] = $keyword_hits_string;

	} # working down GENE OUTPUT file array looking for UID tag


} # gene_hash loop

close OUTPUT;

$total_no_genes = $gene_count;

print "\n\nTotal no of genes: $total_no_genes\n\n";

&print_message("No of hits for individual genes","message");

open(GENE_RESULTS_OUTPUT,">$gene_results_output_file") or die "Can't open file $gene_results_output_file";


#################
# Print headers #
#################
print GENE_RESULTS_OUTPUT "Gene (search gene)\tGene long name\tNumber of hits";

for ($keyword_count = 1; $keyword_count <= $total_no_keywords; $keyword_count++)
{
	print GENE_RESULTS_OUTPUT "\t$keyword_array[$keyword_count]";
}

print GENE_RESULTS_OUTPUT "\tRetnet?\tEpilepsy?\tOMIA?\n";


######################
# Make list of genes #
######################
$retnet_gene_found_count = 0;

for($gene_count = 1; $gene_count <= $total_no_genes; $gene_count++)
{
	if ($found_count_for_this_gene_array[$gene_count] > 0)
	{
		print "$gene_search_short_name_results_array[$gene_count]:\t$gene_search_long_name_results_array[$gene_count]\t$found_count_for_this_gene_array[$gene_count]\n";
		print GENE_RESULTS_OUTPUT "$gene_search_short_name_results_array[$gene_count]\t$gene_search_long_name_results_array[$gene_count]\t$found_count_for_this_gene_array[$gene_count]\t$keyword_hits_string_array[$gene_count]";

		#########################################
		# Check if gene is on the RETNET list   #
		#########################################
		$gene_search_short_name = $gene_search_short_name_results_array[$gene_count];

		if (defined $retnet_gene_hash{$gene_search_short_name})
		{
			$retnet_gene_found_count = $retnet_gene_found_count + 1;
			print "\t Retnet gene $gene_search_short_name_results_array[$gene_count] found\n";
			print GENE_RESULTS_OUTPUT "RET\t";
		}
		else
		{
			print GENE_RESULTS_OUTPUT "\t";
		}


		#########################################
		# Check if gene is on the EPILEPSY list #
		#########################################
		$gene_search_short_name = $gene_search_short_name_results_array[$gene_count];

		if (defined $epilepsy_gene_hash{$gene_search_short_name})
		{
			$epilepsy_gene_found_count = $epilepsy_gene_found_count + 1;
			print "\t Epilepsy gene $gene_search_short_name_results_array[$gene_count] found\n";
			print GENE_RESULTS_OUTPUT "EPI\t";
	}
		else
		{
			print GENE_RESULTS_OUTPUT "\t";
		}


		#########################################
		# Check if gene is on the OMIA list     #
		#########################################
		$gene_search_short_name = $gene_search_short_name_results_array[$gene_count];

		if (defined $OMIA_gene_hash{$gene_search_short_name})
		{
			$OMIA_gene_found_count = $OMIA_gene_found_count + 1;
			print "\t OMIA gene $gene_search_short_name_results_array[$gene_count] found\n";
			print GENE_RESULTS_OUTPUT "$OMIA_gene_hash{$gene_search_short_name}\t";
		}
		else
		{
			print GENE_RESULTS_OUTPUT "\t";
		}

		#########################################
		# Print end of line                     #
		#########################################
		print GENE_RESULTS_OUTPUT "\n";

	}
}

close GENE_RESULTS_OUTPUT;

&print_message("Finished searching list of genes against OMIM","message");

print "Gene list file:       \t$gene_list_file\n";
print "Keywords list file:   \t$keyword_list_file\n\n";

print "Output file:          \t$output_file\n";
print "Output file (genes):  \t$gene_results_output_file\n\n\n";

print "Files for use with varElect (GeneCards): \n\n";
print "  Gene list file:     \t$gene_list_file_cleaned\n";
print "  Keywords list file: \t$keywords_varElect_file\n\n\n";

print "Number of genes on other lists:\n\n";

print "\tRetnet:   \t$retnet_gene_found_count\n";
print "\tEpilepsy: \t$epilepsy_gene_found_count\n";
print "\tOMIA:     \t$OMIA_gene_found_count\n\n";

exit;

############################################################################

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

		if ($gene_list_file eq "test_gene_list.txt"){$_file_name = "test_keywords.txt"}

		# User types 'ls'
		if (($_file_name eq "ls") || ($_file_name eq "")) {print "\n"; system ("ls *"."$_file_type  --ignore='*OMIM_search*'")}

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

}#

#############################################
# Subroutine to execute unix command        #
# Argument 1: command to execute            #
#############################################
sub run_unix_command
{
	my $_unix_command 		= "";
	my $_print_option		= "";

	$_unix_command = $_[0];
	$_print_option = $_[1];

	if ($_print_option ne "no_print") { print("$_unix_command\n\n"); }
	system("$_unix_command");
}



###########################################################
# Subroutine to get filename before file type suffix      #
# (i.e. prefix) (e.g. to get "test" from "test.fasta")    #
# New version using rindex rather than index.  This means #
# it can deal with files like filename.something.vcf      #
###########################################################

sub get_prefix
{
	my $_filename 	= "";
	my $_dot_pos	= 0;

	$_filename = $_[0];
	$_dot_pos = rindex($_filename,".");

	if (rindex($_filename,".") > 0)
	{
		$_filename = substr($_filename, 0, rindex($_filename,"."));
	}
	if (rindex($_filename,".") == -1)
	{
		$_filename = $_filename;
	}

	$_filename = $_filename;

} # get_prefix

##########################################################
# Subroutine to pause until user hits 'return'           #
##########################################################
sub pause
{
	my $_answer = "";
	print "\n Press RETURN to continue\n";
	$_answer=<STDIN>;
}


sub read_retnet_gene_list
{
	$retnet_gene_count = 0;

	system("dos2unix $retnet_file");

	open(RETNET,"$retnet_file") or die "Can't open file $retnet_file";

	while ($retnet_gene = <RETNET>) 
	{
		chomp $retnet_gene;

		if (not defined $retnet_gene_hash{$retnet_gene})
		{
			$retnet_gene_hash{$retnet_gene} = $retnet_gene;
			$retnet_gene_count = $retnet_gene_count + 1;
		}
	}

	close RETNET;

}#read_retnet_gene_list


sub read_epilepsy_gene_list
{
	$epilepsy_gene_count = 0;

	system("dos2unix $epilepsy_file");

	open(EPILEPSY,"$epilepsy_file") or die "Can't open file $epilepsy_file";

	while ($epilepsy_gene = <EPILEPSY>) 
	{
		chomp $epilepsy_gene;

		if (not defined $epilepsy_gene_hash{$epilepsy_gene})
		{
			$epilepsy_gene_hash{$epilepsy_gene} = $epilepsy_gene;
			$epilepsy_gene_count = $epilepsy_gene_count + 1;
		}
	}

	close EPILEPSY;

}#read_epilepsy_gene_list

sub read_OMIA_gene_list
{
	

	$OMIA_gene_count = 0;

	system("dos2unix $OMIA_file");

	open(OMIA,"$OMIA_file") or die "Can't open file $OMIA_file";

	while ($OMIA_data_line = <OMIA>) 
	{
		chomp $OMIA_data_line;
		@item_array = split (/\t/,$OMIA_data_line);
		$OMIA_gene = $item_array[0];
		$OMIA_phenotype = $item_array[1];
		if (not defined $OMIA_gene_hash{$OMIA_gene})
		{
			$OMIA_gene_hash{$OMIA_gene} = $OMIA_phenotype;
			$OMIA_gene_count = $OMIA_gene_count + 1;
		}
	}

	close OMIA;

}#read_OMIA_gene_list
