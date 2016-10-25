#!/usr/bin/perl

#####################################################
#                  ind2phenoQTL                     #
#                                                   #
#  This program gets the phenotypes from the        #
#  example.ind file and makes an example.pheno file #
#  with the phenotypes only.                        #
#                                                   #
#  For QTL data the phenotypes are written with one #
#  phenotype on each line rather than along a       #
#  single line as is the case for binary traits.    #
#                                                   #
#####################################################



##########################################
#  In and out files in the two arguments #
#                                        #
#     >ind2phenoQTL infile outfile       #
#                                        #
##########################################
$in = $ARGV[0]; # .ind file
$out = $ARGV[1]; # .pheno file

open(IN,$in) || die "Cannot open input file: $in";
open(OUT,">$out") || die "Cannot open input file: $out";

print "\nind2phenoQTL -i $in -o $out\n\n";

####################################
# Work down lines in the .ind file #
####################################
while($line = <IN>)
{
	#########################################################
	# if line contains the string 'Case' then print a 2     #
	# if line contains the string 'Control' then print a 1  #
	# Anything else print it as it is.                      #
	#########################################################

	# split line at any number of spaces or tab #
	chomp $line;
	@item=split(/ +|\t/,$line);

	$id = $item[1];
	$gender = $item[2];
	$value = $item[3];

  	if ($line =~ /Case/) 
	{ 
		print OUT ("2.0\n"); 
	}

  	elsif ($line =~ /Control/)
	{ 
		print OUT ("1.0\n"); 
	}

  	else 
	{ 
		print OUT "$value\n";
	}
}

close IN;
close OUT;




