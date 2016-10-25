#!/usr/bin/perl -w

use strict;
use Term::ANSIColor;
#################################
# Tell user about NGS pipelines #
#################################

print "\n\n\n\n";

print color 'bold magenta';

print "\n\n";
print "                 ##################################################\n";
print color 'bold white';
print "                                     NGS Pipelines                 \n";
print " \n";
print "                      NGS data analysis pipelines for samba64\n";
print color 'reset';
print color 'bold magenta';
print "                 ##################################################\n";
print "\n\n";


print color 'bold white';

print "  There are three NGS pipelines available:\n\n\n";

print color 'yellow';

print "    fastq2bam:              \tconverts FASTQ files to BAM files \n\n";
print "    bam2vcf:                \tconverts BAM files to VCF files listing SNPs and Indels\n\n";
print "    fastq2vcf:              \tconverts FASTQ files to VCF files (i.e. the above two combined)\n\n\n\n";

print color 'bold white';

print "  There are some versions for analysing RNA sequence:\n\n\n";

print color 'yellow';

print "    fastq2bam_RNA:          \tconverts FASTQ files to BAM files \n\n";
print "    bam2vcf_RNA:            \tconverts BAM files to VCF files listing SNPs and Indels\n\n\n";



print color 'bold white';
#print color 'bold magenta';
print "Run them by typing the name of the program (or perl /home/genetics/scripts/name.pl)\n\n\n";


print color 'reset';