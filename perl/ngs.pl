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

print "  NGS_pipeline_simple:    \tthis is based on the Broad pipeline but without the frills \n\n";
print "  NGS_pipeline_broad:     \tthis is based on the full Broad pipeline with all the extra steps \n\n";
print "  NGS_pipeline_original:  \tthis is Oliver's original pipeline \n\n\n";


print color 'bold white';

print "  I recommend you start with the simple pipeline.\n\n";
print "  Run it by typing:  perl /home/genetics/scripts/NGS_pipeline_simple.pl\n\n";
print "  Mike\n\n\n";


print color 'reset';