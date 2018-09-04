#!/bin/bash

ENS=93
mkdir -p e${ENS}/download
cd e${ENS}/download

##DOWNLOAD DATA
for i in `seq 1 38`; do wget ftp://ftp.ensembl.org/pub/release-${ENS}/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna_sm.chromosome.$i.fa.gz; done

wget ftp://ftp.ensembl.org/pub/release-${ENS}/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna_sm.chromosome.MT.fa.gz
wget ftp://ftp.ensembl.org/pub/release-${ENS}/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna_sm.chromosome.X.fa.gz
wget ftp://ftp.ensembl.org/pub/release-${ENS}/fasta/canis_familiaris/dna/Canis_familiaris.CanFam3.1.dna_sm.nonchromosomal.fa.gz

wget ftp://ftp.ensembl.org/pub/release-${ENS}/gtf/canis_familiaris/Canis_familiaris.CanFam3.1.${ENS}.gtf.gz

cd ../

##PARSE DATA
#Sequence
for i in `seq 1 38`; do zcat download/Canis_familiaris.CanFam3.1.dna_sm.chromosome.$i.fa >> canfam3.fasta; done
zcat download/Canis_familiaris.CanFam3.1.dna_sm.chromosome.MT.fa >> canfam3.fasta
zcat download/Canis_familiaris.CanFam3.1.dna_sm.chromosome.X.fa >> canfam3.fasta
zcat download/Canis_familiaris.CanFam3.1.dna_sm.nonchromosomal.fa.gz >> canfam3.fasta

samtools faidx canfam3.fasta; 
/opt/bwa/bwa index -a bwtsw canfam3.fasta
#CreateSequenceDictionary REFERENCE=canfam3.fasta OUTPUT=canfam3.fasta.dict
java -jar /opt/picard/CreateSequenceDictionary.jar REFERENCE=canfam3.fasta OUTPUT=canfam3.fasta.dict

#Genes
gunzip download/Canis_familiaris.CanFam3.1.${ENS}.gtf.gz
~/git/Scripts/perl/gtf2bed.pl download/Canis_familiaris.CanFam3.1.${ENS}.gtf > Canis_familiaris.CanFam3.1.${ENS}.bed
gzip download/Canis_familiaris.CanFam3.1.${ENS}.gtf


