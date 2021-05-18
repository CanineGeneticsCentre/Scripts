#!/bin/bash

ENS=102
mkdir -p $PWD/CanFam3.1/e${ENS}/download
cd $PWD/CanFam3.1/e${ENS}/download

##DOWNLOAD DATA
for i in `seq 1 38`; do wget ftp://ftp.ensembl.org/pub/release-${ENS}/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.CanFam3.1.dna_sm.chromosome.$i.fa.gz; done

wget ftp://ftp.ensembl.org/pub/release-${ENS}/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.CanFam3.1.dna_sm.chromosome.MT.fa.gz
wget ftp://ftp.ensembl.org/pub/release-${ENS}/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.CanFam3.1.dna_sm.chromosome.X.fa.gz
wget ftp://ftp.ensembl.org/pub/release-${ENS}/fasta/canis_lupus_familiaris/dna/Canis_lupus_familiaris.CanFam3.1.dna_sm.nonchromosomal.fa.gz

wget ftp://ftp.ensembl.org/pub/release-${ENS}/gtf/canis_lupus_familiaris/Canis_lupus_familiaris.CanFam3.1.${ENS}.gtf.gz

wget ftp://ftp.ensembl.org/pub/release-${ENS}/variation/vcf/canis_lupus_familiaris/canis_lupus_familiaris.vcf.gz
wget ftp://ftp.ensembl.org/pub/release-${ENS}/variation/vcf/canis_lupus_familiaris/canis_lupus_familiaris.vcf.gz.tbi
wget ftp://ftp.ensembl.org/pub/release-${ENS}/variation/vcf/canis_lupus_familiaris/canis_lupus_familiaris_structural_variations.vcf.gz
wget ftp://ftp.ensembl.org/pub/release-${ENS}/variation/vcf/canis_lupus_familiaris/canis_lupus_familiaris_structural_variations.vcf.gz.tbi

cd ../

##PARSE DATA
#Sequence
for i in `seq 1 38`; do zcat download/Canis_lupus_familiaris.CanFam3.1.dna_sm.chromosome.$i.fa.gz >> canfam3.fasta; done
zcat download/Canis_lupus_familiaris.CanFam3.1.dna_sm.chromosome.MT.fa.gz >> canfam3.fasta
zcat download/Canis_lupus_familiaris.CanFam3.1.dna_sm.chromosome.X.fa.gz >> canfam3.fasta
zcat download/Canis_lupus_familiaris.CanFam3.1.dna_sm.nonchromosomal.fa.gz >> canfam3.fasta

module load samtools
samtools faidx canfam3.fasta; 


module load bwa
bwa index -a bwtsw canfam3.fasta

#CreateSequenceDictionary REFERENCE=canfam3.fasta OUTPUT=canfam3.fasta.dict
#java -jar /opt/picard/CreateSequenceDictionary.jar REFERENCE=canfam3.fasta OUTPUT=canfam3.fasta.dict
module load picard/2.9.2
picard_latest CreateSequenceDictionary REFERENCE=canfam3.fasta OUTPUT=canfam3.fasta.dict


#Genes
gunzip download/Canis_lupus_familiaris.CanFam3.1.${ENS}.gtf.gz
~/Git/Scripts/perl/gtf2bed.pl download/Canis_lupus_familiaris.CanFam3.1.${ENS}.gtf > Canis_lupus_familiaris.CanFam3.1.${ENS}.bed
gzip download/Canis_lupus_familiaris.CanFam3.1.${ENS}.gtf

ln -s canfam3.fasta canfam3.fa
ln -s canfam3.fasta.dict canfam3.dict

rm -rf download/

