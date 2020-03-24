#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=END,FAIL
#SBATCH -o logs/job-%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=12
##SBATCH --time=09:00:00

module load apps/bwa

SAMPLE=$1

echo "bwa mem -M -t 12 ${DATA}/${GENOME}/ensembl/${GENOME}.fasta ${SAMPLES}/${SAMPLE}/${SAMPLE}\_R1.fastq.gz ${SAMPLES}/${SAMPLE}/${SAMPLE}\_R2.fastq.gz > ${SAMPLE}\_aligned.sam"
time bwa mem -M -t 12 ${DATA}/${GENOME}/ensembl/${GENOME}.fasta ${SAMPLES}/${SAMPLE}/${SAMPLE}\_R1.fastq.gz ${SAMPLES}/${SAMPLE}/${SAMPLE}\_R2.fastq.gz > ${SAMPLE}\_aligned.sam
