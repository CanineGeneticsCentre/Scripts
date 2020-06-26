#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=END,FAIL
#SBACTH --mem=16
##SBATCH -e outputs/job-%j.error
#SBATCH -o logs/job-%j.out

module load apps/picard

SAMPLE=$1

export PICARD_JAVA_TMPDIR="${HOME}/data/javatempdir"
export PICARD_JAVA_MEM_MX='16g'

picard CollectInsertSizeMetrics INPUT=${SAMPLE}\_final.bam OUTPUT=${SAMPLE}\_insert_size.out HISTOGRAM_FILE=${SAMPLE}\_insert_size.pdf VALIDATION_STRINGENCY=LENIENT
picard CollectWgsMetrics INPUT=${SAMPLE}\_final.bam OUTPUT=${SAMPLE}\_WGS_metrics.out REFERENCE_SEQUENCE=${DATA}/${GENOME}/ensembl/${GENOME}.fasta STOP_AFTER=100000000 VALIDATION_STRINGENCY=LENIENT
