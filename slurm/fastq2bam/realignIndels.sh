#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=END,FAIL
#SBACTH --mem=32G
#SBATCH --ntasks=1
#SBATCH -o logs/job-%j.out

module load apps/gatk

SAMPLE=$1

export GATK_JAVA_TMPDIR="${HOME}/data/javatempdir"
export GATK_JAVA_MEM_MX='32g'

time gatk -T IndelRealigner -R ${DATA}/${GENOME}/ensembl/${GENOME}.fasta -I ${SAMPLE}_dedup.bam -targetIntervals ${SAMPLE}_dedup.bam.intervals -o ${SAMPLE}_clean.dedup.bam -S LENIENT --filter_mismatching_base_and_quals

bam_size=$(wc -c < ${SAMPLE}_clean.dedup.bam)
if [ $bam_size -ge 50000000 ]; then
	rm -rf ${SAMPLE}_dedup.bam
	rm -rf ${SAMPLE}_dedup.bai
	rm -rf ${SAMPLE}_dedup.bam.intervals
fi
