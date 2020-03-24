#!/bin/bash

# options for sbatch
#SBATCH --nodes=1
#SBATCH --export=ALL
#SBATCH -o outputs/job-%j.output

# setting up variables
CF3="${HOME}/data/canfam3/"
RESULTS="${HOME}/data/samples"
SAMPLE=$1

OUTPUT_PATH="${RESULTS}/${SAMPLE}"

#module load apps/samtools
module load apps/picard

#samtools flagstat ${OUTPUT_PATH}/${SAMPLE}_final.bam > ${OUTPUT_PATH}/flagstat_final.out
picard CollectWgsMetrics INPUT=${OUTPUT_PATH}/${SAMPLE}\_final.bam OUTPUT=${OUTPUT_PATH}/${SAMPLE}\_WGS_metrics.out REFERENCE_SEQUENCE=${CF3}/ensembl/canfam3.fasta STOP_AFTER=100000000 VALIDATION_STRINGENCY=LENIENT
picard CollectInsertSizeMetrics INPUT=${OUTPUT_PATH}/${SAMPLE}\_final.bam OUTPUT=${OUTPUT_PATH}/${SAMPLE}\_insert_size.out HISTOGRAM_FILE=${OUTPUT_PATH}/${SAMPLE}\_insert_size.pdf VALIDATION_STRINGENCY=LENIENT

echo ${SAMPLE}
echo "..."

#grep 'in total' ${OUTPUT_PATH}/flagstat_final.out
#grep '%' ${OUTPUT_PATH}/flagstat_final.out | grep -v 'singleton'
#echo "..."

head -8 ${OUTPUT_PATH}/${SAMPLE}_WGS_metrics.out | tail -2 | cut -f 2,3,13-18
echo "..."

head -8 ${OUTPUT_PATH}/${SAMPLE}_insert_size.out | tail -2 | cut -f 5,7
echo "..."
