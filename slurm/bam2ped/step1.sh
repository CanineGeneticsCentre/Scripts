#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL
##SBATCH -e outputs/job-%j.error
#SBATCH -o outputs/job-%j.output

module load apps/gatk

RUN_NAME=$1
SAMPLE=$2

export GATK_JAVA_TMPDIR="${HOME}/data/javatempdir"
export GATK_JAVA_MEM_MX='32g'

OUTPUT_PATH="${RESULTS}/${SAMPLE}/${RUN_NAME}"
mkdir -p "$OUTPUT_PATH"

gatk -T HaplotypeCaller -R ${CF3}/ensembl/canfam3.fasta -L ${CF3}/Canine170K.ens.vcf -I ${RESULTS}/${SAMPLE}/${SAMPLE}\_final.bam --emitRefConfidence GVCF -o ${OUTPUT_PATH}/output.g.vcf
mv ${OUTPUT_PATH}/output.g.vcf ${OUTPUT_PATH}/${SAMPLE}.g.vcf; mv ${OUTPUT_PATH}/output.g.vcf.idx ${OUTPUT_PATH}/${SAMPLE}.g.vcf.idx