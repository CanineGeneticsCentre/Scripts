#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=END,FAIL
#SBACTH --mem=32G
#SBATCH --ntasks=16
#SBATCH -o outputs/job-%j.output
#SBATCH --nodes=1


module load apps/gatk

RUN_NAME=$1
SAMPLE=$2

export GATK_JAVA_TMPDIR="${HOME}/data/javatempdir"
export GATK_JAVA_MEM_MX='32g'

OUTPUT_PATH="${RESULTS}/${SAMPLE}"

rm -rf ${OUTPUT_PATH}/${SAMPLE}.g.vcf*
time gatk -T HaplotypeCaller -R ${CF3}/ensembl/canfam3.fasta -I ${OUTPUT_PATH}/${SAMPLE}\_final.bam -o ${OUTPUT_PATH}/${SAMPLE}.g.vcf.gz --emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -S LENIENT -nct 16
