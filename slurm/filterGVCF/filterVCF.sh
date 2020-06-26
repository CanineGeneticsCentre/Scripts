#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL
##SBATCH -e outputs/job-%j.error
#SBATCH -o outputs/job-%j.output

module load apps/gatk

SAMPLE=$1

export GATK_JAVA_TMPDIR="${HOME}/data/javatempdir"

gatk -T VariantFiltration -R ${CF3}/ensembl/canfam3.fasta -V ${RESULTS}/${SAMPLE}/${SAMPLE}.vcf.gz --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0 || SOR > 3.0 || QUAL < 30" --filterName "basic_filters" -o ${RESULTS}/${SAMPLE}/${SAMPLE}.filtered.vcf.gz

#rm -rf ${RESULTS}/${SAMPLE}/${SAMPLE}.vcf.gz