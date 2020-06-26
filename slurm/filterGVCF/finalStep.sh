#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL
#SBACTH --mem=1
##SBATCH -e outputs/job-%j.error
#SBATCH -o outputs/job-%j.output

SAMPLE=$1

echo
echo ${SAMPLE}
echo "..."

ll ${OUTPUT_PATH}/${SAMPLE}.vcf.gz
zcat ${OUTPUT_PATH}/${SAMPLE}.vcf.gz | grep -v '#' | wc -l
echo "..."

ll ${OUTPUT_PATH}/${SAMPLE}.filtered.vcf.gz
zcat ${OUTPUT_PATH}/${SAMPLE}.filtered.vcf.gz | grep -v '#' | grep -v 'basic_filters' | wc -l
echo