#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL
##SBATCH -e outputs/job-%j.error
#SBATCH -o outputs/job-%j.output
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32000

module load apps/gatk

RUN_NAME=$1
BATCH=$2

export GATK_JAVA_TMPDIR="${HOME}/data/javatempdir"
export GATK_JAVA_MEM_MX='32g'

value=`cat ${RESULTS}/${BATCH}.list`
variants=""
while read SAMPLE; do
	var=`echo -n " --variant ${RESULTS}/$SAMPLE/$SAMPLE.g.vcf.gz"`
    variants="$variants $var"
done < ${RESULTS}/${BATCH}.list

gatk -T CombineGVCFs -R ${CF3}/ensembl/canfam3.fasta -o ${RESULTS}/${BATCH}.g.vcf.gz $variants
