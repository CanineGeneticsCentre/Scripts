#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL
##SBATCH -e outputs/job-%j.error
#SBATCH -o outputs/job-%j.output

module load apps/gatk

export GATK_JAVA_TMPDIR="${HOME}/data/javatempdir"
export GATK_JAVA_MEM_MX='32g'

variants=""
while read BATCH_FILE; do
	var=`echo -n " --variant ${RESULTS}/gVCF_batches/$BATCH_FILE"`
    variants="$variants $var"
done < ${RESULTS}/gVCF_batches/gvcf2vcf.list

gatk -T GenotypeGVCFs -R ${CF3}/ensembl/canfam3.fasta -o ${RESULTS}/merged.vcf.gz $variants