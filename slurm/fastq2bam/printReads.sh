#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=END,FAIL
#SBACTH --mem=32G
#SBATCH --ntasks=16
#SBATCH -o logs/job-%j.out

module load apps/gatk

SAMPLE=$1

export GATK_JAVA_TMPDIR="${HOME}/data/javatempdir"
export GATK_JAVA_MEM_MX='32g'

# time gatk -T PrintReads -R ${DATA}/${GENOME}/ensembl/${GENOME}.fasta -I ${OUTPUT_PATH}/${SAMPLE}_clean.dedup.bam -nct 16 -o ${OUTPUT_PATH}/${SAMPLE}_final.bam -BQSR ${OUTPUT_PATH}/${SAMPLE}.recal.grp -S LENIENT
time gatk -T PrintReads -R ${DATA}/${GENOME}/ensembl/${GENOME}.fasta -I ${SAMPLE}_clean.dedup.bam -nct 16 -o ${SAMPLE}_final.bam -BQSR ${SAMPLE}.recal.grp -S LENIENT --disable_indel_quals

bam_size=$(wc -c < ${SAMPLE}_final.bam)
if [ $bam_size -ge 50000000 ]; then
	rm -rf ${SAMPLE}_clean.dedup.bam
	rm -rf ${SAMPLE}_clean.dedup.bai
	#rm -rf ${OUTPUT_PATH}/${SAMPLE}.recal.grp
fi
