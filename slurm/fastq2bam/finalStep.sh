#!/bin/bash -l

#SBATCH --export=ALL
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL
#SBACTH --mem=1
#SBATCH -e logs/job-%j.err
#SBATCH -o logs/job-%j.out

SAMPLE=$1

NOW=$(date +"%Y-%m-%d")

echo ${SAMPLE}
echo "..."

grep 'in total' flagstat_final.out
grep '%' flagstat_final.out | grep -v 'singleton'
echo "..."

head -8 ${SAMPLE}_WGS_metrics.out | tail -2 | cut -f 2,3,13-18
echo "..."

head -8 ${SAMPLE}_insert_size.out | tail -2 | cut -f 5,7
echo "..."

ls -logh ${SAMPLE}_final.bam | cut -d' ' -f 4,5,6,7
