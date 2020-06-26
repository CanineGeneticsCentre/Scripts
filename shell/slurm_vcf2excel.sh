#!/bin/bash

## Command to Run :  sbatch slurm_vcf2excel.sh <Disease Status File> <Seg Score>

DIS_STATUS=$1
SEG_SCORE=$2

VCF_FILE=${HOME}/data/ens_WGS_186_VEP.vcf.gz
EFFECT_SCORE=1

perl ${HOME}/data/git/Scripts/perl/vcf2excel.pl --default --vcf $VCF_FILE --status_file $DIS_STATUS --effect_score $EFFECT_SCORE --seg_score $SEG_SCORE

