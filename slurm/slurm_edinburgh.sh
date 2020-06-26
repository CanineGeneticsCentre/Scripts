#!/bin/bash

## Command to Run :  sbatch ~/scripts/slurm_edinburgh.sh <SAMPLE [FD06561107]>

#SBATCH --export=ALL
#SBATCH -o outputs/job-%j.output
#SBATCH --ntasks=1
##SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --mail-user=ellen.schofield@aht.org.uk
#SBATCH --mail-type=ALL

SAMPLE=$1

export ASPERA_SCP_PASS=DividePropertyPicture
cd ~/data/downloads/edinburgh12/

AHT_NO=`grep ${SAMPLE} samples.txt | cut -f 1`

echo $AHT_NO;
mkdir $SAMPLE
ln -s $SAMPLE $AHT_NO

ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/${SAMPLE}/${AHT_NO}_R1.fastq.gz ${SAMPLE}/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/${SAMPLE}/${AHT_NO}_R1.fastq.gz.md5 ${SAMPLE}/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/${SAMPLE}/${AHT_NO}_R1_fastqc.zip ${SAMPLE}/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/${SAMPLE}/${AHT_NO}_R1_fastqc.html ${SAMPLE}/

ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/${SAMPLE}/${AHT_NO}_R2_fastq.gz ${SAMPLE}/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/${SAMPLE}/${AHT_NO}_R2_fastq.gz.md5 ${SAMPLE}/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/${SAMPLE}/${AHT_NO}_R2_fastqc.html ${SAMPLE}/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/${SAMPLE}/${AHT_NO}_R2_fastqc.zip ${SAMPLE}/
