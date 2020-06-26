#!/usr/bin/bash

## Command to Run :  bash ~/scripts/edinburgh/edinburgh12.sh <SAMPLE [FD06561107]>

SAMPLE=$1

export ASPERA_SCP_PASS=DividePropertyPicture
cd /users/eschofield/data/downloads/edinburgh12/

if [ -f "samples.txt" ]; then
	AHT_NO=`grep ${SAMPLE} samples.txt | cut -f 1`
	DIR="X19070/2019-07-29"
	
	echo $AHT_NO;
	mkdir $SAMPLE
	ln -s $SAMPLE $AHT_NO

	ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:${DIR}/${SAMPLE}/${AHT_NO}_R1.fastq.gz ${SAMPLE}/
	ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:${DIR}/${SAMPLE}/${AHT_NO}_R1.fastq.gz.md5 ${SAMPLE}/
	ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:${DIR}/${SAMPLE}/${AHT_NO}_R1_fastqc.zip ${SAMPLE}/
	ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:${DIR}/${SAMPLE}/${AHT_NO}_R1_fastqc.html ${SAMPLE}/
	
	ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:${DIR}/${SAMPLE}/${AHT_NO}_R2.fastq.gz ${SAMPLE}/
	ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:${DIR}/${SAMPLE}/${AHT_NO}_R2.fastq.gz.md5 ${SAMPLE}/
	ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:${DIR}/${SAMPLE}/${AHT_NO}_R2_fastqc.html ${SAMPLE}/
	ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:${DIR}/${SAMPLE}/${AHT_NO}_R2_fastqc.zip ${SAMPLE}/
else 
    echo "samples.txt does not exist, please create this file and then re-run this script"
fi

