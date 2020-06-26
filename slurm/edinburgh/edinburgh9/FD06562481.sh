#!/bin/bash

export ASPERA_SCP_PASS=DividePropertyPicture
cd data/downloads/edinburgh9/

mkdir FD06562481/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562481/MWHD_29814_R1.fastq.gz FD06562481/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562481/MWHD_29814_R1.fastq.gz.md5 FD06562481/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562481/MWHD_29814_R1_fastqc.zip FD06562481/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562481/MWHD_29814_R1_fastqc.html FD06562481/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562481/MWHD_29814_R2.fastq.gz FD06562481/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562481/MWHD_29814_R2.fastq.gz.md5 FD06562481/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562481/MWHD_29814_R2_fastqc.zip FD06562481/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562481/MWHD_29814_R2_fastqc.html FD06562481/
