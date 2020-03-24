#!/bin/bash

export ASPERA_SCP_PASS=DividePropertyPicture
cd ~/data/downloads/edinburgh9/

mkdir FD06562758/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562758/SCWT_34784_R1.fastq.gz FD06562758/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562758/SCWT_34784_R1.fastq.gz.md5 FD06562758/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562758/SCWT_34784_R1_fastqc.zip FD06562758/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562758/SCWT_34784_R1_fastqc.html FD06562758/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562758/SCWT_34784_R2.fastq.gz FD06562758/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562758/SCWT_34784_R2.fastq.gz.md5 FD06562758/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562758/SCWT_34784_R2_fastqc.zip FD06562758/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562758/SCWT_34784_R2_fastqc.html FD06562758/

cat SCWT_34784_R1.fastq.gz.md5; echo;
md5sum SCWT_34784_R1.fastq.gz; echo;
cat SCWT_34784_R2.fastq.gz.md5; echo;
md5sum SCWT_34784_R2.fastq.gz
