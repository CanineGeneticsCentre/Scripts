#!/bin/bash

export ASPERA_SCP_PASS=DividePropertyPicture
cd ~/data/downloads/edinburgh9/

mkdir FD06561391/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06561391/WHWT_34443_R1.fastq.gz FD06561391/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06561391/WHWT_34443_R1.fastq.gz.md5 FD06561391/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06561391/WHWT_34443_R1_fastqc.zip FD06561391/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06561391/WHWT_34443_R1_fastqc.html FD06561391/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06561391/WHWT_34443_R2.fastq.gz FD06561391/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06561391/WHWT_34443_R2.fastq.gz.md5 FD06561391/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06561391/WHWT_34443_R2_fastqc.zip FD06561391/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06561391/WHWT_34443_R2_fastqc.html FD06561391/

cat WHWT_34443_R1.fastq.gz.md5;echo;md5sum WHWT_34443_R1.fastq.gz;echo;cat WHWT_34443_R2.fastq.gz.md5;echo;md5sum WHWT_34443_R2.fastq.gz
