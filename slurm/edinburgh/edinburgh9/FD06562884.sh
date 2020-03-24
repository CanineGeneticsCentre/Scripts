#!/bin/bash

export ASPERA_SCP_PASS=DividePropertyPicture
cd ~/data/downloads/edinburgh9/

mkdir FD06562884/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562884/WI_34774_R1.fastq.gz FD06562884/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562884/WI_34774_R1.fastq.gz.md5 FD06562884/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562884/WI_34774_R1_fastqc.zip FD06562884/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562884/WI_34774_R1_fastqc.html FD06562884/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562884/WI_34774_R2.fastq.gz FD06562884/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562884/WI_34774_R2.fastq.gz.md5 FD06562884/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562884/WI_34774_R2_fastqc.zip FD06562884/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18066/2018-06-18/FD06562884/WI_34774_R2_fastqc.html FD06562884/

cat WI_34774_R1.fastq.gz.md5;echo;md5sum WI_34774_R1.fastq.gz;echo;cat WI_34774_R2.fastq.gz.md5;echo;md5sum WI_34774_R2.fastq.gz
