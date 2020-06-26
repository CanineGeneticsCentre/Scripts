#!/usr/bin/bash

export ASPERA_SCP_PASS=DividePropertyPicture
cd /users/eschofield/data/downloads/edinburgh12/

mkdir FD06234274/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/FD06234274/CS_10842_R1.fastq.gz FD06234274/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/FD06234274/CS_10842_R1.fastq.gz.md5 FD06234274/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/FD06234274/CS_10842_R1_fastqc.zip FD06234274/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/FD06234274/CS_10842_R1_fastqc.html FD06234274/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/FD06234274/CS_10842_R2.fastq.gz FD06234274/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/FD06234274/CS_10842_R2.fastq.gz.md5 FD06234274/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/FD06234274/CS_10842_R2_fastqc.zip FD06234274/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X19070/2019-07-29/FD06234274/CS_10842_R2_fastqc.html FD06234274/
