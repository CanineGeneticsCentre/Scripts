#!/bin/bash

export ASPERA_SCP_PASS=DividePropertyPicture
cd /users/eschofield/data/downloads/edinburgh10/

mkdir FD06561210/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06561210/ISP_28899_R1.fastq.gz FD06561210/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06561210/ISP_28899_R1.fastq.gz.md5 FD06561210/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06561210/ISP_28899_R1_fastqc.zip FD06561210/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06561210/ISP_28899_R1_fastqc.html FD06561210/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06561210/ISP_28899_R2.fastq.gz FD06561210/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06561210/ISP_28899_R2.fastq.gz.md5 FD06561210/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06561210/ISP_28899_R2_fastqc.zip FD06561210/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06561210/ISP_28899_R2_fastqc.html FD06561210/
