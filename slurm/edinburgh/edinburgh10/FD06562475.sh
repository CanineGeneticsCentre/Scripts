#!/bin/bash

export ASPERA_SCP_PASS=DividePropertyPicture
cd /users/eschofield/data/downloads/edinburgh10/

mkdir FD06562475/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06562475/FBD_33745_R1.fastq.gz FD06562475/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06562475/FBD_33745_R1.fastq.gz.md5 FD06562475/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06562475/FBD_33745_R1_fastqc.zip FD06562475/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06562475/FBD_33745_R1_fastqc.html FD06562475/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06562475/FBD_33745_R2.fastq.gz FD06562475/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06562475/FBD_33745_R2.fastq.gz.md5 FD06562475/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06562475/FBD_33745_R2_fastqc.zip FD06562475/
ascp -P 33001 -O 33001 -l 20M eschofield@transfer.epcc.ed.ac.uk:X18168/2019-01-22/FD06562475/FBD_33745_R2_fastqc.html FD06562475/
