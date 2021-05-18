#!/bin/bash

ENS=102
mkdir ensembl_v${ENS}
cd ensembl_v${ENS}

git clone https://github.com/Ensembl/ensembl-git-tools.git
export PATH=$PWD/ensembl-git-tools/bin:$PATH

git ensembl --clone api

cd ../
ln -sf ensembl_v${ENS} ENSEMBL
