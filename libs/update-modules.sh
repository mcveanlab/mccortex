#!/bin/bash

set -oeu pipefail
set -o xtrace

# Could also just run:
#   git submodule update --remote
# But that wouldn't run the tests...

# cd xxHash
# git checkout master
# git pull
# make
# cd ..

cd htslib
git checkout develop
git pull
make all
make test >& /dev/null
cd ..

cd bcftools
git checkout develop
git pull
make all
# make test >& /dev/null
cd ..

cd samtools
git checkout develop
git pull
make all
make test >& /dev/null
cd ..

cd bwa
git checkout master
git pull
make
cd ..

cd bfc
git checkout master
git pull
make
cd ..

cd string_buffer
git checkout master
git pull
make all
make test >& /dev/null
cd ..

cd bit_array
git checkout master
git pull
make all
make test >& /dev/null
cd ..

cd seq_file
git checkout master
git pull
make HTSLIB=../htslib all
make test >& /dev/null
cd ..

cd seq-align
git checkout master
git pull
make LIBS_PATH=.. all
make LIBS_PATH=.. test >& /dev/null
cd ..

cd readsim
git checkout master
git pull
make HTSLIB=../htslib
cd ..

cd msg-pool
git checkout master
git pull
make all
make test >& /dev/null
cd ..

cd bioinf-perl
git checkout master
git pull
cd ..

cd sort_r
git checkout master
git pull
make all
make test >& /dev/null
cd ..

cd madcrowlib
git checkout master
git pull
make all
make test >& /dev/null
cd ..

cd biogrok
git checkout master
git pull
cd ..

cd vcf-slim
git checkout master
git pull
make all
make test >& /dev/null
cd ..

cd carrays
git checkout master
git pull
make all
make test >& /dev/null
cd ..

# Now add and commit the updated submodules

# Not external repos
#cd maximal_substrs
# make
#cd misc
# make
