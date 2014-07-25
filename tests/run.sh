#!/bin/bash

set -euo pipefail

#
# Run all the tests!
#   Isaac Turner
#   2014-07-16
#
# cd into each directory and run `make`
#

cwd=`pwd`
echo $cwd

# Get all dependencies used in testing (bioinf-perl, bcftools, samtools etc.)
cd ../libs && make update && cd $cwd

dirs=`ls | grep -v '.*run.sh' | grep -v '^\.' | grep -v old`
echo $dirs

cd .. && make MAXK=31 && make MAXK=63 && cd $cwd

for f in $dirs
do
  echo && echo ===== && echo "Test: $cwd/$f";
  cd $f && make clean && make all && cd ..
done

echo $dirs
echo All tests completed.
