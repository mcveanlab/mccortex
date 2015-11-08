#!/bin/bash

set -euo pipefail

#
# Run all the tests!
#   Isaac Turner
#   2014-07-16
#
# cd into each directory and run `make`
#

if [[ ( $# -gt 1 ) || ( $# -eq 1 && $1 != 'noupdate' && $1 != 'update' ) ]]
then
  echo "./run [update|noupdate]"
  exit -1
fi

cwd=`pwd`
echo $cwd

if [[ $# -eq 0 || $1 == 'update' ]]
then
  # Get all dependencies used in testing (bioinf-perl, bcftools, samtools etc.)
  cd ../libs && make all && cd $cwd
  if [ $? -ne 0 ]; then exit -1; fi
fi

# Run cortex unit tests
cd ..
for k in 31 63 95 127
do
  make test MAXK=$k STRICT=1
done
cd $cwd

# Get list of current tests (all directories except 'old')
dirs=`ls | grep -v '.*run.sh' | grep -v '^\.' | grep -v old`
echo $dirs

cd .. && make MAXK=31 RELEASE=1 && make MAXK=63 && cd $cwd
if [ $? -ne 0 ]; then exit -1; fi

for f in $dirs
do
  echo && echo ===== && echo "Test: $cwd/$f"
  cd $f && make clean && make all && cd ..
  if [ $? -ne 0 ]; then exit -1; fi
done

echo $dirs
echo All tests completed.
