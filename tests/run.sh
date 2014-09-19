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
  cd ../libs && make core common && cd $cwd
fi

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
