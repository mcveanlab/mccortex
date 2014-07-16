#!/bin/bash

set -eu
set -o pipefail

dirs=`ls | grep -v '.*run.sh' | grep -v old`

# remove clean2 for now
dirs=`echo "$dirs" | grep -v clean2`

echo $dirs

cwd=`pwd`
cd .. && make MAXK=31 && make MAXK=63 && cd $cwd

for f in $dirs
do
  echo =====
  echo $f
  cd $f; make clean; make all; cd ..
done

echo $dirs
echo All tests completed.
