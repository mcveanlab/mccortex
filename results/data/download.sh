#!/bin/bash -euo pipefail

for d in chr22 ecoli
do
  cd $d
  make
  cd ..
done
