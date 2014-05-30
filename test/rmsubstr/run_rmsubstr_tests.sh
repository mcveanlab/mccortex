#!/bin/bash

# Exit on error (and if diff returns non zero)
set -ue
set -o pipefail
IFS=$'\n\t'

for i in {0..4}
do
  ../../bin/ctx31 rmsubstr -n 1024 -k 3 test.$i.fa > out.$i.fa 2> /dev/null
  diff -q results.$i.fa out.$i.fa
done

echo "All rmsubstr tests passed."
