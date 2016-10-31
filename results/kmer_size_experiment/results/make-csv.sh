#!/bin/bash
set -eou pipefail

echo "K,NG50,AssemblyErrors"
for f in $@
do
  K=`echo $f | grep -oE 'k[0-9]+' | grep -oE '[0-9]+$'`
  NG50=`grep 'NG50:' $f | grep -oE '[0-9]+$'`
  ERRORS=`grep 'assembly_errors:' $f | grep -oE '[0-9]+$'`
  echo "$K,$NG50,$ERRORS"
done
