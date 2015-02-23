#!/bin/bash

set -euo pipefail

for k in 31 63 95 127
do
  make MAXK=$k $@
done
