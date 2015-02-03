#!/bin/bash

# xtrace prints commands as we run them
set -euo pipefail

for k in 31 63 95 127
do
  make MAXK=$k $@
done
