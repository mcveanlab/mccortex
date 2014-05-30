#!/bin/bash
set -eu
set -o pipefail
IFS=$'\n\t'

gzip -dc ERR039560.breakpoints.txt.gz | ../../../scripts/make-circos.pl my-plot -
cd my-plot
circos
