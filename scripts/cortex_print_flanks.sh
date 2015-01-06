#/bin/bash
set -euo pipefail

gzip -fcd $@ | awk -F '[ \t]' 'm{print $0;m=0;} /^>bubble\..*\.5pflank/{print $1; m=1;}'
