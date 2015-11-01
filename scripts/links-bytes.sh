#!/bin/bash

set -eou pipefail

zcat -fcd $1 | grep '^[FR] ' | awk '{x=x+int(($2+3)/4)}END{print x}'
