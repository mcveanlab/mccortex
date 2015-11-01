#!/bin/bash

set -eou pipefail

zcat -fcd $1 | grep -c '^[FR] '
