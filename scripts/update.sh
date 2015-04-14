#!/bin/bash

set -oeu pipefail

make clean
git pull
git submodule update --init
make
