#!/bin/bash

set -oeu pipefail

cd libs && make clean && cd ..
make clean
git pull
git submodule update --init
make
