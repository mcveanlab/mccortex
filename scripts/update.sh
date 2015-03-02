#!/bin/bash

set -oeu pipefail

make clean && git pull && cd libs && make && cd .. && make
