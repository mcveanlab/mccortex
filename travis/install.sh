#!/bin/bash

set -e

RUN_TRAVIS=`./travis/run.sh`

if [ "$RUN_TRAVIS" == "yes" ]
then
  # Fetch third party code required to compile
  cd libs && make core && cd ..
fi
