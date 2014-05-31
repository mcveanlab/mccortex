#!/bin/sh

set -e

RUN_TRAVIS=`./travis/run.sh`

if [ "$RUN_TRAVIS" == "yes" ]
then
  cd libs && make core && cd ..
fi
