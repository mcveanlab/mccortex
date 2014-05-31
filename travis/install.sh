#!/bin/sh

set -e

if [ "`./travis/run.sh`" == 'yes' ]
then
  cd libs && make core && cd ..
fi
