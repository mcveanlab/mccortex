#!/bin/sh

if [ `./travis/run.sh` == 'yes' ]
then
  cd libs && make core && cd ..
fi
