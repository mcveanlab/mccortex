#!/bin/sh

set -e

if [ "`./travis/run.sh`" == 'yes' ]
then
  if [ "${TRAVIS_BRANCH}" == 'coverity_scan' ]
  then
    ./travis/travisci_build_coverity_scan.sh
  else
    make all RELEASE=1 && make clean && make all && make test
  fi
fi
