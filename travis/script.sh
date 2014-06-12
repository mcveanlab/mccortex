#!/bin/bash

set -e

echo "Branch: ${TRAVIS_BRANCH}"
echo "OS: ${TRAVIS_OS_NAME}"
echo "CC: ${CC}"

RUN_TRAVIS=`./travis/run.sh`

if [ "$RUN_TRAVIS" == "yes" ]
then
  if [ "${TRAVIS_BRANCH}" == 'coverity_scan' ]
  then
    ./travis/travisci_build_coverity_scan.sh
  else
    make all RELEASE=1 && make clean && make all && make test
    make test MAXK=63
    make test MAXK=95
  fi
fi
