#!/bin/bash

set -e

echo "Branch: ${TRAVIS_BRANCH}"
echo "OS: ${TRAVIS_OS_NAME}"
echo "CC: ${CC}"
echo "Perl: ${TRAVIS_PERL_VERSION}"

RUN_TRAVIS=`./travis/run.sh`

if [ "$RUN_TRAVIS" == "yes" ]
then
  if [ "${TRAVIS_BRANCH}" == 'coverity_scan' ]
  then
    ./travis/travisci_build_coverity_scan.sh
  else
    make all RELEASE=1 && make clean
    cd tests && ./run.sh
  fi
fi
