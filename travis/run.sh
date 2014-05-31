#!/bin/bash

set -e

if [ "${TRAVIS_BRANCH}" != 'coverity_scan' ] || \
   ( [ "${CC}" == "gcc" ] && \
     [ "${TRAVIS_OS_NAME}" == "linux" ] );
then
  echo yes
else
  echo no
fi
