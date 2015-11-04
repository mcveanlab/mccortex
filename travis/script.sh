#!/bin/bash

set -e

echo "Branch: ${TRAVIS_BRANCH}"
echo "OS: ${TRAVIS_OS_NAME}"
echo "CC: ${CC}"
echo "Perl: ${TRAVIS_PERL_VERSION}"

# The COVERITY_SCAN_BRANCH environment variable will be set to 1 when the
# Coverity Scan addon is in operation
# Only run if we are not doing Coverity Scan analysis
if [ "${COVERITY_SCAN_BRANCH}" != 1 ]
then
  # Build and run all tests
  cd tests && ./run.sh
fi
