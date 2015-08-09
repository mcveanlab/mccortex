#!/bin/bash

set -e

echo "Branch: ${TRAVIS_BRANCH}"
echo "OS: ${TRAVIS_OS_NAME}"
echo "CC: ${CC}"
echo "Perl: ${TRAVIS_PERL_VERSION}"

# Check we can build a release build
make all RELEASE=1 && make test && make clean

# Build and run all tests
cd tests && ./run.sh
