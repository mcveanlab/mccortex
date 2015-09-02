#!/bin/bash

set -e

RUN_TRAVIS=`./travis/run.sh`

if [ "$RUN_TRAVIS" == "yes" ]
then

  # Install JSON perl package
  # cpanm JSON    <- Only available when language is `perl`
  if [ "$TRAVIS_OS_NAME" == "osx" ]
  then
    # brew update
    # brew outdated <package-name> || brew upgrade <package-name>
    # cpanm JSON
    # sudo cpan JSON
    sudo cpanm JSON
  elif [ "$TRAVIS_OS_NAME" == "linux" ]
  then
    sudo apt-get install libjson-pp-perl
  else
    echo "Unexpected TRAVIS_OX_NAME: $TRAVIS_OS_NAME"
    exit -1
  fi

  # Fetch third party code required to compile
  cd libs && make && cd ..
fi
