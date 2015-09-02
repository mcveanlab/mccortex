#!/bin/bash

set -e

RUN_TRAVIS=`./travis/run.sh`

if [ "$RUN_TRAVIS" == "yes" ]
then

  # Set up installing perl modules locally using local::lib
  # eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"

  # Set up cpanm, install JSON perl package
  # Using default ~/perl5 local directory
  curl -L https://cpanmin.us | perl - App::cpanminus
  # ~/perl5/bin/cpanm local::lib
  ~/perl5/bin/cpanm JSON

  # Fetch third party code required to compile
  cd libs && make && cd ..

fi
