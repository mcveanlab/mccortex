#!/bin/bash

set -e

RUN_TRAVIS=`./travis/run.sh`

if [ "$RUN_TRAVIS" == "yes" ]
then

  # Set up cpanm, install JSON perl package
  # Using default ~/perl5 local directory
  curl -L https://cpanmin.us | perl - App::cpanminus
  ~/perl5/bin/cpanm local::lib
  ~/perl5/bin/cpanm JSON

  # Set up installing perl modules library path with:
  # eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"
  echo '[ $SHLVL -eq 1 ] && eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"' >> ~/.bashrc
  echo '[ $SHLVL -eq 1 ] && eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"' >> ~/.profile

  # Compile third party code
  cd libs && make && cd ..

fi
