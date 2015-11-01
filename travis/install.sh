#!/bin/bash

set -e

# The COVERITY_SCAN_BRANCH environment variable will be set to 1 when the
# Coverity Scan addon is in operation
# Only run if we are not doing Coverity Scan analysis
if [ "${COVERITY_SCAN_BRANCH}" != 1 ]
then

  # Set up cpanm, install JSON perl package
  # Using default ~/perl5 local directory
  curl -L https://cpanmin.us | perl - App::cpanminus
  ~/perl5/bin/cpanm --force --local-lib=~/perl5 local::lib && eval $(perl -I ~/perl5/lib/perl5/ -Mlocal::lib)
  ~/perl5/bin/cpanm --force --local-lib=~/perl5 JSON

  # Set up installing perl modules library path with:
  # eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"
  echo '[ $SHLVL -eq 1 ] && eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"' >> ~/.bashrc
  echo '[ $SHLVL -eq 1 ] && eval "$(perl -I$HOME/perl5/lib/perl5 -Mlocal::lib)"' >> ~/.profile

  # Compile third party code
  cd libs && make && cd ..

fi
