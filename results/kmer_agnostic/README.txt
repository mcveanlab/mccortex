Isaac Turner
2014-10-01

Test effect of using different kmer sizes with perfect coverage

  # In the cortex main directory
  make MAXK=31
  make MAXK=63
  make MAXK=95
  make MAXK=127

  # generate reference
  cd results/chr22/uniq_flanks && make && cd ../../..

  # run experiment
  cd results/kmer_agnostic

  # To produce all the files
  make
  # To run the checks
  make checks

Run time is ~10 mins on my macbook pro. RAM is about 200MB.
