Isaac Turner
2014-09-25

Using 1Mb of chr22 and empirical PhiX Illumina reads to simulate and measure
the power and error rate of read correction using a de Bruijn graph

Requires mccortex/results/data directory

To reproduce:

  cd libs && make core common && cd ..
  make MAXK=31
  cd results/correct
  make

Runtime on my macbook is ~10 minutes
