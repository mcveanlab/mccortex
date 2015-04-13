Record time to read a sequence file buffered vs un-buffered

./file-buffering.sh ../data/chr22/chr22.fa > results20150413mon.csv
../hash_table_benchmark/stats.R results20150413mon.csv > results20150413mon.txt
