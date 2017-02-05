#!/usr/bin/env python
from __future__ import print_function
try: input = raw_input
except: pass
import sys
import csv
import re

#
# Read in distance matrices between raw, clean and perfect graphs for k=...
# output table: kmer,nkmers,raw_nkmers,raw_nreal,clean_nkmers,clean_nreal
#

def usage(argv,err=None):
  if err is not None: print(err,file=sys.stderr)
  print("usage: python",argv[0],"<dist.txt> ...",file=sys.stderr)
  exit(-1)

def load_csv(csvpath):
  m = []
  with open(csvpath) as csvpath:
    csvreader = csv.reader(csvpath, delimiter='\t', quotechar='"')
    next(csvreader) # skip first row (column headers)
    for row in csvreader:
      m.append([ 0 if x == '.' else int(x) for x in row[1:]])
  return m

def main(argv):
  if len(argv) <= 1: usage(argv)
  sep = ','
  print("# Number of kmers in the perfect, raw and cleaned graphs")
  print("# _nreal is the number of real kmers in the raw/cleaned graph")
  print("# raw_errs, clean_errs are the fraction of error kmers in each graph")
  print("# frac_remove_errs is the fraction of kmers removed that were seqn errs")
  print(sep.join(["kmer","nkmers",
                  "raw_nkmers","raw_nreal",
                  "clean_nkmers","clean_nreal",
                  "raw_errs","clean_errs",
                  "frac_remove_errs"]))
  for f in argv[1:]:
    match = re.search('k([0-9]+)', f)
    k = match.group(1)
    m = load_csv(f)
    nkmers = m[2][2]
    raw_nkmers,raw_nreal = m[0][0],m[0][2]
    clean_nkmers,clean_nreal = m[1][1],m[1][2]
    raw_errs = (raw_nkmers-raw_nreal)/float(raw_nkmers)
    clean_errs = (clean_nkmers-clean_nreal)/float(clean_nkmers)
    kmers_removed = raw_nkmers-clean_nkmers
    real_kmers_removed = raw_nreal-clean_nreal
    frac_remove_errs = 1.0 - float(real_kmers_removed)/kmers_removed

    r = [k,m[2][2],m[0][0],m[0][2],m[1][1],m[1][2],
         "%.5f"%raw_errs,"%.5f"%clean_errs,"%.5f"%frac_remove_errs]
    print(sep.join([str(x) for x in r]))

if __name__ == '__main__':
  main(sys.argv)
