#!/usr/bin/env python
from __future__ import print_function
try: input = raw_input
except: pass
import sys

# For various kmer sizes, find the number of sequencing errors that would add a
# new edges between two existing kmers
# Note: there are 3*reflen possible mutations

cov = None
seqn_err_rate = None

def usage():
  print("usage: python count-bad-edges.py [<cov> <seqn-err-rate>]",file=sys.stderr)
  exit(-1)

if len(sys.argv) != 1 and len(sys.argv) != 3: usage()
if len(sys.argv) == 3:
  try: cov,seqn_err_rate = int(sys.argv[1]),float(sys.argv[2])
  except: usage()

s = input().strip().upper()
kmers = [21,31,41,51,61,71,81,91,99]

def est_num_of_added_edges(reflen,nerror_edges,cov,seqn_err_rate):
  nerrors = cov * reflen * seqn_err_rate
  bad_error_rate = nerror_edges / (3.0 * reflen) # errs that create new edges
  return int(nerrors*bad_error_rate)

print("# The number of sequencing errors that would add a new edge between two")
print("# existing kmers. Note: there are 3*reflen possible mutations")

cols = ["kmer","reflen","nkmers","nedges","nerror_edges"]
if cov is not None: cols.extend(["cov","err_rate","est_bad_edges"])
print(",".join([str(x) for x in cols]))

for k in kmers:
  kmers = set()
  edges = set()
  for i in range(len(s)-k):
    kmers.add(s[i:i+k])
    edges.add(s[i:i+k+1])
  kmers.add(s[-k:])
  err_edges = 0
  pk = s[0:k] # first kmer
  for i in range(1,len(s)-k+1):
    nextb = s[i+k-1]
    for c in "ACGT":
      err_edges += (c != nextb and pk[1:]+c in kmers and pk+c not in edges)
    pk = s[i:i+k]
  cols = [k,len(s),len(kmers),len(edges),err_edges]
  if cov is not None:
    cols.extend([cov, seqn_err_rate,
                 est_num_of_added_edges(len(s),err_edges,cov,seqn_err_rate)])
  print(",".join([str(x) for x in cols]))
