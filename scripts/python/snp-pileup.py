#!/usr/bin/env python

from __future__ import print_function

import sys

# Adapted from http://blog.nextgenetics.net/?e=56
# Initially by Damian Kao

# Prints:
#   chr pos ref cov A C G T qA  qC  qG  qT
#   NC_009648.1 21  A 31  0 0 31  0 0 0 38  0
#   NC_009648.1 24  A 32  0 0 32  0 0 0 37  0
#

# Docs:
#   http://samtools.sourceforge.net/samtools.shtml
#   http://samtools.sourceforge.net/pileup.shtml

# samtools mpileup -BQ0 -d10000000 -f /data/ucsc/hg19/hg19.fa AG1.bam | python ../scripts/pileup.py

# inFile = open(sys.argv[1],'r')
inFile = sys.stdin

def num_end(txt,pos):
  j = pos
  while j < len(txt) and ord(txt[j]) >= ord('0') and ord(txt[j]) <= ord('9'):
    j += 1
  return j

def mean_qual(qsum,nq):
  if nq == 0:
    return 0
  else:
    return int((qsum / nq)+0.5);

print('#'+'\t'.join(['chr','pos','ref','cov','A','C','G','T','qA','qC','qG','qT']))

for line in inFile:
  data = line.strip().split('\t')
  chrom = data[0]
  pos = data[1]
  refbase = data[2].upper()
  covg = int(data[3])

  types = {'A':0,'C':0,'G':0,'T':0,'N':0}
  qsum = {'A':0,'C':0,'G':0,'T':0,'N':0}

  # print line

  if(covg > 0):
    bases = data[4].upper()
    quals = data[5]

    i = 0
    q = 0
    while i < len(bases):
      base = bases[i]
      # print("line:"+str(i)+":"+str(q)+": "+base+"; refbase:"+refbase+"; quals[q]:"+quals[q])
      if base == '^':
        i += 2 # ^X  =>  ^ start of read, X is 33+MAPQ
      elif base == '$':
        i += 1 # read end
      elif base == '-' or base == '+':
        i += 1 # indel
        end = num_end(bases, i)
        i = end + int(bases[i:end])
      elif base == '*':
        i += 1 # previous deletion
        q += 1
      elif base == '.' or base == ',':
        types[refbase] += 1
        qsum[refbase] += ord(quals[q])-33;
        q += 1
        i += 1
      elif base in types:
        types[base] += 1
        qsum[base] += ord(quals[q])-33;
        q += 1
        i += 1
      elif base == '<' or base == '>':
        # Reference skip
        i += 1
        q += 1
      else:
        print("q="+str(q)+";i="+str(i),file=sys.stderr)
        print(line,file=sys.stderr)
        print("Something is seriously wrong.",file=sys.stderr)
        sys.exit(-1)

  if q < len(quals):
    raise SystemExit("Didn't finish quality scores: %i < %i: %s" % (q, len(quals), line))

  out = [chrom, pos, refbase, covg,
         types['A'], types['C'], types['G'], types['T'],
         mean_qual(qsum['A'],types['A']),
         mean_qual(qsum['C'],types['C']),
         mean_qual(qsum['G'],types['G']),
         mean_qual(qsum['T'],types['T'])]

  print('\t'.join([str(x) for x in out]))
