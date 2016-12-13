#!/usr/bin/env python

from __future__ import print_function
try: input = raw_input
except: pass
import sys, getopt
import random

# TODO: -f,-F,-o,-P
def usage(argv,err=None):
  if err is not None: print(err,file=sys.stderr)
  print("usage:",argv[0],"<frag-size>","<read-length>",file=sys.stderr)
  print("  -p,--paired <fraglen>   Fragment length PE reads (0 means SE only) [0]",file=sys.stderr)
  print("  -r,--readlen <readlen>  Read length",file=sys.stderr)
  print("  -e,--err <errrate>      Per-base error rate 0<=X<1 [0]",file=sys.stderr)
  print("  -d,--depth <depth>      Per-base depth, if zero perfect coverage [default]",file=sys.stderr)
  print("  -s,--seed <rand>        Seed for random number generator",file=sys.stderr)
  # print("  -f,--fvar <var>         Variance on insert size",file=sys.stderr)
  # print("  -F,--ffrac <frac>       Variance on insert size as a fraction",file=sys.stderr)
  # print("  -o,--out <o.fa>[:<o2.fa>] Write to output file(s)",file=sys.stderr)
  # print("  -P,--profile <in.fq>      Match reads, copy error profile",file=sys.stderr)
  exit(-1)

def reverse_complement(s):
  complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', \
                'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
  # get returns the base if not in dict
  bases = reversed([complement.get(b,b) for b in s])
  return ''.join(bases)

mutate = {'a': "cgt", 'c': "agt", 'g': "act", 't': "acg", \
          'A': "CGT", 'C': "AGT", 'G': "ACT", 'T': "ACG"}

def print_read(r,idx,err=None):
  hdr = ">r"+str(idx)+"\n"
  if err is None: print(hdr+r)
  elif isinstance(err, list):
    e = random.choice(err) # match with a random read
    l = [ random.choice(mutate[c]) if random.random() < e[min(i,len(e)-1)] else c for i,c in enumerate(r) ]
    print(hdr+''.join(l))
  else:
    l = [ random.choice(mutate[c]) if random.random() < err else c for c in r ]
    print(hdr+''.join(l))

def perfect_se_cov(s,readlen,err=None):
  for i in range(len(s)-readlen+1):
    print_read(s[i:i+readlen],i+1,err)

def perfect_pe_cov(s,readlen,fraglen,err=None):
  r = reverse_complement(s)
  for i in range(len(s)-fraglen+1):
    print_read(s[i:i+readlen],str(i+1)+"/1",err)
    j = i+fraglen-readlen
    print_read(r[len(s)-j-readlen:len(s)-j],str(i+1)+"/2",err)

def poiss_se_cov(s,readlen,depth,err=None):
  nreads = int((len(s)*depth)/readlen)
  idx = 1
  for _ in range(nreads):
    i = random.randint(0,len(s)-readlen)
    print_read(s[i:i+readlen],idx,err)
    idx += 1

def poiss_pe_cov(s,readlen,fraglen,depth,err=None):
  npairs = int((len(s)*depth)/(readlen*2))
  r = reverse_complement(s)
  idx = 1
  for _ in range(npairs):
    i = random.randint(0,len(s)-fraglen)
    print_read(s[i:i+readlen],str(idx)+"/1",err)
    j = i+fraglen-readlen
    print_read(r[len(s)-j-readlen:len(s)-j],str(idx)+"/2",err)
    idx += 1

# returns a matrix m[r][i] = error prob of the i-th position of the r-th read
def load_err_profile(file):
  l = []
  # TODO: read form file, covert to prob, add to l
  return l

# TODO:
#  -P,--pvar <fraglen-var>
def main(argv):
  readlen = None
  fraglen = None
  err = None
  seed = random.randrange((1<<32)-1)
  depth = 0 # <= 0 -> perfect coverage
  try:
    opts, args = getopt.getopt(argv[1:],"hp:r:d:e:s:",
                               ["help","paired=","readlen=","depth=","err=","seed="])
  except getopt.GetoptError: usage(argv)
  for opt, arg in opts:
    if opt == '-h': usage(argv)
    elif opt in ("-p", "--paired"):
      try: fraglen = int(arg)
      except: usage(argv,"Bad fragment length argument: "+arg)
    elif opt in ("-r", "--readlen"):
      try: readlen = int(arg)
      except: usage(argv,"Bad read length argument: "+arg)
    elif opt in ("-d", "--depth"):
      try: depth = float(arg)
      except: usage(argv,"Bad depth argument: "+arg)
    elif opt in ("-e", "--err"):
      try: err = float(arg)
      except: usage(argv,"Bad error argument: "+arg)
      if err < 0 or err > 1: usage(argv,"Bad error argument: "+arg)
    elif opt in ("-s", "--seed"):
      try: seed = int(arg); random.seed(seed)
      except: usage(argv,"Bad seed argument: "+arg)
    else:
      usage("Bad option:",opt,arg)
  if len(args) > 0: usage(argv,"Unused args:"+str(args))
  if readlen is None: usage(argv,"Require -r,--readlen")
  print("[generate-reads.py] readlen:",readlen,"fraglen:",fraglen,
        "err:",err,"depth:",depth,"seed:",seed,file=sys.stderr)
  # read input
  s = input().strip()
  paired = (fraglen is not None)
  perfectcov = (depth <= 0)
  if     perfectcov and not paired: perfect_se_cov(s,readlen,err)
  if     perfectcov and     paired: perfect_pe_cov(s,readlen,fraglen,err)
  if not perfectcov and not paired: poiss_se_cov(s,readlen,depth,err)
  if not perfectcov and     paired: poiss_pe_cov(s,readlen,fraglen,depth,err)

if __name__ == '__main__':
  main(sys.argv)
