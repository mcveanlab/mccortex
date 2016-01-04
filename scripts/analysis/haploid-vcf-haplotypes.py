#!/usr/bin/env python3

from __future__ import print_function

import sys
from itertools import groupby
import vcf

def fasta_iter(file_path):
  """
  Given a fasta file. yield tuples of header, sequence
  author: brentp
  url: https://www.biostars.org/p/710/
  """
  with open(file_path) as fh:
    # ditch the boolean (x[0]) and just keep the header or sequence since
    # we know they alternate.
    faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
    for header in faiter:
      # drop the ">"
      header = next(header)[1:].strip()
      # join all sequence lines to one.
      seq = "".join(s.strip() for s in next(faiter))
      yield header, seq

def load_fasta(path):
  # Load all chroms
  chrs = {}
  try:
    g=fasta_iter(path)
    for (n,s) in g:
      chrs[n] = s.upper()
  except FileNotFoundError as fne:
    print("Cannot find file:",path,file=sys.stderr)
    sys.exit(-1)
  return chrs

class BufferedGenerator:
  def __init__(self, iter):
    self.iter = iter
    self.buffer = []

  def __iter__(self):
    return self

  def next(self):
    if self.buffer:
      return self.buffer.pop() # remove last item
    else:
      return self.iter.next()

  def unnext(self,item):
    self.buffer.append(item)


def vcf_is_biallelic(v):
  return len(v.ALT) == 1

def vcf_trim_biallelic(v):
  ref=v.REF
  alt=v.ALT[0]
  rshift=0
  # left trim
  while(len(ref)>0 and len(alt)>0 and ref[0] == alt[0]):
    ref=ref[1:]
    alt=alt[1:]
    rshift+=1
  # right trim
  while(len(ref)>0 and len(alt)>0 and ref[-1] == alt[-1]):
    ref=ref[0:-1]
    alt=alt[0:-1]
  v.REF=ref
  v.ALT[0]=alt
  v.POS += rshift

def print_contigs(chrs, vcfin, k):
  for v in vcfin:
    if vcf_is_biallelic(v):
      vcf_trim_biallelic(v)


def main(args):
  if len(args) != 4:
    print("usage: %s <ref> <vcf> <k>" % (args[0]))
    sys.exit(-1)

  refpath = args[1]
  vcfin = vcf.Reader(filename=args[2])
  k = int(args[3])

  chrs = load_fasta(refpath)

  print_contigs(chrs, vcfin, k)

if __name__ == '__main__':
  main(sys.argv)
