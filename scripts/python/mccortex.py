#!/usr/bin/env python3

from __future__ import print_function

import os, sys, errno
from itertools import groupby

# def faidx_load_all(faidx):
#   chrs = {}
#   for ref in faidx.references:
#     chrs[ref] = faidx.fetch(reference=ref)
#   return(chrs)

# def load_fasta(path):
#   faidx = pysam.FastaFile(path)
#   chrs = faidx_load_all(faidx)
#   faidx.close()
#   return(chrs)

def mkdir_p(path):
  """
  http://stackoverflow.com/a/600612/431087
  """
  try:
    os.makedirs(path)
  except OSError as exc:
    if not (exc.errno == errno.EEXIST and os.path.isdir(path)): raise

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

def quack():
  print("QUACK!")

def reverse_complement(s):
  complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', \
                'a': 't', 'c': 'g', 'g': 'c', 't': 'a'}
  bases = list(s) # get returns the base if not in dict
  bases = reversed([complement.get(base,base) for base in bases])
  return ''.join(bases)

def dna_key(s):
  r = reverse_complement(s)
  return s if s <= r else r

def load_fasta(path):
  """
  Load all chromosomes from a FASTA file.
  """
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

