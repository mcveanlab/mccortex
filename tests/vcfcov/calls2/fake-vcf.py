#!/usr/bin/env python

from __future__ import print_function
from itertools import groupby
import sys

mut = {'A': 'C','C': 'G','G': 'T','T': 'A'}

hdr = ('##fileformat=VCFv4.2\n'
       '##FILTER=<ID=PASS,Description="All filters passed">\n'
       '##fileDate=20151014\n'
       '##reference=ref/ref.fa\n'
       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
       '##FORMAT=<ID=K21R,Number=A,Type=Integer,Description="Coverage on ref (k=21) => sum(kmer_covs)/exp_num_kmers">\n'
       '##FORMAT=<ID=K21A,Number=A,Type=Integer,Description="Coverage on alt (k=21) => sum(kmer_covs)/exp_num_kmers">\n')

col_hdrs=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']

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

def fake_vcf(ref_path):
  # Load all chroms
  chrs = {}

  try:
    g=fasta_iter(ref_path)
    for (n,s) in g:
      chrs[n] = s.upper()
  except FileNotFoundError as fne:
    print("Cannot find file:",ref_path,file=sys.stderr)
    sys.exit(-1)

  # Print header
  print(hdr,end='')
  for name,s in chrs.items():
    print("##contig=<ID=",name,",length=",str(len(s)),">",sep='')
  print('#','\t'.join(col_hdrs),sep='')

  # Generate entries
  for name,s in chrs.items():
    for i in range(0,len(s)):
      print(name,str(i+1),".",s[i],mut[s[i]],".","PASS",".",".",sep='\t');

def main():
  if len(sys.argv) != 2:
      print("usage: %s <ref.fa>" % (sys.argv[0]))
      sys.exit(-1)

  fake_vcf(sys.argv[1])

if __name__ == '__main__':
    main()
