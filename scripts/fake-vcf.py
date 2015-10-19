#!/usr/bin/env python

from __future__ import print_function
from itertools import groupby
import sys

mut = {'A': 'C','C': 'G','G': 'T','T': 'A'}

hdr = ('##fileformat=VCFv4.2\n'
       '##FILTER=<ID=PASS,Description="All filters passed">\n'
       '##fileDate=20151014\n'
       '##reference=ref/ref.fa\n'
       '##contig=<ID=ref,length=599>\n'
       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
       '##FORMAT=<ID=NK21R,Number=A,Type=Integer,Description="Number of exclusive kmers on ref for each allele (k=21)">\n'
       '##FORMAT=<ID=NK21A,Number=A,Type=Integer,Description="Number of exclusive kmers on alt for each allele (k=21)">\n'
       '##FORMAT=<ID=CK21R,Number=A,Type=Integer,Description="Mean ref exclusive kmer coverage (k=21)">\n'
       '##FORMAT=<ID=CK21A,Number=A,Type=Integer,Description="Mean alt exclusive kmer coverage (k=21)">\n')

col_hdr='#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT\n'

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
      chrs[n] = s
  except FileNotFoundError as fne:
    print("Cannot find file:",ref_path,file=sys.stderr)
    sys.exit(-1)

  # Print header
  print(hdr,end='')
  for name,s in chrs.items():
    print("##contig=<ID=",name,",length=",str(len(s)),">",sep='')
  print(col_hdr,end='')

  # Generate entries
  for name,s in chrs.items():
    for i in range(0,len(s)):
      print(name,"\t",str(i+1),"\t.\t",s[i],"\t",mut[s[i]],"\t.\tPASS\t.\t.",sep='');

def main():
  if len(sys.argv) != 2:
      print("usage: %s <ref.fa>" % (sys.argv[0]))
      sys.exit(-1)

  fake_vcf(sys.argv[1])

if __name__ == '__main__':
    main()
