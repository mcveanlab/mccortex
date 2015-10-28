#!/usr/bin/env python

from __future__ import print_function
from itertools import groupby
import sys
import random

mut = {'A': 'C','C': 'G','G': 'T','T': 'A'}

hdr = ('##fileformat=VCFv4.2\n'
       '##FILTER=<ID=PASS,Description="All filters passed">\n'
       '##fileDate=20151014\n'
       '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
       '##FORMAT=<ID=NK21R,Number=A,Type=Integer,Description="Number of exclusive kmers on ref for each allele (k=21)">\n'
       '##FORMAT=<ID=NK21A,Number=A,Type=Integer,Description="Number of exclusive kmers on alt for each allele (k=21)">\n'
       '##FORMAT=<ID=CK21R,Number=A,Type=Integer,Description="Mean ref exclusive kmer coverage (k=21)">\n'
       '##FORMAT=<ID=CK21A,Number=A,Type=Integer,Description="Mean alt exclusive kmer coverage (k=21)">\n')

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

def rand_allele_len():
  x = random.randrange(1,30)
  return max(x,21)-20

# ensure first and last bases of ref and alt alleles do not match
def mut_ref_alt(ref,alt):
  m=dict(mut)
  del m[ref[0]]
  if ref[0] != ref[-1]:
    del m[ref[-1]]
  b=list(m.keys())
  random.shuffle(b)
  alt[0] = b[0]
  alt[-1] = b[-1]
  return alt

def make_alt(ref):
  rlen=len(ref)
  alen=rand_allele_len()
  alt=[]
  for i in range(0,alen):
    alt.append("ACGT"[random.randrange(0,4)])
  # Make sure first and last bases don't match the ref allele
  alt = mut_ref_alt(ref,alt)
  return ''.join(alt)

def get_sample_fields(ks):
  k=str(ks)
  return "NK"+k+"R:NK"+k+"A:CK"+k+"R:CK"+k+"A"

# rcov,acov are coverage numbers for the ref/alt alleles respectively
def make_vcfcov(pos,chrlen,ref,alt,ks,rcov,acov):
  s=max(0,pos-ks+1)
  e=min(pos+len(ref)+ks-1,chrlen)
  nr=(e-s-ks+1) * (rcov >= 1)
  na=(e-s-len(ref)+len(alt)-ks+1) * (acov >= 1)
  return ':'.join([str(i) for i in [nr,na,rcov,acov]])

def spaced_vars(chrs,ks,rcov,acov):
  random.seed()
  gts=get_sample_fields(ks)
  # Generate entries
  for name,s in chrs.items():
    althap=""
    lastp=0
    pos=random.randrange(0,10)
    while pos < len(s):
      rlen=min(rand_allele_len(),len(s)-pos)
      ref=s[pos:(pos+rlen)]
      alt=make_alt(ref)
      althap+=s[lastp:pos]+alt.lower()
      vcfcov=make_vcfcov(pos,len(s),ref,alt,ks,rcov,acov)
      print(name,str(pos+1),".",ref,alt,".","PASS",".",gts,vcfcov,sep="\t");
      lastp=pos+rlen
      pos+=rlen+ks-1
    althap+=s[lastp:len(s)]
    for i in range(0,rcov):
      print(">",name,"_ref\n",s,"\n",sep='',end='',file=sys.stderr)
    for i in range(0,acov):
      print(">",name,"_alt\n",althap,"\n",sep='',end='',file=sys.stderr)

def fake_vcf(ref_path,ks,sample,rcov,acov):
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
  print('##ref=',ref_path,sep='')
  for name,s in chrs.items():
    print("##contig=<ID=",name,",length=",str(len(s)),">",sep='')
  col_hdrs.append(sample)
  print('#','\t'.join(col_hdrs),sep='')

  spaced_vars(chrs,ks,rcov,acov)

def main():
  if len(sys.argv) != 6:
      print("usage: %s <ref.fa> <kmer-size> <sample> <rcov> <acov>" % (sys.argv[0]))
      sys.exit(-1)

  fake_vcf(sys.argv[1], int(sys.argv[2]), sys.argv[3],
           int(sys.argv[4]), int(sys.argv[5]))

if __name__ == '__main__':
    main()
