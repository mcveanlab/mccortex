#!/usr/bin/env python
from __future__ import print_function
try: input = raw_input
except: pass

# usage: python break-contigs-vs-truth.py <k> [input.txt]
# input.txt:
#  <master> :string
#  <query1> :string
#  <query2> :string
#  ...
#
# Output:
#  <contig-id> <ref-strand> <ref-start> <contig-substr>
#  ...
#
# For each query string, report all maximal substring alignments between it and
# the master string. Alignments are of length >= k and can include single base
# mismatches as long as they are flanked by k matching bases.
#
# Next we print coverage of ref using maximal alignments to contigs.
# Finally we report NG50 and number of assembly errors.
#
# Time is on average NlogN where N is number of query bases. Worst case N^2
#  e.g. aaaaaaaaaaaa vs aaaaaaaaaaaa
# Memory usage is linear with length of the master string and number of maximal
# substring matches.
#
# Isaac Turner 2016-09-10
# MIT License


import fileinput
import pyRBT
import sys
from collections import defaultdict

class Alignment:
  def __init__(self,seqid,start1,start2,length):
    self.seqid = seqid
    self.start1,self.start2,self.length = start1,start2,length
  def __cmp__(x,y):
    if x.seqid != y.seqid: return x.seqid - y.seqid
    xoff,yoff = x.start1-x.start2,y.start1-y.start2
    if xoff != yoff: return xoff - yoff
    if x.start1 >= y.start1+y.length: return 1
    if y.start1 >= x.start1+x.length: return -1
    return 0 # one is within the other
  def __lt__(x,y): return x.__cmp__(y)  > 0
  def __ge__(x,y): return x.__cmp__(y) >= 0
  def __eq__(x,y): return x.__cmp__(y) == 0
  def __ne__(x,y): return x.__cmp__(y) != 0
  def __le__(x,y): return x.__cmp__(y) <= 0
  def __lt__(x,y): return x.__cmp__(y)  < 0
  def __str__(self):
    return ("Alignment(id:"+str(self.seqid)+"," +
            "starts:("+str(self.start1)+","+str(self.start2)+")," +
            "len:"+str(self.length)+")")

def dna_reverse_complement(s):
  h = {'a':'t','c':'g','g':'c','t':'a','A':'T','C':'G','G':'C','T':'A'}
  t = [ h.get(c,c) for c in s ]
  return ''.join(t[::-1])

# upper case
def build_occ_hash(occ,k,seq,seqid):
  for i in range(len(seq)-k+1):
    occ[seq[i:i+k].upper()].append((seqid,i))

# Extend a match as far as we have an exact match or a mismatch flanked by k
# exact matches either side. Seeded from an initial matching kmer
def extend_match_kmer_match(a,b,s1,e1,s2,e2,M):
  while True:
    while s1 > 0 and s2 > 0 and a[s1-1].upper() == b[s2-1].upper(): s1,s2 = s1-1,s2-1
    if s1 > M and s2 > M and a[s1-M-1:s1-1].upper() == b[s2-M-1:s2-1].upper():
      s1 -= M+1; s2 -= M+1
    else: break
  while True:
    while e1 < len(a) and e2 < len(b) and a[e1].upper() == b[e2].upper(): e1,e2 = e1+1,e2+1
    if e1+M < len(a) and e2+M < len(b) and a[e1+1:e1+1+M].upper() == b[e2+1:e2+1+M].upper():
      e1 += M+1; e2 += M+1
    else: break
  return (s1,e1,s2,e2)

def count_str_matches(a,b):
  return sum([ i.upper() == j.upper() for i,j in zip(a,b) ])

# Extend a match as far as we have an exact match or N matches within M bases
def extend_match_kmer_match2(a,b,s1,e1,s2,e2,M,N):
  while True:
    while s1 > 0 and s2 > 0 and a[s1-1].upper() == b[s2-1].upper(): s1,s2 = s1-1,s2-1
    if s1 > M and s2 > M and count_str_matches(a[s1-1-M:s1-1], b[s2-1-M:s2-1]) >= N:
      s1 -= M+1; s2 -= M+1
    else: break
  while True:
    while e1 < len(a) and e2 < len(b) and a[e1].upper() == b[e2].upper(): e1,e2 = e1+1,e2+1
    if e1+M < len(a) and e2+M < len(b) and count_str_matches(a[e1+1:e1+1+M], b[e2+1:e2+1+M]) >= N:
      e1 += M+1; e2 += M+1
    else: break
  return (s1,e1,s2,e2)

# Segment tree
# Query if there is an uncovered region in a set of intervals
class GapSegNode:
  def __init__(self,start,end,l=None,r=None):
    (self.start,self.end,self.l,self.r,self.regs) = (start,end,l,r,[])
  def add_interval(self,reg):
    if reg[0] >= self.end or reg[1] <= self.start: return # does not overlap
    if reg[0] <= self.start and reg[1] >= self.end: # contained
      self.regs.append(reg)
    else:
      if self.l is not None: self.l.add_interval(reg)
      if self.r is not None: self.r.add_interval(reg)
  def gap_in_interval(self,reg):
    if reg[0] >= self.end or reg[1] <= self.start: return False # does not overlap
    if len(self.regs) > 0: return False # this node is covered
    if self.l is None and self.r is None: return True # leaf node
    return ((self.l is not None and self.l.gap_in_interval(reg)) or
            (self.r is not None and self.r.gap_in_interval(reg)))
  def __str__(self):
    return "GapSegNode("+str(self.start)+","+str(self.end)+")"
  @staticmethod
  def build_tree(start,end):
    # build tree bottom up
    assert start < end
    l = [ GapSegNode(i,i+1) for i in range(start,end) ]
    while len(l) > 1:
      N = 2*int(len(l)//2)
      m = [ GapSegNode(l[i].start,l[i+1].end,l[i],l[i+1]) for i in range(0,N,2) ]
      if len(l)%2 == 1: m.append(GapSegNode(l[-1].start, l[-1].end, l[-1]))
      l = m
    return l[0]

# `l` is a list of alignments of sequence `seq` against the ref
# for each kmer in seq, keep the longest substring in l that covers it
# discard all substrings that are not kept
def reduce_alignments(seq,l):
  # sort alignments by length (longest to shortest)
  l.sort(key=lambda x: -x.length)
  # iterate over alignments, only taking those that cover uncovered kmers
  gst = GapSegNode.build_tree(0,len(seq))
  keep = []
  for aln in l:
    reg = (aln.start2, aln.start2+aln.length)
    if gst.gap_in_interval(reg):
      keep.append(aln)
      gst.add_interval(reg)
  return keep

# l is a list of (start,length) contig alignments
# removes substrings from l (BEWARE: they get deleted!)
# returns contig_index,contig_length
def ng50_from_coverage(l,reflen):
  # sort by start (ascending) and length (descending) to remove substrings
  # on the ref
  l.sort(key=lambda x: (x[0],-x[1]))
  j = end = 0
  for x in l:
    if x[0]+x[1] > end:
      end = x[0]+x[1]
      l[j] = x
      j += 1
  del(l[j:])
  # sort by length (descending), start position (ascending)
  l.sort(key=lambda x: (-x[1],x[0]))
  halflen = reflen//2
  lensum = n = 0
  while lensum < halflen:
    if n+1 == len(l):
      print("Warning: haven't assembled half of ref, NG50 is underestimate",
            file=sys.stderr)
      break
    lensum += l[n][1]
    n += 1
  return (n,l[n][1])

# for a given alignment, get left hand position on + strand
def convert_ref_strandpos(aln,reflen):
  return aln.start1 if aln.seqid&1 == 0 else reflen-(aln.start1+aln.length)

def main(k,path):
  master = input().strip()
  masterrc = dna_reverse_complement(master)
  print(master)
  print(masterrc)
  occ = defaultdict(list) # [ (strand,pos) ... ]
  build_occ_hash(occ,k,master,0)
  build_occ_hash(occ,k,masterrc,1)
  n_asm_errors = n_no_matches = n_output = 0
  m_cov = []
  strands = ["+","-"]
  print("# Matching contigs sections")
  print("# contig-id contig-substr ref-start ref-strand")
  # Iterate over queries
  qi = 0
  for q in fileinput.input(path):
    q = q.strip()
    rbt = pyRBT.pyRBT()
    l = []
    for i in range(len(q)-k+1):
      kmer = q[i:i+k].upper()
      for (seqid,pos) in occ[kmer]:
        # create alignment
        aln = Alignment(seqid,pos,i,k)
        if aln not in rbt:
          # extend alignment
          s1,s2 = pos,i
          m = master if seqid==0 else masterrc
          (s1,e1,s2,e2) = extend_match_kmer_match(m,q,s1,s1+k,s2,s2+k,5)
          # (s1,e1,s2,e2) = extend_match_kmer_match2(m,q,s1,s1+k,s2,s2+k,10,6)
          aln.start1,aln.start2,aln.length = s1,s2,e1-s1
          # store and print
          rbt.insert(aln)
          l.append(aln)
    for x in l:
      s = convert_ref_strandpos(x,len(master))
      m_cov.append((s,x.length,qi,x.seqid&1,x.start2))
    l = reduce_alignments(q,l)
    for x in l:
      s = convert_ref_strandpos(x,len(master))
      print(qi,q[x.start2:x.start2+x.length],s,strands[x.seqid&1])
    n_asm_errors += max(0,len(l)-1)
    n_no_matches += (len(l) == 0)
    n_output += len(l)
    qi += 1
  (_,ng50) = ng50_from_coverage(m_cov,len(master))
  print() # empty line separates output, now print ref matches
  print("# Ref positions assembled (longest to shortest)")
  print("# start positions are 0-based and indicate the left hand position")
  print("# ref-start length contig-id contig-start contig-strand")
  m_cov.sort() # sort by start, length (both ascending)
  for x in m_cov: print(x[0],x[1],x[2],x[4],strands[x[3]])
  print("contigs_read:",qi,file=sys.stderr)
  print("contigs_printed:",n_output,file=sys.stderr)
  print("assembly_errors:",n_asm_errors,file=sys.stderr)
  print("nomatch_contigs:",n_no_matches,file=sys.stderr)
  # number of unique segments of the reference that were assembled
  print("num_uniq_ref_segs:",len(m_cov),file=sys.stderr)
  print("reflen:",len(master),file=sys.stderr)
  print("NG50:",ng50,file=sys.stderr)

def usage(err=None):
  if err is not None: print(err,file=sys.stderr)
  print("python break-contigs-vs-truth.py <k> [contigs.txt]",file=sys.stderr)
  exit(-1)

if __name__ == '__main__':
  if len(sys.argv) < 2 or len(sys.argv) > 3: usage()
  try: k = int(sys.argv[1])
  except ValueError: usage("Error: invalid kmer value '"+sys.argv[1]+"'")
  path = sys.argv[2] if len(sys.argv) > 2 else "-"
  print("k:",k,"path:",path,file=sys.stderr)
  main(k,path)
