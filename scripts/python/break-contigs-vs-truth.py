from __future__ import print_function
try: input = raw_input
except: pass

# Input:
#  <K> :int
#  <master> :string
#  <query1> :string
#  <query2> :string
#  ...
#
# Output:
#  <qid> <matching-str>
#  ...
#  <nqueries> <nbreaks> <num-zero-match-queries>
# For each query string, report all maximal substrings of length >=k common
# between it and the master string.
# Time is on average NlogN where N is number of query bases. Worst case N^2
#  e.g. aaaaaaaaaaaa vs aaaaaaaaaaaa
# Memory usage is linear with length of the master string and number of maximal
# substring matches.
#
# Isaac Turner 2016-09-08
# MIT License

import fileinput
import pyRBT
from collections import defaultdict

class Alignment:
  def __init__(self,strand,start1,start2,length):
    self.strand = strand
    self.start1,self.start2,self.length = start1,start2,length
  def __cmp__(x,y):
    if x.strand != y.strand: return x.strand - y.strand
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

def dna_reverse_complement(s):
  h = {'a':'t','c':'g','g':'c','t':'a','A':'T','C':'G','G':'C','T':'A'}
  t = [ h[c] if c in h else c for c in s ]
  return ''.join(t[::-1])

k = int(input().strip())
master = input().strip()

occ = defaultdict(list) # [ (strand,pos) ... ]
for i in range(len(master)-k+1):
  occ[master[i:i+k]].append((True,i))
masterrc = dna_reverse_complement(master)
for i in range(len(masterrc)-k+1):
  occ[masterrc[i:i+k]].append((False,i))

def dna_is_substring(k,occ,l,i,s):
  k0 = s[:k]
  for (idx,pos) in occ[k0]:
    if (idx//2 != i and
        (len(l[idx]) > len(s) or
         (len(l[idx]) == len(s) and idx//2 < i))):
      if l[idx][pos:pos+len(s)] == s: return True
  return False

# Remove substrings from set of strings using kmer-based approach
def remove_dna_substrings(k,strs):
  occ = defaultdict(list) # [ (stridx,pos) ... ]
  l = [] # contains forward and reverse strings
  for s in strs:
    for pos in range(len(s)-k+1): occ[s[pos:pos+k]].append((len(l),pos))
    l.append(s)
    s = dna_reverse_complement(s)
    for pos in range(len(s)-k+1): occ[s[pos:pos+k]].append((len(l),pos))
    l.append(s)
  m = []
  for i,s in enumerate(strs):
    if not dna_is_substring(k,occ,l,i,s):
      m.append(s)
  return m

# remove_dna_substrings(3,['acd','acdf','xxy','xxy'])

# Extend a match as far as we have an exact match or a mismatch flanked by k
# exact matches either side. Seeded from an initial matching kmer
def extend_match_kmer_match(a,b,s1,s2,k):
  e1,e2 = s1+k,s2+k
  while True:
    while s1 > 0 and s2 > 0 and a[s1-1] == b[s2-1]: s1,s2 = s1-1,s2-1
    if s1 > k and s2 > k and a[s1-k-1:s1-1] == b[s2-k-1:s2-1]:
      s1 -= k+1; s2 -= k+1
    else: break
  while True:
    while e1 < len(a) and e2 < len(b) and a[e1] == b[e2]: e1,e2 = e1+1,e2+1
    if e1+k < len(a) and e2+k < len(b) and a[e1+1:e1+1+k] == b[e2+1:e2+1+k]:
      e1 += k+1; e2 += k+1
    else: break
  return (s1,e1,s2,e2)

nbreaks = 0
nomatches = 0
noutput = 0

# Iterate over queries
qi = 0
for q in fileinput.input():
  q = q.strip()
  rbt = pyRBT.pyRBT()
  l = []
  for i in range(len(q)-k+1):
    kmer = q[i:i+k]
    for (fwstrand,pos) in occ[kmer]:
      # create alignment
      aln = Alignment(fwstrand,pos,i,k)
      if aln not in rbt:
        # extend alignment
        s1,s2 = pos,i
        m = master if fwstrand else masterrc
        (s1,e1,s2,e2) = extend_match_kmer_match(m,q,s1,s2,k)
        aln.start1,aln.start2,aln.length = s1,s2,e1-s1
        # store and print
        rbt.insert(aln)
        l.append(q[s2:e2])
  l = remove_dna_substrings(k,l)
  for s in l: print(qi,s)
  nbreaks += 1 if len(l) == 0 else len(l)-1
  nomatches += (len(l) == 0)
  noutput += len(l)
  qi += 1

print(qi,noutput,nbreaks,nomatches)
