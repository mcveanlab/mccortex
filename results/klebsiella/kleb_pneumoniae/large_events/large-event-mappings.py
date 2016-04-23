#!/usr/bin/env python2

from __future__ import print_function

import os
import sys
import pysam
import getopt

# Update sys.path to include paths relative to the script dir
def add_rel_path_to_sys_path(rel_inc_path):
  import sys,os,os.path
  scriptdir = os.path.dirname(os.path.realpath(__file__))
  incdir = os.path.realpath(scriptdir+'/'+rel_inc_path)
  sys.path.append(incdir)

# include mccortex/scripts/python/
add_rel_path_to_sys_path('../../../../scripts/python')
import mccortex

def usage(err=None):
  if err is not None:
    print(err,file=sys.stderr)
  print("usage: %s <sam1> <sam2>" % (os.path.basename(__file__)),file=sys.stderr)
  print("  Anaylsis of bubble contigs mapped to two different references", file=sys.stderr)
  sys.exit(-1)

def read_primary(samfh,reads):
  while True:
    try:
      read = next(samfh)
      if (not read.is_secondary) and (not read.is_supplementary):
        reads.append(read)
        break
    except StopIteration:
      break

# Read two entries from 
def read_next_pairs(samfh1,samfh2):
  rds = []
  read_primary(samfh1,rds)
  read_primary(samfh1,rds)
  read_primary(samfh2,rds)
  read_primary(samfh2,rds)

  if (len(rds) != 0) and (len(rds) != 4):
    raise SystemExit("One file is shorter than the other: %i" % (len(rds)))
  else:
    return(rds)

class MappingComparisons:
  # Parameters
  min_len_diff = 100

  # Stats
  nbubbles = 0
  n_multiallelic_skipped = 0
  n_len_diff_small = 0
  pair_non_alternate = 0
  pair_no_mapping = 0
  n_good_pairs = 0
  n_perfect_pairs = 0

  mappings = [] # list of tuples: (m1,m2,m3,m4) m1==m3, m2==m4, m1,m2 map to f1
  maxdiff = 0

# m1, m2 are lists of length 2
# m1 alignments are written to fh1,
# m2 alignments are written to fh2
def write_paired_alignments(fh1,fh2,m1,m2):
    fh1.write(m1[0])
    fh1.write(m1[1])
    fh2.write(m2[0])
    fh2.write(m2[1])

def mapping_get_orig_seq(m):
  return mccortex.reverse_complement(m.seq) if m.is_reverse else m.seq

# return (a,b) or (r(a),r(b)), whichever is less
# where r(x) is reverse_complement()
def mapping_key_pair(a,b):
  ra = mccortex.reverse_complement(a)
  rb = mccortex.reverse_complement(b)
  return (a,b) if (a < ra or (a == ra and b <= rb)) else (ra,rb)

# BAM_CMATCH is 0
def good_mapping(m):
  return(not m.is_unmapped and \
         len(m.cigartuples) == 1 and \
         m.cigartuples[0][0] == 0 and \
         m.cigartuples[0][1] == len(m.seq))

# Which of two mapping map with all bases matching
# 0 => neither, 1 => first, 2 => second, 3 => both
def get_mapping_code(a,b):
  return good_mapping(a) + good_mapping(b)*2

def compare_mappings(m1,m2,stats,bad1,bad2):
  """
  m1, m2 are lists of bubble branches, mapped to two different assemblies
  """
  stats.nbubbles += 1

  # Only deal with biallic
  if len(m1) != 2 or len(m2) != 2:
    stats.n_multiallelic_skipped += 1
    return

  # print(m1[0].query_name, m2[0].query_name, m1[0].seq, m2[0].seq)

  assert(m1[0].query_name == m2[0].query_name and len(m1[0].seq) == len(m2[0].seq))
  assert(m1[1].query_name == m2[1].query_name and len(m1[1].seq) == len(m2[1].seq))

  # Check mappings order
  # if m1[0].query_name > m1[0].query_name:
  #   raise SystemExit("Bad order: "+m1[0].query_name+" vs "+m1[1].query_name)
  # if m2[0].query_name > m2[0].query_name:
  #   raise SystemExit("Bad order: "+m2[0].query_name+" vs "+m2[1].query_name)

  # if (m1[0].is_unmapped == m1[1].is_unmapped) or (m1[0].is_unmapped == m1[1].is_unmapped):
  #   stats.non_alternate_pairs += 1
  #   return

  lendiff = abs(len(m1[0].seq) - len(m1[1].seq))
  stats.maxdiff = max(stats.maxdiff, lendiff)
  if lendiff < stats.min_len_diff:
    stats.n_len_diff_small += 1
    return

  # Must be alternate
  f1code = get_mapping_code(m1[0],m1[1])
  f2code = get_mapping_code(m2[0],m2[1])
  good_mapped_pairs = (f1code & f2code) == 0 and f1code>0 and f2code>0

  if not f1code and not f2code:
    stats.pair_no_mapping += 1
    return

  if not good_mapped_pairs:
    stats.pair_non_alternate += 1
    write_paired_alignments(bad1,bad2,m1,m2)
    return

  # Sort so alignment[0] maps to file[0], and alignment[1] to file[1]
  if f1code == 2 and f2code == 1:
    f1code, f2code = f2code, f1code
    m1[0,1] = m1[1,0]
    m2[0,1] = m2[1,0]

  # First read only maps to file1, second read only maps to file2
  assert(f1code == 1 and f2code == 2)

  # Got a good pair!
  stats.n_good_pairs += 1
  # write_paired_alignments(good1,good2,m1,m2)

  pairs1 = m1[0].get_aligned_pairs(with_seq=True)
  pairs2 = m2[1].get_aligned_pairs(with_seq=True)
  if True not in { str(pair[2]).islower() for pair in (pairs1+pairs2) }:
    stats.n_perfect_pairs += 1

  s1 = mapping_get_orig_seq(m1[0])
  s2 = mapping_get_orig_seq(m1[1])
  s1, s2 = mapping_key_pair(s1,s2)
  stats.mappings.append((m1[0], m1[1], m2[0], m2[1], s1, s2))

  # print(mapping_get_orig_seq(m1[0]), mapping_get_orig_seq(m1[1]))

def main(args):

  if len(args) != 4:
    usage()

  sam1 = args[1]
  sam2 = args[2]
  outdir = args[3]

  print("reading %s and %s" % (sam1, sam2), file=sys.stderr)
  print("writing to dir %s" % (outdir), file=sys.stderr)

  # Make output directory
  # os.makedirs(name, exist_ok=False)
  # os.makedirs(outdir, exist_ok=True)
  mccortex.mkdir_p(outdir)

  stats = MappingComparisons()

  samfhs = []
  samfhs.append(pysam.AlignmentFile(sam1, "r"))
  samfhs.append(pysam.AlignmentFile(sam2, "r"))

  # wh means: open for writing as SAM, write the header
  good1fh = pysam.AlignmentFile(outdir+'/good1.sam', 'wh', template=samfhs[0])
  good2fh = pysam.AlignmentFile(outdir+'/good2.sam', 'wh', template=samfhs[1])
  bad1fh = pysam.AlignmentFile(outdir+'/bad1.sam', 'wh', template=samfhs[0])
  bad2fh = pysam.AlignmentFile(outdir+'/bad2.sam', 'wh', template=samfhs[1])

  # read two reads from each sam file
  while True:
    r = read_next_pairs(samfhs[0], samfhs[1])
    if len(r) == 0: break
    mid = len(r)/2
    compare_mappings(r[0:mid], r[mid:len(r)], stats, bad1fh, bad2fh)

  print("read %i bubbles" % (stats.nbubbles), file=sys.stderr)
  print("skipped %i multiallelic" % (stats.n_multiallelic_skipped), file=sys.stderr)
  print("len diff too small: %i" % (stats.n_len_diff_small), file=sys.stderr)
  print("pair_no_mapping: %i" % (stats.pair_no_mapping), file=sys.stderr)
  print("pair_non_alternate: %i" % (stats.pair_non_alternate), file=sys.stderr)
  print("n_good_pairs: %i" % (stats.n_good_pairs), file=sys.stderr)
  print("n_perfect_pairs: %i" % (stats.n_perfect_pairs), file=sys.stderr)

  print("maxdiff: %i" % (stats.maxdiff), file=sys.stderr)

  events = []

  if len(stats.mappings) > 0:
    # Sort mappings
    stats.mappings.sort(key=lambda t: (t[4], t[5]))

    # remove duplicates
    p = stats.mappings[0]
    events.append(p)
    for i in xrange(1,len(stats.mappings)):
      n = stats.mappings[i]
      if p[4] != n[4] or p[5] != n[5]:
        p = n
        events.append(p)

  # write good events
  for t in events:
    write_paired_alignments(good1fh,good2fh,[t[0],t[1]],[t[2],t[3]])

  # Write histogram of sizes, removing flank sizes
  histfh = open(outdir+'/stats.txt', 'w')
  for t in events:
    fields = t[0].query_name.split(':')
    if len(fields) < 3: raise SystemExit("Name missing flank sizes: "+t[0].query_name)
    flanks = int(fields[1]) + int(fields[2])
    print("%i\t%i" % (len(t[0].seq)-flanks, len(t[1].seq)-flanks), file=histfh)
  histfh.close()

  print("after rmdup, wrote %i events" % (len(events)))

  good1fh.close()
  good2fh.close()
  bad1fh.close()
  bad2fh.close()

  samfhs[0].close()
  samfhs[1].close()

if __name__ == '__main__':
  main(sys.argv)
