#!/usr/bin/env python2

from __future__ import print_function

import sys
import pysam

# Expect '>CHROM:POS:ID:REF:ALT:AIDX:pos:RELPOS:ltrim:LTRIM:rtrim:RTRIM\n'
#   AIDX is allele index (0=>ref, 1=> first ALT, ...)
#   RELPOS is start of ref allele in contig (0-based)
#   LTRIM is number of bases to trim from left side
#   RTRIM is number of bases to trim from right side

def faidx_load_all(faidx):
  chrs = {}
  for ref in faidx.references:
    chrs[ref] = faidx.fetch(reference=ref)
  return(chrs)

def decompose_read_name(read_name):
  title = read_name.split(':')
  ref = title[3]
  alt = title[4]
  relpos = int(title[7])
  ltrim = int(title[9])
  rtrim = int(title[11])
  ref = ref[ltrim:-rtrim]
  alt = alt[ltrim:-rtrim]
  return (ref,alt,relpos)

def validate_var(read):
  (ref,alt,relpos) = decompose_read_name(read.query_name)
  pairs = read.get_aligned_pairs(with_seq=True)
  # dictionary of query position -> ref, ref base
  query = {pair[0]: (pair[1],pair[2]) for pair in pairs}

  # Need a matching base either side and agreement on allele
  for i in xrange(relpos-1, relpos+len(alt)+1):
    if (i not in query) or (query[i][1] is None):
      return False

  # Check distance between flanking bases matches
  dist = query[relpos+len(alt)][0] - query[relpos-1][0]
  if dist != len(alt)+1:
    return False

  # [relpos, relpos+len(alt)) must match
  # mismatches are in lower case
  for i in xrange(relpos, relpos+len(alt)):
    if pairs[i][1].islower():
      return False

  return True

def pretty_frac(nom,denom):
  return "%9d / %d (%.4f%%)" % (nom, denom, (100.0*nom)/denom)

def main(args):
  if len(args) != 2:
    print("usage: %s <sam>" % (args[0]))
    sys.exit(-1)

  sampath = args[1]
  samfile = pysam.AlignmentFile(sampath, "r")

  # refpath = args[2]
  # faidx = pysam.FastaFile(refpath)
  # chrs = faidx_load_all(faidx)
  # faidx.close()
  # chrs[read.reference_name]

  n_mappings = 0
  n_secondary = 0
  n_unmapped = 0
  n_correct = 0
  n_lowmapq = 0

  for read in samfile:
    n_mappings += 1
    if read.is_secondary or read.is_supplementary:
      n_secondary += 1
    elif read.is_unmapped:
      n_unmapped += 1
    elif read.mapping_quality < 20:
      n_lowmapq += 1
    elif validate_var(read):
      n_correct += 1

  n_mapped = n_mappings - n_unmapped - n_secondary - n_lowmapq

  print("n_mappings   = %9d" % (n_mappings))
  print("n_unmapped   =",pretty_frac(n_unmapped, n_mappings))
  print("n_secondary  =",pretty_frac(n_secondary, n_mappings))
  print("n_lowmapq<20 =",pretty_frac(n_lowmapq, n_mappings))
  print("n_mapped     =",pretty_frac(n_mapped, n_mappings))
  print("n_correct    =",pretty_frac(n_correct, n_mapped))

  samfile.close()


if __name__ == '__main__':
  main(sys.argv)
