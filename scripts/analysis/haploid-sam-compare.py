#!/usr/bin/env python2

from __future__ import print_function

import os
import sys
import pysam
import getopt

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
    if str(pairs[i][1]).islower():
      return False

  return True

def pretty_frac(nom,denom):
  return "%9d / %d (%.4f%%)" % (nom, denom, (100.0*nom)/denom)

def haploid_sam_compare(sampath,print_valid):
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
      if print_valid:
        print(read.query_name)

  n_mapped = n_mappings - n_unmapped - n_secondary - n_lowmapq

  print("n_mappings   = %9d" % (n_mappings), file=sys.stderr)
  print("n_unmapped   =",pretty_frac(n_unmapped, n_mappings), file=sys.stderr)
  print("n_secondary  =",pretty_frac(n_secondary, n_mappings), file=sys.stderr)
  print("n_lowmapq<20 =",pretty_frac(n_lowmapq, n_mappings), file=sys.stderr)
  print("n_mapped     =",pretty_frac(n_mapped, n_mappings), file=sys.stderr)
  print("n_correct    =",pretty_frac(n_correct, n_mapped), file=sys.stderr)
  print("n_incorrect  =",pretty_frac(n_mapped-n_correct, n_mapped), file=sys.stderr)

  samfile.close()

def usage(err):
  if(len(err) > 0):
    print(err,file=sys.stderr)
  print("usage: %s [options] <sam>" % (os.path.basename(__file__)),file=sys.stderr)
  print("options:",file=sys.stderr)
  print("  -p, --print-valid     print names of passing contigs",file=sys.stderr)
  sys.exit(-1)

def main(args):
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hp", ["help", "print-valid"])
  except getopt.GetoptError as err:
    # print help information and exit:
    print(err) # will print something like "option -a not recognized"
    usage("")

  print_valid = False

  for o, a in opts:
    if o in ("-h", "--help"):
      usage("")
    elif o in ("-p", "--print-valid"):
      print_valid = True
    else:
      usage("Bad option: %s" % (o))

  if(len(args) > 1):
    usage("Unused arguments")
  if(len(args) < 1):
    usage("")

  sampath = args[0]
  haploid_sam_compare(sampath,print_valid)

if __name__ == '__main__':
  main(sys.argv)
