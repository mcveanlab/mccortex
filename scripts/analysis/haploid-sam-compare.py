#!/usr/bin/env python2

from __future__ import print_function

import os
import sys
import pysam
import getopt

# Expect '>CHROM:POS:ID:REFLEN:ALTLEN:ALTIDX:pos:RELPOS:ltrim:LTRIM:rtrim:RTRIM\n'
#   CHROM,POS,ID are taken from the VCF, POS is 1-based
#   ALTIDX is allele index (0=>ref, 1=> first ALT, ...)
#   RELPOS is start of ref allele in contig (0-based)
#   LTRIM is number of bases trimmed from left side
#   RTRIM is number of bases trimmed from right side

def decompose_read_name(read_name):
  title = read_name.split(':')
  reflen = int(title[3])
  altlen = int(title[4])
  relpos = int(title[7])
  ltrim = int(title[9])
  rtrim = int(title[11])
  reflen -= ltrim + rtrim;
  altlen -= ltrim + rtrim;
  return (reflen,altlen,relpos)

def validate_var(read):
  (reflen,altlen,relpos) = decompose_read_name(read.query_name)
  pairs = read.get_aligned_pairs(with_seq=True)
  # dictionary of query position -> ref, ref base
  query = {pair[0]: (pair[1],pair[2]) for pair in pairs}

  # Need a matching base either side and agreement on allele
  for i in xrange(relpos-1, relpos+altlen+1):
    if (i not in query) or (query[i][1] is None):
      return False

  # Check distance between flanking bases matches
  dist = query[relpos+altlen][0] - query[relpos-1][0]
  if dist != altlen+1:
    return False

  # [relpos, relpos+altlen) must match
  # mismatches are in lower case
  for i in xrange(relpos, relpos+altlen):
    if str(pairs[i][1]).islower():
      return False

  return True

def pretty_frac(nom,denom):
  return "%9d / %d (%.4f%%)" % (nom, denom, (100.0*nom)/denom)

def haploid_sam_compare(sampath,print_valid_fh):
  samfile = pysam.AlignmentFile(sampath, "r")

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
      if print_valid_fh is not None:
        print(read.query_name,file=print_valid_fh)

  n_mapped = n_mappings - n_unmapped - n_secondary - n_lowmapq

  print("n_mappings   = %9d" % (n_mappings))
  print("n_unmapped   =",pretty_frac(n_unmapped, n_mappings))
  print("n_secondary  =",pretty_frac(n_secondary, n_mappings))
  print("n_lowmapq<20 =",pretty_frac(n_lowmapq, n_mappings))
  print("n_mapped     =",pretty_frac(n_mapped, n_mappings))
  print("n_correct    =",pretty_frac(n_correct, n_mapped))
  print("n_incorrect  =",pretty_frac(n_mapped-n_correct, n_mapped))

  samfile.close()

def usage(err):
  if(len(err) > 0):
    print(err,file=sys.stderr)
  print("usage: %s [options] <sam>" % (os.path.basename(__file__)),file=sys.stderr)
  print("options:",file=sys.stderr)
  print("  -p, --print-valid <out.txt>  print names of passing contigs",file=sys.stderr)
  sys.exit(-1)

def main(args):
  try:
    opts, args = getopt.getopt(sys.argv[1:], "hp:", ["help", "print-valid="])
  except getopt.GetoptError as err:
    # print help information and exit:
    print(err) # will print something like "option -a not recognized"
    usage("")

  print_valid=None
  valid_fh=None

  for o, a in opts:
    if o in ("-h", "--help"):
      usage("")
    elif o in ("-p", "--print-valid"):
      print_valid = a
    else:
      usage("Bad option: %s" % (o))

  if(len(args) > 1):
    usage("Unused arguments")
  if(len(args) < 1):
    usage("")

  sampath = args[0]
  if print_valid is not None:
    valid_fh = open(print_valid, 'w')

  haploid_sam_compare(sampath,valid_fh)

  if valid_fh is not None:
    valid_fh.close()

if __name__ == '__main__':
  main(sys.argv)
