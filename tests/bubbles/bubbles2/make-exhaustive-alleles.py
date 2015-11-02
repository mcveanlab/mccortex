#!/usr/bin/env python

from __future__ import print_function
from itertools import groupby
import sys
import random as rand

mut = {'A': 'C','C': 'G','G': 'T','T': 'A'}

# 00 00 00
# 01 00 00
# 11 00 00
# 00 01 00
# 01 01 00
# 11 01 00
# 00 11 00
# 01 11 00
# 11 11 00

# 00 00 01
# 01 00 01
# 11 00 01
# 00 01 01
# 01 01 01
# 11 01 01
# 00 11 01
# 01 11 01
# 11 11 01

# GGCCGATTAGGCTTCCCGGTAGATTTTGAGCGGCAGCCCTTCATCTTGTACAGTCTACGGATA
# AATTCAGTTTGTCGGCTCTAAGCGGCTTACCTGAAGTCCTAATTCTTCTCCGTATACCCCCTG
# CTATGTCAGCGCTGGTGTATACCACTAAAAGAAATATCAGTTTCGCCAGAGAGGATTCTCACA
# TATTTAGTTTCCGACGGCATGCGTCCTCCCTTATGGATCACGATTAATCGACGCATATCAGAT
# AGTGCTTTGACTATATGCGCAACTTAAAATGCCGAAGCATCGCGTGGTAGCCCCGAGCGAGGT
# GGTTTTGTCTCCCAGCCATCAACGCTATTAATTTGTAAAGCAGATGCCCGACTTTTCCAAGCG
# GTTTAGGATCGGTGTACTGGAGCAATTTGCATCGACCCACTGTCTGTGCGGGTGATCCTGCGG
# TACCACACGTCTTGCGTGTTTTGGCACGTCAGCTGTGAACTGACTGCTTCCGGACGCGTCGGA

def rand_seq(n):
    s=[]
    for i in range(0,n):
        s.append("ACGT"[rand.getrandbits(2)])
    return ''.join(s)

def write_seqs(fh,ref,alt,gt,s):
    print('>'+bin(gt)[2:].zfill(6)+'_'+str(s)+'a', file=fh)
    print(alt if gt & (1<<(2*s))   else ref, file=fh)
    print('>'+bin(gt)[2:].zfill(6)+'_'+str(s)+'b', file=fh)
    print(alt if gt & (1<<(2*s+1)) else ref, file=fh)

def generate_genomes(fhs,k):
    for gt in range(0,2**6):
        seq = rand_seq(2*k+1)
        alt = seq[0:k]+mut[seq[k]]+seq[k+1:]
        for s in range(0,3):
            write_seqs(fhs[s], seq, alt, gt, s)

def main():
    if len(sys.argv) != 4:
        print("usage: %s <a.fa> <b.fa> <c.fa>" % (sys.argv[0]))
        sys.exit(-1)

    fhs = []
    for s in range(0,3):
        fh = open(sys.argv[1+s], 'w')
        fhs.append(fh)

    generate_genomes(fhs, 31)

    for fh in fhs:
        fh.close()

if __name__ == '__main__':
    main()
