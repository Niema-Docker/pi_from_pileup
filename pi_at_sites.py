#!/usr/bin/env python3
'''
Given a complete pi_from_pileup output, compute pi from specific positions
'''
from sys import argv, stderr
if __name__ == "__main__":
    if len(argv) != 3:
        print("%s <pi_from_pileup_tsv> <0_based_positions_txt>" % argv[0], file=stderr); exit(1)
    positions = {int(l) for l in open(argv[2])}; pi = 0.; L = 0
    for l in open(argv[1]):
        chrom, pos, D_l = l.split('\t')
        if int(pos) in positions:
            print(l.strip()); pi += float(D_l); L += 1
    if L == 0:
        pi = 0
    else:
        pi /= L
    print("L then PI\t%d\t%s" % (L,pi))
