#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, exit
from os.path import basename 
from Bio import SeqIO
import re

DEFAULT_WINDOW_SIZE = 1000000

if __name__ == '__main__':
    
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-w', '--window_size', default=DEFAULT_WINDOW_SIZE,
            type=int, help='window size of the histogram')
    parser.add_argument('fasta_file', type=str, nargs='+',
            help='fasta file containing one or more sequence records')

    args = parser.parse_args()

    for f in args.fasta_file:
        gid = basename(f).split('.', 1)[0].lower()

        for rec in SeqIO.parse(f, 'fasta'):
            chr1 = rec.id.lower()

            band_id = '.'.join((gid, chr1))

            c = i = 0
            for s in rec.seq:
                if s.upper() not in ('A', 'C', 'G', 'T'):
                    c += 1
                i += 1
                if i % args.window_size == 0:
                    print '%s %s %s %.5f' %(band_id, i-args.window_size, i-1,
                            c/float(args.window_size))  
                    c = 0
            if i % args.window_size != 0:
                l = i % args.window_size
                print '%s %s %s %.5f' %(band_id, i-l, i,
                        c/float(l))  
    
