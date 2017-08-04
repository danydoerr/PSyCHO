#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, stdin, exit
from stat import S_ISFIFO
import os

from Bio.Blast import NCBIStandalone

def parseData(data, out):
    parser = NCBIStandalone.BlastParser()
    it = NCBIStandalone.Iterator(data, parser) 

    for blast_record in it:
        print >> out, blast_record


#matches - Number of matching bases that aren't repeats.
#misMatches - Number of bases that don't match.
#repMatches - Number of matching bases that are part of repeats.
#nCount - Number of 'N' bases.
#qNumInsert - Number of inserts in query.
#qBaseInsert - Number of bases inserted into query.
#tNumInsert - Number of inserts in target.
#tBaseInsert - Number of bases inserted into target.
#strand - defined as + (forward) or - (reverse) for query strand. In mouse, a second '+' or '-' indecates genomic strand.
#qName - Query sequence name.
#qSize - Query sequence size.
#qStart - Alignment start position in query.
#qEnd - Alignment end position in query.
#tName - Target sequence name.
#tSize - Target sequence size.
#tStart - Alignment start position in target.
#tEnd - Alignment end position in target.
#blockCount - Number of blocks in the alignment.
#blockSizes - Comma-separated list of sizes of each block.
#qStarts - Comma-separated list of start position of each block in query.
#tStarts - Comma-separated list of start position of each block in target.



if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
#    parser.add_argument('-o', '--out_filename_prefix',
#            default=DEFAULT_OUT_FILENAME_PREFIX, type=str, 
#            help='prefix of output filename')
    parser.add_argument('-i', '--input', type=str, 
            help='Input BLAST file (by default, input is read from stdin')

    args = parser.parse_args()

    data = stdin.read()

    if not S_ISFIFO(os.fstat(0).st_mode) and args.input:
        data = open(args.input)
    elif not S_ISFIFO(os.fstat(0).st_mode):
        print >> stderr, 'FATAL: Input required\n'
        parser.print_usage()
        exit()

    parseData(data, stdout)
