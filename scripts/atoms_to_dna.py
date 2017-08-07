#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import argv,stdout,stderr,exit
from os.path import basename, join, abspath
from Bio import SeqIO
from Bio.Seq import Seq, Alphabet
from string import maketrans
import csv


TRANS_TABLE=maketrans('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz',\
        'ANCNNNGNNNNNNNNNNNNTNNNNNMancnnngnnnnnnnnnnnntnnnnnn')
DEFAULT_OUT_DIR = '.'

def readSegments(data):
    res = dict()
    for line in csv.reader(data, delimiter='\t'):
        if not res.has_key(line[0]):
            res[line[0]] = list()
        res[line[0]].append(tuple(map(int, line[4:])))
    return res

def paritionRecord(record, segments):
    res = list()
    c = 1
    for start, stop in segments:
        rec = record[start:stop]
        name = record.id.replace('|', '.')
        rec.id = '%s_%s|%s:%s|chromosome|%s' %(name, c, start, stop-1, name)
        rec.description = ''
        rec.seq = Seq(str(rec.seq).translate(TRANS_TABLE),
                Alphabet.generic_dna)
        if rec.seq.count('N') + rec.seq.count('n') != len(rec):
            res.append(rec) 
        c += 1
    return res

if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-o', '--out_dir', default=DEFAULT_OUT_DIR, type=str,
            help='prefix of output filename')
    parser.add_argument('atoms_file', type=str, help='file containing atoms')
    parser.add_argument('fasta_file', type=str, nargs='+', 
            help='fasta files with original genomic sequences')

    args = parser.parse_args()
    
    segments = readSegments(open(args.atoms_file))

    for f in parser.fasta_file:
        outName = join(args.out_dir, basename(f))
        if abspath(outName) == abspath(f):
            print >> stderr, ('Designated output file %s is same as its ' + \
                    'correspondign input file.  Unable to store extracted' + \
                    'sequence data. Exiting.') %abspath(outName)
            exit(1)
        out = open(outName, 'w')
        for record in SeqIO.parse(open(f), 'fasta'):
            if segments.has_key(record.id):
                recs = paritionRecord(record, segments[record.id])
                SeqIO.write(recs, out, 'fasta')
            else:
                print >> stderr, 'No segment found on fasta record %s' %record.id
        out.close()
                 
    
