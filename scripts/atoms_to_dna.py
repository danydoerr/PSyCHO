#!/usr/bin/env python

from sys import argv,stdout,stderr,exit
from os.path import basename, exists
from Bio import SeqIO
from Bio.Seq import Seq, Alphabet
from string import maketrans
import csv


TRANS_TABLE=maketrans('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz',\
        'ANCNNNGNNNNNNNNNNNNTNNNNNMancnnngnnnnnnnnnnnntnnnnnn')

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

    if len(argv) < 3:
        print '\tusage: %s <ATOM FILE> <FASTA FILE 1> ... <FASTA FILE N>' %argv[0]
        exit(1)
    
    segments = readSegments(open(argv[1]))

    for f in argv[2:]:
        outName = basename(f)
        if exists(outName):
            print >> stderr, ('File %s already exists. Unable to store ' + \
                    'extracted sequence data. Exiting.') %outName
            exit(1)
        out = open(outName, 'w')
        for record in SeqIO.parse(open(f), 'fasta'):
            if segments.has_key(record.id):
                recs = paritionRecord(record, segments[record.id])
                SeqIO.write(recs, out, 'fasta')
            else:
                print >> stderr, 'No segment found on fasta record %s' %record.id
        out.close()
                 
    
