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
        res[line[0]].append((int(line[4]), int(line[5]), line[2]))
    return res

def paritionRecord(record, segments, fams):
    res = list()
    c = 1
    for start, stop, f in segments:
        rec = record[start:stop]
        name = record.id.replace('|', '.')
        rec.id = '%s_%s|%s:%s|chromosome|%s' %(name, c, start, stop-1, name)
        rec.description = ''
        rec.seq = Seq(str(rec.seq).translate(TRANS_TABLE),
                Alphabet.generic_dna)
        if rec.seq.count('N') + rec.seq.count('n') != len(rec):
            res.append((star, rec)) 
            if not fams.has_key(f):
                fams[f] = dict()
            if not fams[f].has_key(record.id):
                fams[f][record.id] = list()
            fams[f][record.id].append(rec.id)
        c += 1
    return res

def writePairwiseSimilarities(fams, outDir, genomes, markers):
    import pdb; pdb.set_trace() 
    chrs = dict((k, sorted(v.keys())) for k, v in fams.items())
    gene2fam = dict(chain(*(tuple((g, f) for g in v.values()) for f, v in
        fams.items())))

    genes2pos = dict()
    for Gx in genomes:
        c = 1
        for chrx in chrs[Gx]:
            genes = map(lambda x: x[1], sorted(markers[chrx]))
            genes2pos.update(izip(map(lambda x: x.id, genes), xrange(c,
                c+len(genes))))
            c += len(genes)
    
    for Gx, Gy in combinations(genomes, 2):
        out = open(join(outDir, '%s_%s.sim'))
        for chrx in chrs[Gx]:
            for gx in markers[chrx]:
                f = gene2fam[gx]
                for chry in chrs[Gy]:
                    if fams[f].has_key(chry):
                        for gy in sorted(map(genes2pos.get, fams[f][chry])):
                            print >> out, '\t'.join((chrx, str(gx), chry,
                                str(gy), '*', '1'))


def writeMarkers(fastaFiles, outDir, segments):
    chr2genome = dict()
    fams = dict()
    partitioned_records = dict()
                 
    for f in fastaFiles:
        outName = join(outDir, basename(f))
        if abspath(outName) == abspath(f):
            print >> stderr, ('Designated output file %s is same as its ' + \
                    'correspondign input file.  Unable to store extracted' + \
                    'sequence data. Exiting.') %abspath(outName)
            exit(1)
        out = open(outName, 'w')
        for record in SeqIO.parse(open(f), 'fasta'):
            if chr2genome.has_key(record.id):
                print >> stderr, ('FATAL: record ID %s occurs more than ' + \
                        'once. Make sure that all record IDs are unique.') %(record.id)
                exit(1)
            chr2genome[record.id] = outName

            if segments.has_key(record.id):
                recs = paritionRecord(record, segments[record.id], fams)
                partitioned_records[record.id] = recs
                SeqIO.write(recs, out, 'fasta')
            else:
                print >> stderr, 'No segment found on fasta record %s' %record.id
        out.close()
    return chr2genome, fams, partitioned_records

if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-o', '--out_dir', default=DEFAULT_OUT_DIR, type=str,
            help='prefix of output filename')
    parser.add_argument('atoms_file', type=str, help='file containing atoms')
    parser.add_argument('fasta_file', type=str, nargs='+', 
            help='fasta files with original genomic sequences')

    args = parser.parse_args()

    segments = readSegments(open(args.atoms_file))
    chr2genome, fams, markers = writeMarkers(args.fasta_file, args.out_dir, segments)
    writePairwiseSimilarities(fams, args.out_dir, map(lambda x:
        basename(x).rsplit('.', 1)[0], args.fasta_file), markers)
