#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import argv,stdout,stderr,exit
from os.path import basename, join, abspath
from Bio import SeqIO
from Bio.Seq import Seq, Alphabet
from string import maketrans
import csv

from pairwise_similarities import GM_FILE_KEY, GM_CHR_KEY, GM_ACTV_GNS_KEY, \
        GENOME_MAP_FILE, writeGenomeMap


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

def writePairwiseSimilarities(fams, chr2genome, markers, fastaFiles, outDir):
    import pdb; pdb.set_trace() 
    gene2fam = dict(chain(*(tuple((g, f) for g in v.values()) for f, v in
        fams.items())))
    gMap = dict()

    chrs = dict()
    for chrx, Gx in chr2genome.items():
        if not chrs.has_key(Gx):
            chrs[Gx] = list()
        chrs[Gx].append(chrx)
    for v in chrs.values():
        v.sort()

    genes2pos = dict()
    genomes = list()
    for f in fastaFiles:
        Gx = basename(f).rsplit('.', 1)[0]
        gMap[Gx] = dict()
        gMap[Gx][GM_FILE_KEY] = relpath(f, outDir or '.')
        gMap[Gx][GM_CHR_KEY] = chrs[Gx]
        active_markers = list()
        genomes.append(Gx)

        c = 1
        for chrx in chrs[Gx]:
            genes = map(lambda x: x[1], sorted(markers[chrx]))
            genes2pos.update(izip(map(lambda x: x.id, genes), xrange(c,
                c+len(genes))))
            c += len(genes)
            active_markers.extend(genes)
        gMap[Gx][GM_ACTV_GNS_KEY] = active_markers
    
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
    writeGenomeMap(gMap, genomes, open(join(outDir, GENOME_MAP_FILE)))


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
            chr2genome[record.id] = basename(f).rsplit('.', 1)[0]

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
    writePairwiseSimilarities(fams, chr2genome,  markers, args.fasta_file,args.out_dir)

