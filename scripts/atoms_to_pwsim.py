#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from itertools import chain, combinations 
from sys import argv, stdout, stderr, exit
from os.path import basename, join, abspath, relpath
from Bio import SeqIO
from Bio.Seq import Seq, Alphabet
from string import maketrans
import logging
import csv

from pairwise_similarities import GM_FILE_KEY, GM_CHR_KEY, GM_ACTV_GNS_KEY, \
        GENOME_MAP_FILE, writeGenomeMap, PAT_CHR

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)
LOG_FILENAME = '%s' %basename(argv[0]).rsplit('.py', 1)[0]


TRANS_TABLE=maketrans('ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz',\
        'ANCNNNGNNNNNNNNNNNNTNNNNNMancnnngnnnnnnnnnnnntnnnnnn')
DEFAULT_OUT_DIR = '.'


def readSegments(data):
    res = dict()
    isHeader = True

    fam_members = dict()
    segs = list()
    for line in csv.reader(data, delimiter='\t'):
        if isHeader:
            isHeader = False
            continue
        if not res.has_key(line[0]):
            res[line[0]] = list()
        segs.append((line[0], int(line[4]), int(line[5]), line[2]))
        fam_members[line[2]] = fam_members.get(line[2], 0) + 1

    for chrx, start, stop, f in segs:
        if fam_members[f] < 2:
            continue
        if not res.has_key(chrx):
            res[chrx] = list()
        res[chrx].append((start, stop, f)) 

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
        res.append((start, rec)) 
        if not fams.has_key(f):
            fams[f] = list()
        fams[f].append((record.id, c))
        c += 1
    return map(lambda x: x[1], sorted(res))


def writePairwiseSimilarities(fams, genomes, gMap, outDir):

    chr2genome = dict()
    for Gx, gdata in gMap.items():
        for chrx in gdata[GM_CHR_KEY]:
            if chr2genome.has_key(chrx):
                LOG.fatal(('FASTA record ID %s occurs more than once in ' + \
                        'entire dataset. Make sure that all record IDs are ' + \
                        'unique.') %(chrx))
                exit(1)
            chr2genome[chrx] = Gx

    pwSims = dict(((Gx, Gy), list()) for Gx, Gy in chain(combinations(genomes,
        2), map(lambda x: (x,x), genomes)))

    for f, markers in fams.items():
        # because families are already filtered when reading the input
        # (function readSegments), it is guaranteed that combinations(. 2)
        # won't fail
        for mx, my in combinations(markers, 2):
            Gx = chr2genome[mx[0]]
            Gy = chr2genome[my[0]]
            if Gx == Gy:
                pwSims[(Gx, Gy)].append(tuple(sorted((mx, my))))
            elif pwSims.has_key((Gx, Gy)):
                pwSims[(Gx, Gy)].append((mx, my))
            else:
                pwSims[(Gy, Gx)].append((my, mx))

    for (Gx, Gy), pwData in pwSims.items():
        outFile = join(outDir, '%s_%s.sim' %(Gx, Gy))
        out = open(outFile, 'w')
        for mx, my in sorted(pwData):
            print >> out, '\t'.join((mx[0], str(mx[1]), my[0], str(my[1]), '*', '1'))
        out.close()


def writeMarkers(fastaFiles, outDir, segments):
    gMap = dict()
    gnames = list()
    fams = dict()
    
    for f in fastaFiles:
        outName = join(outDir, basename(f))
        if abspath(outName) == abspath(f):
            LOG.fatal(('Designated output file %s is same as its ' + \
                    'corresponding input file.  Unable to store extracted' + \
                    'sequence data. Exiting.') %abspath(outName))
            exit(1)
        out = open(outName, 'w')

        Gx = basename(f).rsplit('.', 1)[0]
        gnames.append(Gx)
        chrs = list()
        markers = list()
        for record in SeqIO.parse(open(f), 'fasta'):
            if segments.has_key(record.id):
                recs = paritionRecord(record, segments[record.id], fams)
                if recs:
                    markers.extend(recs)
                    SeqIO.write(recs, out, 'fasta')
                    chrs.append(record.id)
            else:
                LOG.info('No segment found on fasta record %s' %record.id)
        out.close()

        gMap[Gx] = dict()
        gMap[Gx][GM_FILE_KEY] = basename(outName)
        gMap[Gx][GM_CHR_KEY] = sorted(chrs)
        gMap[Gx][GM_ACTV_GNS_KEY] = markers
    
    return gnames, gMap, fams


if __name__ == '__main__':

    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-o', '--out_dir', default=DEFAULT_OUT_DIR, type=str,
            help='prefix of output filename')
    parser.add_argument('-s', '--sort_genome_names', action='store_true', 
            help='Rather than using the input argumnt order, Sort genomes ' + \
                    'names in alphabetical order prior to producing pairwise ' + \
                    'similarity file names.')
    parser.add_argument('atoms_file', type=str, help='file containing atoms')
    parser.add_argument('fasta_file', type=str, nargs='+', 
            help='fasta files with original genomic sequences')

    args = parser.parse_args()

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.WARNING)
    ch.setFormatter(logging.Formatter('!! %(message)s'))
    cf = logging.FileHandler('%s.log' %LOG_FILENAME, mode='w', delay=True)
    cf.setLevel(logging.DEBUG)
    cf.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(cf)
    LOG.addHandler(ch)

    segments = readSegments(open(args.atoms_file))
    genomes, gMap, fams = writeMarkers(args.fasta_file, args.out_dir,
            segments)
    if args.sort_genome_names:
        genomes = sorted(genomes)

    writeGenomeMap(gMap, genomes, open(join(args.out_dir, GENOME_MAP_FILE),
        'w'))

    writePairwiseSimilarities(fams, genomes, gMap, args.out_dir)

