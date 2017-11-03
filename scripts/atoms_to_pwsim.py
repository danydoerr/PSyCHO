#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from itertools import chain, combinations, izip, imap
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
DEFAULT_DUP_THRESHOLD = 11


def readChr2GenomeMapping(*files):
    res = dict()
    for f in files:
        gname = basename(f).rsplit('.', 1)[0]
        for line in open(f):
            if line.startswith('>'):
                chrId = line.split()[0][1:]
                res[chrId] = gname
    return res

def readSegments(data, chr2gMap, excludeDups=DEFAULT_DUP_THRESHOLD):

    genomes = sorted(set(chr2gMap.values()))
    g2pos = dict(izip(genomes, xrange(len(genomes))))
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
        segs.append((line[0], line[1], int(line[4]), int(line[5]), line[2]))
        if not fam_members.has_key(line[2]):
            fam_members[line[2]] = [0] * len(genomes)
        fam_members[line[2]][g2pos[chr2gMap[line[0]]]] += 1

    for chrx, sid, start, stop, f in segs:
        if sum(fam_members[f]) < 2 or any(imap(lambda x: x >= excludeDups,
                fam_members[f])):
            continue
        
        if not res.has_key(chrx):
            res[chrx] = list()
        res[chrx].append((start, stop, sid, f)) 

    # make sure segments are sorted
    for v in res.values():
        v.sort()
    return res


def partitionRecord(record, segments, fams):
    # assumes *sorted* segment list (function readSegments ensures that) 
    res = list()
    for start, stop, sid, f in segments:
        rec = record[start:stop]
        name = record.id.replace('|', '.')
        rec.id = '%s_%s|%s:%s|chromosome|%s' %(name, sid, start, stop-1, name)
        rec.description = ''
        rec.seq = Seq(str(rec.seq).translate(TRANS_TABLE),
                Alphabet.generic_dna)
        res.append(rec) 
        if not fams.has_key(f):
            fams[f] = list()
        fams[f].append((record.id, sid))
    return res


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
                pwSims[(Gx, Gy)].append((mx, my))
                pwSims[(Gx, Gy)].append((my, mx))
            elif pwSims.has_key((Gx, Gy)):
                pwSims[(Gx, Gy)].append((mx, my))
            else:
                pwSims[(Gy, Gx)].append((my, mx))

    for (Gx, Gy), pwData in pwSims.items():
        outFile = join(outDir, '%s_%s.sim' %(Gx, Gy))
        out = open(outFile, 'w')
        for mx, my in sorted(pwData, key=lambda x: (x[0][1], x[1][1])):
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
                recs = partitionRecord(record, segments[record.id], fams)
                if recs:
                    markers.extend(map(lambda x: x.id, recs))
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
    parser.add_argument('-e', '--exclude_duplicates',
            default=DEFAULT_DUP_THRESHOLD, type=int,
            help='exclude markers belonging to families that occur <x> ' + \
                    'or more times in any of the genomes.')
    parser.add_argument('-s', '--sort_genome_names', action='store_true', 
            help='rather than using the input argumnt order, sort genomes ' + \
                    'names in alphabetical order prior to producing pairwise ' + \
                    'similarity file names')
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

    chr2gMap = readChr2GenomeMapping(*args.fasta_file)
    segments = readSegments(open(args.atoms_file), chr2gMap,
            args.exclude_duplicates)

    genomes, gMap, fams = writeMarkers(args.fasta_file, args.out_dir,
            segments)
    if args.sort_genome_names:
        genomes = sorted(genomes)

    writeGenomeMap(gMap, genomes, open(join(args.out_dir, GENOME_MAP_FILE),
        'w'))

    writePairwiseSimilarities(fams, genomes, gMap, args.out_dir)

