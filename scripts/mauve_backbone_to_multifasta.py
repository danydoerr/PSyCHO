#!/usr/bin/env python

from sys import stdout,stderr,exit
from optparse import OptionParser
from os.path import basename, join
from itertools import combinations
from Bio import SeqIO
from Bio import pairwise2
from bisect import bisect
import logging
import csv

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)
LOG_FILENAME = '%s.log' %basename(argv[0]).rsplit('.py', 1)[0]

MIN_LENGTH_DEFAULT = 0
QUORUM_DEFAULT = 1

def readFastaFiles(seqFiles, mauveFile):
    seqRecords = dict()
    chr_locations = dict()

    for f in seqFiles:
        # assume that fasta file contains a single entry
        ident= basename(f).rsplit('.', 1)[0]

        LOG.info('reading %s with identifier %s' %(f, ident))
        cur_len = 0
        for record in SeqIO.parse(open(f), 'fasta'):
            if ident not in seqRecords:
                seqRecords[ident] = [f, record]
                chr_locations[ident] = list() 
            else:
                seqRecords[ident][1] += record
            chr_locations[ident].append((cur_len, record.id))
            cur_len += len(record)
    
    # put list into order indicated by mauve backbone filename
    seqList = list()
    chr_loc_list  = list()

    for ident in basename(mauveFile).rsplit('.',1)[0].split('_'):
        if ident not in seqRecords:
            LOG.error(('No corresponding fasta file found for identifier %s'
                + ' of backbone %s. Exiting.') %(ident, basename(mauveFile)))
            exit(1)

        chr_loc_list.append(tuple(chr_locations[ident]))
        seqList.append(tuple(seqRecords[ident]))

    return seqList, chr_loc_list

def parseMauveBackbone(seqData, chr_locations, mauveFile, minLength, quorum):

    genomes = [list() for _ in seqData]
    orthologies = list()

    isHeader = True
    outCount = 0

    for line in csv.reader(open(mauveFile), delimiter='\t'):
        if isHeader:
            isHeader = False
            continue

        cur_orth = list()
        seqRecords = [None] * (len(line)/2)
        for i in range(len(line)/2):
            if line[i*2] == '0':
                continue
            
            ident = basename(seqData[i][0]).rsplit('.', 1)[0]
            orient = int(line[i*2]) >= 0 and '+' or '-'
            start, end = abs(int(line[i*2]))-1, abs(int(line[i*2+1]))

            if end-start <= 0:
                LOG.warning(('sequence length of segment %s[%s:%s] is %s in' + \
                        ' line: \n\t\t%s') %(genomeNames[i], line[i*2], line[i*2+1],
                            end-start-1, '\t'.join(line)))

            if end-start-1 < minLength:
                continue

            seqRecord = seqData[i][1][start:end]
            if orient == '-':
                seqRecord = seqRecord.reverse_complement()

            x = bisect(chr_locations[i], (start, chr(255)))
            chr_start, chr_name = chr_locations[i][x-1]
            gid = '%s_%s' %(ident, outCount)
            g1i = (gid, start+1-chr_start, end-chr_start, chr_name, orient)
            cur_orth.append((i, g1i))
            seqRecord.id = '%s|%s:%s|chromosome|%s|strand|%s' %g1i
            seqRecord.description = ''
            seqRecords[i] = seqRecord
        
        if len(cur_orth) >= quorum:
            for i, g1i in cur_orth:
                genomes[i].append((g1i, seqRecords[i]))
            orthologies.extend(combinations(map(lambda x: x[1], cur_orth), 2))
            outCount += 1

    for i in xrange(len(genomes)):
        genomes[i].sort(key=lambda x: (x[0][3], x[0][1], x[0][2]))
        genomes[i] = map(lambda x: x[1], genomes[i])
    
    LOG.info('Identified %s SBFs larger length %s' %(outCount, minLength))
    return genomes, orthologies


if __name__ == '__main__':

    usage = 'usage: %prog [options] <MAUVE BACKBONE> <FASTA FILE 1> ... <FASTA FILE N>'
    parser = OptionParser(usage=usage)

    parser.add_option('-q', '--quorum', dest='genomeQuorum', help='Minimum' + \
            ' number of genomes that must be covered by each marker ' + \
            '[default=%default]', type=int, default=QUORUM_DEFAULT,
            metavar='INT')

    parser.add_option('-l', '--minlength', dest='minLength',
            help='Miminimum length a segment must have to be included in ' \
                    + 'the output. [default=%default]',
                    type=int, default=MIN_LENGTH_DEFAULT, metavar='INT')

    parser.add_option('-o', '--output_dir', dest='outDir', type=str,
            default='.', help='Output directory [default=%default]')

    (options, args) = parser.parse_args()

    if len(args) < 3:
        parser.print_help()
        exit(1)

    mauveFile = args[0]
    seqFiles = args[1:]

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.ERROR)
    ch.setFormatter(logging.Formatter('!! %(message)s'))
    
    cf = logging.FileHandler('%s.log' %(LOG_FILENAME, mode='w'))
    cf.setLevel(logging.INFO)
    cf.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t++ %(message)s'))

    LOG.addHandler(cf)
    LOG.addHandler(ch)

    #
    # main 
    #


    LOG.info(('start partitioning fasta sequences according to mauve ' \
            + 'backbone data using minimum segment length %s') %options.minLength)
    seqData, chr_locations = readFastaFiles(seqFiles, mauveFile)
    LOG.info('processing backbone file')
    genomes, _ = parseMauveBackbone(seqData, chr_locations, mauveFile,
            options.minLength, options.genomeQuorum)

    # initialize output handles
    for i in xrange(len(seqData)):
        out = open('%s.gos' %join(options.outDir, basename(seqData[i][0]).rsplit('.',
            1)[0]), 'w')
        SeqIO.write(genomes[i], out, 'fasta')

    LOG.info('finished')

