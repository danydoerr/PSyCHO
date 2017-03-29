#!/usr/bin/env python

from sys import stdout,stderr,exit, argv
from optparse import OptionParser
from os.path import basename, join
from itertools import combinations
from Bio import SeqIO
from Bio import pairwise2
from bisect import bisect
import logging
import csv

from pairwise_similarities import PAT_CHR

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
            cur_len += len(record)
            if ident not in seqRecords:
                seqRecords[ident] = [f, record]
                chr_locations[ident] = list() 
            else:
                seqRecords[ident][1] += record

            chr_id = record.id
            m = PAT_CHR.match(record.id)
            if m:
                chr_id = m.group(1)

            chr_locations[ident].append((cur_len, chr_id))
    
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

def parseMauveBackbone(mauveFile, minLength, quorum):

    mauveSegments = [list() for _ in seqData]
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
           
            orient = int(line[i*2]) >= 0 and '+' or '-'
            start, end = abs(int(line[i*2])), abs(int(line[i*2+1]))

            if end-start <= 0:
                LOG.warning(('sequence length of segment seq%s[%s:%s] is %s in' + \
                        ' line: \n\t\t%s') %(i, line[i*2], line[i*2+1],
                            end-start+1, '\t'.join(line)))

            if end-start+1 < minLength:
                continue

            g1i = (outCount, start, end, orient)
            cur_orth.append((i, g1i))
        
        if len(cur_orth) >= quorum:
            for i, g1i in cur_orth:
                mauveSegments[i].append(g1i)
            orthologies.extend(combinations(map(lambda x: x[1], cur_orth), 2))
            outCount += 1

    for i in xrange(len(mauveSegments)):
        mauveSegments[i].sort(key=lambda x: (x[1], x[2]), cmp=lambda y,z:
                cmp(y[0], z[0]) or -cmp(y[1], z[1]))
    
    LOG.info('Identified %s SBFs larger length %s' %(outCount, minLength))
    return mauveSegments, orthologies


def removeDuplAndStrip(mauveSegments):
    res = [mauveSegments[0]]
    for i in xrange(1, len(mauveSegments)):
        # check if disjoint
        if res[-1][2] < mauveSegments[i][1]:
            res.append(mauveSegments[i])
        # check if not contained
        elif res[-1][2] < mauveSegments[i][2]:
            if res[-1][2]-res[-1][1] > mauveSegments[i][2]-mauveSegments[i][1]:
                res[-1] = (res[-1][0], res[-1][1], mauveSegments[i][1]-1,
                        res[-1][3])
                res.append(mauveSegments[i])
            else:
                res.append((mauveSegments[i][0], res[-1][2]+1,
                    mauveSegments[i][2], mauveSegments[i][3]))
        # skip (i.e. do not append) otherwise
    return res


def writeSegments(seqData, chr_locations, mauveSegments, ident, minLength, out, stripNs=False):

    c = i = chr_start = 0
    chr_end, chr_name = chr_locations[c]
    while i < len(mauveSegments):
        iid, start, end, orient = mauveSegments[i]
        gid = '%s_%s' %(ident, iid)

        while c < len(chr_locations)-1 and chr_end <= start-1:
            chr_start = chr_locations[c][0]
            c += 1
            chr_end, chr_name = chr_locations[c]
       
        if end > chr_end:
            if end-chr_end >= minLength and c < len(chr_locations)-1:
                mauveSegments[i] = ('%s\'' %iid, chr_end+1, end, orient)
                i -= 1

            end = chr_end


        while stripNs and end-start+1 >= minLength and seqData[1][start-1].upper() == 'N':
            start += 1
        while stripNs and end-start+1 >= minLength and seqData[1][end-1].upper() == 'N':
            end -= 1
        if end-start+1 >= minLength:
            seqRecord = seqData[1][start-1:end]
            if orient == '-':
                seqRecord = seqRecord.reverse_complement()
            seqRecord.id = '%s|%s:%s|chromosome|%s|strand|%s' %(gid,
                    start-chr_start, end-chr_start, chr_name, orient)
            seqRecord.description = ''
            SeqIO.write(seqRecord, out, 'fasta')
        i += 1


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
    
    cf = logging.FileHandler(LOG_FILENAME, mode='w')
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
    mauveSegments, _ = parseMauveBackbone(mauveFile, options.minLength,
            options.genomeQuorum)

    for x in xrange(len(mauveSegments)):
        ident = basename(seqData[x][0]).rsplit('.', 1)[0]
        out = open('%s.gos' %join(options.outDir, ident), 'w')
        LOG.info('Writing %s' %out.name)
        segments = removeDuplAndStrip(mauveSegments[x])
        writeSegments(seqData[x], chr_locations[x], segments, ident,
                options.minLength, out, True)
        out.close()

    LOG.info('finished')

