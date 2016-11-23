#!/usr/bin/env python

from sys import stdout,stderr,exit
from optparse import OptionParser
from os.path import basename, join
from Bio import SeqIO
from Bio import pairwise2
import logging
import csv

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)

MIN_LENGTH_DEFAULT = 0

def readFastaFiles(seqFiles, mauveFile):
    seqRecords = dict()

    for f in seqFiles:
        # assume that fasta file contains a single entry
        ident= basename(f).rsplit('.', 1)[0]

        LOG.info('reading %s with identifier %s' %(f, ident))
        for record in SeqIO.parse(open(f), 'fasta'):
            if ident not in seqRecords:
                seqRecords[ident] = [f, record]
            else:
                seqRecords[ident][1] += record
    
    # put list into order indicated by mauve backbone filename
    seqList = list()

    for ident in basename(mauveFile).rsplit('.',1)[0].split('_'):
        if ident not in seqRecords:
            LOG.error(('No corresponding fasta file found for identifier %s'
                + ' of backbone %s. Exiting.') %(ident, basename(mauveFile)))
            exit(1)

        seqList.append(tuple(seqRecords[ident]))

    return seqList

def processMauveBackbone(seqData, mauveFile, minLength, unkownSeqPercent,
        outDir):

    isHeader = True

    # initialize output handles
    outFiles = [open('%s.gos' %join(outDir, basename(i).rsplit('.',
        1)[0], 'w')) for i, _ in seqData]

    outCount = 0
    for line in csv.reader(open(mauveFile), delimiter='\t'):
        if isHeader:
            isHeader = False
            continue

        hasOut = 0

        sequences = list()
        for i in range(len(line)/2):
            if line[i*2] == '0':
                continue
            
            ident = basename(seqData[i][0]).rsplit('.', 1)[0]
            orient = int(line[i*2]) >= 0 and '+' or '-'
            start, end = abs(int(line[i*2]))-1, abs(int(line[i*2+1]))

            if end-start <= 0:
                LOG.warning(('sequence length of segment %s[%s:%s] is %s in' + \
                        ' line: \n\t\t%s') %(ident, line[i*2], line[i*2+1],
                            end-start-1, '\t'.join(line)))

            if end-start-1 < minLength:
                continue

            hasOut = 1
            seqRecord = seqData[i][1][start:end]
            if orient == '-':
                seqRecord = seqRecord.reverse_complement()

            seq = seqRecord.seq.lower()

            if unkownSeqPercent > 0 and \
                    len(seq)-(seq.count('a') + seq.count('g') + \
                    seq.count('c') + seq.count('t')) > unkownSeqPercent * len(seq):
                LOG.warning(('sequence %s of segment %s[%s:%s] contains ' + \
                        'more unkown characters than ACGT\'s. Skipping ')%(len(seq) > 50 
                            and (seq[:50] + '...') or seq, ident, line[i*2],
                            line[i*2+1]))
                continue
            seqRecord.id = '%s_%s_%s|strand|%s' %(ident, start+1, end, orient)
            seqRecord.description = ''
            SeqIO.write(seqRecord, outFiles[i], 'fasta')

            sequences.append(seq)

#        if len(sequences) > 5 and len(sequences[0]) + len(sequences[5]) < 600:
#            aln_score = pairwise2.align.globalxx(sequences[0], sequences[5],
#                    score_only=True)
#            print 'aln score %s between seqs of length %s and %s' %(aln_score,
#                    len(sequences[0]), len(sequences[5]))


        outCount += hasOut

    LOG.info('wrote %s SBFs larger length %s to output files' %(outCount, minLength))



if __name__ == '__main__':

    usage = 'usage: %prog [options] <MAUVE BACKBONE> <FASTA FILE 1> ... <FASTA FILE N>'
    parser = OptionParser(usage=usage)

    parser.add_option('-X', '--unkown_seq_percentage', dest='unkownSeqPercent',
            help='Maximum allowed percentage of unknown (XXX\'ed) positions ' \
                    'marker to be included in output. [default=%default]',
                    type=float, default=0, metavar='[0, 1]')

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
    
    cf = logging.FileHandler('%s.log' %(basename(args[0]).rsplit('.', 1)[0]), mode='w')
    cf.setLevel(logging.INFO)
    cf.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t++ %(message)s'))

    LOG.addHandler(cf)
    LOG.addHandler(ch)

    #
    # main 
    #


    LOG.info(('start partitioning fasta sequences according to mauve ' \
            + 'backbone data using minimum segment length %s') %options.minLength)
    seqData = readFastaFiles(seqFiles, mauveFile)
    LOG.info('processing backbone file')
    processMauveBackbone(seqData, mauveFile, options.minLength,
            options.unkownSeqPercent, options.outDir)
    LOG.info('finished')

