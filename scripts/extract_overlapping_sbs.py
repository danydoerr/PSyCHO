#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, exit, maxint, argv
from os.path import join, dirname, basename, isfile
from itertools import izip, chain
from bisect import insort
import logging
import json
import re

from pairwise_similarities import readDists, reverseDistMap, readGenomeMap, \
        PAT_POS, PAT_ID, PAT_CHR, GENOME_MAP_FILE, GM_ACTV_GNS_KEY


LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)
LOG_FILENAME = '%s' %basename(argv[0]).rsplit('.py', 1)[0]

def computeOverlaps(mlists, ciss):
    res = dict()
    for x in xrange(len(mlists)):
        queue = list()
        goss = mlists[x]
        #
        # XXX this whole procedure can be significanlty sped up by sorting start
        # positions before hand and then using a pointer within this list.
        #
        for chrx, start, end, i in goss:
            # only iterate until an interval is found that ends before new
            # interval start
            y = len(queue)-1
            while y >= 0 and queue[y][0] == chrx and queue[y][1] >= start:
                j = queue[y][3]
                if queue[y][1] < end and queue[y][2] < start and \
                        len(ciss[i][0]) == len(ciss[j][0]):
                    ii, jj = i < j and (i, j) or (j, i)
                    if not res.has_key((ii, jj)):
                        res[(ii, jj)] = set()
                    res[(ii, jj)].add(x)
                y -= 1
            
            insort(queue, (chrx, end, start, i))
    return res 


def parseMarkerList(marker_seq_list):
   
    res = list()
    for x in xrange(len(genomes)):
        mlist = list()
        queue = list()
        goss = marker_seq_list[x]
        for i in xrange(len(goss)):
            start = int(PAT_POS.match(goss[i][1]).group(1))
            end = int(PAT_POS.match(goss[i][-2]).group(2))
            chrx = PAT_CHR.match(goss[i][1]).group(1)
            mlist.append((chrx, start, end, i))
        
        mlist.sort()
        res.append(mlist) 

    return res



if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('sb_hierarchy', type=str, 
            help='PSyCHO JSON output file')

    args = parser.parse_args()


    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.ERROR)
    ch.setFormatter(logging.Formatter('!! %(message)s'))
    cf = logging.FileHandler('%s.log' %LOG_FILENAME, mode='w', delay=True)
    cf.setLevel(logging.DEBUG)
    cf.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(cf)
    LOG.addHandler(ch)

    #
    # load hiearchy data
    #
    LOG.info('reading SB hierarchy from file %s' %args.sb_hierarchy)
    jsDict = json.load(open(args.sb_hierarchy))
    ref = jsDict['ref_id']
    genomes = jsDict['genome_names']
    msl = jsDict['marker_seq_list']
    mlists = parseMarkerList(msl)
    ciss = jsDict['raw_sbfs']
    

    LOG.info('identifying overlaps between syntenic blocks (SBs)')
    overlapping_sbs = computeOverlaps(mlists, ciss)

    LOG.info('output synteny blocks...')
    out = stdout 
    for (i, j), o in overlapping_sbs.items():
        for ig, g in enumerate(genomes):
            start_si= PAT_POS.match(msl[ig][i][1]).group(1)
            end_si = PAT_POS.match(msl[ig][i][-2]).group(2)
            start_sj = PAT_POS.match(msl[ig][j][1]).group(1)
            end_sj = PAT_POS.match(msl[ig][j][-2]).group(2)

            if ig:
                out.write('\t')
            out.write('%s:%s-%s, %s-%s' %(g, start_si, end_si, start_sj,
                end_sj))
        out.write('\n')

    LOG.info('done!')
    
