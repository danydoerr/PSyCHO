#!/usr/bin/env python

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter as ADHF
from sys import stdout, stderr, exit, maxint, argv
from os.path import join, dirname, basename, isfile
from itertools import izip, chain
from bisect import insort
import logging
import json
import re

import matplotlib; matplotlib.use('Agg')
from matplotlib import pylab as plt
from matplotlib.patches import Rectangle, Polygon

from psycho import dict2hierarchy 
from pairwise_similarities import readDists, reverseDistMap, readGenomeMap, \
        PAT_POS, PAT_ID, PAT_CHR, GENOME_MAP_FILE, GM_ACTV_GNS_KEY

DRAW_SEG_HEIGHT = 100
DRAW_SCALE = 0.1

DEFAULT_OUT_FILENAME_PREFIX = 'overlapping_SB_'

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
            # only iterate until an interval is found that ends before new interval start
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


def drawSegments(ax, gos, color, scale=1, offset=(0, 0)):

    xoff, yoff = offset
    
    for g in gos:
        start, end = map(int, PAT_POS.match(g).groups()[:2])
        r = Rectangle(((start+xoff)*scale, (yoff-DRAW_SEG_HEIGHT/2.)*scale),
                (end-start)*scale, DRAW_SEG_HEIGHT * scale, fill=True,
                edgecolor='none', facecolor=color)
        ax.add_patch(r)

def drawOverlappingSBFS(genomes, msl, dists, gMap, i, j, source, target, ax):

    start_si= int(PAT_POS.match(msl[source][i][1]).group(1))
    end_si = int(PAT_POS.match(msl[source][i][-2]).group(2))
    start_sj = int(PAT_POS.match(msl[source][j][1]).group(1))
    end_sj = int(PAT_POS.match(msl[source][j][-2]).group(2))

    start_ti= int(PAT_POS.match(msl[target][i][1]).group(1))
    end_ti = int(PAT_POS.match(msl[target][i][-2]).group(2))
    start_tj = int(PAT_POS.match(msl[target][j][1]).group(1))
    end_tj = int(PAT_POS.match(msl[target][j][-2]).group(2))


    
    chrs1 = PAT_CHR.match(msl[source][i][1]).group(1)
    chrs2 = PAT_CHR.match(msl[source][j][1]).group(1)
    chrt1 = PAT_CHR.match(msl[target][i][1]).group(1)
    chrt2 = PAT_CHR.match(msl[target][j][1]).group(1)


    if chrs1 == chrs2:
        OFFSET_SI = OFFSET_SJ = - min(start_si, start_sj) 
    else:
        gap = (end_si-start_si + end_sj-start_sj)/10.
        if chrt1 == chrt2 and start_tj < start_ti:
            OFFSET_SJ = -start_sj-gap
            OFFSET_SI = OFFSET_SJ + end_sj + 2*gap - start_si
        else:
            OFFSET_SI = -start_si-gap
            OFFSET_SJ = OFFSET_SI + end_si + 2*gap - start_sj
    
    if chrt1 == chrt2:
        OFFSET_TI = OFFSET_TJ = (min(end_si, end_sj)+max(start_si, start_sj))/2. \
                + min(OFFSET_SI, OFFSET_SJ) - min(end_ti, end_tj)
    else:
        if start_sj > start_si:
            OFFSET_TI = end_si + min(OFFSET_SI, OFFSET_SJ) - end_ti 
            OFFSET_TJ = OFFSET_TI + end_ti - start_tj 
            gap = end_tj+OFFSET_TJ/10.
        else:
            OFFSET_TJ = end_sj + min(OFFSET_SI, OFFSET_SJ) - end_tj 
            OFFSET_TI = OFFSET_TJ + end_tj - start_ti 
            gap = end_ti+OFFSET_TI/10.
        OFFSET_TI -= gap
        OFFSET_TJ += gap
            

    GAP_Y = 2000

    sgi = set(msl[source][i][1:-1])
    sgj = set(msl[source][j][1:-1])
    tgi = set(msl[target][i][1:-1])
    tgj = set(msl[target][j][1:-1])

    
    r = Rectangle(((start_si+OFFSET_SI)*DRAW_SCALE,
        ((GAP_Y-DRAW_SEG_HEIGHT)/2.)*DRAW_SCALE), (end_si-start_si)*DRAW_SCALE,
        DRAW_SEG_HEIGHT* DRAW_SCALE, fill=True, edgecolor='none', facecolor='m',
        alpha=0.2)
    ax.add_patch(r)
    drawSegments(ax, sgi.difference(sgj), 'm', scale=DRAW_SCALE,
            offset=(OFFSET_SI, GAP_Y/2.))

    r = Rectangle(((start_sj+OFFSET_SJ)*DRAW_SCALE,
        ((GAP_Y-DRAW_SEG_HEIGHT)/2.)*DRAW_SCALE), (end_sj-start_sj)*DRAW_SCALE,
        DRAW_SEG_HEIGHT* DRAW_SCALE, fill=True, edgecolor='none', facecolor='g',
        alpha=0.2)
    ax.add_patch(r)
    drawSegments(ax, sgj.difference(sgi), 'g', scale=DRAW_SCALE,
            offset=(OFFSET_SJ, GAP_Y/2.))
    drawSegments(ax, sgj.intersection(sgi), 'b', scale=DRAW_SCALE,
            offset=(OFFSET_SJ, GAP_Y/2.))

    r = Rectangle(((start_ti+OFFSET_TI) * DRAW_SCALE,
        ((-GAP_Y-DRAW_SEG_HEIGHT)/2.) * DRAW_SCALE),
        (end_ti-start_ti) * DRAW_SCALE, DRAW_SEG_HEIGHT * DRAW_SCALE, fill=True,
        edgecolor='none', facecolor='m', alpha=0.2)
    ax.add_patch(r)
    drawSegments(ax, tgi.difference(tgj), 'm', scale=DRAW_SCALE,
            offset=(OFFSET_TI, -GAP_Y/2.))

    r = Rectangle(((start_tj+OFFSET_TJ)*DRAW_SCALE,
        ((-GAP_Y-DRAW_SEG_HEIGHT)/2.)*DRAW_SCALE), (end_tj-start_tj)*DRAW_SCALE,
        DRAW_SEG_HEIGHT* DRAW_SCALE, fill=True, edgecolor='none', facecolor='g',
        alpha=0.2)
    ax.add_patch(r)
    drawSegments(ax, tgj.difference(tgi), 'g', scale=DRAW_SCALE,
            offset=(OFFSET_TJ, -GAP_Y/2.))
    drawSegments(ax, tgj.intersection(tgi), 'b', scale=DRAW_SCALE,
            offset=(OFFSET_TI, -GAP_Y/2.))

    ax.set_xlim(min((OFFSET_SI+start_si, OFFSET_SJ+start_sj,
        OFFSET_TI+start_ti, OFFSET_TJ+start_tj))*DRAW_SCALE,
        max((end_si+OFFSET_SI, end_sj+OFFSET_SJ, end_ti+OFFSET_TI,
            end_tj+OFFSET_TJ)) * DRAW_SCALE)
    ax.set_ylim((-GAP_Y/2. - DRAW_SEG_HEIGHT) * DRAW_SCALE, (GAP_Y/2. +
        DRAW_SEG_HEIGHT) * DRAW_SCALE)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)

    plt.tick_params(top=False, bottom=True, right=False, labelbottom=True)
    ax.set_yticks([-GAP_Y/2.*DRAW_SCALE, GAP_Y/2.*DRAW_SCALE])
    ax.set_yticklabels([genomes[target], genomes[source]], fontsize=16)
    

    pwDist = dists[(genomes[source], genomes[target])]
    for gs in chain(sgi, sgj):
        sstart, send = map(int, PAT_POS.match(gs).group(1,2))
        chrx, gid = PAT_ID.match(gs).group(1,2)
        gid = int(gid)
        csi = gs in sgi
        csj = gs in sgj

        OFFSET_S = OFFSET_SI
        if csj:
           OFFSET_S = OFFSET_SJ 

        id2pos = dict()

        __f__ = lambda x: (x[0], int(x[1]))
        for gi, g in enumerate(genomes):
            id2pos[gi] = dict((__f__(PAT_ID.match(x).group(1,2)), i)  for i, x \
                    in enumerate(gMap[g][GM_ACTV_GNS_KEY]))
        for gt_id, (_, w) in pwDist[(chrx, gid)].items():
            gt = gMap[genomes[target]][GM_ACTV_GNS_KEY][id2pos[target][gt_id]]
            tstart, tend = map(int, PAT_POS.match(gt).groups()[:2])
            cti = gt in tgi
            ctj = gt in tgj

            OFFSET_T = OFFSET_TI
            if ctj:
               OFFSET_T = OFFSET_TJ

            if (csi and csj) or (cti and ctj):
                OFFSET_S = OFFSET_SI

                xy = [((sstart+OFFSET_S)*DRAW_SCALE, \
                        (GAP_Y-DRAW_SEG_HEIGHT)/2.*DRAW_SCALE), \
                        ((send+OFFSET_S)*DRAW_SCALE, \
                        (GAP_Y-DRAW_SEG_HEIGHT)/2.*DRAW_SCALE), \
                        ((tend+OFFSET_T)*DRAW_SCALE, \
                        (-GAP_Y+DRAW_SEG_HEIGHT)/2.*DRAW_SCALE), \
                        ((tstart+OFFSET_T)*DRAW_SCALE,\
                        (-GAP_Y+DRAW_SEG_HEIGHT)/2.*DRAW_SCALE)]
                p = Polygon(plt.array(xy), closed=True, fill=True, \
                        edgecolor='none', facecolor='b', alpha=w)
                ax.add_patch(p)
            elif csi and cti:
                xy = [((sstart+OFFSET_S)*DRAW_SCALE, \
                        (GAP_Y-DRAW_SEG_HEIGHT)/2.*DRAW_SCALE), \
                        ((send+OFFSET_S)*DRAW_SCALE, \
                        (GAP_Y-DRAW_SEG_HEIGHT)/2.*DRAW_SCALE), \
                        ((tend+OFFSET_T)*DRAW_SCALE, \
                        (-GAP_Y+DRAW_SEG_HEIGHT)/2.*DRAW_SCALE), \
                        ((tstart+OFFSET_T)*DRAW_SCALE,\
                        (-GAP_Y+DRAW_SEG_HEIGHT)/2.*DRAW_SCALE)]
                p = Polygon(plt.array(xy), closed=True, fill=True, \
                        edgecolor='none', facecolor='m', alpha=w)
                ax.add_patch(p)
            elif csj and ctj:
                xy = [((sstart+OFFSET_S)*DRAW_SCALE, \
                        (GAP_Y-DRAW_SEG_HEIGHT)/2.*DRAW_SCALE), \
                        ((send+OFFSET_S)*DRAW_SCALE, \
                        (GAP_Y-DRAW_SEG_HEIGHT)/2.*DRAW_SCALE), \
                        ((tend+OFFSET_T)*DRAW_SCALE, \
                        (-GAP_Y+DRAW_SEG_HEIGHT)/2.*DRAW_SCALE), \
                        ((tstart+OFFSET_T)*DRAW_SCALE,\
                        (-GAP_Y+DRAW_SEG_HEIGHT)/2.*DRAW_SCALE)]
                p = Polygon(plt.array(xy), closed=True, fill=True, \
                        edgecolor='none', facecolor='g', alpha=w)
                ax.add_patch(p)


if __name__ == '__main__':
    parser = ArgumentParser(formatter_class=ADHF)
    parser.add_argument('-o', '--out_filename_prefix',
            default=DEFAULT_OUT_FILENAME_PREFIX, type=str, 
            help='prefix of output filename')
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
    marker_seq_list = jsDict['marker_seq_list']
    mlists = parseMarkerList(marker_seq_list)
    ciss = jsDict['raw_sbfs']
    pwsim_files = jsDict['pwsim_files']
    gMap = readGenomeMap(open(join(dirname(pwsim_files[0]), GENOME_MAP_FILE)))
    
    LOG.info('reading pairwise similarities... ')
    dists = dict()
    for f in pwsim_files:
        LOG.info('  %s' %f)
        _ , dist = readDists(open(f))
        gname1, gname2  = basename(f).split('.', 1)[0].split('_', 1)
        dists[(gname1, gname2)] = dist
        dists[(gname2, gname1)] = reverseDistMap(dist)

    LOG.info('identifying overlaps between syntenic blocks (SBs)')
    overlapping_sbs = computeOverlaps(mlists, ciss)

    LOG.info('drawing figures...')
    for (i, j), o in overlapping_sbs.items():
        out = open('%s%s-%s.pdf' %(args.out_filename_prefix, i, j), 'w')
        LOG.info(('  plotting overlap between SBs %s and %s and storing ' + \
                'figure in %s') %(i, j, out.name))
        source = ref
        f, axs = plt.subplots(len(genomes)-1, 1, sharey=False, sharex=False)

        ax_it = iter(axs) 
        for target in xrange(len(genomes)):
            if source == target:
                continue

            ax = next(ax_it)
            drawOverlappingSBFS(genomes, marker_seq_list, dists, gMap, i, j,
                    source, target, ax)

        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        f.savefig(out, format='pdf')
        out.close()
        #raw_input('Press Enter to continue...')

    LOG.info('done!')
    
