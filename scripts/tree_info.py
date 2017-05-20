#!/usr/bin/env python

from sys import stdout, stderr, exit
from os.path import join, dirname, isfile
from itertools import izip, chain
from optparse import OptionParser
import json
import re

from psycho import dict2hierarchy 
from pairwise_similarities import PAT_POS, PAT_CHR

MIN_SIZE = 10


def compute_coverage(root, marker_seq_list, recovered_markers, ref):
    covered_markers = list(set() for _ in marker_seq_list)

    queue = [root]
    while queue:
        u = queue.pop()
        if u.links:
            for x, start, end in chain(((ref, u.intt[0], u.intt[1]),), u.links):
                gos = marker_seq_list[x][u.id]
                for p in xrange(start, end+1):
                    if gos[p] not in recovered_markers[x]:
                        m = PAT_POS.match(gos[p])
                        if m:
                            covered_markers[x].add(tuple(map(int,
                                m.group(1,2))))
        else:
            for child in u.children:
                if child.children:
                    queue.append(child)

    res = list()
    for x in covered_markers:
        res.append((sum(map(lambda z: z[1]-z[0]+1, x)), len(x)))
    return res

if __name__ == '__main__':
    usage = '%prog [options] <PSYCHO JSON OUTPUT>'
    parser = OptionParser(usage=usage)

    (options, args) = parser.parse_args()

    # check if input file is given, if so, call function with data
    if len(args) != 1:
        parser.print_help()
        exit(1)

    #
    # load hiearchy data
    #
    jsDict = json.load(open(args[0]))
    ref = jsDict['ref_id']
    root = dict2hierarchy(jsDict['sb_hierarchy'])
    genomes = jsDict['genome_names']
    marker_seq_list = jsDict['marker_seq_list']
    recovered_markers = jsDict['recovered_markers']

    queue = [(root, 1)]
    depths = list()
    syn_blocks = 0

    while queue:
        u, u_depth = queue.pop()
        if u.links:
            syn_blocks += 1
        for child in u.children:
            if len(child.children) > 1:
                queue.append((child, u_depth+1))
            else:
                depths.append(u_depth+1)

    print >> stdout, 'max tree depth: %s' %max(depths)
    print >> stdout, 'avg tree depth: %s' %(sum(depths)/float(len(depths)))
    print >> stdout, '# internal nodes w. links: %s' %syn_blocks
    print >> stdout, '# top-level blocks: %s' %sum(1 for x in
            xrange(len(marker_seq_list[ref])) if all(map(lambda y:
                len(marker_seq_list[y][x])-2 > MIN_SIZE,
                xrange(len(genomes)))))

    coverage = compute_coverage(root, marker_seq_list, recovered_markers, ref)
    print >> stdout, 'coverage (#markers):'
    for x in xrange(len(coverage)):
        print >> stdout, '\t%s:\t%s (%s)' %(genomes[x], coverage[x][0], coverage[x][1])


