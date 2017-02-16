#!/usr/bin/env python

from sys import stdout, stderr, exit
from os.path import join, dirname, isfile
from itertools import izip, chain
from optparse import OptionParser
import shelve
import re

from pairwise_similarities import readGenomeMap, GM_ACTV_GNS_KEY, \
        GENOME_MAP_FILE
from psycho import node

PAT_CHR = re.compile('.*\|chromosome\|([^\|]+)(\|.*|$)')
PAT_POS = re.compile('.*\|(\d+):(\d+)(\|.*|$)')

MIN_SIZE = 10


def compute_coverage(root, gene_orders, g2pos, marker_loc, ref):
    covered_markers = list(set() for _ in gene_orders[0])

    queue = [(root, None)]
    while queue:
        u, u_id = queue.pop()
        # remember found u_id in top-down traversal 
        if u.id != None:
            u_id = u.id
        if u_id != None and u.links:
            gos = gene_orders[u_id]
            gos2pos = g2pos[u_id]
            for x, gx1, gx2 in chain(((ref, u.intt[0], u.intt[1]),), u.links):
                covered_markers[x].update(gos[x][gos2pos[x][gx1]:gos2pos[x][gx2]+1])
        else:
            for child in u.children:
                if child.children:
                    queue.append((child, u_id))
    res = list()
    for x in xrange(len(covered_markers)):
        res.append(sum(map(lambda z: marker_loc[x][z][1] -
            marker_loc[x][z][0]+1, covered_markers[x])))
    return res

if __name__ == '__main__':
    usage = '%prog [options] <SHELVE>'
    parser = OptionParser(usage=usage)

    (options, args) = parser.parse_args()

    # check if input file is given, if so, call function with data
    if len(args) != 1:
        parser.print_help()
        exit(1)


    shObj = shelve.open(args[0], flag='r', protocol=-1)
    ref = shObj['ref']

    #
    # load inclusion tree
    #
    root = shObj['inclusion_tree']
    genomes = shObj['genomes']
    gene_orders = shObj['gene_orders']
    
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
    print >> stdout, '# top-level blocks: %s' %sum(1 for x in gene_orders if
            all(map(lambda y: len(y)-2 > MIN_SIZE, x)))

    #
    # read marker positions from marker fasta files
    #
    gMapF = join(dirname(shObj['orig_pw_dists'][0]), GENOME_MAP_FILE)
    if not isfile(gMapF):
        print >> stderr, ('ERROR: unable to find genome map at assumed ' + \
                'location %s. Exiting') %gMapF
        exit(1)

    gMap = readGenomeMap(open(gMapF)) 
    marker_loc = list()
    for Gx in genomes:
        m_dict = dict()
        for x in xrange(len(gMap[Gx][GM_ACTV_GNS_KEY])):
            gx1 = gMap[Gx][GM_ACTV_GNS_KEY][x]
            m_dict[(PAT_CHR.match(gx1).group(1), x+1)] = map(int, \
                    PAT_POS.match(gx1).group(1,2))
        marker_loc.append(m_dict)
           
    g2pos = [[dict(izip(go, xrange(len(go)))) for go in gos] for gos in \
            gene_orders]

    coverage = compute_coverage(root, gene_orders, g2pos, marker_loc, shObj['ref'])
    print >> stdout, 'coverage:'
    for x in xrange(len(coverage)):
        print >> stdout, '\t%s:\t%s' %(genomes[x], coverage[x])


