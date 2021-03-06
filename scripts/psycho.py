#!/usr/bin/env python

from sys import stdout, stderr, argv, exit, maxint, setrecursionlimit
from multiprocessing import Pool, Queue, JoinableQueue, cpu_count
from itertools import izip, combinations, product, chain
from os.path import basename, abspath, relpath, dirname, join, isfile
from bisect import bisect, insort
from optparse import OptionParser
from functools import partial
from collections import deque
from tempfile import mkdtemp
from Queue import Empty
from math import sqrt
import networkx as nx
import logging
import json
import os



#
# enable import from parent/sibling modules
#

from pairwise_similarities import readDists, reverseDistMap, readGenomeMap, \
        DIRECTION_CRICK_STRAND, DIRECTION_WATSON_STRAND, \
        TELOMERE_END, TELOMERE_START, GENOME_MAP_FILE, GM_ACTV_GNS_KEY
from dll import doubly_linked_list as DLL, node as dnode

#setrecursionlimit(100000)

MAX_ITERATION = 10
CONTIG_BOUNDARY = (-1, maxint)
CONTIG_BOUNDARY_KEY = '__CONTIG_BOUNDARY__'

DEFAULT_REF = 'G1'
DEFAULT_COVERAGE_INV = 3

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)
LOG_FILENAME = '%s' %basename(argv[0]).rsplit('.py', 1)[0]

# inverse coverage (1/coverage)

class node(object):
    def __init__(self, left_b, right_b, type=None, parent=None, children=None,
            links=None,id=None):
        self.intt = (left_b, right_b)
        self.id = id
        self.type = type
        self.parent = parent
        self.links = links or list()
        if children:
            self.children = deque(children)
            self.children[0].parent = self
            self.children[0].next_sibling = None
            for i in xrange(1, len(children)):
                self.children[i].next_sibling = self.children[i-1]
                self.children[i].parent = None
        else:
            self.children = deque()


    def __str__(self):
        out = '' 
        if self.children:
            out += '(' + ' '.join(map(str, self.children)) + '):'
        out += '%s{%s-%s}' %(self.type or '?', self.intt[0],
                self.intt[1])
        if self.links:
            out+= '->[%s]' %(' '.join(map(lambda x: '%s:%s-%s' %tuple(x),
                self.links)))
        return out

def hierarchy2dict(root):

    res = list()

    queue = deque(((root, None), ))

    offset = 1
    while queue:
        v, mseq_id = queue.popleft()
        if mseq_id == None and v.id != None:
            mseq_id = v.id
        vDict = {'marker_seq_id': mseq_id,
                'ref_sb': v.intt,
                'linked_sbs': v.links,
                'children': range(offset, offset+len(v.children))}
        res.append(vDict)
        queue.extend(map(lambda x: (x, mseq_id), v.children))
        offset += len(v.children)

    return res


def dict2hierarchy(JSON_list):

    nodes = map(lambda x: node(x['ref_sb'][0], x['ref_sb'][1],
        id=x['marker_seq_id'], links=x['linked_sbs']), JSON_list)
    refd = set(xrange(len(JSON_list)))
    for x in xrange(len(JSON_list)):
        child_ids = JSON_list[x].get('children', ())
        refd.difference_update(child_ids)
        for i in child_ids:
            if nodes[x].children:
                nodes[i].next_sibling = nodes[x].children[-1]
            else:
                nodes[i].parent = nodes[x]
            nodes[x].children.append(nodes[i])

    if len(refd) > 1:
        raise Exception('Restored hierarchy has more than one root!')

    return nodes[next(iter(refd))]


def removeNonUniversalGenes(G, n):

    i = 0
    hasChanged = True
    LOG.info('start masking non-universal genes ...')
    while i < MAX_ITERATION and hasChanged:
        LOG.info('iteration %s, current number of active genes: %s' %(i, len(G)))
        hasChanged = False
        i += 1
        for v in G.nodes():
            if len(set(map(lambda x: x[0], \
                    G.neighbors_iter(v))).difference((v[0], ))) + 1 < n:
                G.remove_node(v)
                hasChanged = True
    LOG.info('Finished')
    return G


def establish_linear_genome_order(genes):

    chromosomes = dict()
    for (chr1, g1i) in genes:
        if not chromosomes.has_key(chr1):
            chromosomes[chr1] = [maxint, 0]
        if g1i < chromosomes[chr1][0]:
            chromosomes[chr1][0] = g1i
        if g1i > chromosomes[chr1][1]:
            chromosomes[chr1][1] = g1i

    chrnames = sorted(chromosomes.keys())
    
    g = list()
    for k in chrnames:
        # stopper for extend()
        g.append(CONTIG_BOUNDARY)
        g.extend([(k, i) for i in xrange(chromosomes[k][0], chromosomes[k][1]+1)])
    # stopper for extend()
    g.append(CONTIG_BOUNDARY)
    return  g


def constructTeamDS(gene_orders, G, ref, delta): 
    L = list()

    for y in xrange(len(gene_orders)):
        cL = DLL()
        c = - delta
        for gi in gene_orders[y]:
            if gi == CONTIG_BOUNDARY:
                c += delta
            elif G.has_node((y, gi)):
                # each element in L is of type (label, character)
                cL.add((c, gi))
            c += 1
        
        L.append(cL)

    # remove all edges from G that are not between reference and non-reference
    for (u, v) in G.edges():
        if (u[0] == ref) == (v[0] == ref):
            G.remove_edge(u, v) 
    return L, G


def smallMax(L, delta):

    head = [x.head for x in L]
    tail = [x.tail for x in L]

    crossed = [head[x] is tail[x] for x in xrange(len(L))]
    B = None 
    while B == None and not all(crossed):
        x = 0
        while B == None and x < len(L):
            if not crossed[x]:
                if head[x].next.data[0]-head[x].data[0] > delta:
                    LOG.debug(('identified (head-) split between %s and %s in ' + \
                            'genome %s') %(head[x].data, head[x].next.data, x))
                    B = (x, L[x].head, head[x])
                else:
                    head[x] = head[x].next
                    crossed[x] |= head[x] is tail[x]
                if tail[x].data[0]-tail[x].prev.data[0] > delta:
                    LOG.debug(('identified (tail-) split between %s and %s in ' + \
                            'genome %s') %(tail[x].prev.data, tail[x].data, x))
                    B = (x, tail[x], L[x].tail)
                else:
                    tail[x] = tail[x].prev
                    crossed[x] |= head[x] is tail[x]
            x += 1

    if B == None:
        B = (0, L[0].head, L[0].tail)
    return B


def division(L, G, B, ref):

    ids = set(x for x in xrange(len(L)) if x != ref)
    x, start, end = B

    # extracted data structures of L and N 
    LX = [None] * len(L)
    GX = nx.Graph()

    Bset = set()
    
    # construct \Sigma(X1)),  where X1 is _always_ a subsequence of the
    # reference genome
    SigmaX1 = set()
    if x != ref:
        cur = start
        SigmaX1.update(map(lambda x: x[1], G.neighbors((x, cur.data[1]))))
        Bset.add(cur.data[1])
        while cur != end:
            cur = cur.next
            SigmaX1.update(map(lambda x: x[1], G.neighbors((x, cur.data[1]))))
            Bset.add(cur.data[1])
    else:
        # construct alphabet of X1 
        cur = start
        SigmaX1.add(cur.data[1])
        while cur != end:
            cur = cur.next
            SigmaX1.add(cur.data[1])


    # extract X1 and X2 outgoing from reference sequence
    X1 = DLL()
    Y1 = L[ref]
    cur = Y1.head
    while cur != None:
        nxt = cur.next

        if cur.data[1] in SigmaX1:
            v = (ref, cur.data[1])
            cids = [u for u, _ in G.neighbors(v)]
            if not ids.issubset(cids):
                LOG.debug(('remove node %s because it is not connected to ' + \
                        'all genomes') %str(v))
                Y1.pop(cur)
                G.remove_node(v)
            elif x != ref and cids.count(x) > 1 and not Bset.issuperset(z[1] \
                    for z in G.neighbors(v) if z[0] == x):
                X1.add(cur.data)
                for u in G.neighbors(v):
                    if u[0] != x:
                        GX.add_edge(v, u)
                    elif u[1] in Bset:
                        # perform split between two gene teams that overlap in
                        # duplicates
                        GX.add_edge(v, u)
                        G.remove_edge(v, u)
            else:
                X1.append_right(Y1.pop(cur))
                GX.add_edges_from((v, u) for u in G.neighbors(v))
                G.remove_node(v)
        cur = nxt
    LX[ref] = X1

    for v in nx.isolates(G):
        G.remove_node(v)
    for v in nx.isolates(GX):
        GX.remove_node(v)


    for y in ids:
        X2 = DLL()
        Y2 = L[y]
        cur = Y2.head
        SigmaX2 = set()
        while cur != None:
            nxt = cur.next
            v = (y, cur.data[1])
            if GX.has_node(v):
                if G.has_node(v): 
                    X2.add(cur.data)
                else:
                    X2.append_right(Y2.pop(cur))
            elif not G.has_node(v):
                Y2.pop(cur)
            cur = nxt
        LX[y] = X2

    if any(x.head == None for x in LX):
        LOG.debug('discarding instance only containing non-universal markers: %s' %map(str, LX))
        LX = [None] * len(LX)
        GX.clear()

    if any(x.head == None for x in L):
        LOG.debug('discarding instance only containing non-universal markers: %s' %map(str, L))
        for x in xrange(len(L)):
            L[x] = None
            G.clear() 
    return LX, GX


def findTeams(L, G, ref, delta):

    res = list()
   
    queue = list()
    queue.append((L, G))
    while queue:
        L, G  = queue.pop()
        B = smallMax(L, delta)

        if L[B[0]].head is B[1] and L[B[0]].tail is B[2]:
            LOG.info('found team %s' %', '.join('G%s:%s-%s (%s)' %(x,
                L[x].head, L[x].tail, len(list(L[x]))) for x in
                xrange(len(L))))
            res.append((L, G)) 
        else:
            LX, GX = division(L, G, B, ref)
            if LX[0] != None and L[0] != None:
                LOG.debug('splitting sequence into %s and %s' %(\
                        ', '.join('G%s:%s-%s (%s)' %(x, L[x].head, L[x].tail, \
                        len(list(L[x]))) for x in xrange(len(L))), \
                        ', '.join('G%s:%s-%s (%s)' %(x, LX[x].head, LX[x].tail, \
                        len(list(LX[x]))) for x in xrange(len(LX)))))
                queue.append((L, G))
                queue.append((LX, GX))
            elif LX[0] != None:
                queue.append((LX, GX))
            elif L[0] != None:
                queue.append((L, G))
    return res


def fixIndels(L, G, gene_orders, g2pos, dists, g_counter, new_markers, genomes,
        ref, delta):
    bounds = getBounds(L, g2pos, delta)
    indels = identifyIndels(L, bounds, gene_orders)

    if sum(map(lambda x: any(x) and 1 or 0, indels)) < 2:
        return L, G, gene_orders, dists, g2pos, g_counter, new_markers

    Rs = computeR(bounds, gene_orders, g2pos, dists, genomes)

    Gindel = nx.Graph()
    for x, y in combinations(xrange(len(genomes)), 2):
        for i, j in product(xrange(len(indels[x])), xrange(len(indels[y]))):
            for gj in indels[y][j]:
                ints = indels[x][i].intersection(dists[(genomes[y], \
                        genomes[x])].get(gj, dict()).keys())
                Gindel.add_edges_from(((y, j, gj), (x, i, gi)) for gi in ints)

    for C in nx.connected_components(Gindel):
        Gs = set(map(lambda x: x[0], C))
        if len(Gs) <= len(genomes)/2.:
            LOG.debug('mark genes %s as insertion.' %(', '.join(map(lambda x: \
                    str(x[2]), C))))
            continue
        
        # it's a deletion! Re-insert the missing marker(s)
        inserted_genes = set()
        for x in sorted(set(xrange(len(genomes))).difference(Gs)):
            for i in xrange(len(bounds[x])):
                has_counterpart = False
                for y, _, gj in C:
                    for gi in dists[(genomes[y], genomes[x])].get(gj, \
                            dict()).keys():
                        p = g2pos[(x, gi)]
                        if p >= bounds[x][i][0] and \
                                p <= bounds[x][i][1]:
                            inserted_genes.add((x, i, gi))
                            inserted = __insert_if_missing__(L[x], x, g2pos, gi)
                            if inserted:
                                LOG.debug(('re-inserting marker %s into ' + \
                                        'genome %s') %(gi, x))
                            has_counterpart = True
                if has_counterpart:
                    continue

                # genes really missing from the syntenic block 
                r, y, k = max((Rs[(x, k)][(y, k)]['weight'], y, k) for \
                        (y, k) in Rs.neighbors((x, i)) if y in Gs)
                dels = sorted(c[2] for c in C if c[:2] == (y, k))
                
                distyx = dists[(genomes[y], genomes[x])]
                # XXX the identification of the insertion point needs much
                # improvement
                for d in dels:
                    left_neighbor = g2pos[(y, d)]-1
                    while left_neighbor >= bounds[y][k][0] and (not \
                            distyx.has_key(gene_orders[y][left_neighbor]) or \
                            not any(g2pos[(x, l)] >= bounds[x][i][0] and \
                            g2pos[(x, l)] <= bounds[x][i][1] for l in \
                            distyx[gene_orders[y][left_neighbor]].keys())):
                        left_neighbor -= 1

                    ins_pos_l = set()
                    ins_pos_r = set()
                    ins_w = dict()
                    if left_neighbor < bounds[y][k][0]:
                        LOG.debug(('unable to find left-bound insertion ' + \
                                'location for gene %s into genome %s') %(d, x))
                    else:
                        for l, (di, w) in distyx[gene_orders[y][ \
                                left_neighbor]].items():
                            if g2pos[(x, l)] >= bounds[x][i][0] and \
                                    g2pos[(x, l)] <= bounds[x][i][1]:
                                p = di == '+' and g2pos[(x, l)]+1 or \
                                        g2pos[(x, l)]
                                ins_pos_l.add(p)
                                ins_w[p] = max(ins_w.get(p, 0), w)

                    right_neighbor = g2pos[(y, d)]+1
                    while right_neighbor <= bounds[y][k][1] and (not \
                            distyx.has_key(gene_orders[y][right_neighbor]) or \
                            not any(g2pos[(x, l)] >= bounds[x][i][0] and \
                            g2pos[(x, l)] <= bounds[x][i][1] for l in \
                            distyx[gene_orders[y][right_neighbor]].keys())):
                        right_neighbor += 1

                    if right_neighbor > bounds[y][k][1]:
                        LOG.debug(('unable to find right-bound insertion ' + \
                                'location for gene %s into genome %s') %(d, x))
                    else:
                        for l, (di, w) in distyx[gene_orders[y][ \
                                right_neighbor]].items():
                            if g2pos[(x, l)] >= bounds[x][i][0] and \
                                    g2pos[(x, l)] <= bounds[x][i][1]:
                                p = di == '+' and g2pos[(x, l)] or \
                                        g2pos[(x, l)] +1
                                ins_pos_r.add(p)
                                ins_w[p] = max(ins_w.get(p, 0), w)
                   
                    if not ins_w:
                        LOG.warning(('unable to find any insertion ' + \
                                'location for gene %s into genome %s') %(d, x))
                        C.clear()
                        inserted_genes.clear()
                        break

                    ins_int = ins_pos_l.intersection(ins_pos_r)
                    ins_pos = None
                    if ins_int:
                        ins_pos = ins_int.pop()
                    else:
                        ins_pos, _ = max(ins_w.items(), key=lambda x: x[1])

                    g_counter[x] += 1
                    ig = (gene_orders[x][ins_pos][0], g_counter[x])
                    new_markers[x].append(ig)
                    gene_orders[x].insert(ins_pos, ig)

                    # update g2pos
                    g2pos[(x, ig)] = ins_pos
                    for p in xrange(ins_pos+1, len(gene_orders[x])):
                        if gene_orders[x][p] != CONTIG_BOUNDARY:
                            g2pos[(x, gene_orders[x][p])] = p
                    # update bounds
                    bounds[x][i]= (bounds[x][i][0], bounds[x][i][1]+1)
                    __insert_if_missing__(L[x], x, g2pos, ig)                        
                    inserted_genes.add((x, i, ig))
                    LOG.info(('inserted missing (deleted) gene %s into ' + \
                            'genome %s') %(ig, x)) 

        for x, i, gi in C:
            inserted = __insert_if_missing__(L[x], x, g2pos, gi)
            if inserted:
                LOG.debug('re-inserting marker %s into genome %s' %(gi, x))
            inserted_genes.add((x, i, gi))

        for _, _, gi in (c for c in inserted_genes if c[0] == ref):
            for x, _, gj in (c for c in inserted_genes if c[0] != ref):
                G.add_edge((ref, gi), (x, gj))

    return L, G, gene_orders, dists, g2pos, g_counter, new_markers

def __insert_if_missing__(Lx, x, g2pos, ig):
    ins_pos = g2pos[(x, ig)]
    n = Lx.head
    while n != None and g2pos[(x, n.data[1])] < ins_pos:
        n = n.next
    if n == None:
        new_n = dnode(data=(Lx.tail.data[0]+1, ig), prev=Lx.tail)
        Lx.head = new_n
        Lx.tail.next = new_n
        Lx.tail = new_n
        LOG.debug(('inserting gene %s at the end of the gene' + \
                ' team data structure') %str(ig))
        return True
    if n.data[1] == ig:
        return False
    p = n.data[0]
    pc = 1.
    new_n = dnode(prev=n.prev)
    new_n.next = n
    n.prev = new_n
    if new_n.prev:
        p += new_n.prev.data[0]
        pc += 1.
        new_n.prev.next = new_n
    else:
        Lx.head = new_n
    new_n.data = (p/pc, ig)

    return True

def identifyIndels(L, bounds, gene_orders):
    indels = [[set() for _ in bx] for bx in bounds]
    for x in xrange(len(L)):
        for k in xrange(len(bounds[x])): 
            start, end = bounds[x][k]
            indels[x][k].update(gene_orders[x][start:end+1])
    
    for x in xrange(len(L)):
        cur = L[x].head
        b = 0
        if cur.data[1] != gene_orders[x][bounds[x][b][0]]:
            LOG.fatal(('gene team bounds %s do not match with L data ' + \
                    'structure of genome %s.') %(str(bounds[x][b])), x)
            indels[x][b].clear()
            continue
        while cur != None:
            indels[x][b].discard(cur.data[1])
            if cur.data[1] == gene_orders[x][bounds[x][b][1]]:
                b += 1
            cur = cur.next

    return indels


def getBounds(L, g2pos, delta):

    bounds = list()

    for x in xrange(len(L)):
        b = list()
        end = start = L[x].head
        while end != None:
            if not end.next or end.next.data[0]-end.data[0]> delta:
                si = g2pos[(x, start.data[1])]
                ei = g2pos[(x, end.data[1])]
                b.append((si, ei))
                start = end.next
            end = end.next
        bounds.append(b) 
    return bounds


def computeR(bounds, gene_orders, g2pos, dists, genomes):

    G = nx.Graph()

    for x, y in combinations(xrange(len(genomes)), 2):
        for i in xrange(len(bounds[x])):
            sx, ex = bounds[x][i]
            for j in xrange(len(bounds[y])):
                sy, ey = bounds[y][j]
                xy = list()
                for cx in xrange(sx, ex+1):
                    for gy in dists[(genomes[x], \
                            genomes[y])].get(gene_orders[x][cx], \
                            dict()).keys():
                        cy = g2pos[(y, gy)]
                        if cy >= sy and cy <= ey:
                            xy.append((cx-sx, cy-sy))
                G.add_edge((x, i), (y, j), weight=r(xy))
    return G


def r(xy):
    n = len(xy)
    sx = sum(map(lambda x: x[0], xy))
    sy = sum(map(lambda x: x[1], xy))
    divisor = (sqrt(n * sum(map(lambda x: x[0]**2, xy)) - sx**2) *\
            sqrt(n*sum(map(lambda x: x[1]**2, xy))-sy**2))
    if divisor == 0.0:
        return 1
    return (n * sum(map(lambda x: x[0]*x[1], xy))-sx*sy)/divisor


def constructCIDS(L, G, ref, delta):

    ids = [x for x in xrange(len(L)) if x != ref]
    bounds = [list() for _ in L]
    g2pos = dict()
    gene_orders = [list() for _ in L]
    pos = dict(((ref, y), list()) for y in ids)
    pos.update(((y, ref), list()) for y in ids)

    p = 0
    prev = -delta
    cur = L[ref].head
    prev_b = None
    while cur != None:
        if cur.data[0] - prev > delta:
            if prev_b != None:
                bounds[ref].append((prev_b, len(gene_orders[ref])-1))
            gene_orders[ref].append(CONTIG_BOUNDARY)
            prev_b = len(gene_orders[ref])
            for y in ids:
                pos[(ref, y)].append(list())
            p += 1
        gene_orders[ref].append(cur.data[1])
        for y in ids:
            pos[(ref, y)].append(list())
        g2pos[(ref, cur.data[1])] = p
        p += 1
        prev = cur.data[0]
        cur = cur.next
    if prev_b != None and prev_b != len(gene_orders[ref]):
        bounds[ref].append((prev_b, len(gene_orders[ref])-1))
    gene_orders[ref].append(CONTIG_BOUNDARY)
    for y in ids:
        pos[(ref, y)].append(list())

    for x in ids:

        Gx = list()
        p = 0
        prev = -delta
        prev_b = None
        cur = L[x].head
        while cur != None:
            if cur.data[0] - prev > delta:
                if prev_b != None:
                    bounds[x].append((prev_b, len(gene_orders[x])-1))
                gene_orders[x].append(CONTIG_BOUNDARY)
                prev_b = len(gene_orders[x])
                pos[(x, ref)].append(list())
                p += 1
            pos[(x, ref)].append(list())
            gene_orders[x].append(cur.data[1])
            for gi in sorted(gy for (y, gy) in G.neighbors((x, cur.data[1])) \
                    if y==ref):
                if g2pos.has_key((ref, gi)):
                    pref = g2pos[(ref, gi)]
                    pos[(x, ref)][-1].append(pref)
                    z = bisect(pos[(ref, x)][pref], p)
                    if not z or pos[(ref, x)][pref][z-1] != p:
                        pos[(ref, x)][pref].insert(z, p)
            p += 1 
            prev = cur.data[0]
            cur = cur.next
        if prev_b != None and prev_b != len(gene_orders[x]):
            bounds[x].append((prev_b, len(gene_orders[x])-1))
        gene_orders[x].append(CONTIG_BOUNDARY)
        pos[(x, ref)].append(list())

    return gene_orders, pos, bounds


def getCharEnvJ(pos2, My, i, j, p, stop_at_left_bound=None):

    k = p 
    l = p+1

    cj = i-1
    DELAYED = set()

    left_bound = max(stop_at_left_bound or 0, 0)
    while k >= left_bound and My[k] < len(pos2[k]) and pos2[k][-1] >= i and \
            pos2[k][My[k]]<= j:
        x = My[k]
        while x < len(pos2[k]) and pos2[k][x] <= j:
            if pos2[k][x] > cj:
                DELAYED.add(pos2[k][x])
            x += 1
        while cj + 1 in DELAYED:
            DELAYED.remove(cj+1)
            cj += 1

        while l < len(pos2) and My[l] < len(pos2[l]) and pos2[l][-1] >= i and \
                pos2[l][My[l]]<= j:
            x = My[l]
            while x < len(pos2[l]) and pos2[l][x] <= j:
                if pos2[l][x] > cj:
                    DELAYED.add(pos2[l][x])
                x += 1
            while cj + 1 in DELAYED:
                DELAYED.remove(cj+1)
                cj += 1
            l += 1
        k -= 1
    return cj


def getCharEnvI(pos2, i, j, p, stop_at_left_bound=None):

    k = p 
    l = p+1

    ci = j+1
    DELAYED = set()

    left_bound = max(stop_at_left_bound or 0, 0)
    while k >= left_bound and len(pos2[k]) and pos2[k][-1] >= i and pos2[k][0] <= j:
        x = 0
        while x < len(pos2[k]) and pos2[k][x] <= j:
            if pos2[k][x] >= i and pos2[k][x] < ci:
                DELAYED.add(pos2[k][x])
            x += 1
        while ci - 1 in DELAYED:
            DELAYED.remove(ci-1)
            ci -= 1

        while l < len(pos2) and len(pos2[l]) and pos2[l][-1] >= i and pos2[l][0] <= j:
            x = 0 
            while x < len(pos2[l]) and pos2[l][x] <= j:
                if pos2[l][x] >= i and pos2[l][x] < ci:
                    DELAYED.add(pos2[l][x])
                x += 1
            while ci - 1 in DELAYED:
                DELAYED.remove(ci-1)
                ci -= 1
            l += 1
        k -= 1
    return ci


def getWitness(pos2, i, j, p, My=None, stop_at_left_bound=None):
    
    CHAR = [0] * (j-i+1)
    c = 0

    k = p 
    l = p+1

    left_bound = max(stop_at_left_bound or 0, 0)
    while k >= left_bound and len(pos2[k]) and pos2[k][0]<= j:
        x = My!=None and My[k] or 0
        while x < len(pos2[k]) and pos2[k][x] <= j:
            if pos2[k][x] >= i and not CHAR[j-pos2[k][x]]:
                CHAR[j-pos2[k][x]] = 1
                c += 1
            x += 1
        if x < 1 or pos2[k][x-1] < i:
            break
        k -= 1

    while l < len(pos2) and len(pos2[l]) and pos2[l][0] <= j:
        x = My!=None and My[l] or 0
        while x < len(pos2[l]) and pos2[l][x] <= j:
            if pos2[l][x] >= i and not CHAR[j-pos2[l][x]]:
                CHAR[j-pos2[l][x]] = 1
                c += 1
            x += 1
        if x < 1 or pos2[l][x-1] < i:
            break
        l += 1

    if c != j-i+1:
        low_i = j
        while low_i >= i and CHAR[j-low_i]:
            low_i -= 1
        high_j = i
        while high_j <= j and CHAR[j-high_j]:
            high_j += 1
        k, l, i, j = 0, -1, low_i+1, high_j-1
        
    return (k, l, i, j)


def getCharsetsJ(pos, M, ids, ref, i, j_max = None, p_set = None, enforce_same_j=False):

    j = j_max == None and maxint-1 or j_max
    prev_j =  j+1
    charsets = None
    while prev_j > j  and i < j:
        prev_j = j
        charsets = [list() for _ in ids]
        jys = [list() for _ in ids]
        for y in xrange(len(ids)):
            if p_set != None:
                ps = p_set[y]
            else:
                ps = pos[(ref, ids[y])][i]
            prev_p = -1
            for p in ps:
                k, l = 0, -1
                if p-1 > prev_p:
                    jp = j
                    while l-k-1 < 1:
                        jp = getCharEnvJ(pos[(ids[y], ref)], M[ids[y]], i, jp,
                                p, stop_at_left_bound=prev_p)
                        k, l, _, jp = getWitness(pos[(ids[y], ref)], i, jp, p,
                                My=M[ids[y]], stop_at_left_bound=prev_p)
                if p-1 != prev_p or k > prev_p:
                    charsets[y].append((p, k+1, l-1, jp))
                    jys[y].append(jp)
                else:
                    _, _k, _l, _jp = charsets[y][-1]
                    charsets[y].append((p, _k, _l, _jp))
                    jys[y].append(jys[y][-1])
                prev_p = p
        
        if enforce_same_j:
            j = min(map(min, jys))
        else:
            j = min(map(max, jys))
    return j, charsets


def getCharsetsI(pos, ids, ref, j, i_min=None, p_set=None,
        enforce_same_i=False):

    i = i_min != None and i_min or 0
    prev_i = i-1
    charsets = None
    while prev_i < i  and i <= j:
        prev_i = i
        charsets = [list() for _ in ids]
        iys = [list() for _ in ids]
        for y in xrange(len(ids)):
            prev_p = -1
            if p_set != None:
                ps = p_set[y]
            else:
                ps = pos[(ref, ids[y])][j]
            for p in ps:
                k, l = 0, -1
                if p +1 > prev_p:
                    ip = i
                    while l-k-1 < 1:
                        ip = getCharEnvI(pos[(ids[y], ref)], ip, j, p,
                                stop_at_left_bound=prev_p)
                        k, l, ip, _ = getWitness(pos[(ids[y], ref)], ip, j, p,
                                stop_at_left_bound=prev_p)
                if p+1 != prev_p or k > prev_p:
                    charsets[y].append((p, k+1, l-1, ip))
                    iys[y].append(ip)
                else:
                    charsets[y].append(charsets[y][-1])
                    iys[y].append(iys[y][-1])
                prev_p = p

        if enforce_same_i:
            i = max(map(max, iys))
        else:
            i = max(map(min, iys))
    return i, charsets


def getIntervals(pos, n, ref, start=None, end=None, M=None):

    res = set() 

    # initialization
    ids = [x for x in xrange(n) if x != ref] 
    M = M or [[0] * len(pos.get((y, ref), ())) for y in xrange(n)]
    i = start != None and start or 0 
    end = end == None and len(pos[(ref, ids[0])]) or end

    # iterating through each position between i and end in the reference 
    while i < end:

        # do not start at an empty position
        if not any(not not pos[(ref, y)][i] for y in xrange(n) if \
                pos.has_key((ref, y))):
            i += 1
            continue

        #
        # enumerate largest intervals containing characters larger equal to i
        #

        ip, i_charsets = getCharsetsI(pos, ids, ref, i)
        j, j_charsets = getCharsetsJ(pos, M, ids, ref, i)

        hits = [[False] * len(x) for x in j_charsets]
        while i < j and not all(map(all, hits)):
            p_set = [[j_charsets[y][p][0] for p in xrange(len(j_charsets[y])) if
                not hits[y][p] and j_charsets[y][p][3] >= j] for y in
                xrange(len(j_charsets))]
          
            jp, jp_charsets = getCharsetsJ(pos, M, ids, ref, i, j_max=j,
                    p_set=p_set, enforce_same_j=True)

            if jp == j:
                prev_i = ip-1
                while prev_i != ip:
                    prev_i = ip
                    ip, ip_charsets = getCharsetsI(pos, ids, ref, i, i_min=ip,
                            p_set=p_set, enforce_same_i=True)
                

                # add identified common intervals to result list
                res_i, res_j = list(), list()
                for y in xrange(len(jp_charsets)):
                    for z in xrange(len(jp_charsets[y])):
                        if z < len(ip_charsets[y])-1 and jp_charsets[y][z] == \
                                jp_charsets[y][z+1] and ip_charsets[y][z] == \
                                ip_charsets[y][z+1]:
                            continue
                        int_i = (ids[y], ip_charsets[y][z][1],
                                ip_charsets[y][z][2])
                        if not res_i or res_i[-1] != int_i:
                            res_i.append(int_i)
                        int_j = (ids[y], jp_charsets[y][z][1],
                                jp_charsets[y][z][2])
                        if not res_j or res_j[-1] != int_j:
                            res_j.append(int_j)
                insort(res_i, (ref, ip, i))
                insort(res_j, (ref, i, j))
                res.update((tuple(res_i), tuple(res_j)))

                # XXX only use largest interval for now
                break

                # update hit data structure
                for y in xrange(len(j_charsets)):
                    x = 0
                    for p, _, _, _ in jp_charsets[y]:
                        while x < len(j_charsets[y]) and j_charsets[y][x][0] != p:
                            x += 1
                        hits[y][x] = True
                j = max([max([j_charsets[y][x] for x in xrange(len(j_charsets[y])) \
                        if not hits[y][x]] + [-1])])
            else:
                j = jp
        for y in xrange(len(M)):
            if pos.has_key((ref, y)):
                for p in pos[(ref, y)][i]:
                    M[y][p] += 1
        i += 1

    return list(res)
        

def getIntervals_parallel(in_queue, out_queue, n, ref, pos):

    # use shared variable (i.e. the list pqr_trees) to communicate the result
    # back to the main thread
    M = [[0] * len(pos.get((y, ref), ())) for y in xrange(n)] 

    prev_i = -1
    while in_queue.qsize() > 0:
        try:
            i = in_queue.get(timeout=1)
        except Empty:
            continue
        for y in xrange(len(M)):
            if pos.has_key((ref, y)):
                for ii in xrange(prev_i+1, i):
                    for p in pos[(ref, y)][ii]:
                        M[y][p] += 1
        ints = getIntervals(pos, n, ref, i, i+1, M)
        for intt in ints:
            out_queue.put(intt)

        in_queue.task_done()
        prev_i = i 


def identifyStrongIntervals(ints, ref):

    # construct 4n array of interval bounds
    bounds = list()

    int_map = dict()

    for x in xrange(len(ints)):
        y = ref
        while ints[x][y][0] <= ref:
            if ints[x][y][0] == ref:
                _, lb, rb = ints[x][y]
                bounds.extend(((lb, '('), (rb, ')')))
                if not int_map.has_key((lb, rb)):
                    int_map[(lb, rb)] = set()
                int_map[(lb, rb)].update((ints[x][z][0],ints[x][z][1], \
                        ints[x][z][2]) for z in xrange(len(ints[x])) if z != y)
            y += 1
    bounds.sort()

    strong_cis = set()
    S = list()
    for a, b in bounds:
        if b == '(':
            S.append(a)
        else:
            strong_cis.add((S.pop(), a))

    # remove singleton strong intervals that have no links
    res = set()
    for (l, r) in strong_cis:
        if l != r or int_map.has_key((l, r)):
            links = sorted(int_map.get((l, r), ()))
            i = 1
            while i < len(links):
                if links[i-1][0] == links[i][0] and \
                        links[i-1][2] >= links[i][2]:
                    del links[i]
                else:
                    i += 1
            res.add(tuple(chain(((ref, l, r), ), links)))
    return res


def constructInclusionTree(strong_cis, pos, gene_orders, bounds, n, ref):

    subtrees = list()
    for l, r in bounds[ref]:
        c_strong_cis = set(c for c in strong_cis if c[0][1] >= l and c[0][2] <= r)
        singletons = set(x[0][1] for x in c_strong_cis if x[0][1] == x[0][2])

        for i in xrange(l, r+1):
            if i not in singletons:
                links = [(ref, i, i)]
                for y in xrange(n):
                    if y != ref and pos.has_key((ref, y)) and pos[(ref, y)]:
                        x = 0
                        # merge intervals of tandem duplicates
                        while x < len(pos[(ref, y)][i]):
                            s = x
                            while x < len(pos[(ref, y)][i])-1  and pos[(ref, \
                                    y)][i][x+1] == pos[(ref, y)][i][x]+1:
                                x += 1
                            links.append((y, pos[(ref, y)][i][s], pos[(ref, y)][i][x]))
                            x += 1
                c_strong_cis.add(tuple(links))


        I = sorted(c_strong_cis, cmp=lambda x, y: -cmp(x[0][2], y[0][2]) or
                cmp(x[0][1], y[0][1]))

        int_map = dict()
        # determine root
        if I[0][0][1] != l or I[0][0][2] !=r:
            I.insert(0, ((ref, l, r), ))
        F = node(I[0][0][1], I[0][0][2], links=tuple(I[0][1:]))
        int_map[F] = 0
        root = F
        k = 1
        while k < len(I):
            if I[k][0][1] >= I[int_map[F]][0][1] and I[k][0][2] <= \
                    I[int_map[F]][0][2]:
                nI = node(I[k][0][1], I[k][0][2], parent=F, \
                        links=tuple(I[k][1:]))
                F.children.append(nI)
                F = nI
                int_map[F] = k
                k += 1
            else:
                F = F.parent
        subtrees.append(root)
    return subtrees


def gos2mseq(goss, gMap, id2genomes):

    res = [list() for _ in id2genomes]

    for gos in goss:
        for y in xrange(len(gos)):
            marker_order = gMap[id2genomes[y]][GM_ACTV_GNS_KEY]
            res[y].append(list())
            for chrx, x in gos[y]:
                if (chrx, x) == CONTIG_BOUNDARY:
                    res[y][-1].append(CONTIG_BOUNDARY_KEY)
                else:
                    res[y][-1].append(marker_order[x-1])
    return res


def parseInput(args, options):
    isMultiChrom = False
    gene_orders_dict = dict()
    dists = dict()

    for f in args:
        isMultiChrom_f, dist = readDists(open(f), 0)
        isMultiChrom = isMultiChrom or isMultiChrom_f

        gname1, gname2  = basename(f).split('.', 1)[0].split('_', 1)
        if not gene_orders_dict.has_key(gname1):
            gene_orders_dict[gname1] = set()
        if not gene_orders_dict.has_key(gname2):
            gene_orders_dict[gname2] = set()

        if not dist:
            LOG.warning('distance table in file %s is empty.' %f)

        revDist = reverseDistMap(dist)

        gene_orders_dict[gname1].update(dist.keys())
        gene_orders_dict[gname2].update(revDist.keys())

        dists[(gname1, gname2)] = dist
        dists[(gname2, gname1)] = revDist

    genomes2id = dict(izip(sorted(gene_orders_dict.keys()),
        xrange(len(gene_orders_dict))))

    gene_orders = [None] * len(genomes2id)
    for gname, gorder in gene_orders_dict.items():
        gene_orders[genomes2id[gname]] = establish_linear_genome_order(gorder)
    
    # discard isMultiChrom information for now
    return genomes2id, gene_orders, dists

if __name__ == '__main__':

    usage = '%prog [options] <PAIRWISE DIST FILE 1> ... <PAIRWISE DIST FILE N>'
    parser = OptionParser(usage=usage)
    parser.add_option('-d', '--delta', dest='delta', default=0, type='int',
            help='Approximate weak common intervals parameter, ' + \
                    'allows for <delta> indels [default: %default]')
    parser.add_option('-r', '--reference', dest='ref', default=DEFAULT_REF,
            type='str', help='Reference genome [default: %default]')
    parser.add_option('-c', '--inv_coverage', dest='cov',
            default=DEFAULT_COVERAGE_INV, type='int',
            help='Minimum inverse covaerge a delta-team must have in order' + \
                    ' to be considered as syntenic context block [default:' + \
                    ' %default]')
    parser.add_option('-o', '--output', dest='outFile', default='',
            type='str', help='Specify output file name. If none is given, ' + \
                    'output is piped to stdout [default: %default]')
    parser.add_option('-t', '--threads', dest='noThreads', default=cpu_count(),
            type='int', help='Max number of parallel threads [default: ' + \
            '%default]')

    (options, args) = parser.parse_args()
    
    if len(args) < 1:
        parser.print_help()
        exit(1)

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.ERROR)
    ch.setFormatter(logging.Formatter('!! %(message)s'))
    cf = logging.FileHandler('%s_d%s.log' %(LOG_FILENAME, options.delta), mode='w', delay=True)
    cf.setLevel(logging.DEBUG)
    cf.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(cf)
    LOG.addHandler(ch)

    genomes2id, gene_orders, dists = parseInput(args, options)
    g2pos = dict()
    for x in xrange(len(gene_orders)):
        for i in xrange(len(gene_orders[x])):
            g2pos[(x, gene_orders[x][i])] = i

    id2genomes = map(lambda x: x[0], sorted(genomes2id.items(), key=lambda x:
        x[1]))

    if options.ref not in genomes2id:
        if options.ref == DEFAULT_REF:
            LOG.fatal('you must specify a reference genome!')
        else:
            LOG.fatal(('reference genome %s is not contained in dataset.' + \
                    'Exiting.') %options.ref)
        parser.print_help()
        exit(1)

    ref = genomes2id[options.ref]

    gMapFile = join(dirname(args[0]), GENOME_MAP_FILE)
    if not isfile(gMapFile):
        LOG.fatal(('Unable to locate genome map file at assumed location ' + \
                '%s. Exiting.') %gMapFile)
        exit(1)
    gMap = readGenomeMap(open(gMapFile))

    G = nx.Graph()
    for Gx, Gy in combinations(id2genomes, 2):
        x = genomes2id[Gx]
        y = genomes2id[Gy]
        for gx, gys in dists[(Gx, Gy)].items():
            for gy in gys.keys():
                G.add_edge((x, gx), (y, gy))

    removeNonUniversalGenes(G, len(id2genomes))
    L, G = constructTeamDS(gene_orders, G, ref, options.delta)
    teams = findTeams(L, G, ref, options.delta)

    CI_instances = list()

    g_counter = [go[-2][1] for go in gene_orders]
    new_markers = [list() for _ in id2genomes]
    root = node(1, len(gene_orders[ref])-2)

    goss = list()
    strong_ciss = list()
    ciss = list()
    for z in xrange(len(teams)):
        L, G = teams[z]
        if all(len(set(x)) > 1 for x in L):
#            L, G, gene_orders, dists, g2pos, g_counter, new_markers = fixIndels(L, G, \
#                    gene_orders, g2pos, dists, g_counter, new_markers, \
#                    id2genomes, ref, options.delta)

            gos, pos, bounds = constructCIDS(L, G, ref, options.delta)
            cov = reduce(lambda x,y: x+y, map(lambda z: [float(gos[z][w[1]][1] \
                    -gos[z][w[0]][1]+1)/(w[1]-w[0]+1) for w in bounds[z]], \
                    xrange(len(bounds))))
            if len(set(map(len, bounds))) != 1 or any(x > options.cov for x
                    in cov):
                LOG.info(('Dismissing marker team due to insufficient ' + \
                        'coverage of markers (1/cov > %s) in at least one ' + \
                        'occurrence: %s. Occurrence lengths: %s') %(options.cov, 
                            str(cov), str(reduce(lambda x,y: x+y, map(lambda z:
                                [w[1]-w[0]+1 for w in z], bounds)))))
                continue

            LOG.info('enumerating intervals for team %s' %', '.join(
                '%s:%s..%s' %(id2genomes[y], gos[y][1],
                    gos[y][-2]) for y in xrange(len(gos))))

            if options.noThreads > 1 and len(gos[ref]) > 100:
                LOG.debug('creating pool with %s threads' %options.noThreads)
                pool = Pool(processes=options.noThreads)
                in_queue = JoinableQueue()
                for i in xrange(1, len(gos[ref])-1):
                    in_queue.put(i)

                out_queue = Queue()
                for _ in xrange(options.noThreads):
                    p = pool.Process(target=getIntervals_parallel,
                            args=(in_queue, out_queue, len(id2genomes), ref,
                                pos))
                    p.start()
                getIntervals_parallel(in_queue, out_queue, len(id2genomes),
                        ref, pos)
                
                in_queue.close()
                in_queue.join()
                
                LOG.debug('closing pool')
                pool.close()

                cis = list() 
                while out_queue.qsize() > 0:
                    try: 
                        cis.append(out_queue.get())
                    except Empty:
                        pass
                
                out_queue.close()
                out_queue.join_thread()
            else:
                cis = getIntervals(pos, len(id2genomes), ref)

            strong_cis = identifyStrongIntervals(cis, ref)
            subtrees = constructInclusionTree(strong_cis, pos, gos, bounds,
                    len(id2genomes), ref)

            for st in subtrees:
                st.id = len(goss)
            goss.append(gos)
            strong_ciss.append(strong_cis)
            ciss.append(cis)

            root.children.extend(subtrees)
            for subtree in subtrees:
                subtree.parent = root

#    covered_markers = set(chain(*(gene_orders[ref][g2pos[ref][u.intt[0]]: \
#        g2pos[ref][u.intt[1]]+1] for u in root.children)))
#    singletons = set(gene_orders[1:-2]).difference(covered_markers)
#
#    i = 0
#    for g1i in singletons:
#        while g2pos[ref][u.children[i].intt[1]] > g2pos[
   
    if len(root.children) == 1 and root.intt == root.children[0].intt:
        root = root.children[0]
        root.parent = None

    LOG.info('DONE! writing hierarchy..')
    outDict = dict()

    outDict['genome_names'] = id2genomes 
    outDict['ref_id'] = ref
    outDict['recovered_markers'] = new_markers
    outDict['marker_seq_list'] = gos2mseq(goss, gMap, id2genomes)
    #outDict['intervals'] = strong_ciss
    outDict['raw_sbfs'] = ciss
    if options.outFile:
        outDict['pwsim_files'] = map(partial(relpath,
            start=dirname(options.outFile)), args)
    else:
        outDict['pwsim_files'] = map(abspath, args)

    outDict['sb_hierarchy'] = hierarchy2dict(root)

    out = stdout
    if options.outFile:
        out = open(options.outFile, 'w')

    print >> out, json.dumps(outDict)


