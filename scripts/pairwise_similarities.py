#!/usr/bin/env python

from itertools import combinations, izip, product, chain, repeat
from os.path import basename, dirname, isfile, join, relpath
from sys import stdout, stderr, exit, maxint, argv
from ConfigParser import ConfigParser
from optparse import OptionParser
from math import floor
import networkx as nx
from os import stat
import logging
import csv
import re

PAT_BLASTTBL = re.compile('^(.*\.\w+)\.blast\w*$', re.IGNORECASE)
PAT_CHR = re.compile('.*\|chromosome\|([^\|]+)(\|.*|$)')
PAT_POS = re.compile('^.*\|([0-9]+):([0-9]+)(\||$)')
PAT_STRAND = re.compile('.*\|strand\|([^\|]+)(\|.*|$)')
PAT_ID = re.compile('^(\w+)_(\d+)\|.*$')

DIRECTION_CRICK_STRAND = '+'
DIRECTION_WATSON_STRAND = '-'

TELOMERE_WEIGHT = 1
TELOMERE_START = 0
TELOMERE_END = maxint

GENOME_MAP_FILE = 'genome_map.cfg'

GM_FILE_KEY = 'fasta_file'
GM_ACTV_GNS_KEY = 'active_genes'
GM_CHR_KEY = 'chromosomes'
GM_CHR_CONFORMATION = 'conformation'
GM_CHR_CIRCULAR = 'circular'
GM_CHR_LINEAR = 'linear'

DEFAULT_MIN_GENOME_NO=1
DEFAULT_MIN_CONTIG_LEN=1
DEFAULT_MAX_ITERATIONS=100

LOG = logging.getLogger(__name__)
LOG.setLevel(logging.DEBUG)
LOG_FILENAME = '%s.log' %basename(argv[0]).rsplit('.py', 1)[0]

def writePairwiseSimilarities(G, gene2genome, gMap, gNames, outDir=None,
        selfComparison=False):

    gPos = dict((gNames[i], i) for i in xrange(len(gNames)))
    activeGenes= dict((k, dict(izip(v[GM_ACTV_GNS_KEY], xrange(1, \
            len(v[GM_ACTV_GNS_KEY])+1)))) for k, v in gMap.items())
    
    dists = dict(((Gx, Gy), list()) for Gx, Gy in combinations(gNames, 2))

    if selfComparison:
        dists.update(((Gx, Gx), list()) for Gx in gNames)

    for g1, g2, data in G.edges_iter(data=True):
        Gx = gene2genome[g1]
        Gy = gene2genome[g2]

        if not selfComparison and Gx == Gy:
            continue

        m = PAT_CHR.match(g1)
        chr1 = 0
        if m:
            chr1 = m.group(1).split('::', 1)[0] 
        m = PAT_CHR.match(g2)
        chr2 = 0
        if m:
            chr2 = m.group(1).split('::', 1)[0]

        bits1self = G.node[g1]['bitscore']
        bits2self = G.node[g2]['bitscore']
        hsp = G[g1][g2]['hsp']
        bits1to2 = hsp[(Gx, Gy)]['bitscore']
        bits2to1 = hsp[(Gy, Gx)]['bitscore']
        strand = hsp[(Gx, Gy)]['orient']
        s = (bits1to2 + bits2to1)/(bits1self + bits2self)
        
        if gPos[Gx] <= gPos[Gy]:
            dists[(Gx, Gy)].append((chr1, activeGenes[Gx][g1], chr2, activeGenes[Gy][g2],
                strand, s))
        if gPos[Gx] >= gPos[Gy]:
            dists[(Gy, Gx)].append((chr2, activeGenes[Gy][g2], chr1, activeGenes[Gx][g1],
                strand, s))
       
    for (Gx, Gy), pwDists in dists.items():
        pwDists.sort()
        out = open(join(outDir or '', '%s_%s.sim' %(Gx, Gy)), 'w')
        for x in pwDists:
            print >> out, '\t'.join(map(str, x))
        out.flush()
        out.close()


def constructGeneGraph(blastfiles, gene2genome, selfComparison=False):
    G = nx.Graph()
    for f in blastfiles:
        LOG.info('Reading BLAST file %s...' %f)
        for line in csv.reader(open(f), delimiter='\t'):
            if gene2genome.has_key(line[0]) and gene2genome.has_key(line[1]):
                if line[0] == line[1]:
                    if not G.has_node(line[0]):
                        G.add_node(line[0])
                    G.node[line[0]]['bitscore'] = max(G.node[line[0]].get('bitscore', 0), float(line[11]))
                    continue

                G0 = gene2genome[line[0]]
                G1 = gene2genome[line[1]]

                if G0 == G1 and not selfComparison:
                    continue

                bitscore = float(line[11])
                if not G.has_edge(line[0], line[1]):
                    G.add_edge(line[0], line[1], hsp=dict())
                if not G[line[0]][line[1]]['hsp'].has_key((G0, G1)) or \
                        G[line[0]][line[1]]['hsp'][(G0, G1)]['bitscore'] < bitscore:
                    hsp = G[line[0]][line[1]]['hsp']
                    if not hsp.has_key((G0, G1)):
                        hsp[(G0, G1)] = dict()
                    hsp[(G0, G1)]['bitscore'] = bitscore

                    m0 = PAT_STRAND.match(line[0])
                    m1 = PAT_STRAND.match(line[1])

                    strand = 1
                    if m0 and m1:
                        strand = m0.group(1) == m1.group(1) and 1 or -1
                    else:
                        strand = int(line[8]) > int(line[9]) and \
                            -1 or 1

                    hsp[(G0, G1)]['orient'] = strand

                # store for each node the best bit score to a gene in ANOTHER
                # genome
                G.node[line[0]][G1] = max(G.node[line[0]].get(G1, 0), bitscore)
                G.node[line[1]][G0] = max(G.node[line[1]].get(G0, 0), bitscore)
    LOG.info('All BLAST files read!')
    for u,v in G.edges():
        if (gene2genome[u] != gene2genome[v] and len(G[u][v]['hsp']) < 2) or \
                not sum(map(lambda x: x['orient'], G[u][v]['hsp'].values())) \
                or not G.node[u].has_key('bitscore') or not \
                G.node[v].has_key('bitscore'):
            G.remove_edge(u, v)
    return G

def getActiveGeneSeqs(genome_map_path):

    from Bio import SeqIO
    gMap_loc = dirname(genome_map_path)
    gMap = readGenomeMap(open(genome_map_path))

    res = dict()

    for Gx, val in gMap.items():
        res[Gx] = list()
        for rec in SeqIO.parse(open(join(gMap_loc, val[GM_FILE_KEY])), 'fasta'):
            if rec.id in val[GM_ACTV_GNS_KEY]:
                res[Gx].append(rec)
    return res
        
def pruneGraph(G, gene2genome, minContigLen, minGenomeNo, stringency, max_iterations):

    i = 0
    hasChanged = True
    LOG.info('Start pruning...')
    while i < max_iterations and hasChanged:
        LOG.info('Iteration %s, current number of genes: %s' %(i, len(G)))
        hasChanged = False
        i += 1
        contigs = dict()
        for v in G.nodes():
            G0 = gene2genome[v]
            if stringency > 0:
                for u in list(G.neighbors_iter(v)):
                    G1 = gene2genome[u]
                    if G[v][u]['hsp'][(G0, G1)]['bitscore'] < stringency * \
                            max((G[vp][u]['hsp'].get((G1, G0), dict()).get('bitscore', \
                            0) for vp in G.neighbors_iter(u))):
                                
                        #G.node[u][G0]:
                        hasChanged = True
                        G.remove_edge(v, u)

            if minGenomeNo > 1 and len(set(map(gene2genome.get, \
                    G.neighbors_iter(v))).difference((gene2genome[v], ))) + 1 < minGenomeNo:
                G.remove_node(v)
                hasChanged = True
            elif minContigLen > 1:
                chrId = PAT_CHR.match(v).group(1)
                if not contigs.has_key((G0, chrId)):
                    contigs[(G0, chrId)] = set()
                contigs[(G0, chrId)].add(v)

        for (G0, chrId), vs in contigs.items():
            if len(vs) < minContigLen:
                hasChanged = True
                for v in vs:
                    G.remove_node(v)
    LOG.info('Finished pruning')
    return G


def constructGenomeMap(G, fastaFiles, outDir):
    res = dict() 
    for f in fastaFiles:
        gName = basename(f).rsplit('.', 1)[0]
        res[gName] = dict()
        res[gName][GM_FILE_KEY] = relpath(f, outDir or '.')
        res[gName][GM_CHR_KEY] = list()
        res[gName][GM_ACTV_GNS_KEY] = list()

        for x in open(join(outDir or '.', res[gName][GM_FILE_KEY])):
            if x.startswith('>'): 
                ident = x[1:x.find(' ')]
                if ident.find(',') >= 0:
                    LOG.fatal(('Erroneous gene identifier %s: ' + \
                            'character \',\'is not allowed in gene ' + \
                            'identifiers. Exiting') %ident)
                    exit(1)
                if G.has_node(ident):
                    res[gName][GM_ACTV_GNS_KEY].append(ident)
                    m = PAT_CHR.match(x)
                    if m: 
                        x_chr = m.group(1).strip()
                        if res[gName][GM_CHR_KEY] and res[gName][GM_CHR_KEY][-1] != x_chr:
                            res[gName][GM_CHR_KEY].append(x_chr)
        # remove unnecessary entry
        if not res[gName][GM_CHR_KEY]:
            del res[gName][GM_CHR_KEY]
    return res 

def writeGenomeMap(gMap, gNames, out):
    config = ConfigParser()
    for Gx in gNames:
        config.add_section(Gx)
        for k, v in gMap[Gx].items():
            if hasattr(v, '__iter__'):
                v = ','.join(map(str, v))
            config.set(Gx, k, v)
    config.write(out)

def readGenomeMap(data):
    config = ConfigParser()
    config.readfp(data)
    
    res = dict()
    for Gx in config.sections():
        res[Gx] = dict()
        for k, v in config.items(Gx):
            if k == GM_CHR_KEY or k == GM_ACTV_GNS_KEY:
                v = map(lambda x: x.strip(), v.split(','))
            res[Gx][k] = v
    return res

def readGenesFromFasta(blastfiles, fastaDir):

    res = dict()
    fastaFiles = list()
    
    for f in blastfiles:
        # read ids of all genes in current genome
        m = PAT_BLASTTBL.match(basename(f))
        if m == None:
            LOG.fatal(('BLAST file %s does not satisfy naming convention ' + \
                    '%s. Exiting') %(f, PAT_BLASTTBL.pattern))
            exit(1)

        fasFileName = join(fastaDir, m.group(1))
        if not isfile(fasFileName):
            LOG.fatal('Could not find required FASTA file %s. Exiting' %fasFileName)
            exit(1)

        fastaFiles.append(fasFileName)

        gName = m.group(1).rsplit('.', 1)[0]
        res[gName] = [x[1:x.find(' ')] for x in open(fasFileName) if
                x.startswith('>')]

    return res, fastaFiles

def readDists(data, edgeThreshold=None):

    res = dict()

    if edgeThreshold == None:
        edgeThreshold = 0

    hasMultipleChromosomes = False

    chr1 = '0'
    chr2 = '0'

    for line in csv.reader(data, delimiter='\t'):
        if not res:
            hasMultipleChromosomes = len(line) == 6

        if hasMultipleChromosomes:
            chr1 = line[0]
            if line[1] == 'TELOMERE_START':
                g1 = TELOMERE_START
            elif line[1] == 'TELOMERE_END':
                g1 = TELOMERE_END
            else:
                g1 = int(line[1])
            chr2 = line[2]
            if line[3] == 'TELOMERE_START':
                g2 = TELOMERE_START
            elif line[3] == 'TELOMERE_END':
                g2 = TELOMERE_END
            else:
                g2 = int(line[3])

            direction = line[4]
            edgeWeight = float(line[5])
        else:
            g1 = int(line[0])
            g2 = int(line[1])
            direction = line[2]
            edgeWeight = float(line[3])

        if edgeWeight < edgeThreshold:
            continue

        l0 = (chr1, g1)
        l1 = (chr2, g2)

        if l0 not in res:
            res[l0] = dict()
        # append mapping pos in mappedGenome and the weight of the corresponding edge
        res[l0][l1] = (direction == '1' and DIRECTION_CRICK_STRAND or \
                DIRECTION_WATSON_STRAND, edgeWeight)

    return hasMultipleChromosomes, res


def readDistsAndOrder(data, g1_chrs=None, g2_chrs=None, edgeThreshold=None, telomereWeight=None, addTelomeres=True):

    hasMultipleChromosomes, dists = readDists(data, edgeThreshold)

    if telomereWeight == None:
        telomereWeight = TELOMERE_WEIGHT

    g1_chromosomes = dict()
    g2_chromosomes = dict()

    for (chr1, g1i), g2ks in dists.items():
        if not g1_chromosomes.has_key(chr1):
            g1_chromosomes[chr1] = set()
        g1_chromosomes[chr1].add(g1i)

        for chr2, g2k in g2ks.keys():
            if not g2_chromosomes.has_key(chr2):
                g2_chromosomes[chr2] = set()
            g2_chromosomes[chr2].add(g2k)

    # construct genome order
    tel1, g1order = __establish_genome_order(g1_chromosomes, g1_chrs or dict())
    tel2, g2order = __establish_genome_order(g2_chromosomes, g2_chrs or dict())

    if addTelomeres:
        # add telomeres
        for (chr1,g1i,orient1), (chr2,g2k,orient2) in product(tel1, tel2):
            t1 = (chr1, g1i)
            t2 = (chr2, g2k)
            if not dists.has_key(t1):
                dists[t1] = dict()
            dists[t1][t2] = (orient1 == orient2 and DIRECTION_CRICK_STRAND or
                    DIRECTION_WATSON_STRAND, 1)

    return hasMultipleChromosomes, g1order, g2order, dists 


def __establish_genome_order(chromosomes, chr_info):
    g = list()
    telomeres = set()
    for k in sorted(chromosomes.keys()):
        # in case of doubt, always establish linear chromosome order
        if not chr_info.has_key(k) or chr_info[k] != GM_CHR_CIRCULAR:
            g.append((k, TELOMERE_START))
            telomeres.add((k, TELOMERE_START, DIRECTION_CRICK_STRAND))
            g.extend([(k, i) for i in sorted(chromosomes[k])])
            g.append((k, TELOMERE_END))
            telomeres.add((k, TELOMERE_END, DIRECTION_WATSON_STRAND))
    return telomeres, g


def reverseDistMap(dists):
    distsrev = dict()
    for g1i, x in dists.items():
        for g2j, d in x.items():
            if g2j not in distsrev:
                distsrev[g2j] = dict()
            distsrev[g2j][g1i] = d
    return distsrev


if __name__ == '__main__':

    usage = 'usage: %prog [options] <BLAST FILE 1> ... <BLAST FILE N>'
    parser = OptionParser(usage=usage)
    parser.add_option('-n', '--min_genome_number', dest='minGenomeNo', type=int,
            help='Minimum number of genomes for which each gene must share ' + \
                    'some similarity in [default: %default]',
                    default=DEFAULT_MIN_GENOME_NO)
    parser.add_option('-l', '--min_contig_len', dest='minContigLen', type=int,
            help='Minimum number of "non-zero similarity" genes contained ' + \
                    'in a contig to be contained in the output [default: ' + \
                    '%default]', default=DEFAULT_MIN_CONTIG_LEN)
    parser.add_option('-s', '--stringency', dest='stringency', default=0,
            type='float', help='Stringency threshold <stringency> \in ' + \
                    '[0, 1]. Edges with stringency lower than stringency' + \
                    'than given threshold value are removed from the ' + \
                    'input graph [default: %default]')
    parser.add_option('-N', '--numeric', dest='sortNumeric', action='store_true',
            help='sort input genomes numerically (rather than alphabetic)',
            default='')
    parser.add_option('-f', '--fasta_dir', dest='fastaDir', type=str,
            help='Directory in which fasta files are located that ' + \
                    'correspond to the input files. If nothing is ' + \
                    'specified, the files are assumed to be in the' + \
                    ' same directory as the blast files. [default:' + \
                    ' \'%default\']', default='')
    parser.add_option('-o', '--out_dir', dest='outDir', type=str,
            help='Output directory. If nothing is specified, output ' + \
                    'files are written to the current working directory.' + \
                    ' [default: \'%default\']', default='')
    parser.add_option('-i', '--max_iterations', dest='maxIterations', type=int,
            help='Maximum number of iterations for pruning the gene ' + \
                    'similarity graph [default: %default]',
                    default=DEFAULT_MAX_ITERATIONS)
    parser.add_option('-S', '--self_comparison', dest='selfComparison',
            action='store_true', default=False, help='Compute also gene ' + \
                    'similarities between the genomes themselves ' + \
                    '[default: %default]')

    (options, args) = parser.parse_args()

    if len(args) < 1 or len(map(isfile, args)) != len(args):
        parser.print_help()
        exit(1)

    # setup logging
    ch = logging.StreamHandler(stderr)
    ch.setLevel(logging.ERROR)
    ch.setFormatter(logging.Formatter('!! %(message)s'))
    cf = logging.FileHandler(LOG_FILENAME, mode='w', delay=True)
    cf.setLevel(logging.INFO)
    cf.setFormatter(logging.Formatter('%(levelname)s\t%(asctime)s\t%(message)s'))
    LOG.addHandler(cf)
    LOG.addHandler(ch)

    blastFiles = args
    if options.sortNumeric:
        _, blastFiles = zip(*sorted(map(lambda x: (int(''.join(re.findall('\d+',
            basename(x)))), x), blastFiles)))
    else:
        blastFiles.sort()

    
    genes, fastaFiles = readGenesFromFasta(blastFiles, options.fastaDir or
            dirname(args[0]))
    gNames = map(lambda x: basename(x).rsplit('.', 1)[0], fastaFiles)
    gene2genome = dict(chain(*(izip(v, repeat(k, len(v))) for k,v in
        genes.items())))
    G = constructGeneGraph(blastFiles, gene2genome, options.selfComparison)
     
    LOG.info('Number of active genes in genomes')
    LOG.info('\n'.join(map(lambda x: '\t%s\t%s'%(x[0], len(x[1])),
        genes.items())))
    if options.minGenomeNo > 1 or options.minContigLen > 1 or \
            options.stringency > 0:
        LOG.info(('Pruning gene similarity graph, min genome occurrence: ' + \
                '%s, min contig len: %s, min stringency: %s') %(options.minGenomeNo, 
                    options.minContigLen, options.stringency))
        pruneGraph(G, gene2genome, options.minContigLen, options.minGenomeNo,
                options.stringency, options.maxIterations)

    gMap = constructGenomeMap(G, fastaFiles, options.outDir)
    LOG.info('Number of active genes')
    LOG.info('\n'.join(map(lambda x: '\t%s\t%s'%(x,
        len(gMap[x][GM_ACTV_GNS_KEY])), gNames)))

    genomeMapOut = open(join(options.outDir, GENOME_MAP_FILE), 'w')
    writeGenomeMap(gMap, gNames, genomeMapOut)
    genomeMapOut.close()

    writePairwiseSimilarities(G, gene2genome, gMap, gNames, options.outDir,
            options.selfComparison)
