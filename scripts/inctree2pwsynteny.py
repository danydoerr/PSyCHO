#!/usr/bin/env python

from sys import stdout, stderr, exit
from optparse import OptionParser
from os.path import dirname, join, isabs
from itertools import izip, chain
import json
import re

from psycho import dict2hierarchy, CONTIG_BOUNDARY_KEY
from pairwise_similarities import readGenomeMap, GENOME_MAP_FILE, \
        GM_ACTV_GNS_KEY, PAT_CHR, PAT_POS

BREWER_COL = ['blues-%s-seq-3', 'bugn-%s-seq-3', 'bupu-%s-seq-3', 'gnbu-%s-seq-3', \
        'greens-%s-seq-3', 'greys-%s-seq-3', 'oranges-%s-seq-3', 'orrd-%s-seq-3', \
        'pubu-%s-seq-3', 'pubugn-%s-seq-3', 'purd-%s-seq-3', 'purples-%s-seq-3', \
        'rdpu-%s-seq-3', 'reds-%s-seq-3', 'ylgn-%s-seq-3', 'ylgnbu-%s-seq-3',\
        'ylorbr-%s-seq-3', 'ylorrd-%s-seq-3']
BREWER_COL_RANGE = range(3,10)


def getSims(ref_int, y_bounds):
    res_x = list()
    res_y = list()

    x_int = [gos[ref][x] for x in xrange(gos2pos[ref][ref_int[0]],
        gos2pos[ref][ref_int[1]]+1)]
    y_int = [gos[y_bounds[0]][x] for x in xrange(gos2pos[y_bounds[0]][y_bounds[1]],
        gos2pos[y_bounds[0]][y_bounds[2]]+1)]

    for gx in x_int:
        res_x.append(list())
        for gy in dists[(G0, genomes[y_bounds[0]])][gx]:
            if gy in y_int:
                res_x[-1].append(gy[1])

    for gy in y_int:
        res_y.append(list())
        for gx in dists[(genomes[y_bounds[0]], G0)][gy]:
            if gx in x_int:
                res_y[-1].append(gx[1])
    return res_x, res_y

if __name__ == '__main__':
    usage = '%prog [options] <PSYCHO JSON OUTPUT> <TARGET GENOME>'
    parser = OptionParser(usage=usage)
    parser.add_option('-c', '--color_links', dest='colorLinks', default=False,
            action='store_true', help='Color links individually')
    parser.add_option('-m', '--min', dest='min', default=3, type='int',
            help='minimum number of positions a wci pair [default: %default]')
    parser.add_option('-s', '--plot_subintervals', dest='subInts', default=False,
            action='store_true', help='Not only draw first level intervals, ' + \
                    'but also subintervals [default: %default]')
    parser.add_option('-l', '--plot_level', dest='level', default=-1,
            type='int', help='Plot a specific level [default: %default]')
    parser.add_option('-k', '--karyotype_dir', dest='karyotypeDir',
            default='.', type='str', help='Directory ' + \
                    'in which the karyotype information is stored ' + \
                    '[default: %default]')
    parser.add_option('-o', '--output_dir', dest='outDir', type=str,
            default='.', help='Output directory [default=%default]')

    (options, args) = parser.parse_args()

    # check if input file is given, if so, call function with data
    if len(args) != 2:
        parser.print_help()
        exit(1)

    if options.level != -1 and options.subInts:
        print >> stderr, 'Options --plot_level and --plot_subintervals ' + \
                'are incompatabile with each other. Exiting'
        parser.print_help()
        exit(1)

    jsDict = json.load(open(args[0])) 
    genomes = jsDict['genome_names']
    genomes2id = dict(zip(genomes, xrange(len(genomes))))

#    # XXX only for debugging purposes
#    from pairwise_similarities import readDists, reverseDistMap
#    from os.path import basename
#
#    dists = dict((tuple(basename(x).split('.')[0].split('_')), \
#            readDists(open(isabs(x) and x or join(dirname(args[0]), x)))[1]) for \
#            x in jsDict['orig_pw_dists'])
#    dists.update(((k[1], k[0]), reverseDistMap(v)) for k,v in dists.items())
#    
#    dpath = dirname(jsDict['orig_pw_dists'][0])
#    genomeMap = readGenomeMap(open(join(isabs(dpath) and dpath or
#        join(dirname(args[0]), dpath), GENOME_MAP_FILE)))

    ref = jsDict['ref_id']
    G0 = genomes[ref]
    G1 = args[1]
    if G1 not in genomes:
        print 'ERROR: target genome %s is not contained in shelve. Exiting.' %G1
        exit(1)

    id1 = genomes2id[G1]

    marker_seq_list = jsDict['marker_seq_list']
    g2pos = [[dict(izip(go, xrange(len(go)))) for go in gos] for gos in \
            marker_seq_list]
    recovered_markers = jsDict['recovered_markers']

    chr1s = sorted(set('%s.%s'%(G0.lower(), PAT_CHR.match(x).group(1).lower()) \
            for x in chain(*map(lambda m: m[ref], marker_seq_list)) if x != \
            CONTIG_BOUNDARY_KEY))
    chr2s = sorted(set('%s.%s'%(G1.lower(), PAT_CHR.match(x).group(1).lower()) \
            for x in chain(*map(lambda m: m[id1], marker_seq_list)) if x != \
            CONTIG_BOUNDARY_KEY))

    #
    # load inclusion tree
    #
    root = dict2hierarchy(jsDict['sb_hierarchy'])
#    print >> stderr, '++ loading inclusion tree of genome %s from shelve' %G0
    
    suffix = 'main'
    if options.subInts:
        suffix = 'sub'
    elif options.level != -1:
        suffix = 'level_%s' %options.level

    #
    # write links
    #
    link_out = open(join(options.outDir, '%s_%s_%s.links' %(G0, G1, suffix)),
            'w')

    queue = [(root, 1)]
    res = set()
    level = dict()
    while queue:
        u, u_depth = queue.pop()

        hasHit = False
        if u.id != None:
            d1 = u.intt[1]-u.intt[0]+1
            if d1 < options.min:
                continue

            gos = marker_seq_list[u.id]
            gos2pos = g2pos[u.id]

            u_start, u_end = u.intt
            while u_start > 0 and gos[ref][u_start] in recovered_markers[ref]:
                u_start -= 1
            while u_end < len(gos[ref])-1 and gos[ref][u_end] in \
                    recovered_markers[ref]:
                u_end += 1

            # if the boundary of the syntenic block has been extended due to
            # recovered markers, take the inner extremity of the neighboring
            # marker, otherwise the outer extremity of the contained marker
            u_start_p = int(PAT_POS.match(gos[ref][u_start]).group(u_start
                == u.intt[0] and 1 or 2))
            u_end_p   = int(PAT_POS.match(gos[ref][u_end]).group(u_end ==
                u.intt[1] and 2 or 1))
            u_chr = PAT_CHR.match(gos[ref][u_start]).group(1)

            if u.links and (not u.parent or u.parent.links != u.links):
                for i in xrange(len(u.links)):
                    y, start, end = u.links[i]

                    d2 = end-start+1
                    if y != id1 or d2 < options.min:
                        continue

                    while start > 0 and gos[y][start] in recovered_markers[y]:
                        start -= 1
                    while end < len(gos[y])-1 and gos[y][end] in \
                            recovered_markers[y]:
                        end += 1

                    # if the boundary of the syntenic block has been extended
                    # due to recovered markers, take the inner extremity of the
                    # neighboring marker, otherwise the outer extremity of the
                    # contained marker
                    v_start_p = int(PAT_POS.match(gos[y][start]).group(\
                            start == u.links[i][1] and 1 or 2))
                    v_end_p = int(PAT_POS.match(gos[y][end]).group(end \
                            == u.links[2] and 2 or 1))
                    v_chr = PAT_CHR.match(gos[y][start]).group(1)
#                   if max(d1,d2)/float(min(d1,d2)) >= 3:
#                       sims = getSims(u.intt, (y, start, end))
#                       import pdb; pdb.set_trace() 
                    link = (G0.lower(), u_chr.lower(), u_start_p, u_end_p,\
                            genomes[y].lower(), v_chr.lower(), v_start_p, \
                            v_end_p)
                    res.add(link)
                    level[link] = level.get(link, u_depth)
                    hasHit = True
        if not hasHit or options.subInts or options.level > u_depth:
            queue.extend(map(lambda y: (y, min(hasHit and u_depth+1 or u_depth,
                5)), filter(lambda x: len(x.children) > 1, u.children)))

    res = sorted(res)
    for r in res:
        if options.level == -1 or level[r] == options.level:
            print >> link_out, '%s.%s %s %s %s.%s %s %s level=%s'%(r + (5-level[r], ))

    link_out.close()

    circos_out = open(join(options.outDir, '%s_%s_%s.circos.conf' %(G0, G1,
        suffix)), 'w')

    print >> circos_out, ('karyotype = %s/karyotype.%s.txt,%s/' + \
            'karyotype.%s.txt') %(options.karyotypeDir, G0,
                    options.karyotypeDir, G1)
    print >> circos_out, 'chromosomes_units = 1000000'
    print >> circos_out, 'chromosomes_display_default = no' 
    print >> circos_out, 'chromosomes                 = %s' %(';'.join(chr1s +
        chr2s))

    print >> circos_out, '<colors>'
    for c in chr1s:
        print >> circos_out, '%s = blue' %c
    for c in chr2s:
        print >> circos_out, '%s = green' %c

    print >> circos_out, '%s.seg = black' %G0.lower()
    print >> circos_out, '%s.no_seg = blue' %G0.lower()
    print >> circos_out, '%s.seg = black' %G1.lower()
    print >> circos_out, '%s.no_seg = green' %G1.lower()

    c_max = len(BREWER_COL_RANGE) * len(BREWER_COL)
    if options.colorLinks:
        x = 0
        for Gx, chrx, startx, endx, _, _, starty, endy, in set(res):
            r = (x % c_max) / len(BREWER_COL)
            i = (x % c_max) % len(BREWER_COL)
            print >> circos_out, '%s.%s.%s.%s.%s.%s.l = %s' %(Gx, chrx, startx, endx,
                    starty, endy, BREWER_COL[i] % BREWER_COL_RANGE[r])
            x += 1
    else:
        c_min = len(chr1s) < len(chr2s) and chr1s or chr2s
        if len(c_min) > c_max:
            print >> stderr, '!! Too many chromosomes. Not enough colors ' + \
                    'available to color all links! Exiting'
            exit(1)
        x = 0
        for c in c_min:
            r = x / len(BREWER_COL)
            i = x % len(BREWER_COL)
            print >> circos_out, '%s.l = %s' %(c, BREWER_COL[i] %BREWER_COL_RANGE[r])
            x += 1

    print >> circos_out, '</colors>'
    print >> circos_out, '<links>\n<link>'
    print >> circos_out, 'file          = %s' %link_out.name
    print >> circos_out, 'radius        = 0.99r'
    print >> circos_out, 'bezier_radius = 0r'
    print >> circos_out, 'ribbon = yes'
    print >> circos_out, 'flat   = yes'
    print >> circos_out, '<rules>\n<rule>'
    print >> circos_out, 'condition     = 1'
    if options.colorLinks: 
        print >> circos_out, 'color         = eval(var(chr1) . \'.\' .' + \
                'var(start1) . \'.\' . var(end1) . \'.\' . var(start2)' + \
                ' . \'.\' . var(end2) . \'.l_a\' . var(level))'
    else:
        print >> circos_out, 'color         = eval(var(chr%s) . \'.l_a\' . var(level))' %(c_min ==
                chr1s and 1 or 2)
    print >> circos_out, 'flow          = continue'
    print >> circos_out, '</rule>\n</rules>'
    print >> circos_out, '</link>\n</links>'
    print >> circos_out, '<ideogram>'
    print >> circos_out, '<spacing>'
    print >> circos_out, 'default = 0.005r'
    print >> circos_out, '</spacing>'
    print >> circos_out, 'radius           = 0.92r'
    print >> circos_out, 'thickness        = 100p'
    print >> circos_out, 'stroke_thickness = 0'
    print >> circos_out, 'show_bands            = yes'
    print >> circos_out, 'fill_bands            = yes'
    print >> circos_out, 'band_stroke_thickness = 0'
    print >> circos_out, 'band_transparency     = 1'
    print >> circos_out, 'show_label       = yes'
    print >> circos_out, 'label_font       = default'
    print >> circos_out, 'label_radius     = 1.045r '
    print >> circos_out, 'label_size       = 60'
    print >> circos_out, 'label_parallel   = yes'
    print >> circos_out, '</ideogram>'
    print >> circos_out, 'show_ticks          = yes'
    print >> circos_out, 'show_tick_labels    = yes'
    print >> circos_out, '<ticks>'
    print >> circos_out, 'radius           = 1r'
    print >> circos_out, 'color            = black'
    print >> circos_out, 'thickness        = 2p'
    print >> circos_out, 'multiplier       = 1e-6'
    print >> circos_out, 'format           = %d'
    print >> circos_out, '<tick>'
    print >> circos_out, 'spacing        = 5u'
    print >> circos_out, 'size           = 10p'
    print >> circos_out, '</tick>'
    print >> circos_out, '<tick>'
    print >> circos_out, 'spacing        = 25u'
    print >> circos_out, 'size           = 15p'
    print >> circos_out, 'show_label     = yes'
    print >> circos_out, 'label_size     = 20p'
    print >> circos_out, 'label_offset   = 10p'
    print >> circos_out, 'format         = %d'
    print >> circos_out, '</tick>'
    print >> circos_out, '</ticks>'
    print >> circos_out, '<image>\nfile  = %s_%s_%s.png' %(G0, G1, suffix)
    print >> circos_out, 'dir   = .\npng   = yes\nradius         = 1500p' 
    print >> circos_out, 'angle_offset      = -90'
    print >> circos_out, 'auto_alpha_colors = yes\nauto_alpha_steps  = 5'
    print >> circos_out, 'background = white\n</image>'
    print >> circos_out, '<<include etc/colors_fonts_patterns.conf>>'
    print >> circos_out, '<<include etc/colors.brewer.conf>>'
    print >> circos_out, '<<include etc/housekeeping.conf>>'
    print >> circos_out, 'data_out_of_range* = warning'

    circos_out.close()

    
