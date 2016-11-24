#!/usr/bin/env python

from sys import stdout, stderr, exit
from optparse import OptionParser
from os.path import dirname, join
import shelve
import re

from psycho import node
from pairwise_similarities import readGenomeMap, GENOME_MAP_FILE, \
        GM_ACTV_GNS_KEY, PAT_CHR

START_END_PAT = re.compile('^.*\|([0-9]+):([0-9]+)(\||$)')

BREWER_COL = ['blues-%s-seq-3', 'bugn-%s-seq-3', 'bupu-%s-seq-3', 'gnbu-%s-seq-3', \
        'greens-%s-seq-3', 'greys-%s-seq-3', 'oranges-%s-seq-3', 'orrd-%s-seq-3', \
        'pubu-%s-seq-3', 'pubugn-%s-seq-3', 'purd-%s-seq-3', 'purples-%s-seq-3', \
        'rdpu-%s-seq-3', 'reds-%s-seq-3', 'ylgn-%s-seq-3', 'ylgnbu-%s-seq-3',\
        'ylorbr-%s-seq-3', 'ylorrd-%s-seq-3']
BREWER_COL_RANGE = range(3,10)


if __name__ == '__main__':
    usage = '%prog [options] <SHELVE> <TARGET GENOME>'
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

    shObj = shelve.open(args[0], flag='r', protocol=-1)
    genomes = shObj['genomes']
    genomes2id = dict(zip(genomes, xrange(len(genomes))))

    genomeMap = readGenomeMap(open(join(dirname(shObj['orig_pw_dists'][0]),
        GENOME_MAP_FILE)))

    ref = shObj['ref']
    G0 = genomes[ref]
    G1 = args[1]
    if G1 not in genomes:
        print 'ERROR: target genome %s is not contained in shelve. Exiting.' %G1
        exit(1)

    id1 = genomes2id[G1]

    gene_orders = shObj['gene_orders']

    g2pos = dict()
    for x in xrange(len(gene_orders)):
        for i in xrange(len(gene_orders[x])):
            g2pos[(x, gene_orders[x][i])] = i

    recovered_markers = shObj['recovered_markers']

    fa1 = [map(int, START_END_PAT.match(sid).group(1,2)) for sid in
            genomeMap[G0][GM_ACTV_GNS_KEY]] 
    fa2 = [map(int, START_END_PAT.match(sid).group(1,2)) for sid in
            genomeMap[G1][GM_ACTV_GNS_KEY]] 
    
    chr1s = sorted(set(map(lambda x: '%s.%s'%(G0.lower(),\
            PAT_CHR.match(x).group(1).lower()), \
            genomeMap[G0][GM_ACTV_GNS_KEY])))
    chr2s = sorted(set(map(lambda x: '%s.%s'%(G1.lower(),\
            PAT_CHR.match(x).group(1).lower()), \
            genomeMap[G1][GM_ACTV_GNS_KEY])))

    #
    # load inclusion tree
    #
    root = shObj['inclusion_tree']
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

    queue = list(map(lambda x: (x, 1), root.children))
    res = set()
    level = dict()
    while queue:
        u, u_depth = queue.pop()

        u_end = u_start = None
        
        if u.intt[0] in recovered_markers:
            p = g2pos[(ref, u.intt[0])]
            while p > 0 and gene_orders[ref][p] in recovered_markers:
                p -= 1
            u_start = fa1[gene_orders[ref][p]][1]

        if u.intt[1] in recovered_markers:
            p = g2pos[(ref, u.intt[1])]
            while p < len(gene_orders[ref])-1 and gene_orders[ref][p] in recovered_markers:
                p += 1
            u_end = fa1[gene_orders[ref][p]][0]

        if u_start == None:
            u_start = fa1[u.intt[0][1]-1][0]

        if u_end == None: 
            u_end = fa1[u.intt[1][1]-1][1]+1
        
        if u.intt[1][1]-u.intt[0][1] < options.min:
            continue

        hasHit = False

        if u.links:
            for (y, start, end) in u.links:
                if y != id1:
                    continue
                v_end = v_start = None
                if start in recovered_markers:
                    p = g2pos[(y, start)]
                    while p > 0 and gene_orders[y][p] in recovered_markers:
                        p -= 1
                    v_start = fa2[gene_orders[y][p]][1]

                if end in recovered_markers:
                    p = g2pos[(ref, end)]
                    while p < len(gene_orders[y])-1 and gene_orders[y][p] in recovered_markers:
                        p += 1
                    v_end = fa2[gene_orders[y][p]][0]

                if v_start == None:
                    v_start = fa2[start[1]-1][0]

                if v_end == None: 
                    v_end = fa2[end[1]-1][1]+1
                
                if v_end -v_start < options.min:
                    continue



                if end[1]-start[1]+1 >= options.min:
#                    if (v_end-v_start)/(u_end-u_start) >= 2:
#                        import pdb; pdb.set_trace() 
                    link = (G0.lower(), u.intt[0][0].lower(), u_start, u_end,\
                            genomes[y].lower(), start[0].lower(), v_start, \
                            v_end)
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
        for Gx, chrx, startx, endx, _, _, starty, endy, in res:
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

    
