configfile: 'config.yaml'

from glob import glob
from os.path import basename, join
from itertools import combinations, repeat, chain

PYEXEC = ' '.join(filter(None, [config['python2_bin'], config['bin']]))
PYSUF = config['pysuf']

GENOMES_DIR = config['genome_data_dir']
GENOMES = sorted(glob('%s/*.fna' %GENOMES_DIR))

MAUVE_CMD = config['mauve_cmd']
MAUVE_OUT = '%s/%s_sw%s' %(config['mauve_out_dir'], basename(MAUVE_CMD),
config['mauve_seed_weight'])

MARKERS_DIR = '%s/%s_marker_l%s' %(config['markers_dir'], basename(MAUVE_OUT),
config['marker_min_length'])
MARKERS_FILES = ['%s/%s.gos' %(MARKERS_DIR, basename(x[:x.rfind('.')])) for x
in GENOMES]

BLAST_DIR = config['blast_dir']
BLAST_PARSTR = config['blast_params'].replace('-', '').replace(' ', '_')
BLAST_OUT = '%s/%s_%s_%s' %(config['blast_out'], basename(MARKERS_DIR), 
config['blast_cmd'], BLAST_PARSTR)

PW_OUT = '%s/%s_pwsim_n2_S_s%s' %(config['pw_out'], basename(BLAST_OUT),
config['pw_stringency'])
PW_SIMS = ['%s/%s_%s.sim' %(PW_OUT, x[0], x[1]) for x in
chain(combinations((basename(y[:y.rfind('.')]) for y in GENOMES), 2),
zip(*repeat([basename(y[:y.rfind('.')]) for y in GENOMES], 2)))]

REF = basename(GENOMES[int(config['psycho_ref'])].rsplit('.', 1)[0])
HIERARCHY_OUT = '%s/%s_psycho_d%s_ref_%s' %(config['psycho_out'],
basename(PW_OUT), config['psycho_delta'], REF)

CIRCOS_CMD = config['circos_cmd']
CIRCOS_PLOT_SUFFIX = config['pwsynteny_plot_params'].find('-s') < 0 and 'main' or 'sub'

#
# MAIN RULE
#

rule all:
    input:
#        '%s/%s_BAPHI_%s.circos.conf' %(HIERARCHY_OUT,  REF,
#        CIRCOS_PLOT_SUFFIX),
#        '%s/hierarchy_d%s.shelve' %(HIERARCHY_OUT,  config['psycho_delta'])
        expand('%s/%s_{target}_%s.png' %(HIERARCHY_OUT,  REF,
        CIRCOS_PLOT_SUFFIX),
        target=[basename(GENOMES[x][:GENOMES[x].rfind('.')]) for x in
        range(len(GENOMES)) if x != config['psycho_ref']])
#        expand('%s/%s_{target}_%s.links' %(HIERARCHY_OUT,  REF,
#        CIRCOS_PLOT_SUFFIX),
#        target=[basename(GENOMES[x][:GENOMES[x].rfind('.')]) for x in
#        range(len(GENOMES)) if x != config['psycho_ref']])

rule run_mauve:
    input:
        GENOMES
    params:
        seed_weight = config['mauve_seed_weight'],
        tmp_dir = config['tmp_dir']
    output:
        xmfa = '%s/alignment.xmfa' %MAUVE_OUT,
        guide_tree = '%s/alignment.tree' %MAUVE_OUT,
        backbone = '%s/%s.backbone' %(MAUVE_OUT,
        '_'.join(basename(x[:x.rfind('.')]) for x in GENOMES))
    log:
        '%s/mauve.log' %MAUVE_OUT
    shell:
        MAUVE_CMD + ' --output={output.xmfa} '
        '--output-guide-tree={output.guide_tree} '
        '--backbone-output={output.backbone} '
        '--seed-weight={params.seed_weight} '
        '--scratch-path-1={params.tmp_dir} {input} > {log}' 

rule mauve_to_markers:
    input:
        backbone = '%s/%s.backbone' %(MAUVE_OUT, '_'.join(
        basename(x[:x.rfind('.')]) for x in GENOMES)),
        genomes = GENOMES
    params:
        min_length = config['marker_min_length'],
        out_dir = MARKERS_DIR
    output:
        MARKERS_FILES 
    shell:
        PYEXEC + 'mauve_backbone_to_multifasta' + PYSUF + 
        ' -l {params.min_length} -o {params.out_dir} '
        '{input.backbone} {input.genomes}'

rule create_blast_db:
    input:
        MARKERS_FILES
    params:
        dbtype = config['blast_db_type'],
        title = 'BLAST database of marker sequences ' + 
        ' '.join(MARKERS_FILES),
        dbname = '%s/%s' %(config['blast_out'], config['blast_db_name'])
    output:
        expand('%s/%s.{dbfile}' %(config['blast_out'], config['blast_db_name']), 
        dbfile=config['blast_db_endings'])
    log:
        '%s/%s.log' %(config['blast_out'], config['blast_db_name'])
    shell:
        join(BLAST_DIR, config['mkblastdb_cmd']) + ' -in \"{input}\" -hash_index '
        '-out {params.dbname} -dbtype {params.dbtype} -title '
        '\"{params.title}\" -logfile {log}'

rule run_blast:
    input:
        markers_file = MARKERS_DIR + '/{markers}.gos',
        blast_db = expand('%s/%s.{dbfile}' %(config['blast_out'],
        config['blast_db_name']), dbfile=config['blast_db_endings'])
    params:
        dbname = '%s/%s' %(config['blast_out'], config['blast_db_name']),
        blast_params = config['blast_params']
    output:
        BLAST_OUT + '/{markers}.gos.blasttbl'
    log:
        BLAST_OUT + '/blastn.log'
    threads: 
        8 
    shell:
        join(BLAST_DIR, config['blast_cmd']) + ' -db {params.dbname} -outfmt 6'
        ' -num_threads ' + '{threads} {params.blast_params} < '
        '{input.markers_file} > {output} 2> {log}'

rule run_pairwise_similarities:
    input:
        expand('%s/{markers}.gos.blasttbl' %BLAST_OUT,
        markers=[basename(x[:x.rfind('.')]) for x in GENOMES])
    params:
        stringency = config['pw_stringency'],
        gos_dir = MARKERS_DIR,
        out_dir = PW_OUT
    output:
        pw_sims = PW_SIMS, 
        genome_map = '%s/genome_map.cfg' %PW_OUT
    log:
        'pairwise_similarities.log'
    shell:
        PYEXEC + 'pairwise_similarities' + PYSUF + ' -n2 -S -s '
        '{params.stringency} -f {params.gos_dir} -o {params.out_dir} {input}'

rule run_psycho:
    input:
        PW_SIMS
    params:
        delta = config['psycho_delta'],
        reference = REF
    threads:
        64
    output:
        '%s/hierarchy_d%s.shelve' %(HIERARCHY_OUT, config['psycho_delta'])
    log:
        'psycho_d%s.log' %config['psycho_delta']
    shell:
        PYEXEC + 'psycho' + PYSUF + ' -d {params.delta} -r '
        '{params.reference} {input} > {output}'

rule create_karyotypes:
    input:
        genome = '%s/{genome}.fna' %GENOMES_DIR,
        markers = '%s/{genome}.gos' %MARKERS_DIR
    output:
        '%s/karyotype.{genome}.txt' %HIERARCHY_OUT
    log:
        '%s/karyotype.{genome}.log' %HIERARCHY_OUT
    shell:
        PYEXEC + 'syn2circos.karyotype' + PYSUF + ' {input.genome} '
        '{input.markers} > {output} 2> {log}'

rule generate_circos_files:
    input:
        shelve = '%s/hierarchy_d%s.shelve' %(HIERARCHY_OUT,
        config['psycho_delta']),
        markers = '%s/{genome}.gos' %MARKERS_DIR
    params:
        karyotype_dir = HIERARCHY_OUT,
        plot_params = config['pwsynteny_plot_params'],
        out_dir = HIERARCHY_OUT
    output:
        '%s/%s_{genome}_%s.circos.conf' %(HIERARCHY_OUT, REF,
        CIRCOS_PLOT_SUFFIX),
        '%s/%s_{genome}_%s.links' %(HIERARCHY_OUT, REF, CIRCOS_PLOT_SUFFIX),
    run:
        shell(PYEXEC + 'inctree2pwsynteny' + PYSUF + ' -k {params.karyotype_dir} ' 
        '-o {params.out_dir} {params.plot_params} {input.shelve} ' + 
        basename(input.markers).rsplit('.')[0])

rule run_circos:
    input:
        circos_conf = '%s/%s_{genome}_%s.circos.conf' %(HIERARCHY_OUT, REF,
        CIRCOS_PLOT_SUFFIX),
        links = '%s/%s_{genome}_%s.links' %(HIERARCHY_OUT, REF,
        CIRCOS_PLOT_SUFFIX),
        ref_karyotype = '%s/karyotype.%s.txt' %(HIERARCHY_OUT, REF),
        target_karyotype = '%s/karyotype.{genome}.txt' %HIERARCHY_OUT 
    params:
        out_dir = HIERARCHY_OUT
    output:
        '%s/%s_{genome}_%s.png' %(HIERARCHY_OUT,  REF, CIRCOS_PLOT_SUFFIX)
    log:
        '%s/%s_{genome}_%s.log' %(HIERARCHY_OUT,  REF, CIRCOS_PLOT_SUFFIX)
    shell:
        CIRCOS_CMD + ' -conf {input.circos_conf} -outputdir {params.out_dir} '
        '2> {log}'
        
