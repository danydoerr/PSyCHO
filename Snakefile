configfile: 'config.yaml'

from glob import glob
from os.path import basename, join
from itertools import combinations, repeat, chain

PYEXEC = ' '.join(filter(None, [config['python2_bin'], config['bin']]))
PYSUF = config['pysuf']

GENOMES_DIR = config['genome_data_dir']
GENOMES = sorted(glob('%s/*.fna' %GENOMES_DIR))

BLAST_DIR = config['blast_dir']
BLAST_PARSTR = config['blast_params'].replace('-', '').replace(' ', '_')
BLAST_OUT = '%s/%s_%s' %(config['blast_out'], config['blast_cmd'], BLAST_PARSTR)
CORES = snakemake.get_argument_parser().parse_args().cores or 1 
BLAST_THREADS = int(max(1, CORES/len(GENOMES)))

SEGMENTATION_OUT = '%s/%s_ai%s_ag%s_al%s_m%s' %(config['segmentation_out'],
        basename(BLAST_OUT), config['sgmtn_alignment_ident'],
        config['sgmtn_alignment_maxgap'], config['sgmtn_alignment_minlen'],
        config['marker_min_length'])
#MARKERS_DIR = '%s/%s_marker_l%s' %(config['markers_dir'], basename(MAUVE_OUT),
#        config['marker_min_length'])
MARKERS_FILES = ['%s/%s' %(SEGMENTATION_OUT, basename(x)) for x in GENOMES]

REF = basename(GENOMES[int(config['psycho_ref'])].rsplit('.', 1)[0])
HIERARCHY_OUT = '%s/%s_psycho_d%s_ref_%s' %(config['psycho_out'],
        basename(SEGMENTATION_OUT), config['psycho_delta'], REF)

CIRCOS_CMD = config['circos_cmd']
CIRCOS_PLOT_SUFFIX = config['pwsynteny_plot_params'].find('-s') < 0 and 'main' or 'sub'

#
# MAIN RULE
#

rule all:
    input:
        MARKERS_FILES
#        expand('%s/%s_{target}_%s.png' %(HIERARCHY_OUT,  REF,
#                CIRCOS_PLOT_SUFFIX), target=[basename(
#                GENOMES[x][:GENOMES[x].rfind('.')]) for x in
#                range(len(GENOMES)) if x != config['psycho_ref']])

rule create_blast_db:
    input:
        GENOMES
    params:
        dbtype = config['blast_db_type'],
        title = 'BLAST database of marker sequences ' + 
        ' '.join(GENOMES),
        dbname = '%s/%s' %(config['blast_out'], config['blast_db_name'])
    output:
        expand('%s/%s.{dbfile}' %(config['blast_out'], config['blast_db_name']), 
        dbfile=config['blast_db_endings'])
    log:
        '%s/%s.log' %(config['blast_out'], config['blast_db_name'])
    shell:
        'mkdir -p "%s";' %config['blast_out'] + 
        join(BLAST_DIR, config['mkblastdb_cmd']) + ' -in \"{input}\" -hash_index '
        '-out {params.dbname} -dbtype {params.dbtype} -title '
        '\"{params.title}\" -logfile {log}'

rule run_blast:
    input:
        markers_file = GENOMES_DIR + '/{genome}.fna',
        blast_db = expand('%s/%s.{dbfile}' %(config['blast_out'],
        config['blast_db_name']), dbfile=config['blast_db_endings'])
    params:
        dbname = '%s/%s' %(config['blast_out'], config['blast_db_name']),
        blast_params = config['blast_params']
    output:
        temp(BLAST_OUT + '/{genome}.psl')
    log:
        BLAST_OUT + '/blastn.log'
    threads: 
        BLAST_THREADS
    shell:
        'mkdir -p "%s";' %BLAST_OUT +
        join(BLAST_DIR, config['blast_cmd']) + ' -db {params.dbname} '
        '-num_threads {threads} {params.blast_params} < {input.markers_file} |'
        + PYEXEC + 'blast2psl' + PYSUF + ' > {output} 2> {log}'

rule concat_psl:
    input:
        expand(BLAST_OUT + '/{genome}.psl', genome=map(lambda x:
                basename(x).rsplit('.', 1)[0], GENOMES))
    output:
        '%s/all.psl' %BLAST_OUT 
    shell:
        'cat {input} > {output};' 

rule atomizer:
    input:
        '%s/all.psl' %BLAST_OUT,
    params:
        out_dir = SEGMENTATION_OUT,
        min_len = config['marker_min_length'],
        al_ident = config['sgmtn_alignment_ident'],
        al_min_len = config['sgmtn_alignment_minlen'],
        al_max_gap = config['sgmtn_alignment_maxgap']
    output:
        '%s/all.atoms' %SEGMENTATION_OUT
    shell:
        config['segmentation_cmd'] + ' {input} --minLength {params.min_len} '
        '--minIdent {params.al_min_len} --minIdent {params.al_ident} --maxGap '
        '{params.al_max_gap} > {output}'

rule atoms_to_markers:
    input:
        atoms_file = '%s/all.atoms' %SEGMENTATION_OUT,
        fasta_files = GENOMES
    params:
        out_dir = SEGMENTATION_OUT
    output:
        MARKERS_FILES
    shell:
        PYEXEC + 'atoms_to_dna' + PYSUF + ' -o {params.out_dir} '
        '{input.atoms_file} {input.fasta_files}'

#rule run_psycho:
#    input:
#        PW_SIMS
#    params:
#        delta = config['psycho_delta'],
#        coverage = config['psycho_cov'],
#        reference = REF
#    threads:
#        64
#    output:
#        '%s/hierarchy_d%s.json' %(HIERARCHY_OUT, config['psycho_delta'])
#    log:
#        'psycho_d%s.log' %config['psycho_delta']
#    shell:
#        PYEXEC + 'psycho' + PYSUF + ' -d {params.delta} -r '
#        '{params.reference} {input} -c {params.coverage} -o {output}'
#
#rule create_karyotypes:
#    input:
#        genome = '%s/{genome}.fna' %GENOMES_DIR,
#        markers = '%s/{genome}.gos' %SEGMENTATION_OUT
#    output:
#        '%s/karyotype.{genome}.txt' %HIERARCHY_OUT
#    log:
#        '%s/karyotype.{genome}.log' %HIERARCHY_OUT
#    shell:
#        PYEXEC + 'syn2circos.karyotype' + PYSUF + ' {input.genome} '
#        '{input.markers} > {output} 2> {log}'
#
#rule generate_circos_files:
#    input:
#        json = '%s/hierarchy_d%s.json' %(HIERARCHY_OUT,
#        config['psycho_delta']),
#        markers = '%s/{genome}.gos' %SEGMENTATION_OUT
#    params:
#        karyotype_dir = HIERARCHY_OUT,
#        plot_params = config['pwsynteny_plot_params'],
#        out_dir = HIERARCHY_OUT
#    output:
#        '%s/%s_{genome}_%s.circos.conf' %(HIERARCHY_OUT, REF,
#        CIRCOS_PLOT_SUFFIX),
#        '%s/%s_{genome}_%s.links' %(HIERARCHY_OUT, REF, CIRCOS_PLOT_SUFFIX),
#    run:
#        shell(PYEXEC + 'inctree2pwsynteny' + PYSUF + ' -k {params.karyotype_dir} ' 
#        '-o {params.out_dir} {params.plot_params} {input.json} ' + 
#        basename(input.markers).rsplit('.')[0])
#
#rule run_circos:
#    input:
#        circos_conf = '%s/%s_{genome}_%s.circos.conf' %(HIERARCHY_OUT, REF,
#        CIRCOS_PLOT_SUFFIX),
#        links = '%s/%s_{genome}_%s.links' %(HIERARCHY_OUT, REF,
#        CIRCOS_PLOT_SUFFIX),
#        ref_karyotype = '%s/karyotype.%s.txt' %(HIERARCHY_OUT, REF),
#        target_karyotype = '%s/karyotype.{genome}.txt' %HIERARCHY_OUT 
#    params:
#        out_dir = HIERARCHY_OUT
#    output:
#        '%s/%s_{genome}_%s.png' %(HIERARCHY_OUT,  REF, CIRCOS_PLOT_SUFFIX)
#    log:
#        '%s/%s_{genome}_%s.log' %(HIERARCHY_OUT,  REF, CIRCOS_PLOT_SUFFIX)
#    shell:
#        CIRCOS_CMD + ' -conf {input.circos_conf} -outputdir {params.out_dir} '
#        '2> {log}'
        
