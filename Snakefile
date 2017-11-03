configfile: 'config.yaml'

from glob import glob
from os.path import basename, join
from itertools import combinations, repeat, chain

PYEXEC = ' '.join(filter(None, [config['python2_bin'], config['bin']]))
PYSUF = config['pysuf']

GENOMES_DIR = config['genome_data_dir']
GENOMES = sorted(filter(lambda x: basename(x).find('_') < 0, 
        glob('%s/*.fna' %GENOMES_DIR)))

BLAST_DIR = config['blast_dir']

BLAST_PARAMS = ' '.join((config['blast_params'], 
        '-perc_identity %s' %config['sgmtn_alignment_ident']))

BLAST_PARSTR = BLAST_PARAMS.replace('-', '').replace(' ', '_')
ATOMS_OUT = join(config['atoms_out'], '%s_%s' %(config['blast_cmd'],
        BLAST_PARSTR))
CORES = snakemake.get_argument_parser().parse_args().cores or 1 
BLAST_THREADS = int(max(1, CORES/len(GENOMES)))

MARKER_OUT = join(config['marker_dir'], '%s_ag%s_al%s_m%s' %(
        basename(ATOMS_OUT), config['sgmtn_alignment_maxgap'],
        config['sgmtn_alignment_minlen'], config['marker_min_length']))
MARKERS_FILES = [join(MARKER_OUT, basename(x)) for x in GENOMES]
PW_SIMS = [join(MARKER_OUT, '%s_%s.sim' %(x,y)) for x,y in
        chain(combinations(map(lambda z: basename(z).rsplit('.', 1)[0], GENOMES), 2), 
        zip(*repeat([basename(z).rsplit('.', 1)[0] for z in GENOMES], 2)))]
GENOME_MAP_FILE = join(MARKER_OUT, 'genome_map.cfg')

REF = basename(GENOMES[int(config['psycho_ref'])].rsplit('.', 1)[0])
HIERARCHY_OUT = join(config['psycho_out'], '%s_psycho_d%s_ref_%s' %(
        basename(MARKER_OUT), config['psycho_delta'], REF))

CIRCOS_CMD = config['circos_cmd']
CIRCOS_PLOT_SUFFIX = config['pwsynteny_plot_params'].find('-s') < 0 and 'main' \
        or 'sub'

#
# MAIN RULE
#

rule all:
    input:
        expand(join(HIERARCHY_OUT, '%s_{target}_%s.png' %(REF,
                CIRCOS_PLOT_SUFFIX)), target=[basename(
                GENOMES[x][:GENOMES[x].rfind('.')]) for x in
                range(len(GENOMES)) if x != config['psycho_ref']])

rule create_blast_db:
    input:
        GENOMES
    params:
        dbtype = config['blast_db_type'],
        title = 'BLAST database of marker sequences ' + 
        ' '.join(GENOMES),
        dbname = '%s/%s' %(config['atoms_out'], config['blast_db_name'])
    output:
        expand('%s/%s.{dbfile}' %(config['atoms_out'], config['blast_db_name']), 
                dbfile=config['blast_db_endings'])
    log:
        '%s/%s.log' %(config['atoms_out'], config['blast_db_name'])
    shell:
        'mkdir -p "%s";' %config['atoms_out'] + 
        join(BLAST_DIR, config['mkblastdb_cmd']) + ' -in \"{input}\" '
        '-hash_index -out {params.dbname} -dbtype {params.dbtype} -title '
        '\"{params.title}\" -logfile {log}'


rule run_blast:
    input:
        markers_file = GENOMES_DIR + '/{genome}.fna',
        blast_db = expand(join(config['atoms_out'], '%s.{dbfile}' %(
                config['blast_db_name'])), dbfile=config['blast_db_endings'])
    params:
        dbname = join(config['atoms_out'], config['blast_db_name']),
        blast_params = BLAST_PARAMS
    output:
        temp(ATOMS_OUT + '_{genome,[^_]+}.psl')
    log:
        join(config['atoms_out'], 'blastn.log')
    threads: 
        BLAST_THREADS
    shell:
        'mkdir -p "%s";' %config['atoms_out']+
        join(BLAST_DIR, config['blast_cmd']) + ' -db {params.dbname} '
        '-num_threads {threads} {params.blast_params} < {input.markers_file} |' +
        PYEXEC + 'blast2psl' + PYSUF +' > {output} 2> {log}'


rule concat_psl:
    input:
        expand(ATOMS_OUT + '_{genome}.psl', genome=map(lambda x:
                basename(x).rsplit('.', 1)[0], GENOMES))
    output:
        ATOMS_OUT + '.psl' 
    shell:
        'cat {input} > {output};' 


rule atomizer:
    input:
        ATOMS_OUT + '.psl'
    params:
        out_dir = MARKER_OUT,
        min_len = config['marker_min_length'],
        al_ident = config['sgmtn_alignment_ident'],
        al_min_len = config['sgmtn_alignment_minlen'],
        al_max_gap = config['sgmtn_alignment_maxgap']
    output:
        ATOMS_OUT + '_ag%s_al%s_m%s.atoms' %(config['sgmtn_alignment_maxgap'],
            config['sgmtn_alignment_minlen'], config['marker_min_length'])
    log:
        ATOMS_OUT + '_ag%s_al%s_m%s.log' %(config['sgmtn_alignment_maxgap'],
            config['sgmtn_alignment_minlen'], config['marker_min_length'])
    shell:
        config['segmentation_cmd'] + ' {input} --minLength {params.min_len} '
        '--minIdent {params.al_min_len} --minIdent {params.al_ident} --maxGap '
        '{params.al_max_gap} > {output} 2> {log}'

#rule atomizer_orig:
#    input:
#        ATOMS_OUT + '.psl'
#    params:
#        out_dir = MARKER_OUT,
#        min_len = config['marker_min_length'],
#        al_ident = config['sgmtn_alignment_ident'],
#        al_min_len = config['sgmtn_alignment_minlen'],
#        al_max_gap = config['sgmtn_alignment_maxgap']
#    output:
#        ATOMS_OUT + '.atoms'
#    log:
#        'atomizer.log'
#    shell:
#        config['orig_segmentation_cmd'] + '--alg IMP  --fasta {input.fasta} --psl {input.psl} --min-len {params.min_len} '
#        '--max-dist 100 --max-gap {params.al_max_gap} --min-sim {params.al_ident} > {output} 2> {log}'


rule atoms_to_markers:
    input:
        atoms_file = ATOMS_OUT + '_ag%s_al%s_m%s.atoms' %(
            config['sgmtn_alignment_maxgap'], config['sgmtn_alignment_minlen'],
            config['marker_min_length']),
        fasta_files = GENOMES
    params:
        out_dir = MARKER_OUT,
        reps  = config['marker_exclude_duplicates']
    output:
        MARKERS_FILES,
        PW_SIMS,
        GENOME_MAP_FILE
    shell:
        PYEXEC + 'atoms_to_pwsim' + PYSUF + ' -o {params.out_dir} -e '
        '{params.reps} {input.atoms_file} {input.fasta_files}'


rule run_psycho:
    input:
        pw_sims = PW_SIMS,
        genome_map = GENOME_MAP_FILE
    params:
        delta = config['psycho_delta'],
        coverage = config['psycho_cov'],
        reference = REF
    threads:
        64
    output:
        join(HIERARCHY_OUT, 'hierarchy_d%s.json' %config['psycho_delta'])
    log:
        'psycho_d%s.log' %config['psycho_delta']
    shell:
        PYEXEC + 'psycho' + PYSUF + ' -d {params.delta} -r '
        '{params.reference} -c {params.coverage} -o {output} {input.pw_sims}'


rule create_karyotypes:
    input:
        genome = join(GENOMES_DIR, '{genome}.fna'),
        markers = join(MARKER_OUT, '{genome}.fna')
    output:
        join(HIERARCHY_OUT, 'karyotype.{genome}.txt')
    log:
        join(HIERARCHY_OUT, 'karyotype.{genome}.log')
    shell:
        PYEXEC + 'syn2circos.karyotype' + PYSUF + ' {input.genome} '
        '{input.markers} > {output} 2> {log}'


rule generate_indeterminate_regions_histogram:
    input:
        GENOMES
    params:
        window_size = config['pwsynteny_indet_regions_window_size']
    output:
        join(HIERARCHY_OUT, 'indet_regions_hist.txt')
    shell:
        PYEXEC + 'syn2circos.indet_regions' + PYSUF + ' -w {params.window_size} '
        '{input} > {output}'
        

rule generate_circos_files:
    input:
        json = join(HIERARCHY_OUT, 'hierarchy_d%s.json' %config['psycho_delta']),
        markers = join(MARKER_OUT, '{genome}.fna'),
        indet_regions_hist = join(HIERARCHY_OUT, 'indet_regions_hist.txt')
    params:
        karyotype_dir = HIERARCHY_OUT,
        plot_params = config['pwsynteny_plot_params'],
        out_dir = HIERARCHY_OUT
    output:
        join(HIERARCHY_OUT, '%s_{genome}_%s.circos.conf' %(REF,
                CIRCOS_PLOT_SUFFIX)),
        join(HIERARCHY_OUT, '%s_{genome}_%s.links' %(REF, CIRCOS_PLOT_SUFFIX)),
    shell:
        PYEXEC + 'inctree2pwsynteny' + PYSUF + ' -k {params.karyotype_dir} ' 
        '-o {params.out_dir} -i {input.indet_regions_hist} {params.plot_params}'
        ' {input.json} {wildcards.genome}' 


rule run_circos:
    input:
        circos_conf = join(HIERARCHY_OUT, '%s_{genome}_%s.circos.conf' %(REF,
                CIRCOS_PLOT_SUFFIX)),
        links = join(HIERARCHY_OUT, '%s_{genome}_%s.links' %(REF,
                CIRCOS_PLOT_SUFFIX)),
        ref_karyotype = join(HIERARCHY_OUT, 'karyotype.%s.txt' %REF),
        target_karyotype = join(HIERARCHY_OUT, 'karyotype.{genome}.txt')
    params:
        out_dir = HIERARCHY_OUT
    output:
        join(HIERARCHY_OUT, '%s_{genome}_%s.png' %(REF, CIRCOS_PLOT_SUFFIX))
    log:
        join(HIERARCHY_OUT, '%s_{genome}_%s.log' %(REF, CIRCOS_PLOT_SUFFIX))
    shell:
        CIRCOS_CMD + ' -conf {input.circos_conf} -outputdir {params.out_dir} '
        '2> {log}'
        
