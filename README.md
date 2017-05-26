CONTENTS OF THIS FILE
------------------------------------------------------------------------------

* Requirements
* Installation
* Tutorial


REQUIREMENTS
------------------------------------------------------------------------------

* PYTHON3 and PYTHON2 (https://www.python.org/)
* PYTHON2 libraries:
    - Biopython
    - NetworkX
* SNAKEMAKE (http://snakemake.bitbucket.org)
* BLAST+ (https://www.ncbi.nlm.nih.gov/blast)
* MAUVE (http://darlinglab.org/mauve/)
* CIRCOS (http://www.circos.ca/)


INSTALLATION
------------------------------------------------------------------------------

1. Clone git repository
    $ git clone http://github.com/danydoerr/PSyCHO
2. Adapt configuration file config.yaml:
    * set blast_dir to the directory containing the 'blastn' and 'makeblastdb'
      binary
    * set mauve_cmd to point to the 'progressiveMauve' binary 
    * set circos_cmd to point to the 'circos' binary 


TUTORIAL
------------------------------------------------------------------------------

1. Prepare your genome sequences such that:
    * All genome sequences of each genome are organized in one multi-FASTA file
    * Sequence records have unique IDs within each FASTA file
    * The filename of each FASTA file must match the regular expression
      '[^_]+.fna', i.e. it _must not_ contain underscores and end with the '.fna'
      suffix
    * Move your genome data into the genomes/ subdirectory

2.  Adjust parameters of the tools in the config.yaml file:
    * blast_params: BLAST parameters for the inference of homology between markers:
    * mauve_seed_weight: minimum alignment score that a seed is required to
      have
    * marker_min_length: minimum length of segment to be considered as genomic
      marker in bp
    * pw_stringency: stringency parameter to discard weak sequence similarities
    * psycho_delta: gap-parameter of delta-teams discovery 
    * psycho_ref: index of reference genome (genomes are alphabetically sorted
      acoording to their file names)
    * pwsynteny_plot_params: parameters of the 'inctree2pwsynteny.py' script
      for pairwise visualization of the synteny hierarchy
3. Run 'snakemake --cores {number of cores you want to allocate for running PSyCHO}
