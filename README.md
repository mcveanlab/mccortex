Population De Novo Assembly and Variant Calling
===============================================

Multi-sample de novo assembly and variant calling using de bruijn graphs.
Variant calling with and without a reference genome. Between closely related
samples or highly diverged ones. From bacterial to mammalian genomes. Minimal
configuration. And it's free.

Isaac Turner's experimental rewrite of cortex_var, to handle larger populations
with better genome assembly.

26 October 2013

Build
-----

Compiles with clang and gcc. Tested on Mac OS X and linux. Requires zlib.
The first compile will take a while since the libraries in libs/ need to be
downloaded and compiled.

To compile for a maximum kmer size of 31:

    make

to compile for a maximum kmer size of 63:

    make MAXK=63

executables appear in the bin/ directory

Commands
--------

    usage: ctx31 <command> [options] <args>
    version: 0.0; zlib: 1.2.5

    Command:  build       FASTA/FASTQ/BAM -> binary graph file
              view        view and check a cortex graph file (.ctx)
              clean       clean errors from a graph
              join        combine graphs, filter graph intersections
              supernodes  pull out supernodes
              subgraph    filter a subgraph
              reads       filter reads against a graph
              extend      extend contigs using a graph
              contigs     pull out contigs for a sample
              inferedges  infer edges before calling `thread`
              thread      thread reads through cleaned population
              pview       view read threading information
              pmerge      merge path files (.ctp)
              call        call variants
              unique      remove duplicated bubbles, produce VCF
              place       place variants and genotype

      Type a command with no arguments to see help.

    Common Options:
      -m <M>       Memory e.g. 1GB [default: 1GB]
      -h <H>       Hash entries [default: 4M (~4 million)]
      -g <G>       Species genome size [default: 3.1Gbp]
      -t <T>       Number of threads [default: 2]
      -k <K>       Kmer size [default: read from binaries]
      -f <file>    Input file
      -p <in.ctp>  Assembly file


Inputs
------

Support reading FASTA, FASTQ, SAM and BAM files. Also supports gzipped files.
File formats are auto-detected.

Cortex graph files can be loaded by specifying a subset of colours:

    in.ctx:0,3-6

will load colours 0, 3, 4, 5 and 6 from graph file in.ctx

Common Pipelines
----------------

Construct graphs for three samples (using 70GB ram):

    ctx31 build -m 70G -k 31 --sample NA12878 --seq input.fq.gz NA12878.ctx
    ctx31 build -m 70G -k 31 --sample Mickey --seq2 reads.1.fq.gz reads.2.fq.gz Mickey.ctx
    ctx31 build -m 70G -k 31 --sample Minnie --seq data.bam Minnie.ctx

Construct graph for the reference (hg19)

    ctx31 build -m 70G -k 31 --sample hg19 --seq hg19.fa.gz hg19.ctx

A) 'Clean' graphs to remove sequencing error (per sample, for high coverage samples)

    ctx31 clean NA12878.clean.ctx NA12878.ctx
    ctx31 clean Mickey.clean.ctx Mickey.ctx
    ctx31 clean Minnie.clean.ctx Minnie.ctx

...then merge graphs into file `refAndSamples.ctx`. Uses 80GB ram:

    ctx31 join -m 80G refAndSamples.clean.ctx hg19.ctx NA12878.clean.ctx Mickey.clean.ctx Minnie.clean.ctx

B) Alternatively merge uncleaned samples then clean on the population (multiple low depth samples). Uses 80GB ram:

    ctx31 join -m 80G refAndSamples.ctx hg19.ctx NA12878.clean.ctx Mickey.clean.ctx Minnie.clean.ctx
    ctx31 clean refAndSamples.clean refAndSamples.ctx

Now we have cortex graphs of the reference and our samples with sequencing error removed.

Coming soon:

1. ctx31 thread
2. ctx31 call
3. ctx31 unique
4. ctx31 place

Getting Helps
-------------

Isaac Turner: turner.isaac@gmail.com

Code And Contributing
------------

Issues can be submitted on github. Pull requests welcome. Please add your name
to the AUTHORS file.

Code is organised as:
* libs/       included library code from other projects
* src/kmer    files that need recompiling based on different max kmer size (MAXK)
* src/basic   files that do not depend on MAXK
* src/tools   tools to do particular jobs
* src/main    files with a main function go in here

Code should compile on mac/linux with clang/gcc without errors or warnings.

License: GPLv2
--------------

Bundled libraries may have different licenses:
* GNU Science Library (GPL)
* seqan (BSD/3-clause)
* city_hash (MIT)
* lookup3 (Public Domain)
* htslib
* seq_file (BSD)
* string_buffer (BSD)
* bit_array (BSD)
* seq-align (GPLv3)

Used in testing:
* bioinf-perl

Citing
------

'Cortex with low memory and read threading' is currently unpublished.  Please
cite previous cortex_var papers:

* De novo assembly and genotyping of variants using colored de Bruijn graphs,
Iqbal(*), Caccamo(*), Turner, Flicek, McVean (Nature Genetics) (2012)
(doi:10.1038/ng.1028)
* High-throughput microbial population genomics using the Cortex variation assembler,
Iqbal, Turner, McVean (Bioinformatics) (Nov 2012)
(doi:10.1093/bioinformatics/bts673)
