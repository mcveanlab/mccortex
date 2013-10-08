Population De Novo Assembly and Variant Calling
===============================================

Multi-sample de novo assembly and variant calling using de bruijn graphs.
Variant calling with and without a reference genome. Between closely related
samples or highly diverged ones. From bacterial to mammalian genomes. Minimal
configuration. And it's free.

Isaac Turner's experimental rewrite of cortex_var, to handle larger populations
with better genome assembly.

7 October 2013

Build
-----

Compiles with clang and gcc. Tested on Mac OS X and linux. Requires zlib.
The first compile will take a while since the libraries in libs/ need to be
downloaded and compiled.

    cd libs; make; cd ..

To compile for a maximum kmer size of 31:

    make

to compile for a maximum kmer size of 63:

    make MAXK=63

executables appear in the bin directory

Commands
--------

    usage: ctx31 <command> [options] <args>
    version: 0.0

    Command:  build       FASTA/FASTQ/BAM -> binary graph file
              view        view and check a cortex graph file (.ctx)
              clean       clean errors from a graph
              join        combine graphs, filter graph intersections
              supernodes  pull out supernodes
              subgraph    filter a subgraph
              reads       filter reads against a graph
              extend      extend contigs using a graph
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

Example pipelines -here- shortly:
1. One high depth sample
2. Multiple low depth samples

Getting Helps
-------------

Isaac Turner: turner.isaac@gmail.com

Code And Contributing
------------

Issues can be submitted on github. Pull requests welcome. Please add your name
to the AUTHORS file.

Code is organised as:
* libs/       included library code from other projects
* src/common  files that need recompiling based on different max kmer size (MAXK)
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
