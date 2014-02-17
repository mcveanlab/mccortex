Population De Novo Assembly and Variant Calling
===============================================

Multi-sample de novo assembly and variant calling using de bruijn graphs.
Variant calling with and without a reference genome. Between closely related
samples or highly diverged ones. From bacterial to mammalian genomes. Minimal
configuration. And it's free.

Isaac Turner's experimental rewrite of cortex_var, to handle larger populations
with better genome assembly.

14 February 2014

[![Build Status](https://magnum.travis-ci.com/noporpoise/ninja-cortex.png?token=HkeonfUv1FrRw6UkpJot&branch=master)](https://magnum.travis-ci.com/noporpoise/ninja-cortex)

Build
-----

Compiles with clang and gcc. Tested on Mac OS X and linux. Requires zlib.
The first compile will take a while since the libraries in libs/ need to be
downloaded and compiled.

To compile for a maximum kmer size of 31:

    make

to compile for a maximum kmer size of 63:

    make MAXK=63

Executables appear in the `bin/` directory. To update the libraries included:

    cd libs; make clean; make

Commands
--------

    usage: ctx31 <command> [options] <args>
    version: 5854d1c; zlib: 1.2.5

    Command:  build       FASTA/FASTQ/BAM -> cortex graph file
              view        view and check a cortex graph file (.ctx)
              healthcheck load and check a cortex graph file (.ctx)
              clean       clean errors from a graph
              join        combine graphs, filter graph intersections
              supernodes  pull out supernodes
              subgraph    filter a subgraph using seed kmers
              reads       filter reads against a graph
              extend      extend contigs using a graph
              contigs     pull out contigs for a sample
              inferedges  infer graph edges before calling `thread`
              thread      thread reads through cleaned population
              pview       view read threading information
              pmerge      merge path files (.ctp)
              call        call variants
              unique      remove duplicated bubbles, produce VCF
              place       place variants and genotype

      Type a command with no arguments to see help.

    Common Options:
      -m --memory <M>      Memory e.g. 1GB [default: 1GB]
      -n --nkmers <H>      Hash entries [default: 4M, ~4 million]
      -c --ncols <C>       Number of graph colours to load at once [default: 1]
      -t --threads <T>     Number of threads [default: 2]
      -k --kmer <K>        Kmer size [default: read from graph files]
      -f --file <file>     Input file
      -o --out <file>      Output file
      -p --paths <in.ctp>  Assembly file

Getting Helps
-------------

Type a command with no arguments to see usage.

Check out the [wiki](https://github.com/noporpoise/ninja-cortex/wiki)

If you find a bug please submitted an [Issue](https://github.com/noporpoise/ninja-cortex/issues) on github.

Cortex mailing list: https://groups.google.com/forum/#!forum/cortex_var

Isaac Turner: turner.isaac@gmail.com

Code And Contributing
------------

Issues can be submitted on github. Pull requests welcome. Please add your name
to the AUTHORS file.

Code should compile on mac/linux with clang/gcc without errors or warnings.

Code is organised as:
* libs/       included library code from other projects
* src/basic   files that do not depend on MAX_KMER_SIZE
* src/kmer    files that need recompiling based on different MAX_KMER_SIZE
* src/tools   one file per cortex command
* src/main    files with a main function go in here

Files only link to files that are above them in the list above. E.g. src/kmer/*
files only include files in src/kmer/, src/basic/ and libs/.

License: GPLv2
--------------

Bundled libraries may have different licenses:
* [GNU Science Library](http://www.gnu.org/software/gsl/) (GPL)
* [CityHash](https://code.google.com/p/cityhash/) (MIT)
* [lookup3](http://burtleburtle.net/bob/c/lookup3.c) (Public Domain)
* [htslib](https://github.com/samtools/htslib) (MIT)
* [bcftools](https://github.com/samtools/bcftools) (MIT)
* [vcflib](https://github.com/ekg/vcflib) (MIT)
* [seq_file](https://github.com/noporpoise/seq_file) (Public Domain)
* [string_buffer](https://github.com/noporpoise/string_buffer) (Public Domain)
* [BitArray](https://github.com/noporpoise/BitArray) (Public Domain)
* [msg-pool](https://github.com/noporpoise/msg-pool) (Public Domain)
* [seq-align](https://github.com/noporpoise/seq-align) (Public Domain)

Used in testing:
* [bioinf-perl](https://github.com/noporpoise/bioinf-perl) (Public Domain)

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
