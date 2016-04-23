McCortex: Population De Novo Assembly and Variant Calling
===============================================

Multi-sample de novo assembly and variant calling using de bruijn graphs.
Variant calling with and without a reference genome. Between closely related
samples or highly diverged ones. From bacterial to mammalian genomes. Minimal
configuration. And it's free.

Isaac Turner's experimental rewrite of cortex_var, to handle larger populations
with better genome assembly. PhD supervisor: Prof Gil McVean. Collaborators: Zam Iqbal, Kiran Garimella. Based at the Wellcome Trust Centre for Human Genetics, University of Oxford.

*Note: Currently under development.* Expect bugs, fixes and vague documentation until we hit our first release. Feel free to try out McCortex and watch this space for the release.

9 Nov 2015

Branch         | Status
---------------|--------
master:        | [![Build Status](https://travis-ci.org/mcveanlab/mccortex.svg?branch=master)](https://travis-ci.org/mcveanlab/mccortex)
develop:       | [![Build Status](https://travis-ci.org/mcveanlab/mccortex.svg?branch=develop)](https://travis-ci.org/mcveanlab/mccortex)
code analysis: | [![Coverity Scan Build Status](https://scan.coverity.com/projects/2329/badge.svg)](https://scan.coverity.com/projects/2329)

Build
-----

McCortex compiles with clang and gcc. Tested on Mac OS X and linux. Requires zlib.
Download with:

    git clone --recursive https://github.com/mcveanlab/mccortex

To compile for a maximum kmer size of 31:

    make all

to compile for a maximum kmer size of 63:

    make MAXK=63 all

Executables appear in the `bin/` directory.


Quickstart: Variant calling
---------------------------

Download and compile McCortex. Can be in any directory, later I'll assume it's in `~/mccortex/`:

    git clone --recursive https://github.com/mcveanlab/mccortex
    cd mccortex
    make all MAXK=31
    make all MAXK=63

Now write a file detailing your samples and their data. Columns are separated by one or more spaces/tabs. File entries are separated by commas. Paired-end read files are separated by a colon ':'. File paths can be relative to the current directory or absolute. Most fileformats are supported:

    cd /path/to/your/data
    echo "#sample_name  SE_files   PE_files                     interleaved_files" >  samples.txt
    echo "Mickey        a.fa,b.fa  reads.1.fq.gz:reads.2.fq.gz  ."                 >> samples.txt
    echo "Minney        .          reads.1.fq.gz:reads.2.fq.gz  in.bam"            >> samples.txt
    echo "Pluto         seq.fq     .                            pluto.cram"        >> samples.txt

Create a job file from your sample file (`samples.txt`). All output will go into the directory we specify (`mc_calls`). We also specify the kmer(s) to use. We'll run at `k=31` and `k=61` and merge the results.

If your data are haploid, we set `--ploidy 1`:

    ~/mccortex/scripts/make-pipeline.pl -r /path/to/ref.fa --ploidy 1 31,61 mc_calls samples.txt > job.k31.k61.mk

If your samples are human, you have a mix of haploid and diploid chromosomes. Therefore you need to specify which samples have only one copy of `chrX` and one of `chrY`. The format is `-P <sample>:<chr>:<ploidy>` where `<sample>` and `<chr>` can be comma-separated lists. Ploidy arguments are read in order.

    ~/mccortex/scripts/make-pipeline.pl -r /path/to/ref.fa --ploidy "-P .:.:2 -P .:chrY:1 -P Mickey:chrX:1" 31,61 mc_calls samples.txt > job.k31.k61.mk

Now you're ready to run. You'll need to pass:
- path to McCortex `CTXDIR=`
- how much memory to use `MEM=`  (2GB for ten E. coli, 100GB for a human)
- number of threads to use `NTHREADS=`

Run the job file:

    make -f job.k31.k61.mk CTXDIR=~/mccortex MEM=100GB NTHREADS=8 \
                           JOINT_CALLING=yes USE_LINKS=no brk-geno-vcf

For a human, running time will be about 8 hours for a single sample and use about 100GB of RAM. 

Job finished? Your results are in: `mc_calls/vcfs/breakpoints.joint.plain.k31.k61.geno.vcf.gz`.

Something go wrong? Take a look at the log file of the last command that ran. You may need to increase memory or compile for a different `MAXK=` value. Once you've fixed the issue, just rerun the `make -f job...` command. Add `--dry-run` to the `make` command to see which commands are going to be run without running them. 


Commands
--------

    usage: mccortex31 <command> [options] <args>
    version: ctx=XXXX zlib=1.2.5 htslib=1.2.1 ASSERTS=ON hash=Lookup3 CHECKS=ON k=3..31
    
    Commands:   breakpoints  use a trusted assembled genome to call large events
                bubbles      find bubbles in graph which are potential variants
                build        construct cortex graph from FASTA/FASTQ/BAM
                calls2vcf    convert bubble/breakpoint calls to VCF
                check        load and check graph (.ctx) and path (.ctp) files
                clean        clean errors from a graph
                contigs      assemble contigs for a sample
                correct      error correct reads
                coverage     print contig coverage
                dist         make colour kmer distance matrix
                index        index a sorted cortex graph file
                inferedges   infer graph edges between kmers before calling `thread`
                join         combine graphs, filter graph intersections
                links        clean and plot link files (.ctp)
                pjoin        merge path files (.ctp)
                popbubbles   pop bubbles in the population graph
                pview        text view of a cortex path file (.ctp)
                reads        filter reads against a graph
                rmsubstr     reduce set of strings to remove substrings
                server       interactively query the graph
                sort         sort the kmers in a graph file
                subgraph     filter a subgraph using seed kmers
                thread       thread reads through cleaned graph to make links
                uniqkmers    generate random unique kmers
                unitigs      pull out unitigs in FASTA, DOT or GFA format
                vcfcov       coverage of a VCF against cortex graphs
                vcfgeno      genotype a VCF after running vcfcov
                view         text view of a cortex graph file (.ctx)
    
      Type a command with no arguments to see help.
    
    Common Options:
      -h, --help            Help message
      -q, --quiet           Silence status output normally printed to STDERR
      -f, --force           Overwrite output files if they already exist
      -m, --memory <M>      Memory e.g. 1GB [default: 1GB]
      -n, --nkmers <H>      Hash entries [default: 4M, ~4 million]
      -t, --threads <T>     Limit on proccessing threads [default: 2]
      -o, --out <file>      Output file
      -p, --paths <in.ctp>  Assembly file to load (can specify multiple times)

Getting Helps
-------------

Type a command with no arguments to see usage. The following may also be useful:
* [wiki](https://github.com/mcveanlab/mccortex/wiki)
* [website](http://mcveanlab.github.io/mccortex)
* [mailing list](https://groups.google.com/forum/#!forum/cortex_var)
* Report a [bug / feature request](https://github.com/mcveanlab/mccortex/issues) on GitHub
* Email me: Isaac Turner <turner.isaac@gmail.com>

Live chat (email me to fix a time):
* [HipChat](http://www.hipchat.com/gbF6Zf4k3) to instant message -- please email me first to arrange a time
* [![Gitter https://gitter.im/mcveanlab/mccortex](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/mcveanlab/mccortex?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

Code And Contributing
---------------------

Issues can be submitted on github. Pull requests welcome. Please add your name
to the AUTHORS file. Code should compile on mac/linux with clang/gcc without errors or warnings.

More on the [wiki](https://github.com/mcveanlab/mccortex/wiki/Contributing)

Unit tests are run with `make test` and integration tests with `cd tests; ./run`. Both of these test suites are run automatically with Travis CI when commits are pushed to GitHub. 

Static analysis can be run with [cppcheck](http://cppcheck.sourceforge.net):

    cppcheck src

or with [clang](http://clang-analyzer.llvm.org):

    rm -rf bin/mccortex31
    scan-build make RECOMPILE=1

Occasionally we also run Coverity Scan. This is done by pushing to the `coverity_scan` branch on github, which triggers Travis CI to upload the latest code to Coverity.

[![Coverity Scan Build Status](https://scan.coverity.com/projects/2329/badge.svg)](https://scan.coverity.com/projects/2329)

    git checkout coverity_scan
    git merge develop
    git checkout --ours .travis.yml
    git checkout --ours configure

License: MIT
------------

Bundled libraries may have different licenses:
* [BitArray](https://github.com/noporpoise/BitArray) (Public Domain)
* [cJSON](http://http://sourceforge.net/projects/cjson/) (MIT)
* [CityHash](https://code.google.com/p/cityhash/) (MIT)
* [htslib](https://github.com/samtools/htslib) (MIT)
* [lookup3](http://burtleburtle.net/bob/c/lookup3.c) (Public Domain)
* [madcrowlib](https://github.com/noporpoise/madcrowlib) (MIT)
* [msg-pool](https://github.com/noporpoise/msg-pool) (Public Domain)
* [seq-align](https://github.com/noporpoise/seq-align) (Public Domain)
* [seq_file](https://github.com/noporpoise/seq_file) (Public Domain)
* [sort_r](https://github.com/noporpoise/sort_r) (Public Domain)
* [carrays](https://github.com/noporpoise/carrays) (Public Domain)
* [string_buffer](https://github.com/noporpoise/string_buffer) (Public Domain)
* [xxHash](https://github.com/Cyan4973/xxHash.git) (BSD)

Used in testing:
* [bcftools](https://github.com/samtools/bcftools) (MIT)
* [bioinf-perl](https://github.com/noporpoise/bioinf-perl) (Public Domain)
* [bwa](https://github.com/lh3/bwa) (MIT)
* [readsim](https://github.com/noporpoise/readsim) (Public Domain)
* [samtools](https://github.com/samtools/samtools) (MIT)

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
