SHELL=/bin/bash

# make          <- compile mccortex dependencies
# make all      <- compile all libraries (deps + analysis)
# make clean    <- clean all libraries

CORE=xxHash htslib string_buffer bit_array seq_file seq-align msg-pool sort_r madcrowlib misc
COMMON=bcftools samtools bwa readsim bioinf-perl maximal_substrs
ANALYSIS=vcflib

ALLTGTS=$(CORE) $(COMMON) $(ANALYSIS)

# Only compile dependencies for compiling ninja-cortex itself
core: $(CORE)

common: $(COMMON)

analysis: $(ANALYSIS)

all: core common analysis

#
# Download
#
gsl-1.16/configure:
	wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
	tar xfz gsl-1.16.tar.gz
	cd gsl-1.16; ./configure; make

# Download of git hosted libraries now done with git submodules

#
# Build
#
gsl: gsl-1.16/configure
	cd gsl-1.16 && make clean && ./configure && make

xxHash: xxHash/Makefile
	cd xxHash && make

htslib: htslib/Makefile
	cd htslib && make

bcftools: htslib bcftools/Makefile
	cd bcftools && make

samtools: htslib samtools/Makefile
	cd samtools && make

bwa: bwa/Makefile
	cd bwa && make

vcflib: vcflib/Makefile
	cd vcflib #&& make

string_buffer: string_buffer/Makefile
	cd string_buffer && make all

bit_array: bit_array/Makefile
	cd bit_array && make all

seq_file: seq_file/Makefile htslib
	cd seq_file && make HTSLIB=../htslib all

seq-align: seq-align/Makefile seq_file bit_array string_buffer sort_r
	cd seq-align && make LIBS_PATH=.. all

readsim: readsim/Makefile htslib
	cd readsim && make HTSLIB=../htslib

msg-pool: msg-pool/Makefile
	cd msg-pool && make all

bioinf-perl: bioinf-perl/README
	cd bioinf-perl

sort_r: sort_r/Makefile
	cd sort_r && make all

madcrowlib: madcrowlib/Makefile
	cd madcrowlib && make all

maximal_substrs:
	cd maximal_substrs && make

misc: misc/Makefile
	cd misc && make

#
# Clean
#
clean:
	cd xxHash && make clean
	cd htslib && make clean
	cd bcftools && make clean
	cd samtools && make clean
	cd vcflib && make clean
	cd string_buffer && make clean
	cd bit_array && make clean
	cd seq_file && make clean
	&& cd seq-align && make clean
	cd readsim && make clean
	&& cd msg-pool && make clean
	cd sort_r && make clean
	cd madcrowlib && make clean
	cd maximal_substrs && make clean
	cd misc && make clean

.PHONY: all clean update core
.PHONY: $(ALLTGTS)
