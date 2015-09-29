SHELL=/bin/bash

# make          <- compile mccortex dependencies
# make all      <- compile all libraries (deps + analysis)
# make clean    <- clean all libraries

CORE=xxHash htslib string_buffer bit_array seq_file seq-align msg-pool sort_r madcrowlib misc
COMMON=bcftools samtools bwa readsim bioinf-perl maximal_substrs vcfhp

ALLTGTS=$(CORE) $(COMMON)

all: core common

# Only compile dependencies for compiling ninja-cortex itself
core: $(CORE)

common: $(COMMON)

# Download of git hosted libraries now done with git submodules

#
# Build
#
xxHash: xxHash/Makefile
	cd xxHash && $(MAKE)

htslib: htslib/Makefile
	cd htslib && $(MAKE)

bcftools: htslib bcftools/Makefile
	cd bcftools && $(MAKE)

samtools: htslib samtools/Makefile
	cd samtools && $(MAKE)

bwa: bwa/Makefile
	cd bwa && $(MAKE)

string_buffer: string_buffer/Makefile
	cd string_buffer && $(MAKE) all

bit_array: bit_array/Makefile
	cd bit_array && $(MAKE) all

seq_file: seq_file/Makefile htslib
	cd seq_file && $(MAKE) HTSLIB=../htslib all

seq-align: seq-align/Makefile seq_file bit_array string_buffer sort_r
	cd seq-align && $(MAKE) LIBS_PATH=.. all

readsim: readsim/Makefile htslib
	cd readsim && $(MAKE) HTSLIB=../htslib

msg-pool: msg-pool/Makefile
	cd msg-pool && $(MAKE) all

bioinf-perl: bioinf-perl/README
	cd bioinf-perl

sort_r: sort_r/Makefile
	cd sort_r && $(MAKE) all

madcrowlib: madcrowlib/Makefile
	cd madcrowlib && $(MAKE) all

maximal_substrs:
	cd maximal_substrs && $(MAKE)

misc: misc/Makefile
	cd misc && $(MAKE)

vcfhp: vcfhp/Makefile htslib
	cd vcfhp && $(MAKE) all test

#
# Clean
#
clean:
	cd xxHash && $(MAKE) clean
	cd htslib && $(MAKE) clean
	cd bcftools && $(MAKE) clean
	cd samtools && $(MAKE) clean
	cd string_buffer && $(MAKE) clean
	cd bit_array && $(MAKE) clean
	cd seq_file && $(MAKE) clean
	cd seq-align && $(MAKE) clean
	cd readsim && $(MAKE) clean
	cd msg-pool && $(MAKE) clean
	cd sort_r && $(MAKE) clean
	cd madcrowlib && $(MAKE) clean
	cd maximal_substrs && $(MAKE) clean
	cd misc && $(MAKE) clean

.PHONY: all clean update core
.PHONY: $(ALLTGTS)
