SHELL=/bin/bash

# make          <- compile mccortex dependencies
# make all      <- compile all libraries (deps + analysis)
# make clean    <- clean all libraries

CORE=xxHash htslib string_buffer bit_array seq_file seq-align msg-pool sort_r madcrowlib carrays misc
OTHER=bcftools samtools bwa bfc readsim bioinf-perl maximal_substrs vcf-slim

ALLTGTS=$(CORE) $(OTHER)

# Only compile dependencies for compiling McCortex itself
core: $(CORE)
other: $(OTHER)

all: core other

# Download of git hosted libraries now done with git submodules

#
# Build
#
xxHash: xxHash/Makefile
	cd xxHash && $(MAKE)

htslib: htslib/Makefile
	# Run configure if config makefile not created
	cd htslib && ( [ -f config.mk ] || \
    ( autoheader && autoconf && ./configure --disable-lzma --disable-bz2 --disable-libcurl )) && \
	  $(MAKE)

bcftools: htslib bcftools/Makefile
	cd bcftools && $(MAKE)

samtools: htslib samtools/Makefile
	cd samtools && $(MAKE)

bwa: bwa/Makefile
	cd bwa && $(MAKE)

bfc: bfc/Makefile
	cd bfc && $(MAKE)

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

carrays: carrays/Makefile
	cd carrays && $(MAKE) all

maximal_substrs: maximal_substrs/Makefile
	cd maximal_substrs && $(MAKE)

misc: misc/Makefile
	cd misc && $(MAKE)

vcf-slim: vcf-slim/Makefile htslib
	cd vcf-slim && $(MAKE) all

#
# Run tests
#
test:
	cd htslib && $(MAKE) test
	cd bit_array && $(MAKE) test
	cd string_buffer && $(MAKE) test
	cd seq-align && $(MAKE) test
	cd msg-pool && $(MAKE) test
	cd sort_r && $(MAKE) test
	cd carrays && $(MAKE) test
	cd vcf-slim && $(MAKE) test
	cd madcrowlib && $(MAKE) test

#
# Clean
#
clean:
	cd xxHash && $(MAKE) clean
	cd htslib && $(MAKE) clean
	cd bcftools && $(MAKE) clean
	cd samtools && $(MAKE) clean
	cd bfc && $(MAKE) clean
	cd string_buffer && $(MAKE) clean
	cd bit_array && $(MAKE) clean
	cd seq_file && $(MAKE) clean
	cd seq-align && $(MAKE) clean
	cd readsim && $(MAKE) clean
	cd msg-pool && $(MAKE) clean
	cd sort_r && $(MAKE) clean
	cd madcrowlib && $(MAKE) clean
	cd carrays && $(MAKE) clean
	cd maximal_substrs && $(MAKE) clean
	cd misc && $(MAKE) clean

.PHONY: all clean test update core
.PHONY: $(ALLTGTS)
