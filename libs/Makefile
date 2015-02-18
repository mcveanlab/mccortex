# make          <- compile all libraries
# make all      <- compile all libraries
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

xxHash/Makefile:
	git clone https://github.com/Cyan4973/xxHash.git

htslib/Makefile:
	git clone --branch master git://github.com/samtools/htslib.git htslib

bcftools/Makefile:
	git clone --branch master https://github.com/samtools/bcftools.git

samtools/Makefile:
	git clone --branch master https://github.com/samtools/samtools.git

bwa/Makefile:
	git clone https://github.com/lh3/bwa.git

vcflib/Makefile:
	git clone --recursive https://github.com/ekg/vcflib.git

string_buffer/Makefile:
	git clone git://github.com/noporpoise/string_buffer.git

bit_array/Makefile:
	git clone git://github.com/noporpoise/BitArray.git bit_array

seq_file/Makefile:
	git clone git://github.com/noporpoise/seq_file.git

seq-align/Makefile:
	git clone git://github.com/noporpoise/seq-align.git

readsim/Makefile:
	git clone git://github.com/noporpoise/readsim.git

msg-pool/Makefile:
	git clone https://github.com/noporpoise/msg-pool.git

bioinf-perl/README:
	git clone https://github.com/noporpoise/bioinf-perl.git

sort_r/Makefile:
	git clone https://github.com/noporpoise/sort_r.git

madcrowlib/Makefile:
	git clone https://github.com/noporpoise/madcrowlib.git

#
# Build
#
gsl: gsl-1.16/configure
	cd gsl-1.16 && make clean && ./configure && make

xxHash: xxHash/Makefile
	cd xxHash && git pull && make

htslib: htslib/Makefile
	cd htslib && git pull && make

bcftools: htslib bcftools/Makefile
	cd bcftools && git pull && make

samtools: htslib samtools/Makefile
	cd samtools && git pull && make

bwa: bwa/Makefile
	cd bwa && git pull && make

vcflib: vcflib/Makefile
	cd vcflib && git pull #&& make

string_buffer: string_buffer/Makefile
	cd string_buffer && git pull && make

bit_array: bit_array/Makefile
	cd bit_array && git pull && make

seq_file: seq_file/Makefile htslib
	cd seq_file && git pull && make HTSLIB=../htslib

seq-align: seq-align/Makefile seq_file bit_array string_buffer sort_r
	cd seq-align && git pull && make LIBS_PATH=../

readsim: readsim/Makefile htslib
	cd readsim && git pull && make HTSLIB=../htslib

msg-pool: msg-pool/Makefile
	cd msg-pool && git pull && make

bioinf-perl: bioinf-perl/README
	cd bioinf-perl && git pull

sort_r: sort_r/Makefile
	cd sort_r && git pull && make

madcrowlib: madcrowlib/Makefile
	cd madcrowlib && git pull && make

maximal_substrs:
	cd maximal_substrs && make

misc: misc/Makefile
	cd misc && make

#
# Clean
#
clean:
	[ -d gsl-1.16 ] && cd gsl-1.16 && make clean
	[ -d xxHash ] && cd xxHash && make clean
	[ -d htslib ] && cd htslib && make clean
	[ -d bcftools ] && cd bcftools && make clean
	[ -d samtools ] && cd samtools && make clean
	[ -d vcflib ] && cd vcflib && make clean
	[ -d string_buffer ] && cd string_buffer && make clean
	[ -d bit_array ] && cd bit_array && make clean
	[ -d seq_file ] && cd seq_file && make clean
	[ -d seq-align ] && cd seq-align && make clean
	[ -d readsim ] && cd readsim && make clean
	[ -d msg-pool ] && cd msg-pool && make clean
	[ -d sort_r ] && cd sort_r && make clean
	[ -d madcrowlib ] && cd madcrowlib && make clean
	[ -d maximal_substrs ] && cd maximal_substrs && make clean
	[ -d misc ] && cd misc && make clean

.PHONY: all clean update core
.PHONY: $(ALLTGTS)
