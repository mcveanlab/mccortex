# make          <- compile all libraries
# make all      <- compile all libraries
# make clean    <- clean all libraries

UPDATE=htslib bcftools vcflib string_buffer bit_array seq_file seq-align \
       readsim msg-pool sort_r misc bioinf-perl

ALLTGTS=$(UPDATE) gsl

# Only compile dependencies for compiling ninja-cortex itself
core: htslib string_buffer bit_array seq_file seq-align msg-pool sort_r misc

update: $(UPDATE)
all: $(ALLTGTS)

#
# Download
#
gsl-1.16/configure:
	wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
	tar xfz gsl-1.16.tar.gz
	cd gsl-1.16; ./configure; make

htslib/Makefile:
	git clone --branch 'bcftools+calling' git://github.com/samtools/htslib.git htslib

bcftools/Makefile:
	git clone https://github.com/samtools/bcftools.git

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

#
# Build
#
gsl: gsl-1.16/configure
	cd gsl-1.16 && make clean && ./configure && make

htslib: htslib/Makefile
	cd htslib && git pull && make

bcftools: htslib bcftools/Makefile
	cd bcftools && git pull && make

vcflib: vcflib/Makefile
	cd vcflib && git pull && make

string_buffer: string_buffer/Makefile
	cd string_buffer && git pull && make

bit_array: bit_array/Makefile
	cd bit_array && git pull && make

seq_file: seq_file/Makefile htslib
	cd seq_file && git pull && make HTSLIB=../htslib

seq-align: seq-align/Makefile bit_array string_buffer
	cd seq-align && git pull && make LIBS_PATH=../

readsim: readsim/Makefile htslib
	cd readsim && git pull && make HTSLIB=../htslib

msg-pool: msg-pool/Makefile
	cd msg-pool && git pull && make

bioinf-perl: bioinf-perl/README
	cd bioinf-perl && git pull

sort_r: sort_r/Makefile
	cd sort_r && git pull && make

misc: misc/Makefile
	cd misc && make

#
# Clean
#
clean:
	# [ -d gsl-1.16 ] && cd gsl-1.16 && make clean
	[ -d htslib ] && cd htslib && make clean
	[ -d bcftools ] && cd bcftools && make clean
	[ -d vcflib ] && cd vcflib && make clean
	[ -d string_buffer ] && cd string_buffer && make clean
	[ -d bit_array ] && cd bit_array && make clean
	[ -d seq_file ] && cd seq_file && make clean
	[ -d seq-align ] && cd seq-align && make clean
	[ -d readsim ] && cd readsim && make clean
	[ -d msg-pool ] && cd msg-pool && make clean
	[ -d sort_r ] && cd sort_r && make clean
	[ -d misc ] && cd misc && make clean

.PHONY: all clean update core
.PHONY: $(ALLTGTS)
