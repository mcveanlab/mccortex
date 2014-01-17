# make          <- compile all libraries
# make all      <- compile all libraries
# make clean    <- clean all libraries

UPDATE=htslib bcftools vcflib string_buffer bit_array seq_file seq-align readsim msg-pool misc
ALLTGTS=$(UPDATE) gsl

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
	# nothing to do for bioinf-perl

misc: misc/Makefile
	cd misc && make

#
# Clean
#
clean:
	# cd gsl-1.16 && make clean
	cd htslib && make clean
	cd bcftools && make clean
	cd vcflib && make clean
	cd string_buffer && make clean
	cd bit_array && make clean
	cd seq_file && make clean
	cd seq-align && make clean
	cd readsim && make clean
	cd msg-pool && make clean
	cd misc && make clean

.PHONY: all clean update
.PHONY: $(ALLTGTS)
