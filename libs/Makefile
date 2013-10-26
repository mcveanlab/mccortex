# make          <- compile all libraries
# make all      <- compile all libraries
# make pull     <- git pull all repos
# make common   <- compile all repos (except gsl)
# make update   <- pull + common

all: gsl htslib string_buffer bit_array seq_file seq-align readsim

gsl-1.16/configure:
	wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.16.tar.gz
	tar xfz gsl-1.16.tar.gz
	cd gsl-1.16; ./configure; make

htslib/Makefile:
	git clone git://github.com/samtools/htslib.git

string_buffer/Makefile:
	git clone git://github.com/noporpoise/string_buffer.git

bit_array/Makefile:
	git clone git://github.com/noporpoise/BitArray.git
	mv BitArray bit_array

seq_file/Makefile:
	git clone git://github.com/noporpoise/seq_file.git

seq-align/Makefile:
	git clone git://github.com/noporpoise/seq-align.git

readsim/Makefile:
	git clone git://github.com/noporpoise/readsim.git

fetch: gsl-1.16/configure htslib/Makefile \
	     string_buffer/Makefile bit_array/Makefile \
	     seq_file/Makefile seq-align/Makefile readsim/Makefile

pull:
	cd htslib; git pull
	cd string_buffer; git pull
	cd bit_array; git pull
	cd seq_file; git pull
	cd seq-align; git pull
	cd readsim; git pull
	cd bioinf-perl; git pull

update: pull common

common: htslib string_buffer bit_array seq_file seq-align readsim

gsl: gsl-1.16/configure
	cd gsl-1.16; make clean; ./configure; make

htslib: htslib/Makefile
	cd htslib; make

string_buffer: string_buffer/Makefile
	cd string_buffer; make

bit_array: bit_array/Makefile
	cd bit_array; make

seq_file: seq_file/Makefile
	cd seq_file; make HTSLIB=../htslib

seq-align: seq-align/Makefile
	cd seq-align; make LIBS_PATH=../

readsim: readsim/Makefile
	cd readsim; make HTSLIB=../htslib

clean:
	cd gsl-1.16; make clean
	cd htslib; make clean
	cd string_buffer; make clean
	cd bit_array; make clean
	cd seq_file; make clean
	cd seq-align; make clean
	cd readsim; make clean

.PHONY: all clean fetch pull common update
.PHONY: gsl htslib string_buffer bit_array seq_file seq-align readsim
