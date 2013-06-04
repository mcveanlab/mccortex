# make          <- compile all libraries
# make all      <- compile all libraries
# make pull     <- git pull all repos
# make common   <- compile all repos (except gsl)
# make update   <- pull + common

all: htslib string_buffer bit_array seq_file seq-align gsl

gsl-1.15/configure:
	wget ftp://ftp.gnu.org/gnu/gsl/gsl-1.15.tar.gz
	tar xfvz gsl-1.15.tar.gz

htslib/Makefile:
	git clone git://github.com/samtools/htslib.git

bit_array/Makefile:
	git clone git://github.com/noporpoise/BitArray.git
	mv BitArray bit_array

string_buffer/Makefile:
	git clone git://github.com/noporpoise/string_buffer.git

seq_file/Makefile:
	git clone git://github.com/noporpoise/seq_file.git

seq-align/Makefile:
	git clone git://github.com/noporpoise/seq-align.git

cortex_binary_reader/Makefile:
	git clone git://github.com/noporpoise/cortex_binary_reader.git

# seqan-svn:
# 	git clone git://github.com/seqan/seqan-svn.git

fetch: gsl-1.15/configure htslib/Makefile bit_array/Makefile \
	     string_buffer/Makefile seq_file/Makefile \
	     seq-align/Makefile cortex_binary_reader/Makefile

pull:
	cd htslib; git pull
	cd string_buffer; git pull
	cd bit_array; git pull
	cd seq_file; git pull
	cd seq-align; git pull
	cd cortex_binary_reader; git pull
	cd bioinf-perl; git pull
	# cd seqan-svn; git pull

update: pull common

common: htslib string_buffer bit_array seq_file seq-align cortex_binary_reader

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

gsl: gsl-1.15/Makefile
	cd gsl-1.15; make clean; ./configure; make

cortex_binary_reader: cortex_binary_reader/Makefile
	cd cortex_binary_reader; make;

# seqan: seqan-svn
#.PHONY: seqan

clean:
	cd htslib; make clean
	cd string_buffer; make clean
	cd bit_array; make clean
	cd seq_file; make clean
	cd seq-align; make clean
	cd gsl-1.15; make clean
	cd cortex_binary_reader; make clean

.PHONY: all clean fetch pull common update
.PHONY: htslib string_buffer bit_array seq_file seq-align cortex_binary_reader
