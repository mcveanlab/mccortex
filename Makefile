# ifndef CC
  CC = gcc
# endif

#Arguments (all are optional):
# MAXK=31
# NUM_COLS=12
# LIB_PATH=/usr/local/lib/ or where ever your system libraries are
# DEBUG=1
# CITY_HASH=1   (experimental -- use CityHash hash function)

# e.g. for WTCHG cluster3
# a) pass LIB_PATH=/usr/local/lib/ to compile
# b) also set LD_LIBRARY_PATH=/usr/local/lib/:$(LD_LIBRARY_PATH) before running

# Targets:
#  make clean
#  make                       (build cortex_var)
#  make cortex_var
#  make run_basic_tests
#  make run_cortex_var_tests
#  make run_hash_table_tests
#  make all                   (build all targets)

# The following does not compile and is not part of `make all`
# make run_cmdline_tests

BIN = bin
TEMP_TEST_DIR = data/tempfiles_can_be_deleted

# Use bash as shell
SHELL := /bin/bash

BITFIELDS=$(shell echo $$((($(MAXK)+31)/32)))

ifndef MAXK
  BITFIELDS = 1
  MAXK = 31
else
  BITFIELDS=$(shell echo $$((($(MAXK)+31)/32)))
endif

ifeq ($(BITFIELDS),0)
  $(error Invalid MAXK value '$(MAXK)'. Please choose from 31,63,95,..(32*n-1) [default: 31])
endif

MAX_KMER_SIZE=$(shell echo $$(($(BITFIELDS)*32-1)))
MIN_KMER_SIZE=$(shell echo $$(($(MAX_KMER_SIZE)-30)))

ifeq ($(MIN_KMER_SIZE),1)
	MIN_KMER_SIZE=3
endif

ifndef NUM_COLS
  NUM_COLS = 1
endif

# Use City hash instead of lookup3?
ifdef CITY_HASH
	HASH_KEY_FLAGS = -DUSE_CITY_HASH=1
endif

# Test if running on a mac
UNAME=$(shell uname)
ifeq ($(UNAME),Darwin)
  MAC = 1
endif

# Library paths
IDIR_GSL_HEADERS = libs/gsl-1.15
IDIR_HTS = libs/htslib/htslib
IDIR_STRS = libs/string_buffer
IDIR_SEQ = libs/seq_file
IDIR_ALIGN = libs/seq-align/src
IDIR_SEQAN = libs/seqan
IDIR_CITY = libs/city_hash
IDIR_LOOKUP3 = libs/lookup3_hash

LIB_GSL = libs/gsl-1.15/.libs/libgsl.a
LIB_HTS = libs/htslib/htslib/libhts.a
LIB_ALIGN = libs/seq-align/src/libalign.a

# Various places to look for Cunit (first found is the one used)
IDIR_CUNIT1=/opt/local/include/CUnit*/
IDIR_CUNIT2=~/lib/CUnit*/doc/headers/

LDIR_CUNIT1=/opt/local/lib
LDIR_CUNIT2=~/lib/CUnit*/CUnit/Sources/.libs/

IDIR_CUNIT=$(shell if ls $(IDIR_CUNIT1) &> /dev/null; then echo $(IDIR_CUNIT1); \
                 elif ls $(IDIR_CUNIT2) &> /dev/null; then echo $(IDIR_CUNIT2); fi)

LDIR_CUNIT=$(shell if ls $(LDIR_CUNIT1) &> /dev/null; then echo $(LDIR_CUNIT1); \
                 elif ls $(LDIR_CUNIT2) &> /dev/null; then echo $(LDIR_CUNIT2); fi)

INCLUDES = -I src/common/ \
           -I $(IDIR_CITY) -I $(IDIR_LOOKUP3) \
           -I $(IDIR_GSL_HEADERS) -I $(IDIR_HTS) \
           -I $(IDIR_ALIGN) -I $(IDIR_SEQ) -I $(IDIR_STRS)

ifdef MAC
	MACFLAG = -fnested-functions
endif

# Comment out this line to turn off adding the commit version
# (it already checks if hg is installed)
#VERSION_STR=$(shell if [ `command -v hg` ]; then echo ' (commit' `hg id --num --id`')'; else echo; fi)

WARNS = -Wall -Wextra -Winit-self -Wmissing-include-dirs \
        -Wstrict-aliasing -Wdiv-by-zero \
        -Wcast-qual -Wcast-align -Wmissing-noreturn \
        -Wwrite-strings -Waggregate-return \
        -Wwrite-strings -Wundef -std=c99

# -Wshadow -Wconversion -Wshorten-64-to-32 -Woverlength-strings -Wunreachable-code
# -Wenum-compare -Wlogical-op -Wfloat-equal -Wbad-function-cast

CFLAGS = $(ARCH) $(WARNS) $(MACFLAG) -DVERSION_STR='"$(VERSION_STR)"' \
         -DMAX_KMER_SIZE=$(MAX_KMER_SIZE) -DMIN_KMER_SIZE=$(MIN_KMER_SIZE) \
         -DNUM_BITFIELDS_IN_BKMER=$(BITFIELDS) -DNUM_OF_COLOURS=$(NUM_COLS) \
         $(HASH_KEY_FLAGS)

# Optimisations tags for testing
OPT_TESTS = -O0 -Wstack-protector -fstack-protector

# If not debugging, add optimisations and -DNDEBUG=1 to turn off assert() calls
# only applies to cortex_var
ifdef DEBUG
	OPT = $(OPT_TESTS)
	DEBUG_ARGS = -g -ggdb -DDEBUG=1
	DEBUG_LIBS = 
else
	OPT = -O3 -DNDEBUG=1
	DEBUG_ARGS = 
	DEBUG_LIBS =
endif

# Resolve some issues linking libz:
# 1) pass LIB_PATH=/usr/local/lib/ to compile on WTCHG cluster3
# 2) also set LD_LIBRARY_PATH=/usr/local/lib/:$(LD_LIBRARY_PATH)
ifdef LIB_PATH
	INCLUDES := -I $(LIB_PATH) -L $(LIB_PATH) $(INCLUDES)
endif

INCLUDES_TESTS = $(INCLUDES) -I $(IDIR_BASIC_TESTS) \
                -I $(IDIR_HASH_TABLE_TESTS) -I $(IDIR_CORTEX_VAR_TESTS) \
                -I$(IDIR_CUNIT) -L$(LDIR_CUNIT)

# Library linking
LIB_OBJS = $(LIB_GSL) $(LIB_HTS) $(LIB_ALIGN)
LINK_LIBS = -lpthread -lz -lm $(DEBUG_LIBS)

COMMON_SRCS = $(wildcard src/common/*.c) \
              libs/city_hash/city.c libs/lookup3_hash/lookup3.c \
              libs/string_buffer/string_buffer.c

COMMON_HDRS = $(wildcard src/common/*.h) \
              libs/city_hash/city.h libs/lookup3_hash/lookup3.h \
              libs/string_buffer/string_buffer.h

PROC_CALLS_SRCS = src/tools/proc_calls.c $(COMMON_SRCS)
PLACE_CALLS_SRCS = src/tools/place_calls.c src/common/call_seqan.o $(COMMON_SRCS)
FILTER_SUBGRAPH_SRCS = src/tools/filter_subgraph.c $(COMMON_SRCS)
FILTER_READS_SRCS = src/tools/filter_reads.c $(COMMON_SRCS)
CTX_CALL_SRCS = src/tools/ctx_call.c $(COMMON_SRCS)

PROC_CALLS_HDRS = $(COMMON_HDRS)
PLACE_CALLS_HDRS = $(COMMON_HDRS)
FILTER_SUBGRAPH_HDRS = $(COMMON_HDRS)
FILTER_READS_HDRS = $(COMMON_HDRS)
CTX_CALL_HDRS = $(COMMON_HDRS)

MAXK_NUMCOLS = $(MAXK)_c$(NUM_COLS)

FILTERGRAPH_BIN=$(BIN)/filter_subgraph_$(MAXK_NUMCOLS)
FILTERREADS_BIN=$(BIN)/filter_reads_$(MAXK_NUMCOLS)

# DEPS are common dependencies that do not need to be re-built per target
DEPS=$(LIB_OBJS) $(BIN)

all: filter_reads filter_subgraph proc_calls place_calls ctx_call

proc_calls: bin/proc_calls
bin/proc_calls:  $(PROC_CALLS_SRCS) $(PROC_CALLS_HDRS) Makefile | $(DEPS)
	$(CC) -o bin/proc_calls $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(PROC_CALLS_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled proc_calls

place_calls: bin/place_calls
bin/place_calls: $(PLACE_CALLS_SRCS) $(PLACE_CALLS_HDRS) Makefile | $(DEPS)
	$(CC) -o bin/place_calls $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(PLACE_CALLS_SRCS) $(LIB_OBJS) $(LINK_LIBS) -lstdc++
	@echo Sucessfully compiled place_calls

filter_subgraph: $(FILTERGRAPH_BIN)
$(FILTERGRAPH_BIN): $(FILTER_SUBGRAPH_SRCS) $(FILTER_SUBGRAPH_HDRS) Makefile | $(DEPS)
	$(CC) -o $(FILTERGRAPH_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(FILTER_SUBGRAPH_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled filter_subgraph

filter_reads: $(FILTERREADS_BIN)
$(FILTERREADS_BIN): $(FILTER_READS_SRCS) $(FILTER_READS_HDRS) Makefile | $(DEPS)
	$(CC) -o $(FILTERREADS_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(FILTER_READS_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled filter_reads

ctx_call: bin/ctx_call
bin/ctx_call: $(CTX_CALL_SRCS) $(CTX_CALL_HDRS) Makefile | $(DEPS)
	$(CC) -o bin/ctx_call $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_CALL_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled ctx_call

src/common/call_seqan.o: src/common/call_seqan.cpp src/common/call_seqan.h
	$(CXX) -Wall -Wextra -I $(IDIR_SEQAN) -c -o src/common/call_seqan.o src/common/call_seqan.cpp

$(BIN):
	mkdir -p $(BIN)

$(TEMP_TEST_DIR):
	mkdir -p $(TEMP_TEST_DIR)

$(LIB_OBJS):
	cd libs; make

cunit:
	@if [[ ! -e '$(IDIR_CUNIT)' || ! -e '$(LDIR_CUNIT)' ]]; \
	then echo "Error: Cannot find CUnit"; exit 1; fi

clean:
	rm -rf $(BIN)/* $(TEMP_TEST_DIR)/* src/common/*.o

.PHONY: all clean cunit
.PHONY: cortex_var proc_calls place_calls filter_reads filter_subgraph ctx_call

