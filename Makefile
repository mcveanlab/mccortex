ifndef CC
  CC = gcc
endif

#Arguments (all are optional):
# MAXK=31
# NUM_COLS=12
# DEBUG=1
# CITY_HASH=1   (use CityHash hash function)

# e.g. for WTCHG cluster3
# 1) pass LIB_PATH=/usr/local/lib/ to compile on WTCHG cluster3
# 2) set LD_LIBRARY_PATH=/usr/local/lib/:$(LD_LIBRARY_PATH) before running

# Targets:
#  make
#  make clean
#  make all

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

# Library paths
IDIR_GSL_HEADERS = libs/gsl-1.15
IDIR_HTS = libs/htslib/htslib
IDIR_STRS = libs/string_buffer
IDIR_SEQ = libs/seq_file
IDIR_ALIGN = libs/seq-align/src
IDIR_SEQAN = libs/seqan-library-1.4.0/include/
IDIR_CITY = libs/city_hash
IDIR_LOOKUP3 = libs/lookup3_hash

LIB_GSL = libs/gsl-1.15/.libs/libgsl.a
LIB_HTS = libs/htslib/htslib/libhts.a
LIB_ALIGN = libs/seq-align/src/libalign.a

INCLUDES = -I src/common/ \
           -I $(IDIR_CITY) -I $(IDIR_LOOKUP3) \
           -I $(IDIR_GSL_HEADERS) -I $(IDIR_HTS) \
           -I $(IDIR_ALIGN) -I $(IDIR_SEQ) -I $(IDIR_STRS)

WARNS = -Wall -Wextra -Winit-self -Wmissing-include-dirs \
        -Wstrict-aliasing -Wdiv-by-zero \
        -Wcast-qual -Wcast-align -Wmissing-noreturn \
        -Wwrite-strings -Waggregate-return -Wundef

# -Wshadow -Wconversion -Wshorten-64-to-32 -Woverlength-strings -Wunreachable-code
# -Wenum-compare -Wlogical-op -Wfloat-equal -Wbad-function-cast

CFLAGS = -std=c99 $(WARNS) \
         -DMAX_KMER_SIZE=$(MAX_KMER_SIZE) -DMIN_KMER_SIZE=$(MIN_KMER_SIZE) \
         -DNUM_BITFIELDS_IN_BKMER=$(BITFIELDS) -DNUM_OF_COLOURS=$(NUM_COLS) \
         $(HASH_KEY_FLAGS)

# Optimisations tags for testing
OPT_TESTS = -O0 -Wstack-protector -fstack-protector

# If not debugging, add optimisations and -DNDEBUG=1 to turn off assert() calls
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

CTX_SUBGRAPH_SRCS = src/tools/ctx_subgraph.c $(COMMON_SRCS)
CTX_READS_SRCS = src/tools/ctx_reads.c $(COMMON_SRCS)
CTX_CALL_SRCS = src/tools/ctx_call.c $(COMMON_SRCS)
CTX_UNIQUE_SRCS = src/tools/ctx_unique.c $(COMMON_SRCS)
CTX_PLACE_SRCS = src/tools/ctx_place.c src/common/call_seqan.o $(COMMON_SRCS)

CTX_SUBGRAPH_HDRS = $(COMMON_HDRS)
CTX_READS_HDRS = $(COMMON_HDRS)
CTX_CALL_HDRS = $(COMMON_HDRS)
CTX_UNIQUE_HDRS = $(COMMON_HDRS)
CTX_PLACE_HDRS = $(COMMON_HDRS)

MAXK_NUMCOLS = k$(MAXK)_c$(NUM_COLS)

CTX_SUBGRAPH_BIN=bin/ctx_subgraph_$(MAXK_NUMCOLS)
CTX_READS_BIN=bin/ctx_reads_$(MAXK_NUMCOLS)
CTX_CALL_BIN=bin/ctx_call_$(MAXK_NUMCOLS)

# DEPS are common dependencies that do not need to be re-built per target
DEPS=$(LIB_OBJS) bin

all: ctx_reads ctx_subgraph ctx_unique ctx_place ctx_call

ctx_subgraph: $(CTX_SUBGRAPH_BIN)
$(CTX_SUBGRAPH_BIN): $(CTX_SUBGRAPH_SRCS) $(CTX_SUBGRAPH_HDRS) Makefile | $(DEPS)
	$(CC) -o $(CTX_SUBGRAPH_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_SUBGRAPH_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled $(CTX_SUBGRAPH_BIN)

ctx_reads: $(CTX_READS_BIN)
$(CTX_READS_BIN): $(CTX_READS_SRCS) $(CTX_READS_HDRS) Makefile | $(DEPS)
	$(CC) -o $(CTX_READS_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_READS_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled $(CTX_READS_BIN)

ctx_call: $(CTX_CALL_BIN)
$(CTX_CALL_BIN): $(CTX_CALL_SRCS) $(CTX_CALL_HDRS) Makefile | $(DEPS)
	$(CC) -o $(CTX_CALL_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_CALL_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled $(CTX_CALL_BIN)

ctx_unique: bin/ctx_unique
bin/ctx_unique:  $(CTX_UNIQUE_SRCS) $(CTX_UNIQUE_HDRS) Makefile | $(DEPS)
	$(CC) -o bin/ctx_unique $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_UNIQUE_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled ctx_unique

ctx_place: bin/ctx_place
bin/ctx_place: $(CTX_PLACE_SRCS) $(CTX_PLACE_HDRS) Makefile | $(DEPS)
	$(CC) -o bin/ctx_place $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_PLACE_SRCS) $(LIB_OBJS) $(LINK_LIBS) -lstdc++
	@echo Sucessfully compiled ctx_place

src/common/call_seqan.o: src/common/call_seqan.cpp src/common/call_seqan.h
	$(CXX) -Wall -Wextra -I $(IDIR_SEQAN) -c -o src/common/call_seqan.o src/common/call_seqan.cpp

bin:
	mkdir -p bin

$(TEMP_TEST_DIR):
	mkdir -p $(TEMP_TEST_DIR)

$(LIB_OBJS):
	cd libs; make

cunit:
	@if [[ ! -e '$(IDIR_CUNIT)' || ! -e '$(LDIR_CUNIT)' ]]; \
	then echo "Error: Cannot find CUnit"; exit 1; fi

clean:
	rm -rf bin/* src/common/*.o

.PHONY: all clean cunit
.PHONY: ctx_call ctx_unique ctx_place ctx_reads ctx_subgraph

