ifndef CC
  CC = gcc
endif

#Arguments (all are optional):
# MAXK=31
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

# Resolve some issues linking libz:
# 1) pass LIB_PATH=/usr/local/lib/ to compile on WTCHG cluster3
# 2) also set LD_LIBRARY_PATH=/usr/local/lib/:$(LD_LIBRARY_PATH)
ifdef LIB_PATH
	EXTRA_INCLUDES := -I $(LIB_PATH) -L $(LIB_PATH)
endif

INCLUDES = $(EXTRA_INCLUDES) -I src/common/ \
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
         -DNUM_BITFIELDS_IN_BKMER=$(BITFIELDS) $(HASH_KEY_FLAGS)

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

# Library linking
LIB_OBJS = $(LIB_GSL) $(LIB_HTS) $(LIB_ALIGN)
LINK_LIBS = -lpthread -lz -lm $(DEBUG_LIBS)

COMMON_SRCS = $(wildcard src/common/*.c) \
              libs/city_hash/city.c libs/lookup3_hash/lookup3.c \
              libs/string_buffer/string_buffer.c

COMMON_HDRS = $(wildcard src/common/*.h) \
              libs/city_hash/city.h libs/lookup3_hash/lookup3.h \
              libs/string_buffer/string_buffer.h

CTX_BUILD_SRCS = src/tools/ctx_build.c $(COMMON_SRCS)
CTX_CLEAN_SRCS = src/tools/ctx_clean.c $(COMMON_SRCS)
CTX_SUBGRAPH_SRCS = src/tools/ctx_subgraph.c $(COMMON_SRCS)
CTX_READS_SRCS = src/tools/ctx_reads.c $(COMMON_SRCS)
CTX_THREAD_SRCS = src/tools/ctx_thread.c $(COMMON_SRCS)
CTP_VIEW_SRCS = src/tools/ctP_VIEW.c $(COMMON_SRCS)
CTX_CALL_SRCS = src/tools/ctx_call.c $(COMMON_SRCS)
CTX_UNIQUE_SRCS = src/tools/ctx_unique.c $(COMMON_SRCS)
CTX_PLACE_SRCS = src/tools/ctx_place.c src/common/call_seqan.o $(COMMON_SRCS)
CTX_COVG_SRCS = src/tools/ctx_covg.c $(COMMON_SRCS)

CTX_BUILD_BIN=bin/ctx_build_k$(MAXK)
CTX_CLEAN_BIN=bin/ctx_clean_k$(MAXK)
CTX_SUBGRAPH_BIN=bin/ctx_subgraph_k$(MAXK)
CTX_READS_BIN=bin/ctx_reads_k$(MAXK)
CTX_THREAD_BIN=bin/ctx_thread_k$(MAXK)
CTP_VIEW_BIN=bin/ctp_view_k$(MAXK)
CTX_CALL_BIN=bin/ctx_call_k$(MAXK)
CTX_COVG_BIN=bin/ctx_covg_k$(MAXK)

# DEPS are common dependencies that do not need to be re-built per target
DEPS=$(LIB_OBJS) bin

TOOLS=ctx_build ctx_clean ctx_reads ctx_subgraph \
      ctx_thread ctp_view ctx_call ctx_unique ctx_place ctx_covg

all: $(TOOLS)

ctx_build: $(CTX_BUILD_BIN)
$(CTX_BUILD_BIN): $(CTX_BUILD_SRCS) $(COMMON_HDRS) Makefile | $(DEPS)
	$(CC) -o $(CTX_BUILD_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_BUILD_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled $(CTX_BUILD_BIN)

ctx_clean: $(CTX_CLEAN_BIN)
$(CTX_CLEAN_BIN): $(CTX_CLEAN_SRCS) $(COMMON_HDRS) Makefile | $(DEPS)
	$(CC) -o $(CTX_CLEAN_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_CLEAN_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled $(CTX_CLEAN_BIN)

ctx_subgraph: $(CTX_SUBGRAPH_BIN)
$(CTX_SUBGRAPH_BIN): $(CTX_SUBGRAPH_SRCS) $(COMMON_HDRS) Makefile | $(DEPS)
	$(CC) -o $(CTX_SUBGRAPH_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_SUBGRAPH_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled $(CTX_SUBGRAPH_BIN)

ctx_reads: $(CTX_READS_BIN)
$(CTX_READS_BIN): $(CTX_READS_SRCS) $(COMMON_HDRS) Makefile | $(DEPS)
	$(CC) -o $(CTX_READS_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_READS_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled $(CTX_READS_BIN)

ctx_thread: $(CTX_THREAD_BIN)
$(CTX_THREAD_BIN): $(CTX_THREAD_SRCS) $(COMMON_HDRS) Makefile | $(DEPS)
	$(CC) -o $(CTX_THREAD_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_THREAD_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled $(CTX_THREAD_BIN)

ctp_view: $(CTP_VIEW_BIN)
$(CTP_VIEW_BIN): $(CTP_VIEW_SRCS) $(COMMON_HDRS) Makefile | $(DEPS)
	$(CC) -o $(CTP_VIEW_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTP_VIEW_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled $(CTP_VIEW_BIN)

ctx_call: $(CTX_CALL_BIN)
$(CTX_CALL_BIN): $(CTX_CALL_SRCS) $(COMMON_HDRS) Makefile | $(DEPS)
	$(CC) -o $(CTX_CALL_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_CALL_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled $(CTX_CALL_BIN)

ctx_unique: bin/ctx_unique
bin/ctx_unique:  $(CTX_UNIQUE_SRCS) $(COMMON_HDRS) Makefile | $(DEPS)
	$(CC) -o bin/ctx_unique $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_UNIQUE_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled ctx_unique

ctx_covg: $(CTX_COVG_BIN)
$(CTX_COVG_BIN): $(CTX_COVG_SRCS) $(COMMON_HDRS) Makefile | $(DEPS)
	$(CC) -o $(CTX_COVG_BIN) $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_COVG_SRCS) $(LIB_OBJS) $(LINK_LIBS)
	@echo Sucessfully compiled $(CTX_COVG_BIN)

ctx_place: bin/ctx_place
bin/ctx_place: $(CTX_PLACE_SRCS) $(COMMON_HDRS) Makefile | $(DEPS)
	$(CC) -o bin/ctx_place $(DEBUG_ARGS) $(OPT) $(CFLAGS) $(INCLUDES) $(CTX_PLACE_SRCS) $(LIB_OBJS) $(LINK_LIBS) -lstdc++
	@echo Sucessfully compiled ctx_place

src/common/call_seqan.o: src/common/call_seqan.cpp src/common/call_seqan.h
	$(CXX) -Wall -Wextra -I $(IDIR_SEQAN) -c -o src/common/call_seqan.o src/common/call_seqan.cpp

$(LIB_OBJS):
libs/string_buffer/string_buffer.h:
	cd libs; make

bin:
	mkdir -p bin

$(TEMP_TEST_DIR):
	mkdir -p $(TEMP_TEST_DIR)

clean:
	rm -rf bin/* src/common/*.o

.PHONY: all clean $(TOOLS)
