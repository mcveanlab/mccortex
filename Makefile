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
IDIR_HASH = libs/hash_functions

LIB_GSL=libs/gsl-1.15/.libs/libgsl.a
LIB_HTS=libs/htslib/htslib/libhts.a
LIB_ALIGN=libs/seq-align/src/libalign.a
LIB_STRS=libs/string_buffer/libstrbuf.a

# Resolve some issues linking libz:
# 1) pass LIB_PATH=/usr/local/lib/ to compile on WTCHG cluster3
# 2) also set LD_LIBRARY_PATH=/usr/local/lib/:$(LD_LIBRARY_PATH)
ifdef LIB_PATH
	EXTRA_INCS := -I $(LIB_PATH) -L $(LIB_PATH)
endif

INCS=-I src/basic/ -I src/common/ -I $(IDIR_HASH) -I $(IDIR_STRS) -I $(IDIR_HTS) \
     -I $(IDIR_SEQ) -I $(IDIR_ALIGN) -I $(IDIR_GSL_HEADERS) $(EXTRA_INCS)

# Library linking
LIB_OBJS=$(LIB_GSL) $(LIB_STRS) $(LIB_HTS) $(LIB_ALIGN) $(wildcard libs/hash_functions/*.o)
LINK=-lstdc++ -lpthread -lz -lm

# -Winit-self -Wmissing-include-dirs
# -Wstrict-aliasing -Wdiv-by-zero -Wunreachable-code
# -Wcast-qual -Wcast-align -Wmissing-noreturn
# -Wwrite-strings -Waggregate-return -Wundef
# -Wshadow -Wconversion -Wshorten-64-to-32 -Woverlength-strings
# -Wenum-compare -Wlogical-op -Wfloat-equal -Wbad-function-cast

# Optimisations tags for testing
OPT_TESTS = -O0 -Wstack-protector -fstack-protector

# If not debugging, add optimisations and -DNDEBUG=1 to turn off assert() calls
ifdef DEBUG
	OPT = $(OPT_TESTS)
	DEBUG_ARGS = -g -ggdb -DDEBUG=1
else
	OPT = -O3 -DNDEBUG=1
	DEBUG_ARGS = 
endif

CFLAGS = -std=c99 -Wall -Wextra \
         -DMAX_KMER_SIZE=$(MAX_KMER_SIZE) -DMIN_KMER_SIZE=$(MIN_KMER_SIZE) \
         -DNUM_BITFIELDS_IN_BKMER=$(BITFIELDS) $(HASH_KEY_FLAGS) $(OPT) $(DEBUG_ARGS)

# basic objects compile without MAXK
# common and tool objects require MAXK
BASIC_OBJDIR=build/basic
COMMON_OBJDIR=build/common$(MAXK)
TOOLS_OBJDIR=build/tools$(MAXK)

BASIC_FILES=$(notdir $(wildcard src/basic/*.c)) call_seqan.o
BASIC_OBJS=$(addprefix $(BASIC_OBJDIR)/, $(BASIC_FILES:.c=.o))
BASIC_HDRS=$(wildcard src/basic/*.h)

COMMON_FILES=$(notdir $(wildcard src/common/*.c))
COMMON_OBJS=$(addprefix $(COMMON_OBJDIR)/, $(COMMON_FILES:.c=.o))
COMMON_HDRS=$(wildcard src/common/*.h)

TOOLS_FILES=$(notdir $(wildcard src/tools/*.c))
TOOLS_OBJS=$(addprefix $(TOOLS_OBJDIR)/, $(TOOLS_FILES:.c=.o))
TOOLS_HDRS=$(wildcard src/toools/*.h)

# DEPS are common dependencies that do not need to be re-built per target
DEPS=Makefile bin build/basic $(COMMON_OBJDIR) $(LIB_OBJS)

TOOLS=ctx_build ctx_clean ctx_reads ctx_subgraph ctx_intersect ctx_join \
      ctx_thread ctp_view ctx_call ctx_covg ctx_extend ctx_contigs ctx_diverge \
      ctx_unique ctx_place

# all: $(TOOLS)
all: ctx

.SUFFIXES: .c .o _k$(MAXK)
.SECONDARY:

build/basic/call_seqan.o: src/basic/call_seqan.cpp src/basic/call_seqan.h | $(DEPS)
	$(CXX) -Wall -Wextra -I $(IDIR_SEQAN) -c src/basic/call_seqan.cpp -o $@

$(BASIC_OBJDIR)/%.o: src/basic/%.c $(BASIC_HDRS) | build/basic
	$(CC) $(DEBUG_ARGS) $(CFLAGS) $(INCS) -c $< -o $@

$(COMMON_OBJDIR)/%.o: src/common/%.c $(BASIC_HDRS) $(COMMON_HDRS) | $(COMMON_OBJDIR)
	$(CC) $(DEBUG_ARGS) $(CFLAGS) $(INCS) -c $< -o $@

$(TOOLS_OBJDIR)/%.o: src/tools/%.c $(BASIC_HDRS) $(COMMON_HDRS) $(TOOLS_HDRS) | $(TOOLS_OBJDIR)
	$(CC) $(DEBUG_ARGS) $(CFLAGS) $(INCS) -c $< -o $@

ctx: bin/ctx$(MAXK)
bin/ctx$(MAXK): $(TOOLS_OBJS) $(COMMON_OBJS) $(BASIC_OBJS) $(LIB_OBJS)
	$(CC) -o $@ $(CFLAGS) $(INCS) $(TOOLS_OBJS) $(COMMON_OBJS) $(BASIC_OBJS) $(LIB_OBJS) $(LINK)

# directories
bin:
	mkdir -p bin

$(BASIC_OBJDIR):
	mkdir -p $(BASIC_OBJDIR)

$(COMMON_OBJDIR):
	mkdir -p $(COMMON_OBJDIR)

$(TOOLS_OBJDIR):
	mkdir -p $(TOOLS_OBJDIR)

# libraries
$(LIB_OBJS):
libs/string_buffer/string_buffer.h:
	cd libs; make

clean:
	rm -rf bin build

.PHONY: all clean $(TOOLS)
