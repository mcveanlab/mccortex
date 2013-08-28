ifndef CC
  CC = gcc
endif

#Arguments (all are optional):
# MAXK=31
# DEBUG=1
# CITY_HASH=1   (use CityHash hash function)
# RECOMPILE=1   (recompile all from source)

# Resolve some issues linking libz:
# e.g. for WTCHG cluster3
# 1) pass LIB_PATH=/usr/local/lib/ to compile on WTCHG cluster3
# 2) set LD_LIBRARY_PATH=/usr/local/lib/:$(LD_LIBRARY_PATH) before running

# Targets:
#  make
#  make clean
#  make all

# For release
# RECOMPILE=1

# Use bash as shell
SHELL := /bin/bash

NUM_BKMER_WORDS=$(shell echo $$((($(MAXK)+31)/32)))

ifndef MAXK
  NUM_BKMER_WORDS = 1
  MAXK = 31
else
  NUM_BKMER_WORDS=$(shell echo $$((($(MAXK)+31)/32)))
endif

ifeq ($(NUM_BKMER_WORDS),0)
  $(error Invalid MAXK value '$(MAXK)'. Please choose from 31,63,95,..(32*n-1) [default: 31])
endif

MAX_KMER_SIZE=$(shell echo $$(($(NUM_BKMER_WORDS)*32-1)))
MIN_KMER_SIZE=$(shell echo $$(($(MAX_KMER_SIZE)-30)))

ifeq ($(MIN_KMER_SIZE),1)
	MIN_KMER_SIZE=3
endif

# Use City hash instead of lookup3?
ifdef CITY_HASH
	HASH_KEY_FLAGS = -DUSE_CITY_HASH=1
endif

# Library paths
IDIR_GSL_HEADERS = libs/gsl-1.16
IDIR_HTS = libs/htslib/htslib
IDIR_STRS = libs/string_buffer
IDIR_SEQ = libs/seq_file
IDIR_ALIGN = libs/seq-align/src
IDIR_SEQAN = libs/seqan-library-1.4.0/include/
IDIR_HASH = libs/hash_functions

LIB_GSL=libs/gsl-1.16/.libs/libgsl.a
LIB_HTS=libs/htslib/htslib/libhts.a
LIB_ALIGN=libs/seq-align/src/libalign.a
# LIB_STRS=libs/string_buffer/libstrbuf.a
LIB_STRS=libs/string_buffer/string_buffer.c

ifdef LIB_PATH
	EXTRA_INCS := -I $(LIB_PATH) -L $(LIB_PATH)
endif

INCS=-I src/basic/ -I src/common/ -I src/tools/ \
     -I $(IDIR_HASH) -I $(IDIR_STRS) -I $(IDIR_HTS) \
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
	ifdef PROFILE
	#-DNDEBUG=1
		OPT = -O3
		DEBUG_ARGS = -g -ggdb
	else
	#-DNDEBUG=1
		OPT = -O3
		DEBUG_ARGS = 
	endif
endif

CFLAGS = -std=c99 -Wall -Wextra $(OPT) $(DEBUG_ARGS) \
         -DMAX_KMER_SIZE=$(MAX_KMER_SIZE) -DMIN_KMER_SIZE=$(MIN_KMER_SIZE) \
         -DNUM_BKMER_WORDS=$(NUM_BKMER_WORDS) $(HASH_KEY_FLAGS)

# basic objects compile without MAXK
# common and tool objects require MAXK
BASIC_OBJDIR=build/basic
COMMON_OBJDIR=build/common$(MAXK)
TOOLS_OBJDIR=build/tools$(MAXK)

BASIC_SRCS=$(wildcard src/basic/*.c)
BASIC_HDRS=$(wildcard src/basic/*.h)
BASIC_FILES=$(notdir $(BASIC_SRCS)) call_seqan.o
BASIC_OBJS=$(addprefix $(BASIC_OBJDIR)/, $(BASIC_FILES:.c=.o))

COMMON_SRCS=$(wildcard src/common/*.c)
COMMON_HDRS=$(wildcard src/common/*.h)
COMMON_FILES=$(notdir $(COMMON_SRCS))
COMMON_OBJS=$(addprefix $(COMMON_OBJDIR)/, $(COMMON_FILES:.c=.o))

TOOLS_SRCS=$(wildcard src/tools/*.c)
TOOLS_HDRS=$(wildcard src/tools/*.h)
TOOLS_FILES=$(notdir $(TOOLS_SRCS))
TOOLS_OBJS=$(addprefix $(TOOLS_OBJDIR)/, $(TOOLS_FILES:.c=.o))

HDRS=$(TOOLS_HDRS) $(COMMON_HDRS) $(BASIC_HDRS)

# DEPS are common dependencies that do not need to be re-built per target
DEPS=Makefile bin $(BASIC_OBJDIR) $(COMMON_OBJDIR) $(TOOLS_OBJDIR) $(LIB_OBJS)

# RECOMPILE=1 to recompile all from source
ifdef RECOMPILE
	OBJS=$(TOOLS_SRCS) $(COMMON_SRCS) $(BASIC_SRCS) build/basic/call_seqan.o $(LIB_OBJS)
else
	OBJS=$(TOOLS_OBJS) $(COMMON_OBJS) $(BASIC_OBJS) $(LIB_OBJS)
endif

all: ctx

build/basic/call_seqan.o: src/basic/call_seqan.cpp src/basic/call_seqan.h | $(DEPS)
	$(CXX) -Wall -Wextra -I $(IDIR_SEQAN) -c src/basic/call_seqan.cpp -o $@

$(BASIC_OBJDIR)/%.o: src/basic/%.c $(BASIC_HDRS) | $(DEPS)
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

$(COMMON_OBJDIR)/%.o: src/common/%.c $(BASIC_HDRS) $(COMMON_HDRS) | $(DEPS)
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

$(TOOLS_OBJDIR)/%.o: src/tools/%.c $(BASIC_HDRS) $(COMMON_HDRS) $(TOOLS_HDRS) | $(DEPS)
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

ctx: bin/ctx$(MAXK)
bin/ctx$(MAXK): src/main/ctx.c $(OBJS) $(HDRS) | bin
	$(CC) -o $@ $(CFLAGS) $(INCS) src/main/ctx.c $(OBJS) $(LINK)

traversal: bin/traversal$(MAXK)
bin/traversal$(MAXK): src/main/traversal.c $(OBJS) $(HDRS) | bin
	$(CC) -o $@ $(CFLAGS) $(INCS) src/main/traversal.c $(OBJS) $(LINK)

hashtest: bin/hashtest$(MAXK)
bin/hashtest$(MAXK): src/main/hashtest.c $(OBJS) $(HDRS) | bin
	$(CC) -o $@ $(CFLAGS) $(INCS) src/main/hashtest.c $(OBJS) $(LINK)

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

.PHONY: all clean ctx traversal
