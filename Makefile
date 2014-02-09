ifndef CC
  CC = gcc
endif

#Arguments (all are optional):
# MAXK=31
# DEBUG=1       (compile with debug symbols)
# VERBOSE=1     (print all the things!)
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

CTX_VERSION=0.0.1

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
	HASH_KEY_FLAGS=-DUSE_CITY_HASH=1
endif

# Library paths
# IDIR_GSL_HEADERS=libs/gsl-1.16
IDIR_HTS=libs/htslib/htslib
IDIR_STRS=libs/string_buffer
IDIR_SEQ=libs/seq_file
IDIR_ALIGN=libs/seq-align/src
IDIR_BITARR=libs/bit_array
IDIR_MSGPOOL=libs/msg-pool
IDIR_MISC=libs/misc

# LIB_GSL=libs/gsl-1.16/.libs/libgsl.a
LIB_HTS=libs/htslib/libhts.a
LIB_ALIGN=libs/seq-align/src/libalign.a
# LIB_STRS=libs/string_buffer/libstrbuf.a
LIB_STRS=libs/string_buffer/string_buffer.c
LIB_MISC=$(wildcard libs/misc/*.o)

ifdef LIB_PATH
	EXTRA_INCS := -I $(LIB_PATH) -L $(LIB_PATH)
endif

INCS=-I $(IDIR_MISC) -I $(IDIR_BITARR) -I $(IDIR_STRS) -I $(IDIR_HTS) \
     -I $(IDIR_SEQ) -I $(IDIR_ALIGN) -I $(IDIR_MSGPOOL) $(EXTRA_INCS)

# INCS=-I src/basic/ -I src/kmer/ -I src/tools/ $(INCS_EXTERNAL)
# -I $(IDIR_GSL_HEADERS)

# Library linking
LIB_OBJS=$(LIB_MISC) $(LIB_STRS) $(LIB_HTS) $(LIB_ALIGN)
LINK=-lpthread -lz -lm
# $(LIB_GSL)

# -Winit-self -Wmissing-include-dirs
# -Wstrict-aliasing -Wdiv-by-zero -Wunreachable-code
# -Wcast-qual -Wcast-align -Wmissing-noreturn
# -Wwrite-strings -Waggregate-return -Wundef
# -Wshadow -Wconversion -Wshorten-64-to-32 -Woverlength-strings
# -Wenum-compare -Wlogical-op -Wfloat-equal -Wbad-function-cast

# Optimisations tags for testing
# -Wstack-protector -fstack-protector
# -fsanitize=thread

# -fno-strict-aliasing
USEFUL_CFLAGS=-Wshadow -Wstrict-aliasing=2

IGNORE_CFLAGS=-Wno-cast-align -Wno-shorten-64-to-32 \
              -Wno-aggregate-return -Wno-conversion

OVERKILL_CFLAGS = -Winit-self -Wmissing-include-dirs \
                 -Wstrict-aliasing -Wdiv-by-zero -Wunreachable-code \
                 -Wcast-qual -Wmissing-noreturn \
                 -Wwrite-strings -Wundef \
                 -Wshadow -Woverlength-strings \
                 -Wenum-compare -Wfloat-equal -Wbad-function-cast \
                 -fstack-protector-all -D_FORTIFY_SOURCE=2

#CLANG_ONLY=-fsanitize-undefined-trap-on-error

# If not debugging, add optimisations and -DNDEBUG=1 to turn off assert() calls
ifdef DEBUG
	OPT = -O0 $(OVERKILL_CFLAGS) $(USEFUL_CFLAGS) $(IGNORE_CFLAGS)
	DEBUG_ARGS = -g -ggdb -gdwarf-2 -g3
else
	ifdef PROFILE
		#-DNDEBUG=1
		OPT = -O4 $(OVERKILL_CFLAGS) $(USEFUL_CFLAGS) $(IGNORE_CFLAGS)
		DEBUG_ARGS = -g -ggdb
	else
		#-DNDEBUG=1
		OPT = -O4 $(OVERKILL_CFLAGS) $(USEFUL_CFLAGS) $(IGNORE_CFLAGS)
		DEBUG_ARGS = 
	endif
endif

ifdef VERBOSE
	OUTPUT_ARGS=-DCTXVERBOSE=1
endif

CFLAGS = -std=c99 -Wall -Wextra $(OPT) $(DEBUG_ARGS) $(OUTPUT_ARGS) \
          $(HASH_KEY_FLAGS) -D_USESAM=1 -DCTXCHECKS=1

KMERARGS=-DMAX_KMER_SIZE=$(MAX_KMER_SIZE) -DMIN_KMER_SIZE=$(MIN_KMER_SIZE) \
         -DNUM_BKMER_WORDS=$(NUM_BKMER_WORDS)

# basic objects compile without MAXK
# kmer and tool objects require MAXK
BASIC_OBJDIR=build/basic
BASIC_SRCS=$(wildcard src/basic/*.c)
BASIC_HDRS=$(wildcard src/basic/*.h) src/basic/version.h
BASIC_FILES=$(notdir $(BASIC_SRCS))
BASIC_OBJS=$(addprefix $(BASIC_OBJDIR)/, $(BASIC_FILES:.c=.o))

KMER_OBJDIR=build/kmer$(MAXK)
KMER_SRCS=$(wildcard src/kmer/*.c)
KMER_HDRS=$(wildcard src/kmer/*.h)
KMER_FILES=$(notdir $(KMER_SRCS))
KMER_OBJS=$(addprefix $(KMER_OBJDIR)/, $(KMER_FILES:.c=.o))

TOOLS_OBJDIR=build/tools$(MAXK)
TOOLS_SRCS=$(wildcard src/tools/*.c)
TOOLS_HDRS=$(wildcard src/tools/*.h)
TOOLS_FILES=$(notdir $(TOOLS_SRCS))
TOOLS_OBJS=$(addprefix $(TOOLS_OBJDIR)/, $(TOOLS_FILES:.c=.o))

HDRS=$(TOOLS_HDRS) $(KMER_HDRS) $(BASIC_HDRS)

DIRS=bin $(BASIC_OBJDIR) $(KMER_OBJDIR) $(TOOLS_OBJDIR)

# DEPS are kmer dependencies that do not need to be re-built per target
DEPS=Makefile $(DIRS) $(LIB_OBJS)

# RECOMPILE=1 to recompile all from source
ifdef RECOMPILE
	OBJS=$(TOOLS_SRCS) $(KMER_SRCS) $(BASIC_SRCS) $(LIB_OBJS)
else
	OBJS=$(TOOLS_OBJS) $(KMER_OBJS) $(BASIC_OBJS) $(LIB_OBJS)
endif

.DEFAULT_GOAL := ctx

all: ctx tests hashtest

# This Makefile mastery borrowed from htslib [https://github.com/samtools/htslib]
# If git repo, grab commit hash to use in version
# Force version.h to be remade if $(CTX_VERSION) has changed.
ifneq "$(wildcard .git)" ""
#CTX_VERSION := $(shell git describe --always --dirty)
CTX_VERSION := $(shell git log --pretty=format:'%h' -n 1 --tags)
src/basic/version.h: $(if $(wildcard src/basic/version.h),$(if $(findstring "$(CTX_VERSION)",$(shell cat src/basic/version.h)),,force))
endif

src/basic/version.h:
	echo '#define CTXVERSIONSTR "$(CTX_VERSION)"' > $@

$(BASIC_OBJDIR)/%.o: src/basic/%.c $(BASIC_HDRS) | $(DEPS)
	$(CC) $(CFLAGS) -I src/basic/ $(INCS) -c $< -o $@

$(KMER_OBJDIR)/%.o: src/kmer/%.c $(BASIC_HDRS) $(KMER_HDRS) | $(DEPS)
	$(CC) $(CFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ $(INCS) -c $< -o $@

$(TOOLS_OBJDIR)/%.o: src/tools/%.c $(BASIC_HDRS) $(KMER_HDRS) $(TOOLS_HDRS) | $(DEPS)
	$(CC) $(CFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ -I src/tools/ $(INCS) -c $< -o $@

ctx: bin/ctx$(MAXK)
bin/ctx$(MAXK): src/main/ctx.c $(OBJS) $(HDRS) | bin
	$(CC) -o $@ $(CFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ -I src/tools/ src/main/ctx.c $(INCS) $(OBJS) $(LINK)

tests: bin/tests$(MAXK)
bin/tests$(MAXK): src/main/tests.c src/tests/* $(OBJS) $(HDRS) | bin
	$(CC) -o $@ $(CFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ -I src/tools/ -I src/tests/ src/tests/*.c src/main/tests.c $(INCS) $(OBJS) $(LINK)

hashtest: bin/hashtest$(MAXK)
bin/hashtest$(MAXK): src/main/hashtest.c $(OBJS) $(HDRS) | bin
	$(CC) -o $@ $(CFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ src/main/hashtest.c $(INCS) $(OBJS) $(LINK)

# directories
$(DIRS):
	mkdir -p $@

# libraries
# This triggers the compiling of library dependencies for first install
$(LIB_OBJS): libs/string_buffer/string_buffer.h
libs/string_buffer/string_buffer.h:
	cd libs; make

clean:
	rm -rf bin build

force:

.PHONY: all clean ctx force
