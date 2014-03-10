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
#  make test

# For release
# RECOMPILE=1

# Use bash as shell
SHELL := /bin/bash
CC ?= gcc

## Toggle Release Version
#
# RELEASE=1
# CTX_VERSION=0.1
# RECOMPILE=1
#
##

MAXK ?= 31
LOWK = $(word 1, $(sort $(MAXK) 31))
TEST_KMER=$(shell echo $$[(($(MAXK)+31)/32)*32 - 1])

ifneq ($(MAXK),$(TEST_KMER))
  KMERERROR=1
endif
ifneq ($(LOWK),31)
  KMERERROR=1
endif
ifdef KMERERROR
  $(error Invalid MAXK value '$(MAXK)'. Please choose from 31,63,95,..(32*n-1) [default: 31])
endif

MAX_KMER_SIZE=$(MAXK)
MIN_KMER_SIZE=$(shell echo $$[$(MAX_KMER_SIZE)-30] | sed 's/^1$$/3/g')

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
LIB_STRS=libs/string_buffer/string_buffer.o

MISC_SRCS=$(wildcard libs/misc/*.c)
MISC_OBJS=$(MISC_SRCS:.c=.o)

ifdef RELEASE
	LIB_MISC=$(MISC_SRCS)
else
	LIB_MISC=$(MISC_OBJS)
endif

ifdef LIB_PATH
	EXTRA_INCS := -I $(LIB_PATH) -L $(LIB_PATH)
endif

INCS=-I $(IDIR_MISC) -I $(IDIR_BITARR) -I $(IDIR_STRS) -I $(IDIR_HTS) \
     -I $(IDIR_SEQ) -I $(IDIR_ALIGN) -I $(IDIR_MSGPOOL) $(EXTRA_INCS)

# Library linking
LIB_OBJS=$(LIB_MISC) $(LIB_STRS) $(LIB_HTS) $(LIB_ALIGN)
LINK=-lpthread -lz -lm

CFLAGS = -std=c99 -Wall -Wextra
CPPFLAGS=$(HASH_KEY_FLAGS) -D_USESAM=1
KMERARGS=-DMIN_KMER_SIZE=$(MIN_KMER_SIZE) -DMAX_KMER_SIZE=$(MAX_KMER_SIZE)

# -fno-strict-aliasing
USEFUL_CFLAGS=-Wshadow -Wstrict-aliasing=2

IGNORE_CFLAGS=-Wno-cast-align -Wno-aggregate-return -Wno-conversion

OVERKILL_CFLAGS = -Winit-self -Wmissing-include-dirs \
                  -Wstrict-aliasing -Wdiv-by-zero \
                  -Wcast-qual -Wmissing-noreturn \
                  -Wwrite-strings -Wundef \
                  -Wshadow -Wfloat-equal -Wbad-function-cast \
                  -fstack-protector-all -D_FORTIFY_SOURCE=2

CLANG_CFLAGS=-fsanitize-undefined-trap-on-error -Wno-shorten-64-to-32

CFLAGS := $(CFLAGS) $(OVERKILL_CFLAGS) $(USEFUL_CFLAGS) $(IGNORE_CFLAGS)

PLATFORM := $(shell uname)
COMPILER := $(shell ($(CC) -v 2>&1) | tr A-Z a-z )

ifneq (,$(findstring clang,$(COMPILER)))
	CFLAGS := $(CFLAGS) $(CLANG_CFLAGS)
endif

# If not debugging, add optimisations and -DNDEBUG=1 to turn off assert() calls
ifdef DEBUG
	OPT = -O0 -g -ggdb -gdwarf-2 -g3
else
	# Could add -DNDEBUG=1 here to turn off asserts
	ifneq (,$(findstring gcc,$(COMPILER)))
		OPT = -O4
		TGTFLAGS = -fwhole-program
	else
		OPT = -O3
	endif

	ifneq (,$(findstring lto,$(COMPILER)))
		OPT := -flto $(OPT)
	endif
endif

CFLAGS := $(OPT) $(CFLAGS)

ifdef VERBOSE
	CPPFLAGS := $(CPPFLAGS) -DCTXVERBOSE=1
endif

ifndef RELEASE
	CPPFLAGS := $(CPPFLAGS) -DCTXCHECKS=1
endif

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

CMDS_OBJDIR=build/commands$(MAXK)
CMDS_SRCS=$(wildcard src/commands/*.c)
CMDS_HDRS=$(wildcard src/commands/*.h)
CMDS_FILES=$(notdir $(CMDS_SRCS))
CMDS_OBJS=$(addprefix $(CMDS_OBJDIR)/, $(CMDS_FILES:.c=.o))

TESTS_OBJDIR=build/tests$(MAXK)
TESTS_SRCS=$(wildcard src/tests/*.c)
TESTS_HDRS=$(wildcard src/tests/*.h)
TESTS_FILES=$(notdir $(TESTS_SRCS))
TESTS_OBJS=$(addprefix $(TESTS_OBJDIR)/, $(TESTS_FILES:.c=.o))

HDRS=$(CMDS_HDRS) $(KMER_HDRS) $(BASIC_HDRS)

DIRS=bin $(BASIC_OBJDIR) $(KMER_OBJDIR) $(TOOLS_OBJDIR) $(CMDS_OBJDIR) $(TESTS_OBJDIR)

# DEPS are kmer dependencies that do not need to be re-built per target
DEPS=Makefile $(DIRS) $(LIB_OBJS)

# RECOMPILE=1 to recompile all from source
ifdef RECOMPILE
	OBJS=$(CMDS_SRCS) $(TOOLS_SRCS) $(KMER_SRCS) $(BASIC_SRCS) $(LIB_OBJS)
	TESTS_OBJS=$(TESTS_SRCS)
else
	OBJS=$(CMDS_OBJS) $(TOOLS_OBJS) $(KMER_OBJS) $(BASIC_OBJS) $(LIB_OBJS)
endif

.DEFAULT_GOAL := ctx

all: ctx tests hashtest tables

# Run tests
test: tests
	./bin/tests$(MAXK)

# This Makefile mastery borrowed from htslib [https://github.com/samtools/htslib]
# If git repo, grab commit hash to use in version
# Force version.h to be remade if $(CTX_VERSION) has changed.
ifndef CTX_VERSION
ifneq "$(wildcard .git)" ""
#CTX_VERSION := $(shell git describe --always --dirty)
CTX_VERSION := $(shell git log --pretty=format:'%h' -n 1 --tags)
src/basic/version.h: $(if $(wildcard src/basic/version.h),$(if $(findstring "$(CTX_VERSION)",$(shell cat src/basic/version.h)),,force))
endif
endif

src/basic/version.h:
	echo '#define CTXVERSIONSTR "$(CTX_VERSION)"' > $@

$(BASIC_OBJDIR)/%.o: src/basic/%.c $(BASIC_HDRS) | $(DEPS)
	$(CC) -o $@ $(OBJFLAGS) $(CFLAGS) $(CPPFLAGS) -I src/basic/ $(INCS) -c $<

$(KMER_OBJDIR)/%.o: src/kmer/%.c $(BASIC_HDRS) $(KMER_HDRS) | $(DEPS)
	$(CC) -o $@ $(OBJFLAGS) $(CFLAGS) $(CPPFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ $(INCS) -c $<

$(TOOLS_OBJDIR)/%.o: src/tools/%.c $(BASIC_HDRS) $(KMER_HDRS) $(TOOLS_HDRS) | $(DEPS)
	$(CC) -o $@ $(OBJFLAGS) $(CFLAGS) $(CPPFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ -I src/tools/ $(INCS) -c $<

$(CMDS_OBJDIR)/%.o: src/commands/%.c $(BASIC_HDRS) $(KMER_HDRS) $(TOOLS_HDRS) $(CMDS_HDRS) | $(DEPS)
	$(CC) -o $@ $(OBJFLAGS) $(CFLAGS) $(CPPFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ -I src/tools/ -I src/commands/ $(INCS) -c $<

$(TESTS_OBJDIR)/%.o: src/tests/%.c $(BASIC_HDRS) $(KMER_HDRS) $(TOOLS_HDRS) | $(DEPS)
	$(CC) -o $@ $(OBJFLAGS) $(CFLAGS) $(CPPFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ -I src/tools/ $(INCS) -c $<

# Misc library code
libs/misc/%.o: libs/misc/%.c libs/misc/%.h
	$(CC) -o libs/misc/$*.o $(OBJFLAGS) $(CFLAGS) $(CPPFLAGS) -c libs/misc/$*.c

ctx: bin/ctx$(MAXK)
bin/ctx$(MAXK): src/main/ctx.c $(OBJS) $(HDRS) | bin
	$(CC) -o $@ $(TGTFLAGS) $(CFLAGS) $(CPPFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ -I src/tools/ -I src/commands/ $(INCS) src/main/ctx.c $(OBJS) $(LINK)

tests: bin/tests$(MAXK)
bin/tests$(MAXK): src/main/tests.c $(TESTS_OBJS) $(OBJS) $(HDRS) | bin
	$(CC) -o $@ $(TGTFLAGS) $(CFLAGS) $(CPPFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ -I src/tools/ -I src/commands/ -I src/tests/ $(INCS) src/main/tests.c $(TESTS_OBJS) $(OBJS) $(LINK)

hashtest: bin/hashtest$(MAXK)
bin/hashtest$(MAXK): src/main/hashtest.c $(OBJS) $(HDRS) | bin
	$(CC) -o $@ $(TGTFLAGS) $(CFLAGS) $(CPPFLAGS) $(KMERARGS) -I src/basic/ -I src/kmer/ -I src/tools/ $(INCS) src/main/hashtest.c $(OBJS) $(LINK)

tables: bin/tables
bin/tables: src/main/tables.c | bin
	$(CC) -o $@ $(TGTFLAGS) $(CFLAGS) $<

# directories
$(DIRS):
	mkdir -p $@

# libraries
# This triggers the compiling of library dependencies for first install
$(LIB_OBJS): libs/string_buffer/string_buffer.h
libs/string_buffer/string_buffer.h:
	cd libs; make

clean:
	rm -rf bin build lib/misc/*.o

force:

.PHONY: all clean ctx test force
