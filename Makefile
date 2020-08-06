# Makefile for the Cansam library, which provides tools for SAM/BAM files.
#
#    Copyright (C) 2010-2012 Genome Research Ltd.
#
#    Author: John Marshall <jm18@sanger.ac.uk>
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
#    this list of conditions and the following disclaimer in the documentation
#    and/or other materials provided with the distribution.
# 3. Neither the names Genome Research Ltd and Wellcome Trust Sanger Institute
#    nor the names of its contributors may be used to endorse or promote
#    products derived from this software without specific prior written
#    permission.
#
# THIS SOFTWARE IS PROVIDED BY GENOME RESEARCH LTD AND ITS CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL GENOME RESEARCH LTD OR ITS CONTRIBUTORS
# BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
# THE POSSIBILITY OF SUCH DAMAGE.

srcdir = .
# SRC will be $(srcdir)/ when building in a separate directory
SRC =

CXX      = g++
CXXFLAGS = -Wall -Wextra -g -O2 -I$(srcdir)
LDFLAGS  =
LDLIBS   = -lz
AR       = ar
RANLIB   = ranlib

TOOLS   = samcat samcount samgroupbyname samhead samsort samsplit
OUTPUTS = libcansam.a $(TOOLS) test/runtests
all: $(OUTPUTS)

lib: libcansam.a

LIBOBJS = lib/alignment.o lib/collection.o lib/header.o lib/sambamio.o \
	  lib/samstream.o lib/ostream.o lib/rawfilebuf.o \
	  lib/interval.o lib/intervalmap.o \
	  lib/exception.o lib/system.o lib/utilities.o lib/version.o

libcansam.a: $(LIBOBJS)
	-rm -f $@
	$(AR) rc $@ $(LIBOBJS)
	-$(RANLIB) $@


sam_alignment_h = cansam/sam/alignment.h cansam/types.h cansam/sam/header.h
sam_header_h    = cansam/sam/header.h cansam/types.h
sam_interval_h  = cansam/interval.h cansam/types.h
sam_intervalmap_h=cansam/intervalmap.h cansam/interval.h cansam/types.h
lib_sambamio_h  = lib/sambamio.h cansam/sam/stream.h
lib_utilities_h = lib/utilities.h cansam/types.h

lib/alignment.o: lib/alignment.cpp $(sam_alignment_h) cansam/exception.h \
		 $(sam_header_h) $(lib_utilities_h) lib/wire.h
lib/collection.o: lib/collection.cpp $(sam_header_h) cansam/exception.h
lib/exception.o: lib/exception.cpp cansam/exception.h
lib/header.o: lib/header.cpp $(sam_header_h) cansam/exception.h $(lib_utilities_h)
lib/interval.o: lib/interval.cpp $(sam_interval_h) cansam/exception.h \
		$(sam_alignment_h)
lib/intervalmap.o: lib/intervalmap.cpp $(sam_intervalmap_h)
lib/ostream.o: lib/ostream.cpp $(sam_alignment_h) $(sam_header_h) \
	       $(lib_utilities_h)
lib/rawfilebuf.o: lib/rawfilebuf.cpp cansam/streambuf.h cansam/exception.h
lib/sambamio.o: lib/sambamio.cpp $(lib_sambamio_h) $(sam_alignment_h) \
		cansam/exception.h cansam/sam/stream.h $(lib_utilities_h) lib/wire.h
lib/samstream.o: lib/samstream.cpp cansam/sam/stream.h $(sam_alignment_h) \
		 cansam/exception.h cansam/streambuf.h $(lib_sambamio_h)
lib/system.o: lib/system.cpp
lib/utilities.o: lib/utilities.cpp lib/utilities.h
lib/version.o: lib/version.cpp cansam/version.h


MISC_OBJS = tools/samcat.o tools/samcount.o tools/samgroupbyname.o \
	    tools/samhead.o tools/samsort.o tools/samsplit.o \
	    tools/utilities.o examples/simplecat.o

samcat: tools/samcat.o tools/utilities.o libcansam.a
	$(CXX) $(LDFLAGS) -o $@ tools/samcat.o tools/utilities.o libcansam.a $(LDLIBS)

samcount: tools/samcount.o tools/utilities.o libcansam.a
	$(CXX) $(LDFLAGS) -o $@ tools/samcount.o tools/utilities.o libcansam.a $(LDLIBS)

samgroupbyname: tools/samgroupbyname.o tools/utilities.o libcansam.a
	$(CXX) $(LDFLAGS) -o $@ tools/samgroupbyname.o tools/utilities.o libcansam.a $(LDLIBS)

samhead: tools/samhead.o tools/utilities.o libcansam.a
	$(CXX) $(LDFLAGS) -o $@ tools/samhead.o tools/utilities.o libcansam.a $(LDLIBS)

samsort: tools/samsort.o libcansam.a
	$(CXX) $(LDFLAGS) -o $@ tools/samsort.o libcansam.a $(LDLIBS)

samsplit: tools/samsplit.o tools/utilities.o libcansam.a
	$(CXX) $(LDFLAGS) -o $@ tools/samsplit.o tools/utilities.o libcansam.a $(LDLIBS)

simplecat: examples/simplecat.o libcansam.a
	$(CXX) $(LDFLAGS) -o $@ examples/simplecat.o libcansam.a $(LDLIBS)


tools/samcat.o: tools/samcat.cpp $(sam_alignment_h) $(sam_header_h) \
		cansam/sam/stream.h tools/utilities.h
tools/samcount.o: tools/samcount.cpp $(sam_alignment_h) $(sam_header_h) \
		  cansam/sam/stream.h tools/utilities.h
tools/samgroupbyname.o: tools/samgroupbyname.cpp cansam/sam/algorithm.h \
			$(sam_alignment_h) cansam/exception.h $(sam_header_h) \
			cansam/sam/stream.h tools/utilities.h
tools/samhead.o: tools/samhead.cpp cansam/exception.h $(sam_header_h) \
		 cansam/sam/stream.h tools/utilities.h
tools/samsort.o: tools/samsort.cpp tools/samsort.h $(sam_alignment_h)
tools/samsplit.o: tools/samsplit.cpp $(sam_alignment_h) cansam/exception.h \
		  $(sam_header_h) cansam/sam/stream.h $(lib_utilities_h) \
		  tools/utilities.h
tools/utilities.o: tools/utilities.cpp tools/utilities.h cansam/version.h
examples/simplecat.o: examples/simplecat.cpp cansam/sam/header.h cansam/sam/alignment.h


test: test/runtests
	test/runtests test $(SRC)test

TEST_OBJS = test/runtests.o test/alignment.o test/header.o test/sam.o \
	    test/wire.o test/interval.o

test/runtests: $(TEST_OBJS) libcansam.a
	$(CXX) $(LDFLAGS) -o $@ $(TEST_OBJS) libcansam.a $(LDLIBS)

test/runtests.o: test/runtests.cpp test/test.h cansam/exception.h
test/alignment.o: test/alignment.cpp test/test.h $(sam_alignment_h)
test/header.o: test/header.cpp test/test.h $(sam_header_h)
test/interval.o: test/interval.cpp test/test.h $(sam_intervalmap_h)
test/sam.o: test/sam.cpp test/test.h $(sam_alignment_h) cansam/sam/stream.h
test/wire.o: test/wire.cpp test/test.h lib/wire.h

.PHONY: all clean doc docclean install lib tags test testclean uninstall

prefix      = /usr
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin
includedir  = $(prefix)/include
libdir      = $(exec_prefix)/lib
mandir      = $(prefix)/share/man
man1dir     = $(mandir)/man1
man3dir     = $(mandir)/man3

MKDIR_P = mkdir -p
INSTALL = install -p
INSTALL_PROGRAM = $(INSTALL)
INSTALL_DATA    = $(INSTALL) -m 644

install: libcansam.a $(TOOLS)
	$(MKDIR_P) $(DESTDIR)$(includedir)/cansam/sam $(DESTDIR)$(libdir)
	$(INSTALL_DATA) $(SRC)cansam/*.h $(DESTDIR)$(includedir)/cansam
	$(INSTALL_DATA) $(SRC)cansam/sam/*.h $(DESTDIR)$(includedir)/cansam/sam
	$(INSTALL_DATA) libcansam.a $(DESTDIR)$(libdir)/libcansam.a
	$(MKDIR_P) $(DESTDIR)$(bindir)
	$(INSTALL_PROGRAM) $(TOOLS) $(DESTDIR)$(bindir)
	$(MKDIR_P) $(DESTDIR)$(man1dir) $(DESTDIR)$(man3dir)
	$(INSTALL_DATA) $(SRC)tools/*.1 $(DESTDIR)$(man1dir)
	$(INSTALL_DATA) $(SRC)tools/cansam.3 $(DESTDIR)$(man3dir)/cansam.3

uninstall:
	-rm -rf $(DESTDIR)$(includedir)/cansam
	-rm -f $(DESTDIR)$(libdir)/libcansam.a
	-rm -f $(DESTDIR)$(bindir)/samcat $(DESTDIR)$(man1dir)/samcat.1
	-rm -f $(DESTDIR)$(bindir)/samcount
	-rm -f $(DESTDIR)$(bindir)/samgroupbyname $(DESTDIR)$(man1dir)/samgroupbyname.1
	-rm -f $(DESTDIR)$(bindir)/samhead $(DESTDIR)$(man1dir)/samhead.1
	-rm -f $(DESTDIR)$(bindir)/samsort $(DESTDIR)$(man1dir)/samsort.1
	-rm -f $(DESTDIR)$(bindir)/samsplit $(DESTDIR)$(man1dir)/samsplit.1
	-rm -f $(DESTDIR)$(man3dir)/cansam.3

doc:
	doxygen doc/Doxyfile
	@echo "Remove <p></p> from enum (see namespacesam.html)" >&2

tags:
	ctags -f TAGS [clt]*/*.h [cl]*/*/*.h [clt]*/*.cpp

clean:
	-rm -f $(OUTPUTS) $(LIBOBJS) $(MISC_OBJS) $(TEST_OBJS) TAGS

docclean:
	-rm -rf doc/html doc/latex

testclean:
	-rm -f test/*-out.[bs]am
