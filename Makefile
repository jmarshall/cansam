CXXFLAGS = -Wall -Wextra -O2 -I. -ansi -pedantic
LDFLAGS  = -L.

OUTPUTS = libcansam.a samcat samsort simplecat test/runtests
all: $(OUTPUTS)

lib: libcansam.a

LIBOBJS = lib/alignment.o lib/collection.o lib/header.o lib/bam.o lib/sam.o \
	  lib/samstream.o lib/istream.o lib/ostream.o lib/rawfilebuf.o \
	  lib/utilities.o lib/zio.o

libcansam.a: $(LIBOBJS)
	$(AR) cr $@ $(LIBOBJS)
	ranlib $@

lib/alignment.o: lib/alignment.cpp sam/alignment.h
lib/bam.o: lib/bam.cpp sam/stream.h lib/wire.h
lib/collection.o: lib/collection.cpp sam/collection.h
lib/header.o: lib/header.cpp sam/header.h
lib/istream.o: lib/istream.cpp sam/header.h sam/alignment.h
lib/ostream.o: lib/ostream.cpp sam/header.h sam/alignment.h
lib/rawfilebuf.o: lib/rawfilebuf.cpp sam/iobuffer.h
lib/sam.o: lib/sam.cpp
lib/samstream.o: lib/samstream.cpp sam/stream.h
lib/utilities.o: lib/utilities.cpp lib/utilities.h
lib/zio.o: lib/zio.cpp lib/zio.h


MISC_OBJS = utilities/samcat.o utilities/samsort.o examples/simplecat.o

samcat: utilities/samcat.o libcansam.a
	$(CXX) $(LDFLAGS) -o $@ utilities/samcat.o -lcansam

samsort: utilities/samsort.o libcansam.a
	$(CXX) $(LDFLAGS) -o $@ utilities/samsort.o -lcansam

simplecat: examples/simplecat.o libcansam.a
	$(CXX) $(LDFLAGS) -o $@ examples/simplecat.o -lcansam

utilities/samcat.o: utilities/samcat.cpp sam/header.h sam/alignment.h
utilities/samsort.o: utilities/samsort.cpp utilities/samsort.h
examples/simplecat.o: examples/simplecat.cpp sam/header.h sam/alignment.h


test: test/runtests
	test/runtests

TEST_OBJS = test/runtests.o test/alignment.o test/header.o test/wire.o

test/runtests: $(TEST_OBJS) libcansam.a
	$(CXX) $(LDFLAGS) -o $@ $(TEST_OBJS) -lcansam

test/runtests.o: test/runtests.cpp test/test.h
test/alignment.o: test/alignment.cpp test/test.h sam/alignment.h
test/header.o: test/header.cpp test/test.h sam/header.h
test/wire.o: test/wire.cpp test/test.h lib/wire.h

.PHONY: all clean doc docclean lib tags test

doc:
	doxygen doc/Doxyfile
	@echo "Remove <p></p> from enum (see namespacesam.html)" >&2

tags:
	@#ctags -o TAGS */*.h */*.cpp
	ctags -o TAGS {examples,lib,sam,test,utilities}/*.[ch]*

clean:
	-rm -f $(OUTPUTS) $(LIBOBJS) $(MISC_OBJS) $(TEST_OBJS) TAGS

docclean:
	-rm -rf doc/html doc/latex
