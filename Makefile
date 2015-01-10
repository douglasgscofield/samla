# Makefile uses recursive $(MAKE) to build separate versions.
#
#   all (default):   debug version unoptimised (O_FLAG=-O0)
#   opt:             debug version optimised (O_FLAG=-O2)
#   profile:         default with gprof profiling (P_FLAG=-pg)
#   release:         no debug optimised (D_FLAG= O_FLAG=-O2)
#   release-profile: no debug optimised (D_FLAG= O_FLAG=-O2 P_FLAG=-pg)

CXX=		g++
O_FLAG=     -O0
D_FLAG=     -D_WITH_DEBUG -ggdb -g3 -fvar-tracking-assignments -fno-inline -fno-inline-small-functions -fno-eliminate-unused-debug-types
P_FLAG=         
# My fork of vcflib is now included as a git submodule
VCFLIB=     vcflib
CXXFLAGS=	-I$(VCFLIB) -Wall -D_FILE_OFFSET_BITS=64 $(O_FLAG) $(D_FLAG) $(P_FLAG)

RELEASE_O_FLAG = -O2
RELEASE_D_FLAG = -D_WITH_DEBUG
RELEASE_P_FLAG =

PROG=		samla

LIBS=		-lm -lz

VCFLIB_LIBS=-L$(VCFLIB)/tabixpp/ -ltabix

OBJS=		samla.o

HEAD_COMM=  version.h \
			SimpleOpt.h \
			$(VCFLIB)/src/Variant.h \
			$(VCFLIB)/src/split.h \
			$(VCFLIB)/src/join.h

HEAD=		$(HEAD_COMM)

VCFLIB_AR=	vcflib.a

VCFLIB_OBJS=$(VCFLIB)/src/Variant.o \
			$(VCFLIB)/src/split.o
VCFLIB_OTHER_OBJS=$(VCFLIB)/smithwaterman/SmithWatermanGotoh.o \
			$(VCFLIB)/smithwaterman/Repeats.o \
			$(VCFLIB)/smithwaterman/LeftAlign.o \
			$(VCFLIB)/smithwaterman/IndelAllele.o \
			$(VCFLIB)/src/ssw.o \
			$(VCFLIB)/src/ssw_cpp.o \
			$(VCFLIB)/fastahack/Fasta.o \
			$(VCFLIB)/fsom/fsom.o \
			$(VCFLIB)/tabixpp/tabix.o \
			$(VCFLIB)/tabixpp/bgzf.o
VCFLIB_DISORDER=$(VCFLIB)/smithwaterman/disorder.c
# NOTE: disorder.c is compiled fresh each time which seems wasteful, but it is C
# source, and is built in its own program as pure C, whereas it is build with
# other programs including vcflib as if it were C++.  If we always built as pure
# C then we would have to situationally have extern "C" around it.  It is smoother
# to simply build it new each time it is needed.  This is an issue inherited from
# vcflib.

#---------------------------  Main program


all: $(PROG)

samla: $(OBJS) $(VCFLIB_AR)
	$(CXX) $^ -o $@ $(VCFLIB_DISORDER) $(LIBS) $(VCFLIB_LIBS) $(CXXFLAGS) $(PROG_O_FLAG)

version.h: .FORCE
	./git-getversion.sh > version.h
	echo "#define CXX_VERSION \""`$(CXX) --version | head -n 1`"\"" >> version.h
	echo "#define CXXFLAGS \"$(CXXFLAGS)\"" >> version.h

.FORCE:

opt: samla-opt
	$(MAKE) O_FLAG=-O2
	mv samla $^

samla-opt:

profile: samla-profile
	$(MAKE) P_FLAG=-pg
	mv samla $^

samla-profile:

release: samla-release
	$(MAKE) O_FLAG=$(RELEASE_O_FLAG) D_FLAG=$(RELEASE_D_FLAG) P_FLAG=$(RELEASE_P_FLAG)
	mv samla $^

samla-release:

release-profile: samla-release-profile
	$(MAKE) O_FLAG=-O2 D_FLAG="-D_WITH_DEBUG" P_FLAG=-pg
	mv samla $^

samla-release-profile:


#---------------------------  Individual object files

$(VCFLIB_AR): $(VCFLIB_OBJS) $(VCFLIB_OTHER_OBJS)
	$(AR) -cru $@ $(VCFLIB_OBJS) $(VCFLIB_OTHER_OBJS)
	ranlib $@

# SimpeOpt.h is from http://code.jellycan.com/simpleopt and processes command-line args

# rebuild everything if the common headers change
$(OBJS): $(HEAD_COMM)

# rebuild the main file if any header changes
samla.o: $(HEAD)

$(VCFLIB_OBJS):
	cd $(VCFLIB) && make


#---------------------------  Other targets


clean:
	rm -f gmon.out *.o $(PROG) $(VCFLIB_AR)

clean-all: clean


#---------------------------  Obsolete and/or waiting for cleanup/reuse


