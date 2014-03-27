
CXX=		g++
#CXX=		/opt/local/bin/gcc
#CXXFLAGS=	-Wall -O3 -D_WITH_DEBUG 
CXXFLAGS=	-I../vcflib -Wall -D_FILE_OFFSET_BITS=64 -D_WITH_DEBUG -ggdb -g3 -fvar-tracking-assignments -fno-inline -fno-inline-small-functions -O0 -fno-eliminate-unused-debug-types
#CXXFLAGS=	-Wall -pg -g -D_WITH_DEBUG 

PROG=		samla

LIBS=		-lm -lz

#VCFLIB_LIBS=-L../vcflib/ -L../vcflib/tabixpp/ -ltabix
VCFLIB_LIBS=-L../vcflib/tabixpp/ -ltabix

OBJS=		samla.o

HEAD_COMM=  SimpleOpt.h \
			../vcflib/src/Variant.h \
			../vcflib/src/split.h \
			../vcflib/src/join.h

HEAD=		$(HEAD_COMM)

VCFLIB_AR=	vcflib.a

VCFLIB_OBJS=../vcflib/src/Variant.o \
			../vcflib/src/split.o
#VCFLIB_OTHER_OBJS=../vcflib/smithwaterman/SmithWatermanGotoh.o \
#			../vcflib/smithwaterman/Repeats.o \
#			../vcflib/smithwaterman/disorder.o \
#			../vcflib/smithwaterman/LeftAlign.o \
#			../vcflib/smithwaterman/IndelAllele.o \
#			../vcflib/src/ssw.o \
#			../vcflib/src/ssw_cpp.o \
#			../vcflib/fastahack/Fasta.o \
#			../vcflib/fsom/fsom.o \
#			../vcflib/tabixpp/tabix.o \
#			../vcflib/tabixpp/bgzf.o
VCFLIB_OTHER_OBJS=../vcflib/smithwaterman/SmithWatermanGotoh.o \
			../vcflib/smithwaterman/Repeats.o \
			../vcflib/smithwaterman/LeftAlign.o \
			../vcflib/smithwaterman/IndelAllele.o \
			../vcflib/src/ssw.o \
			../vcflib/src/ssw_cpp.o \
			../vcflib/fastahack/Fasta.o \
			../vcflib/fsom/fsom.o \
			../vcflib/tabixpp/tabix.o \
			../vcflib/tabixpp/bgzf.o
VCFLIB_DISORDER=../vcflib/smithwaterman/disorder.c
# TODO: disorder.c is compiled fresh each time which is wasteful, should probably
# have extern "C" around it within vcflib... but that breaks its own compilation
# within other parts of vcflib.

#---------------------------  Main program


all: $(PROG)

samla: $(OBJS) $(VCFLIB_AR)
	$(CXX) $^ -o $@ $(VCFLIB_DISORDER) $(LIBS) $(VCFLIB_LIBS) $(CXXFLAGS)


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
	cd ../vcflib && make


#---------------------------  Other targets


clean:
	rm -f gmon.out *.o $(PROG) $(VCFLIB_AR)

clean-all: clean


#---------------------------  Obsolete and/or waiting for cleanup/reuse


