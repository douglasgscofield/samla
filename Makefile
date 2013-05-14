
CXX=		/usr/bin/g++
#CXX=		/opt/local/bin/gcc
#CXXFLAGS=	-Wall -O3 -D_WITH_DEBUG 
CXXFLAGS=	-Wall -g -D_FILE_OFFSET_BITS=64 -D_WITH_DEBUG -fno-inline -I../vcflib
#CXXFLAGS=	-Wall -pg -g -D_WITH_DEBUG 

PROG=		samla

LIBS=		-lz

OBJS=		samla.o

HEAD_COMM=  SimpleOpt.h

HEAD=		$(HEAD_COMM)


#---------------------------  Main program


all: $(PROG)

smorgas: $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS) $(LIBS)


#---------------------------  Individual object files


# SimpeOpt.h is from http://code.jellycan.com/simpleopt and processes command-line args

# rebuild everything if the common headers change
$(OBJS): $(HEAD_COMM)

# rebuild the main file if any header changes
samla.o: $(HEAD)


#---------------------------  Other targets


clean:
	rm -f gmon.out *.o $(PROG)

clean-all: clean


#---------------------------  Obsolete and/or waiting for cleanup/reuse


