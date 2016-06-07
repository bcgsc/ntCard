CXX=g++
OPTFLAGS=-O3
LIBPATH=-Ilib -ldl

all:ntcard 

SRCS=ntcard.cpp lib/Uncompress.cpp  lib/Fcontrol.cpp lib/Sequence.cpp lib/Options.cpp lib/FastaReader.cpp 

ntcard:$(SRCS)
	$(CXX) $(OPTFLAGS) $(LIBPATH) -o $@ $^

clean:
	rm ntcard

