CXX=g++
OPTFLAGS=-O3 -fopenmp
LIBPATH=-Ilib -ldl

all:ntcard 

SRCS=ntcard.cpp lib/Uncompress.cpp lib/SignalHandler.cpp lib/Fcontrol.cpp 

ntcard:$(SRCS)
	$(CXX) $(OPTFLAGS) $(LIBPATH) -o $@ $^

clean:
	rm ntcard

