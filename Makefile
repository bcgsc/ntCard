CXX=g++ -Wall
OPTFLAGS=-O3 -fopenmp
LIBPATH=-Ilib -ldl

all:ntcard

COMMON_SRC=lib/Uncompress.cpp lib/SignalHandler.cpp lib/Fcontrol.cpp

SRCS=ntcard.cpp $(COMMON_SRC)

ntcard:$(SRCS)
	$(CXX) $(OPTFLAGS) $(LIBPATH) -o $@ $^

nthll:nthll.cpp $(COMMON_SRC)
	$(CXX) $(OPTFLAGS) $(LIBPATH) -o $@ $^ 

clean:
	rm ntcard

