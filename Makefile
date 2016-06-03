CXX=g++
OPTFLAGS=-O3
LIBPATH=-Ilib

EXEC:ntcard 

SRCS=ntcard.cpp 

ntcard:$(SRCS)
	$(CXX) $(OPTFLAGS) $(LIBPATH) -o $@ $^

clean:
	rm $(EXEC)

