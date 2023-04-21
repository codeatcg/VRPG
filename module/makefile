
CXX = g++
CPPFLAG = -D PYMODULE
PYSO = minipg$(shell python3-config --extension-suffix)
INCLUDE =  $(shell python3 -m pybind11 --includes)
CXXFLAG = -O3 -Wall -std=c++11 -lpthread

all : gfa2view $(PYSO)

gfa2view : gfa2view.o minipg.o
	$(CXX) $(CXXFLAG) -o gfa2view gfa2view.o minipg.o
%.o : %.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@

$(PYSO) : minipg.cpp
	$(CXX) $(CPPFLAG) $(CXXFLAG) -shared -fPIC $(INCLUDE) minipg.cpp  -o $(PYSO)

.PHONY : clean

clean:
	rm *.o

