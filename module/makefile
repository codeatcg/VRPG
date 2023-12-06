
CXX = g++
CPPFLAG = -D PYMODULE
PYSO = minipg$(shell python3-config --extension-suffix)
INCLUDE =  $(shell python3 -m pybind11 --includes)
CXXFLAG = -O3 -Wall -std=c++11 -pthread -lz

OBJ1 = gfa2view.o minipg.o gz.o
OBJ2 = anno.o refgene.o ghAnno.o gz.o

all : gfa2view GraphAnno $(PYSO)

gfa2view : $(OBJ1)
	$(CXX) $(CXXFLAG) $(OBJ1) -o $@

GraphAnno : $(OBJ2)
	$(CXX) $(CXXFLAG) $(OBJ2) -o $@

%.o : %.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@

$(PYSO) : minipg.cpp
	$(CXX) $(CPPFLAG) $(CXXFLAG) -shared -fPIC $(INCLUDE) minipg.cpp  -o $(PYSO)

.PHONY : clean

clean:
	rm *.o

