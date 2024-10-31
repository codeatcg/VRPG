
CXX = g++
CPPFLAG = -D PYMODULE
PYSO = minipg$(shell python3-config --extension-suffix)
INCLUDE =  $(shell python3 -m pybind11 --includes)
CXXFLAG = -O3 -Wall -std=c++11 -pthread
LIBS = -lz

OBJ1 = gfa2view.o minipg.o gz.o
OBJ2 = anno.o refgene.o ghAnno.o gz.o
OBJ2 = anno.o refgene.o refbed.o ghAnno.o gz.o 
OBJ3 = gaf2rbed.o
OBJ4 = nodeSeq.o

all : gfa2view GraphAnno gaf2rbed nodeSeq $(PYSO)

gfa2view : $(OBJ1)
	$(CXX) $(CXXFLAG) $(OBJ1) -o $@ $(LIBS)

GraphAnno : $(OBJ2)
	$(CXX) $(CXXFLAG) $(OBJ2) -o $@ $(LIBS)

gaf2rbed : $(OBJ3)
	$(CXX) $(CXXFLAG) $(OBJ3) -o $@

nodeSeq : $(OBJ4)
	$(CXX) $(CXXFLAG) $(OBJ4) -o $@
    
%.o : %.cpp
	$(CXX) $(CXXFLAG) -c $< -o $@ $(LIBS)

$(PYSO) : minipg.cpp
	$(CXX) $(CPPFLAG) $(CXXFLAG) -shared -fPIC $(INCLUDE) minipg.cpp  -o $(PYSO)

.PHONY : clean

clean:
	rm *.o

