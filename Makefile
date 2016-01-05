FLAGS = -std=c++11 -lm
CXX = g++
IDIR = -I/home/sdorsher  -I/home/sdorsher/tnt -I/home/sdorsher/jama

wave: main.o ReferenceElement.o
	$(CXX) $(FLAGS) $(IDIR) main.o ReferenceElement.o -o wave

main.o: main.cpp ReferenceElement.h
	$(CXX) $(FLAGS) $(IDIR) -c main.cpp

ReferenceElement.o: ReferenceElement.cpp
	$(CXX) $(FLAGS) $(IDIR) -c ReferenceElement.cpp TNT2.h tnt_array1D_extn.h tnt_array2D_extn.h


