CXX=g++
CXXFLAGS=-g -std=c++11 -Wall -pedantic -O1 -larmadillo
BIN=fem1

SRC=$(wildcard *.cpp)
OBJ=$(SRC:%.cpp=%.o)

all:$(OBJ)
	$(CXX) 	-o 	$(BIN)	$^

%.o:%.c
	$(CXX) 	$@ 	-c 	$<

clean:
	rm -f *.o *.out
	rm $(BIN)