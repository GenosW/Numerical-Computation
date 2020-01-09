# requirements on ubuntu/linux mint
# sudo apt-get install build-essential libopenblas-base libopenblas-dev libeigen3-dev

CXX=g++
# Single thread -> no omp
#CXXFLAGS=-std=c++11 -Wall -pedantic -O3 -ffast-math -march=native
# Linux
CXXFLAGS= -v -std=c++14 -O1 -Wall -pedantic -march=native -ffast-math #-fopenmp
# Mac
#CXXFLAGS= -v -std=c++14 -O3 -Wall -pedantic -march=native -ffast-math -Xpreprocessor -fopenmp -lomp

SRCS = main.cpp

.DEFAULT_GOAL := all

all: strassen
 
strassen: $(SRCS) Makefile
	$(CXX) $(SRCS) -o strassen $(CXXFLAGS)

.PHONY: clean
clean:
	rm strassen