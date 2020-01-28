CXX=g++
# Linux
CXXFLAGS= -v -std=c++14 -O1 -Wall -pedantic #-march=native -ffast-math

PROG = strassenTest
HDRS = strassen.hpp
SRCS = main.cpp strassen.cpp
OBJS = $(SRCS: .cpp = .0)

.DEFAULT_GOAL := all

all: $(PROG)
 
$(PROG): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $^
$(OBJS): $(HDRS)
# My notes:
# $@ ... file of the target rule (should be strassenTest)
# $^ ... The names of all the prerequisites, with spaces between them. (should be als OBJS here)

# Lecture notes:
# $@ ... contains the filename of the current rule’s target.
# $^ ... lists all source files

.PHONY: clean
clean:
	rm $(PROG)