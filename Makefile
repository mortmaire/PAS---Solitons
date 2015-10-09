CC=g++
CFLAGS=-O2 -Wall -static
LIBS=-lpng -lfftw3 -lz

# main: main.cpp
# 	g++ main.cpp -o m $(LFLAGS) $(CFLAGS)-O2
# 	g++ main.cpp -o m $(LFLAGS) $(CFLAGS)
	
# SRC=$(wildcard *.cpp)
SRC=main.cpp
# 
main: $(SRC)
	g++ -o $@ $^ $(CFLAGS) $(LIBS)

# SRCS = $(wildcard *.cpp)
# SRCS = main2.cpp

# PROGS = $(patsubst %.cpp,%,$(SRCS))
