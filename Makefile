CC=g++
CFLAGS=-O2 -Wall -static
LFLAGS=-lpng -lfftw3 -lz

main: main.cpp
# 	g++ main.cpp -o m $(LFLAGS) $(CFLAGS)-O2
	g++ main.cpp -o m $(LFLAGS) $(CFLAGS)