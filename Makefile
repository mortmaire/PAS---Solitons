CC=g++
CFLAGS=-O2 -Wall
LFLAGS=-lpng -lfftw3

main: main.cpp
# 	g++ main.cpp -o m $(LFLAGS) $(CFLAGS)-O2
	g++ main.cpp -o m $(LFLAGS) $(CFLAGS)