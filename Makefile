VERSION = 1.0
SHELL = /bin/sh

CC = g++
CFLAGS = -std=c++17 -O3 -funroll-loops -fopenmp -ffast-math
#CFLAGS = -std=c++11 -g -Wall -Wextra -pedantic -funroll-loops -DNDEBUG
OBJS = gsa-is/gsacak.o sais-lite-2.4.1/sais.o utils.o sacamats.o main.o 

EXEC = clean matchstat

.cpp.o:
	$(CC) $(CFLAGS) -c $<

all: $(EXEC)

matchstat: $(OBJS)
	$(CC) $(CFLAGS) -o matchstat $(OBJS) $(LIBS)

clean:
	/bin/rm -f matchstat *.o
