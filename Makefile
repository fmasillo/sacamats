VERSION = 1.0
SHELL = /bin/sh

CC = g++
CFLAGS = -std=c++20 -O3 -funroll-loops -fopenmp -march=native -DNDEBUG -DLIBSAIS_OPENMP
#CFLAGS = -std=c++11 -g -Wall -Wextra -pedantic -funroll-loops -DNDEBUG
OBJS = utils.o induce_parallel.o libsais/src/libsais.o sacamats.o main.o #malloc_count/malloc_count.o
LIBS = -ltbb

EXEC = clean sacamats

.cpp.o:
	$(CC) $(CFLAGS) -c $<

all: $(EXEC)

sacamats: $(OBJS)
	$(CC) $(CFLAGS) -o sacamats $(OBJS) $(LIBS)

clean:
	/bin/rm -f *.o libsais/src/*.o
	
#sacamats *.o
