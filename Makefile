VERSION = 1.0
SHELL = /bin/sh

CC = g++
CFLAGS = -std=c++17 -O3 -funroll-loops -fopenmp -ffast-math
#CFLAGS = -std=c++11 -g -Wall -Wextra -pedantic -funroll-loops -DNDEBUG
OBJS = gsa-is/gsacak.o sais-lite-2.4.1/sais.o utils.o sacamats.o main.o 

EXEC = clean sacamats

.cpp.o:
	$(CC) $(CFLAGS) -c $<

all: $(EXEC)

sacamats: $(OBJS)
	$(CC) $(CFLAGS) -o sacamats $(OBJS) $(LIBS)

clean:
	/bin/rm -f sacamats *.o gsa-is/*.o sais-lite-2.4.1/*.o
