# makefile for APC524 HW1
CC=gcc
CFLAGS= -lm -O3 -Wall -pedantic -std=c99 -fopenmp
# CFLAGS= -lm -g -Wall -pedantic -std=c99
EXEC = interpolate beautify test trim binary interp2
DEPS = cube.o

.PHONY: all

all: $(EXEC)

$(EXEC): $(DEPS)
	$(CC) $@.c $(DEPS) $(CFLAGS) -o $@

%.o:%.c cube.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)
