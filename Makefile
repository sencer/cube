# makefile for APC524 HW1
CC=gcc
CFLAGS= -lm -Wall -pedantic -std=c99
EXEC = interpolate beautify test trim binary interp2
DEPS = cube.o

.PHONY: all debug

all: CFLAGS += -O3 -fopenmp
all: $(EXEC)

debug: CFLAGS += -g
debug: $(EXEC)

$(EXEC): $(DEPS)
	$(CC) $@.c $(DEPS) $(CFLAGS) -o $@

%.o:%.c cube.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)
