# makefile for APC524 HW1
CC=gcc
CFLAGS= -lm -O3 -Wall -pedantic -std=c99
# CFLAGS= -lm -g -Wall -pedantic -std=c99
EXEC = interpolate beautify
DEPS = cube.o

.PHONY: all

all: interpolate beautify

beautify: $(DEPS)
interpolate: $(DEPS)

%.o:%.c cube.h layer.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)
