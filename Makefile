CC=gcc
CFLAGS= -lm -Wall -pedantic -std=c99
EXEC = beautify qedat interpolate
DEPS = cube.o

.PHONY: all debug

all: CFLAGS += -O3 -fopenmp
all: $(EXEC)

debug: CFLAGS += -g
debug: $(EXEC)

$(EXEC): % : $(DEPS) %.c
	$(CC) $@.c $(DEPS) $(CFLAGS) -o $@

%.o:%.c cube.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)

%.c: ;
%.h: ;
Makefile: ;
