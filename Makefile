CC=gcc
CFLAGS= -lm -Wall -pedantic -std=c99
EXEC = beautify qedat interpolate
DEPS = cube.o 3d.o files.o

.PHONY: all debug clean

all: CFLAGS += -fPIC -O3 -fopenmp
all: $(EXEC)

debug: CFLAGS += -g
debug: $(EXEC)

$(EXEC): % : %.c libcube.a
	$(CC) $@.c -L. -lcube $(CFLAGS) -o $@

libcube.a:$(DEPS)
	ar rvs libcube.a $(DEPS)
	rm -f *.o

%.o:%.c cube.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)

%.c: ;
%.h: ;
Makefile: ;
