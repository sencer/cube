# CC             = clang-3.6
CC             = gcc
CFLAGS         = -std=c99
EXEC           = beautify qedat #interpolate spherical average beautify2 testc test1
DEPS           = cube.o 3d.o files.o

.PHONY: all debug clean
all:   CFLAGS += -O3 -Wno-unused-result -fopenmp
debug: CFLAGS += -g -Wall -pedantic -fno-omit-frame-pointer -fsanitize=address

all: $(EXEC)
debug: $(EXEC)

$(EXEC): % : %.c $(DEPS)
	$(CC) $(CFLAGS) $@.c $(DEPS) -lm -o $@

%.o:%.c cube.h
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f *.o $(EXEC)  cube.py* cube_wrap.* _cube.so

%.c: ;
%.h: ;
Makefile: ;
