# Gaussian Cube File C library

This is a small, untested C library w/o any error-checking features to
manipulate Gaussian Cube Files -which is a chemical data format. Use at your
own risk.

While code is not documented, yet, `cube.h` file may give you some idea of what
is there, and what is not.

Also `beautify.c`, which tries to wrap back badly separated portions of, say,
electron densities of a periodic image; and `interpolate.c`, which interpolates
several new images between two cube files can be used as examples.
