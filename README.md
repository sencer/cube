# Gaussian Cube File C library

This is a small, untested C library w/o any error-checking features to
manipulate Gaussian Cube Files -which is a chemical data format. Use at your
own risk.

While code is not documented, yet, `cube.h` file may give you some idea of what
is there, and what is not.

There are also several examples using the library:

beautify.c : which tries to wrap back badly separated portions of, say,
electron densities of a periodic image

interpolate.c : interpolates several images between two cube files, obtained
for example from a MD simulation. Files needs to be in cube folder, and
named as %d.cube.

qedat.c : reads Quantum Espresso cp.x dat files and writes cube files. Needs dat
files named as dat/%d.dat with %d = nfi; cp.pos and an input file named cube.in

    nat  ntyp
    num_of_atoms_type_i Z_i
    Lattice parameters matrix, 3x3
