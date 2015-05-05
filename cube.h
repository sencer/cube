#ifndef Cube_H_XTY81TBG
#define Cube_H_XTY81TBG 1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "3d.h"

#define L_LENGTH 110

typedef struct atom {
  int Z;                           // atomic number
  double coor[3];                  // coordinates
} Atom;

typedef struct cube_data{
  int nat,                         // number of atoms
      ngrid[3];                    // number of grid points
  double vsize[3][3],              // voxel widths
         invvsize[3][3],
         origin[3];                // location of origin (wrt coordinate system
  double *data;                    // of the atoms)
  char comment[2][255];
  Atom *atoms;                     // atoms as in Atom.
} Cube;

Cube *CubeInit(int nat, int ngrid[3]);         // allocate space for a cube file in memory as Cube
Cube *CubeInitFrom(Cube *cube);              // allocate space same size with another Cube
Cube *CubeCopy(Cube *cube);                  // completely copy some Cube into a new Cube
Cube *CubeRead(char *filename);              // read a cube file into memory as "Cube"
void CubeDelete(Cube *cube);                 // free the space used for Cube data in the memory
void CubeWrite(Cube *cube, char *filename);  // write Cube into a file


void CubeCopyData(Cube *dest, Cube *source);          // copy volumetric data from a Cube to another
void CubeCopyAtoms(Cube *dest, Cube *source);         // copy atoms from a Cube to another
void CubeSetVoxels(Cube *dest, double vsize[3][3]);   // set voxel sizes of a Cube
void CubeSetOrigin(Cube *dest, double origin[3]);       // set origin of Cube data

int CubeDataSize(Cube *cube);                    // return the number of data points of Cube
double CubeVVolume(Cube *cube);                  // return volume of a single Voxes in Bohr^3
double CubeVolume(Cube *cube);                   // return total volume in Bohr^3

int *CubeRegionIndices(Cube *cube, int p[3], int r[3]);
Cube *CubeGetRegion(Cube *cube, int p[3], int r[3]);
void CubePutRegion(Cube *dest, Cube *source, int p[3]);

int *CubeLayerIndices(Cube *cube, int dir, int n);    // get indices of the nth Layer normal to x/y/z direction
Cube *CubeGetLayer(Cube *cube, int dir, int n);
void CubePutLayer(Cube *dest, Cube *source, int dir, int n);
void CubeRotateLayers(Cube *cube, int dir, int n); // rotate all the layers in direction dir n times
                                                     // returns the amount origin is shifted to keep everything
                                                     // in place wrt Atoms

void CubeMoveAtoms(Cube *cube, int dir, double r);  // move all atoms in direction dir r-much.

void CubeBeautify(Cube *c, double thr);
void CubeTrim(Cube **c, double thr);
Cube **CubeInterpolate(Cube *first, Cube *last, int n);
double CubeLayerMax(Cube *c, int dir, int layer);

void CubeWrapAtoms(Cube *c);
#endif /* end of include guard: Cube_H_XTY81TBG */

// vim: foldmethod=syntax
