#ifndef Cube_H_XTY81TBG
#define Cube_H_XTY81TBG 1
#include <stdio.h>

#define L_LENGTH 110

typedef struct atom {
  int Z;                           // atomic number
  double coor[3];                  // coordinates
} Atom;

typedef struct cube_data{
  int nat,                         // number of atoms
      ngrid[3];                    // number of grid points
  double vsize[3],                 // voxel widths
         origin[3];                // location of origin (wrt coordinate system
  double *data;                    // of the atoms)
  Atom *atoms;                     // atoms as in Atom.
} Cube;

typedef struct layer_data{
  int dim[2],         // dimensions
      size;           // total size of Layer = dim[0] * dim[1]
  double **data;      // A 2D array
} Layer;


Cube *CubeInit(int nat, int* ngrid);         // allocate space for a cube file in memory as Cube
Cube *CubeInitFrom(Cube *cube);              // allocate space same size with another Cube
Cube *CubeCopy(Cube *cube);                  // completely copy some Cube into a new Cube
Cube *CubeRead(char *filename);              // read a cube file into memory as "Cube"
void CubeDelete(Cube *cube);                 // free the space used for Cube data in the memory
void CubeWrite(Cube *cube, char *filename);  // write Cube into a file


void CubeCopyData(Cube *dest, Cube *source);     // copy volumetric data from a Cube to another
void CubeCopyAtoms(Cube *dest, Cube *source);    // copy atoms from a Cube to another
void CubeSetVoxels(Cube *dest, double *vsize);   // set voxel sizes of a Cube, only orthorhombic implemented
void CubeSetOrigin(Cube *dest, double *origin);  // set origin of Cube data

int CubeDataSize(Cube *cube);                    // return the number of data points of Cube
double CubeVVolume(Cube *cube);                  // return volume of a single Voxes in Bohr^3
double CubeVolume(Cube *cube);                   // return total volume in Bohr^3

Layer *LayerInit(int d1, int d2);                 // allocate a Layer
Layer *CubeGetLayer(Cube *cube, int dir, int n);  // get nth Layer of a Cube, normal to 0:x 1:y 2:z direction
int *LayerIndices(Cube *cube, int dir, int n);    // get indices of the nth Layer normal to x/y/z direction
void CubePutLayer(Cube *cube, Layer *layer, int dir, int n); // put a layer into specified location in Cube
double CubeRotateLayers(Cube *cube, int dir, int n); // rotate all the layers in direction dir n times
                                                     // returns the amount origin is shifted to keep everything
                                                     // in place wrt Atoms

void LayerDelete(Layer *layer);                      // free the memory allocated for a Layer

void CubeMoveAtoms(Cube *cube, int dir, double r);  // move all atoms in direction dir r-much.

#endif /* end of include guard: Cube_H_XTY81TBG */

// vim: foldmethod=syntax
