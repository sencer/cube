#include "cube.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


Cube *CubeInit(int nat, int* ngrid)
{
  Cube *cube = malloc(sizeof(Cube));
  cube->atoms = malloc(nat * sizeof(Atom));
  cube->data = malloc(ngrid[0] * ngrid[1] * ngrid[2] * sizeof(double));

  cube->nat = nat;
  memcpy(cube->ngrid, ngrid, 3 * sizeof(int));

  return cube;
}

Cube *CubeInitFrom(Cube *cube)
{
  Cube *newCube = CubeInit(cube->nat, cube->ngrid);
  CubeSetVoxels(newCube, cube->vsize);
  memcpy(newCube->origin, cube->origin, 3 * sizeof(double));
  return newCube;
}

Cube *CubeCopy(Cube *cube)
{
  Cube *newCube = CubeInitFrom(cube);
  CubeCopyAtoms(newCube, cube);
  CubeCopyData(newCube, cube);
  return newCube;
}

void CubeDelete(Cube *cube)
{
  free(cube->atoms);
  free(cube->data);
  free(cube);
}

void CubeCopyAtoms(Cube *dest, Cube *source)
{
  memcpy(dest->atoms, source->atoms, source->nat * sizeof(Atom));
}

void CubeCopyData(Cube *dest, Cube *source)
{
  memcpy(dest->data, source->data, CubeDataSize(source) * sizeof(double));
}

void CubeSetVoxels(Cube *dest, double *vsize)
{
  memcpy(dest->vsize, vsize, 3 * sizeof(double));
}

void CubeSetOrigin(Cube *dest, double *origin)
{
  memcpy(dest->origin, origin, 3 * sizeof(double));
}

int CubeDataSize(Cube *cube)
{
  return cube->ngrid[0] * cube->ngrid[1] * cube->ngrid[2];
}

double CubeVVolume(Cube *cube)
{
  return cube->vsize[0] * cube->vsize[1] * cube->vsize[2];
}

double CubeVolume(Cube *cube)
{
  return CubeDataSize(cube) * CubeVVolume(cube);
}

int *CubeRegionIndices(Cube *cube, int *p, int *r)
{
  int c[] = { cube->ngrid[0], cube->ngrid[1], cube->ngrid[2] },
      d[] = { r[0] - p[0] + 1,  r[1] - p[1] + 1 ,  r[2] - p[2] + 1  },
     *indices = malloc( d[0] * d[1] * d[2] * sizeof(int)),
     x, y, z;

  for( x = 0; x < d[0] ; x++ ){
    for( y = 0; y < d[1] ; y++ ){
      for( z = 0; z < d[2] ; z++ ){
        indices[d[2]*d[1]*x+d[2]*y+z] = c[2] * c[1] * (x + p[0]) + c[2] * (y + p[1]) + z + p[2];
      }
    }
  }

  return indices;
}

Cube *CubeGetRegion(Cube *cube, int *p, int *r)
{
  int *indices = CubeRegionIndices(cube, p, r),
      d[3] = { r[0] - p[0] + 1, r[1] - p[1] + 1, r[2] - p[2] + 1 },
      dim = d[0] * d[1] * d[2],
      i;
  Cube *c = CubeInit(0, d);
  for(i = 0; i < 3; i++) {
    c->origin[i] = cube->origin[i] + p[i] * cube->vsize[i];
    c->vsize[i] = cube->vsize[i];
  }
  for(i=0; i<dim; i++)
  {
    c->data[i] = cube->data[indices[i]];
  }
  free(indices);
  return c;
}

void CubePutRegion(Cube *dest, Cube *source, int *p)
{
  int dim = CubeDataSize(source),
      r[3] = { p[0]+source->ngrid[0]-1,
               p[1]+source->ngrid[1]-1,
               p[2]+source->ngrid[2]-1 },
      *indices = CubeRegionIndices(dest, p, r), i;
  for( i = 0; i < dim; i++ )
  {
    dest->data[indices[i]] = source->data[i];
  }
  free(indices);
}


int *LayerIndices(Cube *cube, int dir, int n)
{
  int p[3] = {0, 0, 0},
      r[3] = {cube->ngrid[0]-1, cube->ngrid[1]-1, cube->ngrid[2]-1};
  p[dir] = n;
  r[dir] = n;
  return CubeRegionIndices(cube, p, r);
}

Cube *CubeGetLayer(Cube *cube, int dir, int n)
{
  int p[3] = {0, 0, 0},
      r[3] = {cube->ngrid[0]-1, cube->ngrid[1]-1, cube->ngrid[2]-1};
  p[dir] = n;
  r[dir] = n;
  return CubeGetRegion(cube, p, r);
}

void CubePutLayer(Cube *dest, Cube *source, int dir, int n)
{
  int p[3] = {0, 0, 0};
  p[dir] = n;
  CubePutRegion(dest, source, p);
}

double CubeRotateLayers(Cube *cube, int dir, int n)
{
  double shift, *data = malloc(CubeDataSize(cube) * sizeof(double));
  int i, j, k,
      dim = CubeDataSize(cube)/cube->ngrid[dir],
      *d_index,
      *c_index;
  for (i = 0; i < cube->ngrid[dir]; ++i)
  {
    k = i + n;
    if(k < 0)
    {
      k += cube->ngrid[dir];
    }
    else if(k >= cube->ngrid[dir])
    {
      k -= cube->ngrid[dir];
    }
    c_index = LayerIndices(cube, dir, i);
    d_index = LayerIndices(cube, dir, k);
    for (j = 0; j < dim; ++j)
    {
      data[d_index[j]] = cube->data[c_index[j]];
    }
  }
  free(c_index);
  free(d_index);
  free(cube->data);
  cube->data = data;
  shift = n * cube->vsize[dir];
  cube->origin[dir] -= shift;

  return shift;
}

void CubeMoveAtoms(Cube *cube, int dir, double r)
{
  int i;
  for (i = 0; i < cube->nat; ++i)
  {
    cube->atoms[i].coor[dir] += r;
  }
}

Cube *CubeRead(char* filename)
{
  FILE *f  = fopen(filename, "r");
  char line[L_LENGTH];
  int i, nat, ngrid[3];
  double origin[3], vsize[3];

  //Header data: # of grid points, atoms, species, as well as v_size
  fgets(line, L_LENGTH, f); //
  fgets(line, L_LENGTH, f); // Dismiss Title Lines

  fgets(line, L_LENGTH, f);
  sscanf(line, "%d %lf %lf %lf", &nat, origin, origin + 1, origin + 2);

  // assume orthorhombic for now
  fgets(line, L_LENGTH, f);
  sscanf(line,"%d %lf %*f %*f", ngrid, vsize);
  fgets(line, L_LENGTH, f);
  sscanf(line,"%d %*f %lf %*f", ngrid + 1, vsize + 1);
  fgets(line, L_LENGTH, f);
  sscanf(line,"%d %*f %*f %lf", ngrid + 2, vsize + 2);

  Cube *cube   = CubeInit(nat, ngrid);
  memcpy(cube->origin, origin, sizeof(origin));
  memcpy(cube->vsize, vsize, sizeof(vsize));


  //read nuclear positions of all atoms
  for (i = 0; i < cube->nat; ++i){
    fgets(line, L_LENGTH, f);
    sscanf(line, "%d %*f %lf %lf %lf", &(cube->atoms[i].Z)
                                     , &(cube->atoms[i].coor[0])
                                     , &(cube->atoms[i].coor[1])
                                     , &(cube->atoms[i].coor[2]));
  }
  //read grid data
  i = 0;
  while(fgets(line, L_LENGTH, f)){
    sscanf(line, "%lf %lf %lf %lf %lf %lf", &cube->data[i+0],
                                            &cube->data[i+1],
                                            &cube->data[i+2],
                                            &cube->data[i+3],
                                            &cube->data[i+4],
                                            &cube->data[i+5]);
    i += 6;
  }

  fclose(f);
  return cube;
}

void CubeWrite(Cube *cube, char *filename)
{
  FILE *f = fopen(filename, "w");

  fprintf(f, "\n\n%-5d %12.6f %12.6f %12.6f\n", cube->nat,
                                                cube->origin[0],
                                                cube->origin[1],
                                                cube->origin[2]);

  fprintf(f, "%-5d %12.6f %12.6f %12.6f\n", cube->ngrid[0],
                                            cube->vsize[0], 0., 0.);
  fprintf(f, "%-5d %12.6f %12.6f %12.6f\n", cube->ngrid[1], 0.,
                                            cube->vsize[1], 0.);
  fprintf(f, "%-5d %12.6f %12.6f %12.6f\n", cube->ngrid[2], 0., 0.,
                                            cube->vsize[2]);

  for (int i = 0; i < cube->nat; ++i)
  {
    fprintf(f, "%-3d %8.4f %14.10f %14.10f %14.10f\n", cube->atoms[i].Z,
                                                   (float) cube->atoms[i].Z,
                                                   cube->atoms[i].coor[0],
                                                   cube->atoms[i].coor[1],
                                                   cube->atoms[i].coor[2]);
  }
  int size = cube->ngrid[0] * cube->ngrid[1] * cube->ngrid[2];
  for (int i = 0; i < size; ++i)
  {
    fprintf(f, "%14.6E ", cube->data[i]);
    if (i % 6 == 5)
    {
      fprintf(f, "\n");
    }
  }

  fclose(f);
}

// vim: foldmethod=syntax
