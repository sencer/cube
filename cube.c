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

Layer *LayerInit(int d0, int d1)
{
  int i;
  Layer *layer = malloc(sizeof(Layer));
  layer->size = d0 * d1;
  layer->dim[0] = d0;
  layer->dim[1] = d1;
  layer->data = malloc(d0 * sizeof(double*));
  for (i = 0; i < d0; ++i)
  {
    layer->data[i] = malloc(d1 * sizeof(double));
  }
  return layer;
}

Layer *CubeGetLayer(Cube *cube, int dir, int n)
{
/* Somewhat unnaturally
 * dir 0 gets x-planes as y-z
 * dir 1 gets y-planes as x-z
 * dir 2 gets z-planes as x-y
 * n for n-th plane
 */
  int i, j;
  Layer *layer;
  int d[] = { cube->ngrid[0], cube->ngrid[1], cube->ngrid[2] };
  switch (dir) {
    case 0:
      layer = LayerInit(d[1], d[2]);

      for (i = 0; i < d[1]; ++i)
      {
        memcpy(layer->data[i], &(cube->data[i*d[2]+n*d[2]*d[1]]), d[2]*sizeof(double));
      }

      break;

    case 1:

      layer = LayerInit(d[0], d[2]);

      for (i = 0; i < d[0]; ++i)
      {
        memcpy(layer->data[i], &(cube->data[n*d[2]+i*d[1]*d[2]]), d[2]*sizeof(double));
      }

      break;

    case 2:

      layer = LayerInit(d[0], d[1]);

      for (i = 0; i < d[0]; ++i)
      {
        for (j = 0; j < d[1]; ++j)
        {
          layer->data[i][j] = cube->data[n + d[2] * j + i * d[1] * d[2]];
        }
      }

      break;
  }
  return layer;
}

void CubePutLayer(Cube *cube, Layer *layer, int dir, int n)
{
  int i, j;
  int d[] = { cube->ngrid[0], cube->ngrid[1], cube->ngrid[2] };
  switch (dir) {
    case 0:
      for (i = 0; i < d[1]; ++i)
      {
        memcpy(&(cube->data[i*d[2]+n*d[2]*d[1]]), layer->data[i], d[2]*sizeof(double));
      }
      break;
    case 1:
      for (i = 0; i < d[0]; ++i)
      {
        memcpy(&(cube->data[n*d[2]+i*d[1]*d[2]]), layer->data[i], d[2]*sizeof(double));
      }
      break;
    case 2:
      for (i = 0; i < d[0]; ++i)
      {
        for (j = 0; j < d[1]; ++j)
        {
          cube->data[n + d[2] * j + i * d[1] * d[2]] = layer->data[i][j];
        }
      }
      break;
  }
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

void LayerDelete(Layer *layer)
{
  int i;
  for (i = 0; i < layer->dim[0]; ++i)
  {
    free(layer->data[i]);
  }
  free(layer->data);
  free(layer);
}

int *LayerIndices(Cube *cube, int dir, int n)
{
  int i, j;
  int d[] = { cube->ngrid[0], cube->ngrid[1], cube->ngrid[2] };
  int *indices = malloc(d[0] * d[1] * d[2] / d[dir] * sizeof(int));
  switch (dir) {
    case 0:
      for (i = 0; i < d[1]; ++i)
      {
        for (j = 0; j < d[2]; ++j)
        {
          indices[i * d[2] + j] = i * d[2] + n * d[2] * d[1] + j;
        }
      }
      break;
    case 1:
      for (i = 0; i < d[0]; ++i)
      {
        for (j = 0; j < d[2]; ++j)
        {
          indices[i * d[2] + j] = n * d[2] + i * d[1] * d[2] + j;
        }
      }
      break;
    case 2:
      for (i = 0; i < d[0]; ++i)
      {
        for (j = 0; j < d[1]; ++j)
        {
          indices[i * d[1] + j] = n + d[2] * j + i * d[1] * d[2];
        }
      }
      break;
  }
  return indices;
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
