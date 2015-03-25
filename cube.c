#include "cube.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


Cube *CubeInit(int nat, int* ngrid)
{
  /* @brief Allocates the memory needed for a Cube
   * @p int nat number of atoms
   * @p int[3] ngrid, number of grid points in 0,1,2=x,y,z directions
   * @return Cube with no data & atoms. Neither they are set to 0.
   */
  Cube *cube = malloc(sizeof(Cube));
  cube->atoms = malloc(nat * sizeof(Atom));
  cube->data = malloc(ngrid[0] * ngrid[1] * ngrid[2] * sizeof(double));

  cube->nat = nat;
  memcpy(cube->ngrid, ngrid, 3 * sizeof(int));

  memset(cube->vsize, 0, sizeof(cube->vsize));
  memset(cube->origin, 0, sizeof(cube->origin));

  sprintf(cube->comment[0], "");
  sprintf(cube->comment[1], "");

  return cube;
}

Cube *CubeInitFrom(Cube *cube)
{
  /* @brief Initialize a Cube with same nat, ngrid, vsize, origin  with @p cube
   * @p Cube *cube the Cube to copy nat etc from
   * @return Cube new cube with empty atoms list and no data
   * data isn't set to 0 either!
   */
  Cube *newCube = CubeInit(cube->nat, cube->ngrid);
  CubeSetVoxels(newCube, cube->vsize);
  CubeSetOrigin(newCube, cube->origin);
  return newCube;
}

Cube *CubeCopy(Cube *cube)
{
  /* @brief Copy @p cube exactly
   * @p Cube *cube Cube to copy
   * @return Cube new Cube, same with @p cube
   */
  Cube *newCube = CubeInitFrom(cube);
  CubeCopyAtoms(newCube, cube);
  CubeCopyData(newCube, cube);
  return newCube;
}

void CubeDelete(Cube *cube)
{
  /* @brief free the memory allocated for @p cube
   * @p Cube *cube Cube to be deleted
   */
  free(cube->atoms);
  free(cube->data);
  free(cube);
}

void CubeCopyAtoms(Cube *dest, Cube *source)
{
  /* @brief Copy atoms from a Cube to another
   * @p Cube dest   The cube to be copied to
   * @p Cube source The cube to be copied from
   */
  memcpy(dest->atoms, source->atoms, source->nat * sizeof(Atom));
}

void CubeCopyData(Cube *dest, Cube *source)
{
  /* @brief Copy grid data from a Cube to another
   * @p Cube dest   The cube to be copied to
   * @p Cube source The cube to be copied from
   */
  memcpy(dest->data, source->data, CubeDataSize(source) * sizeof(double));
}

void CubeSetVoxels(Cube *dest, double vsize[3][3])
{
  /* @brief Set voxel size of a cube
   * @p Cube dest   The cube to apply new voxel sizes
   * @p double[3][3] vsize an array of length = 3, containing new voxel sizes
   */
  memcpy(dest->vsize, vsize, sizeof(dest->vsize));
}

void CubeSetOrigin(Cube *dest, double origin[3])
{
  /* @brief Set origin for the data in the Cube file wrt origin for the atoms
   * @p Cube dest   The cube to apply new coordinates for the
   * @p double[3] origin an array of length = 3, containing the new origin
   */
  memcpy(dest->origin, origin, 3 * sizeof(double));
}

int CubeDataSize(Cube *cube)
{
  /* @brief a helper method that returns the number of data points kept
   * in the Cube. To find the actual size in memory, multiply with
   * sizeof(double)
   * @p Cube *cube
   * @return int size
   */
  return cube->ngrid[0] * cube->ngrid[1] * cube->ngrid[2];
}

// TODO Won't work if one of the lattice parameters does not contain two 
// zero elements. Although that is unusual.
double CubeVVolume(Cube *cube)
{
  /* @brief a helper method that returns the physical volume of a voxel
   * in the Cube.
   * @p Cube *cube
   * @return double volume
   */
  return cube->vsize[0][0] * cube->vsize[1][1] * cube->vsize[2][2];
}

double CubeVolume(Cube *cube)
{
  /* @brief a helper method that returns the physical volume of the
   * all data in the Cube.
   * @p Cube *cube
   * @return double volume
   */
  return CubeDataSize(cube) * CubeVVolume(cube);
}

int *CubeRegionIndices(Cube *cube, int *p, int *r)
{
  /* @brief a helper method. returns indices of data points in Cube *cube in
   * the parallel-piped region from p to q
   * @p Cube *cube
   * @p int *p, p[3] 3D indices of one corner of the region
   * @p int *r, r[3] 3D indices of one corner of the region
   * @return int* a 1D array of integers containing 1D indices
   */
  int c[] = { cube->ngrid[0], cube->ngrid[1], cube->ngrid[2] },
      /* TODO I probably want the absolute values of r[0] - p[0]*/
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
  /* @brief returns the Cube containing only the data contained in
   * between p[3] and r[3], and all the atoms in Cube *cube
   * @p Cube *cube, the cube to get a region from
   * @p int *p, p[3] 3D indices of one corner of the region
   * @p int *r, r[3] 3D indices of one corner of the region
   * @return Cube, with same atoms as *cube and a subset of the data.
   */
  int *indices = CubeRegionIndices(cube, p, r),
      d[3] = { r[0] - p[0] + 1, r[1] - p[1] + 1, r[2] - p[2] + 1 },
      dim = d[0] * d[1] * d[2],
      i, j;
  Cube *c = CubeInit(cube->nat, d);
  CubeCopyAtoms(c, cube);
  CubeSetVoxels(c, cube->vsize);
  CubeSetOrigin(c, cube->origin);
  for(i = 0; i < 3; i++)
  {
    for(j = 0; j < 3; j++)
    {
      c->origin[j] += p[i] * cube->vsize[i][j];
    }
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
  /* @brief overwrites the data in *dest "above" p[3], with the data *source
   * @p Cube *dest
   * @p Cube *source
   * @p int *p, p[3] 3D indices of one corner of the region
   * @return void
   */
  // TODO need to check if there is enough space in dest-above-p
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
  /* @brief get 1D indices of the nth layer from the cube in direction dir
   * @p Cube *cube source cube
   * @p int dir direction 0 is x, 1 is y, 2 is z
   * @p int n nth layer
   * @return int* a 1D array of 1D integers
   */
  // TODO needs to check if n is out of boundaries
  int p[3] = {0, 0, 0},
      r[3] = {cube->ngrid[0]-1, cube->ngrid[1]-1, cube->ngrid[2]-1};
  p[dir] = n;
  r[dir] = n;
  return CubeRegionIndices(cube, p, r);
}

Cube *CubeGetLayer(Cube *cube, int dir, int n)
{
  /* @brief get nth layer from the cube in direction dir
   * @p Cube *cube source cube
   * @p int dir direction 0 is x, 1 is y, 2 is z
   * @p int n nth layer
   * @return int* a 1D array of 1D integers
   */
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

void CubeRotateLayers(Cube *cube, int dir, int n)
{
  double *data = malloc(CubeDataSize(cube) * sizeof(double));
  int i, j, k,
      dim = CubeDataSize(cube) / cube->ngrid[dir],
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
  for (j=0; j < 3; j++)
  {
    cube->origin[j] -= n * cube->vsize[dir][j];
  }
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
  double origin[3], vsize[3][3];

  //Header data: # of grid points, atoms, species, as well as v_size
  fgets(line, L_LENGTH, f); //
  fgets(line, L_LENGTH, f); // Dismiss Title Lines

  fgets(line, L_LENGTH, f);
  sscanf(line, "%d %lf %lf %lf", &nat, origin, origin + 1, origin + 2);

  fgets(line, L_LENGTH, f);
  sscanf(line,"%d %lf %lf %lf", ngrid, vsize[0], vsize[0] + 1, vsize[0] + 2);
  fgets(line, L_LENGTH, f);
  sscanf(line,"%d %lf %lf %lf", ngrid + 1, vsize[1], vsize[1] + 1, vsize[1] + 2);
  fgets(line, L_LENGTH, f);
  sscanf(line,"%d %lf %lf %lf", ngrid + 2, vsize[2], vsize[2] + 1, vsize[2] + 2);

  Cube *cube   = CubeInit(nat, ngrid);
  memcpy(cube->origin, origin, sizeof(origin));
  memcpy(cube->vsize, vsize, sizeof(vsize));


  //read nuclear positions of all atoms
  for (i = 0; i < cube->nat; ++i){
    fgets(line, L_LENGTH, f);
    sscanf(line, "%d %*f %lf %lf %lf", &(cube->atoms[i].Z),
                                      &(cube->atoms[i].coor[0]),
                                      &(cube->atoms[i].coor[1]),
                                      &(cube->atoms[i].coor[2]));
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

  fprintf(f, "%s\n", cube->comment[0]);
  fprintf(f, "%s\n", cube->comment[1]);

  fprintf(f, "%-5d %12.6f %12.6f %12.6f\n", cube->nat,
                                                cube->origin[0],
                                                cube->origin[1],
                                                cube->origin[2]);

  fprintf(f, "%-5d %12.6f %12.6f %12.6f\n", cube->ngrid[0],
                                            cube->vsize[0][0],
                                            cube->vsize[0][1],
                                            cube->vsize[0][2]);
  fprintf(f, "%-5d %12.6f %12.6f %12.6f\n", cube->ngrid[1],
                                            cube->vsize[1][0],
                                            cube->vsize[1][1],
                                            cube->vsize[1][2]);
  fprintf(f, "%-5d %12.6f %12.6f %12.6f\n", cube->ngrid[2],
                                            cube->vsize[2][0],
                                            cube->vsize[2][1],
                                            cube->vsize[2][2]);
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

double GetLayerMax(Cube *c, int dir, int layer, int dim)
{
  int *indices = LayerIndices(c, dir, layer);
  double max = -100;
  for(int i = 0; i < dim; ++i)
  {
    if(c->data[indices[i]] > max)
      max = c->data[indices[i]];
  }
  free(indices);
  return max;
}

void CubeBeautify(Cube *c, double thr)
{
  int dim;
  for(int i = 0; i < 3; ++i)
  {
    dim = CubeDataSize(c) / c->ngrid[i];
    for(int j = 0; j < c->ngrid[i] / 2; ++j)
    {
      if(GetLayerMax(c, i, j, dim) < thr)
      {
        CubeRotateLayers(c, i, -j);
        break;
      }
      else if(GetLayerMax(c, i, c->ngrid[i] - 1 - j, dim) < thr)
      {
        CubeRotateLayers(c, i, j);
        break;
      }
    }
  }
}

void CubeTrim(Cube **c, double thr)
{
  Cube *tmp = *c;
  int dim, p[6] = { 0, 0, 0, tmp->ngrid[0] - 1, tmp->ngrid[1] - 1, tmp->ngrid[2] - 1 };
  for(int i = 0; i < 3; ++i)
  {
    dim = CubeDataSize(tmp) / tmp->ngrid[i];
    for(int j = 0; j < tmp->ngrid[i] / 2; ++j)
    {
      if(GetLayerMax(tmp, i, j, dim) > thr) break;
      p[i] += 1;
    }
    for(int j = tmp->ngrid[i] - 1; j > tmp->ngrid[i] / 2 - 1; --j)
    {
      if(GetLayerMax(tmp, i, j, dim) > thr) break;
      p[i + 3] -= 1;
    }
  }
  *c = CubeGetRegion(tmp, p, p + 3);
  CubeDelete(tmp);
}

double WeightedAverage(double a, double b, double weight)
{
  return (1 - weight) * a + weight * b;
}

Cube **CubeInterpolate(Cube *first, Cube *last, int n)
{
  double w;
  int dim;
  Cube **cubes = malloc(n * sizeof(Cube));
  for(int i = 0; i < n; ++i)
  {
    w = 1.0 * (i + 1) / (n + 1);
    cubes[i] = CubeInitFrom(first);
    for (int j = 0; j < first->nat; ++j)
    {
      cubes[i]->atoms[j].Z = first->atoms[j].Z;
      for (int k = 0; k < 3; ++k)
      {
        cubes[i]->atoms[j].coor[k] = WeightedAverage(first->atoms[j].coor[k], last->atoms[j].coor[k], w);
      }
    }
    dim = CubeDataSize(first);
    for(int j = 0; j < dim; ++j)
    {
      cubes[i]->data[j] = WeightedAverage(first->data[j], last->data[j], w);
    }
  }
  return cubes;
}

// vim: foldmethod=syntax
