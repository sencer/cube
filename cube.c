#include "cube.h"

// TODO CubeGetRegion etc should return only data

Cube *CubeInit(int nat, int ngrid[3])
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
  memset(cube->invvsize, 0, sizeof(cube->vsize));
  memset(cube->origin, 0, sizeof(cube->origin));

  cube->comment[0][0] = '\0';
  cube->comment[1][0] = '\0';

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
  Inv3D(vsize, dest->invvsize);
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

int *CubeRegionIndices(Cube *cube, int p[3], int r[3])
{
  /* @brief a helper method. returns indices of data points in Cube *cube in
   * the parallel-piped region from p to q
   * @p Cube *cube
   * @p int p[3] 3D indices of one corner of the region
   * @p int r[3] 3D indices of one corner of the region
   * @return int* a 1D array of integers containing 1D indices
   */
  int c[] = { cube->ngrid[0], cube->ngrid[1], cube->ngrid[2] },
      /* TODO I probably want the absolute values of r[0] - p[0]*/
      d[] = { r[0] - p[0] + 1,  r[1] - p[1] + 1 ,  r[2] - p[2] + 1  },
      *indices = malloc( d[0] * d[1] * d[2] * sizeof(int)),
      x, y, z, xp, yp, zp;

  for( x = 0; x < d[0] ; x++ ){
    xp = x + p[0];
    xp = xp - cube->ngrid[0] * floor((double)xp / cube->ngrid[0]);
    for( y = 0; y < d[1] ; y++ ){
      yp = y + p[1];
      yp = yp - cube->ngrid[1] * floor((double)yp / cube->ngrid[1]);
      for( z = 0; z < d[2] ; z++ ){
        zp = z + p[2];
        zp = zp - cube->ngrid[2] * floor((double)zp / cube->ngrid[2]);
        indices[d[2]*d[1]*x+d[2]*y+z] = c[2] * c[1] * xp + c[2] * yp + zp;
      }
    }
  }

  return indices;
}

Cube *CubeSphericalRegion(Cube *c, int pos[3], double rad)
{
  int r[3], p[3],
      d[3] = {
        ceil(rad/c->vsize[0][0]),
        ceil(rad/c->vsize[1][1]),
        ceil(rad/c->vsize[2][2])
      };

  double v[3];

  for (int i = 0; i < 3; ++i)
  {
    p[i] = pos[i] - d[i] + 1;
    r[i] = pos[i] + d[i] - 1;
  }

  Cube *nc = CubeGetRegion(c, p, r);

  for (int i = 0; i < d[0] * 2 - 1; ++i)
  {
    for (int j = 0; j < d[1] * 2 - 1; ++j)
    {
      for (int k = 0; k < d[2] * 2 - 1; ++k)
      {
        for (int l = 0; l < 3; ++l)
        {
          v[l] = (d[0] - i - 1) * c->vsize[0][l]
               + (d[1] - j - 1) * c->vsize[1][l]
               + (d[2] - k - 1) * c->vsize[2][l];
        }
        if(VecLen(v) > rad)
        {
          nc->data[nc->ngrid[2]*nc->ngrid[1]*i +nc->ngrid[2]*j + k] = 0;
        }
      }
    }
  }
  return nc;
}

Cube *CubeGetRegion(Cube *cube, int p[3], int r[3])
{
  /* @brief returns the Cube containing only the data contained in
   * between p[3] and r[3], and all the atoms in Cube *cube
   * @p Cube *cube, the cube to get a region from
   * @p int p[3] 3D indices of one corner of the region
   * @p int r[3] 3D indices of one corner of the region
   * @return Cube with a subset of the data.
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

void CubePutRegion(Cube *dest, Cube *source, int p[3])
{
  /* @brief overwrites the data in *dest "above" p[3], with the data *source
   * @p Cube *dest
   * @p Cube *source
   * @p int p[3] 3D indices of one corner of the region
   * @return void
   */
  int dim = CubeDataSize(source),
      r[3] = {
        p[0]+source->ngrid[0]-1,
        p[1]+source->ngrid[1]-1,
        p[2]+source->ngrid[2]-1
      },
      *indices = CubeRegionIndices(dest, p, r), i;
  for( i = 0; i < dim; i++ )
  {
    dest->data[indices[i]] = source->data[i];
  }
  free(indices);
}

int *CubeLayerIndices(Cube *cube, int dir, int n)
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
      *d_index = NULL,
      *c_index = NULL;
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
    c_index = CubeLayerIndices(cube, dir, i);
    d_index = CubeLayerIndices(cube, dir, k);
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
  CubeSetVoxels(cube, vsize);


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
    i += sscanf(line, "%lf %lf %lf %lf %lf %lf", &cube->data[i+0],
        &cube->data[i+1],
        &cube->data[i+2],
        &cube->data[i+3],
        &cube->data[i+4],
        &cube->data[i+5]);
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
    fprintf(f, "%13.5E", cube->data[i]);
    if (! (((i+1) % cube->ngrid[2] ) % 6))
    {
      fprintf(f, "\n");
    }
  }

  fclose(f);
}

double CubeLayerMax(Cube *c, int dir, int layer)
{
  int *indices = CubeLayerIndices(c, dir, layer),
      dim = CubeDataSize(c) / c->ngrid[dir];
  double max = -100;
  for(int i = 0; i < dim; ++i)
  {
    if(fabs(c->data[indices[i]]) > max)
      max = fabs(c->data[indices[i]]);
  }
  free(indices);
  return max;
}

void CubeBeautify(Cube *c, double thr)
{
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < c->ngrid[i] / 2; ++j)
    {
      if(CubeLayerMax(c, i, j) < thr)
      {
        CubeRotateLayers(c, i, -j);
        break;
      }
      else if(CubeLayerMax(c, i, c->ngrid[i] - 1 - j) < thr)
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
  int p[6] = { 0, 0, 0, tmp->ngrid[0] - 1, tmp->ngrid[1] - 1, tmp->ngrid[2] - 1 };
  for(int i = 0; i < 3; ++i)
  {
    for(int j = 0; j < tmp->ngrid[i]-1; ++j)
    {
      if(CubeLayerMax(tmp, i, j) > thr) break;
      p[i] += 1;
    }
    for(int j = tmp->ngrid[i] - 1; j > 0; --j)
    {
      if(CubeLayerMax(tmp, i, j) > thr) break;
      p[i + 3] -= 1;
    }
  }
  *c = CubeGetRegion(tmp, p, p + 3);
  CubeDelete(tmp);
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

void CubeR2BoxI(Cube *c, double r[3], int box[3])
{

  double coor[3];

  CubeR2Box(c, r, coor);

  for (int i = 0; i < 3; ++i)
  {
    box[i] = round(coor[i]);
  }
}

void CubeR2Box(Cube *c, double r[3], double box[3])
{

  double coor[3];

  for (int i = 0; i < 3; ++i)
  {
    coor[i] = r[i] - c->origin[i];
  }

  for(int i = 0; i < 3; ++i)
  {
    box[i] = VecDot(coor, c->invvsize[i]);
  }
}

void CubeBox2R(Cube *c, double box[3], double r[3])
{
  for (int i = 0; i < 3; ++i)
  {
    VecShift(r, c->vsize[i], box[i]);
  }
  VecShift(r, c->origin, 1);
}

void CubeWrapAtoms(Cube *c)
{
  for (int i = 0; i < c->nat; ++i)
  {
    for (int j = 0; j < 3; ++j)
    {
      VecShift(c->atoms[i].coor, c->vsize[j], -c->ngrid[j] *
          floor(VecDot(c->invvsize[j], c->atoms[i].coor)/c->ngrid[j]));
    }
  }
}

// vim: foldmethod=syntax
