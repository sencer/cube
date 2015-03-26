#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include "cube.h"

#define THR 6.E-3
#define MIN(a, b) ((a<b)?a:b)
#define MAX(a, b) ((a>b)?a:b)
#define DEG 57.2957795786
#define B2A 0.529177249

int is_dat(char const *name)
{
  size_t len = strlen(name);
  return len > 4 && strcmp(name + len - 4, ".dat") == 0;
}

int FilesList(char files[5000][11])
{
  DIR *directory = opendir("dat");
  struct dirent *ent;
  int nfile = 0;

  while ((ent = readdir(directory)) != NULL)
  {
    if(is_dat(ent->d_name))
    {
      strcpy(files[nfile], ent->d_name);
      nfile++;
    }
  }
  return nfile;
}

void ReadInput(int *nat, int *ntyp, int *num, int *z, double celldm[3][3])
{
  FILE *f = fopen("cube.in", "r");
  char *line;
  size_t len = 0;
  getline(&line, &len, f);
  sscanf(line, "%d %d", nat, ntyp);
  for (int i = 0; i < *ntyp; ++i)
  {
    getline(&line, &len, f);
    sscanf(line, "%d %d", num+i, z+i);
    if(i>0) num[i]+=num[i-1];
  }
  getline(&line, &len, f);
  sscanf(line, "%lf %lf %lf", celldm[0], celldm[1], celldm[2]);
  getline(&line, &len, f);
  sscanf(line, "%lf %lf %lf", celldm[0]+1, celldm[1]+1, celldm[2]+1);
  getline(&line, &len, f);
  sscanf(line, "%lf %lf %lf", celldm[0]+2, celldm[1]+2, celldm[2]+2);
}

void GetAtoms(int index, Cube *c, int nat, int *z, int *num)
{
  // read atoms
  FILE *f = fopen("cp.pos", "r");
  char nfi[12];
  char *line;
  size_t len = 0;
  sprintf(nfi, " %d ", index);
  do
  {
    getline(&line, &len, f);
  }
  while(strstr(line, nfi) == NULL);
  int type = 0;
  for (int i = 0; i < nat; ++i)
  {
    getline(&line, &len, f);
    sscanf(line, "%lf %lf %lf", &c->atoms[i].coor[0], &c->atoms[i].coor[1],
        &c->atoms[i].coor[2]);
    c->atoms[i].Z = z[type];
    if(i+1 == num[type]) type++;
  }
}

FILE *ReadDatHeader(char *fname, int *ngrid)
{
  size_t len = 0;
  char *line;
  FILE *f = fopen(fname, "rb");

  getline(&line, &len, f);
  while(strstr(line, "INFO nr1") == NULL)
  {
    getline(&line, &len, f);
  }
  sscanf(line, "    <INFO nr1=\"%d\" nr2=\"%d\" nr3=\"%d\"/>",
      &ngrid[0], &ngrid[1], &ngrid[2]);
  getline(&line, &len, f);
  return f;
}

void ReadDatData(FILE *f, int ngrid[3], int dim, Cube *c)
{
  size_t len = 0;
  char *line;
  double *data = (double *) malloc(dim * sizeof(double));
  for (int i = 0; i < ngrid[2]; ++i)
  {
    while(strstr(line, " type=\"real\" size=\"") == NULL)
    {
      getline(&line, &len, f);
    }
    sprintf(line, "xxx");
    fread(data, 12, 1, f);
    fread(data, sizeof(double), dim, f);
    for(int j = 0; j < ngrid[1]; j++)
    {
      for(int k = 0; k < ngrid[0]; k++)
      {
        c->data[ngrid[1]*ngrid[2]*k + ngrid[2] * j + i] = data[j * ngrid[0] + k];
      }
    }
  }
  free(data);
  fclose(f);
}

double pAverageNN(Cube *c, int coor[3])
{
  int p[3] = { MAX(coor[0] - 1, 0), MAX(coor[1] - 1, 0), MAX(coor[2] - 1, 0) },
      r[3] = { MIN(coor[0] + 1, c->ngrid[0]-1), MIN(coor[1] + 1, c->ngrid[1]-1), MIN(coor[2] + 1, c->ngrid[2]-1) },
      dim, *ind = CubeRegionIndices(c, p, r);

  double av = 0;
  dim = (r[0]-p[0]+1)*(r[1]-p[1]+1)*(r[2]-p[2]+1);

  for(int i = 0; i < dim; ++i)
  {
    av += fabs(c->data[ind[i]]);
  }

  free(ind);

  return av / dim;
}

void Inverse(double celldm[3][3], double inv[3][3])
{
  double det = 0;

  for(int i=0;i<3;++i)
      det = det + (celldm[0][i]*(celldm[1][(i+1)%3]*celldm[2][(i+2)%3] - celldm[1][(i+2)%3]*celldm[2][(i+1)%3]));

   for(int i=1;i<4;++i)
   {
      for(int j=1;j<4;++j)
      {
        inv[i-1][j-1] = (celldm[i%3][j%3]*celldm[(i+1)%3][(j+1)%3] - celldm[i%3][(j+1)%3]*celldm[(i+1)%3][j%3])/ det;
      }
   }
}

double VecMult(double a[3], double b[3])
{
  return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

double VecLen(double a[3])
{
  return sqrt(VecMult(a, a));
}

void VecShift(double pnt[3], double vec[3], double c)
{
  for(int i = 0; i < 3; ++i)
  {
    pnt[i] += vec[i] * c;
  }
}


int Coor2DataPnt(int atom, Cube *c, double celldm[3][3], double inv[3][3])
{
  int coor[3];

  for(int i = 0; i < 3; ++i)
  {
    // TODO Fix for the cases where origin of the cube is not 0

    VecShift(c->atoms[atom].coor, celldm[i],
        -1 * floor(VecMult(c->atoms[atom].coor, inv[i])));

    coor[i] = VecMult(c->atoms[atom].coor, inv[i]) * celldm[i][i]/c->vsize[i][i];
  }

  return pAverageNN(c, coor) > THR ? 1: 0;
}

int main()
{
  char files[5000][11];
  double celldm[3][3], inv[3][3], *norm;
  int nat = 0,
      ntyp = 0,
      num[30],
      z[30],
      ngrid[3],
      dim,
      dummy[3],
      nfile = FilesList(files);
  char fname[20];

  if(nfile < 1) return -1;

  ReadInput(&nat, &ntyp, num, z, celldm);
  Inverse(celldm, inv);
  norm = (double [3]) {
    VecLen(celldm[0]),
    VecLen(celldm[1]),
    VecLen(celldm[2])
  };
  sprintf(fname, "dat/%s", files[0]);
  fclose(ReadDatHeader(fname, ngrid));
  dim = ngrid[0] * ngrid[1];

#pragma omp parallel for private(fname) shared(dummy)
  for(int file = 0; file < nfile; ++file)
  {
    FILE *f;
    char com[255] = "";
    int pos = 0;
    Cube *c = CubeInit(nat, ngrid);

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        c->vsize[i][j] = celldm[i][j] / ngrid[i];
      }
    }

    sprintf(fname, "dat/%s", files[file]);
    f = ReadDatHeader(fname, dummy);
    ReadDatData(f, ngrid, dim, c);
    GetAtoms(atoi(files[file]), c, nat, z, num);
    for (int i = 0; i < c->nat; i++)
    {
      if(c->atoms[i].Z != 22) break;
      if(Coor2DataPnt(i, c, celldm, inv))
      {
        pos += sprintf(com + pos, "%d ", i);
      }
    }
    CubeRotateLayers(c, 0, -37);
    CubeRotateLayers(c, 1, 6);
    CubeRotateLayers(c, 2, -112);
    CubeBeautify(c, THR);
    CubeTrim(&c, THR);
    sprintf(fname, "%d.cube", atoi(files[file]));

    sprintf(c->comment[0], "%-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f",
        norm[0]*B2A, norm[1]*B2A, norm[2]*B2A,
        DEG*acos(VecMult(celldm[0], celldm[2])/(norm[2]*norm[0])),
        DEG*acos(VecMult(celldm[1], celldm[2])/(norm[1]*norm[2])),
        DEG*acos(VecMult(celldm[0], celldm[1])/(norm[0]*norm[1])));

    strcpy(c->comment[1], com);
    CubeWrite(c, fname);
    CubeDelete(c);
  }
  return 0;
}

// vi: foldmethod=syntax
