#define _GNU_SOURCE

#include "cube.h"
#include "files.h"
#include "3d.h"

#define THR 6.E-3

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

int Coor2DataPnt(int atom, Cube *c, double celldm[3][3], double inv[3][3])
{
  int coor[3];

  for(int i = 0; i < 3; ++i)
  {
    // TODO Fix for the cases where origin of the cube is not 0


    VecShift(c->atoms[atom].coor, celldm[i],
        -1 * floor(VecDot(c->atoms[atom].coor, inv[i])));

    coor[i] = VecDot(c->atoms[atom].coor, inv[i]) * celldm[i][i]/c->vsize[i][i];
  }

  return pAverageNN(c, coor) > THR ? 1: 0;
}

int main()
{
  char files[5000][11];
  double celldm[3][3], inv[3][3], *norm, vsize[3][3];
  int nat = 0,
      ntyp = 0,
      num[30],
      z[30],
      ngrid[3],
      dim,
      dummy[3],
      nfile = FilesList(files, "dat");
  char fname[20];

  if(nfile < 1) return -1;

  ReadInput(&nat, &ntyp, num, z, celldm);
  Inv3D(celldm, inv);
  norm = (double [3]) {
    VecLen(celldm[0]),
    VecLen(celldm[1]),
    VecLen(celldm[2])
  };
  sprintf(fname, "dat/%s", files[0]);
  fclose(ReadDatHeader(fname, ngrid));
  dim = ngrid[0] * ngrid[1];

  memcpy(vsize, celldm, sizeof(celldm));
  for(int i = 0; i < 3; ++i)
  {
    VecScale(vsize[i], 1/ngrid[i]);
  }

#pragma omp parallel for private(fname) shared(dummy)
  for(int file = 0; file < nfile; ++file)
  {
    FILE *f;
    char com[255] = "";
    int pos = 0;
    Cube *c = CubeInit(nat, ngrid);
    CubeSetVoxels(c, vsize);

    sprintf(fname, "dat/%s", files[file]);
    f = ReadDatHeader(fname, dummy);
    ReadDatData(f, ngrid, dim, c);
    GetAtoms(atoi(files[file]), c, nat, z, num);
    for (int i = 0; i < c->nat; i++)
    {
      /* if(c->atoms[i].Z != 22) break; */
      if(Coor2DataPnt(i, c, celldm, inv))
      {
        pos += sprintf(com + pos, "%d ", i);
      }
    }
    CubeBeautify(c, THR);
    CubeTrim(&c, THR);
    sprintf(fname, "cube/%d.cube", atoi(files[file]));

    sprintf(c->comment[0], "%-10.6f %-10.6f %-10.6f %-10.6f %-10.6f %-10.6f",
        norm[0]*B2A, norm[1]*B2A, norm[2]*B2A,
        DEG*acos(VecDot(celldm[0], celldm[2])/(norm[2]*norm[0])),
        DEG*acos(VecDot(celldm[1], celldm[2])/(norm[1]*norm[2])),
        DEG*acos(VecDot(celldm[0], celldm[1])/(norm[0]*norm[1])));

    strcpy(c->comment[1], com);
    CubeWrite(c, fname);
    CubeDelete(c);
  }
  return 0;
}

// vi: foldmethod=syntax
