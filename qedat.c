#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <dirent.h>
#include "cube.h"
#define THR 7.E-3

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

void ReadInput(int *nat, int *ntyp, int *num, int *z, double *celldm)
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
  sscanf(line, "%lf %lf %lf", celldm + 0, celldm + 3, celldm + 6);
  getline(&line, &len, f);
  sscanf(line, "%lf %lf %lf", celldm + 1, celldm + 4, celldm + 7);
  getline(&line, &len, f);
  sscanf(line, "%lf %lf %lf", celldm + 2, celldm + 5, celldm + 8);
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

int main()
{
  char files[5000][11];
  double celldm[9];
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
  sprintf(fname, "dat/%s", files[0]);
  fclose(ReadDatHeader(fname, ngrid));
  dim = ngrid[0] * ngrid[1];

#pragma omp parallel for private(fname) shared(dummy)
  for(int file = 0; file < nfile; ++file)
  {
    FILE *f;
    Cube *c = CubeInit(nat, ngrid);
    for (int i = 0; i < 3; ++i)
      for (int j = 0; j < 3; ++j)
        c->vsize[i][j] = celldm[i*3+j] / ngrid[i];
    sprintf(fname, "dat/%s", files[file]);
    f = ReadDatHeader(fname, dummy);
    ReadDatData(f, ngrid, dim, c);
    GetAtoms(atoi(files[file]), c, nat, z, num);
    CubeBeautify(c, THR);
    CubeTrim(&c, THR);
    sprintf(fname, "%d.cube", atoi(files[file]));
    CubeWrite(c, fname);
    CubeDelete(c);
  }
  return 0;
}
