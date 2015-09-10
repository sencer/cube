#define _GNU_SOURCE
#define THR 6.E-3

#include "cube.h"
#include "files.h"
#include "3d.h"
#include <unistd.h>


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
  char *line = malloc(100 * sizeof(char));
  FILE *f = fopen(fname, "rb");

  getline(&line, &len, f);
  while(strstr(line, "INFO nr1") == NULL)
  {
    getline(&line, &len, f);
  }
  sscanf(line, "    <INFO nr1=\"%d\" nr2=\"%d\" nr3=\"%d\"/>",
      &ngrid[0], &ngrid[1], &ngrid[2]);
  getline(&line, &len, f);
  free(line);
  return f;
}

void ReadDatData(FILE *f, int ngrid[3], int dim, Cube *c)
{
  size_t len = 0;
  char *line = malloc(100 * sizeof(char));
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
  free(line);
  free(data);
  fclose(f);
}

double Surround(int atom, Cube *c, double r, int flag)
{
  int p[3], size;
  CubeR2BoxI(c, c->atoms[atom].coor, p);
  Cube *tmp = CubeSphericalRegion(c, p, r);
  size = CubeDataSize(tmp);
  double total = 0;
  for (int i = 0; i < size; ++i)
  {
    total += fabs(tmp->data[i]);
  }
  if (flag)
  {
    return total / size;
  }
  else
  {
    return total;
  }
}

int main()
{
  char files[5000][11], fname[20], cname[20];
  double celldm[3][3], inv[3][3], *norm, vsize[3][3];
  int nat = 0,
      ntyp = 0,
      num[30],
      z[30],
      ngrid[3],
      dim,
      dummy[3],
      nfile = FilesList(files, "dat");

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
    VecScale(vsize[i], 1.0/ngrid[i]);
  }

/* #pragma omp parallel for private(fname,cname) shared(dummy) */
  for(int file = 0; file < nfile; ++file)
  {
    sprintf(cname, "cube/%d.cube", atoi(files[file]));
    if (access(cname, F_OK) == -1)
    {
      FILE *f;
      char com[255] = "";
      int pos = 0;
      double rho, total = 0;
      Cube *c = CubeInit(nat, ngrid);
      CubeSetVoxels(c, vsize);

      sprintf(fname, "dat/%s", files[file]);
      f = ReadDatHeader(fname, dummy);
      ReadDatData(f, ngrid, dim, c);
      GetAtoms(atoi(files[file]), c, nat, z, num);
      /* for (int i = 0; i < CubeDataSize(c); ++i) */
      /* { */
      /*   total += fabs(c->data[i]); */
      /* } */
      for (int i = 0; i < c->nat; i++)
      {
        /* if(c->atoms[i].Z != 22) break; */
        if(Surround(i, c, 0.4, 1) > THR)
        {
          rho = Surround(i, c, 2.1, 0) ;
          /* printf("%d %lf %lf\n", atoi(files[file]), c->atoms[i].coor[2], rho /total); */
          printf("%d %lf %lf\n", atoi(files[file]), c->atoms[i].coor[2], rho);
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
      CubeWrite(c, cname);
      CubeDelete(c);
    }
  }
  return 0;
}

// vi: foldmethod=syntax
