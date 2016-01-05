#define _GNU_SOURCE
#define THR 6.E-3

#include "cube.h"
#include "files.h"
#include "3d.h"
#include <unistd.h>


void ReadInput(int *nat, int *ntyp, int *num, int *z, double celldm[3][3])
{
  // read cube.in for system information.
  // file should be in the following format
  // nat ntyp
  // n_typ_1 z_typ_1
  // n_typ_2 z_typ_2
  // ..
  // n_typ_i z_typ_i
  // 3-lines of cell dimension information from cp.cel
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
  fclose(f);
  free(line);
}

void GetAtoms(int index, Cube *c, int nat, int *z, int *num)
{
  // Get ionic coordinates for frame *index* from cp.pos
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
  free(line);
}

void ReadDatHeader(FILE *f, int *ngrid)
{
  size_t len = 0;
  char *line = malloc(100 * sizeof(char));
  /* FILE *f = fopen(fname, "rb"); */

  getline(&line, &len, f);
  while(strstr(line, "INFO nr1") == NULL)
  {
    getline(&line, &len, f);
  }
  sscanf(line, "    <INFO nr1=\"%d\" nr2=\"%d\" nr3=\"%d\"/>",
      &ngrid[0], &ngrid[1], &ngrid[2]);
  getline(&line, &len, f);
  free(line);
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
  // get list of dat/*.dat, and sort them
  char files[5000][11];
  int nfile = FilesList(files, "dat");
  qsort(&files[0], nfile, sizeof(files[0]), CompStrInt);
  if(nfile < 1) return -1;  // exit if no files

  // read following system information from cube.in
  // number of atoms -> nat
  // number of atomic types -> ntyp
  // number of each atom in system -> num[]
  // proton number of each atom in the same order with num[] -> z[]
  // 3x3 matrix of the cell dimensions -> celldm[3][3]
  int nat = 0, ntyp = 0, num[30], z[30];
  double celldm[3][3];
  ReadInput(&nat, &ntyp, num, z, celldm);

  // we will need
  // the inverse of celldm matrix -> inv[3][3]
  // and lengths of cell vectors -> norm[3]
  double inv[3][3];
  double norm[3] = {
    VecLen(celldm[0]),
    VecLen(celldm[1]),
    VecLen(celldm[2])
  };
  Inv3D(celldm, inv);

  // a string for
  // the name of the dat file to read -> fname[20]
  // the name of the cube file to write -> cname[20]
  // a FILE for the file i/o
  char fname[20], cname[20];
  FILE *f;

  // read the grid information from the files[0]
  int ngrid[3];
  sprintf(fname, "dat/%s", files[0]);
  f = fopen(fname, "r");
  ReadDatHeader(f, ngrid);
  fclose(f);

  // let's calculate ngrid[0] * ngrid[1] once and for all.
  int dim = ngrid[0] * ngrid[1];

  // calculate the voxel sizes
  double vsize[3][3];
  memcpy(vsize, celldm, sizeof(celldm));
  for(int i = 0; i < 3; ++i)
  {
    VecScale(vsize[i], 1.0/ngrid[i]);
  }

  int dummy[3];

#pragma omp parallel for private(f, fname, cname) shared(dummy)
  for(int file = 0; file < nfile; ++file)
  {
    // set cube file name
    sprintf(cname, "cube/%d.cube", atoi(files[file]));
    // if the file exists, skip
    if (access(cname, F_OK) == -1)
    {
      // create a Cube to read dat file into
      Cube *c = CubeInit(nat, ngrid);
      CubeSetVoxels(c, vsize);

      // Read the next dat file (header info is read to skip this part)
      sprintf(fname, "dat/%s", files[file]);
      f = fopen(fname, "r");
      ReadDatHeader(f, dummy);
      ReadDatData(f, ngrid, dim, c);
      fclose(f);

      // get ionic positions for the corresponding frame from cp.pos
      GetAtoms(atoi(files[file]), c, nat, z, num);

      int pos = 0;
      char com[255] = "";
      double rho;

      for (int i = 0; i < c->nat; i++)
      {
        if(Surround(i, c, 0.4, 1) > THR && c->atoms[i].Z == 22)
        {
          rho = Surround(i, c, 2.1, 0) ;
          printf("%d %lf %lf\n", atoi(files[file]), c->atoms[i].coor[2], rho);
          pos += sprintf(com + pos, "%d ", i);
        }
        fflush(stdout);
      }

      CubeBeautify(c, THR);
      CubeTrim(&c, THR);

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
