#include "cube.h"
#include <stdlib.h>
#include <math.h>

void usage()
{
  printf("Usage:\n");
  printf(" beautify infile outfile [threshold = 4 bohr]\n");
  printf(" threshold: the maximum distance to shift the data");
}
int main(int argc, char *argv[])
{
  int i, j, k, min_i[3] = {0, 0, 0};
  double sum, thr, min[3] = {1000, 1000, 1000};
  int *ind, dim;

  if(argc < 3 || argc > 4)
  {
    usage();
    return -1;
  }
  else if(argc == 4)
  {
    thr = atof(argv[3]);
  }
  else
  {
    thr = 4;
  }

  Cube *c = CubeRead(argv[1]);

  // looking one way, for now (+x, +y, +z, not -x, -y, -z)
  for (i = 0; i < 3; ++i)
  {
    dim = CubeDataSize(c)/c->ngrid[i];
    for (j = 0; j < c->ngrid[i]; ++j)
    {
      ind = LayerIndices(c, i, j);
      sum = 0;
      for (k = 0; k < dim; ++k)
      {
        sum += fabs(c->data[ind[k]]);
      }
      if(sum < min[i]) {
        min[i] = sum;
        min_i[i] = j;
      }
      if (j * c->vsize[i] > thr )
      {
        break;
      }
    }
    CubeRotateLayers(c, i, -1 * min_i[i]);
  }
  CubeWrite(c, argv[2]);
  return 0;
}
