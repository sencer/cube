#include "cube.h"
#include <stdlib.h>
#include <math.h>

int main(int argc, char *argv[])
{
  Cube *c = CubeRead(argv[1]);
  double max, thr = atof(argv[3]);
  int *c_i, p[6], i, k, j, n, dim;

  for (i = 0; i < 6; ++i)
  {
    dim = CubeDataSize(c)/c->ngrid[i%3];
    for (j = 0; j < c->ngrid[i%3]; ++j)
    {
      n = i < 3 ? j : c->ngrid[i-3] - 1 - j;
      c_i = LayerIndices(c, i % 3, n);
      max = 0;
      for (k = 0; k < dim; ++k)
      {
        if(max < fabs(c->data[c_i[k]]))
          max = fabs(c->data[c_i[k]]);
      }
      free(c_i);
      if (max > thr)
      {
        p[i] = i<3?j-1:c->ngrid[i-3]-j;
        break;
      }
    }
  }
  printf("%d %d %d %d %d %d\n", *p, *(p+1), *(p+2), *(p+3), *(p+4), *(p+5));
  Cube *nc = CubeGetRegion(c, p, p+3);
  Atom *ptr = nc->atoms;
  nc->nat = c->nat;
  nc->atoms = c->atoms;
  CubeWrite(nc, argv[2]);
  nc->atoms = ptr;
  CubeDelete(nc);
  CubeDelete(c);
  return 0;
}
