#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cube.h"


double linear(double a, double b, double weight)
{
  return (1 - weight) * a + weight * b;
}

int main(int argc, char *argv[])
{
  if (argc != 4)
  {
    printf("Usage:\n\n");
    printf("interpolate num_images file1 file2");
  }
  else
  {
    int nimage = atoi(argv[1]),
        f = atoi(argv[2]),
        l = atoi(argv[3]),
        k, j = 0, i;
    char fname[20];

    sprintf(fname, "%d.cube", f);
    Cube *first = CubeRead(fname);
    sprintf(fname, "%d.cube", l);
    Cube *last  = CubeRead(fname);
    Cube *image = CubeInitFrom(first);

    for (k = 1; k <= nimage; ++k)
    {
      double w = 1.0 * k / (nimage + 1);

      for(i = 0; i < first->nat; i++)
      {
        image->atoms[i].Z = first->atoms[i].Z;
        for(j=0; j<3; j++)
        {
          image->atoms[i].coor[j] = linear(first->atoms[i].coor[j],
                                           last->atoms[i].coor[j], w);
        }
      }

      int ndata = CubeDataSize(first);
      for (i = 0; i < ndata; ++i)
      {
        image->data[i] = linear(first->data[i], last->data[i], w);
      }

      sprintf(fname, "%d.cube", (int) (f + (l-f) * w));
      CubeWrite(image, fname);
    }
    CubeDelete(image);
    CubeDelete(first);
    CubeDelete(last);
  }
  return 0;
}
