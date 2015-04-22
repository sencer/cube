#include "cube.h"
#include <stdlib.h>
#include <string.h>

void usage()
{
  printf("Usage:\n");
  printf(" beautify [-t threshold] [-trim] infile outfile\n");
  printf(" -t threshold: threshold for trimming and rotating\n");
  printf(" -trim do not trim\n");
}

int main(int argc, char *argv[])
{
  int i, trim = 1;
  double thr = 7.E-3;

  if(argc < 3 || argc > 6)
  {
    usage();
    return -1;
  }
  else
  {
    for (i = 1; i < argc; ++i)
    {
      if(strcmp(argv[i], "-trim") == 0)
      {
        trim = 0;
      }
      else if(strcmp(argv[i], "-t") == 0)
      {
        thr = atof(argv[i+1]);
      }
    }
  }
  Cube *c = CubeRead(argv[argc - 2]);
  CubeBeautify(c, thr);
  if(trim) CubeTrim(&c, thr);
  CubeWrite(c, argv[argc - 1]);
  CubeDelete(c);
  return 0;
}
