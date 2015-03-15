#include "cube.h"
#include <stdlib.h>
#include <string.h>
#include <dirent.h>
#include <math.h>


void SwapPtr(char **arg1, char **arg2)
{
    char *tmp = *arg1;
    *arg1 = *arg2;
    *arg2 = tmp;
}

void QuickSort(char *args[], unsigned int len)
{
    unsigned int i, pvt=0;

    if (len <= 1)
        return;

    // swap a randomly selected value to the last node
    SwapPtr(args+((unsigned int)rand() % len), args+len-1);

    // reset the pivot index to zero, then scan
    for (i=0;i<len-1;++i)
    {
        if (atoi(args[i]) < atoi(args[len-1]))
            SwapPtr(args+i, args+pvt++);
    }

    // move the pivot value into its place
    SwapPtr(args+pvt, args+len-1);

    // and invoke on the subsequences. does NOT include the pivot-slot
    QuickSort(args, pvt++);
    QuickSort(args+pvt, len - pvt);
}

int is_extension(char const *ext, char const *name)
{
  size_t len = strlen(name), elen = strlen(ext);
  return len > elen && strcmp(name + len - elen, ext) == 0;
}

int FilesList(char *files[])
{
  DIR *directory = opendir("cube");
  struct dirent *ent;
  int nfile = 0;

  while ((ent = readdir(directory)) != NULL)
  {
    if(is_extension(".cube", ent->d_name))
    {
      files[nfile] = ent->d_name;
      nfile++;
    }
  }
  QuickSort(files, nfile);
  return nfile;
}

int main()
{
  char *files[5000];
  int n = FilesList(files), fi, li, step = 10, im;
  char fname[20], lname[20];
#pragma omp parallel for shared(step, files) private(fname, lname, fi, li, im)
  for(int i = 0; i < n - 1; ++i)
  {
    sprintf(fname, "cube/%s", files[i]);
    sprintf(lname, "cube/%s", files[i+1]);
    Cube *first = CubeRead(fname);
    Cube *last  = CubeRead(lname);
    fi = atoi(files[i]);
    li = atoi(files[i+1]);
    im = (li - fi - 1) / step;
    Cube **cubes = CubeInterpolate(first, last, im);
    for(int j = 0; j < im; ++j)
    {
      sprintf(fname, "%d.cube", fi + ( j + 1 ) * step);
      CubeWrite(cubes[j], fname);
    }
    free(cubes);
  }
  return 0;
}
