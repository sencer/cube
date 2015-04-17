#include "files.h"

int hasExtension(char const *name, char const *ext)
{
  char dotext[8];
  sprintf(dotext, ".%s", ext);

  size_t ln = strlen(name),
         le = strlen(dotext);
  return ln > le && strcmp(name + ln - le, dotext) == 0;
}

int FilesList(char files[5000][11], const char *ext)
{
  DIR *directory = opendir(ext);
  struct dirent *ent;
  int nfile = 0;

  while ((ent = readdir(directory)) != NULL)
  {
    if(hasExtension(ent->d_name, ext))
    {
      strcpy(files[nfile], ent->d_name);
      nfile++;
    }
  }
  return nfile;
}
