#ifndef FILES_HELPERS
#define FILES_HELPERS 5678

#include <dirent.h>
#include <string.h>
#include <stdio.h>

int hasExtension(char const *name, char const *ext);
int FilesList(char files[5000][11], const char *ext);

#endif

