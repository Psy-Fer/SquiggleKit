#include <stdlib.h>


typedef struct Path
{
  int k;
  int *px;
  int *py;
} Path;


double std(double *x, double *y, int n, int m, double *cost, int squared);
int path(double *cost, int n, int m, int startx, int starty, Path *p);
void subsequence(double *x, double *y, int n, int m, double *cost);
int subsequence_path(double *cost, int n, int m, int starty, Path *p);
