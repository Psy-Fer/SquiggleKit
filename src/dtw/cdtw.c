/*  
    This code is written by Davide Albanese <davide.albanese@gmail.com>.
    (C) mlpy Developers.

    This program is free software: you can redistribute it and/or modify
    it underthe terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include "cdtw.h"

double
min3(double a, double b, double c)
{
  double min;
  
  min = a;
  if (b < min)
    min = b;
  if (c < min)
    min = c;
  return min;
}

// Paliwal adjustment window used for restricting the warping function
// r: window length
int
paliwal_window(int i, int j, int n, int m, int r)
{
  double s, f;
  
  s = ((double) m) / n;
  f = fabs(i - (((double) j) / s));
  
  if (f <= r)
    return 1;
  else
    return 0;
}

// euclidean distance
double e_dist(double x, double y)
{
  return fabs(x - y);
}

// squared euclidean distance
double se_dist(double x, double y)
{
  return pow(x - y, 2);
}


// Returns the (unnormalized) minimum-distance warp path 
// between time series x and y and the cost matrix C
double 
std(double *x, double *y, int n, int m, double *cost, int squared)
{
  int i, j;
  double (*dist)(double, double);
  
  if (squared == 0)
    dist = &e_dist;
  else
    dist = &se_dist;
  
  cost[0] = (*dist)(x[0], y[0]);
  
  for (i=1; i<n; i++)
    cost[i*m] = (*dist)(x[i], y[0]) + cost[(i-1)*m];
  
  for (j=1; j<m; j++)
    cost[j] = (*dist)(x[0], y[j]) + cost[(j-1)];
  
  for (i=1; i<n; i++)
    for (j=1; j<m; j++)
      cost[i*m+j] = (*dist)(x[i], y[j]) + 
	min3(cost[(i-1)*m+j], cost[(i-1)*m+(j-1)], cost[i*m+(j-1)]);
  
  return cost[n*m-1];
}

// Compute the warp path starting at cost[startx, starty]
// If startx = -1 -> startx = n-1; if starty = -1 -> starty = m-1
int
path(double *cost, int n, int m, int startx, int starty, Path *p)
{
  int i, j, k, z1, z2;
  int *px;
  int *py;
  double min_cost;
  
  if ((startx >= n) || (starty >= m))
    return 0;
  
  if (startx < 0)
    startx = n - 1;
  
  if (starty < 0)
    starty = m - 1;
      
  i = startx;
  j = starty;
  k = 1;
  
  // allocate path for the worst case
  px = (int *) malloc ((startx+1) * (starty+1) * sizeof(int));
  py = (int *) malloc ((startx+1) * (starty+1) * sizeof(int));
  
  px[0] = i;
  py[0] = j;
  
  while ((i > 0) || (j > 0))
    {
      if (i == 0)
	j--;
      else if (j == 0)
	i--;
      else
	{
	  min_cost = min3(cost[(i-1)*m+j],
			  cost[(i-1)*m+(j-1)], 
			  cost[i*m+(j-1)]);
	  
	  if (cost[(i-1)*m+(j-1)] == min_cost)
	    {
	      i--;
	      j--;
	    }
	  else if (cost[i*m+(j-1)] == min_cost)
	    j--;
	  else
	    i--;
	}
      
      px[k] = i;
      py[k] = j;
      k++;      
    }
  
  p->px = (int *) malloc (k * sizeof(int));
  p->py = (int *) malloc (k * sizeof(int));
  for (z1=0, z2=k-1; z1<k; z1++, z2--)
    {
      p->px[z1] = px[z2];
      p->py[z1] = py[z2];
    }
  p->k = k;
  
  free(px);
  free(py);
  
  return 1;
}


//
void
subsequence(double *x, double *y, int n, int m, double *cost)
{
  int i, j;
    
  cost[0] = fabs(x[0]-y[0]);
  
  for (i=1; i<n; i++)
    cost[i*m] = fabs(x[i]-y[0]) + cost[(i-1)*m];
  
  for (j=1; j<m; j++)
    cost[j] = fabs(x[0]-y[j]); // subsequence variation: D(0,j) := c(x0, yj)
  
  for (i=1; i<n; i++)
    for (j=1; j<m; j++)
      cost[i*m+j] = fabs(x[i]-y[j]) +
	min3(cost[(i-1)*m+j], cost[(i-1)*m+(j-1)], cost[i*m+(j-1)]);

}

  
int 
subsequence_path(double *cost, int n, int m, int starty, Path *p)
{
  int i, z1, z2;
  int a_star;
  int *tmpx, *tmpy;

  // find path
  if (!path(cost, n, m, -1, starty, p))
    return 0;
  
  // find a_star
  a_star = 0;
  for (i=1; i<p->k; i++)
    if (p->px[i] == 0)
      a_star++;
    else
      break;
  
  // rebuild path
  tmpx = p->px;
  tmpy = p->py;
  p->px = (int *) malloc ((p->k-a_star) * sizeof(int));
  p->py = (int *) malloc ((p->k-a_star) * sizeof(int));
  for (z1=0, z2=a_star; z2<p->k; z1++, z2++)
    {
      p->px[z1] = tmpx[z2];
      p->py[z1] = tmpy[z2];
    }
  p->k = p->k-a_star;
  
  free(tmpx);
  free(tmpy);
  
  return 1;
}
