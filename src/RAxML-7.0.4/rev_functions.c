/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 *
 *  This program is free software; you may redistribute it and/or modify its
 *  under the terms of the GNU General Public License as published by the Free
 *  Software Foundation; either version 2 of the License, or (at your option)
 *  any later version.
 *
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
 *  for more details.
 * 
 *
 *  For any other enquiries send an Email to Alexandros Stamatakis
 *  Alexandros.Stamatakis@epfl.ch
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "axml.h"


void tred2 (double *a, const int n, const int np, double *d, double *e)
{
#define a(i,j) a[(j-1)*np + (i-1)]
#define e(i)   e[i-1]
#define d(i)   d[i-1]
  int i, j, k, l;
  double f, g, h, hh, scale;
  for (i = n; i > 1; i--) {
    l = i-1;
    h = 0;
    scale = 0;
    if ( l > 1 ) {
      for ( k = 1; k <= l; k++ )
	scale += fabs(a(i,k));
      if (scale == 0) 
	e(i) = a(i,l);
      else {
	for (k = 1; k <= l; k++) {
	  a(i,k) /= scale;
	  h += a(i,k) * a(i,k);
	}
	f = a(i,l);
	g = -sqrt(h);
	if (f < 0) g = -g;
	e(i) = scale *g;
	h -= f*g;
	a(i,l) = f-g;
	f = 0;
	for (j = 1; j <=l ; j++) {
	  a(j,i) = a(i,j) / h;
	  g = 0;
	  for (k = 1; k <= j; k++)
	    g += a(j,k)*a(i,k);
	  for (k = j+1; k <= l; k++)
	    g += a(k,j)*a(i,k);
	  e(j) = g/h;
	  f += e(j)*a(i,j);
	}
	hh = f/(h+h);
	for (j = 1; j <= l; j++) {
	  f = a(i,j);
	  g = e(j) - hh * f;
	  e(j) = g;
	  for (k = 1; k <= j; k++) 
	    a(j,k) -= f*e(k) + g*a(i,k);
	}
      }
    } else 
      e(i) = a(i,l);
    d(i) = h;
  }
  d(1) = 0;
  e(1) = 0;
  for (i = 1; i <= n; i++) {
    l = i-1;
    if (d(i) != 0) {
      for (j = 1; j <=l; j++) {
	g = 0;
	for (k = 1; k <= l; k++)
	  g += a(i,k)*a(k,j);
	for (k=1; k <=l; k++)
	  a(k,j) -= g * a(k,i);
      }
    }
    d(i) = a(i,i);
    a(i,i) = 1;
    for (j=1; j<=l; j++)
      a(i,j) = a(j,i) = 0;
  }

  return;
#undef a
#undef e
#undef d
}


double pythag(double a, double b) {
  double absa = fabs(a), absb = fabs(b);
  return (absa > absb) ?
       absa * sqrt(1+ (absb/absa)*(absb/absa)) :
    absb == 0 ?
       0 :
       absb * sqrt(1+ (absa/absb)*(absa/absb));
}

void tqli(double *d, double *e, int n, int np, double *z) 
{
#define z(i,j) z[(j-1)*np + (i-1)]
#define e(i)   e[i-1]
#define d(i)   d[i-1]
  
  int i = 0, iter = 0, k = 0, l = 0, m = 0;
  double b = 0, c = 0, dd = 0, f = 0, g = 0, p = 0, r = 0, s = 0;
 
  for(i=2; i<=n; i++)
    e(i-1) = e(i);
  e(n) = 0;

  for (l = 1; l <= n; l++) 
    {
      iter = 0;
    labelExtra:
     
      for (m = l; (m < n); m++) 
	{
	  dd = fabs(d(m))+fabs(d(m+1));
	 
	  if (fabs(e(m))+dd == dd) 
	    break;
	}
     
      if (m != l) 
	{
	  assert(iter < 30); 
	   
	  iter++;
	  g = (d(l+1)-d(l))/(2*e(l));
	  r = pythag(g,1.);
	  g = d(m)-d(l)+e(l)/(g+(g<0?-r:r));
	  s = 1; 
	  c = 1;
	  p = 0;
	 
	  for (i = m-1; i>=l; i--) 
	    {
	      f = s*e(i);
	      b = c*e(i);
	      r = pythag(f,g);
	     
	      e(i+1) = r;
	      if (r == 0) 
		{
		  d (i+1) -= p;
		  e (m) = 0;
		  
		  goto labelExtra;
		}
	      s = f/r;
	      c = g/r;
	      g = d(i+1)-p;
	      r = (d(i)-g)*s + 2*c*b;
	      p = s*r;
	      d(i+1) = g + p;
	      g = c*r - b;
	      for (k=1; k <= n; k++) 
		{
		  f = z(k,i+1);
		  z(k,i+1) = s * z(k,i) + c*f;
		  z(k,i) = c * z(k,i) - s*f;
		}
	    }
	  d(l) -= p;
	  e(l) = g;
	  e(m) = 0;
	  
	  goto labelExtra;
	}
    }
 
  return;
#undef z
#undef e
#undef d
  
}
