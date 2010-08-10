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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with thousands 
 * of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef WIN32
#include <unistd.h>
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "axml.h"

static const double MNBRAK_GOLD =    1.618034;
static const double MNBRAK_TINY =      1.e-20;
static const double MNBRAK_GLIMIT =     100.0;
static const double BRENT_ZEPS  =      1.e-5;
static const double BRENT_CGOLD =   0.3819660;

extern int optimizeRatesInvocations;
extern int optimizeRateCategoryInvocations;
extern int optimizeAlphaInvocations;
extern int optimizeInvarInvocations;
extern double masterTime;
extern char ratesFileName[1024];
extern char workdir[1024];
extern char run_id[128];
extern char lengthFileName[1024];
extern char lengthFileNameModel[1024];



#ifdef _USE_PTHREADS
extern int rateCategoryJobs;
#endif

/*********************FUNCTIONS FOOR EXACT MODEL OPTIMIZATION UNDER GTRGAMMA ***************************************/


static double evaluateInvar(tree *tr, double invar, int model)
{
  double result;
  
  tr->invariants[model] = invar;      
#ifdef _LOCAL_DATA
  masterBarrier(THREAD_COPY_INVARIANTS ,tr);
#endif  
  result = evaluatePartitionGeneric(tr, tr->start, model);     
  return result;   
}

static int brakInvar(double *param, double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double lim_inf, 
		     double lim_sup, tree *tr, int model)
{
   double ulim,u,r,q,fu,dum;

   u = 0.0;

   *param = *ax;
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fa = -evaluateInvar(tr, *param, model);

   *param = *bx;
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;  
   *fb = -evaluateInvar(tr, *param, model);

   if (*fb > *fa) 
     {
       SHFT(dum,*ax,*bx,dum)
       SHFT(dum,*fb,*fa,dum)
     }

   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   *param = *cx;
   if(*param > lim_sup) *param = *cx = lim_sup;
   if(*param < lim_inf) *param = *cx = lim_inf;
   *fc = -evaluateInvar(tr, *param, model); 

   while (*fb > *fc) 
     {        
       if(*ax > lim_sup) *ax = lim_sup;
       if(*ax < lim_inf) *ax = lim_inf;
       if(*bx > lim_sup) *bx = lim_sup;
       if(*bx < lim_inf) *bx = lim_inf;
       if(*cx > lim_sup) *cx = lim_sup;
       if(*cx < lim_inf) *cx = lim_inf;      

       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(fabs(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if(u > lim_sup) u = lim_sup;
       if(u < lim_inf) u = lim_inf;
       if(ulim > lim_sup) ulim = lim_sup;
       if(ulim < lim_inf) ulim = lim_inf;

       if ((*bx-u)*(u-*cx) > 0.0) 
	 {
	   *param = u;	   
	   fu = -evaluateInvar(tr, *param, model);
	   if (fu < *fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       *fa=(*fb);
	       *fb=fu;	      
	       return(0);
	     } 
	   else if (fu > *fb) 
	     {
	       *cx=u;
	       *fc=fu;		       
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   *param = u;
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}	
	   fu= -evaluateInvar(tr, *param, model);
	 } 
       else if ((*cx-u)*(u-ulim) > 0.0) 
	 {
	   *param = u;	   
	   fu = -evaluateInvar(tr, *param, model);
	   if (fu < *fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))	        	      
	       SHFT(*fb,*fc,fu, -evaluateInvar(tr, *param, model))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= 0.0) 
	 {
	   u = ulim;
	   *param = u;	    
	   fu = -evaluateInvar(tr, *param, model);
	 } 
       else 
	 {
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   *param = u;
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu = -evaluateInvar(tr, *param, model);
	 }
       SHFT(*ax,*bx,*cx,u)
       SHFT(*fa,*fb,*fc,fu)


     }
   
   return(0);
}


static double brentInvar(double ax, double bx, double cx, double fb, double tol, double *xmin, int model, tree *tr)
{
  int iter;
  double a,b,d = 0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x = w = v = bx;
  fw = fv = fx = fb;

  for(iter = 1; iter <= ITMAX; iter++)
    {
      xm = 0.5 * (a + b);
      tol2 = 2.0 * (tol1 = tol * fabs(x) + BRENT_ZEPS);
      if(fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
	{
	  *xmin = x;
	  return -fx;
	}
      if(fabs(e) > tol1)
	{
	  r = (x - w) * (fx - fv);
	  q = (x - v) * (fx - fw);
	  p = (x -v ) * q - (x - w) * r;
	  q = 2.0 * (q - r);
	  if(q > 0.0)
	    p = -p;
	  q = fabs(q);
	  etemp = e;
	  e = d;
	  if((fabs(p) >= fabs(0.5 * q * etemp)) || (p <= q * (a-x)) || (p >= q * (b - x)))
	    d = BRENT_CGOLD * (e = (x >= xm ? a - x : b - x));
	  else
	    {
	      d = p / q;
	      u = x + d;
	      if( u - a < tol2 || b - u < tol2)
		d = SIGN(tol1, xm - x);
	    }
	}
      else
	{
	  d = BRENT_CGOLD * (e = (x >= xm ? a - x: b - x));
	}
      u = (fabs(d) >= tol1 ? x + d: x +SIGN(tol1, d));
      
      fu = -evaluateInvar(tr, u, model);

      if(fu <= fx)
	{
	  if(u >= x)
	    a = x;
	  else
	    b = x;
	  SHFT(v,w,x,u)
	    SHFT(fv,fw,fx,fu)
	    }
      else
	{
	  if(u < x)
	    a = u;
	  else
	    b = u;
	  if(fu <= fw || w == x)
	    {
	      v = w;
	      w = u;
	      fv = fw;
	      fw = fu;
	    }
	  else
	    {
	      if(fu <= fv || v == x || v == w)
		{
		  v = u;
		  fv = fu;
		}
	    }	    
	}
    }

  printf("\n. Too many iterations in BRENT !");
  exit(-1);
  return(-1);

}


static void optInvar(tree *tr, analdef *adef, double modelEpsilon, int model)
{
  double param, a, b, c, fa, fb, fc, x;
  double lim_inf = INVAR_MIN;
  double lim_sup = INVAR_MAX;
  double startLH;
  double startInvar = tr->invariants[model];
  double endInvar;  

  if(adef->useMultipleModel)         
    startLH = evaluatePartitionGeneric(tr, tr->start, model);             
  else
    startLH = tr->likelihood;

  a = tr->invariants[model] + 0.1;
  b = tr->invariants[model] - 0.1; 
  
  if(b < lim_inf) b = lim_inf;
  
  brakInvar(&param, &a, &b, &c, &fa, &fb, &fc, lim_inf, lim_sup, tr, model);
          
  endInvar = brentInvar(a, b, c, fb, modelEpsilon, &x, model, tr);  

  if(startLH > endInvar)
    {    
      tr->invariants[model] = startInvar;
#ifdef _LOCAL_DATA
      masterBarrier(THREAD_COPY_INVARIANTS ,tr);
#endif 
      evaluateInvar(tr, tr->invariants[model], model);      
    } 
}






/**********************************************************************************************************/
/* ALPHA PARAM ********************************************************************************************/

static double evaluateAlpha(tree *tr, double alpha, int model)
{
  double result;

  tr->alphas[model] = alpha;    
  
  makeGammaCats(model, tr->alphas, tr->gammaRates); 
#ifdef _LOCAL_DATA
  masterBarrier(THREAD_COPY_GAMMA_RATES, tr);
#endif

  result = evaluateGenericInitravPartition(tr, tr->start, model);     
  return result;
}


static int brakAlpha(double *param, double *ax, double *bx, double *cx, double *fa, double *fb, double *fc, double lim_inf, double lim_sup, tree *tr, int model)
{
   double ulim,u,r,q,fu,dum;

   u = 0.0;

   *param = *ax;
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   *fa = -evaluateAlpha(tr, *param, model);

   *param = *bx;
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;  
   *fb = -evaluateAlpha(tr, *param, model);

   if (*fb > *fa) 
     {
       SHFT(dum,*ax,*bx,dum)
       SHFT(dum,*fb,*fa,dum)
     }

   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   *param = *cx;
   if(*param > lim_sup) *param = *cx = lim_sup;
   if(*param < lim_inf) *param = *cx = lim_inf;
   *fc = -evaluateAlpha(tr, *param, model); 

   while (*fb > *fc) 
     {        
       if(*ax > lim_sup) *ax = lim_sup;
       if(*ax < lim_inf) *ax = lim_inf;
       if(*bx > lim_sup) *bx = lim_sup;
       if(*bx < lim_inf) *bx = lim_inf;
       if(*cx > lim_sup) *cx = lim_sup;
       if(*cx < lim_inf) *cx = lim_inf;      

       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(fabs(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if(u > lim_sup) u = lim_sup;
       if(u < lim_inf) u = lim_inf;
       if(ulim > lim_sup) ulim = lim_sup;
       if(ulim < lim_inf) ulim = lim_inf;

       if ((*bx-u)*(u-*cx) > 0.0) 
	 {
	   *param = u;	   
	   fu = -evaluateAlpha(tr, *param, model);
	   if (fu < *fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       *fa=(*fb);
	       *fb=fu;	      
	       return(0);
	     } 
	   else if (fu > *fb) 
	     {
	       *cx=u;
	       *fc=fu;		       
	       return(0);
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   *param = u;
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}	
	   fu= -evaluateAlpha(tr, *param, model);
	 } 
       else if ((*cx-u)*(u-ulim) > 0.0) 
	 {
	   *param = u;	   
	   fu = -evaluateAlpha(tr, *param, model);
	   if (fu < *fc) 
	     {
	       SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))	        	      
	       SHFT(*fb,*fc,fu, -evaluateAlpha(tr, *param, model))
	     }
	 } 
       else if ((u-ulim)*(ulim-*cx) >= 0.0) 
	 {
	   u = ulim;
	   *param = u;	    
	   fu = -evaluateAlpha(tr, *param, model);
	 } 
       else 
	 {
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   *param = u;
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}
	   fu = -evaluateAlpha(tr, *param, model);
	 }
       SHFT(*ax,*bx,*cx,u)
       SHFT(*fa,*fb,*fc,fu)


     }
   
   return(0);
}


static double brentAlpha(double ax, double bx, double cx, double fb, double tol, double *xmin, int model, tree *tr)
{
  int iter;
  double a,b,d = 0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx); 

  x = w = v = bx;
  fw = fv = fx = fb; 

  for(iter = 1; iter <= ITMAX; iter++)
    {
      xm = 0.5 * (a + b);
      tol2 = 2.0 * (tol1 = tol * fabs(x) + BRENT_ZEPS);

      if(fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
	{
	  *xmin = x;
	  return -fx;
	}

      if(fabs(e) > tol1)
	{
	  r = (x - w) * (fx - fv);
	  q = (x - v) * (fx - fw);
	  p = (x - v) * q - (x - w) * r;	 

	  q = 2.0 * (q - r);

	  if(q > 0.0)
	    p = -p;

	  q = fabs(q);
	  etemp = e;
	  e = d;

	  if(fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a-x) || p >= q * (b - x))
	    {	     
	      d = BRENT_CGOLD * (e = (x >= xm ? a - x : b - x));
	    }
	  else
	    {	     
	      d = p / q;
	      u = x + d;
	      if( u - a < tol2 || b - u < tol2)
		d = SIGN(tol1, xm - x);
	    }
	}
      else
	{	 
	  d = BRENT_CGOLD * (e = (x >= xm ? a - x: b - x));
	}
      
      u = (fabs(d) >= tol1 ? x + d: x + SIGN(tol1, d));
          

      fu = -evaluateAlpha(tr, u, model);

      if(fu <= fx)
	{
	  if(u >= x)
	    a = x;
	  else
	    b = x;

	  SHFT(v,w,x,u)
	  SHFT(fv,fw,fx,fu)
	}
      else
	{
	  if(u < x)
	    a = u;
	  else
	    b = u;
	  if(fu <= fw || w == x)
	    {
	      v = w;
	      w = u;
	      fv = fw;
	      fw = fu;
	    }
	  else
	    {
	      if(fu <= fv || v == x || v == w)
		{
		  v = u;
		  fv = fu;
		}
	    }	    
	}
    }

  printf("\n. Too many iterations in BRENT !");
  exit(-1);
  return(-1);

}




static void optAlpha(tree *tr, analdef *adef, double modelEpsilon, int model)
{
  double param, a, b, c, fa, fb, fc, x;
  double lim_inf = ALPHA_MIN;
  double lim_sup = ALPHA_MAX;
  double startLH;
  double startAlpha = tr->alphas[model];
  double endAlpha;

  if(adef->useMultipleModel)
    startLH = evaluateGenericInitravPartition(tr, tr->start, model);
  else
    startLH = tr->likelihood;

  

  a = tr->alphas[model] + 0.1;
  b = tr->alphas[model] - 0.1; 
  
  if(b < lim_inf) b = lim_inf;
  
  brakAlpha(&param, &a, &b, &c, &fa, &fb, &fc, lim_inf, lim_sup, tr, model);
             
  endAlpha = brentAlpha(a, b, c, fb, modelEpsilon, &x, model, tr);  

  if(startLH > endAlpha)
    {    
      tr->alphas[model] = startAlpha;
      evaluateAlpha(tr, tr->alphas[model], model);      
    }  
  
}


/*******************************************************************************************************************/
/*******************RATES ******************************************************************************************/

static void setRateModel(tree *tr, int model, double rate, int position)
{
  int dataType;

  /*if(adef->model == M_PROTCAT || adef->model == M_PROTGAMMA)
    dataType = AA_DATA;
  else
  dataType = DNA_DATA;*/

  /*dataType = tr->dataType[model]; */
  dataType = tr->partitionData[model].dataType;

  switch(dataType)
    {
    case AA_DATA:
      assert(position >= 0 && position < AA_RATES);
      tr->initialRates_AA[model * AA_RATES  + position] = rate;
      break;
    case DNA_DATA:
      assert(position >= 0 && position < DNA_RATES);
      tr->initialRates_DNA[model * DNA_RATES + position] = rate;
      break;
    default:
      assert(0);
    }
}



static double getRateModel(tree *tr, int model, int position)
{
  int dataType;

  dataType = tr->partitionData[model].dataType;

  switch(dataType)
    {
    case AA_DATA:
      assert(position >= 0 && position < AA_RATES);
      return tr->initialRates_AA[model * AA_RATES  + position];
      break;
    case DNA_DATA:
      assert(position >= 0 && position < DNA_RATES);
      return tr->initialRates_DNA[model * DNA_RATES + position];
      break;
    default:
      assert(0);
    }

  return 0.0;
}







static double evaluateRate(tree *tr, int i, double rate, analdef *adef, int model)
{
  double result;

  setRateModel(tr, model, rate, i);  
  initReversibleGTR(tr, adef, model);
#ifdef _LOCAL_DATA 
  masterBarrier(THREAD_COPY_REVERSIBLE, tr);  
#endif
  
  result = evaluateGenericInitravPartition(tr, tr->start, model);   

  return result; 
}

static double brentRates(double ax, double bx, double cx, double fb, double tol, double *xmin, int model, tree *tr, analdef *adef, int i)
{
  int iter;
  double a,b,d = 0.0, etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
  double e=0.0;

  a=((ax < cx) ? ax : cx);
  b=((ax > cx) ? ax : cx);
  x = w = v = bx;
  fw = fv = fx = fb;

  for(iter = 1; iter <= ITMAX; iter++)
    {
      xm = 0.5 * (a + b);
      tol2 = 2.0 * (tol1 = tol * fabs(x) + BRENT_ZEPS);
      if(fabs(x - xm) <= (tol2 - 0.5 * (b - a)))
	{
	  *xmin = x;
	  return -fx;
	}
      if(fabs(e) > tol1)
	{
	  r = (x - w) * (fx - fv);
	  q = (x - v) * (fx - fw);
	  p = (x -v ) * q - (x - w) * r;
	  q = 2.0 * (q - r);
	  if(q > 0.0)
	    p = -p;
	  q = fabs(q);
	  etemp = e;
	  e = d;
	  if((fabs(p) >= fabs(0.5 * q * etemp)) || (p <= q * (a-x)) || (p >= q * (b - x)))
	    d = BRENT_CGOLD * (e = (x >= xm ? a - x : b - x));
	  else
	    {
	      d = p / q;
	      u = x + d;
	      if( u - a < tol2 || b - u < tol2)
		d = SIGN(tol1, xm - x);
	    }
	}
      else
	{
	  d = BRENT_CGOLD * (e = (x >= xm ? a - x: b - x));
	}
      u = (fabs(d) >= tol1 ? x + d: x +SIGN(tol1, d));
      
      fu = (-evaluateRate(tr, i, u, adef, model));    
    
      if(fu <= fx)
	{
	  if(u >= x)
	    a = x;
	  else
	    b = x;
	  SHFT(v,w,x,u)
	    SHFT(fv,fw,fx,fu)
	    }
      else
	{
	  if(u < x)
	    a = u;
	  else
	    b = u;
	  if(fu <= fw || w == x)
	    {
	      v = w;
	      w = u;
	      fv = fw;
	      fw = fu;
	    }
	  else
	    {
	      if(fu <= fv || v == x || v == w)
		{
		  v = u;
		  fv = fu;
		}
	    }	    
	}
    }

  printf("\n. Too many iterations in BRENT !");
  exit(-1);
  return(-1);
}


static int brakRates(double *param, double *ax, double *bx, double *cx, double *fa, double *fb, 
		     double *fc, double lim_inf, double lim_sup, 
		     tree *tr, int i, analdef *adef, int model)
{
   double ulim,u,r,q,fu,dum;

   u = 0.0;

   *param = *ax;

   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;
   
   *fa = -evaluateRate(tr, i, *param, adef, model);

   *param = *bx;
   if(*param > lim_sup) *param = lim_sup;
   if(*param < lim_inf) *param = lim_inf;

   *fb = -evaluateRate(tr, i, *param, adef, model);

   if (*fb > *fa) 
     {
       SHFT(dum,*ax,*bx,dum)
       SHFT(dum,*fa,*fb,dum)
     }

   *cx=(*bx)+MNBRAK_GOLD*(*bx-*ax);
   *param = *cx;
   if(*param > lim_sup) *param = *cx = lim_sup;
   if(*param < lim_inf) *param = *cx = lim_inf;
   *fc =  -evaluateRate(tr, i, *param, adef, model);

   while (*fb > *fc) 
     {                     
       if(*ax > lim_sup) *ax = lim_sup;
       if(*ax < lim_inf) *ax = lim_inf;
       if(*bx > lim_sup) *bx = lim_sup;
       if(*bx < lim_inf) *bx = lim_inf;
       if(*cx > lim_sup) *cx = lim_sup;
       if(*cx < lim_inf) *cx = lim_inf;
       
       r=(*bx-*ax)*(*fb-*fc);
       q=(*bx-*cx)*(*fb-*fa);
       u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
               (2.0*SIGN(MAX(fabs(q-r),MNBRAK_TINY),q-r));
       ulim=(*bx)+MNBRAK_GLIMIT*(*cx-*bx);
       
       if(u > lim_sup) u = lim_sup;
       if(u < lim_inf) u = lim_inf;
       if(ulim > lim_sup) ulim = lim_sup;
       if(ulim < lim_inf) ulim = lim_inf;

       if ((*bx-u)*(u-*cx) > 0.0)
	 {
	   *param = u;	 
	   fu = -evaluateRate(tr, i, *param, adef, model);
	   if (fu < *fc) 
	     {
	       *ax=(*bx);
	       *bx=u;
	       *fa=(*fb);
	       *fb=fu;	       
	       return(0);
	     } 
	   else 
	     {
	       if (fu > *fb) 
		 {
		   *cx=u;
		   *fc=fu;			  
		   return(0);
		 }
	     }
	   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
	   *param = u;
	   if(*param > lim_sup) {*param = u = lim_sup;}
	   if(*param < lim_inf) {*param = u = lim_inf;}	   
	   fu= -evaluateRate(tr, i, *param, adef, model);
	 } 
       else 
	 {
	   if ((*cx-u)*(u-ulim) > 0.0) 
	     {
	       *param = u;	       
	       fu = -evaluateRate(tr, i, *param, adef, model);
	       if (fu < *fc) 
		 {
		   SHFT(*bx,*cx,u,*cx+MNBRAK_GOLD*(*cx-*bx))		  
		   SHFT(*fb,*fc,fu, -evaluateRate(tr, i, *param, adef, model))
		     }
	     } 
	   else
	     {
	       if ((u-ulim)*(ulim-*cx) >= 0.0) 
		 {
		   u = ulim;
		   *param = u;		  
		   fu = -evaluateRate(tr, i, *param, adef, model);		   
		 } 
	       else 
		 {
		   u=(*cx)+MNBRAK_GOLD*(*cx-*bx);
		   *param = u;
		   if(*param > lim_sup) {*param = u = lim_sup;}
		   if(*param < lim_inf) {*param = u = lim_inf;}		   
		   fu = -evaluateRate(tr, i, *param, adef, model);
		 }
	     }
	 }     
       SHFT(*ax,*bx,*cx,u)
	 SHFT(*fa,*fb,*fc,fu)
	 }

   
   return(0);
}


static double optRates(tree *tr, analdef *adef, double modelEpsilon, int model)
{
  int i, dataType, numberOfRates;
  double param, a, b, c, fa, fb, fc, x;
  double lim_inf = RATE_MIN;
  double lim_sup = RATE_MAX;
  double *startRates, *initialRates;
  double startLH;
  double endLH = unlikely;   
  
  /*dataType = tr->dataType[model];*/
  dataType = tr->partitionData[model].dataType;

  switch(dataType)
    {
    case AA_DATA:
      assert(0);
      /* should not get here */
      initialRates  = tr->initialRates_AA;
      numberOfRates = AA_RATES;
      break;
    case DNA_DATA:
      initialRates  = tr->initialRates_DNA;
      numberOfRates = DNA_RATES;
      break;
    default:
      assert(0);
    }

  startRates = (double *)malloc(sizeof(double) * numberOfRates);

  if(adef->useMultipleModel)
    startLH = evaluateGenericInitravPartition(tr, tr->start, model);
  else
    startLH = tr->likelihood;  

 

  for(i = 0; i < numberOfRates; i++)
    startRates[i] = initialRates[model * numberOfRates  + i];

  for(i = 0; i < numberOfRates; i++)
    {           
      a = initialRates[model * numberOfRates + i] + 0.1;
      b = initialRates[model * numberOfRates + i] - 0.1;
      
      if(a < lim_inf) a = lim_inf;
      if(a > lim_sup) a = lim_sup;

      if(b < lim_inf) b = lim_inf;
      if(b > lim_sup) b = lim_sup;    
      
      brakRates(&param, &a, &b, &c, &fa, &fb, &fc, lim_inf, lim_sup, tr, i, adef, model);                
      
      endLH = brentRates(a, b, c, fb, modelEpsilon, &x,  model, tr, adef, i);                   

      if(startLH > endLH)
	{       		
	  initialRates[model * numberOfRates + i] = startRates[i];       
	  
	  initReversibleGTR(tr, adef, model);
#ifdef _LOCAL_DATA
	  masterBarrier(THREAD_COPY_REVERSIBLE, tr);
#endif

	  if(adef->useMultipleModel)
	    startLH = evaluateGenericInitravPartition(tr, tr->start, model);
	  else
	    {
	      evaluateGenericInitrav(tr, tr->start);
	      startLH = tr->likelihood;
	    }
	}
      else	
	startLH = endLH;                     	
    } 

  free(startRates);
  return endLH;
}


/*****************************************************************************************************/

void resetBranches(tree *tr)
{
  nodeptr  p, q;
  int  nodes, i;
  
  nodes = tr->mxtips  +  3 * (tr->mxtips - 2);
  p = tr->nodep[1];
  while (nodes-- > 0) 
    {   
      for(i = 0; i < tr->numBranches; i++)
	p->z[i] = defaultz;
	
      q = p->next;
      while(q != p)
	{	
	  for(i = 0; i < tr->numBranches; i++)
	    q->z[i] = defaultz;	    
	  q = q->next;
	}
      p++;
    }
}

/*
static void cacheZs(tree *tr, double *zs)
{ 
  nodeptr  p;
  int  nodes, i;
  
  nodes = tr->mxtips  +  3 * (tr->mxtips - 2);
  p = tr->nodep[1];
  while (nodes-- > 0) 
    {
      for(i = 0; i <  tr->numBranches; i++)	
	zs[nodes * tr->numBranches + i] = p->z[i]; 	
      p++;
    }
} 


static void restoreZs(tree *tr, double *zs)
{ 
  nodeptr  p;
  int  nodes, i;
  
  nodes = tr->mxtips  +  3 * (tr->mxtips - 2);
  p = tr->nodep[1];
  while (nodes-- > 0) 
    {
      for(i = 0; i <  tr->numBranches; i++)
	{
	  p->z[i]       = zs[nodes * tr->numBranches + i]; 
	  p->back->z[i] = zs[nodes * tr->numBranches + i]; 
	}
      p++;
    }
    }*/



void modOpt(tree *tr, analdef *adef)
{
  double currentLikelihood;
  int i; 
  double modelEpsilon = MODEL_EPSILON;
  int model;  

  assert(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I); 

  for(model = 0; model < tr->NumberOfModels; model++)
    {     
      for(i = 0; i <  DNA_RATES; i++)
	tr->initialRates_DNA[model * DNA_RATES + i] = 0.5;
      for(i = 0; i <  AA_RATES; i++)
	tr->initialRates_AA[model * AA_RATES + i] = 0.5;

      if(adef->useInvariant)
	{
	  int lower, upper;
	  int count = 0;
	  int total = 0;
	  	  	 
	  lower = tr->partitionData[model].lower;
	  upper = tr->partitionData[model].upper;

	  
	  for(i = lower; i < upper; i++)
	    {
	      if(tr->invariant[i] < 4) 		
		count += tr->cdta->aliaswgt[i];		  		
	      total += tr->cdta->aliaswgt[i];
	    }
	  tr->invariants[model] = ((double)count)/((double) total);
#ifdef _LOCAL_DATA
	  masterBarrier(THREAD_COPY_INVARIANTS ,tr);
#endif 
	}   
       
      tr->alphas[model] = 1.0;     

      initReversibleGTR(tr, adef, model);      
#ifdef _LOCAL_DATA
      masterBarrier(THREAD_COPY_REVERSIBLE, tr);
#endif
      makeGammaCats(model, tr->alphas, tr->gammaRates); 
#ifdef _LOCAL_DATA
      masterBarrier(THREAD_COPY_GAMMA_RATES, tr);
#endif
    }

  resetBranches(tr);

  if(adef->mode != MEHRING_ALGO && !adef->rapidML_Addition)
    tr->start = tr->nodep[1];

  /* no need for individual models here, just an init on params equal for all partitions*/
  
  
  evaluateGenericInitrav(tr, tr->start);

  /*printf("1 %f\n", tr->likelihood);*/

  treeEvaluate(tr, 1.0);    

#ifdef _DEBUG_AA
  printf("I %f\n", tr->likelihood);
#endif
  /*printf("2 %f\n", tr->likelihood); */

  do
    {       
      currentLikelihood = tr->likelihood;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	switch(tr->partitionData[model].dataType)	  
	  {
	  case AA_DATA:
	    break;
	  case DNA_DATA:	   
	    optRates(tr, adef, modelEpsilon, model); 	    
	    break;
	  default:
	    assert(0);
	  }
	}           
      
   
      for(model = 0; model < tr->NumberOfModels; model++)
	optAlpha(tr, adef, modelEpsilon, model); 
       
      if(adef->useInvariant)
	{	      
	  /* Need to execute initrav() here because it is not executed by evalInvar */
	  onlyInitrav(tr, tr->start);
	  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    {	      
	      optInvar(tr, adef, modelEpsilon, model); 	      
	    }
	}   

      onlyInitrav(tr, tr->start);
      /* printf("4 %f\n", tr->likelihood); */
      treeEvaluate(tr, 0.25);    
      /* printf("5 %f\n", tr->likelihood); */
#ifdef _DEBUG_AA           
      printf("O %f\n", tr->likelihood);
#endif

      modelEpsilon /= 10.0;
      if(modelEpsilon < LIKELIHOOD_EPSILON)
	modelEpsilon = LIKELIHOOD_EPSILON;           
    }
  while(fabs(currentLikelihood - tr->likelihood) > adef->likelihoodEpsilon);      

}

static void modOptSpecial(tree *tr, analdef *adef)
{
  double currentLikelihood;
  int i; 
  double modelEpsilon = MODEL_EPSILON;
  int model;  

  /* TODO-LOCAL, not sure that this works */

  /* TODO add inits for CAT here just to be safe */

  for(model = 0; model < tr->NumberOfModels; model++)
    {      
      for(i = 0; i <  DNA_RATES; i++)
	tr->initialRates_DNA[model * DNA_RATES + i] = 0.5;
      for(i = 0; i <  AA_RATES; i++)
	tr->initialRates_AA[model * AA_RATES + i] = 0.5;

      if(adef->useInvariant)
	{
	  int lower, upper;
	  int count = 0;
	  int total = 0;

	  /* TODO-MIX simplify */
	  
	  if(tr->NumberOfModels == 1)
	    {
	      lower = 0;
	      upper = tr->cdta->endsite;	      
	    }
	  else
	    {
	      /*lower = tr->modelIndices[model][0];
		upper = tr->modelIndices[model][1];*/
	      lower = tr->partitionData[model].lower;
	      upper = tr->partitionData[model].upper;
	    }
	  
	  for(i = lower; i < upper; i++)
	    {
	      if(tr->invariant[i] < 4) 		
		count += tr->cdta->aliaswgt[i];		  		
	      total += tr->cdta->aliaswgt[i];
	    }
	  tr->invariants[model] = ((double)count)/((double) total);
#ifdef _LOCAL_DATA
	  masterBarrier(THREAD_COPY_INVARIANTS ,tr);
#endif 	
	}   
       
      tr->alphas[model] = 1.0;

      initReversibleGTR(tr, adef, model);
#ifdef _LOCAL_DATA
      masterBarrier(THREAD_COPY_REVERSIBLE, tr);
#endif
      makeGammaCats(model,  tr->alphas, tr->gammaRates); 
#ifdef _LOCAL_DATA
      masterBarrier(THREAD_COPY_GAMMA_RATES, tr);
#endif
    }

  resetBranches(tr);
  tr->start = tr->nodep[1];

  /* no need for individual models here, just an init on params equal for all partitions*/

  evaluateGenericInitrav(tr, tr->start);
 

  /*printf("INIT %1.20f\n", tr->likelihood);*/

  treeEvaluate(tr, 1.0);

  /*  printf("TREE EVAL %1.20f\n", tr->likelihood); */
 
  do
    {
      currentLikelihood = tr->likelihood;
                  
      if(adef->model == M_GTRGAMMA)
	{	       	           
	  for(model = 0; model < tr->NumberOfModels; model++)
	    optRates(tr, adef, modelEpsilon, model);   			
	
	  for(model = 0; model < tr->NumberOfModels; model++)
	    optAlpha(tr, adef, modelEpsilon, model);       	 	

	  if(adef->useInvariant)
	    {	      
	      /* Need to execute initrav() here because it is not executed by evalInvar */

	      onlyInitrav(tr, tr->start);

	      for(model = 0; model < tr->NumberOfModels; model++)
		{		  
		  optInvar(tr, adef, modelEpsilon, model); 	      
		}
	    }
	}
      else
	{
	  if(adef->model == M_PROTGAMMA)
	    {
	      if(adef->proteinMatrix == GTR)
		{
		  for(model = 0; model < tr->NumberOfModels; model++)
		    optRates(tr, adef, modelEpsilon, model); 
		}	      
	      	      
	      for(model = 0; model < tr->NumberOfModels; model++)				 
		{		  
		  optAlpha(tr, adef, modelEpsilon, model);    	      		 
		  /*printf("ALPHA[%d] %f", model, tr->likelihood);*/
		}
	      
	      if(adef->useInvariant)
		{		 	
		  /* Need to execute initrav() here because it is not executed by evalInvar */

		  onlyInitrav(tr, tr->start);

		  for(model = 0; model < tr->NumberOfModels; model++)	      
		    {
		      optInvar(tr, adef, modelEpsilon, model); 	     		    		    
		    }
		}
	    }
	  else
	    {	     
	      assert(adef->model == M_GTRCAT || adef->model == M_PROTCAT);		  		  
	      
	      if((adef->model == M_PROTCAT && adef->proteinMatrix == GTR) || adef->model == M_GTRCAT)
		{		     
		  for(model = 0; model < tr->NumberOfModels; model++)
		    optRates(tr, adef, modelEpsilon, model); 		    
		}
	      
	      optimizeAllRateCategories(tr);		 		   	   
	    } 
	}
      

      onlyInitrav(tr, tr->start);
      
      treeEvaluate(tr, 0.25);      

      /*printf("ModOPT: %1.20f\n", tr->likelihood);*/

      modelEpsilon /= 10.0;
      if(modelEpsilon < LIKELIHOOD_EPSILON)
	modelEpsilon = LIKELIHOOD_EPSILON;           
    }
  while(fabs(currentLikelihood - tr->likelihood) > adef->likelihoodEpsilon);      
}


void modOptModel(tree *tr, analdef *adef, int model)
{ 
  double currentLikelihood, globalLikelihood;
  int i; 
  double modelEpsilon = MODEL_EPSILON;
  
  /* TODO-MIX, not sure that this works here  */

  resetBranches(tr);        
  
  for(i = 0; i <  DNA_RATES; i++)
    tr->initialRates_DNA[model * DNA_RATES + i] = 0.5;
  for(i = 0; i <  AA_RATES; i++)
    tr->initialRates_AA[model * AA_RATES + i] = 0.5;
        
  if(adef->useInvariant)
    {
      int lower, upper;
      int count = 0;
      int total = 0;	  	       

      lower = tr->partitionData[model].lower;
      upper = tr->partitionData[model].upper;

	  
      for(i = lower; i < upper; i++)
	{
	  if(tr->invariant[i] < 4) 		
	    count += tr->cdta->aliaswgt[i];		  		
	  total += tr->cdta->aliaswgt[i];
	}
      tr->invariants[model] = ((double)count)/((double) total);	
#ifdef _LOCAL_DATA
      masterBarrier(THREAD_COPY_INVARIANTS ,tr);
#endif 
    }   	
      
  tr->alphas[model] = 1.0;

  initReversibleGTR(tr, adef, model);
#ifdef _LOCAL_DATA
  masterBarrier(THREAD_COPY_REVERSIBLE, tr);
#endif
  makeGammaCats(model,  tr->alphas, tr->gammaRates);  
#ifdef _LOCAL_DATA
  masterBarrier(THREAD_COPY_GAMMA_RATES, tr);
#endif
  
  tr->start = tr->nodep[1];  


  evaluateGenericInitravPartition(tr, tr->start, model);


  treeEvaluatePartition(tr, 1.0, model);
  /*printf("Partition %d likelihood %f\n", model, tr->likelihood);      */
  globalLikelihood = tr->likelihood;  
  
  do
    {	  
      currentLikelihood = globalLikelihood;	
      globalLikelihood = 0;

      evaluateGenericInitravPartition(tr, tr->start, model);  
	  	  
      if(adef->model == M_GTRGAMMA)
	{	       	           	     
	  optRates(tr, adef, modelEpsilon, model);   				      	     
	  optAlpha(tr, adef, modelEpsilon, model);       	 	
	      
	  if(adef->useInvariant)
	    {	      		  		  
	      initravPartition(tr, tr->start, model);
	      initravPartition(tr, tr->start->back, model);		  	       
	      optInvar(tr, adef, modelEpsilon, model); 	   
	    }
	}	    
      else
	{
	  if(adef->model == M_PROTGAMMA)
	    {
	      if(adef->proteinMatrix == GTR)		    
		optRates(tr, adef, modelEpsilon, model); 		  
		  
	      optAlpha(tr, adef, modelEpsilon, model);    	      		 
		  
	      if(adef->useInvariant)
		{		 			     		      
		  initravPartition(tr, tr->start, model);
		  initravPartition(tr, tr->start->back, model);		      			
		  optInvar(tr, adef, modelEpsilon, model); 	
		}
	    }		
	  else
	    {		  
	      /* don't remove this, it's for testing ! 
		 
	      if(adef->model == M_GTRCAT)
	      {
	      for(model = 0; model < tr->NumberOfModels; model++)
	      optRates(tr, adef, modelEpsilon, model);		   
	      }
	      */
	      
	      assert(0);	     	  
	    } 
	}
	  
      initravPartition(tr, tr->start, model);
      initravPartition(tr, tr->start->back, model);                        	
      treeEvaluatePartition(tr, 0.25, model);          
      globalLikelihood = tr->likelihood;
    
      /*printf("AFTER ALL OPT %f -> %f\n", currentLikelihood, globalLikelihood);*/
      
      modelEpsilon /= 10.0;
      if(modelEpsilon < LIKELIHOOD_EPSILON)
	modelEpsilon = LIKELIHOOD_EPSILON;           
    }
  while(fabs(currentLikelihood - globalLikelihood) > adef->likelihoodEpsilon);     
  

}





/*********************FUNCTIONS FOOR EXACT MODEL OPTIMIZATION UNDER GTRGAMMA ***************************************/


/*********************FUNCTIONS FOR APPROXIMATE MODEL OPTIMIZATION ***************************************/

static void optimizeAlphaMULT(tree *tr, int model, analdef *adef)
{
  int k; 
  double currentTT, startLikelihood, maxLikelihoodMinus,
    maxLikelihoodPlus, maxTTPlus, maxTTMinus, spacing,
    tree_likelihood;
  boolean finish = FALSE;

  spacing = 0.5/((double)optimizeAlphaInvocations);

  currentTT = tr->alphas[model];
  
  /* TODO-MIX this might be an initrav too much in the single model case */

  tree_likelihood = evaluateGenericInitravPartition(tr, tr->start, model);
 
  maxTTPlus = currentTT;
  maxTTMinus = currentTT;

  startLikelihood = tree_likelihood;
  maxLikelihoodMinus = tree_likelihood;
  maxLikelihoodPlus = tree_likelihood;
      
  k = 1;       

  while(!finish && ((currentTT - spacing * k) > ALPHA_MIN)) 
    {                           
      tree_likelihood = evaluateAlpha(tr, currentTT - spacing * k, model);

      if(tree_likelihood > maxLikelihoodMinus)
	{
	  finish = (fabs(tree_likelihood - maxLikelihoodMinus) < adef->likelihoodEpsilon);

	  maxLikelihoodMinus = tree_likelihood;
	  maxTTMinus = currentTT - spacing * k;
	}    
      else
	finish = TRUE;

      k++;
    }

  finish = FALSE;
  k = 1;
  tree_likelihood = startLikelihood;

  while(!finish && ((currentTT + spacing * k) < ALPHA_MAX)) 
    {            
      tree_likelihood = evaluateAlpha(tr, currentTT + spacing * k, model);

      if(tree_likelihood > maxLikelihoodPlus)
	{
	  finish = (fabs(tree_likelihood - maxLikelihoodPlus) < adef->likelihoodEpsilon);

	  maxLikelihoodPlus = tree_likelihood;
	  maxTTPlus = currentTT + spacing * k;
	}	 
      else
	finish = TRUE;

      k++;
    }  

  if(maxLikelihoodPlus > startLikelihood || maxLikelihoodMinus > startLikelihood)
    {
      if(maxLikelihoodPlus > maxLikelihoodMinus)	
	tr->alphas[model] = maxTTPlus;	
      else	
	tr->alphas[model] = maxTTMinus;	    
      
      makeGammaCats(model, tr->alphas, tr->gammaRates); 
#ifdef _LOCAL_DATA
      masterBarrier(THREAD_COPY_GAMMA_RATES, tr);
#endif  
    }
  else
    { 
      tr->alphas[model] = currentTT;
      
      makeGammaCats(model, tr->alphas, tr->gammaRates);
#ifdef _LOCAL_DATA
      masterBarrier(THREAD_COPY_GAMMA_RATES, tr);
#endif
    }    
}



static double alterRatesMULT(tree *tr, int k, analdef *adef, int model)
{
  int i;
  double granularity = 0.1/((double)optimizeRatesInvocations);
  double bestLikelihood, maxLikelihoodMinus, maxLikelihoodPlus,
    treeLikelihood, originalRate, maxRateMinus, maxRatePlus;
  boolean finish = FALSE;
  
  treeLikelihood =  evaluateGenericInitravPartition(tr, tr->start, model);  

  bestLikelihood = treeLikelihood;
  maxLikelihoodMinus = treeLikelihood;
  maxLikelihoodPlus = treeLikelihood;
  
  i = 1;

  originalRate = getRateModel(tr, model, k);
  
  maxRateMinus = maxRatePlus = originalRate;

  while(!finish && ((originalRate - granularity * i) > RATE_MIN))
    {      
      treeLikelihood = evaluateRate(tr, k, originalRate - granularity * i, adef, model);

      if(treeLikelihood > maxLikelihoodMinus)
	{
	  finish = (fabs(treeLikelihood - maxLikelihoodMinus) < adef->likelihoodEpsilon);
	  
	  maxLikelihoodMinus = treeLikelihood;
	  maxRateMinus = originalRate - granularity * i;
	}      
      else
	finish = TRUE;

      i++;
    }

  i = 1;
  treeLikelihood = bestLikelihood;

  setRateModel(tr, model, originalRate, k); 
  finish = FALSE;

  while(!finish && ((originalRate + i * granularity) < RATE_MAX))
    {     
      treeLikelihood = evaluateRate(tr, k, originalRate + granularity * i, adef, model);

      if(treeLikelihood > maxLikelihoodPlus)
	{
	  finish = (fabs(treeLikelihood - maxLikelihoodPlus) < adef->likelihoodEpsilon);
	  
	  maxLikelihoodPlus = treeLikelihood;
	  maxRatePlus = originalRate + granularity * i;
	}      
      else
	finish = TRUE;

      i++;
    }

  if(maxLikelihoodPlus > bestLikelihood || maxLikelihoodMinus > bestLikelihood)
    {
      if(maxLikelihoodPlus > maxLikelihoodMinus)
	{	  
	  setRateModel(tr, model, maxRatePlus, k);	 
	  initReversibleGTR(tr, adef, model);
#ifdef _LOCAL_DATA
	  masterBarrier(THREAD_COPY_REVERSIBLE, tr);
#endif
	  return maxLikelihoodPlus;
	}
      else
	{
	  setRateModel(tr, model, maxRateMinus, k);	 
	  initReversibleGTR(tr, adef, model);
#ifdef _LOCAL_DATA
	  masterBarrier(THREAD_COPY_REVERSIBLE, tr);
#endif
	  return maxLikelihoodMinus;
	}              
    }
  else
    {      
      setRateModel(tr, model, originalRate, k);     
      initReversibleGTR(tr, adef, model);
#ifdef _LOCAL_DATA
      masterBarrier(THREAD_COPY_REVERSIBLE, tr);
#endif
      return bestLikelihood;
    }
}

static void optimizeRates(tree *tr, analdef *adef)
{
  int i, model, numberOfRates;

  for(model = 0; model < tr->NumberOfModels; model++)         	
    {       
      switch(tr->partitionData[model].dataType)
	{
	case AA_DATA:
	  /* 
	     numberOfRates = AA_RATES; 
	  */
	  
	  /* no rate opt for AA TODO-MIX remove GTR for prots!*/

	  /* 
	     for(i = 0; i < numberOfRates; i++)		      	    	   
	     alterRatesMULT(tr, i, adef, model);   
	  */
	  break;
	case DNA_DATA:
	  numberOfRates = DNA_RATES;
	  for(i = 0; i < numberOfRates; i++)		      	    	   
	    alterRatesMULT(tr, i, adef, model);   
	  break;
	default:
	  assert(0);
	}           
    }

  evaluateGenericInitrav(tr, tr->start);
         
  optimizeRatesInvocations++;
}





static int catCompare(const void *p1, const void *p2)
{
 rateCategorize *rc1 = (rateCategorize *)p1;
 rateCategorize *rc2 = (rateCategorize *)p2;

  double i = rc1->accumulatedSiteLikelihood;
  double j = rc2->accumulatedSiteLikelihood;
  
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}


static void categorize(tree *tr, rateCategorize *rc)
{
  int i, k, found;
  double temp, diff, min;

  for (i = 0; i < tr->cdta->endsite; i++) 
      {
	temp = tr->cdta->patrat[i];
	found = 0;
	for(k = 0; k < tr->NumberOfCategories; k++)
	  {
	    if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
	      {
		found = 1;
		tr->cdta->rateCategory[i] = k;				
		break;
	      }
	  }
	if(!found)
	  {
	    min = fabs(temp - rc[0].rate);
	    tr->cdta->rateCategory[i] = 0;

	    for(k = 1; k < tr->NumberOfCategories; k++)
	    {
	      diff = fabs(temp - rc[k].rate);
	      if(diff < min)
		{
		  min = diff;
		  tr->cdta->rateCategory[i] = k;
		}
	    }
	  }
      }

  for(k = 0; k < tr->NumberOfCategories; k++)
    tr->cdta->patrat[k] = rc[k].rate; 

}







void optRateCat(tree *tr, int i, double lower_spacing, double upper_spacing, double *lhs)
{
  double initialRate, initialLikelihood, 
    leftLH, rightLH, leftRate, rightRate, v;
  const double epsilon = 0.00001;
  int k;

  tr->cdta->patrat[i] = tr->cdta->patratStored[i];     
  initialRate = tr->cdta->patrat[i];
      
  initialLikelihood = evaluatePartialGeneric(tr, i, initialRate); 
    

  leftLH = rightLH = initialLikelihood;
  leftRate = rightRate = initialRate;
  
  k = 1;
  
  while((initialRate - k * lower_spacing > 0.0001) && 
	((v = evaluatePartialGeneric(tr, i, initialRate - k * lower_spacing)) 
	 > leftLH) && 
	(fabs(leftLH - v) > epsilon))  
    {	  
#ifndef WIN32
      if(isnan(v))
	assert(0);
#endif
      
      leftLH = v;
      leftRate = initialRate - k * lower_spacing;
      k++;	  
    }      
  
  k = 1;
  
  while(((v = evaluatePartialGeneric(tr, i, initialRate + k * upper_spacing)) > rightLH) &&
	(fabs(rightLH - v) > epsilon))    	
    {
#ifndef WIN32
      if(isnan(v))
	assert(0);
#endif     
      rightLH = v;
      rightRate = initialRate + k * upper_spacing;	 
      k++;
    }           
  
  if(rightLH > initialLikelihood || leftLH > initialLikelihood)
    {
      if(rightLH > leftLH)	    
	{	     
	  tr->cdta->patrat[i] = rightRate;
	  lhs[i] = rightLH;
	}
      else
	{	      
	  tr->cdta->patrat[i] = leftRate;
	  lhs[i] = leftLH;
	}
    }
  else
    lhs[i] = initialLikelihood;
  
  tr->cdta->patratStored[i] = tr->cdta->patrat[i];     
}

#ifdef _LOCAL_DATA

void optRateCat_LOCAL(tree *localTree, int lower, int upper, double lower_spacing, double upper_spacing, double *lhs)
{
  double initialRate, initialLikelihood, 
    leftLH, rightLH, leftRate, rightRate, v;
  const double epsilon = 0.00001;
  int k, i; 

  /* double sum = 0.0; */

  for(i = lower; i < upper; i++)
    {
      localTree->strided_patrat[i] = localTree->strided_patratStored[i];     
      initialRate = localTree->strided_patrat[i];
      
      initialLikelihood = evaluatePartialGeneric(localTree, i, initialRate); 
            
      leftLH = rightLH = initialLikelihood;
      leftRate = rightRate = initialRate;
      
      k = 1;
      
      while((initialRate - k * lower_spacing > 0.0001) && 
	    ((v = evaluatePartialGeneric(localTree, i, initialRate - k * lower_spacing)) 
	     > leftLH) && 
	    (fabs(leftLH - v) > epsilon))  
	{	  
#ifndef WIN32
	  if(isnan(v))
	    assert(0);
#endif
	  
	  leftLH = v;
	  leftRate = initialRate - k * lower_spacing;
	  k++;	  
	}      
      
      k = 1;
      
      while(((v = evaluatePartialGeneric(localTree, i, initialRate + k * upper_spacing)) > rightLH) &&
	    (fabs(rightLH - v) > epsilon))    	
	{
#ifndef WIN32
	  if(isnan(v))
	    assert(0);
#endif     
	  rightLH = v;
	  rightRate = initialRate + k * upper_spacing;	 
	  k++;
	}           
      
      if(rightLH > initialLikelihood || leftLH > initialLikelihood)
	{
	  if(rightLH > leftLH)	    
	    {	     
	      localTree->strided_patrat[i] = rightRate;
	      lhs[i] = rightLH;
	    }
	  else
	    {	      
	      localTree->strided_patrat[i] = leftRate;
	      lhs[i] = leftLH;
	    }
	}
      else
	lhs[i] = initialLikelihood;            

      /* sum += lhs[i]; */

      localTree->strided_patratStored[i] = localTree->strided_patrat[i];     
    }

  /* printf("%f \n", sum); */
}



#endif


void optimizeRateCategories(tree *tr, int _maxCategories)
{
  int i, k;
  double  temp, wtemp;   
  double lower_spacing, upper_spacing;
  int maxCategories = _maxCategories;
  double initialLH = tr->likelihood;
  double *oldRat =    (double *)malloc(sizeof(double) * tr->cdta->endsite);
  double *ratStored = (double *)malloc(sizeof(double) * tr->cdta->endsite);
  double *oldwr =     (double *)malloc(sizeof(double) * tr->cdta->endsite);
  double *oldwr2 =    (double *)malloc(sizeof(double) * tr->cdta->endsite);
  double *lhs =       (double *)malloc(sizeof(double) * tr->cdta->endsite);
  int *oldCategory =  (int *)malloc(sizeof(int) * tr->cdta->endsite);  
  int oldNumber;   
  
  assert(isTip(tr->start->number, tr->rdta->numsp));   
  determineFullTraversal(tr->start, tr);

  if(optimizeRateCategoryInvocations == 1)
    {
      lower_spacing = 0.5 / ((double)optimizeRateCategoryInvocations);
      upper_spacing = 1.0 / ((double)optimizeRateCategoryInvocations);
    }
  else
    {
      lower_spacing = 0.05 / ((double)optimizeRateCategoryInvocations);
      upper_spacing = 0.1 / ((double)optimizeRateCategoryInvocations);
    }

  if(lower_spacing < 0.001)
    lower_spacing = 0.001;

  if(upper_spacing < 0.001)
    upper_spacing = 0.001;

  optimizeRateCategoryInvocations++;
  
  oldNumber = tr->NumberOfCategories;

  for(i = 0; i < tr->cdta->endsite; i++)
    {    
      oldCategory[i] = tr->cdta->rateCategory[i];
      ratStored[i] = tr->cdta->patratStored[i];    
      oldRat[i] = tr->cdta->patrat[i];
      oldwr[i] =  tr->cdta->wr[i];
      oldwr2[i] =  tr->cdta->wr2[i];
    }

#ifdef _USE_PTHREADS
  tr->lhs = lhs;
  tr->lower_spacing = lower_spacing;
  tr->upper_spacing = upper_spacing; 
  masterBarrier(THREAD_RATE_CATS, tr);
#else
  for(i = 0; i < tr->cdta->endsite; i++)      
    optRateCat(tr, i, lower_spacing, upper_spacing, lhs);
#endif

  /* {
    double sum = 0.0;

    for (i = 0; i < tr->cdta->endsite; i++) 
      {
	printf("%d %f %f %f\n", i, lhs[i],  tr->cdta->patrat[i], tr->cdta->patratStored[i]);
	sum += lhs[i];
      }
    printf("LH %f\n", sum);
    }*/
  
       
  {     
    rateCategorize *rc = (rateCategorize *)malloc(sizeof(rateCategorize) * tr->cdta->endsite);
    int where;
    int found = 0;
    for (i = 0; i < tr->cdta->endsite; i++)
      {
	rc[i].accumulatedSiteLikelihood = 0;
	rc[i].rate = 0;
      }
      
    where = 1;   
    rc[0].accumulatedSiteLikelihood = lhs[0];
    rc[0].rate = tr->cdta->patrat[0];
    tr->cdta->rateCategory[0] = 0;
    
    for (i = 1; i < tr->cdta->endsite; i++) 
      {
	temp = tr->cdta->patrat[i];
	found = 0;
	for(k = 0; k < where; k++)
	  {
	    if(temp == rc[k].rate || (fabs(temp - rc[k].rate) < 0.001))
	      {
		found = 1;						
		rc[k].accumulatedSiteLikelihood += lhs[i];	
		break;
	      }
	  }
	if(!found)
	  {	    
	    rc[where].rate = temp;	    
	    rc[where].accumulatedSiteLikelihood += lhs[i];	    
	    where++;
	  }
	}

    qsort(rc, where, sizeof(rateCategorize), catCompare);
  
    if(where < maxCategories)
      {
	tr->NumberOfCategories = where;
	categorize(tr, rc);
      }
    else
      {
	tr->NumberOfCategories = maxCategories;	
	categorize(tr, rc);
      }
      
    free(rc);
  
    for (i = 0; i < tr->cdta->endsite; i++) 
      {	
	temp = tr->cdta->patrat[tr->cdta->rateCategory[i]];

	tr->cdta->wr[i]  = wtemp = temp * tr->cdta->aliaswgt[i];
	tr->cdta->wr2[i] = temp * wtemp;
      }     


   

#ifdef _LOCAL_DATA
    /* TODO this could actually be merged with evaluateGenericInitrav */
    /* will help to reduce sync points */
    /* same holds for many other model opt procedures of type change model 
       param and then execute initrav */
    
    masterBarrier(THREAD_COPY_RATE_CATS, tr);
#endif


    evaluateGenericInitrav(tr, tr->start);

    /* printf("%f\n", tr->likelihood); */

    if(tr->likelihood < initialLH)
      {	
	tr->NumberOfCategories = oldNumber;
	for (i = 0; i < tr->cdta->endsite; i++)
	    {
	      tr->cdta->patratStored[i] = ratStored[i]; 
	      tr->cdta->rateCategory[i] = oldCategory[i];
	      tr->cdta->patrat[i] = oldRat[i];	    
	      tr->cdta->wr[i]  = oldwr[i];
	      tr->cdta->wr2[i] = oldwr2[i];
	    }       

#ifdef _LOCAL_DATA
	masterBarrier(THREAD_COPY_RATE_CATS, tr);
#endif


	evaluateGenericInitrav(tr, tr->start);
      }
    }
  
  free(oldCategory);
  free(oldRat);
  free(ratStored);
  free(oldwr);
  free(oldwr2); 
  free(lhs); 
}
  




static void optimizeAlphas(tree *tr, analdef *adef)
{
  int model;   

  for(model = 0; model < tr->NumberOfModels; model++)       
    optimizeAlphaMULT(tr, model, adef);	  
	    

  evaluateGenericInitrav(tr, tr->start);    

  optimizeAlphaInvocations++;
}



static void optimizeInvariantMULT(tree *tr, int model, analdef *adef)
{
  int k; 
  double currentTT, startLikelihood, maxLikelihoodMinus,
    maxLikelihoodPlus, maxTTPlus, maxTTMinus, spacing,
    tree_likelihood;
  boolean finish = FALSE;

  spacing = 0.1/((double)optimizeInvarInvocations);

  currentTT = tr->invariants[model];
  
  tree_likelihood = evaluatePartitionGeneric(tr, tr->start, model);
 
  maxTTPlus = currentTT;
  maxTTMinus = currentTT;

  startLikelihood = tree_likelihood;
  maxLikelihoodMinus = tree_likelihood;
  maxLikelihoodPlus = tree_likelihood;
      
  k = 1;       

  while(!finish && ((currentTT - spacing * k) > INVAR_MIN)) 
    {                            
      tree_likelihood =  evaluateInvar(tr, currentTT - spacing * k, model);

      if(tree_likelihood > maxLikelihoodMinus)
	{
	  finish = (fabs(tree_likelihood - maxLikelihoodMinus) < adef->likelihoodEpsilon);

	  maxLikelihoodMinus = tree_likelihood;
	  maxTTMinus = currentTT - spacing * k;
	}    
      else
	finish = TRUE;

      k++;
    }

  finish = FALSE;
  k = 1;
  tree_likelihood = startLikelihood;

  while(!finish && ((currentTT + spacing * k) < INVAR_MAX)) 
    {           
      tree_likelihood =  evaluateInvar(tr, currentTT + spacing * k, model);

      if(tree_likelihood > maxLikelihoodPlus)
	{
	  finish = (fabs(tree_likelihood - maxLikelihoodPlus) < adef->likelihoodEpsilon);

	  maxLikelihoodPlus = tree_likelihood;
	  maxTTPlus = currentTT + spacing * k;
	}	 
      else
	finish = TRUE;

      k++;
    }  

  if(maxLikelihoodPlus > startLikelihood || maxLikelihoodMinus > startLikelihood)
    {
      if(maxLikelihoodPlus > maxLikelihoodMinus)	
	tr->invariants[model] = maxTTPlus;	
      else	
	tr->invariants[model] = maxTTMinus;
#ifdef _LOCAL_DATA
      masterBarrier(THREAD_COPY_INVARIANTS ,tr);
#endif 
    }
  else
    { 
      tr->invariants[model] = currentTT;
#ifdef _LOCAL_DATA
      masterBarrier(THREAD_COPY_INVARIANTS ,tr);
#endif  
    }    
}




static void optimizeInvariants(tree *tr, analdef *adef)
{
  int model;
  
  /* Don't need full tree traversal for invar ! just do it once at this point*/

  onlyInitrav(tr, tr->start);
     
  for(model = 0; model < tr->NumberOfModels; model++) 
    {
      optimizeInvariantMULT(tr, model, adef);
    }

  evaluateGeneric(tr, tr->start);        

  optimizeInvarInvocations++;
}




/*static void rateCategorizeExpectation(tree *tr, double *vector, int n)
{
  int i, j;
  double e;
  double allSum;

  for(j = 0; j < 4; j++)
    printf("Rate %d %f\n", j, tr->gammaRates[j]);
  printf("\n");

  for(i = 0; i < n; i++)
    {
      allSum = 0.25 * (vector[i * 4] + vector[i * 4 + 1] + vector[i * 4 + 2] + vector[i * 4 + 3]);

      e = 0.0;
      for(j = 0; j < 4; j++)	
	e += (tr->gammaRates[j] * ((0.25 * vector[i * 4 + j]) / allSum));
      printf("Site %d E %f\n", i, e);
    }
    }*/



void quickAndDirtyOptimization(tree *tr, analdef *adef)
{
  int oldInv;
  double initialLH;     
   
  /* printf("ENTER %1.80f\n", tr->likelihood); */

  /* 
     TODO-MIX for AA data only there is one useless call to initrav from within 
     optimizeRates 
  */
  
  oldInv = optimizeRatesInvocations;    

  if(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
    {
      do
	{	    
	  initialLH = tr->likelihood;
	  optimizeRates(tr, adef);       
	}
      while((fabs(tr->likelihood - initialLH) > adef->likelihoodEpsilon) && optimizeRatesInvocations < oldInv + 10); 
    }
  else
    {
      /* TODO-MIX this is extremely ugly, for the time being just for backward compatibility */     

      do
	{	    
	  initialLH = tr->likelihood;
	  optimizeRates(tr, adef);       
	}
      while(tr->likelihood > initialLH && optimizeRatesInvocations < oldInv + 10);
    }

  /* printf("ONE %f \n", tr->likelihood); */


  if(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
    {           
      oldInv = optimizeAlphaInvocations;	
      do
	{    	   
	  initialLH = tr->likelihood;	     
	  optimizeAlphas(tr, adef);		   	  
	}
      while((fabs(tr->likelihood - initialLH) > adef->likelihoodEpsilon) && optimizeAlphaInvocations < oldInv + 10);
    }

  /* printf("%f \n", tr->likelihood); */

  if(tr->rateHetModel == GAMMA_I)
    {            
      oldInv = optimizeInvarInvocations;
      
      do
	{    	   
	  initialLH = tr->likelihood;
	  optimizeInvariants(tr, adef);		   	  
	}
      while((fabs(tr->likelihood - initialLH) > adef->likelihoodEpsilon) && optimizeInvarInvocations < oldInv + 10);
    }   
  
  /* printf("%f \n", tr->likelihood); */

  if(tr->rateHetModel == CAT)   
    optimizeRateCategories(tr, adef->categories);             

  /* printf("EXIT %1.80f\n", tr->likelihood);*/
}

int optimizeModel (tree *tr, analdef *adef)
{
  double startLH = tr->likelihood;

  if(!adef->categorizeGamma)
    {
      quickAndDirtyOptimization(tr, adef);
    }
  else
    {
      if(tr->mixedData)
	{
	  assert(0);
	}
      else
	{
	  switch(adef->model)
	    {          
	    case M_PROTCAT: 
	    case M_GTRCAT:
	      {
		int j;
		double catLikelihood;
		
		int    *catVector        = (int *)malloc(sizeof(int) * tr->cdta->endsite);	 
		double *wr               = (double *)malloc(sizeof(double) * tr->cdta->endsite);
		double *wr2              = (double *)malloc(sizeof(double) * tr->cdta->endsite);
		double *patrat           = (double *)malloc(sizeof(double) * tr->cdta->endsite);
		double *initialRates_DNA = (double *)malloc(tr->NumberOfModels * DNA_RATES * sizeof(double));
		double *initialRates_AA  = (double *)malloc(tr->NumberOfModels * AA_RATES * sizeof(double));
		
		memcpy(catVector,    tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);	  
		memcpy(initialRates_DNA, tr->initialRates_DNA, tr->NumberOfModels * DNA_RATES * sizeof(double));
		memcpy(initialRates_AA,  tr->initialRates_AA,  tr->NumberOfModels * AA_RATES * sizeof(double));
		
		memcpy(wr, tr->cdta->wr, sizeof(double) * tr->cdta->endsite);
		memcpy(wr2, tr->cdta->wr2, sizeof(double) * tr->cdta->endsite);
		memcpy(patrat, tr->cdta->patrat, sizeof(double) * tr->cdta->endsite);
		
		catLikelihood = tr->likelihood;				    

		catToGamma(tr, adef);

		
		for(j = 0; j < tr->cdta->endsite; j++)
		  tr->cdta->wr[j] = tr->cdta->aliaswgt[j];
		

		evaluateGenericInitrav(tr, tr->start);	 
		
		/*printf("GAMMA-ENTRY %f\n", tr->likelihood);*/
		
		quickAndDirtyOptimization(tr, adef);
		
		/*printf("GAMMA exit %f\n", tr->likelihood);*/
		
		categorizeGeneric(tr, tr->start);
		
		
		gammaToCat(tr, adef);
			
		evaluateGenericInitrav(tr, tr->start);	
	     
		/*printf("CAT exit %f\n", tr->likelihood);*/
		
		if(catLikelihood > tr->likelihood)
		  {	  
		    int model;
		    /*printf("Need roll-back Rates + CATs %f %f !\n", catLikelihood, tr->likelihood);*/
		    memcpy(tr->cdta->rateCategory, catVector, sizeof(int) * tr->cdta->endsite);
		    
		    memcpy(tr->initialRates_DNA, initialRates_DNA, tr->NumberOfModels * DNA_RATES * sizeof(double));
		    memcpy(tr->initialRates_AA, initialRates_AA, tr->NumberOfModels * AA_RATES * sizeof(double));
		    
		    memcpy(tr->cdta->wr, wr, sizeof(double) * tr->cdta->endsite);
		    memcpy(tr->cdta->wr2,wr2, sizeof(double) * tr->cdta->endsite);
		    memcpy(tr->cdta->patrat, patrat, sizeof(double) * tr->cdta->endsite);
		    
		    for(model = 0;  model < tr->NumberOfModels; model++)		      
		      initReversibleGTR(tr, adef, model);		     
#ifdef _LOCAL_DATA
		    masterBarrier(THREAD_COPY_REVERSIBLE, tr);
#endif
		    

		    evaluateGenericInitrav(tr, tr->start);     
	      
		    /*printf("Roll Back %f\n", tr->likelihood);*/
		  }

		free(catVector);
		free(initialRates_DNA);
		free(initialRates_AA);
		free(wr);
		free(wr2);
		free(patrat);	 
	      }      
	      break;         	   
	    default:
	      assert(0);
	    }
	}
    }

  if(optimizeRatesInvocations > 90)
    optimizeRatesInvocations = 90;  
  if(optimizeRateCategoryInvocations > 90)
    optimizeRateCategoryInvocations = 90;
  if(optimizeAlphaInvocations > 90)
    optimizeAlphaInvocations = 90;
  if(optimizeInvarInvocations > 90)
    optimizeInvarInvocations = 90;

  if(startLH > tr->likelihood) return 0;
  else return 1;
}




void optimizeAllRateCategories(tree *tr)
{
  int i;
  double temp, wtemp;   
  double lower_spacing, upper_spacing; 
  double initialLH = tr->likelihood;
  double *oldRat = (double *)malloc(sizeof(double) * tr->cdta->endsite);
  double *ratStored = (double *)malloc(sizeof(double) * tr->cdta->endsite);
  double *oldwr = (double *)malloc(sizeof(double) * tr->cdta->endsite);
  double *oldwr2 = (double *)malloc(sizeof(double) * tr->cdta->endsite);
  double *lhs = (double *)malloc(sizeof(double) * tr->cdta->endsite);
  int *oldCategory = (int *)malloc(sizeof(int) * tr->cdta->endsite);  
  int oldNumber;  
  
  assert(isTip(tr->start->number, tr->rdta->numsp));     
  determineFullTraversal(tr->start, tr); 

  if(optimizeRateCategoryInvocations == 1)
    {
      lower_spacing = 0.5 / ((double)optimizeRateCategoryInvocations);
      upper_spacing = 1.0 / ((double)optimizeRateCategoryInvocations);
    }
  else
    {
      lower_spacing = 0.05 / ((double)optimizeRateCategoryInvocations);
      upper_spacing = 0.1 / ((double)optimizeRateCategoryInvocations);
    }

  if(lower_spacing < 0.001)
    lower_spacing = 0.001;

  if(upper_spacing < 0.001)
    upper_spacing = 0.001;

  optimizeRateCategoryInvocations++;
  
  oldNumber = tr->NumberOfCategories;

  for(i = 0; i < tr->cdta->endsite; i++)
    {    
      oldCategory[i] = tr->cdta->rateCategory[i];
      ratStored[i] = tr->cdta->patratStored[i];    
      oldRat[i] = tr->cdta->patrat[i];
      oldwr[i] =  tr->cdta->wr[i];
      oldwr2[i] =  tr->cdta->wr2[i];
    }

#ifdef _USE_PTHREADS
  tr->lhs = lhs;
  tr->lower_spacing = lower_spacing;
  tr->upper_spacing = upper_spacing;  
  masterBarrier(THREAD_RATE_CATS, tr);
#else
  for(i = 0; i < tr->cdta->endsite; i++)      
    optRateCat(tr, i, lower_spacing, upper_spacing, lhs);
#endif

     
  tr->NumberOfCategories = tr->cdta->endsite;

  for (i = 0; i < tr->cdta->endsite; i++) 
    {	
      temp = tr->cdta->patrat[i];
      tr->cdta->rateCategory[i] = i;
      tr->cdta->wr[i]  = wtemp = temp * tr->cdta->aliaswgt[i];
      tr->cdta->wr2[i] = temp * wtemp;
    }     

#ifdef _LOCAL_DATA
  /* TODO this could actually be merged with evaluateGenericInitrav below */
  /* will help to reduce sync points */
  /* same holds for many other model opt procedures of type change model 
       param and then execute initrav */
  masterBarrier(THREAD_COPY_RATE_CATS, tr);
#endif


  evaluateGenericInitrav(tr, tr->start);

  
  if(tr->likelihood < initialLH)
    {	
      tr->NumberOfCategories = oldNumber;
      for (i = 0; i < tr->cdta->endsite; i++)
	{
	  tr->cdta->patratStored[i] = ratStored[i]; 
	  tr->cdta->rateCategory[i] = oldCategory[i];
	  tr->cdta->patrat[i] = oldRat[i];	    
	  tr->cdta->wr[i]  = oldwr[i];
	  tr->cdta->wr2[i] = oldwr2[i];
	  
	}       
      
#ifdef _LOCAL_DATA
      /* TODO this could actually be merged with evaluateGenericInitrav */
      /* will help to reduce sync points */
      /* same holds for many other model opt procedures of type change model 
	 param and then execute initrav */
      masterBarrier(THREAD_COPY_RATE_CATS, tr);
#endif


      evaluateGenericInitrav(tr, tr->start);         
    }
  
  free(oldCategory);
  free(oldRat);
  free(ratStored);
  free(oldwr);
  free(oldwr2); 
  free(lhs); 
}
  

static void cacheZ(tree *tr, double *zs, int model)
{ 
  nodeptr  p;
  int  nodes;
  
  nodes = tr->mxtips  +  3 * (tr->mxtips - 2);
  p = tr->nodep[1];
  while (nodes-- > 0) 
    {
      zs[nodes] = p->z[model]; 
      p++;
    }
} 


static void restoreZ(tree *tr, double *zs, int model)
{ 
  nodeptr  p;
  int  nodes;
  
  nodes = tr->mxtips  +  3 * (tr->mxtips - 2);
  p = tr->nodep[1];
  while (nodes-- > 0) 
    {
      p->z[model]       = zs[nodes]; 
      p->back->z[model] = zs[nodes]; 
      p++;
    }
}



static double treeLengthRec(nodeptr p, tree *tr, int model)
{  
  if(/*p->tip*/ isTip(p->number, tr->rdta->numsp))  
    return 0.0;    
  else
    {
      double x, acc = 0;
      nodeptr q;

      x = p->z[model];
      if (x < zmin) 
	x = zmin;      	 
      x = -log(x) * tr->fracchange;
      q = p->next;
      while(q != p)
	{
	  acc +=  treeLengthRec(q->back, tr, model);
	  q = q->next;
	}

      return acc + x;
    }
}

double treeLength(tree *tr, int model)
{ 
  return (treeLengthRec(tr->start, tr, model) + treeLengthRec(tr->start->back, tr, model));
}


/*static double gappyness(tree *tr, int lower, int upper)
{
  int i, j, count = 0, total = 0;
  double result;

  for(i = 1; i <= tr->mxtips; i++)
    {
    char *y = tr->yVector[i];

      for(j = lower; j < upper; j++)
	{
	  if(y[j] == 15)
	    count++;	  
	  total++;
	}
    }

  result = ((double)count) / ((double)total);

  return result;
  }*/

void optimizeRatesOnly(tree *tr, analdef *adef)
{
  int i = 0;  
  char temporaryFileName[1024] = "";
  FILE *out;
    

  assert(0);
  /* TODO-MIX needs to be checked */

  getStartingTree(tr, adef); 

  strcpy(temporaryFileName, ratesFileName);     
  out = fopen(temporaryFileName, "w");         

  adef->likelihoodEpsilon = 0.5;

  if(adef->model == M_GTRGAMMA || adef->model == M_PROTGAMMA)
    {   
      modOpt(tr, adef);
      /*printf("GAMMA model OPT %f\n", tr->likelihood);*/
          
      categorizeGeneric(tr, tr->start);           

      for(i = 0; i < tr->cdta->endsite; i++)
	fprintf(out, "%d %f\n", i, tr->gammaRates[tr->cdta->rateCategory[i]]);


      if(adef->treeLength)
	{	 
	  FILE *tlf;

	  int lower, upper;
	  int windowSize = 100;
	  int increment  = 20;
	  double *branches = malloc(4 * tr->mxtips * sizeof(double)), tl;
	  
	  
	  /*printf("FileNames %s and %s\n", lengthFileName, lengthFileNameModel);*/
	  
	  for(i = 0; i < tr->cdta->endsite; i++)
	    tr->cdta->wr[i] = tr->cdta->aliaswgt[i];


	  evaluateGenericInitrav(tr, tr->start);
	  
	  /*printf("GAMMA %f\n", tr->likelihood);*/
	  cacheZ(tr, branches, 0);
	 
	  
	  tlf = fopen(lengthFileName, "w"); 
	  
	  for(upper = windowSize, lower = 0; upper < tr->cdta->endsite && upper > lower; 
	      lower += increment, upper = MIN(upper + increment, tr->cdta->endsite))
	    {
	      /*resetBranches(tr);*/
	      restoreZ(tr, branches, 0);
	      /*tr->modelIndices[0][0] = lower;
		tr->modelIndices[0][1] = upper;*/
	      tr->partitionData[0].lower = lower;
	      tr->partitionData[0].upper = upper;
	      tr->start = tr->nodep[1];


	      evaluateGenericInitravPartition(tr, tr->start, 0);   	      

	      /*printf("Winodw %d to %d : %f ->", lower, upper, tr->likelihood);*/
	      treeEvaluatePartition(tr, 1.0, 0);
	      tl = treeLength(tr, 0);
	      /*printf(" %f Length %f\n", tr->likelihood, tl);*/
	      fprintf(tlf, "%d %f\n", lower + (int)((upper - lower) / 2), tl);
	    }      
	  
	  fclose(tlf);
	  
	  modOpt(tr, adef);
	  /*printf("GAMMA model OPT %f\n", tr->likelihood);*/
	  
	  tlf = fopen(lengthFileNameModel, "w"); 
	  
	  for(upper = windowSize, lower = 0; upper < tr->cdta->endsite && upper > lower; 
	      lower += increment, upper = MIN(upper + increment, tr->cdta->endsite))
	    {	  
	      /*tr->modelIndices[0][0] = lower;
		tr->modelIndices[0][1] = upper;*/
	      tr->partitionData[0].lower = lower;
	      tr->partitionData[0].upper = upper;

	      tr->start = tr->nodep[1];


	      evaluateGenericInitravPartition(tr, tr->start, 0);
      
	      /*printf("W-M-OPT %d to %d : %f ->", lower, upper, tr->likelihood);*/
	      modOptModel(tr, adef, 0);
	      tl = treeLength(tr, 0);
	      /*printf(" %f Length %f alpha %f gappyness %f\n", tr->likelihood, tl, 
		tr->alphas[0], gappyness(tr, lower, upper));*/
	      fprintf(tlf, "%d %f\n", lower + (int)((upper - lower) / 2), tl);
	    }      
	  
	  fclose(tlf);
	  
	  free(branches);
	}                 
    }
  else
    {       
      modOptSpecial(tr, adef);       

      for(i = 0; i < tr->cdta->endsite; i++)
	/*fprintf(out, "%d %f\n", i, log(tr->cdta->patrat[i]));*/
	fprintf(out, "%d %f\n", i, tr->cdta->patrat[i]);	     
    }

  fclose(out);
}



/************************* ARNDT_MODE *************************************/

void optimizeArndt(tree *tr, analdef *adef)
{  
  int model;

  modOpt(tr, adef);
    
  printf("\nFinal ML Optimization Likelihood: %f\n", tr->likelihood); 
  printf("\nModel Information:\n\n");
 

  for(model = 0; model < tr->NumberOfModels; model++)		    		    
    {
      double tl;
      char typeOfData[1024];
      
      switch(tr->partitionData[model].dataType)
	{
	case AA_DATA:
	  strcpy(typeOfData,"AA");
	  break;
	case DNA_DATA:
	  strcpy(typeOfData,"DNA");
	  break;
	default:
	  assert(0);
	}     

      printf("Model Parameters of Partition %d, Name: %s, Type of Data: %s\n", 
	     model, tr->partitionData[model].partitionName, typeOfData);
      printf("alpha: %f\n", tr->alphas[model]);
      
      if(adef->useInvariant)      
	printf("invar: %f\n", tr->invariants[model]);    	
                 
      if(adef->perGeneBranchLengths)
	tl = treeLength(tr, model);
      else
	tl = treeLength(tr, 0);
     
      printf("Tree-Length: %f\n", tl);       
      
      switch(tr->partitionData[model].dataType)
	{
	case AA_DATA:
	  break;
	case DNA_DATA:
	  {
	    int k;
	    char *names[6] = {"a<->c", "a<->g", "a<->t", "c<->g", "c<->t", "g<->t"};	 
	    for(k = 0; k < DNA_RATES; k++)			    	     		    
	      printf("rate %s: %f\n", names[k], tr->initialRates_DNA[model * DNA_RATES + k]);	      
		   
	    printf("rate %s: %f\n", names[5], 1.0);
	  }      
	  break;
	default:
	  assert(0);
	}
     
      printf("\n");
    }		    		    
  

  {
    FILE *of;
    char outName[1024];
    int count = 0;
    int billig = (int)(-1.0 * tr->likelihood);
    strcpy(outName, workdir);
    strcat(outName, "RAxML_arndt.");
    strcat(outName,  run_id);
    

    of = fopen(outName, "w");

    if(adef->useInvariant)
      {
	double invar = INVAR_MIN;

	while(invar <= 0.9999)
	  {

	    double alpha =  ALPHA_MIN;
	    tr->invariants[0] = invar;

	    while(alpha <= 5.0)
	      {      	
		double lh, rho, eta;
		
		lh = evaluateAlpha(tr, alpha, 0);
			       
		/* fprintf(of, "%1.10f %1.10f %f\n", invar, alpha, lh);	 */

		if(billig == (int)(-1.0 * lh))
		  printf("%1.10f %1.10f\n", invar, alpha);

		eta = 1.0 / (1.0 + alpha);

		rho = invar + eta - eta * invar;

		fprintf(of, "%1.10f %f\n", rho, lh);

		alpha += 0.01;
		/*if(count % 10000 == 0)
		  printf("%f %f\n", alpha, invar);*/
		count++;
	      }
	    invar += 0.01;
	  }
      }
    else
      {
	double alpha =  ALPHA_MIN;
	
	while(alpha < 5.0)
	  {      	
	    double lh, eta;
	    
	    lh = evaluateAlpha(tr, alpha, 0);

	    eta = 1.0 / (1.0 + alpha);

	    /*fprintf(of, "%1.10f %f\n", alpha, lh);*/
	    fprintf(of, "%1.10f %f\n", eta, lh);

	    if(billig == (int)(-1.0 * lh))
	      printf("0.0 %1.10f\n", alpha);

	    alpha += 0.01;

	    /*if(count % 10000 == 0)
	      printf("%f\n", alpha);*/
	    count++;
	  }	


      }

    printf("Result written to: %s\n", outName);

    fclose(of);
    exit(0);
  }
}

