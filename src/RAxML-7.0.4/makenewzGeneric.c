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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses with 
 *  thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#ifndef WIN32 
#include <unistd.h>
#endif

#ifdef _USE_OMP
extern volatile int             NumberOfThreads;
#include <omp.h>
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "axml.h"

#ifdef _USE_PTHREADS
extern double *reductionBuffer;
extern double *reductionBufferTwo;
extern int NumberOfThreads;
#endif




/*******************/

static void sumCAT(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
		   char *tipX1, char *tipX2, int lower, int n)
{
  int i;
  double *x1, *x2;  
  
  switch(tipCase)
    {
    case TIP_TIP:
      for (i = lower; i < n; i++) 
	{    
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);
	  
	  sum[i * 4]     = x1[0] * x2[0];
	  sum[i * 4 + 1] = x1[1] * x2[1];
	  sum[i * 4 + 2] = x1[2] * x2[2];
	  sum[i * 4 + 3] = x1[3] * x2[3];		    	    	  
	}    
      break;
    case TIP_INNER:         
      for (i = lower; i < n; i++) 
	{    
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];
	  
	  sum[i * 4]     = x1[0] * x2[0];
	  sum[i * 4 + 1] = x1[1] * x2[1];
	  sum[i * 4 + 2] = x1[2] * x2[2];
	  sum[i * 4 + 3] = x1[3] * x2[3];				 
	}   
      break;
    case INNER_INNER:      
      for (i = lower; i < n; i++) 
	{     
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  
	  sum[i * 4]     = x1[0] * x2[0];
	  sum[i * 4 + 1] = x1[1] * x2[1];
	  sum[i * 4 + 2] = x1[2] * x2[2];
	  sum[i * 4 + 3] = x1[3] * x2[3];				   
	}
      break;  
    default:
      assert(0);
    }
}

static void coreGTRCAT(int lower, int upper, int numberOfCategories, double *sum, 
		       double *d1, double *d2, double *wrptr, double *wr2ptr, 
		       double *rptr, double *EIGN, int *cptr, double lz)
{
  int i;
  double 
    *d, *d_start,
    tmp_0, tmp_1, tmp_2, inv_Li, dlnLidlz, d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;
  double e[6];
  double dd1, dd2, dd3;
 
  e[0] = EIGN[0];
  e[1] = EIGN[0] * EIGN[0];
  e[2] = EIGN[1];
  e[3] = EIGN[1] * EIGN[1];
  e[4] = EIGN[2];
  e[5] = EIGN[2] * EIGN[2];
          
  d = d_start = (double *)malloc(numberOfCategories * 3 * sizeof(double));   

  dd1 = e[0] * lz;
  dd2 = e[2] * lz;
  dd3 = e[4] * lz;	  
   
  for(i = 0; i < numberOfCategories; i++)
    {        
      *d++ = exp(dd1 * rptr[i]); 
      *d++ = exp(dd2 * rptr[i]);
      *d++ = exp(dd3 * rptr[i]);	    	    	  	   
    }       

  for (i = lower; i < upper; i++) 	
    {	    	 	   	    
      d = &d_start[3 * cptr[i]];
		
      inv_Li = sum[4 * i];	   
      inv_Li += (tmp_0 = d[0] * sum[4 * i + 1]);
      inv_Li += (tmp_1 = d[1] * sum[4 * i + 2]);
      inv_Li += (tmp_2 = d[2] * sum[4 * i + 3]);
		
      inv_Li = 1.0/inv_Li;
		
      dlnLidlz   = tmp_0 * e[0];
      d2lnLidlz2 = tmp_0 * e[1];
		
      dlnLidlz   += tmp_1 * e[2];
      d2lnLidlz2 += tmp_1 * e[3];
      
      dlnLidlz   += tmp_2 * e[4];	    	   	   	   
      d2lnLidlz2 += tmp_2 * e[5];
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      
      dlnLdlz   += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wr2ptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  free(d_start);
}



static void sumCATMULT(int tipCase, double *sum, double *x1_start, double *x2_start, double *tipVector,
		       char *tipX1, char *tipX2, int *modelptr, int loopStart, int n)
{
  double *x1, *x2;
  int i, model;

  switch(tipCase)
    {
    case TIP_TIP:
      for (i = loopStart; i < n; i++) 
	{    
	  model = modelptr[i];
	  x1 = &(tipVector[model * 64 + tipX1[i] * 4]);
	  x2 = &(tipVector[model * 64 + tipX2[i] * 4]);
	  
	  sum[i * 4]     = x1[0] * x2[0];
	  sum[i * 4 + 1] = x1[1] * x2[1];
	  sum[i * 4 + 2] = x1[2] * x2[2];
	  sum[i * 4 + 3] = x1[3] * x2[3];	    	    	  
	}    
      break;
    case TIP_INNER:
      for (i = loopStart; i < n; i++) 
	{
	  model = modelptr[i];
	  x1 = &(tipVector[model * 64 + 4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];
	  
	  sum[i * 4]     = x1[0] * x2[0];
	  sum[i * 4 + 1] = x1[1] * x2[1];
	  sum[i * 4 + 2] = x1[2] * x2[2];
	  sum[i * 4 + 3] = x1[3] * x2[3];				 
	}   
      break;
    case INNER_INNER:      
      for (i = loopStart; i < n; i++) 
	{     
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  
	  sum[i * 4]     = x1[0] * x2[0];
	  sum[i * 4 + 1] = x1[1] * x2[1];
	  sum[i * 4 + 2] = x1[2] * x2[2];
	  sum[i * 4 + 3] = x1[3] * x2[3];				   
	}
      break;
    default:
      assert(0);
    }  
}


static void coreGTRCATMULT(int lower, int upper, int numberOfCategories, int numberOfModels, double *sum, 
			   double *d1, double *d2, double *wrptr, double *wr2ptr, 
			   double *rptr, double *EIGN, int *cptr, int *modelptr, double lz)
{
  int i, model;
  double 
    *d, *d_start,
    tmp_0, tmp_1, tmp_2, inv_Li, dlnLidlz, d2lnLidlz2,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;
  double *e1, *e2, *e3, *s1, *s2, *s3;
   
  d = d_start = (double *)malloc(numberOfModels * numberOfCategories * 3 * sizeof(double));                 
  e1 = (double *)malloc(numberOfModels * sizeof(double));
  e2 = (double *)malloc(numberOfModels * sizeof(double));
  e3 = (double *)malloc(numberOfModels * sizeof(double));
  s1 = (double *)malloc(numberOfModels * sizeof(double));
  s2 = (double *)malloc(numberOfModels * sizeof(double));
  s3 = (double *)malloc(numberOfModels * sizeof(double));
 	  
  for(model = 0; model < numberOfModels; model++)
    {	       		
      double dd1, dd2, dd3;

      e1[model] = EIGN[model * 3]     * EIGN[model * 3];
      e2[model] = EIGN[model * 3 + 1] * EIGN[model * 3 + 1];
      e3[model] = EIGN[model * 3 + 2] * EIGN[model * 3 + 2];
      s1[model] = EIGN[model * 3];
      s2[model] = EIGN[model * 3 + 1];
      s3[model] = EIGN[model * 3 + 2];
	      
      dd1 =  EIGN[model * 3] * lz;
      dd2 =  EIGN[model * 3 + 1] * lz;
      dd3 =  EIGN[model * 3 + 2] * lz; 
      
      for(i = 0; i < numberOfCategories; i++)
	{	
	  *d++ = exp(dd1 * rptr[i]); 
	  *d++ = exp(dd2 * rptr[i]);
	  *d++ = exp(dd3 * rptr[i]);	    	    	  	   
	}       
    }	    	    

  for (i = lower; i < upper; i++) 	
    {	    	 	   	    
     model = modelptr[i];
	      
     d = &d_start[model * numberOfCategories * 3 + 3 * cptr[i]];
	      
     inv_Li = sum[4 * i];	   
     inv_Li += (tmp_0 = d[0] * sum[4 * i + 1]);
     inv_Li += (tmp_1 = d[1] * sum[4 * i + 2]);
     inv_Li += (tmp_2 = d[2] * sum[4 * i + 3]);
	      
     inv_Li = 1.0/inv_Li;	    	  	   	    
     
     dlnLidlz   = tmp_0 * s1[model];
     d2lnLidlz2 = tmp_0 * e1[model];
     
     dlnLidlz   += tmp_1 * s2[model];
     d2lnLidlz2 += tmp_1 * e2[model];
     
     dlnLidlz   += tmp_2 * s3[model];	    	   	   	   
     d2lnLidlz2 += tmp_2 * e3[model];
     
     dlnLidlz   *= inv_Li;
     d2lnLidlz2 *= inv_Li;
     
     dlnLdlz  += wrptr[i] * dlnLidlz;
     d2lnLdlz2 += wr2ptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  free(d_start);
  free(e1);
  free(e2);
  free(e3);
  free(s1);
  free(s2);
  free(s3);
}



static void sumGTRCATPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			  char *tipX1, char *tipX2, int lower, int n)
{
  int i;
  double *sum, *left, *right;

  switch(tipCase)
    {
    case TIP_TIP:          
      for (i = lower; i < n; i++) 
	{    
	  left  = &(tipVector[20 * tipX1[i]]);
	  right = &(tipVector[20 * tipX2[i]]);	
	  sum = &sumtable[20 * i];
	  
	  sum[0] = left[0] * right[0];			
	  sum[1] = left[1] * right[1];	
	  sum[2] = left[2] * right[2];			
	  sum[3] = left[3] * right[3];
	  sum[4] = left[4] * right[4];			
	  sum[5] = left[5] * right[5];	
	  sum[6] = left[6] * right[6];			
	  sum[7] = left[7] * right[7];
	  sum[8] = left[8] * right[8];			
	  sum[9] = left[9] * right[9];	
	  sum[10] = left[10] * right[10];			
	  sum[11] = left[11] * right[11];
	  sum[12] = left[12] * right[12];			
	  sum[13] = left[13] * right[13];	
	  sum[14] = left[14] * right[14];			
	  sum[15] = left[15] * right[15];
	  sum[16] = left[16] * right[16];			
	  sum[17] = left[17] * right[17];	
	  sum[18] = left[18] * right[18];			
	  sum[19] = left[19] * right[19];	
	}    
      break;
    case TIP_INNER:     
      for (i = lower; i < n; i++) 
	{    
	  left = &(tipVector[20 * tipX1[i]]);			       
	  right = &x2[20 * i];
	  sum = &sumtable[20 * i];
	  
	  sum[0] = left[0] * right[0];			
	  sum[1] = left[1] * right[1];	
	  sum[2] = left[2] * right[2];			
	  sum[3] = left[3] * right[3];
	  sum[4] = left[4] * right[4];			
	  sum[5] = left[5] * right[5];	
	  sum[6] = left[6] * right[6];			
	  sum[7] = left[7] * right[7];
	  sum[8] = left[8] * right[8];			
	  sum[9] = left[9] * right[9];	
	  sum[10] = left[10] * right[10];			
	  sum[11] = left[11] * right[11];
	  sum[12] = left[12] * right[12];			
	  sum[13] = left[13] * right[13];	
	  sum[14] = left[14] * right[14];			
	  sum[15] = left[15] * right[15];
	  sum[16] = left[16] * right[16];			
	  sum[17] = left[17] * right[17];	
	  sum[18] = left[18] * right[18];			
	  sum[19] = left[19] * right[19];				      
	}   
      break;
    case INNER_INNER:	     
      for (i = lower; i < n; i++) 
	{     
	  left  = &x1[20 * i];
	  right = &x2[20 * i];
	  sum = &sumtable[20 * i];
	  
	  sum[0] = left[0] * right[0];			
	  sum[1] = left[1] * right[1];	
	  sum[2] = left[2] * right[2];			
	  sum[3] = left[3] * right[3];
	  sum[4] = left[4] * right[4];			
	  sum[5] = left[5] * right[5];	
	  sum[6] = left[6] * right[6];			
	  sum[7] = left[7] * right[7];
	  sum[8] = left[8] * right[8];			
	  sum[9] = left[9] * right[9];	
	  sum[10] = left[10] * right[10];			
	  sum[11] = left[11] * right[11];
	  sum[12] = left[12] * right[12];			
	  sum[13] = left[13] * right[13];	
	  sum[14] = left[14] * right[14];			
	  sum[15] = left[15] * right[15];
	  sum[16] = left[16] * right[16];			
	  sum[17] = left[17] * right[17];	
	  sum[18] = left[18] * right[18];			
	  sum[19] = left[19] * right[19];			
	}   
      break;
    default:
      assert(0);
    }
}

static void coreGTRCATPROT(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int lower, int upper,
			   double *wrptr, double *wr2ptr, double *ext_dlnLdlz, double *ext_d2lnLdlz2, double *sumtable)
{
  int i, l;
  double *d1, *d_start, *sum; 
  double e[19], s[19]; 
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;
  double dd[19];
  double  dlnLdlz = 0.0;                 
  double  d2lnLdlz2 = 0.0; 

  for(l = 0; l < 19; l++)
    {
      e[l] = EIGN[l] * EIGN[l];
      s[l] = EIGN[l];
    }        
  
  d1 = d_start = (double *)malloc(numberOfCategories * 19 * sizeof(double));   
	
  for(l = 0; l < 19; l++)
    dd[l] = s[l] * lz; 

  for(i = 0; i < numberOfCategories; i++)
    {     
      for(l = 0; l < 19; l++)
	*d1++ = exp(dd[l] * rptr[i]); 	      	    	  	   
    }       

  for (i = lower; i < upper; i++) 	
    {	   	 	    
      d1 = &d_start[19 * cptr[i]];
      sum = &sumtable[20 * i];
	      
      inv_Li = sum[0];	   
      inv_Li += (tmp = d1[0] * sum[1]);
      dlnLidlz   = tmp * s[0];	   	   	    
      d2lnLidlz2 = tmp * e[0];
	      
      for(l = 1; l < 19; l++)	      
	{
	  inv_Li     += (tmp = d1[l] * sum[l + 1]);
	  dlnLidlz   += tmp * s[l];
	  d2lnLidlz2 += tmp * e[l];
	}	   	    
	      
      inv_Li = 1.0/inv_Li;  	    	    	    	    	    
	      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;

      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wr2ptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  free(d_start);
}


static void sumGTRCATPROTMULT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			      char *tipX1, char *tipX2, int *modelptr, int lower, int n)
{
  int i, model;
  double *sum, *left, *right;

  switch(tipCase)
    {
    case TIP_TIP:          
      for (i = lower; i < n; i++) 
	{    
	  model = modelptr[i];
	  left  = &(tipVector[model * 460 + 20 * tipX1[i]]);
	  right = &(tipVector[model * 460 + 20 * tipX2[i]]);	
	  sum = &sumtable[20 * i];
	  
	  sum[0] = left[0] * right[0];			
	  sum[1] = left[1] * right[1];	
	  sum[2] = left[2] * right[2];			
	  sum[3] = left[3] * right[3];
	  sum[4] = left[4] * right[4];			
	  sum[5] = left[5] * right[5];	
	  sum[6] = left[6] * right[6];			
	  sum[7] = left[7] * right[7];
	  sum[8] = left[8] * right[8];			
	  sum[9] = left[9] * right[9];	
	  sum[10] = left[10] * right[10];			
	  sum[11] = left[11] * right[11];
	  sum[12] = left[12] * right[12];			
	  sum[13] = left[13] * right[13];	
	  sum[14] = left[14] * right[14];			
	  sum[15] = left[15] * right[15];
	  sum[16] = left[16] * right[16];			
	  sum[17] = left[17] * right[17];	
	  sum[18] = left[18] * right[18];			
	  sum[19] = left[19] * right[19];	
	}    
      break;
    case TIP_INNER:     
      for (i = lower; i < n; i++) 
	{    
	  model = modelptr[i];
	  left = &(tipVector[model * 460 + 20 * tipX1[i]]);			       
	  right = &x2[20 * i];
	  sum = &sumtable[20 * i];
	  
	  sum[0] = left[0] * right[0];			
	  sum[1] = left[1] * right[1];	
	  sum[2] = left[2] * right[2];			
	  sum[3] = left[3] * right[3];
	  sum[4] = left[4] * right[4];			
	  sum[5] = left[5] * right[5];	
	  sum[6] = left[6] * right[6];			
	  sum[7] = left[7] * right[7];
	  sum[8] = left[8] * right[8];			
	  sum[9] = left[9] * right[9];	
	  sum[10] = left[10] * right[10];			
	  sum[11] = left[11] * right[11];
	  sum[12] = left[12] * right[12];			
	  sum[13] = left[13] * right[13];	
	  sum[14] = left[14] * right[14];			
	  sum[15] = left[15] * right[15];
	  sum[16] = left[16] * right[16];			
	  sum[17] = left[17] * right[17];	
	  sum[18] = left[18] * right[18];			
	  sum[19] = left[19] * right[19];				      
	}   
      break;
    case INNER_INNER:	     
      for (i = lower; i < n; i++) 
	{     
	  left  = &x1[20 * i];
	  right = &x2[20 * i];
	  sum = &sumtable[20 * i];
	  
	  sum[0] = left[0] * right[0];			
	  sum[1] = left[1] * right[1];	
	  sum[2] = left[2] * right[2];			
	  sum[3] = left[3] * right[3];
	  sum[4] = left[4] * right[4];			
	  sum[5] = left[5] * right[5];	
	  sum[6] = left[6] * right[6];			
	  sum[7] = left[7] * right[7];
	  sum[8] = left[8] * right[8];			
	  sum[9] = left[9] * right[9];	
	  sum[10] = left[10] * right[10];			
	  sum[11] = left[11] * right[11];
	  sum[12] = left[12] * right[12];			
	  sum[13] = left[13] * right[13];	
	  sum[14] = left[14] * right[14];			
	  sum[15] = left[15] * right[15];
	  sum[16] = left[16] * right[16];			
	  sum[17] = left[17] * right[17];	
	  sum[18] = left[18] * right[18];			
	  sum[19] = left[19] * right[19];			
	}   
      break;
    default:
      assert(0);
    }
}


static void coreGTRCATPROTMULT(double *EIGN, double lz, int numberOfCategories, double *rptr, int *cptr, int lower, int upper,
			       double *wrptr, double *wr2ptr, double *ext_dlnLdlz, double *ext_d2lnLdlz2, double *sumtable,
			       int numberOfModels, int *modelptr)
{
  double  *sum;
  int     i, l;
  double  dlnLdlz = 0.0;                 
  double  d2lnLdlz2 = 0.0;                  
  double *d1, *d_start; 
  double *e, *s, *e_start, *s_start; 
  int model; 
  double dd[19];	
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  e  = e_start =    (double *)malloc(19 * numberOfModels * sizeof(double));
  s  = s_start =    (double *)malloc(19 * numberOfModels * sizeof(double));            
  d1 = d_start =   (double *)malloc(numberOfModels * numberOfCategories * 19 * sizeof(double));  

  for(model = 0; model < numberOfModels; model++)
    {	    
      for(l = 0; l < 19; l++)
	{
	  e[model * 19 + l] = EIGN[model * 19 + l] * EIGN[model * 19 + l];
	  s[model * 19 + l] = EIGN[model * 19 + l];							       	
	  dd[l] = s[model * 19 + l] * lz;
	}		  	     
      
      for(i = 0; i < numberOfCategories; i++)
	{	 
	  for(l = 0; l < 19; l++)			
	    *d1++ = exp(dd[l] * rptr[i]); 	      	    	  	   		    
	}
    }
  
  
  for (i = lower; i < upper; i++) 	
    {	    	 	
      model = modelptr[i];	   
      sum = &sumtable[20 * i];
      
      d1 = &(d_start[model * numberOfCategories * 19 + 19 * cptr[i]]);	
      e  = &(e_start[model * 19]);
      s  = &(s_start[model * 19]);
      
      inv_Li = sum[0];	   
      inv_Li += (tmp = d1[0] * sum[1]);
      dlnLidlz   = tmp * s[0];	   	   	    
      d2lnLidlz2 = tmp * e[0];
      
      for(l = 1; l < 19; l++)	      
	{
	  inv_Li     += (tmp = d1[l] * sum[l + 1]);
	  dlnLidlz   += tmp * s[l];
	  d2lnLidlz2 += tmp * e[l];
	}	   	    
      
      inv_Li = 1.0/inv_Li;  	    	    	    	    	    
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wr2ptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }	
  
  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
  
  free(d_start);
  free(s_start);  
  free(e_start);
}


static void coreGTRGAMMAINVAR(double *invariants, double *frequencies, double *gammaRates, double *EIGN,
			      double *sumtable, double *ext_dlnLdlz, double *ext_d2lnLdlz2,
			      int *iptr, int *wrptr, int lower, int upper, double lz)
{   
  double  *sum, *diagptable, *diagptable_start;
  int     i;
    double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;    
  double freqs[4]; 
  double scaler =  0.25 * (1.0 - invariants[0]);   
  double tmp_1, tmp_2, tmp_3;
  double inv_Li, dlnLidlz, d2lnLidlz2;
 
 

  freqs[0] = frequencies[0] * invariants[0];
  freqs[1] = frequencies[1] * invariants[0];
  freqs[2] = frequencies[2] * invariants[0];
  freqs[3] = frequencies[3] * invariants[0];       
   
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 36);
  
  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];
      kisqr = ki * ki;
	      
      *diagptable++ = exp (EIGN[0] * ki * lz);
      *diagptable++ = exp (EIGN[1] * ki * lz);
      *diagptable++ = exp (EIGN[2] * ki * lz);
      
      *diagptable++ = EIGN[0] * ki;
      *diagptable++ = EIGN[0] * EIGN[0] * kisqr;
      
      *diagptable++ = EIGN[1] * ki;
      *diagptable++ = EIGN[1] * EIGN[1] * kisqr;
      
      *diagptable++ = EIGN[2] * ki;
      *diagptable++ = EIGN[2] * EIGN[2] * kisqr;
    }

  for (i = lower; i < upper; i++) 
    {	    	    	   
      diagptable = diagptable_start;
      sum = &(sumtable[i * 16]);

      inv_Li  = *sum++;
      inv_Li += (tmp_1 = *diagptable++ * *sum++);
      inv_Li += (tmp_2 = *diagptable++ * *sum++);
      inv_Li += (tmp_3 = *diagptable++ * *sum++);
      
      dlnLidlz   = tmp_1 * *diagptable++;
      d2lnLidlz2 = tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;	    
      
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;	    	    
      

      inv_Li += *sum++;
      inv_Li += (tmp_1 = *sum++ *  *diagptable++);
      inv_Li += (tmp_2 = *sum++ *  *diagptable++);
      inv_Li += (tmp_3 = *sum++ *  *diagptable++);	    	   	   	     	  	   	   
      
      dlnLidlz   += tmp_1 * *diagptable++;
      d2lnLidlz2 += tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;
      
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;
          
      inv_Li += *sum++;
      inv_Li += (tmp_1 = *sum++ *  *diagptable++);
      inv_Li += (tmp_2 = *sum++ *  *diagptable++);
      inv_Li += (tmp_3 = *sum++ *  *diagptable++);	   
      
      dlnLidlz   += tmp_1 * *diagptable++;
      d2lnLidlz2 += tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;
            
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;
      
      inv_Li += *sum++;
      inv_Li += (tmp_1 = *sum++ * *diagptable++);
      inv_Li += (tmp_2 = *sum++ * *diagptable++);
      inv_Li += (tmp_3 = *sum++ * *diagptable++);	          
      
      dlnLidlz   += tmp_1 * *diagptable++;
      d2lnLidlz2 += tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;
      
      
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;
     
      /* end */	    	
      
      inv_Li *= scaler;	   
      
      if(iptr[i] < 4)	      
	inv_Li += freqs[iptr[i]];
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLidlz *= scaler;
      d2lnLidlz2 *= scaler;
      
      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;
	  
  free(diagptable_start);
}




static void sumGAMMA(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
		     char *tipX1, char *tipX2, int lower, const int n)
{
  double *x1, *x2, *sum;
  int i;

  switch(tipCase)
    {
    case TIP_TIP:   
#ifdef _USE_OMP 
#pragma omp parallel for private(x1, x2, sum)      
      for (i = 0; i < n; i++) 
#else
      for (i = lower; i < n; i++) 
#endif
	{     
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);
	  sum = &sumtable[i * 16];
	  
	  sum[0] = x1[0] * x2[0];
	  sum[1] = x1[1] * x2[1];
	  sum[2] = x1[2] * x2[2];
	  sum[3] = x1[3] * x2[3];		    
	  
	  sum[4] = x1[0] * x2[0];
	  sum[5] = x1[1] * x2[1];
	  sum[6] = x1[2] * x2[2];
	  sum[7] = x1[3] * x2[3];
	  
	  sum[8]  = x1[0] * x2[0];
	  sum[9]  = x1[1] * x2[1];
	  sum[10] = x1[2] * x2[2];
	  sum[11] = x1[3] * x2[3];
	  
	  sum[12] = x1[0] * x2[0];
	  sum[13] = x1[1] * x2[1];
	  sum[14] = x1[2] * x2[2];
	  sum[15] = x1[3] * x2[3];	         
	}    
      break;
    case TIP_INNER: 
#ifdef _USE_OMP 
#pragma omp parallel for private(x1, x2, sum)
      for (i = 0; i < n; i++) 
#else
      for (i = lower; i < n; i++) 
#endif
	{     
	  x1  = &(tipVector[4 * tipX1[i]]);
	  x2  = &x2_start[16 * i];
	  sum = &sumtable[16 * i];
	  
	  sum[0] = x1[0] * x2[0];
	  sum[1] = x1[1] * x2[1];
	  sum[2] = x1[2] * x2[2];
	  sum[3] = x1[3] * x2[3];		    
	  
	  sum[4] = x1[0] * x2[4];
	  sum[5] = x1[1] * x2[5];
	  sum[6] = x1[2] * x2[6];
	  sum[7] = x1[3] * x2[7];
	  
	  sum[8]  = x1[0] * x2[8];
	  sum[9]  = x1[1] * x2[9];
	  sum[10] = x1[2] * x2[10];
	  sum[11] = x1[3] * x2[11];
	  
	  sum[12] = x1[0] * x2[12];
	  sum[13] = x1[1] * x2[13];
	  sum[14] = x1[2] * x2[14];
	  sum[15] = x1[3] * x2[15];	  	
	}   
      break;
    case INNER_INNER:
#ifdef _USE_OMP      
#pragma omp parallel for private(x1, x2, sum)
      for (i = 0; i < n; i++) 
#else
      for (i = lower; i < n; i++) 
#endif
	{     	      
	  x1  = &x1_start[16 * i];
	  x2  = &x2_start[16 * i];
	  sum = &sumtable[16 * i];
	  
	  sum[0] = x1[0] * x2[0];
	  sum[1] = x1[1] * x2[1];
	  sum[2] = x1[2] * x2[2];
	  sum[3] = x1[3] * x2[3];		    
	  
	  sum[4] = x1[4] * x2[4];
	  sum[5] = x1[5] * x2[5];
	  sum[6] = x1[6] * x2[6];
	  sum[7] = x1[7] * x2[7];
	  
	  sum[8]  = x1[8] * x2[8];
	  sum[9]  = x1[9] * x2[9];
	  sum[10] = x1[10] * x2[10];
	  sum[11] = x1[11] * x2[11];
	  
	  sum[12] = x1[12] * x2[12];
	  sum[13] = x1[13] * x2[13];
	  sum[14] = x1[14] * x2[14];
	  sum[15] = x1[15] * x2[15];	  	  
	}
      break;
    default:
      assert(0);
    }
}


static void coreGTRGAMMA(int lower, const int upper, double *sumtable, 
		  double *d1,  double *d2, double *EIGN, double *gammaRates, double lz, int *wrptr)
{
  int i;
  double 
    *diagptable, *diagptable_start, *sum, 
    tmp_1, tmp_2, tmp_3, inv_Li, dlnLidlz, d2lnLidlz2, ki, kisqr,
    dlnLdlz = 0.0,
    d2lnLdlz2 = 0.0;



  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 36);    
		 
  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];	 
      kisqr = ki * ki;
	      
      *diagptable++ = exp (EIGN[0] * ki * lz);
      *diagptable++ = exp (EIGN[1] * ki * lz);
      *diagptable++ = exp (EIGN[2] * ki * lz);
      
      *diagptable++ = EIGN[0] * ki;
      *diagptable++ = EIGN[0] * EIGN[0] * kisqr;
      
      *diagptable++ = EIGN[1] * ki;
      *diagptable++ = EIGN[1] * EIGN[1] * kisqr;
      
      *diagptable++ = EIGN[2] * ki;
      *diagptable++ = EIGN[2] * EIGN[2] * kisqr;
    }

#ifdef _USE_OMP
  {   
     double *di  = (double*)malloc(sizeof(double) * NumberOfThreads);
     double *di2 = (double*)malloc(sizeof(double) * NumberOfThreads);
#pragma omp parallel private(diagptable, sum, inv_Li, tmp_1, tmp_2, tmp_3, dlnLidlz, d2lnLidlz2) reduction(+ : dlnLdlz) reduction( + : d2lnLdlz2)
     {
       double private_di  = 0.0;
       double private_di2 = 0.0;
       int    ref = omp_get_thread_num();
#pragma omp for
       for (i = 0; i < upper; i++)
#else
       for (i = lower; i < upper; i++) 
#endif
	 {	   	    	   	    
	   diagptable = diagptable_start;
	   sum = &(sumtable[i * 16]);
	   
	   inv_Li  = *sum++;
	   inv_Li += (tmp_1 = *diagptable++ * *sum++);
	   inv_Li += (tmp_2 = *diagptable++ * *sum++);
	   inv_Li += (tmp_3 = *diagptable++ * *sum++);
	   
	   dlnLidlz   = tmp_1 * *diagptable++;
	   d2lnLidlz2 = tmp_1 * *diagptable++;
	   
	   dlnLidlz   += tmp_2 * *diagptable++;
	   d2lnLidlz2 += tmp_2 * *diagptable++;	    
	   
	   dlnLidlz   += tmp_3 * *diagptable++;
	   d2lnLidlz2 += tmp_3 * *diagptable++;	    	    
	   
	   
	   inv_Li += *sum++;
	   inv_Li += (tmp_1 = *sum++ *  *diagptable++);
	   inv_Li += (tmp_2 = *sum++ *  *diagptable++);
	   inv_Li += (tmp_3 = *sum++ *  *diagptable++);	    	   	   	     	  	   	   
	   
	   dlnLidlz   += tmp_1 * *diagptable++;
	   d2lnLidlz2 += tmp_1 * *diagptable++;
	   
	   dlnLidlz   += tmp_2 * *diagptable++;
	   d2lnLidlz2 += tmp_2 * *diagptable++;
	   
	   dlnLidlz   += tmp_3 * *diagptable++;
	   d2lnLidlz2 += tmp_3 * *diagptable++;
	   
	   inv_Li += *sum++;
	   inv_Li += (tmp_1 = *sum++ *  *diagptable++);
	   inv_Li += (tmp_2 = *sum++ *  *diagptable++);
	   inv_Li += (tmp_3 = *sum++ *  *diagptable++);	   
	   
	   dlnLidlz   += tmp_1 * *diagptable++;
	   d2lnLidlz2 += tmp_1 * *diagptable++;
	   
	   dlnLidlz   += tmp_2 * *diagptable++;
	   d2lnLidlz2 += tmp_2 * *diagptable++;
	   
	   dlnLidlz   += tmp_3 * *diagptable++;
	   d2lnLidlz2 += tmp_3 * *diagptable++;
	   
	   inv_Li += *sum++;
	   inv_Li += (tmp_1 = *sum++ * *diagptable++);
	   inv_Li += (tmp_2 = *sum++ * *diagptable++);
	   inv_Li += (tmp_3 = *sum++ * *diagptable++);	          
	   
	   dlnLidlz   += tmp_1 * *diagptable++;
	   d2lnLidlz2 += tmp_1 * *diagptable++;
	   
	   dlnLidlz   += tmp_2 * *diagptable++;
	   d2lnLidlz2 += tmp_2 * *diagptable++;
	   	   
	   dlnLidlz   += tmp_3 * *diagptable++;
	   d2lnLidlz2 += tmp_3 * *diagptable++;
	   	   	   
	   inv_Li = 1.0 / inv_Li;
      
	   dlnLidlz   *= inv_Li;
	   d2lnLidlz2 *= inv_Li;
#ifndef _USE_OMP	   
	   dlnLdlz  += wrptr[i] * dlnLidlz;
	   d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
	 }
#else
           private_di  += wrptr[i] * dlnLidlz;
	   private_di2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
         }
		
       di[ref]  = private_di;
       di2[ref] = private_di2;
    }

    for(i = 0; i < NumberOfThreads; i++)
      {
        dlnLdlz   += di[i];
        d2lnLdlz2 += di2[i];
      }
    free(di);
    free(di2);
  }
#endif     
  

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  free(diagptable_start);
}




/*************************************************************************************/

static void sumGAMMAMULT(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
			 char *tipX1, char *tipX2, int *modelptr, int lower, int n)
{
  double *x1, *x2, *sum;
  int i;

  switch(tipCase)
    {
    case TIP_TIP:         
      for (i = lower; i < n; i++) 
	{     
	  x1 = &(tipVector[64 * modelptr[i] + 4 * tipX1[i]]);
	  x2 = &(tipVector[64 * modelptr[i] + 4 * tipX2[i]]);
	  sum = &sumtable[i * 16];
	  
	  sum[0] = x1[0] * x2[0];
	  sum[1] = x1[1] * x2[1];
	  sum[2] = x1[2] * x2[2];
	  sum[3] = x1[3] * x2[3];		    
	  
	  sum[4] = x1[0] * x2[0];
	  sum[5] = x1[1] * x2[1];
	  sum[6] = x1[2] * x2[2];
	  sum[7] = x1[3] * x2[3];
	  
	  sum[8]  = x1[0] * x2[0];
	  sum[9]  = x1[1] * x2[1];
	  sum[10] = x1[2] * x2[2];
	  sum[11] = x1[3] * x2[3];
	  
	  sum[12] = x1[0] * x2[0];
	  sum[13] = x1[1] * x2[1];
	  sum[14] = x1[2] * x2[2];
	  sum[15] = x1[3] * x2[3];	         
	}    
      break;
    case TIP_INNER:          
      for (i = lower; i < n; i++) 
	{     
	  x1  = &(tipVector[64 * modelptr[i] + 4 * tipX1[i]]);
	  x2  = &x2_start[16 * i];
	  sum = &sumtable[16 * i];
	  
	  sum[0] = x1[0] * x2[0];
	  sum[1] = x1[1] * x2[1];
	  sum[2] = x1[2] * x2[2];
	  sum[3] = x1[3] * x2[3];		    
	  
	  sum[4] = x1[0] * x2[4];
	  sum[5] = x1[1] * x2[5];
	  sum[6] = x1[2] * x2[6];
	  sum[7] = x1[3] * x2[7];
	  
	  sum[8]  = x1[0] * x2[8];
	  sum[9]  = x1[1] * x2[9];
	  sum[10] = x1[2] * x2[10];
	  sum[11] = x1[3] * x2[11];
	  
	  sum[12] = x1[0] * x2[12];
	  sum[13] = x1[1] * x2[13];
	  sum[14] = x1[2] * x2[14];
	  sum[15] = x1[3] * x2[15];	  	
	}   
      break;
    case INNER_INNER:	  
      for (i = lower; i < n; i++) 
	{     	      
	  x1  = &x1_start[16 * i];
	  x2  = &x2_start[16 * i];
	  sum = &sumtable[16 * i];
	  
	  sum[0] = x1[0] * x2[0];
	  sum[1] = x1[1] * x2[1];
	  sum[2] = x1[2] * x2[2];
	  sum[3] = x1[3] * x2[3];		    
	  
	  sum[4] = x1[4] * x2[4];
	  sum[5] = x1[5] * x2[5];
	  sum[6] = x1[6] * x2[6];
	  sum[7] = x1[7] * x2[7];
	  
	  sum[8]  = x1[8] * x2[8];
	  sum[9]  = x1[9] * x2[9];
	  sum[10] = x1[10] * x2[10];
	  sum[11] = x1[11] * x2[11];
	  
	  sum[12] = x1[12] * x2[12];
	  sum[13] = x1[13] * x2[13];
	  sum[14] = x1[14] * x2[14];
	  sum[15] = x1[15] * x2[15];	  	  
	}
      break;
    default:
      assert(0);
    }  
}

static void coreGTRGAMMAMULT(int lower, int upper, double *sumtable, 
			     double *ext_dlnLdlz,  double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wrptr,
			     int numberOfModels, int *modelptr)
{
  double  *sum, *diagptable, *diagptable_start;
  int     i;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double tmp_1, tmp_2, tmp_3;
  double inv_Li, dlnLidlz, d2lnLidlz2;
  double ki, kisqr;    
  int model;
    
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 36 * numberOfModels);	    	  
	    
  for(model = 0; model < numberOfModels; model++)
    {	      		
      diagptable = &diagptable_start[model * 36];
      
      for(i = 0; i < 4; i++)
	{
	  ki = gammaRates[model * 4 + i];	 
	  kisqr = ki * ki;
	  
	  *diagptable++ = exp (EIGN[model * 3]     * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 1] * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 2] * ki * lz);
	  
	  *diagptable++ = EIGN[model * 3] * ki;
	  *diagptable++ = EIGN[model * 3] * EIGN[model * 3] * kisqr;
	  
	  *diagptable++ = EIGN[model * 3 + 1] * ki;
	  *diagptable++ = EIGN[model * 3 + 1] * EIGN[model * 3 + 1] * kisqr;
	  
	  *diagptable++ = EIGN[model * 3 + 2] * ki;
	  *diagptable++ = EIGN[model * 3 + 2] * EIGN[model * 3 + 2] * kisqr;
	}
    }
  
  for (i = lower; i < upper; i++) 
    {	   	    	   
      diagptable = &diagptable_start[36 * modelptr[i]];
      sum = &(sumtable[i * 16]);
      
      inv_Li  = *sum++;
      inv_Li += (tmp_1 = *diagptable++ * *sum++);
      inv_Li += (tmp_2 = *diagptable++ * *sum++);
      inv_Li += (tmp_3 = *diagptable++ * *sum++);
      
      dlnLidlz   = tmp_1 * *diagptable++;
      d2lnLidlz2 = tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;	    
      
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;	    	    
      

      inv_Li += *sum++;
      inv_Li += (tmp_1 = *sum++ *  *diagptable++);
      inv_Li += (tmp_2 = *sum++ *  *diagptable++);
      inv_Li += (tmp_3 = *sum++ *  *diagptable++);	    	   	   	     	  	   	   
      
      dlnLidlz   += tmp_1 * *diagptable++;
      d2lnLidlz2 += tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;
      
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;
          
      inv_Li += *sum++;
      inv_Li += (tmp_1 = *sum++ *  *diagptable++);
      inv_Li += (tmp_2 = *sum++ *  *diagptable++);
      inv_Li += (tmp_3 = *sum++ *  *diagptable++);	   
      
      dlnLidlz   += tmp_1 * *diagptable++;
      d2lnLidlz2 += tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;
            
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;
      
      inv_Li += *sum++;
      inv_Li += (tmp_1 = *sum++ * *diagptable++);
      inv_Li += (tmp_2 = *sum++ * *diagptable++);
      inv_Li += (tmp_3 = *sum++ * *diagptable++);	          
      
      dlnLidlz   += tmp_1 * *diagptable++;
      d2lnLidlz2 += tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;
      
      
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  free(diagptable_start);	  
}



static void coreGTRGAMMAMULTINVAR(int lower, int upper, double *sumtable, 
				  double *ext_dlnLdlz,  double *ext_d2lnLdlz2, double *EIGN, double *gammaRates, double lz, int *wrptr,
				  int numberOfModels, int *modelptr, double *frequencies, double *invariants, int *iptr)
{
  double  *sum, *diagptable, *diagptable_start;
  int     i;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;     
  int model;
  double *scalers;
  double *freqs; 
  double tmp_1, tmp_2, tmp_3;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  diagptable_start = (double *)malloc(sizeof(double) * 36 * numberOfModels);
  freqs            = (double *)malloc(4 * numberOfModels * sizeof(double));
  scalers          = (double *)malloc(numberOfModels * sizeof(double));

  for(model = 0; model < numberOfModels; model++)
    {
      scalers[model] = 0.25 * (1.0 - invariants[model]); 
      freqs[4 * model]     = frequencies[4 * model]     * invariants[model]; 
      freqs[4 * model + 1] = frequencies[4 * model + 1] * invariants[model];
      freqs[4 * model + 2] = frequencies[4 * model + 2] * invariants[model];
      freqs[4 * model + 3] = frequencies[4 * model + 3] * invariants[model];
    }

  for(model = 0; model < numberOfModels; model++)
    {	    	      
      diagptable = &diagptable_start[model * 36];
      
      for(i = 0; i < 4; i++)
	{
	  ki = gammaRates[model * 4 + i];	 
	  kisqr = ki * ki;
		  
	  *diagptable++ = exp (EIGN[model * 3]     * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 1] * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 2] * ki * lz);
		  
	  *diagptable++ = EIGN[model * 3] * ki;
	  *diagptable++ = EIGN[model * 3] * EIGN[model * 3] * kisqr;
	  
	  *diagptable++ = EIGN[model * 3 + 1] * ki;
	  *diagptable++ = EIGN[model * 3 + 1] * EIGN[model * 3 + 1] * kisqr;
	  
	  *diagptable++ = EIGN[model * 3 + 2] * ki;
	  *diagptable++ = EIGN[model * 3 + 2] * EIGN[model * 3 + 2] * kisqr;
	}
    }	    	    

  for (i = lower; i < upper; i++) 
    {	   	    	   
      diagptable = &diagptable_start[36 * modelptr[i]];
      sum = &(sumtable[i * 16]);

      inv_Li  = *sum++;
      inv_Li += (tmp_1 = *diagptable++ * *sum++);
      inv_Li += (tmp_2 = *diagptable++ * *sum++);
      inv_Li += (tmp_3 = *diagptable++ * *sum++);
      
      dlnLidlz   = tmp_1 * *diagptable++;
      d2lnLidlz2 = tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;	    
      
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;	    	    
      

      inv_Li += *sum++;
      inv_Li += (tmp_1 = *sum++ *  *diagptable++);
      inv_Li += (tmp_2 = *sum++ *  *diagptable++);
      inv_Li += (tmp_3 = *sum++ *  *diagptable++);	    	   	   	     	  	   	   
      
      dlnLidlz   += tmp_1 * *diagptable++;
      d2lnLidlz2 += tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;
      
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;
          
      inv_Li += *sum++;
      inv_Li += (tmp_1 = *sum++ *  *diagptable++);
      inv_Li += (tmp_2 = *sum++ *  *diagptable++);
      inv_Li += (tmp_3 = *sum++ *  *diagptable++);	   
      
      dlnLidlz   += tmp_1 * *diagptable++;
      d2lnLidlz2 += tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;
            
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;
      
      inv_Li += *sum++;
      inv_Li += (tmp_1 = *sum++ * *diagptable++);
      inv_Li += (tmp_2 = *sum++ * *diagptable++);
      inv_Li += (tmp_3 = *sum++ * *diagptable++);	          
      
      dlnLidlz   += tmp_1 * *diagptable++;
      d2lnLidlz2 += tmp_1 * *diagptable++;
      
      dlnLidlz   += tmp_2 * *diagptable++;
      d2lnLidlz2 += tmp_2 * *diagptable++;
      
      
      dlnLidlz   += tmp_3 * *diagptable++;
      d2lnLidlz2 += tmp_3 * *diagptable++;
     
      /* end */	    	
      
      inv_Li *= scalers[modelptr[i]];
      
      if(iptr[i] < 4)	      
	inv_Li += freqs[4 * modelptr[i] + iptr[i]];
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLidlz   *= scalers[modelptr[i]];
      d2lnLidlz2 *= scalers[modelptr[i]];
      
      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  free(diagptable_start);
  free(freqs);
  free(scalers);
}



static void sumGAMMAPROT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
		     char *tipX1, char *tipX2, int lower, int n)
{
  int i, l;
  double *left, *right, *sum; 

  switch(tipCase)
    {
    case TIP_TIP:    
      for(i = lower; i < n; i++) 
	{     
	  left  = &(tipVector[20 * tipX1[i]]);
	  right = &(tipVector[20 * tipX2[i]]);      
	  sum   = &(sumtable[i * 80]);	 

	  for(l = 0; l < 4; l++)
	    {	     
	      *sum++ = left[0] * right[0];
	      *sum++ = left[1] * right[1];
	      *sum++ = left[2] * right[2];
	      *sum++ = left[3] * right[3];
	      *sum++ = left[4] * right[4];
	      *sum++ = left[5] * right[5];
	      *sum++ = left[6] * right[6];
	      *sum++ = left[7] * right[7];
	      *sum++ = left[8] * right[8];
	      *sum++ = left[9] * right[9];
	      *sum++ = left[10] * right[10];
	      *sum++ = left[11] * right[11];
	      *sum++ = left[12] * right[12];
	      *sum++ = left[13] * right[13];
	      *sum++ = left[14] * right[14];
	      *sum++ = left[15] * right[15];
	      *sum++ = left[16] * right[16];
	      *sum++ = left[17] * right[17];
	      *sum++ = left[18] * right[18];
	      *sum++ = left[19] * right[19];
	    }	    		    	 
	}    
      break;
    case TIP_INNER:  
      for(i = lower; i < n; i++) 
	{     
	  left = &(tipVector[20 * tipX1[i]]);		    	  
	  sum  = &(sumtable[i * 80]);	 

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[80 * i + l * 20]);

	      *sum++ = left[0] * right[0];
	      *sum++ = left[1] * right[1];
	      *sum++ = left[2] * right[2];
	      *sum++ = left[3] * right[3];
	      *sum++ = left[4] * right[4];
	      *sum++ = left[5] * right[5];
	      *sum++ = left[6] * right[6];
	      *sum++ = left[7] * right[7];
	      *sum++ = left[8] * right[8];
	      *sum++ = left[9] * right[9];
	      *sum++ = left[10] * right[10];
	      *sum++ = left[11] * right[11];
	      *sum++ = left[12] * right[12];
	      *sum++ = left[13] * right[13];
	      *sum++ = left[14] * right[14];
	      *sum++ = left[15] * right[15];
	      *sum++ = left[16] * right[16];
	      *sum++ = left[17] * right[17];
	      *sum++ = left[18] * right[18];
	      *sum++ = left[19] * right[19];
	    }	    	      	
	}    
      break;
    case INNER_INNER:
      for(i = lower; i < n; i++) 
	{     		    
	  sum = &sumtable[i * 80];	
	  
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[80 * i + l * 20]);
	      right = &(x2[80 * i + l * 20]);
	      
	      *sum++ = left[0] * right[0];	
	      *sum++ = left[1] * right[1];
	      *sum++ = left[2] * right[2];
	      *sum++ = left[3] * right[3];
	      *sum++ = left[4] * right[4];
	      *sum++ = left[5] * right[5];
	      *sum++ = left[6] * right[6];
	      *sum++ = left[7] * right[7];
	      *sum++ = left[8] * right[8];
	      *sum++ = left[9] * right[9];
	      *sum++ = left[10] * right[10];
	      *sum++ = left[11] * right[11];
	      *sum++ = left[12] * right[12];
	      *sum++ = left[13] * right[13];
	      *sum++ = left[14] * right[14];
	      *sum++ = left[15] * right[15];
	      *sum++ = left[16] * right[16];
	      *sum++ = left[17] * right[17];
	      *sum++ = left[18] * right[18];
	      *sum++ = left[19] * right[19];
	    }	  	
	}
      break;
    default:
      assert(0);
    }
}


static void coreGTRGAMMAPROT(double *gammaRates, double *EIGN, double *sumtable, int lower, int upper, int *wrptr,
			     double *ext_dlnLdlz, double *ext_d2lnLdlz2, double lz)
{
  double  *sum, *diagptable, *diagptable_start;
  int     i, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr; 
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 228);

  for(i = 0; i < 4; i++)
    {	   
      ki = gammaRates[i];	 
      kisqr = ki * ki;
      
      for(l = 0; l < 19; l++)
	{	      
	  *diagptable++ = exp (EIGN[l] * ki * lz);
	  *diagptable++ = EIGN[l] * ki;
	  *diagptable++ = EIGN[l] * EIGN[l] * kisqr;
	}			    	      
    }

  for (i = lower; i < upper; i++) 
    {	  	    	   
      diagptable = diagptable_start;
      sum = &sumtable[i * 80];    
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;
        
      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}

      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}
      
      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}
      	   	   	 
      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    } 

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  free(diagptable_start);
}


static void coreGTRGAMMAPROTINVAR(double *gammaRates, double *EIGN, double *sumtable, int lower, int upper, int *wrptr,
				  double *ext_dlnLdlz, double *ext_d2lnLdlz2, double lz, double *frequencies,
				  double *invariants, int *iptr)
{
  double  *sum, *diagptable, *diagptable_start;
  int     i, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;    
  double freqs[20]; 
  double scaler =  0.25 * (1.0 - invariants[0]);
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;

  for(i = 0; i < 20; i++)    
    freqs[i] = frequencies[i] * invariants[0];
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 228);

   for(i = 0; i < 4; i++)
     {	   
       ki = gammaRates[i];	 
       kisqr = ki * ki;
       
       for(l = 0; l < 19; l++)
	 {	      
	   *diagptable++ = exp (EIGN[l] * ki * lz);
	   *diagptable++ = EIGN[l] * ki;
	   *diagptable++ = EIGN[l] * EIGN[l] * kisqr;
	 }			    	      
     }               
   
   for(i = lower; i < upper; i++) 
     {	  	    	   
       diagptable = diagptable_start;
       sum = &sumtable[i * 80];    
       inv_Li   = 0.0;
       dlnLidlz = 0.0;
       d2lnLidlz2 = 0.0;
       
       inv_Li += *sum++;           
       for(l = 1; l < 20; l++)
	 {		
	   inv_Li     += (tmp = *diagptable++  * *sum++);
	   dlnLidlz   += tmp * *diagptable++;
	   d2lnLidlz2 += tmp * *diagptable++;				
	 }
       
       inv_Li += *sum++;           
       for(l = 1; l < 20; l++)
	 {		
	   inv_Li     += (tmp = *diagptable++  * *sum++);
	   dlnLidlz   += tmp * *diagptable++;
	   d2lnLidlz2 += tmp * *diagptable++;				
	 }
       
       inv_Li += *sum++;           
       for(l = 1; l < 20; l++)
	 {		
	   inv_Li     += (tmp = *diagptable++  * *sum++);
	   dlnLidlz   += tmp * *diagptable++;
	   d2lnLidlz2 += tmp * *diagptable++;				
	 }
       
       inv_Li += *sum++;           
       for(l = 1; l < 20; l++)
	 {		
	   inv_Li     += (tmp = *diagptable++  * *sum++);
	   dlnLidlz   += tmp * *diagptable++;
	   d2lnLidlz2 += tmp * *diagptable++;				
	 }       
       
       inv_Li *= scaler;
       
       if(iptr[i] < 20)	      
	 inv_Li += freqs[iptr[i]];
       
       inv_Li = 1.0 / inv_Li;
       
       dlnLidlz   *= inv_Li;
       d2lnLidlz2 *= inv_Li;
       
       dlnLidlz *= scaler;
       d2lnLidlz2 *= scaler;
       
       dlnLdlz  += wrptr[i] * dlnLidlz;
       d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
     }

   *ext_dlnLdlz   = dlnLdlz;
   *ext_d2lnLdlz2 = d2lnLdlz2; 
   
   free(diagptable_start);	  
}



static void sumGAMMAPROTMULT(int tipCase, double *sumtable, double *x1, double *x2, double *tipVector,
			     char *tipX1, char *tipX2, int *modelptr, int lower, int n)
{
  int i, l;
  double *left, *right, *sum;

  switch(tipCase)
    {
    case TIP_TIP:    
      for(i = lower; i < n; i++) 
	{     
	  left  = &(tipVector[460 * modelptr[i] + 20 * tipX1[i]]);
	  right = &(tipVector[460 * modelptr[i] + 20 * tipX2[i]]);      
	  sum   = &(sumtable[i * 80]);	 

	  for(l = 0; l < 4; l++)
	    {	     
	      *sum++ = left[0] * right[0];
	      *sum++ = left[1] * right[1];
	      *sum++ = left[2] * right[2];
	      *sum++ = left[3] * right[3];
	      *sum++ = left[4] * right[4];
	      *sum++ = left[5] * right[5];
	      *sum++ = left[6] * right[6];
	      *sum++ = left[7] * right[7];
	      *sum++ = left[8] * right[8];
	      *sum++ = left[9] * right[9];
	      *sum++ = left[10] * right[10];
	      *sum++ = left[11] * right[11];
	      *sum++ = left[12] * right[12];
	      *sum++ = left[13] * right[13];
	      *sum++ = left[14] * right[14];
	      *sum++ = left[15] * right[15];
	      *sum++ = left[16] * right[16];
	      *sum++ = left[17] * right[17];
	      *sum++ = left[18] * right[18];
	      *sum++ = left[19] * right[19];
	    }	    		    	 
	}    
      break;
    case TIP_INNER:  
      for(i = lower; i < n; i++) 
	{     
	  left = &(tipVector[460 * modelptr[i] + 20 * tipX1[i]]);		    	  
	  sum  = &(sumtable[i * 80]);	 

	  for(l = 0; l < 4; l++)
	    {
	      right = &(x2[80 * i + l * 20]);

	      *sum++ = left[0] * right[0];
	      *sum++ = left[1] * right[1];
	      *sum++ = left[2] * right[2];
	      *sum++ = left[3] * right[3];
	      *sum++ = left[4] * right[4];
	      *sum++ = left[5] * right[5];
	      *sum++ = left[6] * right[6];
	      *sum++ = left[7] * right[7];
	      *sum++ = left[8] * right[8];
	      *sum++ = left[9] * right[9];
	      *sum++ = left[10] * right[10];
	      *sum++ = left[11] * right[11];
	      *sum++ = left[12] * right[12];
	      *sum++ = left[13] * right[13];
	      *sum++ = left[14] * right[14];
	      *sum++ = left[15] * right[15];
	      *sum++ = left[16] * right[16];
	      *sum++ = left[17] * right[17];
	      *sum++ = left[18] * right[18];
	      *sum++ = left[19] * right[19];
	    }	    	      	
	}    
      break;
    case INNER_INNER:
      for(i = lower; i < n; i++) 
	{     		    
	  sum = &sumtable[i * 80];	
	  
	  for(l = 0; l < 4; l++)
	    {
	      left  = &(x1[80 * i + l * 20]);
	      right = &(x2[80 * i + l * 20]);
	      
	      *sum++ = left[0] * right[0];	
	      *sum++ = left[1] * right[1];
	      *sum++ = left[2] * right[2];
	      *sum++ = left[3] * right[3];
	      *sum++ = left[4] * right[4];
	      *sum++ = left[5] * right[5];
	      *sum++ = left[6] * right[6];
	      *sum++ = left[7] * right[7];
	      *sum++ = left[8] * right[8];
	      *sum++ = left[9] * right[9];
	      *sum++ = left[10] * right[10];
	      *sum++ = left[11] * right[11];
	      *sum++ = left[12] * right[12];
	      *sum++ = left[13] * right[13];
	      *sum++ = left[14] * right[14];
	      *sum++ = left[15] * right[15];
	      *sum++ = left[16] * right[16];
	      *sum++ = left[17] * right[17];
	      *sum++ = left[18] * right[18];
	      *sum++ = left[19] * right[19];
	    }	  	
	}
      break;
    default:
      assert(0);
    } 
}

static void coreGTRGAMMAPROTMULT(double *gammaRates, double *EIGN, double *sumtable, int lower, int upper, int *wrptr,
				 double *ext_dlnLdlz, double *ext_d2lnLdlz2, double lz, int numberOfModels, int *modelptr)
{
  double  *sum, *diagptable, *diagptable_start;
  int     i, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr; 
  double tmp;
  double inv_Li, dlnLidlz, d2lnLidlz2;
  int model;            
        
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 228 * numberOfModels); 

  for(model = 0; model < numberOfModels; model++)
    {	    	      
      diagptable = &diagptable_start[model * 228];
      
      for(i = 0; i < 4; i++)
	{	   
	  ki = gammaRates[model * 4 + i];	 
	  kisqr = ki * ki;
	  
	  for(l = 0; l < 19; l++)
	    {	      
	      *diagptable++ = exp (EIGN[model * 19 + l] * ki * lz);
	      *diagptable++ = EIGN[model * 19 + l] * ki;
	      *diagptable++ = EIGN[model * 19 + l] * EIGN[model * 19 + l] * kisqr;
	    }			    	      
	}
    }	  	 
	  
  for (i = lower; i < upper; i++) 
    {	  		     	      
      diagptable = &diagptable_start[228 * modelptr[i]];
      sum = &sumtable[i * 80];    
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;
        
      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}

      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}
      
      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}
      	   	   	 
      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }

  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  free(diagptable_start);
}


static void coreGTRGAMMAPROTMULTINVAR(double *gammaRates, double *EIGN, double *sumtable, int lower, int upper, int *wrptr,
				      double *ext_dlnLdlz, double *ext_d2lnLdlz2, double lz, int numberOfModels, int *modelptr,
				      double *invariants, double *frequencies, int *iptr)
{
  double  *sum, *diagptable, *diagptable_start;
  int     i, l;
  double  dlnLdlz = 0;
  double d2lnLdlz2 = 0;
  double ki, kisqr;    
  int model;
  double *scalers;
  double *freqs;   
  double tmp;  
  double inv_Li, dlnLidlz, d2lnLidlz2;

  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 228 * numberOfModels);
  freqs   = (double *)malloc(20 * numberOfModels * sizeof(double));
  scalers = (double *)malloc(numberOfModels * sizeof(double));
      
  for(model = 0; model < numberOfModels; model++)
    {
      scalers[model] = 0.25 * (1.0 - invariants[model]); 
      for(i = 0; i < 20; i++)
	freqs[20 * model + i]     = frequencies[20 * model + i] * invariants[model];       
    }
  
  for(model = 0; model < numberOfModels; model++)
    {	     	      
      diagptable = &diagptable_start[model * 228];
      
      for(i = 0; i < 4; i++)
	{	   
	  ki = gammaRates[model * 4 + i];	 
	  kisqr = ki * ki;
	  
	  for(l = 0; l < 19; l++)
	    {	      
	      *diagptable++ = exp (EIGN[model * 19 + l] * ki * lz);
	      *diagptable++ = EIGN[model * 19 + l] * ki;
	      *diagptable++ = EIGN[model * 19 + l] * EIGN[model * 19 + l] * kisqr;
	    }			    	      
	}
    }	  	 
  
  for(i = lower; i < upper; i++) 
    {	  
      model = modelptr[i];      
      diagptable = &diagptable_start[228 * model];
      sum = &sumtable[i * 80];    
      inv_Li   = 0.0;
      dlnLidlz = 0.0;
      d2lnLidlz2 = 0.0;
        
      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}

      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}
      
      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}
      	   	   	 
      inv_Li += *sum++;           
      for(l = 1; l < 20; l++)
	{		
	  inv_Li     += (tmp = *diagptable++  * *sum++);
	  dlnLidlz   += tmp * *diagptable++;
	  d2lnLidlz2 += tmp * *diagptable++;				
	}       	 
      
      inv_Li *= scalers[model];
      
      if(iptr[i] < 20)	      
	inv_Li += freqs[20 * model + iptr[i]];
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLidlz   *= scalers[model];
      d2lnLidlz2 *= scalers[model];
      
      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }
  
  *ext_dlnLdlz   = dlnLdlz;
  *ext_d2lnLdlz2 = d2lnLdlz2;

  free(diagptable_start);
  free(freqs);
  free(scalers);
}


static void makenewzMixedData(int model, int pNumber, int qNumber, tree *tr)
{
  double   
    *x1_start = (double*)NULL, 
    *x2_start = (double*)NULL;
  char 
    *tipX1 = (char*)NULL,
    *tipX2 = (char*)NULL;      
  int tipCase;

  /*int l = tr->modelIndices[model][0];
    int u = tr->modelIndices[model][1];	  	*/
  int l = tr->partitionData[model].lower;
  int u = tr->partitionData[model].upper;

  int width  = u - l;
  /*int offset = tr->modelOffsets[model]; */
  int offset = tr->partitionData[model].modelOffset;

  if(isTip(pNumber, tr->rdta->numsp) || isTip(qNumber, tr->rdta->numsp))
    {	    	    	   	               
      if(!( isTip(pNumber, tr->rdta->numsp) && isTip(qNumber, tr->rdta->numsp)) )
	{	 
	  tipCase = TIP_INNER;
	  if(isTip(qNumber, tr->rdta->numsp))
	    {	
	      tipX1 = tr->yVector[qNumber];
	      tipX1 = &tipX1[l];
	      
	      x2_start = getLikelihoodArray(pNumber,  tr->mxtips, tr->xVector);
	      x2_start = &x2_start[offset];
	    }	    
	  else
	    {
	      tipX1 = tr->yVector[pNumber];
	      tipX1 = &tipX1[l];

	      x2_start = getLikelihoodArray(qNumber,  tr->mxtips, tr->xVector);
	      x2_start = &x2_start[offset];
	    }
	}
      else
	{
	  tipCase = TIP_TIP;
	  tipX1 = tr->yVector[pNumber];
	  tipX1 = &tipX1[l];

	  tipX2 = tr->yVector[qNumber];
	  tipX2 = &tipX2[l];
	}
    }
  else
    {
      tipCase = INNER_INNER;
            
      x1_start = getLikelihoodArray(pNumber, tr->mxtips, tr->xVector);
      x1_start = &x1_start[offset];

      x2_start = getLikelihoodArray(qNumber, tr->mxtips, tr->xVector);
      x2_start = &x2_start[offset];
    }

  switch(tr->partitionData[model].dataType)
    { 
    case DNA_DATA:
      switch(tr->rateHetModel)
	{
	case CAT:
	  sumCAT(tipCase, &(tr->sumBuffer[offset]), x1_start, 
		 x2_start, &(tr->tipVectorDNA[model * 64]), tipX1, tipX2, 
		 0, width);		 
	  break;	
	case GAMMA:
	case GAMMA_I:
	  sumGAMMA(tipCase, &(tr->sumBuffer[offset]), 
		   x1_start, x2_start, 
		   &(tr->tipVectorDNA[model *64]), tipX1, tipX2,
		   0, width);
	  break;
	default:
	  assert(0);
	}
      break;
    case AA_DATA:
      switch(tr->rateHetModel)
	{
	case CAT:
	  sumGTRCATPROT(tipCase, &(tr->sumBuffer[offset]), x1_start, x2_start,  
			&(tr->tipVectorAA[model * 460]), 
			tipX1, tipX2, 0, width);		  	
	  break;
	case GAMMA:
	case GAMMA_I:
	  sumGAMMAPROT(tipCase,  &tr->sumBuffer[offset], x1_start, x2_start,  
		       &(tr->tipVectorAA[model * 460]), 
		       tipX1, tipX2, 0, width);
	  break;
	default:
	  assert(0);
	}
      break;
    default:
      assert(0);
    }	  
}     


#ifdef _LOCAL_DATA

void makenewzIterative(tree *localTree, int startIndex, int endIndex)
{  
  int pNumber, qNumber; 
       
  pNumber = localTree->td[0].ti[0].pNumber;
  qNumber = localTree->td[0].ti[0].qNumber;

  newviewIterative(localTree, startIndex, endIndex);
 
  if(localTree->mixedData)
    {
      int model;
      
      assert(0);

      for(model = 0; model < localTree->NumberOfModels; model++)
	makenewzMixedData(model, pNumber, qNumber, localTree);
    }
  else
    {
       double   
	 *x1_start = (double*)NULL, 
	 *x2_start = (double*)NULL;
       char 
	 *tipX1 = (char*)NULL,
	 *tipX2 = (char*)NULL;      
       int tipCase;
                
       if(isTip(pNumber, localTree->mxtips) || isTip(qNumber, localTree->mxtips))
	 {	    	    	   	               
	   if(!( isTip(pNumber, localTree->mxtips) && isTip(qNumber, localTree->mxtips)) )
	     {	 
	       tipCase = TIP_INNER;
	       if(isTip(qNumber, localTree->mxtips))
		 {	
		   tipX1 = &(localTree->strided_yVector[qNumber][startIndex]);	      
		   x2_start = getLikelihoodArray(pNumber, localTree->mxtips, localTree->xVector); 
		   
		 }	    
	       else
		 {
		   tipX1 = &(localTree->strided_yVector[pNumber][startIndex]);	     
		   x2_start = getLikelihoodArray(qNumber, localTree->mxtips, localTree->xVector); 
		 }
	     }
	   else
	     {
	       tipCase = TIP_TIP;
	       tipX1 = &(localTree->strided_yVector[pNumber][startIndex]);
	       tipX2 = &(localTree->strided_yVector[qNumber][startIndex]);
	     }
	 }
       else
	 {
	   tipCase = INNER_INNER;
	   
	   x1_start = getLikelihoodArray(pNumber, localTree->mxtips, localTree->xVector);
	   x2_start = getLikelihoodArray(qNumber, localTree->mxtips, localTree->xVector);
	 }
       
       switch(localTree->likelihoodFunction)
	{
	case GTRCAT:	 	 
	  sumCAT(tipCase, &(localTree->sumBuffer[startIndex * 4]), x1_start, x2_start, &(localTree->tipVectorDNA[0]), tipX1, tipX2, 
		 0, (endIndex - startIndex));
	  break;
	case GTRGAMMA:
	case GTRGAMMAI:
	  sumGAMMA(tipCase, &(localTree->sumBuffer[startIndex * 16]), x1_start, x2_start, &(localTree->tipVectorDNA[0]), tipX1, tipX2, 
		   0, (endIndex - startIndex));
	  break;
	case GTRGAMMAMULT:
	case GTRGAMMAMULTI:
	  if(localTree->multiBranch)
	    sumGAMMAMULT(tipCase, &(localTree->sumBuffer[startIndex * 16]), x1_start, x2_start, 
			 &(localTree->tipVectorDNA[0]), tipX1, tipX2,  &(localTree->strided_model[startIndex]),
			 0, (endIndex - startIndex));
	  else
	    sumGAMMAMULT(tipCase, &(localTree->sumBuffer[startIndex * 16]), x1_start, x2_start, 
			 &(localTree->tipVectorDNA[0]), tipX1, tipX2,  
			 &(localTree->strided_model[startIndex]), 0, (endIndex - startIndex));
	  break;
	case GTRCATMULT: 	 	  
	  if(localTree->multiBranch)
	    sumCATMULT(tipCase, &(localTree->sumBuffer[startIndex * 4]), x1_start, x2_start, &(localTree->tipVectorDNA[0]),
		       tipX1, tipX2, &(localTree->strided_model[startIndex]), 0, (endIndex - startIndex));
	  else
	    sumCATMULT(tipCase, &(localTree->sumBuffer[startIndex * 4]), x1_start, x2_start, &(localTree->tipVectorDNA[0]),
		       tipX1, tipX2, &(localTree->strided_model[startIndex]), startIndex, endIndex);
	  break;
	case PROTCAT:
	  sumGTRCATPROT(tipCase, &(localTree->sumBuffer[startIndex * 20]), x1_start, x2_start,  &(localTree->tipVectorAA[0]), 
			tipX1, tipX2, 0, (endIndex - startIndex));
	  break;
	case PROTCATMULT:
	  if(localTree->multiBranch)
	    sumGTRCATPROTMULT(tipCase, &(localTree->sumBuffer[startIndex * 20]), x1_start, x2_start,  &(localTree->tipVectorAA[0]), 
			      tipX1, tipX2, &(localTree->strided_model[startIndex]), 0, (endIndex - startIndex));
	  else
	    sumGTRCATPROTMULT(tipCase, &(localTree->sumBuffer[startIndex * 20]), x1_start, x2_start,  &(localTree->tipVectorAA[0]), 
			      tipX1, tipX2, &(localTree->strided_model[startIndex]), 0, (endIndex - startIndex));
	  break;
	case PROTGAMMA:
	case PROTGAMMAI:
	  sumGAMMAPROT(tipCase,  &(localTree->sumBuffer[startIndex * 80]), x1_start, x2_start,  &(localTree->tipVectorAA[0]), 
		       tipX1, tipX2, 0, (endIndex - startIndex));
	  break;
	case PROTGAMMAMULT:
	case PROTGAMMAMULTI:
	  if(localTree->multiBranch)
	    sumGAMMAPROTMULT(tipCase,  &(localTree->sumBuffer[startIndex * 80]), x1_start, x2_start,  &(localTree->tipVectorAA[0]), 
			     tipX1, tipX2, &(localTree->strided_model[startIndex]), 0, (endIndex - startIndex));
	  else
	    sumGAMMAPROTMULT(tipCase,  &(localTree->sumBuffer[startIndex * 80]), x1_start, x2_start,  &(localTree->tipVectorAA[0]), 
			     tipX1, tipX2, &(localTree->strided_model[startIndex]), 0, (endIndex - startIndex));
	  break;
	default:
	  assert(0);
	}
    }
}


#else

void makenewzIterative(tree *tr, int startIndex, int endIndex)
{  
  int pNumber, qNumber; 
       
  pNumber = tr->td[0].ti[0].pNumber;
  qNumber = tr->td[0].ti[0].qNumber;

  newviewIterative(tr, startIndex, endIndex);
 
  if(tr->mixedData)
    {
      int model;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	makenewzMixedData(model, pNumber, qNumber, tr);
    }
  else
    {
       double   
	 *x1_start = (double*)NULL, 
	 *x2_start = (double*)NULL;
       char 
	 *tipX1 = (char*)NULL,
	 *tipX2 = (char*)NULL;      
       int tipCase;
                
       if(isTip(pNumber, tr->rdta->numsp) || isTip(qNumber, tr->rdta->numsp))
	 {	    	    	   	               
	   if(!( isTip(pNumber, tr->rdta->numsp) && isTip(qNumber, tr->rdta->numsp)) )
	     {	 
	       tipCase = TIP_INNER;
	       if(isTip(qNumber, tr->rdta->numsp))
		 {	
		   tipX1 = tr->yVector[qNumber];	      
		   x2_start = getLikelihoodArray(pNumber,  tr->mxtips, tr->xVector);
		   
		 }	    
	       else
		 {
		   tipX1 = tr->yVector[pNumber];	     
		   x2_start = getLikelihoodArray(qNumber,  tr->mxtips, tr->xVector);
		 }
	     }
	   else
	     {
	       tipCase = TIP_TIP;
	       tipX1 = tr->yVector[pNumber];
	       tipX2 = tr->yVector[qNumber];
	     }
	 }
       else
	 {
	   tipCase = INNER_INNER;
	   
	   x1_start = getLikelihoodArray(pNumber, tr->mxtips, tr->xVector);
	   x2_start = getLikelihoodArray(qNumber, tr->mxtips, tr->xVector);
	 }
       
       switch(tr->likelihoodFunction)
	{
	case GTRCAT:
	  sumCAT(tipCase, tr->sumBuffer, x1_start, x2_start, &(tr->tipVectorDNA[0]), tipX1, tipX2, startIndex, 
		 endIndex);
	  break;
	case GTRGAMMA:
	case GTRGAMMAI:
	  sumGAMMA(tipCase, tr->sumBuffer, x1_start, x2_start, &(tr->tipVectorDNA[0]), tipX1, tipX2, 
		   startIndex, endIndex);
	  break;
	case GTRGAMMAMULT:
	case GTRGAMMAMULTI:
	  sumGAMMAMULT(tipCase, tr->sumBuffer, x1_start, x2_start, &(tr->tipVectorDNA[0]), tipX1, tipX2,  tr->model,
		       startIndex, endIndex);
	  break;
	case GTRCATMULT: 
	  sumCATMULT(tipCase, tr->sumBuffer, x1_start, x2_start, &(tr->tipVectorDNA[0]),
		     tipX1, tipX2, tr->model, startIndex, endIndex);
	  break;
	case PROTCAT:
	  sumGTRCATPROT(tipCase, tr->sumBuffer, x1_start, x2_start,  &(tr->tipVectorAA[0]), 
			tipX1, tipX2, startIndex, endIndex);
	  break;
	case PROTCATMULT:
	  sumGTRCATPROTMULT(tipCase, tr->sumBuffer, x1_start, x2_start,  &(tr->tipVectorAA[0]), 
			    tipX1, tipX2, tr->model, startIndex, endIndex);
	  break;
	case PROTGAMMA:
	case PROTGAMMAI:
	  sumGAMMAPROT(tipCase,  tr->sumBuffer, x1_start, x2_start,  &(tr->tipVectorAA[0]), 
		       tipX1, tipX2, startIndex, endIndex);
	  break;
	case PROTGAMMAMULT:
	case PROTGAMMAMULTI:
	  sumGAMMAPROTMULT(tipCase,  tr->sumBuffer, x1_start, x2_start,  &(tr->tipVectorAA[0]), 
			   tipX1, tipX2, tr->model, startIndex, endIndex);
	  break;
	default:
	  assert(0);
	}
    }
}

#endif


static void coreMixedData(tree *tr, int model, double *dlnLdlz, double *d2lnLdlz2, double lz)
{
  /*int l = tr->modelIndices[model][0];
    int u = tr->modelIndices[model][1];	  	*/

  int l = tr->partitionData[model].lower;
  int u = tr->partitionData[model].upper;

  int width  = u - l;
  /*int offset = tr->modelOffsets[model];*/
  int offset = tr->partitionData[model].modelOffset;

  switch(/*tr->dataType[model]*/ tr->partitionData[model].dataType)
    { 
    case DNA_DATA:
      switch(tr->rateHetModel)
	{
	case CAT:
	  coreGTRCAT(0, width, tr->NumberOfCategories, &(tr->sumBuffer[offset]), 
		     dlnLdlz, d2lnLdlz2, &(tr->cdta->wr[l]), &(tr->cdta->wr2[l]),
		     &(tr->cdta->patrat[0]), &(tr->EIGN_DNA[model * 3]),  &(tr->cdta->rateCategory[l]), lz);
	  break;
	case GAMMA:
	  coreGTRGAMMA(0, width, &tr->sumBuffer[offset], 
		       dlnLdlz, d2lnLdlz2, &(tr->EIGN_DNA[model * 3]), &(tr->gammaRates[model * 4]), lz, 
		       &(tr->cdta->aliaswgt[l]));
	  break;
	case GAMMA_I:
	  coreGTRGAMMAINVAR(&(tr->invariants[model]), &(tr->frequencies_DNA[model * 4]), 
			    &(tr->gammaRates[model * 4]), &(tr->EIGN_DNA[model * 3]),
			    &(tr->sumBuffer[offset]), dlnLdlz, d2lnLdlz2,
			    &(tr->invariant[l]), &(tr->cdta->aliaswgt[l]), 0, width, lz);
	  break;
	default:
	  assert(0);
	}
      break;
    case AA_DATA:
      switch(tr->rateHetModel)
	{
	case CAT: 
	  coreGTRCATPROT(&(tr->EIGN_AA[model * 19]), lz, tr->NumberOfCategories,  &(tr->cdta->patrat[0]), 
			 &(tr->cdta->rateCategory[l]), 0, width,
			 &(tr->cdta->wr[l]), &(tr->cdta->wr2[l]), dlnLdlz, d2lnLdlz2, 
			 &(tr->sumBuffer[offset]));		 		  	
	  break;	  
	case GAMMA:
	  coreGTRGAMMAPROT(&(tr->gammaRates[model * 4]), &(tr->EIGN_AA[model * 19]), 
			   &(tr->sumBuffer[offset]), 0, width, &(tr->cdta->aliaswgt[l]),
			   dlnLdlz, d2lnLdlz2, lz);
	  break;
	case GAMMA_I:
	  coreGTRGAMMAPROTINVAR(&(tr->gammaRates[model * 4]), &(tr->EIGN_AA[model * 19]), 
				&(tr->sumBuffer[offset]), 0, width, &(tr->cdta->aliaswgt[l]),
				dlnLdlz, d2lnLdlz2, lz, &(tr->frequencies_AA[model * 20]),
				&(tr->invariants[model]), &(tr->invariant[l]));
	  break;
	default:
	  assert(0);
	}
      break;
    default:
      assert(0);
    }		   
}

#ifdef _LOCAL_DATA

void execCore(tree *localTree, double *dlnLdlz, double *d2lnLdlz2, int lower, int upper, int model)
{
  double lz = localTree->coreLZ;

  if(localTree->mixedData)
    {         
      assert(0);
     
      if(localTree->multiBranch)	
	coreMixedData(localTree, model, dlnLdlz, d2lnLdlz2, lz);	      
      else
	{
	  int i;

	  double 
	    local_dlnLdlz, 
	    local_d2lnLdlz2;

	  *dlnLdlz   = 0.0;
	  *d2lnLdlz2 = 0.0;

	  for(i = 0; i < localTree->NumberOfModels; i++)
	    {
	      coreMixedData(localTree, i, &local_dlnLdlz, &local_d2lnLdlz2, lz);
	     
	      *dlnLdlz   += local_dlnLdlz;
	      *d2lnLdlz2 += local_d2lnLdlz2;
	    }
	}
    }
  else
    {
      switch(localTree->likelihoodFunction)
	{
	case GTRCAT:
	  coreGTRCAT(lower, upper, localTree->NumberOfCategories, localTree->sumBuffer, 
		     dlnLdlz, d2lnLdlz2, &(localTree->strided_wr[0]), &(localTree->strided_wr2[0]),
		     &(localTree->strided_patrat[0]), &(localTree->EIGN_DNA[0]),  &(localTree->strided_rateCategory[0]), lz);
	  break;
	case GTRCATMULT:      	
	  if(localTree->multiBranch)
	    coreGTRCAT(lower, upper, localTree->NumberOfCategories, localTree->sumBuffer, 
		       dlnLdlz, d2lnLdlz2, &(localTree->strided_wr[0]), &(localTree->strided_wr2[0]),
		       &(localTree->strided_patrat[0]), &(localTree->EIGN_DNA[model * 3]),  &(localTree->strided_rateCategory[0]), lz);
	  else
	    coreGTRCATMULT(lower, upper, localTree->NumberOfCategories, localTree->NumberOfModels, localTree->sumBuffer, 
			   dlnLdlz, d2lnLdlz2, &(localTree->strided_wr[0]), &(localTree->strided_wr2[0]),
			   &(localTree->strided_patrat[0]), localTree->EIGN_DNA,  &(localTree->strided_rateCategory[0]), 
			   localTree->strided_model, lz);
	  break;
	case PROTCAT: 
	  coreGTRCATPROT(localTree->EIGN_AA, lz, localTree->NumberOfCategories,  &(localTree->strided_patrat[0]), 
			 &(localTree->strided_rateCategory[0]), lower, upper,
			 &(localTree->strided_wr[0]), &(localTree->strided_wr2[0]), dlnLdlz, d2lnLdlz2, 
			 localTree->sumBuffer);
	  break;
	case PROTCATMULT:
	  if(localTree->multiBranch)
	    coreGTRCATPROT(&(localTree->EIGN_AA[model * 19]), lz, localTree->NumberOfCategories,  &(localTree->strided_patrat[0]), 
			   &(localTree->strided_rateCategory[0]), lower, upper,
			   &(localTree->strided_wr[0]), &(localTree->strided_wr2[0]), dlnLdlz, d2lnLdlz2, 
			   localTree->sumBuffer);
	  else
	    coreGTRCATPROTMULT(localTree->EIGN_AA, lz, localTree->NumberOfCategories,  &(localTree->strided_patrat[0]), 
			       &(localTree->strided_rateCategory[0]), lower, upper,
			       &(localTree->strided_wr[0]), &(localTree->strided_wr2[0]), dlnLdlz, d2lnLdlz2, 
			       localTree->sumBuffer, localTree->NumberOfModels, localTree->strided_model);
	  break;
	case GTRGAMMA:
	  coreGTRGAMMA(lower, upper, localTree->sumBuffer, 
		       dlnLdlz, d2lnLdlz2, localTree->EIGN_DNA, localTree->gammaRates, lz, 
		       localTree->strided_aliaswgt);
	  break;
	case GTRGAMMAI:
	  coreGTRGAMMAINVAR(localTree->invariants, localTree->frequencies_DNA, localTree->gammaRates, localTree->EIGN_DNA,
			    localTree->sumBuffer, dlnLdlz, d2lnLdlz2,
			    localTree->strided_invariant, localTree->strided_aliaswgt, lower, upper, lz);
	  break; 
	case GTRGAMMAMULT:
	  if(localTree->multiBranch)
	    coreGTRGAMMA(lower, upper, localTree->sumBuffer, 
			 dlnLdlz, d2lnLdlz2, &(localTree->EIGN_DNA[model * 3]), &(localTree->gammaRates[model * 4]), lz, 
			 localTree->strided_aliaswgt);
	  else
	    coreGTRGAMMAMULT(lower, upper, localTree->sumBuffer, 
			     dlnLdlz,  d2lnLdlz2, localTree->EIGN_DNA, localTree->gammaRates, 
			     lz, localTree->strided_aliaswgt,
			     localTree->NumberOfModels, localTree->strided_model);
	  
	  break;
	case GTRGAMMAMULTI:
	  if(localTree->multiBranch)
	    coreGTRGAMMAINVAR(&(localTree->invariants[model]), &(localTree->frequencies_DNA[model * 4]), 
			      &(localTree->gammaRates[model * 4]), &(localTree->EIGN_DNA[model * 3]),
			      localTree->sumBuffer, dlnLdlz, d2lnLdlz2,
			      localTree->strided_invariant, localTree->strided_aliaswgt, lower, upper, lz);
	  else
	    coreGTRGAMMAMULTINVAR(lower, upper, localTree->sumBuffer, 
				  dlnLdlz, d2lnLdlz2, localTree->EIGN_DNA, localTree->gammaRates, lz, 
				  localTree->strided_aliaswgt,
				  localTree->NumberOfModels, localTree->strided_model, localTree->frequencies_DNA, localTree->invariants, 
				  localTree->strided_invariant);
	  break;
	case PROTGAMMA:
	  coreGTRGAMMAPROT(localTree->gammaRates, localTree->EIGN_AA, localTree->sumBuffer, lower, upper, localTree->strided_aliaswgt,
			   dlnLdlz, d2lnLdlz2, lz);
	  break;
	case PROTGAMMAI:
	  coreGTRGAMMAPROTINVAR(localTree->gammaRates, localTree->EIGN_AA, localTree->sumBuffer, lower, upper, localTree->strided_aliaswgt,
				dlnLdlz, d2lnLdlz2, lz, localTree->frequencies_AA,
				localTree->invariants, localTree->strided_invariant);
	  break;
	case PROTGAMMAMULT:
	  if(localTree->multiBranch)
	    coreGTRGAMMAPROT(&(localTree->gammaRates[model * 4]), &(localTree->EIGN_AA[model * 19]), 
			     localTree->sumBuffer, lower, upper, localTree->strided_aliaswgt,
			     dlnLdlz, d2lnLdlz2, lz);
	  else
	    coreGTRGAMMAPROTMULT(localTree->gammaRates, localTree->EIGN_AA, localTree->sumBuffer, lower, upper, localTree->strided_aliaswgt,
				 dlnLdlz, d2lnLdlz2, lz, localTree->NumberOfModels, localTree->strided_model);
	  break;
	case PROTGAMMAMULTI:
	  if(localTree->multiBranch)
	    coreGTRGAMMAPROTINVAR(&(localTree->gammaRates[model * 4]), &(localTree->EIGN_AA[model * 19]), 
				  localTree->sumBuffer, lower, upper, localTree->strided_aliaswgt,
				  dlnLdlz, d2lnLdlz2, lz, &(localTree->frequencies_AA[model * 20]),
				  &(localTree->invariants[model]), localTree->strided_invariant);
	  else
	    coreGTRGAMMAPROTMULTINVAR(localTree->gammaRates, localTree->EIGN_AA, localTree->sumBuffer, lower, upper, 
				      localTree->strided_aliaswgt, dlnLdlz, d2lnLdlz2, lz, 
				      localTree->NumberOfModels, localTree->strided_model,
				      localTree->invariants, localTree->frequencies_AA, localTree->strided_invariant);      
	  break;
	default:
	  assert(0);
	}
    }

}

#else


void execCore(tree *tr, double *dlnLdlz, double *d2lnLdlz2, int lower, int upper, int model)
{
  double lz = tr->coreLZ;

  if(tr->mixedData)
    {              
      if(tr->multiBranch)	
	coreMixedData(tr, model, dlnLdlz, d2lnLdlz2, lz);	      
      else
	{
	  int i;

	  double 
	    local_dlnLdlz, 
	    local_d2lnLdlz2;

	  *dlnLdlz   = 0.0;
	  *d2lnLdlz2 = 0.0;

	  for(i = 0; i < tr->NumberOfModels; i++)
	    {
	      coreMixedData(tr, i, &local_dlnLdlz, &local_d2lnLdlz2, lz);
	     
	      *dlnLdlz   += local_dlnLdlz;
	      *d2lnLdlz2 += local_d2lnLdlz2;
	    }
	}
    }
  else
    {
      switch(tr->likelihoodFunction)
	{
	case GTRCAT:
	  coreGTRCAT(lower, upper, tr->NumberOfCategories, tr->sumBuffer, 
		     dlnLdlz, d2lnLdlz2, &(tr->cdta->wr[0]), &(tr->cdta->wr2[0]),
		     &(tr->cdta->patrat[0]), &(tr->EIGN_DNA[0]),  &(tr->cdta->rateCategory[0]), lz);
	  break;
	case GTRCATMULT:      
	  if(tr->multiBranch)
	    coreGTRCAT(lower, upper, tr->NumberOfCategories, tr->sumBuffer, 
		       dlnLdlz, d2lnLdlz2, &(tr->cdta->wr[0]), &(tr->cdta->wr2[0]),
		       &(tr->cdta->patrat[0]), &(tr->EIGN_DNA[model * 3]),  &(tr->cdta->rateCategory[0]), lz);
	  else
	    coreGTRCATMULT(lower, upper, tr->NumberOfCategories, tr->NumberOfModels, tr->sumBuffer, 
			   dlnLdlz, d2lnLdlz2, &(tr->cdta->wr[0]), &(tr->cdta->wr2[0]),
			   &(tr->cdta->patrat[0]), tr->EIGN_DNA,  &(tr->cdta->rateCategory[0]), 
			   tr->model, lz);
	  break;
	case PROTCAT: 
	  coreGTRCATPROT(tr->EIGN_AA, lz, tr->NumberOfCategories,  &(tr->cdta->patrat[0]), 
			 &(tr->cdta->rateCategory[0]), lower, upper,
			 &(tr->cdta->wr[0]), &(tr->cdta->wr2[0]), dlnLdlz, d2lnLdlz2, 
			 tr->sumBuffer);
	  break;
	case PROTCATMULT:
	  if(tr->multiBranch)
	    coreGTRCATPROT(&(tr->EIGN_AA[model * 19]), lz, tr->NumberOfCategories,  &(tr->cdta->patrat[0]), 
			   &(tr->cdta->rateCategory[0]), lower, upper,
			   &(tr->cdta->wr[0]), &(tr->cdta->wr2[0]), dlnLdlz, d2lnLdlz2, 
			   tr->sumBuffer);
	  else
	    coreGTRCATPROTMULT(tr->EIGN_AA, lz, tr->NumberOfCategories,  &(tr->cdta->patrat[0]), 
			       &(tr->cdta->rateCategory[0]), lower, upper,
			       &(tr->cdta->wr[0]), &(tr->cdta->wr2[0]), dlnLdlz, d2lnLdlz2, 
			       tr->sumBuffer, tr->NumberOfModels, tr->model);
	  break;
	case GTRGAMMA:
	  coreGTRGAMMA(lower, upper, tr->sumBuffer, 
		       dlnLdlz, d2lnLdlz2, tr->EIGN_DNA, tr->gammaRates, lz, 
		       tr->cdta->aliaswgt);
	  break;
	case GTRGAMMAI:
	  coreGTRGAMMAINVAR(tr->invariants, tr->frequencies_DNA, tr->gammaRates, tr->EIGN_DNA,
			    tr->sumBuffer, dlnLdlz, d2lnLdlz2,
			    tr->invariant, tr->cdta->aliaswgt, lower, upper, lz);
	  break; 
	case GTRGAMMAMULT:
	  if(tr->multiBranch)
	    coreGTRGAMMA(lower, upper, tr->sumBuffer, 
			 dlnLdlz, d2lnLdlz2, &(tr->EIGN_DNA[model * 3]), &(tr->gammaRates[model * 4]), lz, 
			 tr->cdta->aliaswgt);
	  else
	    coreGTRGAMMAMULT(lower, upper, tr->sumBuffer, 
			     dlnLdlz,  d2lnLdlz2, tr->EIGN_DNA, tr->gammaRates, 
			     lz, tr->cdta->aliaswgt,
			     tr->NumberOfModels, tr->model);
	  
	  break;
	case GTRGAMMAMULTI:
	  if(tr->multiBranch)
	    coreGTRGAMMAINVAR(&(tr->invariants[model]), &(tr->frequencies_DNA[model * 4]), 
			      &(tr->gammaRates[model * 4]), &(tr->EIGN_DNA[model * 3]),
			      tr->sumBuffer, dlnLdlz, d2lnLdlz2,
			      tr->invariant, tr->cdta->aliaswgt, lower, upper, lz);
	  else
	    coreGTRGAMMAMULTINVAR(lower, upper, tr->sumBuffer, 
				  dlnLdlz, d2lnLdlz2, tr->EIGN_DNA, tr->gammaRates, lz, 
				  tr->cdta->aliaswgt,
				  tr->NumberOfModels, tr->model, tr->frequencies_DNA, tr->invariants, 
				  tr->invariant);
	  break;
	case PROTGAMMA:
	  coreGTRGAMMAPROT(tr->gammaRates, tr->EIGN_AA, tr->sumBuffer, lower, upper, tr->cdta->aliaswgt,
			   dlnLdlz, d2lnLdlz2, lz);
	  break;
	case PROTGAMMAI:
	  coreGTRGAMMAPROTINVAR(tr->gammaRates, tr->EIGN_AA, tr->sumBuffer, lower, upper, tr->cdta->aliaswgt,
				dlnLdlz, d2lnLdlz2, lz, tr->frequencies_AA,
				tr->invariants, tr->invariant);
	  break;
	case PROTGAMMAMULT:
	  if(tr->multiBranch)
	    coreGTRGAMMAPROT(&(tr->gammaRates[model * 4]), &(tr->EIGN_AA[model * 19]), 
			     tr->sumBuffer, lower, upper, tr->cdta->aliaswgt,
			     dlnLdlz, d2lnLdlz2, lz);
	  else
	    coreGTRGAMMAPROTMULT(tr->gammaRates, tr->EIGN_AA, tr->sumBuffer, lower, upper, tr->cdta->aliaswgt,
				 dlnLdlz, d2lnLdlz2, lz, tr->NumberOfModels, tr->model);
	  break;
	case PROTGAMMAMULTI:
	  if(tr->multiBranch)
	    coreGTRGAMMAPROTINVAR(&(tr->gammaRates[model * 4]), &(tr->EIGN_AA[model * 19]), 
				  tr->sumBuffer, lower, upper, tr->cdta->aliaswgt,
				  dlnLdlz, d2lnLdlz2, lz, &(tr->frequencies_AA[model * 20]),
				  &(tr->invariants[model]), tr->invariant);
	  else
	    coreGTRGAMMAPROTMULTINVAR(tr->gammaRates, tr->EIGN_AA, tr->sumBuffer, lower, upper, 
				      tr->cdta->aliaswgt, dlnLdlz, d2lnLdlz2, lz, 
				      tr->NumberOfModels, tr->model,
				      tr->invariants, tr->frequencies_AA, tr->invariant);      
	  break;
	default:
	  assert(0);
	}
    }

}

#endif

static void topLevelMakenewz(tree *tr, int lower, int upper, int model, double z0, int maxiter, double *result)
{
  double   z, zprev, zstep;    
  double  dlnLdlz, d2lnLdlz2;
              
  z = z0;
  
  do 
    {
      int curvatOK = FALSE;
      
      zprev = z;
	  
      zstep = (1.0 - zmax) * z + zmin;
	  
      do 
	{		
	  double lz;	  	  	    	               
	    	    
	  if (z < zmin) z = zmin;
	  else if (z > zmax) z = zmax;
	  lz    = log(z);
	    	 	 	    	  
	  tr->coreLZ = lz;

	  
#ifdef _USE_PTHREADS
	  {
	    int i;

	    dlnLdlz   = 0.0;
	    d2lnLdlz2 = 0.0;
	    
	    tr->modelNumber = model;

	    masterBarrier(THREAD_MAKENEWZ, tr);

	    for(i = 0; i < NumberOfThreads; i++)
	      {
		dlnLdlz   += reductionBuffer[i];
		d2lnLdlz2 += reductionBufferTwo[i]; 
	      }
	  }
#else
	  execCore(tr, &dlnLdlz, &d2lnLdlz2, lower, upper, model);
#endif

	    	
	  if ((d2lnLdlz2 >= 0.0) && (z < zmax))
	    zprev = z = 0.37 * z + 0.63;  /*  Bad curvature, shorten branch */
	  else
	    curvatOK = TRUE;
	    
	} 
      while (! curvatOK);

      if (d2lnLdlz2 < 0.0) 
	{
	  double tantmp = -dlnLdlz / d2lnLdlz2;  /* prevent overflow */
	  if (tantmp < 100) 
	    {
	      z *= exp(tantmp);
	      if (z < zmin) 
		{
		  z = zmin;	    
		}
	      if (z > 0.25 * zprev + 0.75)    /*  Limit steps toward z = 1.0 */
		z = 0.25 * zprev + 0.75;
	    } 
	  else 
	    {
	      z = 0.25 * zprev + 0.75;
	    }
	  }
      if (z > zmax) z = zmax;
	  
    } 
  while ((--maxiter > 0) && (ABS(z - zprev) > zstep));	
   
  result[0] = z;
}
			     


static void switchAndAssign(tree *tr)
{
#ifdef _USE_PTHREADS        
  masterBarrier(THREAD_SUM_MAKENEWZ, tr);         
#else  
  makenewzIterative(tr, 0, tr->cdta->endsite);	
#endif

}


void makenewzGeneric(tree *tr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result)
{      
  int i, l, u, modelCounter;
  
  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;        
  for(i = 0; i < tr->numBranches; i++)      
    {
      tr->td[0].ti[0].qz[i] =  z0[i];      
    }

  tr->td[0].count = 1;


  if(!p->x)
    computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp, tr->numBranches);
  if(!q->x)
    computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp, tr->numBranches);        

  if(tr->mixedData)
    { 
      switch(tr->rateHetModel)
	{
	case CAT:
	  switchAndAssign(tr);
	  break;
	case GAMMA:
	case GAMMA_I:
	  switchAndAssign(tr);
	  break;
	default:
	  assert(0);
	}
    }
  else
    {
      switch(tr->likelihoodFunction)
	{
	case GTRCAT:
	  switchAndAssign(tr);                   
	  break;
	case GTRCATMULT: 
	  switchAndAssign(tr);           
	  break;
	case PROTCAT:
	  switchAndAssign(tr);             
	  break;
	case PROTCATMULT:
	  switchAndAssign(tr);        
	  break;
	case GTRGAMMA:
	  switchAndAssign(tr);         
	  break;
	case GTRGAMMAI:
	  switchAndAssign(tr);               
	  break;
	case GTRGAMMAMULT:
	  switchAndAssign(tr);            
	  break;
	case GTRGAMMAMULTI:
	  switchAndAssign(tr);     
	  break;
	case PROTGAMMA:
	  switchAndAssign(tr);     
	  break;
	case PROTGAMMAI:
	  switchAndAssign(tr);     
	  break;
	case PROTGAMMAMULT:
	  switchAndAssign(tr);      
	  break;
	case PROTGAMMAMULTI:
	  switchAndAssign(tr);     
	  break;
	default:
	  assert(0);
	}   
    }


  if(tr->multiBranch)         	  			
    for(modelCounter = 0; modelCounter < tr->NumberOfModels; modelCounter++)
      {
	l = tr->partitionData[modelCounter].lower;
	u = tr->partitionData[modelCounter].upper;

	
	topLevelMakenewz(tr, l, u, modelCounter, z0[modelCounter], maxiter, &(result[modelCounter]));
      }	  	     
  else
    topLevelMakenewz(tr, 0, tr->cdta->endsite, 0, z0[0], maxiter, &(result[0]));
}

void makenewzGenericDistance(tree *tr, int maxiter, double *z0, double *result, int taxon1, int taxon2)
{      
  int i, l, u, modelCounter;
  
  assert(taxon1 != taxon2);
  assert(0 < taxon1 && taxon1 <= tr->mxtips);
  assert(0 < taxon2 && taxon2 <= tr->mxtips);

  tr->td[0].ti[0].pNumber = taxon1;
  tr->td[0].ti[0].qNumber = taxon2;   
  tr->td[0].ti[0].tipCase = TIP_TIP;
   
  for(i = 0; i < tr->numBranches; i++)      
    {
      tr->td[0].ti[0].qz[i] =  defaultz;      
    }

  tr->td[0].count = 1;        

  if(tr->mixedData)
    { 
      switch(tr->rateHetModel)
	{	
	case GAMMA:
	case GAMMA_I:
	  switchAndAssign(tr);
	  break;
	default:
	  assert(0);
	}
    }
  else
    {
      switch(tr->likelihoodFunction)
	{	
	case GTRGAMMA:
	  switchAndAssign(tr);         
	  break;
	case GTRGAMMAI:
	  switchAndAssign(tr);               
	  break;
	case GTRGAMMAMULT:
	  switchAndAssign(tr);            
	  break;
	case GTRGAMMAMULTI:
	  switchAndAssign(tr);     
	  break;
	case PROTGAMMA:
	  switchAndAssign(tr);     
	  break;
	case PROTGAMMAI:
	  switchAndAssign(tr);     
	  break;
	case PROTGAMMAMULT:
	  switchAndAssign(tr);      
	  break;
	case PROTGAMMAMULTI:
	  switchAndAssign(tr);     
	  break;
	default:
	  assert(0);
	}   
    }


  if(tr->multiBranch)         	  			
    for(modelCounter = 0; modelCounter < tr->NumberOfModels; modelCounter++)
      {
	l = tr->partitionData[modelCounter].lower;
	u = tr->partitionData[modelCounter].upper;
	
	topLevelMakenewz(tr, l, u, modelCounter, z0[modelCounter], maxiter, &(result[modelCounter]));
      }	  	     
  else
    topLevelMakenewz(tr, 0, tr->cdta->endsite, 0, z0[0], maxiter, &(result[0]));
}


/********************************************************************************************************/

void execCorePartition(tree *tr, double *dlnLdlz, double *d2lnLdlz2, int lower, int upper, int model)
{
  double lz = tr->coreLZ;

  switch(tr->likelihoodFunction)
    {       
    case GTRGAMMA:  /* needed for rate opt*/
    case GTRGAMMAMULT:      
       coreGTRGAMMA(lower, upper, tr->sumBuffer, 
		    dlnLdlz, d2lnLdlz2, &(tr->EIGN_DNA[model * 3]), &(tr->gammaRates[model * 4]), lz, 
		    tr->cdta->aliaswgt);            
      break;
    case GTRGAMMAI:  /* needed for rate opt*/
    case GTRGAMMAMULTI:      
	coreGTRGAMMAINVAR(&(tr->invariants[model]), &(tr->frequencies_DNA[model * 4]), 
			  &(tr->gammaRates[model * 4]), &(tr->EIGN_DNA[model * 3]),
			  tr->sumBuffer, dlnLdlz, d2lnLdlz2,
			  tr->invariant, tr->cdta->aliaswgt, lower, upper, lz);      
      break;   
    case PROTGAMMA:  /* needed for rate opt*/
    case PROTGAMMAMULT:     
	coreGTRGAMMAPROT(&(tr->gammaRates[model * 4]), &(tr->EIGN_AA[model * 19]), 
			 tr->sumBuffer, lower, upper, tr->cdta->aliaswgt,
			 dlnLdlz, d2lnLdlz2, lz);      
      break;
    case PROTGAMMAI: /* needed for rate opt*/
    case PROTGAMMAMULTI:      
      coreGTRGAMMAPROTINVAR(&(tr->gammaRates[model * 4]), &(tr->EIGN_AA[model * 19]), 
			    tr->sumBuffer, lower, upper, tr->cdta->aliaswgt,
			    dlnLdlz, d2lnLdlz2, lz, &(tr->frequencies_AA[model * 20]),
			    &(tr->invariants[model]), tr->invariant);          
      break;
    default:
      assert(0);
    }

}


void makenewzIterativePartition(tree *tr, int startIndex, int endIndex, int model)
{
  double   
    *x1_start = (double*)NULL, 
    *x2_start = (double*)NULL;
   char 
    *tipX1 = (char*)NULL,
    *tipX2 = (char*)NULL;      
  int tipCase;
  int pNumber, qNumber; 

#ifdef _MULTI_GENE
  if(tr->doMulti)
    {
      pNumber = tr->td[model].ti[0].pNumber;
      qNumber = tr->td[model].ti[0].qNumber; 
    }
  else
#endif
    {
      pNumber = tr->td[0].ti[0].pNumber;
      qNumber = tr->td[0].ti[0].qNumber;  
    }

  newviewIterativePartition(tr, startIndex, endIndex, model);

  if(isTip(pNumber, tr->rdta->numsp) || isTip(qNumber, tr->rdta->numsp))
    {	    	    	   	               
      if(!(isTip(pNumber, tr->rdta->numsp) && isTip(qNumber, tr->rdta->numsp)))
	{	 
	  tipCase = TIP_INNER;
	  if(isTip(qNumber, tr->rdta->numsp))
	    {	
	      tipX1 = tr->yVector[qNumber];
	      x2_start =  getLikelihoodArray(pNumber, tr->mxtips, tr->xVector);
	    }	    
	  else
	    {
	      tipX1 = tr->yVector[pNumber];
	      x2_start = getLikelihoodArray(qNumber, tr->mxtips, tr->xVector);
	    }
	}
      else
	{
	  tipCase = TIP_TIP;
	  tipX1 = tr->yVector[pNumber];
	  tipX2 = tr->yVector[qNumber];
	}
    }
  else
    {
      tipCase = INNER_INNER;            

      x1_start = getLikelihoodArray(pNumber, tr->mxtips, tr->xVector);
      x2_start = getLikelihoodArray(qNumber, tr->mxtips, tr->xVector);
    }

  switch(tr->likelihoodFunction)
    {    
    case GTRGAMMA: /* needed for rate opt*/
    case GTRGAMMAI: /* needed for rate opt*/
    case GTRGAMMAMULT:
    case GTRGAMMAMULTI:
      sumGAMMAMULT(tipCase, tr->sumBuffer, x1_start, x2_start, &(tr->tipVectorDNA[0]), tipX1, tipX2,  tr->model,
		   startIndex, endIndex);
      break;    
    case PROTGAMMA: /* needed for rate opt*/
    case PROTGAMMAI:  /* needed for rate opt*/
    case PROTGAMMAMULT:
    case PROTGAMMAMULTI:
       sumGAMMAPROTMULT(tipCase,  tr->sumBuffer, x1_start, x2_start,  &(tr->tipVectorAA[0]), 
			tipX1, tipX2, tr->model, startIndex, endIndex);
      break;
    default:
      assert(0);
    }
}


static void switchAndAssignPartition(tree *tr, int lower, int upper, int model)
{
#ifdef _USE_PTHREADS    
  tr->modelNumber = model;
  masterBarrier(THREAD_SUM_MAKENEWZ_PARTITION, tr);         
#else
  makenewzIterativePartition(tr, lower, upper, model);	
#endif
}

static void topLevelMakenewzPartition(tree *tr, int lower, int upper, int model, double z0, int maxiter, double *result)
{
  double   z, zprev, zstep;    
  double  dlnLdlz, d2lnLdlz2;
              
  z = z0;
  

  do 
    {
      int curvatOK = FALSE;
      
      zprev = z;
	  
      zstep = (1.0 - zmax) * z + zmin;
	  
      do 
	{		
	  double lz;	  	  	    	               
	    	    
	  if (z < zmin) z = zmin;
	  else if (z > zmax) z = zmax;
	  lz    = log(z);
	    	 	 	    	  
	  tr->coreLZ = lz;

	  
#ifdef _USE_PTHREADS
	  {
	    int i;

	    dlnLdlz   = 0.0;
	    d2lnLdlz2 = 0.0;
	    
	    tr->modelNumber = model;

	    masterBarrier(THREAD_MAKENEWZ_PARTITION, tr);

	    for(i = 0; i < NumberOfThreads; i++)
	      {
		dlnLdlz   += reductionBuffer[i];
		d2lnLdlz2 += reductionBufferTwo[i]; 
	      }
	  }
#else
	  execCorePartition(tr, &dlnLdlz, &d2lnLdlz2, lower, upper, model);	  
#endif

	    	
	  if ((d2lnLdlz2 >= 0.0) && (z < zmax))
	    zprev = z = 0.37 * z + 0.63;  /*  Bad curvature, shorten branch */
	  else
	    curvatOK = TRUE;
	    
	} 
      while (! curvatOK);

      if (d2lnLdlz2 < 0.0) 
	{
	  double tantmp = -dlnLdlz / d2lnLdlz2;  /* prevent overflow */
	  if (tantmp < 100) 
	    {
	      z *= exp(tantmp);
	      if (z < zmin) 
		{
		  z = zmin;	    
		}
	      if (z > 0.25 * zprev + 0.75)    /*  Limit steps toward z = 1.0 */
		z = 0.25 * zprev + 0.75;
	    } 
	  else 
	    {
	      z = 0.25 * zprev + 0.75;
	    }
	  }
      if (z > zmax) z = zmax;
	  
    } 
  while ((--maxiter > 0) && (ABS(z - zprev) > zstep));	
   
  *result = z;
}
			  


double makenewzPartitionGeneric(tree *tr, nodeptr p, nodeptr q, double z0, int maxiter, int model)
{      
  int i, l, u;
  double result;
  
  assert(!tr->mixedData);  

  l = tr->partitionData[model].lower;
  u = tr->partitionData[model].upper;   

  assert(!tr->mixedData);
  /* should not be called when using mixed DNA/AA data */

#ifdef _MULTI_GENE
  if(tr->doMulti)
    {
      tr->td[model].ti[0].pNumber = p->number;
      tr->td[model].ti[0].qNumber = q->number;	  
      tr->td[model].ti[0].qz[model] =  q->z[model];	  
      tr->td[model].count = 1;	       

      if(!p->xs[model])
	computeMultiTraversalInfo(p, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->rdta->numsp, model);  
      if(!q->xs[model])
	computeMultiTraversalInfo(q, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->rdta->numsp, model);
      
      printf("%d\n", tr->td[model].count);

	/*
	  if(!p->xs[model])
	  computeFullMultiTraversalInfo(p, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->rdta->numsp, model);
	  if(!q->xs[model])
	  computeFullMultiTraversalInfo(q, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->rdta->numsp, model);
	*/

      /*printf("MAMA\n");*/

      if(!isTip(p->number, tr->mxtips))
	assert(p->xs[model]);
      
      if(!isTip(q->number, tr->mxtips))    
	assert(q->xs[model]);          
      
      /*printf("PAPA\n");*/
    }
  else
#endif
    {
      tr->td[0].ti[0].pNumber = p->number;
      tr->td[0].ti[0].qNumber = q->number;

      for(i = 0; i < tr->numBranches; i++)    
	tr->td[0].ti[0].qz[i] =  q->z[i]; 

      tr->td[0].count = 1;
      
      if(!p->x)
	computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp, tr->numBranches);
      if(!q->x)
	computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp, tr->numBranches);       
    }
 
  switchAndAssignPartition(tr, l, u, model);            
        	
  topLevelMakenewzPartition(tr, l, u, model, z0, maxiter, &result);
      
  /*free(tr->sumBuffer);*/

  return result;
}
