/*  RAxML-VI-HPC (version 2.2) a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright August 2006 by Alexandros Stamatakis
 *
 *  Partially derived from
 *  fastDNAml, a program for estimation of phylogenetic trees from sequences by Gary J. Olsen
 *  
 *  and 
 *
 *  Programs of the PHYLIP package by Joe Felsenstein.
 
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


extern int NumberOfThreads;

static void categorizeGTRGAMMA(char *tipX1, double *x2_start, int *ex2, double *EIGN, double *pz,
			       double *gammaRates, double *tipVector, int *rateCategory, int lower, int upper)
{
  double   z, lz, ki, term[4];    
  int     i, index = 0;
  double  *diagptable, *diagptable_start, allSum, min; 
  double *x1, *x2;
       
  z = pz[0];
  
  if (z < zmin) z = zmin;
  lz = log(z);
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 4 * 3);

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];	 
      
      *diagptable++ = exp (EIGN[0] * ki * lz);
      *diagptable++ = exp (EIGN[1] * ki * lz);
      *diagptable++ = exp (EIGN[2] * ki * lz);
    }          
        
  for (i = lower; i < upper; i++) 
    {
      x1 = &(tipVector[4 * tipX1[i]]);
      x2 = &x2_start[16 * i];
      diagptable = diagptable_start;
      
      /* cat 0 */
      
      term[0]  =  x1[0] * x2[0];
      term[0] += x1[1] * x2[1] * *diagptable++;
      term[0] += x1[2] * x2[2] * *diagptable++;
      term[0] += x1[3] * x2[3] * *diagptable++;     
      
      term[0] = log(term[0]) + (ex2[i])*log(minlikelihood);
      
      /* cat 1 */
      
      term[1]  = x1[0] * x2[4];
      term[1] += x1[1] * x2[5] * *diagptable++;
      term[1] += x1[2] * x2[6] * *diagptable++;
      term[1] += x1[3] * x2[7] * *diagptable++;    
      
      term[1] = log(term[1]) + (ex2[i])*log(minlikelihood);
      
      /* cat 2 */
      
      term[2]  = x1[0] * x2[8];
      term[2] += x1[1] * x2[9] * *diagptable++;
      term[2] += x1[2] * x2[10] * *diagptable++;
      term[2] += x1[3] * x2[11] * *diagptable++;     
      
      term[2] = log(term[2]) + (ex2[i])*log(minlikelihood);	 
      
      /* cat 3 */
      
      term[3]  = x1[0] * x2[12];
      term[3] += x1[1] * x2[13] * *diagptable++;
      term[3] += x1[2] * x2[14] * *diagptable++;
      term[3] += x1[3] * x2[15] * *diagptable++;     
      
      term[3] = log(term[3]) + (ex2[i])*log(minlikelihood);
      
      allSum = 0.25 * (term[0] + term[1] + term[2] + term[3]);
      
      allSum = 1.0 / allSum;
      
      term[0] = 0.25 * term[0] * allSum;
      term[1] = 0.25 * term[1] * allSum;
      term[2] = 0.25 * term[2] * allSum;
      term[3] = 0.25 * term[3] * allSum;
      
      min = largeDouble;
      if(term[0] < min)
	{
	  min   = term[0];
	  index = 0;
	}	
      if(term[1] < min)
	{
	  min   = term[1];
	  index = 1;
	}
      if(term[2] < min)
	{
	  min   = term[2];
	  index = 2;
	}
      if(term[3] < min)
	{
	  min   = term[3];
	  index = 3;
	}
      
      rateCategory[i] = index;              
    }
  
  free(diagptable_start);    
  return;
}

static void categorizeGTRGAMMAINVAR(char *tipX1, double *x2_start, int *ex2, double *pz, double *gammaRates, double *EIGN, 
				    double *tipVector, int *rateCategory, int lower, int upper, 
				    int *iptr, double *invariants, double *frequencies)
{
  double   z, lz, ki, term[4];      
  int     i, index = 0;
  double  *diagptable, *diagptable_start, allSum, min;
  double *x1, *x2; 
  double freqs[4];
  double scaler = (1.0 - invariants[0]);

  freqs[0] = frequencies[0] * invariants[0]; 
  freqs[1] = frequencies[1] * invariants[0];
  freqs[2] = frequencies[2] * invariants[0];
  freqs[3] = frequencies[3] * invariants[0];    
   
 
  z = pz[0];
  
  if (z < zmin) z = zmin;
  lz = log(z);
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 4 * 3);

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];	 
      
      *diagptable++ = exp (EIGN[0] * ki * lz);
      *diagptable++ = exp (EIGN[1] * ki * lz);
      *diagptable++ = exp (EIGN[2] * ki * lz);
    }      


 
  for (i = lower; i < upper; i++) 
    {
      x1 = &(tipVector[4 * tipX1[i]]);
      x2 = &x2_start[16 * i];
      diagptable = diagptable_start;
      
      /* cat 0 */
      
      term[0]  =  x1[0] * x2[0];
      term[0] += x1[1] * x2[1] * *diagptable++;
      term[0] += x1[2] * x2[2] * *diagptable++;
      term[0] += x1[3] * x2[3] * *diagptable++;     

      if(iptr[i] < 4)
	term[0] = log((scaler * term[0] + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
      else
	term[0] = log(scaler * term[0]) + (ex2[i]) * log(minlikelihood);
      
      /* cat 1 */
      
      term[1]  = x1[0] * x2[4];
      term[1] += x1[1] * x2[5] * *diagptable++;
      term[1] += x1[2] * x2[6] * *diagptable++;
      term[1] += x1[3] * x2[7] * *diagptable++;             
      
      if(iptr[i] < 4)
	term[1] = log((scaler * term[1] + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
      else
	term[1] = log(scaler * term[1]) + (ex2[i]) * log(minlikelihood);

      /* cat 2 */
      
      term[2]  = x1[0] * x2[8];
      term[2] += x1[1] * x2[9] * *diagptable++;
      term[2] += x1[2] * x2[10] * *diagptable++;
      term[2] += x1[3] * x2[11] * *diagptable++;     
          	 
      
      if(iptr[i] < 4)
	term[2] = log((scaler * term[2] + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
      else
	term[2] = log(scaler * term[2]) + (ex2[i]) * log(minlikelihood);

      /* cat 3 */
      
      term[3]  = x1[0] * x2[12];
      term[3] += x1[1] * x2[13] * *diagptable++;
      term[3] += x1[2] * x2[14] * *diagptable++;
      term[3] += x1[3] * x2[15] * *diagptable++;     
      
      if(iptr[i] < 4)
	term[3] = log((scaler * term[3] + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
      else
	term[3] = log(scaler * term[3]) + (ex2[i]) * log(minlikelihood);
      
      allSum = 0.25 * (term[0] + term[1] + term[2] + term[3]);
      
      allSum = 1.0 / allSum;
      
      term[0] = 0.25 * term[0] * allSum;
      term[1] = 0.25 * term[1] * allSum;
      term[2] = 0.25 * term[2] * allSum;
      term[3] = 0.25 * term[3] * allSum;
      
      min = largeDouble;
      if(term[0] < min)
	{
	  min   = term[0];
	  index = 0;
	}	
      if(term[1] < min)
	{
	  min   = term[1];
	  index = 1;
	}
      if(term[2] < min)
	{
	  min   = term[2];
	  index = 2;
	}
      if(term[3] < min)
	{
	  min   = term[3];
	  index = 3;
	}
      
      rateCategory[i] = index;     
    }
  
  free(diagptable_start);    
  return;
}

static void categorizeGTRGAMMAMULT (char *tipX1, double *x2_start, int *ex2, double *pz, double *gammaRates, double *EIGN, 
				    double *tipVector, int *rateCategory, int lower, int upper, int *modelptr, int NumberOfModels,
				    int multiBranch) 
{
  double   z, lz = 0.0, term[4], ki;     
  int     i, index = 0;  
  double  *diagptable, *diagptable_start, allSum, min;
  double *x1, *x2;
  int model;
    
  if(!multiBranch)
    {
      z = pz[0]; 
      if (z < zmin) z = zmin;
      lz = log(z);
    }
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 12 * NumberOfModels);

  for(model = 0; model < NumberOfModels; model++)
    {
      if(multiBranch)
	{
	  z = pz[model]; 
	  if (z < zmin) z = zmin;
	  lz = log(z);
	}
      
      diagptable = &diagptable_start[12 * model];
      
      for(i = 0; i < 4; i++)
	{
	  ki = gammaRates[model * 4 + i];	 
	  
	  *diagptable++ = exp (EIGN[model * 3]     * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 1] * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 2] * ki * lz);
	}     
    }
  

  for (i = lower; i < upper; i++) 
    {
      model = modelptr[i];

      x1 = &(tipVector[64 * model + 4 * tipX1[i]]);	 
      x2 = &x2_start[16 * i];
	   
      diagptable = &diagptable_start[12 * model];
	  
      /* cat 0 */
      
      term[0] =  x1[0] * x2[0];	 
      term[0] += x1[1] * x2[1] * *diagptable++;	
      term[0] += x1[2] * x2[2] * *diagptable++;	  
      term[0] += x1[3] * x2[3] * *diagptable++;     

      term[0] = log(term[0]) + ex2[i] *log(minlikelihood);
	
      /* cat 1 */
	  
      term[1] = x1[0] * x2[4];
      term[1] += x1[1] * x2[5] * *diagptable++;
      term[1] += x1[2] * x2[6] * *diagptable++;
      term[1] += x1[3] * x2[7] * *diagptable++;     

      term[1] = log(term[1]) + ex2[i] * log(minlikelihood);
	  
      /* cat 2 */
	  
      term[2]  = x1[0] * x2[8];
      term[2] += x1[1] * x2[9] * *diagptable++;
      term[2] += x1[2] * x2[10] * *diagptable++;
      term[2] += x1[3] * x2[11] * *diagptable++;     
	  
      term[2] = log(term[2]) + ex2[i] *log(minlikelihood);

      /* cat 3 */
	  
      term[3]  = x1[0] * x2[12];
      term[3] += x1[1] * x2[13] * *diagptable++;
      term[3] += x1[2] * x2[14] * *diagptable++;
      term[3] += x1[3] * x2[15] * *diagptable++;     
	  
      term[3] = log(term[3]) + (ex2[i])*log(minlikelihood);

      allSum = 0.25 * (term[0] + term[1] + term[2] + term[3]);
      
      allSum = 1.0 / allSum;
      
      term[0] = 0.25 * term[0] * allSum;
      term[1] = 0.25 * term[1] * allSum;
      term[2] = 0.25 * term[2] * allSum;
      term[3] = 0.25 * term[3] * allSum;
      
      min = largeDouble;
      if(term[0] < min)
	{
	  min   = term[0];
	  index = 0;
	}	
      if(term[1] < min)
	{
	  min   = term[1];
	  index = 1;
	}
      if(term[2] < min)
	{
	  min   = term[2];
	  index = 2;
	}
      if(term[3] < min)
	{
	  min   = term[3];
	  index = 3;
	}
      
      rateCategory[i] = model * 4 + index;    
    }
  
  free(diagptable_start);       
  
  return;
}
  
static void categorizeGTRGAMMAMULTINVAR (char *tipX1, double *x2_start, int *ex2, double *pz, double *gammaRates, double *EIGN,
					 double *tipVector, int *rateCategory, int lower, int upper, int *modelptr, int NumberOfModels,
					 int *iptr, double *invariants, double *frequencies, int multiBranch) 
{
  double   z, lz = 0.0, term[4], ki;    
  int     i, index = 0;  
  double  *diagptable, *diagptable_start, allSum, min;
  double *x1, *x2;
  int model;
  double *scalers;
  double *freqs;
   
  if(!multiBranch)
    {
      z = pz[0]; 
      if (z < zmin) z = zmin;
      lz = log(z);
    }
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 12 * NumberOfModels);
  freqs      = (double *)malloc(4 * NumberOfModels * sizeof(double));
  scalers    = (double *)malloc(NumberOfModels * sizeof(double)); 

  for(model = 0; model < NumberOfModels; model++)
    {
       if(multiBranch)
	 {
	   z = pz[model]; 
	   if (z < zmin) z = zmin;
	   lz = log(z);
	 }

      scalers[model] = (1.0 - invariants[model]); 
      freqs[4 * model]     = frequencies[4 * model]     * invariants[model]; 
      freqs[4 * model + 1] = frequencies[4 * model + 1] * invariants[model];
      freqs[4 * model + 2] = frequencies[4 * model + 2] * invariants[model];
      freqs[4 * model + 3] = frequencies[4 * model + 3] * invariants[model];
      
      diagptable = &diagptable_start[12 * model];
      
      for(i = 0; i < 4; i++)
	{
	  ki = gammaRates[model * 4 + i];	 
	  
	  *diagptable++ = exp (EIGN[model * 3] * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 1] * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 2] * ki * lz);
	}

     
    }
  

  for (i = lower; i < upper; i++) 
    {
      model = modelptr[i];

      x1 = &(tipVector[64 * model + 4 * tipX1[i]]);	 
      x2 = &x2_start[16 * i];
	   
      diagptable = &diagptable_start[12 * model];
	  
      /* cat 0 */
      
      term[0] =  x1[0] * x2[0];	 
      term[0] += x1[1] * x2[1] * *diagptable++;	
      term[0] += x1[2] * x2[2] * *diagptable++;	  
      term[0] += x1[3] * x2[3] * *diagptable++;     

       if(iptr[i] < 4)	  
	 term[0] = log((scalers[model] * term[0] + freqs[model * 4 + iptr[i]]) * pow(minlikelihood, ex2[i]));
       else
	 term[0] = log(scalers[model] * term[0]) + ex2[i] * log(minlikelihood);
     	
      /* cat 1 */
	  
      term[1] = x1[0] * x2[4];
      term[1] += x1[1] * x2[5] * *diagptable++;
      term[1] += x1[2] * x2[6] * *diagptable++;
      term[1] += x1[3] * x2[7] * *diagptable++;           
	        
      if(iptr[i] < 4)	  
	term[1] = log((scalers[model] * term[1] + freqs[model * 4 + iptr[i]]) * pow(minlikelihood, ex2[i]));
      else
	term[1] = log(scalers[model] * term[1]) + ex2[i] * log(minlikelihood);

      /* cat 2 */
	  
      term[2]  = x1[0] * x2[8];
      term[2] += x1[1] * x2[9] * *diagptable++;
      term[2] += x1[2] * x2[10] * *diagptable++;
      term[2] += x1[3] * x2[11] * *diagptable++;     	       
      
      if(iptr[i] < 4)	  
	term[2] = log((scalers[model] * term[2] + freqs[model * 4 + iptr[i]]) * pow(minlikelihood, ex2[i]));
      else
	term[2] = log(scalers[model] * term[2]) + ex2[i] * log(minlikelihood);

      /* cat 3 */
	  
      term[3]  = x1[0] * x2[12];
      term[3] += x1[1] * x2[13] * *diagptable++;
      term[3] += x1[2] * x2[14] * *diagptable++;
      term[3] += x1[3] * x2[15] * *diagptable++;     
	        
      if(iptr[i] < 4)	  
	term[3] = log((scalers[model] * term[3] + freqs[model * 4 + iptr[i]]) * pow(minlikelihood, ex2[i]));
      else
	term[3] = log(scalers[model] * term[3]) + ex2[i] * log(minlikelihood);      

      allSum = 0.25 * (term[0] + term[1] + term[2] + term[3]);
      
      allSum = 1.0 / allSum;
      
      term[0] = 0.25 * term[0] * allSum;
      term[1] = 0.25 * term[1] * allSum;
      term[2] = 0.25 * term[2] * allSum;
      term[3] = 0.25 * term[3] * allSum;
      
      min = largeDouble;
      if(term[0] < min)
	{
	  min   = term[0];
	  index = 0;
	}	
      if(term[1] < min)
	{
	  min   = term[1];
	  index = 1;
	}
      if(term[2] < min)
	{
	  min   = term[2];
	  index = 2;
	}
      if(term[3] < min)
	{
	  min   = term[3];
	  index = 3;
	}
      
      rateCategory[i] = model * 4 + index;
      
    }
  
  free(diagptable_start);       
  free(scalers);
  free(freqs);

  return;
}
  
static void categorizeGTRGAMMAPROT (char *tipX1, double *x2, int *ex2, double *pz, double *gammaRates, double *EIGN, 
				    double *tipVector, int *rateCategory, int lower, int upper)
{
  double   z, lz, term[4], ki, min, allSum;    
  int     i, j, index = 0;     
  double  *diagptable, *diagptable_start;
  double  *left, *right;           
    
  z = pz[0];
   
  if (z < zmin) z = zmin;
  lz = log(z);
       
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 4 * 19);

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];	 

      for(j = 0; j < 19; j++)
	*diagptable++ = exp (EIGN[j] * ki * lz);	
    }
	       
  for (i = lower; i < upper; i++) 
    {
      left = &(tipVector[20 * tipX1[i]]);
      diagptable = diagptable_start;
	    	 	
      for(j = 0; j < 4; j++)
	{
	  right = &(x2[80 * i + 20 * j]);
	  term[j]  =  0.0;
	  term[j] +=  left[0] * right[0];
	  term[j] +=  left[1] * right[1] * *diagptable++;
	  term[j] +=  left[2] * right[2] * *diagptable++;
	  term[j] +=  left[3] * right[3] * *diagptable++;
	  term[j] +=  left[4] * right[4] * *diagptable++;
	  term[j] +=  left[5] * right[5] * *diagptable++;
	  term[j] +=  left[6] * right[6] * *diagptable++;	
	  term[j] +=  left[7] * right[7] * *diagptable++;
	  term[j] +=  left[8] * right[8] * *diagptable++;
	  term[j] +=  left[9] * right[9] * *diagptable++;
	  term[j] +=  left[10] * right[10] * *diagptable++;
	  term[j] +=  left[11] * right[11] * *diagptable++;
	  term[j] +=  left[12] * right[12] * *diagptable++;
	  term[j] +=  left[13] * right[13] * *diagptable++;
	  term[j] +=  left[14] * right[14] * *diagptable++;
	  term[j] +=  left[15] * right[15] * *diagptable++;
	  term[j] +=  left[16] * right[16] * *diagptable++;
	  term[j] +=  left[17] * right[17] * *diagptable++;
	  term[j] +=  left[18] * right[18] * *diagptable++;	
	  term[j] +=  left[19] * right[19] * *diagptable++;
	  
	  term[j] = log(term[j]) +  ex2[i] * log(minlikelihood);
	}	  
	    
      allSum = 0.25 * (term[0] + term[1] + term[2] + term[3]);
      
      allSum = 1.0 / allSum;
	
      term[0] = 0.25 * term[0] * allSum;
      term[1] = 0.25 * term[1] * allSum;
      term[2] = 0.25 * term[2] * allSum;
      term[3] = 0.25 * term[3] * allSum;
    
      min = largeDouble;
      
      if(term[0] < min)
	{
	  min   = term[0];
	  index = 0;
	}	
      if(term[1] < min)
	{
	  min   = term[1];
	  index = 1;
	}
      if(term[2] < min)
	{
	  min   = term[2];
	  index = 2;
	}
      if(term[3] < min)
	{
	  min   = term[3];
	  index = 3;
	}
	
      rateCategory[i] = index;
    }
  
  free(diagptable_start); 
  return;
}	

static void categorizeGTRGAMMAPROTINVAR (char *tipX1, double *x2, int *ex2, double *pz, double *gammaRates, double *EIGN, 
					 double *tipVector, int *rateCategory, int lower, int upper,
					 int *iptr, double *invariants, double *frequencies)
  {
    double   z, lz, term[4], ki, min, allSum;    
    int     i, j, index = 0;     
    double  *diagptable, *diagptable_start;
    double  *left, *right;   
    double scaler = (1.0 - invariants[0]);
    double freqs[20];   
    
    for(i = 0; i < 20; i++)
      freqs[i] = frequencies[i] * invariants[0];
            
    z = pz[0];
   
    if (z < zmin) z = zmin;
    lz = log(z);
       
    diagptable = diagptable_start = (double *)malloc(sizeof(double) * 4 * 19);

    for(i = 0; i < 4; i++)
      {
	ki = gammaRates[i];	 

	for(j = 0; j < 19; j++)
	  *diagptable++ = exp (EIGN[j] * ki * lz);	
      }
	    
  
    for (i = lower; i < upper; i++) 
      {
	left = &(tipVector[20 * tipX1[i]]);
	diagptable = diagptable_start;
	    	 	
	for(j = 0; j < 4; j++)
	  {
	    right = &(x2[80 * i + 20 * j]);
	    term[j] = 0;
	    term[j] +=  left[0] * right[0];
	    term[j] +=  left[1] * right[1] * *diagptable++;
	    term[j] +=  left[2] * right[2] * *diagptable++;
	    term[j] +=  left[3] * right[3] * *diagptable++;
	    term[j] +=  left[4] * right[4] * *diagptable++;
	    term[j] +=  left[5] * right[5] * *diagptable++;
	    term[j] +=  left[6] * right[6] * *diagptable++;	
	    term[j] +=  left[7] * right[7] * *diagptable++;
	    term[j] +=  left[8] * right[8] * *diagptable++;
	    term[j] +=  left[9] * right[9] * *diagptable++;
	    term[j] +=  left[10] * right[10] * *diagptable++;
	    term[j] +=  left[11] * right[11] * *diagptable++;
	    term[j] +=  left[12] * right[12] * *diagptable++;
	    term[j] +=  left[13] * right[13] * *diagptable++;
	    term[j] +=  left[14] * right[14] * *diagptable++;
	    term[j] +=  left[15] * right[15] * *diagptable++;
	    term[j] +=  left[16] * right[16] * *diagptable++;
	    term[j] +=  left[17] * right[17] * *diagptable++;
	    term[j] +=  left[18] * right[18] * *diagptable++;	
	    term[j] +=  left[19] * right[19] * *diagptable++;

	    if(iptr[i] < 20)	   
	      term[j] = log((scaler * term[j] + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term[j] = log(scaler * term[j]) + ex2[i] * log(minlikelihood);	 
	  }	  
	    
	allSum = 0.25 * (term[0] + term[1] + term[2] + term[3]);
      
	allSum = 1.0 / allSum;
	
	term[0] = 0.25 * term[0] * allSum;
	term[1] = 0.25 * term[1] * allSum;
	term[2] = 0.25 * term[2] * allSum;
	term[3] = 0.25 * term[3] * allSum;
    
	min = largeDouble;
	if(term[0] < min)
	  {
	    min   = term[0];
	    index = 0;
	  }	
	if(term[1] < min)
	  {
	    min   = term[1];
	    index = 1;
	  }
	if(term[2] < min)
	  {
	    min   = term[2];
	    index = 2;
	  }
	if(term[3] < min)
	  {
	    min   = term[3];
	    index = 3;
	  }
	
	rateCategory[i] = index;
      }
    
    free(diagptable_start); 
    return;
  }
	
static void categorizeGTRGAMMAPROTMULT (char *tipX1, double *x2, int *ex2, double *pz, double *gammaRates, double *EIGN, 
					double *tipVector, int *rateCategory, int lower, int upper, int *modelptr, int NumberOfModels, 
					int multiBranch)
  {
    double   z, lz = 0.0, term[4], ki, min, allSum;    
    int     i, j, index = 0;     
    double  *diagptable, *diagptable_start;
    double  *left, *right;  
    int model;
            
    if(!multiBranch)
      {
	z = pz[0];   
	if (z < zmin) z = zmin;
	lz = log(z);
      }
       
    diagptable = diagptable_start = (double *)malloc(sizeof(double) * 4 * 19 * NumberOfModels);

    for(model = 0; model < NumberOfModels; model++)
      {
	if(multiBranch)
	  {
	    z = pz[model];   
	    if (z < zmin) z = zmin;
	    lz = log(z);
	  }
       
	diagptable = &diagptable_start[4 * 19 * model];
	
	for(i = 0; i < 4; i++)
	  {
	    ki = gammaRates[model * 4 + i];	 
	    
	    for(j = 0; j < 19; j++)
	      *diagptable++ = exp (EIGN[model * 19 + j] * ki * lz);	
	  }
	

    }
    	            	   
    for (i = lower; i < upper; i++) 
      {
	model = modelptr[i];

	left = &(tipVector[460 * model + 20 * tipX1[i]]);
	diagptable = &diagptable_start[76 * model];
	    	 	
	for(j = 0; j < 4; j++)
	  {
	    right = &(x2[80 * i + 20 * j]);
	    term[j] = 0;
	    term[j] +=  left[0] * right[0];
	    term[j] +=  left[1] * right[1] * *diagptable++;
	    term[j] +=  left[2] * right[2] * *diagptable++;
	    term[j] +=  left[3] * right[3] * *diagptable++;
	    term[j] +=  left[4] * right[4] * *diagptable++;
	    term[j] +=  left[5] * right[5] * *diagptable++;
	    term[j] +=  left[6] * right[6] * *diagptable++;	
	    term[j] +=  left[7] * right[7] * *diagptable++;
	    term[j] +=  left[8] * right[8] * *diagptable++;
	    term[j] +=  left[9] * right[9] * *diagptable++;
	    term[j] +=  left[10] * right[10] * *diagptable++;
	    term[j] +=  left[11] * right[11] * *diagptable++;
	    term[j] +=  left[12] * right[12] * *diagptable++;
	    term[j] +=  left[13] * right[13] * *diagptable++;
	    term[j] +=  left[14] * right[14] * *diagptable++;
	    term[j] +=  left[15] * right[15] * *diagptable++;
	    term[j] +=  left[16] * right[16] * *diagptable++;
	    term[j] +=  left[17] * right[17] * *diagptable++;
	    term[j] +=  left[18] * right[18] * *diagptable++;	
	    term[j] +=  left[19] * right[19] * *diagptable++;

	    term[j] = log(term[j]) +  ex2[i] * log(minlikelihood);
	  }	  
	    
	allSum = 0.25 * (term[0] + term[1] + term[2] + term[3]);
      
	allSum = 1.0 / allSum;
	
	term[0] = 0.25 * term[0] * allSum;
	term[1] = 0.25 * term[1] * allSum;
	term[2] = 0.25 * term[2] * allSum;
	term[3] = 0.25 * term[3] * allSum;
    
	min = largeDouble;
	if(term[0] < min)
	  {
	    min   = term[0];
	    index = 0;
	  }	
	if(term[1] < min)
	  {
	    min   = term[1];
	    index = 1;
	  }
	if(term[2] < min)
	  {
	    min   = term[2];
	    index = 2;
	  }
	if(term[3] < min)
	  {
	    min   = term[3];
	    index = 3;
	  }
	
	rateCategory[i] = index;
      }
    
    free(diagptable_start); 
    return;
  }
	

static void categorizeGTRGAMMAPROTMULTINVAR (char *tipX1, double *x2, int *ex2, double *pz, double *gammaRates, double *EIGN,
					     double *tipVector, int *rateCategory, int lower, int upper, int *modelptr, 
					     int NumberOfModels, int *iptr, double *invariants, double *frequencies, int multiBranch)
  {
    double   z, lz = 0.0, term[4], ki, min, allSum;    
    int     i, j, index = 0;     
    double  *diagptable, *diagptable_start;
    double  *left, *right;   
    int model;
    double *scalers;
    double *freqs;   
            
    if(!multiBranch)
      {
	z = pz[0];   
	if (z < zmin) z = zmin;
	lz = log(z);
      }
       
    diagptable = diagptable_start = (double *)malloc(sizeof(double) * 4 * 19 * NumberOfModels);
    freqs      = (double *)malloc(20 * NumberOfModels * sizeof(double));
    scalers    = (double *)malloc(NumberOfModels * sizeof(double));

    for(model = 0; model < NumberOfModels; model++)
      {
	if(multiBranch)
	  {
	    z = pz[model];   
	    if (z < zmin) z = zmin;
	    lz = log(z);
	  }
       
	diagptable = &diagptable_start[4 * 19 * model];
	scalers[model] = (1.0 - invariants[model]); 

	for(i = 0; i < 20; i++)
	  freqs[20 * model + i]     = frequencies[20 * model + i] * invariants[model];
	
	for(i = 0; i < 4; i++)
	  {
	    ki = gammaRates[model * 4 + i];	 
	    
	    for(j = 0; j < 19; j++)
	      *diagptable++ = exp (EIGN[model * 19 + j] * ki * lz);	
	  }	
    }
 
    for (i = lower; i < upper; i++) 
      {
	model = modelptr[i];

	left = &(tipVector[460 * model + 20 * tipX1[i]]);
	diagptable = &diagptable_start[76 * model];
	    	 	
	for(j = 0; j < 4; j++)
	  {
	    right = &(x2[80 * i + 20 * j]);
	    term[j] = 0;
	    term[j] +=  left[0] * right[0];
	    term[j] +=  left[1] * right[1] * *diagptable++;
	    term[j] +=  left[2] * right[2] * *diagptable++;
	    term[j] +=  left[3] * right[3] * *diagptable++;
	    term[j] +=  left[4] * right[4] * *diagptable++;
	    term[j] +=  left[5] * right[5] * *diagptable++;
	    term[j] +=  left[6] * right[6] * *diagptable++;	
	    term[j] +=  left[7] * right[7] * *diagptable++;
	    term[j] +=  left[8] * right[8] * *diagptable++;
	    term[j] +=  left[9] * right[9] * *diagptable++;
	    term[j] +=  left[10] * right[10] * *diagptable++;
	    term[j] +=  left[11] * right[11] * *diagptable++;
	    term[j] +=  left[12] * right[12] * *diagptable++;
	    term[j] +=  left[13] * right[13] * *diagptable++;
	    term[j] +=  left[14] * right[14] * *diagptable++;
	    term[j] +=  left[15] * right[15] * *diagptable++;
	    term[j] +=  left[16] * right[16] * *diagptable++;
	    term[j] +=  left[17] * right[17] * *diagptable++;
	    term[j] +=  left[18] * right[18] * *diagptable++;	
	    term[j] +=  left[19] * right[19] * *diagptable++;

	    if(iptr[i] < 20)	  
	      term[j] = log((scalers[model] * term[j] + freqs[model * 20 + iptr[i]]) * pow(minlikelihood, ex2[i]));
	    else
	      term[j] = log(scalers[model] * term[j]) + ex2[i] * log(minlikelihood);	   
	  }
	    
	allSum = 0.25 * (term[0] + term[1] + term[2] + term[3]);
      
	allSum = 1.0 / allSum;
	
	term[0] = 0.25 * term[0] * allSum;
	term[1] = 0.25 * term[1] * allSum;
	term[2] = 0.25 * term[2] * allSum;
	term[3] = 0.25 * term[3] * allSum;
    
	min = largeDouble;
	if(term[0] < min)
	  {
	    min   = term[0];
	    index = 0;
	  }	
	if(term[1] < min)
	  {
	    min   = term[1];
	    index = 1;
	  }
	if(term[2] < min)
	  {
	    min   = term[2];
	    index = 2;
	  }
	if(term[3] < min)
	  {
	    min   = term[3];
	    index = 3;
	  }
	
	rateCategory[i] = index;	
      }
    
    free(diagptable_start); 
    free(freqs);
    free(scalers);   

    return;
  }
       

#ifdef _LOCAL_DATA

void categorizeIterative(tree *localTree, int startIndex, int endIndex)
{
  double *x2_start, *pz;
  char *tipX1;
  int *ex2;
  int pNumber, qNumber;

  pNumber = localTree->td[0].ti[0].pNumber;
  qNumber = localTree->td[0].ti[0].qNumber;
  pz      = localTree->td[0].ti[0].qz;
   
  if(localTree->td[0].count > 1)    
    newviewIterative(localTree, 0, (endIndex - startIndex));   

  x2_start = getLikelihoodArray(qNumber, localTree->mxtips, localTree->xVector);
  ex2      = getScalingArray(qNumber,    localTree->mySpan, localTree->mxtips, localTree->expArray);  

  tipX1   = &localTree->strided_yVector[pNumber][startIndex];

  if(localTree->mixedData)
    {
      assert(0);
    }
  else
    {
      switch(localTree->likelihoodFunction)
	{   
	case GTRGAMMA:
	  categorizeGTRGAMMA(tipX1, x2_start, ex2, localTree->EIGN_DNA, pz, localTree->gammaRates, 
			     localTree->tipVectorDNA, &(localTree->strided_rateCategory[startIndex]), 
			     0, (endIndex - startIndex)); 
	  break;
	case GTRGAMMAI:
	  categorizeGTRGAMMAINVAR(tipX1, x2_start, ex2, pz, localTree->gammaRates, localTree->EIGN_DNA, 
				  localTree->tipVectorDNA, &(localTree->strided_rateCategory[startIndex]),
				  0, (endIndex - startIndex), &(localTree->invariant[startIndex]), 
				  localTree->invariants, localTree->frequencies_DNA); 
	  break;
	case GTRGAMMAMULT:
	  categorizeGTRGAMMAMULT(tipX1, x2_start, ex2, pz, localTree->gammaRates, localTree->EIGN_DNA, 
				 localTree->tipVectorDNA, &(localTree->strided_rateCategory[startIndex]),
				 0, (endIndex - startIndex), &(localTree->strided_model[startIndex]), 
				 localTree->NumberOfModels, localTree->multiBranch); 
	  break;
	case GTRGAMMAMULTI:
	  categorizeGTRGAMMAMULTINVAR(tipX1, x2_start, ex2, pz, localTree->gammaRates, localTree->EIGN_DNA, 
				      localTree->tipVectorDNA, &(localTree->strided_rateCategory[startIndex]),
				      0, (endIndex - startIndex), &(localTree->strided_model[startIndex]), 
				      localTree->NumberOfModels, 
				      &(localTree->strided_invariant[startIndex]), localTree->invariants, 
				      localTree->frequencies_DNA, localTree->multiBranch); 
	  break;
	case PROTGAMMA:	  
	  categorizeGTRGAMMAPROT(tipX1, x2_start, ex2, pz, localTree->gammaRates, localTree->EIGN_AA, 
				 localTree->tipVectorAA, 
				 &(localTree->strided_rateCategory[startIndex]),
				 0, (endIndex - startIndex)); 
	  break;
	case PROTGAMMAI:
	  categorizeGTRGAMMAPROTINVAR(tipX1, x2_start, ex2, pz, localTree->gammaRates, localTree->EIGN_AA, 
				      localTree->tipVectorAA, &(localTree->strided_rateCategory[startIndex]),
				      0, (endIndex - startIndex), &(localTree->strided_invariant[startIndex]), 
				      localTree->invariants, 
				      localTree->frequencies_AA); 
	  break;
	case PROTGAMMAMULT:
	  categorizeGTRGAMMAPROTMULT(tipX1, x2_start, ex2, pz, localTree->gammaRates, localTree->EIGN_AA, 
				     localTree->tipVectorAA, &(localTree->strided_rateCategory[startIndex]),
				     0, (endIndex - startIndex), &(localTree->strided_model[startIndex]), 
				     localTree->NumberOfModels, localTree->multiBranch); 
	  break;
	case PROTGAMMAMULTI:
	  categorizeGTRGAMMAPROTMULTINVAR(tipX1, x2_start, ex2, pz, localTree->gammaRates, localTree->EIGN_AA, 
					  localTree->tipVectorAA, &(localTree->strided_rateCategory[startIndex]),
					  0, (endIndex - startIndex), &(localTree->strided_model[startIndex]), 
					  localTree->NumberOfModels, 
					  &(localTree->strided_invariant[startIndex]), 
					  localTree->invariants, localTree->frequencies_AA, localTree->multiBranch); 
	  break;
	default:
	  assert(0);
	}
    }
}

#else

void categorizeIterative(tree *tr, int startIndex, int endIndex)
{
  double *x2_start, *pz;
  char *tipX1;
  int *ex2;
  int pNumber, qNumber;

  pNumber = tr->td[0].ti[0].pNumber;
  qNumber = tr->td[0].ti[0].qNumber;
  pz      = tr->td[0].ti[0].qz;
   
  x2_start = getLikelihoodArray(qNumber, tr->mxtips, tr->xVector);
  ex2      = getScalingArray(qNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);  

  tipX1   = tr->yVector[pNumber];

  if(tr->mixedData)
    {
      int model, branchIndex, width, l, u, offset;

      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(tr->multiBranch)
	    branchIndex = model;
	  else
	    branchIndex = 0;

	  l = tr->partitionData[model].lower;
	  u = tr->partitionData[model].upper;
	  
	  offset =  tr->partitionData[model].modelOffset;

	  width = u - l;

	  switch(tr->partitionData[model].dataType)
	    {
	    case DNA_DATA:
	      switch(tr->rateHetModel)
		{
		case GAMMA:
		  categorizeGTRGAMMA(&(tipX1[l]), &(x2_start[offset]), 
				     &(ex2[l]), &(tr->EIGN_DNA[model * 3]), 
				     &(pz[branchIndex]), &(tr->gammaRates[model * 4]), 
				     &(tr->tipVectorDNA[model * 64]), &(tr->cdta->rateCategory[l]), 
				     0, width);
		  break;
		case GAMMA_I:
		  categorizeGTRGAMMAINVAR(&(tipX1[l]), &(x2_start[offset]), 
					  &(ex2[l]), &(pz[branchIndex]), 
					  &(tr->gammaRates[model * 4]), &(tr->EIGN_DNA[model * 3]), 
					  &(tr->tipVectorDNA[model * 64]), 
					  &(tr->cdta->rateCategory[l]),
					  0, width, &(tr->invariant[l]), &(tr->invariants[model]), 
					  &(tr->frequencies_DNA[model * 4])); 
		  break;
		default:
		  assert(0);
		}
	      break;
	    case AA_DATA:
	      switch(tr->rateHetModel)
		{
		case GAMMA:
		  categorizeGTRGAMMAPROT(&(tipX1[l]), 
					 &(x2_start[offset]), 
					 &(ex2[l]), &(pz[branchIndex]), 
					 &(tr->gammaRates[model * 4]), 
					 &(tr->EIGN_AA[model * 19]), 
					 &(tr->tipVectorAA[model * 460]), 
					 &(tr->cdta->rateCategory[l]),
					 0, width); 
		  break;
		case GAMMA_I:
		  categorizeGTRGAMMAPROTINVAR(&(tipX1[l]), &(x2_start[offset]), 
					      &(ex2[l]), &(pz[branchIndex]), 
					      &(tr->gammaRates[model * 4]), 
					      &(tr->EIGN_AA[model * 19]), 
					      &(tr->tipVectorAA[model * 460]), 
					      &(tr->cdta->rateCategory[l]),
					      0, width, &(tr->invariant[l]), 
					      &(tr->invariants[model]), &(tr->frequencies_AA[model * 20])); 
		  break;
		default:
		  assert(0);
		}
	      break;
	    default:
	      assert(0);
	    }
	}
    }
  else
    {
      switch(tr->likelihoodFunction)
	{   
	case GTRGAMMA:
	  categorizeGTRGAMMA(tipX1, x2_start, ex2, tr->EIGN_DNA, pz, tr->gammaRates, tr->tipVectorDNA, tr->cdta->rateCategory, 
			     startIndex, endIndex); 
	  break;
	case GTRGAMMAI:
	  categorizeGTRGAMMAINVAR(tipX1, x2_start, ex2, pz, tr->gammaRates, tr->EIGN_DNA, tr->tipVectorDNA, tr->cdta->rateCategory,
				  startIndex, endIndex, tr->invariant, tr->invariants, tr->frequencies_DNA); 
	  break;
	case GTRGAMMAMULT:
	  categorizeGTRGAMMAMULT(tipX1, x2_start, ex2, pz, tr->gammaRates, tr->EIGN_DNA, tr->tipVectorDNA, tr->cdta->rateCategory,
				 startIndex, endIndex, tr->model, tr->NumberOfModels, tr->multiBranch); 
	  break;
	case GTRGAMMAMULTI:
	  categorizeGTRGAMMAMULTINVAR(tipX1, x2_start, ex2, pz, tr->gammaRates, tr->EIGN_DNA, tr->tipVectorDNA, tr->cdta->rateCategory,
				      startIndex, endIndex, tr->model, tr->NumberOfModels, tr->invariant, tr->invariants, 
				      tr->frequencies_DNA, tr->multiBranch); 
	  break;
	case PROTGAMMA:
	  categorizeGTRGAMMAPROT(tipX1, x2_start, ex2, pz, tr->gammaRates, tr->EIGN_AA, tr->tipVectorAA, tr->cdta->rateCategory,
				 startIndex, endIndex); 
	  break;
	case PROTGAMMAI:
	  categorizeGTRGAMMAPROTINVAR(tipX1, x2_start, ex2, pz, tr->gammaRates, tr->EIGN_AA, tr->tipVectorAA, tr->cdta->rateCategory,
				      startIndex, endIndex, tr->invariant, tr->invariants, tr->frequencies_AA); 
	  break;
	case PROTGAMMAMULT:
	  categorizeGTRGAMMAPROTMULT(tipX1, x2_start, ex2, pz, tr->gammaRates, tr->EIGN_AA, tr->tipVectorAA, tr->cdta->rateCategory,
				     startIndex, endIndex, tr->model, tr->NumberOfModels, tr->multiBranch); 
	  break;
	case PROTGAMMAMULTI:
	  categorizeGTRGAMMAPROTMULTINVAR(tipX1, x2_start, ex2, pz, tr->gammaRates, tr->EIGN_AA, tr->tipVectorAA, tr->cdta->rateCategory,
					  startIndex, endIndex, tr->model, tr->NumberOfModels, 
					  tr->invariant, tr->invariants, tr->frequencies_AA, tr->multiBranch); 
	  break;
	default:
	  assert(0);
	}
    }
}

#endif

void categorizeGeneric (tree *tr, nodeptr p)
{
  nodeptr q = p->back;
  int i;

  assert(isTip(p->number, tr->rdta->numsp)  || isTip(q->number, tr->rdta->numsp));
  
  if(isTip(q->number, tr->rdta->numsp))
    {
      nodeptr tmp = q;
      q = p;
      p = tmp;
    }         
  
  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  for(i = 0; i < tr->numBranches; i++)    
    tr->td[0].ti[0].qz[i] =  q->z[i];

  tr->td[0].count = 1;

  if(!q->x)
    computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp, tr->numBranches);

#ifdef _LOCAL_DATA
  masterBarrier(THREAD_CATEGORIZE, tr); 
#else

  if(tr->td[0].count > 1)
    {    
#ifdef _USE_PTHREADS  
      masterBarrier(THREAD_NEWVIEW, tr);       
#else
      newviewIterative(tr, 0, tr->cdta->endsite);   
#endif
    }
 
  for(i = 0; i < tr->NumberOfModels; i++)
    {      
      tr->cdta->patrat[i * 4]     = tr->gammaRates[i * 4];
      tr->cdta->patrat[i * 4 + 1] = tr->gammaRates[i * 4 + 1];
      tr->cdta->patrat[i * 4 + 2] = tr->gammaRates[i * 4 + 2];
      tr->cdta->patrat[i * 4 + 3] = tr->gammaRates[i * 4 + 3];
    }
  
  tr->NumberOfCategories = 4 * tr->NumberOfModels;

#ifdef _USE_PTHREADS 
  masterBarrier(THREAD_CATEGORIZE, tr); 
#else
  categorizeIterative(tr, 0,  tr->cdta->endsite);  
#endif  
 
  for(i = 0; i < tr->cdta->endsite; i++)
    {
      double temp, wtemp;
      temp = tr->gammaRates[tr->cdta->rateCategory[i]];     
      tr->cdta->wr[i]  = wtemp = temp * tr->cdta->aliaswgt[i];
      tr->cdta->wr2[i] = temp * wtemp;
    }
#endif 
}
