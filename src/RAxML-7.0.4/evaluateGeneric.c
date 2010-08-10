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


#ifdef _USE_PTHREADS
extern double *reductionBuffer;
extern int NumberOfThreads;
#endif

#ifdef _USE_OMP
extern volatile int             NumberOfThreads;
#include <omp.h>
#endif



static double evaluateGTRCATPROTMULT (int *ex1, int *ex2, int *cptr, int *modelptr, int *wptr,
				      double *x1, double *x2, double *EIGN, double *rptr, double *tipVector, double *pz, 
				      char *tipX1, int lower, int n, int numberOfCategories, int numberOfModels, int multiBranch)
{
  double   
    sum, z, lz = 0.0, ki, lza[19], term;      
  double  
    *diagptable, *diagptable_start, *left, *right;  
  int 
    model, modelCounter, i, l;           
  
  if(!multiBranch)
    {
      z = pz[0];
      
      if (z < zmin) z = zmin;
      lz = log(z); 
    }
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * numberOfCategories * 19 * numberOfModels); 
  
  for(modelCounter = 0; modelCounter < numberOfModels; modelCounter++)
    {      
      if(multiBranch)
	{
	  z = pz[modelCounter];	   
	  if (z < zmin) z = zmin;
	  lz = log(z); 
	}

      for(l = 0; l < 19; l++)      
	lza[l] = EIGN[modelCounter * 19 + l] * lz;             
      
      for(i = 0; i <  numberOfCategories; i++)
	{	
	  ki = rptr[i];	 
	  
	  *diagptable++ = exp (ki * lza[0]);	   
	  *diagptable++ = exp (ki * lza[1]);
	  *diagptable++ = exp (ki * lza[2]);	   
	  *diagptable++ = exp (ki * lza[3]);
	  *diagptable++ = exp (ki * lza[4]);	   
	  *diagptable++ = exp (ki * lza[5]);
	  *diagptable++ = exp (ki * lza[6]);	   
	  *diagptable++ = exp (ki * lza[7]);
	  *diagptable++ = exp (ki * lza[8]);	   
	  *diagptable++ = exp (ki * lza[9]);
	  *diagptable++ = exp (ki * lza[10]);	   
	  *diagptable++ = exp (ki * lza[11]);
	  *diagptable++ = exp (ki * lza[12]);	   
	  *diagptable++ = exp (ki * lza[13]);
	  *diagptable++ = exp (ki * lza[14]);	   
	  *diagptable++ = exp (ki * lza[15]);
	  *diagptable++ = exp (ki * lza[16]);	   
	  *diagptable++ = exp (ki * lza[17]);
	  *diagptable++ = exp (ki * lza[18]);	          
	}
    }
  
  sum = 0.0;   
  
  if(tipX1)
    {          	  
      for (i = lower; i < n; i++) 
	{	    		    
	  model = modelptr[i];
	  left  = &(tipVector[model * 460 + 20 * tipX1[i]]);
	  right = &(x2[20 * i]);
	  
	  diagptable = &diagptable_start[model * 19 * numberOfCategories + 19 * cptr[i]];	           
	  
	  term =  left[0] * right[0];	    
	  term += left[1] * right[1] * *diagptable++; 
	  term += left[2] * right[2] * *diagptable++; 
	  term += left[3] * right[3] * *diagptable++; 
	  term += left[4] * right[4] * *diagptable++;	   
	  term += left[5] * right[5] * *diagptable++; 
	  term += left[6] * right[6] * *diagptable++; 
	  term += left[7] * right[7] * *diagptable++; 
	  term += left[8] * right[8] * *diagptable++;
	  term += left[9] * right[9] * *diagptable++; 	    
	  term += left[10] * right[10] * *diagptable++; 
	  term += left[11] * right[11] * *diagptable++; 
	  term += left[12] * right[12] * *diagptable++; 
	  term += left[13] * right[13] * *diagptable++;
	  term += left[14] * right[14] * *diagptable++; 	    
	  term += left[15] * right[15] * *diagptable++; 
	  term += left[16] * right[16] * *diagptable++; 
	  term += left[17] * right[17] * *diagptable++; 
	  term += left[18] * right[18] * *diagptable++;
	  term += left[19] * right[19] * *diagptable++; 
	  
	  term = (log(term)) + (ex2[i] * log(minlikelihood));	          
	  
	  sum += wptr[i] * term;
	}
              
      free(diagptable_start);               
      return  sum;
    }                      


  for (i = lower; i < n; i++) 
    {		       	
      model      = modelptr[i];
      diagptable = &diagptable_start[model * 19 * numberOfCategories + 19 * cptr[i]];	
      
      left  = &x1[20 * i];
      right = &x2[20 * i];
      
      term =  left[0] * right[0];
      term += left[1] * right[1] * *diagptable++; 
      term += left[2] * right[2] * *diagptable++; 
      term += left[3] * right[3] * *diagptable++; 
      term += left[4] * right[4] * *diagptable++;	   
      term += left[5] * right[5] * *diagptable++; 
      term += left[6] * right[6] * *diagptable++; 
      term += left[7] * right[7] * *diagptable++; 
      term += left[8] * right[8] * *diagptable++;
      term += left[9] * right[9] * *diagptable++; 	    
      term += left[10] * right[10] * *diagptable++; 
      term += left[11] * right[11] * *diagptable++; 
      term += left[12] * right[12] * *diagptable++; 
      term += left[13] * right[13] * *diagptable++;
      term += left[14] * right[14] * *diagptable++; 	    
      term += left[15] * right[15] * *diagptable++; 
      term += left[16] * right[16] * *diagptable++; 
      term += left[17] * right[17] * *diagptable++; 
      term += left[18] * right[18] * *diagptable++;
      term += left[19] * right[19] * *diagptable++; 
      
      term = log(term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
      sum += wptr[i] * term;
    }
  
  free(diagptable_start);        
  return  sum;         
}



static double evaluateGTRCATPROT (int *ex1, int *ex2, int *cptr, int *wptr,
				  double *x1, double *x2, double *EIGN, double *rptr, double *tipVector, double pz, 
				  char *tipX1, int lower, int n, int numberOfCategories)
{
  double   sum, z, lz, ki, term, lza[19];
  double  *diagptable, *diagptable_start,  *left, *right;
  int     i, l;  
               
  z = pz;
  
  if (z < zmin) z = zmin;
  lz = log(z);
  
  for(l = 0; l < 19; l++)      
    lza[l] = EIGN[l] * lz;     
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * numberOfCategories * 19);    
  
  for(i = 0; i <  numberOfCategories; i++)
    {	
      ki = rptr[i];	 
      
      *diagptable++ = exp (ki * lza[0]);	   
      *diagptable++ = exp (ki * lza[1]);
      *diagptable++ = exp (ki * lza[2]);	   
      *diagptable++ = exp (ki * lza[3]);
      *diagptable++ = exp (ki * lza[4]);	   
      *diagptable++ = exp (ki * lza[5]);
      *diagptable++ = exp (ki * lza[6]);	   
      *diagptable++ = exp (ki * lza[7]);
      *diagptable++ = exp (ki * lza[8]);	   
      *diagptable++ = exp (ki * lza[9]);
      *diagptable++ = exp (ki * lza[10]);	   
      *diagptable++ = exp (ki * lza[11]);
      *diagptable++ = exp (ki * lza[12]);	   
      *diagptable++ = exp (ki * lza[13]);
      *diagptable++ = exp (ki * lza[14]);	   
      *diagptable++ = exp (ki * lza[15]);
      *diagptable++ = exp (ki * lza[16]);	   
      *diagptable++ = exp (ki * lza[17]);
      *diagptable++ = exp (ki * lza[18]);	          
    }
  
  sum = 0.0;   
  
  if(tipX1)
    {      
      
      for (i = lower; i < n; i++) 
	{	       	
	  left = &(tipVector[20 * tipX1[i]]);
	  right = &(x2[20 * i]);
	  
	  diagptable = &diagptable_start[19 * cptr[i]];	           
	  
	  term =  left[0] * right[0];	    
	  term += left[1] * right[1] * *diagptable++; 
	  term += left[2] * right[2] * *diagptable++; 
	  term += left[3] * right[3] * *diagptable++; 
	  term += left[4] * right[4] * *diagptable++;	   
	  term += left[5] * right[5] * *diagptable++; 
	  term += left[6] * right[6] * *diagptable++; 
	  term += left[7] * right[7] * *diagptable++; 
	  term += left[8] * right[8] * *diagptable++;
	  term += left[9] * right[9] * *diagptable++; 	    
	  term += left[10] * right[10] * *diagptable++; 
	  term += left[11] * right[11] * *diagptable++; 
	  term += left[12] * right[12] * *diagptable++; 
	  term += left[13] * right[13] * *diagptable++;
	  term += left[14] * right[14] * *diagptable++; 	    
	  term += left[15] * right[15] * *diagptable++; 
	  term += left[16] * right[16] * *diagptable++; 
	  term += left[17] * right[17] * *diagptable++; 
	  term += left[18] * right[18] * *diagptable++;
	  term += left[19] * right[19] * *diagptable++; 
	  
	  
	  term = log(term) + (ex2[i] * log(minlikelihood));	          
	  sum += wptr[i] * term;
	}

      free(diagptable_start);              
      return  sum;
    }                   
    
  for (i = lower; i < n; i++) 
    {		       	
      diagptable = &diagptable_start[19 * cptr[i]];	
      
      left  = &x1[20 * i];
      right = &x2[20 * i];
      
      term =  left[0] * right[0];
      term += left[1] * right[1] * *diagptable++; 
      term += left[2] * right[2] * *diagptable++; 
      term += left[3] * right[3] * *diagptable++; 
      term += left[4] * right[4] * *diagptable++;	   
      term += left[5] * right[5] * *diagptable++; 
      term += left[6] * right[6] * *diagptable++; 
      term += left[7] * right[7] * *diagptable++; 
      term += left[8] * right[8] * *diagptable++;
      term += left[9] * right[9] * *diagptable++; 	    
      term += left[10] * right[10] * *diagptable++; 
      term += left[11] * right[11] * *diagptable++; 
      term += left[12] * right[12] * *diagptable++; 
      term += left[13] * right[13] * *diagptable++;
      term += left[14] * right[14] * *diagptable++; 	    
      term += left[15] * right[15] * *diagptable++; 
      term += left[16] * right[16] * *diagptable++; 
      term += left[17] * right[17] * *diagptable++; 
      term += left[18] * right[18] * *diagptable++;
      term += left[19] * right[19] * *diagptable++; 
      
      term = log(term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
      
      sum += wptr[i] * term;

      
    }
    
  free(diagptable_start);         
  return  sum;         
} 


static double evaluateGTRCATMULT (int *ex1, int *ex2, int *cptr, int *modelptr, int *wptr,
				  double *x1_start, double *x2_start, double *EIGN, double *rptr, 
				  double *tipVector, double *pz, 
				  char *tipX1, int lower, int n, int numberOfCategories, int numberOfModels, int multiBranch)
{
  double   sum, z, lz = 0.0, ki, lz1, lz2, lz3, term;      
  double  *diagptable, *diagptable_start, *x1, *x2;
  int model, modelCounter, i;         

  if(!multiBranch)
    {
      z = pz[0];      
      if (z < zmin) z = zmin;
      lz = log(z);
    }  
    
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * numberOfCategories * 3 * numberOfModels);
    
  for(modelCounter = 0; modelCounter < numberOfModels; modelCounter++)
    {
      if(multiBranch)
	{	
	  z = pz[modelCounter];      
	  if (z < zmin) z = zmin;
	  lz = log(z);
	}

      lz1 = EIGN[modelCounter * 3] * lz;
      lz2 = EIGN[modelCounter * 3 + 1] * lz;
      lz3 = EIGN[modelCounter * 3 + 2] * lz;          
      
      for(i = 0; i <  numberOfCategories; i++)
	{	
	  ki = rptr[i];	 
	  *diagptable++ = exp (ki * lz1);
	  *diagptable++ = exp (ki * lz2);
	  *diagptable++ = exp (ki * lz3);	
	}
    }
	    
  sum = 0.0;   
    
  if(tipX1)
    {                    
      for (i = lower; i < n; i++) 
	{	    		
	 

	  model = modelptr[i];
	  x1 = &(tipVector[model * 64 + 4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];
	  
	  diagptable = &diagptable_start[model * 3 *  numberOfCategories + 3 * cptr[i]];	    	    
	  
	  term =  x1[0] * x2[0];	   
	  term += x1[1] * x2[1] * *diagptable++;	  
	  term += x1[2] * x2[2] * *diagptable++;	    
	  term += x1[3] * x2[3] * *diagptable++; 
	  
	  term = log(term) + (ex2[i] * log(minlikelihood));	   	    	   	  	 	  	  

	  sum += wptr[i] * term;	
	}

      
      free(diagptable_start); 
      
      return  sum;
    }                   
  

  for (i = lower; i < n; i++) 
    {		             
      model = modelptr[i];
      x1 = &x1_start[4 * i];
      x2 = &x2_start[4 * i];
      
      diagptable = &diagptable_start[model * 3 *  numberOfCategories + 3 * cptr[i]];	
      
      term =  x1[0] * x2[0];
      term += x1[1] * x2[1] * *diagptable++;
      term += x1[2] * x2[2] * *diagptable++;
      term += x1[3] * x2[3] * *diagptable++;
      
      term = log(term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
      sum += wptr[i] * term;
    }
       
  free(diagptable_start); 
   
  return  sum;         
}
    

static double evaluateGTRCAT (int *ex1, int *ex2, int *cptr, int *wptr,
			      double *x1_start, double *x2_start, double *EIGN, double *rptr, double *tipVector, 
			      double pz, 
			      char *tipX1, int lower, int n, int numberOfCategories)
{
  double  sum, z, lz, ki, lz1, lz2, lz3, term;       
  int     i;  
  double  *diagptable, *diagptable_start, *x1, *x2;                   
   
  z = pz;
  
  if (z < zmin) z = zmin;
  lz = log(z);
  
  lz1 = EIGN[0] * lz;
  lz2 = EIGN[1] * lz;
  lz3 = EIGN[2] * lz;
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * numberOfCategories * 3);    
 
  for(i = 0; i <  numberOfCategories; i++)
    {	
      ki = rptr[i];	 
      *diagptable++ = exp(ki * lz1);
      *diagptable++ = exp(ki * lz2);
      *diagptable++ = exp(ki * lz3);	
    }
	    
  sum = 0.0;   
 

  if(tipX1)
    {      
      for (i = lower; i < n; i++) 
	{	    		   	  
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[4 * i];
	  
	  diagptable = &diagptable_start[3 * cptr[i]];	    	    
	  
	  term =  x1[0] * x2[0];	   
	  term += x1[1] * x2[1] * *diagptable++;	  
	  term += x1[2] * x2[2] * *diagptable++;	    
	  term += x1[3] * x2[3] * *diagptable++; 
	  
	  term = log(term) + (ex2[i] * log(minlikelihood));	   	    	   	 	  	  	 
		
	  sum += wptr[i] * term;
	}
	
     
      free(diagptable_start); 
      
      return  sum;
    }               


  for (i = lower; i < n; i++) 
    {		          	
      x1 = &x1_start[4 * i];
      x2 = &x2_start[4 * i];
      
      diagptable = &diagptable_start[3 * cptr[i]];	
      
      term =  x1[0] * x2[0];
      term += x1[1] * x2[1] * *diagptable++;
      term += x1[2] * x2[2] * *diagptable++;
      term += x1[3] * x2[3] * *diagptable++;
      
      term = log(term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
      
      sum += wptr[i] * term;
    }
	   
  free(diagptable_start); 
       
  return  sum;         
} 


static double evaluateGTRGAMMA(int *ex1, int *ex2, int *wptr,
			       double *x1_start, double *x2_start, double *EIGN, double *gammaRates, 
			       double *tipVector, double pz, 
			       char *tipX1, int lower, const int n)
{
  double   sum = 0.0, z, lz, term, ki;    
  int     i;
  double  *diagptable, *x1, *x2; 
           
  z = pz; 

  if (z < zmin) z = zmin;
  lz = log(z);
  
  diagptable = (double *)malloc(sizeof(double) * 16);

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];	 
      diagptable[i * 4]     = 1.0;
      diagptable[i * 4 + 1] = exp (EIGN[0] * ki * lz);
      diagptable[i * 4 + 2] = exp (EIGN[1] * ki * lz);
      diagptable[i * 4 + 3] = exp (EIGN[2] * ki * lz);
    }      

  if(tipX1)
    {    
#ifdef _USE_OMP
      {
	double *sums = (double*)malloc(sizeof(double) * NumberOfThreads);
#pragma omp parallel private(x1, x2, term)      
	{
	  double local_sum = 0.0;
	  int ref = omp_get_thread_num();
#pragma omp for	
	  for (i = 0; i < n; i++)
#else
	  for (i = lower; i < n; i++) 
#endif
	    {
	      x1 = &(tipVector[4 * tipX1[i]]);	 
	      x2 = &x2_start[16 * i];	          
	  
	      /* cat 0 */
	    
	      term =  x1[0] * x2[0] * diagptable[0];	 
	      term += x1[1] * x2[1] * diagptable[1];	
	      term += x1[2] * x2[2] * diagptable[2];	  
	      term += x1[3] * x2[3] * diagptable[3];     
	      
	      /* cat 1 */
	      
	      term += x1[0] * x2[4] * diagptable[4];
	      term += x1[1] * x2[5] * diagptable[5];
	      term += x1[2] * x2[6] * diagptable[6];
	      term += x1[3] * x2[7] * diagptable[7];     
	      
	      /* cat 2 */
	      
	      term += x1[0] * x2[8] * diagptable[8];
	      term += x1[1] * x2[9] * diagptable[9];
	      term += x1[2] * x2[10] * diagptable[10];
	      term += x1[3] * x2[11] * diagptable[11];     
	      
	      /* cat 3 */
	      
	      term += x1[0] * x2[12] * diagptable[12];
	      term += x1[1] * x2[13] * diagptable[13];
	      term += x1[2] * x2[14] * diagptable[14];
	      term += x1[3] * x2[15] * diagptable[15];     	  	  
	      
	      term = log(0.25 * term) + ex2[i] * log(minlikelihood);	 

#ifndef _USE_OMP
	      sum += wptr[i] * term;
	    }
#else
              local_sum += wptr[i] * term;
	    }
	
	sums[ref] = local_sum;
      }
      for(i = 0; i < NumberOfThreads; i++)
	sum += sums[i];
      free(sums);
    }
#endif
	            
      free(diagptable); 
            
      return  sum;
    }
  else
    {   
#ifdef _USE_OMP      
      {
	double *sums = (double*)malloc(sizeof(double) * NumberOfThreads);
#pragma omp parallel private(x1, x2, term)      
	{
	  double local_sum = 0.0;
	  int ref = omp_get_thread_num();
#pragma omp for	
	  for (i = 0; i < n; i++)  
      
#else
	  for (i = lower; i < n; i++) 
#endif
	    {	  	 	  	  
	      x1 = &x1_start[16 * i];
	      x2 = &x2_start[16 * i];
	      
	      /* cat 0 */
	      
	      term =  x1[0] * x2[0] * diagptable[0];
	      term += x1[1] * x2[1] * diagptable[1];
	      term += x1[2] * x2[2] * diagptable[2];
	      term += x1[3] * x2[3] * diagptable[3];     
	      
	      /* cat 1 */
	      
	      term += x1[4] * x2[4] * diagptable[4];
	      term += x1[5] * x2[5] * diagptable[5];
	      term += x1[6] * x2[6] * diagptable[6];
	      term += x1[7] * x2[7] * diagptable[7];     
	      
	      /* cat 2 */
	      
	      term += x1[8] * x2[8] * diagptable[8];
	      term += x1[9] * x2[9] * diagptable[9];
	      term += x1[10] * x2[10] * diagptable[10];
	      term += x1[11] * x2[11] * diagptable[11];     
	      
	      /* cat 3 */
	      
	      term += x1[12] * x2[12] * diagptable[12];
	      term += x1[13] * x2[13] * diagptable[13];
	      term += x1[14] * x2[14] * diagptable[14];
	      term += x1[15] * x2[15] * diagptable[15];     
	      
	      term = log(0.25 * term) + (ex1[i] + ex2[i]) * log(minlikelihood);
#ifndef _USE_OMP
	      sum += wptr[i] * term;
	    }
#else
              local_sum += wptr[i] * term;
	    }
	
	sums[ref] = local_sum;
      }
      for(i = 0; i < NumberOfThreads; i++)
	sum += sums[i];
      free(sums);
    }
#endif
            
      free(diagptable); 
    	
      return  sum;
    }
} 






static double evaluateGTRGAMMAINVAR (int *ex1, int *ex2, int *wptr, int *iptr,
				     double *x1_start, double *x2_start, double *EIGN, double *gammaRates, 
				     double *tipVector, double *tFreqs, double invariants, double pz, 
				     char *tipX1, int lower, int n)
{ 
  int     i;
  double  *diagptable, *x1, *x2; 
  double 
    freqs[4], 
    scaler = 0.25 * (1.0 - invariants),
    sum = 0.0, 
    z, lz, term, ki; 

  freqs[0] = tFreqs[0] * invariants; 
  freqs[1] = tFreqs[1] * invariants;
  freqs[2] = tFreqs[2] * invariants;
  freqs[3] = tFreqs[3] * invariants;
  
  z = pz;
  
  if (z < zmin) z = zmin;
 
  lz = log(z);
  
  diagptable = (double *)malloc(sizeof(double) * 16); 

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];	 
      
      diagptable[i * 4]     = 1.0;
      diagptable[i * 4 + 1] = exp (EIGN[0] * ki * lz);
      diagptable[i * 4 + 2] = exp (EIGN[1] * ki * lz);
      diagptable[i * 4 + 3] = exp (EIGN[2] * ki * lz);
    }      

  if(tipX1)
    {   
      for (i = lower; i < n; i++) 
	{
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &x2_start[16 * i];
	  
	  /* cat 0 */
	    
	  term =  x1[0] * x2[0] * diagptable[0];	 
	  term += x1[1] * x2[1] * diagptable[1];	
	  term += x1[2] * x2[2] * diagptable[2];	  
	  term += x1[3] * x2[3] * diagptable[3];     
	
	  /* cat 1 */
	  
	  term += x1[0] * x2[4] * diagptable[4];
	  term += x1[1] * x2[5] * diagptable[5];
	  term += x1[2] * x2[6] * diagptable[6];
	  term += x1[3] * x2[7] * diagptable[7];     
	  
	  /* cat 2 */
	  
	  term += x1[0] * x2[8] * diagptable[8];
	  term += x1[1] * x2[9] * diagptable[9];
	  term += x1[2] * x2[10] * diagptable[10];
	  term += x1[3] * x2[11] * diagptable[11];     
	  
	  /* cat 3 */
	  
	  term += x1[0] * x2[12] * diagptable[12];
	  term += x1[1] * x2[13] * diagptable[13];
	  term += x1[2] * x2[14] * diagptable[14];
	  term += x1[3] * x2[15] * diagptable[15];	 

	  if(iptr[i] < 4)	   
	    term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	  else
	    term = log(scaler * term) + (ex2[i] * log(minlikelihood));	 

	    sum += wptr[i] * term;
	}
	    	  	  	            
      free(diagptable); 
      
      return  sum;
    }
  else
    {           		

      for (i = lower; i < n; i++) 
	{	  	 	  	
	  x1 = &x1_start[16 * i];
	  x2 = &x2_start[16 * i];
	  
	  /* cat 0 */
	    
	  term =  x1[0] * x2[0] * diagptable[0];	 
	  term += x1[1] * x2[1] * diagptable[1];	
	  term += x1[2] * x2[2] * diagptable[2];	  
	  term += x1[3] * x2[3] * diagptable[3];     
	
	  /* cat 1 */
	  
	  term += x1[4] * x2[4] * diagptable[4];
	  term += x1[5] * x2[5] * diagptable[5];
	  term += x1[6] * x2[6] * diagptable[6];
	  term += x1[7] * x2[7] * diagptable[7];     
	  
	  /* cat 2 */
	  
	  term += x1[8] * x2[8] * diagptable[8];
	  term += x1[9] * x2[9] * diagptable[9];
	  term += x1[10] * x2[10] * diagptable[10];
	  term += x1[11] * x2[11] * diagptable[11];     
	  
	  /* cat 3 */
	  
	  term += x1[12] * x2[12] * diagptable[12];
	  term += x1[13] * x2[13] * diagptable[13];
	  term += x1[14] * x2[14] * diagptable[14];
	  term += x1[15] * x2[15] * diagptable[15];	      
	  
	  if(iptr[i] < 4)	   
	    term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i] + ex1[i]));
	  else
	    term = log(scaler * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));

	    sum += wptr[i] * term;
	}
	  	            
      free(diagptable);            
      return  sum;
    }
} 


static double evaluateGTRGAMMAMULT (int *ex1, int *ex2, int *modelptr, int *wptr,
				    double *x1_start, double *x2_start, double *EIGN, double *gammaRates, 
				    double *tipVector, double *pz, 
				    char *tipX1, int lower, int n, int numberOfModels, int multiBranch)
{
  double   sum = 0.0, z, lz = 0.0, term, ki;      
  int     i, model; 
  double  *diagptable, *diagptable_start, *x1, *x2; 
    
  if(!multiBranch)
    {
      z = pz[0]; 
      if (z < zmin) z = zmin;
      lz = log(z);
    }
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 12 * numberOfModels);

  for(model = 0; model < numberOfModels; model++)
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
	  
	  *diagptable++ = exp (EIGN[model * 3] * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 1] * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 2] * ki * lz);
	}
    }
	  
  if(tipX1)
    {         
      for (i = lower; i < n; i++) 
	{
	  model = modelptr[i];

	  x1 = &(tipVector[64 * model + 4 * tipX1[i]]);	 
	  x2 = &x2_start[16 * i];
	   
	  diagptable = &diagptable_start[12 * model];
	  
	  /* cat 0 */
	    
	  term =  x1[0] * x2[0];	 
	  term += x1[1] * x2[1] * *diagptable++;	
	  term += x1[2] * x2[2] * *diagptable++;	  
	  term += x1[3] * x2[3] * *diagptable++;     
	
	  /* cat 1 */
	  
	  term += x1[0] * x2[4];
	  term += x1[1] * x2[5] * *diagptable++;
	  term += x1[2] * x2[6] * *diagptable++;
	  term += x1[3] * x2[7] * *diagptable++;     
	  
	  /* cat 2 */
	  
	  term += x1[0] * x2[8];
	  term += x1[1] * x2[9] * *diagptable++;
	  term += x1[2] * x2[10] * *diagptable++;
	  term += x1[3] * x2[11] * *diagptable++;     
	  
	  /* cat 3 */
	  
	  term += x1[0] * x2[12];
	  term += x1[1] * x2[13] * *diagptable++;
	  term += x1[2] * x2[14] * *diagptable++;
	  term += x1[3] * x2[15] * *diagptable++;     
	  
	  term = log(0.25 * term) + (ex2[i] * log(minlikelihood));

	  sum += wptr[i] * term;
	}
	            
      free(diagptable_start); 
            
      return  sum;
    }
  else
    {                 
      for (i = lower; i < n; i++) 
	{	
	  model = modelptr[i];

	  diagptable = &diagptable_start[model * 12];

	  x1 = &x1_start[16 * i];
	  x2 = &x2_start[16 * i];
	    
	  /* cat 0 */
	  
	  term =  x1[0] * x2[0];
	  term += x1[1] * x2[1] * *diagptable++;
	  term += x1[2] * x2[2] * *diagptable++;
	  term += x1[3] * x2[3] * *diagptable++;     
	  
	  /* cat 1 */
	  
	  term += x1[4] * x2[4];
	  term += x1[5] * x2[5] * *diagptable++;
	  term += x1[6] * x2[6] * *diagptable++;
	  term += x1[7] * x2[7] * *diagptable++;     
	  
	  /* cat 2 */
	  
	  term += x1[8] * x2[8];
	  term += x1[9] * x2[9] * *diagptable++;
	  term += x1[10] * x2[10] * *diagptable++;
	  term += x1[11] * x2[11] * *diagptable++;     
	  
	  /* cat 3 */
	  
	  term +=  x1[12] * x2[12];
	  term += x1[13] * x2[13] * *diagptable++;
	  term += x1[14] * x2[14] * *diagptable++;
	  term += x1[15] * x2[15] * *diagptable++;     
	    
	  term = log(0.25 * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));

	  sum += wptr[i] * term;
	}

      free(diagptable_start); 
          
      return  sum;
    }
}

static double evaluateGTRGAMMAMULTINVAR(int *ex1, int *ex2, int *modelptr, int *wptr, int *iptr,
					double *x1_start, double *x2_start, double *EIGN, double *gammaRates, 
					double *tipVector, double *tFreqs, double *invariants, double *pz, 
					char *tipX1, int lower, int n, int numberOfModels, int multiBranch)
{
  double   sum = 0.0, z, lz = 0.0, term, ki;      
  int     i, model;  
  double  *diagptable, *diagptable_start, *x1, *x2, *scalers, *freqs; 
    
  if(!multiBranch)
    {
      z = pz[0]; 
      if (z < zmin) z = zmin;
      lz = log(z);
    }
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 12 * numberOfModels);
  freqs      = (double *)malloc(4 * numberOfModels * sizeof(double));
  scalers    = (double *)malloc(numberOfModels * sizeof(double)); 

  for(model = 0; model < numberOfModels; model++)
    {
      if(multiBranch)
	{
	  z = pz[model]; 
	  if (z < zmin) z = zmin;
	  lz = log(z);
	}

      scalers[model] = 0.25 * (1.0 - invariants[model]); 
      freqs[4 * model]     = tFreqs[4 * model]     * invariants[model]; 
      freqs[4 * model + 1] = tFreqs[4 * model + 1] * invariants[model];
      freqs[4 * model + 2] = tFreqs[4 * model + 2] * invariants[model];
      freqs[4 * model + 3] = tFreqs[4 * model + 3] * invariants[model];
   
      diagptable = &diagptable_start[12 * model];
      
      for(i = 0; i < 4; i++)
	{
	  ki = gammaRates[model * 4 + i];	 
	  
	  *diagptable++ = exp (EIGN[model * 3] * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 1] * ki * lz);
	  *diagptable++ = exp (EIGN[model * 3 + 2] * ki * lz);
	}
    }

	  
  if(tipX1)
    {           
      for (i = lower; i < n; i++) 
	{
	  model = modelptr[i];

	  x1 = &(tipVector[64 * model + 4 * tipX1[i]]);	 
	  x2 = &x2_start[16 * i];
	   
	  diagptable = &diagptable_start[12 * model];
	  
	  /* cat 0 */
	    
	  term =  x1[0] * x2[0];	 
	  term += x1[1] * x2[1] * *diagptable++;	
	  term += x1[2] * x2[2] * *diagptable++;	  
	  term += x1[3] * x2[3] * *diagptable++;     
	
	  /* cat 1 */
	  
	  term += x1[0] * x2[4];
	  term += x1[1] * x2[5] * *diagptable++;
	  term += x1[2] * x2[6] * *diagptable++;
	  term += x1[3] * x2[7] * *diagptable++;     
	  
	  /* cat 2 */
	  
	  term += x1[0] * x2[8];
	  term += x1[1] * x2[9] * *diagptable++;
	  term += x1[2] * x2[10] * *diagptable++;
	  term += x1[3] * x2[11] * *diagptable++;     
	  
	  /* cat 3 */
	  
	  term += x1[0] * x2[12];
	  term += x1[1] * x2[13] * *diagptable++;
	  term += x1[2] * x2[14] * *diagptable++;
	  term += x1[3] * x2[15] * *diagptable++;     	  

	  if(iptr[i] < 4)	  
	    term = log(((scalers[model] * term) + freqs[model * 4 + iptr[i]]) * pow(minlikelihood, ex2[i]));
	  else
	    term = log(scalers[model] * term) + (ex2[i] * log(minlikelihood));

	    sum += wptr[i] * term;
	}
  	  	                 
      free(diagptable_start); 
      free(scalers);
      free(freqs);
           
      return  sum;
    }
  else
    {            
      for (i = lower; i < n; i++) 
	{	
	  model = modelptr[i];

	  diagptable = &diagptable_start[model * 12];

	  x1 = &x1_start[16 * i];
	  x2 = &x2_start[16 * i];
	    
	  /* cat 0 */
	  
	  term =  x1[0] * x2[0];
	  term += x1[1] * x2[1] * *diagptable++;
	  term += x1[2] * x2[2] * *diagptable++;
	  term += x1[3] * x2[3] * *diagptable++;     
	  
	  /* cat 1 */
	  
	  term += x1[4] * x2[4];
	  term += x1[5] * x2[5] * *diagptable++;
	  term += x1[6] * x2[6] * *diagptable++;
	  term += x1[7] * x2[7] * *diagptable++;     
	  
	  /* cat 2 */
	  
	  term += x1[8] * x2[8];
	  term += x1[9] * x2[9] * *diagptable++;
	  term += x1[10] * x2[10] * *diagptable++;
	  term += x1[11] * x2[11] * *diagptable++;     
	  
	  /* cat 3 */
	  
	  term +=  x1[12] * x2[12];
	  term += x1[13] * x2[13] * *diagptable++;
	  term += x1[14] * x2[14] * *diagptable++;
	  term += x1[15] * x2[15] * *diagptable++;     	    	  

	  if(iptr[i] < 4)	    
	    term = log(((scalers[model] * term) + freqs[model * 4 + iptr[i]]) * pow(minlikelihood, ex2[i] + ex1[i]));
	  else
	    term = log(scalers[model] * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
	    sum += wptr[i] * term;
	}
      
      free(diagptable_start); 
      free(scalers);
      free(freqs);
        
      return  sum;
    }
}

static double evaluateGTRGAMMAPROT (int *ex1, int *ex2, int *wptr,
				    double *x1, double *x2, double *EIGN, double *gammaRates, 
				    double *tipVector, double pz, 
				    char *tipX1, int lower, int n)
{
  double   sum, z, lz, term, ki;        
  int     i, j;   
  double  *diagptable, *diagptable_start, *left, *right;           
  
  z = pz;
  
  if (z < zmin) z = zmin;
  lz = log(z);
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 4 * 19);
  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];	 
      
      for(j = 0; j < 19; j++)
	*diagptable++ = exp (EIGN[j] * ki * lz);	
    }
  
  sum = 0.0;
  
  if(tipX1)
    {            
      for (i = lower; i < n; i++) 
	{
	  left = &(tipVector[20 * tipX1[i]]);
	  diagptable = diagptable_start;
	  
	  term = 0;
	  
	  for(j = 0; j < 4; j++)
	    {
	      right = &(x2[80 * i + 20 * j]);
	      
	      term +=  left[0] * right[0];
	      term +=  left[1] * right[1] * *diagptable++;
	      term +=  left[2] * right[2] * *diagptable++;
	      term +=  left[3] * right[3] * *diagptable++;
	      term +=  left[4] * right[4] * *diagptable++;
	      term +=  left[5] * right[5] * *diagptable++;
	      term +=  left[6] * right[6] * *diagptable++;	
	      term +=  left[7] * right[7] * *diagptable++;
	      term +=  left[8] * right[8] * *diagptable++;
	      term +=  left[9] * right[9] * *diagptable++;
	      term +=  left[10] * right[10] * *diagptable++;
	      term +=  left[11] * right[11] * *diagptable++;
	      term +=  left[12] * right[12] * *diagptable++;
	      term +=  left[13] * right[13] * *diagptable++;
	      term +=  left[14] * right[14] * *diagptable++;
	      term +=  left[15] * right[15] * *diagptable++;
	      term +=  left[16] * right[16] * *diagptable++;
	      term +=  left[17] * right[17] * *diagptable++;
	      term +=  left[18] * right[18] * *diagptable++;	
	      term +=  left[19] * right[19] * *diagptable++;
	    }	  
	  
	  term = log(0.25 * term) + (ex2[i] * log(minlikelihood));	   
	  
	  sum += wptr[i] * term;
	}
                
      free(diagptable_start); 
      
      return  sum;
    }              

  for (i = lower; i < n; i++) 
    {	  	 	  
      diagptable = diagptable_start;
      
      term = 0;
      
      for(j = 0; j < 4; j++)
	{
	  left  = &(x1[80 * i + 20 * j]);
	  right = &(x2[80 * i + 20 * j]);	    
	  term +=  left[0] * right[0];
	  term +=  left[1] * right[1] * *diagptable++;
	  term +=  left[2] * right[2] * *diagptable++;
	  term +=  left[3] * right[3] * *diagptable++;
	  term +=  left[4] * right[4] * *diagptable++;
	  term +=  left[5] * right[5] * *diagptable++;
	  term +=  left[6] * right[6] * *diagptable++;	
	  term +=  left[7] * right[7] * *diagptable++;
	  term +=  left[8] * right[8] * *diagptable++;
	  term +=  left[9] * right[9] * *diagptable++;
	  term +=  left[10] * right[10] * *diagptable++;
	  term +=  left[11] * right[11] * *diagptable++;
	  term +=  left[12] * right[12] * *diagptable++;
	  term +=  left[13] * right[13] * *diagptable++;
	  term +=  left[14] * right[14] * *diagptable++;
	  term +=  left[15] * right[15] * *diagptable++;
	  term +=  left[16] * right[16] * *diagptable++;
	  term +=  left[17] * right[17] * *diagptable++;
	  term +=  left[18] * right[18] * *diagptable++;	
	  term +=  left[19] * right[19] * *diagptable++;	
	}
      
      term = log(0.25 * term) + ((ex1[i] + ex2[i])*log(minlikelihood));
      
      sum += wptr[i] * term;
    }
            
  free(diagptable_start); 
       
  return  sum;
}

static double evaluateGTRGAMMAPROTINVAR (int *ex1, int *ex2, int *wptr, int *iptr,
					 double *x1, double *x2, double *EIGN, double *gammaRates, 
					 double *tipVector,double *tFreqs, double invariants, double pz, 
					 char *tipX1, int lower, int n)
{
  double   
    sum, z, lz, term, ki, freqs[20],
    scaler = 0.25 * (1.0 - invariants);        
  int     i, j;     
  double  
    *diagptable, *diagptable_start, *left, *right;   
    
  for(i = 0; i < 20; i++)
    freqs[i] = tFreqs[i] * invariants;
            
  z = pz;
  
  if (z < zmin) z = zmin;
  lz = log(z);
  
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 4 * 19);

  for(i = 0; i < 4; i++)
    {
      ki = gammaRates[i];	 
      
      for(j = 0; j < 19; j++)
	*diagptable++ = exp (EIGN[j] * ki * lz);	
    }
	    
  sum = 0.0;
  
  if(tipX1)
    {          
      for (i = lower; i < n; i++) 
	{
	  left = &(tipVector[20 * tipX1[i]]);
	  diagptable = diagptable_start;
	  
	  term = 0;
	  
	  for(j = 0; j < 4; j++)
	    {
	      right = &(x2[80 * i + 20 * j]);
	      
	      term +=  left[0] * right[0];
	      term +=  left[1] * right[1] * *diagptable++;
	      term +=  left[2] * right[2] * *diagptable++;
	      term +=  left[3] * right[3] * *diagptable++;
	      term +=  left[4] * right[4] * *diagptable++;
	      term +=  left[5] * right[5] * *diagptable++;
	      term +=  left[6] * right[6] * *diagptable++;	
	      term +=  left[7] * right[7] * *diagptable++;
	      term +=  left[8] * right[8] * *diagptable++;
	      term +=  left[9] * right[9] * *diagptable++;
	      term +=  left[10] * right[10] * *diagptable++;
	      term +=  left[11] * right[11] * *diagptable++;
	      term +=  left[12] * right[12] * *diagptable++;
	      term +=  left[13] * right[13] * *diagptable++;
	      term +=  left[14] * right[14] * *diagptable++;
	      term +=  left[15] * right[15] * *diagptable++;
	      term +=  left[16] * right[16] * *diagptable++;
	      term +=  left[17] * right[17] * *diagptable++;
	      term +=  left[18] * right[18] * *diagptable++;	
	      term +=  left[19] * right[19] * *diagptable++;
	    }	  
	  
	  if(iptr[i] < 20)	   
	    term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
	  else
	    term = log(scaler * term) + (ex2[i] * log(minlikelihood));
	  
	  sum += wptr[i] * term;
	}
          
      free(diagptable_start); 
      
      return  sum;
    }                
    
  for (i = lower; i < n; i++) 
    {	  	 	  
      diagptable = diagptable_start;
      
      term = 0;
      
      for(j = 0; j < 4; j++)
	{
	  left  = &(x1[80 * i + 20 * j]);
	  right = &(x2[80 * i + 20 * j]);	    
	  term +=  left[0] * right[0];
	  term +=  left[1] * right[1] * *diagptable++;
	  term +=  left[2] * right[2] * *diagptable++;
	  term +=  left[3] * right[3] * *diagptable++;
	  term +=  left[4] * right[4] * *diagptable++;
	  term +=  left[5] * right[5] * *diagptable++;
	  term +=  left[6] * right[6] * *diagptable++;	
	  term +=  left[7] * right[7] * *diagptable++;
	  term +=  left[8] * right[8] * *diagptable++;
	  term +=  left[9] * right[9] * *diagptable++;
	  term +=  left[10] * right[10] * *diagptable++;
	  term +=  left[11] * right[11] * *diagptable++;
	  term +=  left[12] * right[12] * *diagptable++;
	  term +=  left[13] * right[13] * *diagptable++;
	  term +=  left[14] * right[14] * *diagptable++;
	  term +=  left[15] * right[15] * *diagptable++;
	  term +=  left[16] * right[16] * *diagptable++;
	  term +=  left[17] * right[17] * *diagptable++;
	  term +=  left[18] * right[18] * *diagptable++;	
	  term +=  left[19] * right[19] * *diagptable++;	
	}
      
      if(iptr[i] < 20)	   
	term = log(((scaler * term) + freqs[iptr[i]]) * pow(minlikelihood, ex2[i] + ex1[i]));
      else
	term = log(scaler * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
      sum += wptr[i] * term;
    }
            
  free(diagptable_start); 
       
  return  sum;
}


static double evaluateGTRGAMMAPROTMULT (int *ex1, int *ex2, int *modelptr, int *wptr,
					double *x1, double *x2, double *EIGN, double *gammaRates, 
					double *tipVector, double *pz, 
					char *tipX1, int lower, int n, int numberOfModels, int multiBranch)
{
  double   sum, z, lz = 0.0, term, ki;        
  int     i, j, model;
  double  *diagptable, *diagptable_start, *left, *right;  
                
  if(!multiBranch)
    {
      z = pz[0];
      if (z < zmin) z = zmin;
      lz = log(z);
    }
       
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 4 * 19 * numberOfModels);

  for(model = 0; model < numberOfModels; model++)
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
	    
  sum = 0.0;
      
  if(tipX1)
    {       	
      for (i = lower; i < n; i++) 
	{	    
	  model = modelptr[i];
	  
	  left = &(tipVector[460 * model + 20 * tipX1[i]]);
	  diagptable = &diagptable_start[76 * model];	   
	  
	  term = 0;
	  
	  for(j = 0; j < 4; j++)
	    {
	      right = &(x2[80 * i + 20 * j]);
	      
	      term +=  left[0] * right[0];
	      term +=  left[1] * right[1] * *diagptable++;
	      term +=  left[2] * right[2] * *diagptable++;
	      term +=  left[3] * right[3] * *diagptable++;
	      term +=  left[4] * right[4] * *diagptable++;
	      term +=  left[5] * right[5] * *diagptable++;
	      term +=  left[6] * right[6] * *diagptable++;	
	      term +=  left[7] * right[7] * *diagptable++;
	      term +=  left[8] * right[8] * *diagptable++;
	      term +=  left[9] * right[9] * *diagptable++;
	      term +=  left[10] * right[10] * *diagptable++;
	      term +=  left[11] * right[11] * *diagptable++;
	      term +=  left[12] * right[12] * *diagptable++;
	      term +=  left[13] * right[13] * *diagptable++;
	      term +=  left[14] * right[14] * *diagptable++;
	      term +=  left[15] * right[15] * *diagptable++;
	      term +=  left[16] * right[16] * *diagptable++;
	      term +=  left[17] * right[17] * *diagptable++;
	      term +=  left[18] * right[18] * *diagptable++;	
	      term +=  left[19] * right[19] * *diagptable++;
	    }	  
	  
	  term = log(0.25 * term) + (ex2[i] * log(minlikelihood));	   
	  
	  sum += wptr[i] * term;
	}
      
      free(diagptable_start); 
      
      return  sum;
    }               
  
  for (i = lower; i < n; i++) 
    {	  	      	  
      diagptable = &diagptable_start[76 * modelptr[i]];
      
      term = 0;
      
      for(j = 0; j < 4; j++)
	{
	  left  = &(x1[80 * i + 20 * j]);
	  right = &(x2[80 * i + 20 * j]);	    
	  term +=  left[0] * right[0];
	  term +=  left[1] * right[1] * *diagptable++;
	  term +=  left[2] * right[2] * *diagptable++;
	  term +=  left[3] * right[3] * *diagptable++;
	  term +=  left[4] * right[4] * *diagptable++;
	  term +=  left[5] * right[5] * *diagptable++;
	  term +=  left[6] * right[6] * *diagptable++;	
	  term +=  left[7] * right[7] * *diagptable++;
	  term +=  left[8] * right[8] * *diagptable++;
	  term +=  left[9] * right[9] * *diagptable++;
	  term +=  left[10] * right[10] * *diagptable++;
	  term +=  left[11] * right[11] * *diagptable++;
	  term +=  left[12] * right[12] * *diagptable++;
	  term +=  left[13] * right[13] * *diagptable++;
	  term +=  left[14] * right[14] * *diagptable++;
	  term +=  left[15] * right[15] * *diagptable++;
	  term +=  left[16] * right[16] * *diagptable++;
	  term +=  left[17] * right[17] * *diagptable++;
	  term +=  left[18] * right[18] * *diagptable++;	
	  term +=  left[19] * right[19] * *diagptable++;	
	}
      
      term = log(0.25 * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
      
      sum += wptr[i] * term;
    }
  
  free(diagptable_start); 
      
  return  sum;
}






static double evaluateGTRGAMMAPROTMULTINVAR(int *ex1, int *ex2, int *modelptr, int *wptr, int *iptr,
					    double *x1, double *x2, double *EIGN, double *gammaRates, 
					    double *tipVector, double *tFreqs, double *invariants, double *pz, 
					    char *tipX1, int lower, int n, int numberOfModels, int multiBranch)
{
  double   sum, z, lz = 0.0, term, ki;        
  int     i, j, model;  
  double  *diagptable, *diagptable_start, *left, *right, *scalers, *freqs;    
    
  if(!multiBranch)
    {
      z = pz[0];   
      if (z < zmin) z = zmin;
      lz = log(z);
    }
       
  diagptable = diagptable_start = (double *)malloc(sizeof(double) * 4 * 19 * numberOfModels);
  freqs      = (double *)malloc(20 * numberOfModels * sizeof(double));
  scalers    = (double *)malloc(numberOfModels * sizeof(double));

  for(model = 0; model < numberOfModels; model++)
    {
      if(multiBranch)
	{
	  z = pz[model];   
	  if (z < zmin) z = zmin;
	  lz = log(z);
	}
     
      diagptable = &diagptable_start[4 * 19 * model];
      scalers[model] = 0.25 * (1.0 - invariants[model]); 
	
      for(i = 0; i < 20; i++)
	freqs[20 * model + i] = tFreqs[20 * model + i] * invariants[model]; 

      for(i = 0; i < 4; i++)
	{
	  ki = gammaRates[model * 4 + i];	 
	  
	  for(j = 0; j < 19; j++)
	    *diagptable++ = exp (EIGN[model * 19 + j] * ki * lz);	
	}
    }	    
	    
  sum = 0.0;
      
  if(tipX1)
    {      	  
      for (i = lower; i < n; i++) 
	{	    
	  model = modelptr[i];
	  
	  left = &(tipVector[460 * model + 20 * tipX1[i]]);
	  diagptable = &diagptable_start[76 * model];	   
	  
	  term = 0;
	  
	  for(j = 0; j < 4; j++)
	    {
	      right = &(x2[80 * i + 20 * j]);
	      
	      term +=  left[0] * right[0];
	      term +=  left[1] * right[1] * *diagptable++;
	      term +=  left[2] * right[2] * *diagptable++;
	      term +=  left[3] * right[3] * *diagptable++;
	      term +=  left[4] * right[4] * *diagptable++;
	      term +=  left[5] * right[5] * *diagptable++;
	      term +=  left[6] * right[6] * *diagptable++;	
	      term +=  left[7] * right[7] * *diagptable++;
	      term +=  left[8] * right[8] * *diagptable++;
	      term +=  left[9] * right[9] * *diagptable++;
	      term +=  left[10] * right[10] * *diagptable++;
	      term +=  left[11] * right[11] * *diagptable++;
	      term +=  left[12] * right[12] * *diagptable++;
	      term +=  left[13] * right[13] * *diagptable++;
	      term +=  left[14] * right[14] * *diagptable++;
	      term +=  left[15] * right[15] * *diagptable++;
	      term +=  left[16] * right[16] * *diagptable++;
	      term +=  left[17] * right[17] * *diagptable++;
	      term +=  left[18] * right[18] * *diagptable++;	
	      term +=  left[19] * right[19] * *diagptable++;
	    }	  	        	  	   
	  
	  if(iptr[i] < 20)	  
	    term = log(((scalers[model] * term) + freqs[model * 20 + iptr[i]]) * pow(minlikelihood, ex2[i]));
	  else
	    term = log(scalers[model] * term) + (ex2[i] * log(minlikelihood));
	  
	  sum += wptr[i] * term;
	}
	   
      free(diagptable_start); 
      free(freqs);
      free(scalers);
      
      return  sum;
    }      
 
  for (i = lower; i < n; i++) 
    {	  	   
      model = modelptr[i];
      diagptable = &diagptable_start[76 * model];
      
      term = 0;
      
      for(j = 0; j < 4; j++)
	{
	  left  = &(x1[80 * i + 20 * j]);
	  right = &(x2[80 * i + 20 * j]);	    
	  term +=  left[0] * right[0];
	  term +=  left[1] * right[1] * *diagptable++;
	  term +=  left[2] * right[2] * *diagptable++;
	  term +=  left[3] * right[3] * *diagptable++;
	  term +=  left[4] * right[4] * *diagptable++;
	  term +=  left[5] * right[5] * *diagptable++;
	  term +=  left[6] * right[6] * *diagptable++;	
	  term +=  left[7] * right[7] * *diagptable++;
	  term +=  left[8] * right[8] * *diagptable++;
	  term +=  left[9] * right[9] * *diagptable++;
	  term +=  left[10] * right[10] * *diagptable++;
	  term +=  left[11] * right[11] * *diagptable++;
	  term +=  left[12] * right[12] * *diagptable++;
	  term +=  left[13] * right[13] * *diagptable++;
	  term +=  left[14] * right[14] * *diagptable++;
	  term +=  left[15] * right[15] * *diagptable++;
	  term +=  left[16] * right[16] * *diagptable++;
	  term +=  left[17] * right[17] * *diagptable++;
	  term +=  left[18] * right[18] * *diagptable++;	
	  term +=  left[19] * right[19] * *diagptable++;	
	}
      
      if(iptr[i] < 20)	    
	term = log(((scalers[model] * term) + freqs[model * 20 + iptr[i]]) * pow(minlikelihood, ex2[i] + ex1[i]));
      else
	term = log(scalers[model] * term) + ((ex1[i] + ex2[i]) * log(minlikelihood));
      
      sum += wptr[i] * term;
    }
	                
  free(diagptable_start); 
  free(freqs);
  free(scalers);
  
  return  sum;
}






/*********************************************************************************************/

static double evaluateMixedData(int model, tree *tr, char *tip, int *ex1, int *ex2, 
				double *x1_start, double *x2_start, double pz)
{
  double result = 0.0;

  /*int l = tr->modelIndices[model][0];
    int u = tr->modelIndices[model][1];	  	*/

  int l = tr->partitionData[model].lower;
  int u = tr->partitionData[model].upper;

  int width  = u - l;
  /*int offset = tr->modelOffsets[model];   */
  int offset = tr->partitionData[model].modelOffset;

  switch(/*tr->dataType[model]*/ tr->partitionData[model].dataType)
    { 
    case DNA_DATA:
      switch(tr->rateHetModel)
	{
	case CAT:
	  if(tip)
	    result = evaluateGTRCAT(ex1, &(ex2[l]), &(tr->cdta->rateCategory[l]), &(tr->cdta->aliaswgt[l]),
				    x1_start, &(x2_start[offset]), &(tr->EIGN_DNA[model * 3]), 
				    &(tr->cdta->patrat[0]), &(tr->tipVectorDNA[model * 64]), pz, 
				    &(tip[l]), 0, width, tr->NumberOfCategories); 
	  else
	    result = evaluateGTRCAT(&(ex1[l]), &(ex2[l]), &(tr->cdta->rateCategory[l]), &(tr->cdta->aliaswgt[l]),
				    &(x1_start[offset]), &(x2_start[offset]), &(tr->EIGN_DNA[model * 3]), 
				    &(tr->cdta->patrat[0]), &(tr->tipVectorDNA[model * 64]), pz, 
				    tip, 0, width, tr->NumberOfCategories);
	  
	  break;
	case GAMMA:
	  if(tip)
	    result = evaluateGTRGAMMA(ex1, &(ex2[l]), &(tr->cdta->aliaswgt[l]),
				      x1_start, &(x2_start[offset]), 
				      &(tr->EIGN_DNA[model * 3]), &(tr->gammaRates[model * 4]), 
				      &(tr->tipVectorDNA[model * 64]), pz, 
				      &(tip[l]), 0, width); 
	  else
	    result = evaluateGTRGAMMA(&(ex1[l]), &(ex2[l]), &(tr->cdta->aliaswgt[l]),
				      &(x1_start[offset]), &(x2_start[offset]), 
				      &(tr->EIGN_DNA[model * 3]), &(tr->gammaRates[model * 4]), 
				      &(tr->tipVectorDNA[model * 64]), pz, 
				      tip, 0, width); 
	  
	  break;
	case GAMMA_I:
	  if(tip)
	    result = evaluateGTRGAMMAINVAR(ex1, &(ex2[l]), &(tr->cdta->aliaswgt[l]), 
					   &(tr->invariant[l]),
					   x1_start, &(x2_start[offset]), 
					   &(tr->EIGN_DNA[model * 3]), 
					   &(tr->gammaRates[model * 4]), 
					   &(tr->tipVectorDNA[model * 64]), 
					   &(tr->frequencies_DNA[model * 4]), 
					   tr->invariants[model], pz, 
					   &(tip[l]), 0, width);
	  else
	    result = evaluateGTRGAMMAINVAR(&(ex1[l]), &(ex2[l]), &(tr->cdta->aliaswgt[l]), 
					   &(tr->invariant[l]),
					   &(x1_start[offset]), &(x2_start[offset]), 
					   &(tr->EIGN_DNA[model * 3]), 
					   &(tr->gammaRates[model * 4]), 
					   &(tr->tipVectorDNA[model * 64]), 
					   &(tr->frequencies_DNA[model * 4]), 
					   tr->invariants[model], pz, 
					   tip, 0, width);	  	  
	  break;
	default:
	  assert(0);
	}
      break;
    case AA_DATA:
      switch(tr->rateHetModel)
	{
	case CAT:
	  if(tip)
	    result = evaluateGTRCATPROT(ex1, &(ex2[l]), &(tr->cdta->rateCategory[l]), &(tr->cdta->aliaswgt[l]),
					x1_start, &(x2_start[offset]), 
					&(tr->EIGN_AA[19 * model]), &(tr->cdta->patrat[0]), &(tr->tipVectorAA[model * 460]), 
					pz, &(tip[l]), 0, width, tr->NumberOfCategories);
	  else
	    result = evaluateGTRCATPROT(&(ex1[l]), &(ex2[l]), &(tr->cdta->rateCategory[l]), &(tr->cdta->aliaswgt[l]),
					&(x1_start[offset]), &(x2_start[offset]), 
					&(tr->EIGN_AA[19 * model]), &(tr->cdta->patrat[0]), &(tr->tipVectorAA[model * 460]), 
					pz, tip, 0, width, tr->NumberOfCategories);
	  
	  break;
	case GAMMA:
	  if(tip)
	    result = evaluateGTRGAMMAPROT(ex1, &(ex2[l]), &(tr->cdta->aliaswgt[l]),
					  x1_start, &(x2_start[offset]), 
					  &(tr->EIGN_AA[19 * model]), &(tr->gammaRates[model * 4]), 
					  &(tr->tipVectorAA[model * 460]), pz, 
					  &(tip[l]), 0, width);
	  else
	    result = evaluateGTRGAMMAPROT(&(ex1[l]), &(ex2[l]), &(tr->cdta->aliaswgt[l]),
					  &(x1_start[offset]), &(x2_start[offset]), 
					  &(tr->EIGN_AA[19 * model]), &(tr->gammaRates[model * 4]), 
					  &(tr->tipVectorAA[model * 460]), pz, 
					  tip, 0, width);
	  break;
	case GAMMA_I:
	  if(tip)
	    result = evaluateGTRGAMMAPROTINVAR(ex1, &(ex2[l]), &(tr->cdta->aliaswgt[l]), &(tr->invariant[l]),
					       x1_start, &(x2_start[offset]), 
					       &(tr->EIGN_AA[model * 19]), 
					       &(tr->gammaRates[model * 4]), 
					       &(tr->tipVectorAA[model * 460]), 
					       &(tr->frequencies_AA[model * 20]), 
					       tr->invariants[model], 
					       pz, 
					       &(tip[l]), 0, width);
	  else
	    result = evaluateGTRGAMMAPROTINVAR(&(ex1[l]), &(ex2[l]), &(tr->cdta->aliaswgt[l]), &(tr->invariant[l]),
					       &(x1_start[offset]), &(x2_start[offset]), 
					       &(tr->EIGN_AA[model * 19]), 
					       &(tr->gammaRates[model * 4]), 
					       &(tr->tipVectorAA[model * 460]), 
					       &(tr->frequencies_AA[model * 20]), 
					       tr->invariants[model], 
					       pz, 
					       tip, 0, width);
	  break;
	default:
	  assert(0);
	}
      break;
    default:
      assert(0);
    }	  
  return result;
}


#ifdef _LOCAL_DATA

double evaluateIterative(tree *localTree, int startIndex, int endIndex)
{
  double 
    result = 0.0,
    *x1_start = (double*)NULL, 
    *x2_start = (double*)NULL;
  int    
    *ex1 = (int*)NULL, 
    *ex2 = (int*)NULL;
  char *tip = (char*)NULL;   
  int pNumber, qNumber;
  double *pz;

  pNumber = localTree->td[0].ti[0].pNumber;
  qNumber = localTree->td[0].ti[0].qNumber;
  pz      = localTree->td[0].ti[0].qz;

 
  newviewIterative(localTree, startIndex, endIndex);

  if(isTip(pNumber, localTree->mxtips) || isTip(qNumber, localTree->mxtips))
    {	        	    
      if(isTip(qNumber, localTree->mxtips))
	{	
	  
	  x2_start = getLikelihoodArray(pNumber, localTree->mxtips, localTree->xVector);
	  ex2      = getScalingArray(pNumber, localTree->mySpan, localTree->mxtips, localTree->expArray);

	  tip = localTree->strided_yVector[qNumber];	 	      
	}           
      else
	{
	  x2_start = getLikelihoodArray(qNumber, localTree->mxtips, localTree->xVector);
	  ex2      = getScalingArray(qNumber, localTree->mySpan, localTree->mxtips, localTree->expArray); 

	  tip = localTree->strided_yVector[pNumber];
	}
    }
  else
    {                 
      x1_start = getLikelihoodArray(pNumber, localTree->mxtips, localTree->xVector);
      x2_start = getLikelihoodArray(qNumber, localTree->mxtips, localTree->xVector);
      ex1      = getScalingArray(pNumber, localTree->mySpan, localTree->mxtips, localTree->expArray);
      ex2      = getScalingArray(qNumber, localTree->mySpan, localTree->mxtips, localTree->expArray);
    }  

  

  if(localTree->mixedData)
    {                           
      assert(0);

      /*
	for(model = 0; model < localTree->NumberOfModels; model++)
	{
	if(multiBranch)
	branchIndex = model;
	else
	branchIndex = 0;
	result += evaluateMixedData(model, tr, tip, ex1, ex2, x1_start, x2_start, pz[branchIndex]);	  	          
	}
      */
    }
  else
    {
      switch(localTree->likelihoodFunction)
	{
	case GTRCAT:
	  result =  evaluateGTRCAT(ex1, ex2, &(localTree->strided_rateCategory[startIndex]), 
				   &(localTree->strided_aliaswgt[startIndex]),
				   x1_start, x2_start, 
				   localTree->EIGN_DNA, localTree->strided_patrat, localTree->tipVectorDNA, pz[0], 
				   tip, 0, (endIndex - startIndex), localTree->NumberOfCategories);	 
	  break;
	case GTRCATMULT:    
	  result = evaluateGTRCATMULT(ex1, ex2, &(localTree->strided_rateCategory[startIndex]), 
				      &(localTree->strided_model[startIndex]), 
				      &(localTree->strided_aliaswgt[startIndex]),
				      x1_start, x2_start, localTree->EIGN_DNA, 
				      localTree->strided_patrat, localTree->tipVectorDNA, pz, 
				      tip, 0, (endIndex - startIndex), localTree->NumberOfCategories, localTree->NumberOfModels, 
				      localTree->multiBranch);       
	  break;
	case PROTCAT:
	  result = evaluateGTRCATPROT(ex1, ex2, &(localTree->strided_rateCategory[startIndex]), &(localTree->strided_aliaswgt[startIndex]),
				      x1_start, x2_start, localTree->EIGN_AA, localTree->strided_patrat, localTree->tipVectorAA, pz[0], 
				      tip, 0, (endIndex - startIndex), localTree->NumberOfCategories); 
	  break;
	case PROTCATMULT:
	  result = evaluateGTRCATPROTMULT(ex1, ex2, &(localTree->strided_rateCategory[startIndex]), 
					  &(localTree->strided_model[startIndex]), &(localTree->strided_aliaswgt[startIndex]),
					  x1_start, x2_start, localTree->EIGN_AA, localTree->strided_patrat, localTree->tipVectorAA, pz, 
					  tip, 0, (endIndex - startIndex), localTree->NumberOfCategories, localTree->NumberOfModels,
					  localTree->multiBranch); 
	  break;
	case GTRGAMMA:
	  result = evaluateGTRGAMMA(ex1, ex2, &(localTree->strided_aliaswgt[startIndex]),
				    x1_start, x2_start, localTree->EIGN_DNA, localTree->gammaRates, localTree->tipVectorDNA, pz[0], 
				    tip, 0, (endIndex - startIndex)); 
	  break;
	case GTRGAMMAI:
	  result = evaluateGTRGAMMAINVAR(ex1, ex2, &(localTree->strided_aliaswgt[startIndex]), &(localTree->strided_invariant[startIndex]),
					 x1_start, x2_start, localTree->EIGN_DNA, localTree->gammaRates, 
					 localTree->tipVectorDNA, localTree->frequencies_DNA, localTree->invariants[0], pz[0], 
					 tip, 0, (endIndex - startIndex)); 
	  break;
	case GTRGAMMAMULT:
	  result = evaluateGTRGAMMAMULT(ex1, ex2, &(localTree->strided_model[startIndex]), &(localTree->strided_aliaswgt[startIndex]),
					x1_start, x2_start, localTree->EIGN_DNA, localTree->gammaRates, localTree->tipVectorDNA, pz, 
					tip, 0, (endIndex - startIndex), localTree->NumberOfModels, localTree->multiBranch); 
	  break;
	case GTRGAMMAMULTI:
	  result = evaluateGTRGAMMAMULTINVAR(ex1, ex2, &(localTree->strided_model[startIndex]), &(localTree->strided_aliaswgt[startIndex]), 
					     &(localTree->strided_invariant[startIndex]),
					     x1_start, x2_start, localTree->EIGN_DNA, localTree->gammaRates, 
					     localTree->tipVectorDNA, localTree->frequencies_DNA, localTree->invariants, pz, 
					     tip, 0, (endIndex - startIndex), localTree->NumberOfModels, localTree->multiBranch); 
	  break;
	case PROTGAMMA:
	  result = evaluateGTRGAMMAPROT(ex1, ex2, &(localTree->strided_aliaswgt[startIndex]),
					x1_start, x2_start, localTree->EIGN_AA, localTree->gammaRates, localTree->tipVectorAA, pz[0], 
					tip, 0, (endIndex - startIndex)); 
	  break;
	case PROTGAMMAI:
	  result = evaluateGTRGAMMAPROTINVAR(ex1, ex2, &(localTree->strided_aliaswgt[startIndex]), 
					     &(localTree->strided_invariant[startIndex]),
					     x1_start, x2_start, localTree->EIGN_AA, localTree->gammaRates, 
					     localTree->tipVectorAA, localTree->frequencies_AA, localTree->invariants[0], pz[0], 
					     tip, 0, (endIndex - startIndex)); 
	  break;
	case PROTGAMMAMULT:
	  result = evaluateGTRGAMMAPROTMULT(ex1, ex2, &(localTree->strided_model[startIndex]), &(localTree->strided_aliaswgt[startIndex]),
					    x1_start, x2_start, localTree->EIGN_AA, localTree->gammaRates, localTree->tipVectorAA, pz, 
					    tip, 0, (endIndex - startIndex), localTree->NumberOfModels, localTree->multiBranch); 
	  break;
	case PROTGAMMAMULTI:     
	  result = evaluateGTRGAMMAPROTMULTINVAR(ex1, ex2, &(localTree->strided_model[startIndex]), 
						 &(localTree->strided_aliaswgt[startIndex]), &(localTree->strided_invariant[startIndex]),
						 x1_start, x2_start, localTree->EIGN_AA, localTree->gammaRates, 
						 localTree->tipVectorAA, localTree->frequencies_AA, localTree->invariants, pz, 
						 tip, 0, (endIndex - startIndex), localTree->NumberOfModels, localTree->multiBranch
						 ); 
	  break;
	default:
	  assert(0);
	}
    }

  return result;
}

#else

double evaluateIterative(tree *tr, int startIndex, int endIndex)
{
  double 
    result = 0.0,
    *x1_start = (double*)NULL, 
    *x2_start = (double*)NULL;
  int    
    *ex1 = (int*)NULL, 
    *ex2 = (int*)NULL;
  char *tip = (char*)NULL;   
  int pNumber, qNumber;
  double *pz;

  pNumber = tr->td[0].ti[0].pNumber;
  qNumber = tr->td[0].ti[0].qNumber;
  pz      = tr->td[0].ti[0].qz;

  newviewIterative(tr, startIndex, endIndex);

  if(isTip(pNumber, tr->rdta->numsp) || isTip(qNumber, tr->rdta->numsp))
    {	        	    
      if(isTip(qNumber, tr->rdta->numsp))
	{	
	  
	  x2_start = getLikelihoodArray(pNumber, tr->mxtips, tr->xVector);
	  ex2      = getScalingArray(pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);	 

	  tip = tr->yVector[qNumber];	 	      
	}           
      else
	{
	  x2_start = getLikelihoodArray(qNumber, tr->mxtips, tr->xVector);
	  ex2      = getScalingArray(qNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);	 

	  tip = tr->yVector[pNumber];
	}
    }
  else
    {                 
      x1_start = getLikelihoodArray(pNumber, tr->mxtips, tr->xVector);
      x2_start = getLikelihoodArray(qNumber, tr->mxtips, tr->xVector);
      ex1      = getScalingArray(pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
      ex2      = getScalingArray(qNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);     
    }


  if(tr->mixedData)
    {
      int model, branchIndex;           
            
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  if(tr->multiBranch)
	    branchIndex = model;
	  else
	    branchIndex = 0;
	  result += evaluateMixedData(model, tr, tip, ex1, ex2, x1_start, x2_start, pz[branchIndex]);	  	          
	}
    }
  else
    {
      switch(tr->likelihoodFunction)
	{
	case GTRCAT:
	  result =  evaluateGTRCAT(ex1, ex2, tr->cdta->rateCategory, tr->cdta->aliaswgt,
				   x1_start, x2_start, tr->EIGN_DNA, tr->cdta->patrat, tr->tipVectorDNA, pz[0], 
				   tip, startIndex, endIndex, tr->NumberOfCategories);   
	  break;
	case GTRCATMULT:    
	  result = evaluateGTRCATMULT(ex1, ex2, tr->cdta->rateCategory, tr->model, tr->cdta->aliaswgt,
				      x1_start, x2_start, tr->EIGN_DNA, tr->cdta->patrat, tr->tipVectorDNA, pz, 
				      tip, startIndex, endIndex, tr->NumberOfCategories, tr->NumberOfModels, tr->multiBranch);       
	  break;
	case PROTCAT:
	  result = evaluateGTRCATPROT(ex1, ex2, tr->cdta->rateCategory, tr->cdta->aliaswgt,
				      x1_start, x2_start, tr->EIGN_AA, tr->cdta->patrat, tr->tipVectorAA, pz[0], 
				      tip, startIndex, endIndex, tr->NumberOfCategories); 
	  break;
	case PROTCATMULT:
	  result = evaluateGTRCATPROTMULT(ex1, ex2, tr->cdta->rateCategory, tr->model, tr->cdta->aliaswgt,
					  x1_start, x2_start, tr->EIGN_AA, tr->cdta->patrat, tr->tipVectorAA, pz, 
					  tip, startIndex, endIndex, tr->NumberOfCategories, tr->NumberOfModels, tr->multiBranch); 
	  break;
	case GTRGAMMA:
	  result = evaluateGTRGAMMA(ex1, ex2, tr->cdta->aliaswgt,
				    x1_start, x2_start, tr->EIGN_DNA, tr->gammaRates, tr->tipVectorDNA, pz[0], 
				    tip, startIndex, endIndex); 
	  break;
	case GTRGAMMAI:
	  result = evaluateGTRGAMMAINVAR(ex1, ex2, tr->cdta->aliaswgt, tr->invariant,
					 x1_start, x2_start, tr->EIGN_DNA, tr->gammaRates, 
					 tr->tipVectorDNA, tr->frequencies_DNA, tr->invariants[0], pz[0], 
					 tip, startIndex, endIndex); 
	  break;
	case GTRGAMMAMULT:
	  result = evaluateGTRGAMMAMULT(ex1, ex2, tr->model, tr->cdta->aliaswgt,
					x1_start, x2_start, tr->EIGN_DNA, tr->gammaRates, tr->tipVectorDNA, pz, 
					tip, startIndex, endIndex, tr->NumberOfModels, tr->multiBranch); 
	  break;
	case GTRGAMMAMULTI:
	  result = evaluateGTRGAMMAMULTINVAR(ex1, ex2, tr->model, tr->cdta->aliaswgt, tr->invariant,
					     x1_start, x2_start, tr->EIGN_DNA, tr->gammaRates, 
					     tr->tipVectorDNA, tr->frequencies_DNA, tr->invariants, pz, 
					     tip, startIndex, endIndex, tr->NumberOfModels, tr->multiBranch); 
	  break;
	case PROTGAMMA:
	  result = evaluateGTRGAMMAPROT(ex1, ex2, tr->cdta->aliaswgt,
					x1_start, x2_start, tr->EIGN_AA, tr->gammaRates, tr->tipVectorAA, pz[0], 
					tip, startIndex, endIndex); 
	  break;
	case PROTGAMMAI:
	  result = evaluateGTRGAMMAPROTINVAR(ex1, ex2, tr->cdta->aliaswgt, tr->invariant,
					     x1_start, x2_start, tr->EIGN_AA, tr->gammaRates, 
					     tr->tipVectorAA, tr->frequencies_AA, tr->invariants[0], pz[0], 
					     tip, startIndex, endIndex); 
	  break;
	case PROTGAMMAMULT:
	  result = evaluateGTRGAMMAPROTMULT(ex1, ex2, tr->model, tr->cdta->aliaswgt,
					    x1_start, x2_start, tr->EIGN_AA, tr->gammaRates, tr->tipVectorAA, pz, 
					    tip, startIndex, endIndex, tr->NumberOfModels, tr->multiBranch); 
	  break;
	case PROTGAMMAMULTI:     
	  result = evaluateGTRGAMMAPROTMULTINVAR(ex1, ex2, tr->model, tr->cdta->aliaswgt, tr->invariant,
						 x1_start, x2_start, tr->EIGN_AA, tr->gammaRates, 
						 tr->tipVectorAA, tr->frequencies_AA, tr->invariants, pz, 
						 tip, startIndex, endIndex, tr->NumberOfModels, tr->multiBranch
						 ); 
	  break;
	default:
	  assert(0);
	}
    }

  return result;
}

#endif




double evaluateGeneric (tree *tr, nodeptr p)
{
  double result;
  nodeptr q = p->back; 
  int i;

#ifdef _MULTI_GENE
  if(tr->doMulti)
    {
      int model;
      
      result = 0.0;
      
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  double erg;	  	 	  	 

	  tr->td[model].ti[0].pNumber = p->number;
	  tr->td[model].ti[0].qNumber = q->number;
	  
	  tr->td[model].ti[0].qz[model] =  q->z[model];
	  
	  tr->td[model].count = 1;
	      
	  if(!p->xs[model])
	    computeMultiTraversalInfo(p, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->rdta->numsp, model);
	  if(!q->xs[model])
	    computeMultiTraversalInfo(q, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->rdta->numsp, model);
	  
	  printf("Model %d Count %d %d %d\n", model, tr->td[model].count, (p->backs[model])?1:0, (q->backs[model])?1:0);

	  /*
	    if(!p->backs[model])
	    {	     
	    assert(q->backs[model] != p);
	    result += tr->perPartitionLH[model];
	    
	    newviewIterativePartition(tr, tr->partitionData[model].lower, tr->partitionData[model].upper, model);
	    }
	    else
	  */
	  if(p->backs[model])
	    {	    	   	
	      tr->td[model].ti[0].qNumber = p->backs[model]->number;	  	      
   		
	      /* printf("%d %d\n", (p->backs[model])?1:0, (q->backs[model])?1:0); */
	      result += (erg = evaluateIterativePartition(tr, tr->partitionData[model].lower, 
							  tr->partitionData[model].upper, model));	      	     
	    }
	  else
	    {
	      result += tr->perPartitionLH[model];
	    
	      newviewIterativePartition(tr, tr->partitionData[model].lower, tr->partitionData[model].upper, model);
	    }
	}     
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
      
#ifdef _USE_PTHREADS 
      masterBarrier(THREAD_EVALUATE, tr); 
      for(i = 0, result = 0.0; i < NumberOfThreads; i++)	  
	result += reductionBuffer[i];  	           
#else
      result = evaluateIterative(tr, 0, tr->cdta->endsite);
#endif
    }
  
  tr->likelihood = result;    
 
  return result;
}

double evaluateGenericInitrav (tree *tr, nodeptr p)
{
  double result;   

#ifdef _MULTI_GENE
  if(tr->doMulti)
    {	
      int model;
      result = 0.0;
      determineFullMultiTraversal(tr);
      for(model = 0; model < tr->NumberOfModels; model++)          
	{
	  tr->perPartitionLH[model] = evaluateIterativePartition(tr, tr->partitionData[model].lower, 
								 tr->partitionData[model].upper, model);
	  result += tr->perPartitionLH[model];
	}
    }
  else
#endif
    {
      determineFullTraversal(p, tr);
      
#ifdef _USE_PTHREADS 
      {
	int i;
	
	masterBarrier(THREAD_EVALUATE, tr); 
	for(i = 0, result = 0.0; i < NumberOfThreads; i++)    
	  result += reductionBuffer[i];  	  
      }
#else
      result = evaluateIterative(tr, 0, tr->cdta->endsite);
#endif
    }
     
  tr->likelihood = result;    
      
  return result;
}

void onlyInitravPartition(tree *tr, nodeptr p, int model)
{   
  int lower, upper;

  determineFullTraversal(p, tr); 

  lower = tr->partitionData[model].lower;
  upper = tr->partitionData[model].upper;

  tr->modelNumber = model;
#ifdef _USE_PTHREADS  
  masterBarrier(THREAD_NEWVIEW_PARTITION, tr);  	 
#else
  newviewIterativePartition(tr, lower, upper, model);   
#endif     
}

void onlyInitrav(tree *tr, nodeptr p)
{   
  determineFullTraversal(p, tr);  

#ifdef _USE_PTHREADS  
  masterBarrier(THREAD_NEWVIEW, tr);  	 
#else
  newviewIterative(tr, 0, tr->cdta->endsite);   
#endif     
}



double evaluateGenericInitravPartition(tree *tr, nodeptr p, int model)
{
  double result;   
  int lower, upper;

  determineFullTraversal(p, tr); 

  lower = tr->partitionData[model].lower;
  upper = tr->partitionData[model].upper;

#ifdef _USE_PTHREADS 
  {
    int i;
    tr->modelNumber = model;    
    masterBarrier(THREAD_EVALUATE_PARTITION, tr); 
    for(i = 0, result = 0.0; i < NumberOfThreads; i++)      
      result += reductionBuffer[i];               
  }
#else
  result = evaluateIterativePartition(tr, lower, upper, model);
#endif
  tr->likelihood = result;    

  return result;
}


#ifdef _LOCAL_DATA

double evaluateIterativePartition(tree *localTree, int lower, int upper, int model)
{
  double 
    result = 0.0,
    *x1_start = (double*)NULL, 
    *x2_start = (double*)NULL;
  int    
    *ex1 = (int*)NULL, 
    *ex2 = (int*)NULL;
  char *tip = (char*)NULL;   
  int pNumber, qNumber, branchIndex;
  double *pz;

  if(localTree->multiBranch)
    branchIndex = model;
  else
    branchIndex = 0;

  pNumber = localTree->td[0].ti[0].pNumber;
  qNumber = localTree->td[0].ti[0].qNumber;
  pz      = localTree->td[0].ti[0].qz;  

  newviewIterativePartition(localTree, lower, upper, model);

   if(isTip(pNumber, localTree->mxtips) || isTip(qNumber, localTree->mxtips))
    {	        	    
      if(isTip(qNumber, localTree->mxtips))
	{	
	  
	  x2_start = getLikelihoodArray(pNumber, localTree->mxtips, localTree->xVector);
	  ex2      = getScalingArray(pNumber, localTree->mySpan, localTree->mxtips, localTree->expArray);

	  tip = localTree->strided_yVector[qNumber];	 	      
	}           
      else
	{
	  x2_start = getLikelihoodArray(qNumber, localTree->mxtips, localTree->xVector);
	  ex2      = getScalingArray(qNumber, localTree->mySpan, localTree->mxtips, localTree->expArray); 

	  tip = localTree->strided_yVector[pNumber];
	}
    }
  else
    {                 
      x1_start = getLikelihoodArray(pNumber, localTree->mxtips, localTree->xVector);
      x2_start = getLikelihoodArray(qNumber, localTree->mxtips, localTree->xVector);
      ex1      = getScalingArray(pNumber, localTree->mySpan, localTree->mxtips, localTree->expArray);
      ex2      = getScalingArray(qNumber, localTree->mySpan, localTree->mxtips, localTree->expArray);
    }
  

  if(localTree->mixedData)
    {         
      assert(0);
      /*result = evaluateMixedData(model, tr, tip, ex1, ex2, x1_start, x2_start, pz[branchIndex]);*/
    }
  else
    {
      switch(localTree->likelihoodFunction)
	{
	case GTRCAT:
	case GTRCATMULT:    
	  result = evaluateGTRCAT(ex1, ex2, localTree->strided_rateCategory, 
				  localTree->strided_aliaswgt,
				  x1_start, x2_start, 
				  &(localTree->EIGN_DNA[model * 3]), 
				  localTree->strided_patrat, 
				  &(localTree->tipVectorDNA[model * 64]), 
				  pz[branchIndex], 
				  tip, lower, upper, localTree->NumberOfCategories);       
	  break; 
	case GTRGAMMA: 
	case GTRGAMMAMULT:
	  result = evaluateGTRGAMMA(ex1, ex2, localTree->strided_aliaswgt,
				    x1_start, x2_start, 
				    &(localTree->EIGN_DNA[model * 3]), 
				    &(localTree->gammaRates[model * 4]), 
				    &(localTree->tipVectorDNA[model * 64]), 
				    pz[branchIndex], 
				    tip, lower, upper); 
	  break;
	case GTRGAMMAI: /* needed for rate opt*/
	case GTRGAMMAMULTI:
	  result = evaluateGTRGAMMAINVAR(ex1, ex2, localTree->strided_aliaswgt, localTree->strided_invariant,
					 x1_start, x2_start, 
					 &(localTree->EIGN_DNA[model * 3]), 
					 &(localTree->gammaRates[model * 4]), 
					 &(localTree->tipVectorDNA[model * 64]), 
					 &(localTree->frequencies_DNA[model * 4]), 
					 localTree->invariants[model], 
					 pz[branchIndex], 
					 tip, lower, upper); 
	  break;
	case PROTGAMMA:  /* needed for rate opt*/
	case PROTGAMMAMULT:
	  result = evaluateGTRGAMMAPROT(ex1, ex2, localTree->strided_aliaswgt,
					x1_start, x2_start, 
					&(localTree->EIGN_AA[model * 19]), 
					&(localTree->gammaRates[model * 4]), 
					&(localTree->tipVectorAA[model * 460]), 
					pz[branchIndex], 
					tip, lower, upper); 
	  break;
	case PROTGAMMAI: /* needed for rate opt*/
	case PROTGAMMAMULTI:
	  result = evaluateGTRGAMMAPROTINVAR(ex1, ex2, localTree->strided_aliaswgt, 
					     localTree->strided_invariant,
					     x1_start, x2_start, 
					     &(localTree->EIGN_AA[model * 19]), 
					     &(localTree->gammaRates[model * 4]), 
					     &(localTree->tipVectorAA[model * 460]), 
					     &(localTree->frequencies_AA[model * 20]), 
					     localTree->invariants[model], 
					     pz[branchIndex], 
					     tip, lower, upper);         
	  break;
	default:
	  assert(0);
	}
    }

  return result;
}


#else


double evaluateIterativePartition(tree *tr, int lower, int upper, int model)
{
  double 
    result = 0.0,
    *x1_start = (double*)NULL, 
    *x2_start = (double*)NULL;
  int    
    *ex1 = (int*)NULL, 
    *ex2 = (int*)NULL;
  char *tip = (char*)NULL;   
  int pNumber, qNumber, branchIndex;
  double *pz;

  if(tr->multiBranch)
    branchIndex = model;
  else
    branchIndex = 0;


#ifdef _MULTI_GENE
  if(tr->doMulti)
    {
      pNumber = tr->td[model].ti[0].pNumber;
      qNumber = tr->td[model].ti[0].qNumber;
      pz      = tr->td[model].ti[0].qz;
    }
  else
    {
      pNumber = tr->td[0].ti[0].pNumber;
      qNumber = tr->td[0].ti[0].qNumber;
      pz      = tr->td[0].ti[0].qz;
    }
#else
  pNumber = tr->td[0].ti[0].pNumber;
  qNumber = tr->td[0].ti[0].qNumber;
  pz      = tr->td[0].ti[0].qz;
#endif

  newviewIterativePartition(tr, lower, upper, model);

  if(isTip(pNumber, tr->rdta->numsp) || isTip(qNumber, tr->rdta->numsp))
    {	        	    
      if(isTip(qNumber, tr->rdta->numsp))
	{	
	  x2_start = getLikelihoodArray(pNumber, tr->mxtips, tr->xVector);
	  ex2     = getScalingArray(pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);	 

	  tip = tr->yVector[qNumber];	 	      
	}           
      else
	{
	  x2_start = getLikelihoodArray(qNumber, tr->mxtips, tr->xVector);
	  ex2      = getScalingArray(qNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);

	  tip = tr->yVector[pNumber];
	}
    }
  else
    {           
      x1_start = getLikelihoodArray(pNumber, tr->mxtips, tr->xVector);
      x2_start = getLikelihoodArray(qNumber, tr->mxtips, tr->xVector);

      ex1      = getScalingArray(pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
      ex2      = getScalingArray(qNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
      
     
    }

  if(tr->mixedData)
    {         
      result = evaluateMixedData(model, tr, tip, ex1, ex2, x1_start, x2_start, pz[branchIndex]);     
    }
  else
    {
      switch(tr->likelihoodFunction)
	{
	case GTRCAT:
	case GTRCATMULT:    
	  result = evaluateGTRCAT(ex1, ex2, tr->cdta->rateCategory, 
				  tr->cdta->aliaswgt,
				  x1_start, x2_start, 
				  &(tr->EIGN_DNA[model * 3]), 
				  tr->cdta->patrat, 
				  &(tr->tipVectorDNA[model * 64]), 
				  pz[branchIndex], 
				  tip, lower, upper, tr->NumberOfCategories);       
	  break; 
	case GTRGAMMA:  /* needed for rate opt*/
	case GTRGAMMAMULT:
	  result = evaluateGTRGAMMA(ex1, ex2, tr->cdta->aliaswgt,
				    x1_start, x2_start, 
				    &(tr->EIGN_DNA[model * 3]), 
				    &(tr->gammaRates[model * 4]), 
				    &(tr->tipVectorDNA[model * 64]), 
				    pz[branchIndex], 
				    tip, lower, upper); 
	  break;
	case GTRGAMMAI: /* needed for rate opt*/
	case GTRGAMMAMULTI:
	  result = evaluateGTRGAMMAINVAR(ex1, ex2, tr->cdta->aliaswgt, tr->invariant,
					 x1_start, x2_start, 
					 &(tr->EIGN_DNA[model * 3]), 
					 &(tr->gammaRates[model * 4]), 
					 &(tr->tipVectorDNA[model * 64]), 
					 &(tr->frequencies_DNA[model * 4]), 
					 tr->invariants[model], 
					 pz[branchIndex], 
					 tip, lower, upper); 
	  break;
	case PROTGAMMA:  /* needed for rate opt*/
	case PROTGAMMAMULT:
	  result = evaluateGTRGAMMAPROT(ex1, ex2, tr->cdta->aliaswgt,
					x1_start, x2_start, 
					&(tr->EIGN_AA[model * 19]), 
					&(tr->gammaRates[model * 4]), 
					&(tr->tipVectorAA[model * 460]), 
					pz[branchIndex], 
					tip, lower, upper); 
	  break;
	case PROTGAMMAI: /* needed for rate opt*/
	case PROTGAMMAMULTI:
	  result = evaluateGTRGAMMAPROTINVAR(ex1, ex2, tr->cdta->aliaswgt, tr->invariant,
					     x1_start, x2_start, 
					     &(tr->EIGN_AA[model * 19]), 
					     &(tr->gammaRates[model * 4]), 
					     &(tr->tipVectorAA[model * 460]), 
					     &(tr->frequencies_AA[model * 20]), 
					     tr->invariants[model], 
					     pz[branchIndex], 
					     tip, lower, upper);         
	  break;
	default:
	  assert(0);
	}
    }

  return result;
}

#endif

double evaluatePartitionGeneric (tree *tr, nodeptr p, int model)
{
  double result;
  nodeptr q = p->back; 
  int i, lower, upper; 

  lower = tr->partitionData[model].lower;
  upper = tr->partitionData[model].upper;


  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  for(i = 0; i < tr->numBranches; i++)    
    tr->td[0].ti[0].qz[i] =  q->z[i];

  tr->td[0].count = 1;
  if(!p->x)
    computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp, tr->numBranches);
  if(!q->x)
    computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp, tr->numBranches);  
   
#ifdef _USE_PTHREADS
  tr->modelNumber = model;
  masterBarrier(THREAD_EVALUATE_PARTITION, tr); 
  for(i = 0, result = 0.0; i < NumberOfThreads; i++)
    result += reductionBuffer[i];      
#else
  result = evaluateIterativePartition(tr, lower, upper, model);
#endif
  tr->likelihood = result;       
 
  return result;
}


