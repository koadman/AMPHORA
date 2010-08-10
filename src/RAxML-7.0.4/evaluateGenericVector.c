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





static void evaluateGTRCAT_VECTOR (int *ex2, int *cptr,
				   double *x2_start, double *v, double *EIGN, double *rptr, double *tipVector, 
				   double pz, 
				   char *tipX1, int lower, int n, int numberOfCategories)
{
  double   z, lz, ki, lz1, lz2, lz3, term;    
  int     i; 
  double  *diagptable, *diagptable_start;
  double *x1, *x2;              
  
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
      v[i] = term;      
    }
  free(diagptable_start);        
}         
  
static void evaluateGTRCATMULT_VECTOR (int *ex2, int *cptr, int *modelptr, 
				       double *x2_start, double *v, double *EIGN, double *rptr, 
				       double *tipVector, double *pz, 
				       char *tipX1, int lower, int n, int numberOfCategories, int numberOfModels, int multiBranch)
{
  double   z, lz = 0.0, ki, lz1, lz2, lz3, term;    
  int     i;
  double  *diagptable, *diagptable_start; 
  double *x1, *x2;
  int model, modelCounter;        
  
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
      
      term = (log(term)) + (ex2[i] * log(minlikelihood));	   	    	   
	
      v[i] = term; 			  	    	 
    }
  free(diagptable_start); 
}
    
static void evaluateGTRCATPROT_VECTOR (int *ex2, int *cptr,
				       double *x2, double *v, double *EIGN, double *rptr, double *tipVector, double pz, 
				       char *tipX1, int lower, int n, int numberOfCategories)
{
  double   z, lz, ki, lza[19];       
  int     i, l; 
  double  *diagptable, *diagptable_start; 
  double term;
  double *left, *right;            
   
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
            
      term = log(term) + ex2[i] * log(minlikelihood);	          

      v[i] = term;	  
    }

  free(diagptable_start);          
}
      
   
static void evaluateGTRCATPROTMULT_VECTOR (int *ex2, int *cptr, int *modelptr, 
					   double *x2,  double *v, double *EIGN, double *rptr, double *tipVector, 
					   double *pz, 
					   char *tipX1, int lower, int n, int numberOfCategories, int numberOfModels, int multiBranch)
{
  double   z, lz = 0.0, ki, lza[19];       
  int     i, l; 
  double  *diagptable, *diagptable_start;   
  double term;
  double *left, *right;   
  int model, modelCounter;     
    
  if(!multiBranch)
    {
      z = pz[0];   
      if (z < zmin) z = zmin;
      lz = log(z);  
    }

  diagptable = diagptable_start = (double *)malloc(sizeof(double) * numberOfCategories * 19 * numberOfModels); 
 
  for(modelCounter = 0; modelCounter < numberOfModels; modelCounter++)
    {
      if(!multiBranch)
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
               

  for(i = lower; i < n; i++) 
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
      
      term = log(term) + ex2[i] * log(minlikelihood);	          
      
      v[i] = term;	   
    }
  
  free(diagptable_start);              
}         
      
static void evaluateGTRGAMMA_VECTOR (int *ex2,
				     double *x2_start, double *v, double *EIGN, double *gammaRates, 
				     double *tipVector, double pz, 
				     char *tipX1, int lower, int n)
{
  double   z, lz, term, ki;    
  int     i;  
  double  *diagptable, *diagptable_start;
  double *x1, *x2;  
  
  z = pz;
  
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
  
   
  for (i = lower; i < n; i++) 
    {
      x1 = &(tipVector[4 * tipX1[i]]);
      x2 = &x2_start[16 * i];
      diagptable = diagptable_start;
      
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
      
      term = log(0.25 * term) + ex2[i] * log(minlikelihood);
      v[i] = term;     
    }
  
  free(diagptable_start);   
} 

static void evaluateGTRGAMMAMULT_VECTOR(int *ex2, int *modelptr, 
					double *x2_start, double *v, double *EIGN, double *gammaRates, 
					double *tipVector, double *pz, 
					char *tipX1, int lower, int n, int numberOfModels, int multiBranch)
{
  double   z, lz = 0.0, term, ki;    
  int     i;
  double  *diagptable, *diagptable_start; 
  double *x1, *x2;
  int model;
  
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
      
      term = log(0.25 * term) + ex2[i] * log(minlikelihood);
      v[i] = term;          
    }
          
  free(diagptable_start); 
}
      
static void evaluateGTRGAMMAMULTINVAR_VECTOR(int *ex2, int *modelptr, int *iptr,
					     double *x2_start, double *v, double *EIGN, double *gammaRates, 
					     double *tipVector, double *tFreqs, double *invariants, double *pz, 
					     char *tipX1, int lower, int n, int numberOfModels, int multiBranch)
{
  double   z, lz = 0.0, term, ki;    
  int     i, model; 
  double  *diagptable, *diagptable_start;
  double *x1, *x2;  
  double *scalers;
  double *freqs;    
  
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
	term = log((scalers[model] * term + freqs[model * 4 + iptr[i]]) * pow(minlikelihood, ex2[i]));
      else
	term = log(scalers[model] * term) + ex2[i] * log(minlikelihood);
	
      v[i] = term;      
    }      
          
  free(diagptable_start); 
  free(scalers);
  free(freqs);
}
      
static void evaluateGTRGAMMAINVAR_VECTOR(int *ex2, int *iptr,
					 double *x2_start, double *v, double *EIGN, double *gammaRates, 
					 double *tipVector, double *tFreqs, double invariants, double pz, 
					 char *tipX1, int lower, int n)
{
  double   z, lz, term, ki;    
  int     i; 
  double  *diagptable;  
  double freqs[4];  
  double *x1, *x2;
  double scaler = 0.25 * (1.0 - invariants);   

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
	term = log((scaler * term + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
      else
	term = log(scaler * term) + ex2[i] * log(minlikelihood);	 
      
      v[i] = term;		  	   	
    }
          
  free(diagptable);            
}
 
static void evaluateGTRGAMMAPROT_VECTOR(int *ex2, 
					double *x2,  double *v, double *EIGN, double *gammaRates, 
					double *tipVector, double pz, 
					char *tipX1, int lower, int n)
{
  double   z, lz, term, ki;    
  int     i, j;
  double  *diagptable, *diagptable_start;
  double  *left, *right;  
          
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
      
      term = log(0.25 * term) + ex2[i] * log(minlikelihood);	   
      v[i] = term;		  	   	 
    }
    
  free(diagptable_start); 
}
    
static void evaluateGTRGAMMAPROTINVAR_VECTOR(int *ex2, int *iptr,
					     double *x2,  double *v, double *EIGN, double *gammaRates, 
					     double *tipVector,double *tFreqs, double invariants, double pz, 
					     char *tipX1, int lower, int n)
{
  double   z, lz, term, ki;    
  int     i, j;    
  double  *diagptable, *diagptable_start;
  double *left, *right;
  double scaler = 0.25 * (1.0 - invariants);
  double freqs[20]; 
  
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
	term = log((scaler * term + freqs[iptr[i]]) * pow(minlikelihood, ex2[i]));
      else
	term = log(scaler * term) + ex2[i] * log(minlikelihood);
      
      v[i] = term;		  	   	 
    }
    

  free(diagptable_start); 
}
    	 
static void evaluateGTRGAMMAPROTMULT_VECTOR (int *ex2, int *modelptr,
					     double *x2,  double *v, double *EIGN, double *gammaRates, 
					     double *tipVector, double *pz, 
					     char *tipX1, int lower, int n, int numberOfModels, int multiBranch)
{
  double   z, lz = 0.0, term, ki;    
  int     i, j;  
  double  *diagptable, *diagptable_start;
  double *left, *right;
  int model;
      
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
      
      term = log(0.25 * term) + ex2[i] * log(minlikelihood);	   
      v[i] = term;		  	   	 
    }
    
  free(diagptable_start); 
}
	
static void evaluateGTRGAMMAPROTMULTINVAR_VECTOR(int *ex2, int *modelptr, int *iptr,
						 double *x2,  double *v, double *EIGN, double *gammaRates, 
						 double *tipVector, double *tFreqs, double *invariants, double *pz, 
						 char *tipX1, int lower, int n, int numberOfModels, int multiBranch)
{
  double   z, lz = 0.0, term, ki;    
  int     i, j;  
  double  *diagptable, *diagptable_start;
  double *left, *right;
  int model;
  double *scalers;
  double *freqs;  
      
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
	term = log((scalers[model] * term + freqs[model * 20 + iptr[i]]) * pow(minlikelihood, ex2[i]));
      else
	term = log(scalers[model] * term) + ex2[i] * log(minlikelihood);
      
      v[i] = term;		  	   	 
    }
    

  free(diagptable_start); 
  free(freqs);
  free(scalers);
}
       
    

#ifdef _LOCAL_DATA
  
void evaluateGenericVectorIterative(tree *localTree, int startIndex, int endIndex)    
{
  double *x2_start, *pz, *v;
  char *tip;
  int *ex2;
  int pNumber, qNumber;

  v = localTree->strided_siteLL_Vector;

  pNumber = localTree->td[0].ti[0].pNumber;
  qNumber = localTree->td[0].ti[0].qNumber;
  pz      = localTree->td[0].ti[0].qz;

  newviewIterative(localTree, startIndex, endIndex);


  x2_start = getLikelihoodArray(qNumber,  localTree->mxtips, localTree->xVector); 
  ex2      = getScalingArray(qNumber, localTree->cdta->endsite, localTree->mxtips, localTree->expArray); 

  tip   = localTree->strided_yVector[pNumber];


  if(localTree->mixedData)
    {
      assert(0);
    }
  else
    {
      switch(localTree->likelihoodFunction)
	{
	case GTRCAT:
	  evaluateGTRCAT_VECTOR(ex2, localTree->strided_rateCategory, 
				x2_start, v, localTree->EIGN_DNA, localTree->strided_patrat, localTree->tipVectorDNA, pz[0], 
				tip, startIndex, endIndex, localTree->NumberOfCategories);    
	  break;
	case GTRCATMULT:
	  evaluateGTRCATMULT_VECTOR(ex2, localTree->strided_rateCategory, localTree->strided_model,
				    x2_start, v, localTree->EIGN_DNA, localTree->strided_patrat, localTree->tipVectorDNA, pz, 
				    tip, startIndex, endIndex, localTree->NumberOfCategories, localTree->NumberOfModels, 
				    localTree->multiBranch);
	  break;
	case PROTCAT:
	  evaluateGTRCATPROT_VECTOR(ex2, localTree->strided_rateCategory,
				    x2_start, v, localTree->EIGN_AA, localTree->strided_patrat, localTree->tipVectorAA, pz[0], 
				    tip, startIndex, endIndex, localTree->NumberOfCategories);
	  break;
	case PROTCATMULT:
	  evaluateGTRCATPROTMULT_VECTOR(ex2, localTree->strided_rateCategory, localTree->strided_model, 
					x2_start, v, localTree->EIGN_AA, localTree->strided_patrat, localTree->tipVectorAA, pz, 
					tip, startIndex, endIndex, localTree->NumberOfCategories, localTree->NumberOfModels, localTree->multiBranch);
	  break;
	case GTRGAMMA:
	  evaluateGTRGAMMA_VECTOR(ex2, 
				  x2_start, v, localTree->EIGN_DNA, localTree->gammaRates, localTree->tipVectorDNA, pz[0], 
				  tip, startIndex, endIndex); 
	  break;
	case GTRGAMMAI:
	  evaluateGTRGAMMAINVAR_VECTOR(ex2, localTree->strided_invariant,
				       x2_start, v, localTree->EIGN_DNA, localTree->gammaRates, 
				       localTree->tipVectorDNA, localTree->frequencies_DNA, localTree->invariants[0], pz[0], 
				       tip, startIndex, endIndex);
	  break;
	case GTRGAMMAMULT:
	  evaluateGTRGAMMAMULT_VECTOR(ex2, localTree->strided_model,
				      x2_start, v, localTree->EIGN_DNA, localTree->gammaRates, localTree->tipVectorDNA, pz, 
				      tip, startIndex, endIndex, localTree->NumberOfModels, localTree->multiBranch);
	  break;
	case GTRGAMMAMULTI:
	  evaluateGTRGAMMAMULTINVAR_VECTOR(ex2, localTree->strided_model,  localTree->strided_invariant,
					   x2_start, v, localTree->EIGN_DNA, localTree->gammaRates, 
					   localTree->tipVectorDNA, localTree->frequencies_DNA, localTree->invariants, pz, 
					   tip, startIndex, endIndex, localTree->NumberOfModels, localTree->multiBranch); 
	  break;
	case PROTGAMMA:
	  evaluateGTRGAMMAPROT_VECTOR(ex2,
				      x2_start, v, localTree->EIGN_AA, localTree->gammaRates, localTree->tipVectorAA, pz[0], 
				      tip, startIndex, endIndex);
	  break;
	case PROTGAMMAI:
	  evaluateGTRGAMMAPROTINVAR_VECTOR(ex2,  localTree->strided_invariant,
					   x2_start, v, localTree->EIGN_AA, localTree->gammaRates, 
					   localTree->tipVectorAA, localTree->frequencies_AA, localTree->invariants[0], pz[0], 
					   tip, startIndex, endIndex);
	  break;
	case PROTGAMMAMULT:
	  evaluateGTRGAMMAPROTMULT_VECTOR(ex2, localTree->strided_model,
					  x2_start, v, localTree->EIGN_AA, localTree->gammaRates, localTree->tipVectorAA, pz, 
					  tip, startIndex, endIndex, localTree->NumberOfModels, localTree->multiBranch);
	  break;
	case PROTGAMMAMULTI:
	  evaluateGTRGAMMAPROTMULTINVAR_VECTOR(ex2, localTree->strided_model,  localTree->strided_invariant,
					       x2_start, v, localTree->EIGN_AA, localTree->gammaRates, 
					       localTree->tipVectorAA, localTree->frequencies_AA, localTree->invariants, pz, 
					       tip, startIndex, endIndex, localTree->NumberOfModels, localTree->multiBranch
					       );
	  break;
	default:
	  assert(0);
	} 
    }

}
  
#else

void evaluateGenericVectorIterative(tree *tr, int startIndex, int endIndex)    
{
  double *x2_start, *pz, *v;
  char *tip;
  int *ex2;
  int pNumber, qNumber;

  v = tr->siteLL_Vector;

  pNumber = tr->td[0].ti[0].pNumber;
  qNumber = tr->td[0].ti[0].qNumber;
  pz      = tr->td[0].ti[0].qz;

  newviewIterative(tr, 0, tr->cdta->endsite);


  x2_start = getLikelihoodArray(qNumber,  tr->mxtips, tr->xVector); 
  ex2      = getScalingArray(qNumber, tr->cdta->endsite, tr->mxtips, tr->expArray); 

  tip   = tr->yVector[pNumber];


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
		case CAT:
		  evaluateGTRCAT_VECTOR(&(ex2[l]), &(tr->cdta->rateCategory[l]), 
					&(x2_start[offset]), &(v[l]), &(tr->EIGN_DNA[model * 3]), 
					tr->cdta->patrat, &(tr->tipVectorDNA[64 * model]), pz[branchIndex], 
					&(tip[l]), 0, width, tr->NumberOfCategories);
		  break;
		case GAMMA:		  
		  evaluateGTRGAMMA_VECTOR(&(ex2[l]), 
					  &(x2_start[offset]), &(v[l]), 
					  &(tr->EIGN_DNA[model * 3]), 
					  &(tr->gammaRates[model * 4]), 
					  &(tr->tipVectorDNA[model * 64]), pz[branchIndex], 
					  &(tip[l]), 0, width); 
		  break;
		case GAMMA_I:
		  evaluateGTRGAMMAINVAR_VECTOR(&(ex2[l]), &(tr->invariant[l]),
					       &(x2_start[offset]), &(v[l]), 
					       &(tr->EIGN_DNA[model * 3]), 
					       &(tr->gammaRates[model * 4]), 
					       &(tr->tipVectorDNA[model * 64]), 
					       &(tr->frequencies_DNA[model * 4]), 
					       tr->invariants[model], pz[branchIndex], 
					       &(tip[l]), 0, width);
		  break;
		default:
		  assert(0);
		}
	      break;
	    case AA_DATA:
	      switch(tr->rateHetModel)
		{
		case CAT:
		  evaluateGTRCATPROT_VECTOR(&(ex2[l]), &(tr->cdta->rateCategory[l]),
					    &(x2_start[offset]), &(v[l]),
					    &(tr->EIGN_AA[model * 19]), 
					    tr->cdta->patrat, &(tr->tipVectorAA[model * 460]), pz[branchIndex], 
					    &(tip[l]), 0, width, tr->NumberOfCategories);
		  break;
		case GAMMA:		  
		  evaluateGTRGAMMAPROT_VECTOR(&(ex2[l]),
					      &(x2_start[offset]), &(v[l]), 
					      &(tr->EIGN_AA[model * 19]), 
					      &(tr->gammaRates[model * 4]), 
					      &(tr->tipVectorAA[model * 460]), pz[branchIndex], 
					      &(tip[l]), 0, width);
		  break;
		case GAMMA_I:		 
		  evaluateGTRGAMMAPROTINVAR_VECTOR(&(ex2[l]),  &(tr->invariant[l]),
						   &(x2_start[offset]), &(v[l]), 
						   &(tr->EIGN_AA[model * 19]), 
						   &(tr->gammaRates[model * 4]), 
						   &(tr->tipVectorAA[model * 460]), 
						   &(tr->frequencies_AA[model * 20]) , 
						   tr->invariants[model], pz[branchIndex], 
						   &(tip[l]), 0, width);
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
	case GTRCAT:
	  evaluateGTRCAT_VECTOR(ex2, tr->cdta->rateCategory, 
				x2_start, v, tr->EIGN_DNA, tr->cdta->patrat, tr->tipVectorDNA, pz[0], 
				tip, startIndex, endIndex, tr->NumberOfCategories);    
	  break;
	case GTRCATMULT:
	  evaluateGTRCATMULT_VECTOR(ex2, tr->cdta->rateCategory, tr->model,
				    x2_start, v, tr->EIGN_DNA, tr->cdta->patrat, tr->tipVectorDNA, pz, 
				    tip, startIndex, endIndex, tr->NumberOfCategories, tr->NumberOfModels, tr->multiBranch);
	  break;
	case PROTCAT:
	  evaluateGTRCATPROT_VECTOR(ex2, tr->cdta->rateCategory,
				    x2_start, v, tr->EIGN_AA, tr->cdta->patrat, tr->tipVectorAA, pz[0], 
				    tip, startIndex, endIndex, tr->NumberOfCategories);
	  break;
	case PROTCATMULT:
	  evaluateGTRCATPROTMULT_VECTOR(ex2, tr->cdta->rateCategory, tr->model, 
					x2_start, v, tr->EIGN_AA, tr->cdta->patrat, tr->tipVectorAA, pz, 
					tip, startIndex, endIndex, tr->NumberOfCategories, tr->NumberOfModels, tr->multiBranch);
	  break;
	case GTRGAMMA:
	  evaluateGTRGAMMA_VECTOR(ex2, 
				  x2_start, v, tr->EIGN_DNA, tr->gammaRates, tr->tipVectorDNA, pz[0], 
				  tip, startIndex, endIndex); 
	  break;
	case GTRGAMMAI:
	  evaluateGTRGAMMAINVAR_VECTOR(ex2, tr->invariant,
				       x2_start, v, tr->EIGN_DNA, tr->gammaRates, 
				       tr->tipVectorDNA, tr->frequencies_DNA, tr->invariants[0], pz[0], 
				       tip, startIndex, endIndex);
	  break;
	case GTRGAMMAMULT:
	  evaluateGTRGAMMAMULT_VECTOR(ex2, tr->model,
				      x2_start, v, tr->EIGN_DNA, tr->gammaRates, tr->tipVectorDNA, pz, 
				      tip, startIndex, endIndex, tr->NumberOfModels, tr->multiBranch);
	  break;
	case GTRGAMMAMULTI:
	  evaluateGTRGAMMAMULTINVAR_VECTOR(ex2, tr->model,  tr->invariant,
					   x2_start, v, tr->EIGN_DNA, tr->gammaRates, 
					   tr->tipVectorDNA, tr->frequencies_DNA, tr->invariants, pz, 
					   tip, startIndex, endIndex, tr->NumberOfModels, tr->multiBranch); 
	  break;
	case PROTGAMMA:
	  evaluateGTRGAMMAPROT_VECTOR(ex2,
				      x2_start, v, tr->EIGN_AA, tr->gammaRates, tr->tipVectorAA, pz[0], 
				      tip, startIndex, endIndex);
	  break;
	case PROTGAMMAI:
	  evaluateGTRGAMMAPROTINVAR_VECTOR(ex2,  tr->invariant,
					   x2_start, v, tr->EIGN_AA, tr->gammaRates, 
					   tr->tipVectorAA, tr->frequencies_AA, tr->invariants[0], pz[0], 
					   tip, startIndex, endIndex);
	  break;
	case PROTGAMMAMULT:
	  evaluateGTRGAMMAPROTMULT_VECTOR(ex2, tr->model,
					  x2_start, v, tr->EIGN_AA, tr->gammaRates, tr->tipVectorAA, pz, 
					  tip, startIndex, endIndex, tr->NumberOfModels, tr->multiBranch);
	  break;
	case PROTGAMMAMULTI:
	  evaluateGTRGAMMAPROTMULTINVAR_VECTOR(ex2, tr->model,  tr->invariant,
					       x2_start, v, tr->EIGN_AA, tr->gammaRates, 
					       tr->tipVectorAA, tr->frequencies_AA, tr->invariants, pz, 
					       tip, startIndex, endIndex, tr->NumberOfModels, tr->multiBranch
					       );
	  break;
	default:
	  assert(0);
	} 
    }

}

#endif

void evaluateGenericVector (tree *tr, nodeptr p, double *v)
{
  nodeptr q = p->back;
  int i;  
    
  assert(isTip(p->number, tr->rdta->numsp)  || isTip(q->number, tr->rdta->numsp));
      
  if(isTip(q->number, tr->rdta->numsp))
    {
      nodeptr tmp = p;
      p = q;
      q = tmp;	      
    }   
  
  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  for(i = 0; i < tr->numBranches; i++)    
    tr->td[0].ti[0].qz[i] =  q->z[i];

  tr->siteLL_Vector = v;

  tr->td[0].count = 1;
  if(!q->x)
    computeTraversalInfo(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp, tr->numBranches); 
 
#ifdef _USE_PTHREADS
  masterBarrier(THREAD_EVALUATE_VECTOR, tr);
#else
  evaluateGenericVectorIterative(tr, 0, tr->cdta->endsite); 
#endif
}
