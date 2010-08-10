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
#include <pthread.h>
extern int NumberOfThreads;
extern pthread_mutex_t jobMutex;
extern pthread_cond_t  jobCond;
#endif



static void newviewGTRCAT( traversalInfo *ti,  double *EV,  double *EI,  double *EIGN, 
			      double *rptr,  int *cptr, 
                              double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
			      int    *ex1,  int *ex2,  int *ex3, char *tipX1, char *tipX2,
			      int lower, int n,  int numberOfCategories, double z1, double z2)
{         
  double  
    *left, *left_start,
    *x1, *x2, *x3;
  double
    d1c, d1g, d1t, d2c, d2g, d2t,
    ump_x1_1, ump_x1_2, ump_x1_3, ump_x1_0, 
    ump_x2_0, ump_x2_1, ump_x2_2, ump_x2_3, x1px2, ki,
    lz10, lz11, lz12, lz20, lz21, lz22;  
  int i;  
  
  left_start = left = (double *)malloc(24 * numberOfCategories * sizeof(double));                           
  
  lz10 = EIGN[0] * z1;
  lz11 = EIGN[1] * z1;
  lz12 = EIGN[2] * z1;
  
  lz20 = EIGN[0] * z2;
  lz21 = EIGN[1] * z2;
  lz22 = EIGN[2] * z2;
  
  for(i = 0; i < numberOfCategories; i++)
    {		   
      ki = rptr[i];	             
      
      d1c = exp(ki * lz10);
      d1g = exp(ki * lz11);
      d1t = exp(ki * lz12);	  
      
      *left++ = d1c * EI[0];
      *left++ = d1g * EI[1];
      *left++ = d1t * EI[2];
      
      *left++ = d1c * EI[3];
      *left++ = d1g * EI[4];
      *left++ = d1t * EI[5];
      
      *left++ = d1c * EI[6];
      *left++ = d1g * EI[7];
      *left++ = d1t * EI[8];
      
      *left++ = d1c * EI[9];
      *left++ = d1g * EI[10];
      *left++ = d1t * EI[11];
            
      d2c = exp(ki * lz20);
      d2g = exp(ki * lz21);
      d2t = exp(ki * lz22);	
      
      *left++ = d2c * EI[0];
      *left++ = d2g * EI[1];
      *left++ = d2t * EI[2];
      
      *left++ = d2c * EI[3];
      *left++ = d2g * EI[4];
      *left++ = d2t * EI[5];
      
      *left++ = d2c * EI[6];
      *left++ = d2g * EI[7];
      *left++ = d2t * EI[8];
      
      *left++ = d2c * EI[9];
      *left++ = d2g * EI[10];
      *left++ = d2t * EI[11];
    }         
        


  switch(ti->tipCase)
    {
    case TIP_TIP:
      {	  	      	         	      	

      for (i = lower; i < n; i++) 
	{		   	    	   	   	     	   
	  x1 = &(tipVector[4 * tipX1[i]]);
	  x2 = &(tipVector[4 * tipX2[i]]);	  
	  x3 = &x3_start[4 * i];	 
	  
	  left = &left_start[cptr[i] * 24];
	 
	  ump_x1_0 =  x1[0];
	  ump_x1_0 += x1[1] * *left++;
	  ump_x1_0 += x1[2] * *left++;
	  ump_x1_0 += x1[3] * *left++;	      		   	    
	  
	  ump_x1_1 =  x1[0];
	  ump_x1_1 += x1[1] * *left++;
	  ump_x1_1 += x1[2] * *left++;
	  ump_x1_1 += x1[3] * *left++;	      		  	 
	  
	  ump_x1_2 =  x1[0];
	  ump_x1_2 += x1[1] * *left++;
	  ump_x1_2 += x1[2] * *left++;
	  ump_x1_2 += x1[3] * *left++;	      		   	 
	  
	  ump_x1_3 =  x1[0];
	  ump_x1_3 += x1[1] * *left++;
	  ump_x1_3 += x1[2] * *left++;
	  ump_x1_3 += x1[3] * *left++;	      		   	  
	  
	  ump_x2_0 =  x2[0];
	  ump_x2_0 += x2[1] * *left++;
	  ump_x2_0 += x2[2] * *left++;
	  ump_x2_0 += x2[3] * *left++;		     	  	
	  
	  ump_x2_1 =  x2[0];
	  ump_x2_1 += x2[1] * *left++;
	  ump_x2_1 += x2[2] * *left++;
	  ump_x2_1 += x2[3] * *left++;		     	 	     
	  
	  ump_x2_2 =  x2[0];
	  ump_x2_2 += x2[1] * *left++;
	  ump_x2_2 += x2[2] * *left++;
	  ump_x2_2 += x2[3] * *left++;		     	   	      
	  
	  ump_x2_3 =  x2[0];
	  ump_x2_3 += x2[1] * *left++;
	  ump_x2_3 += x2[2] * *left++;
	  ump_x2_3 += x2[3] * *left++;		     	  	       	   
	  
	  x1px2 = ump_x1_0 * ump_x2_0;	      
	  x3[0] = x1px2 *  EV[0];
	  x3[1] = x1px2 *  EV[1];
	  x3[2] = x1px2 *  EV[2];
	  x3[3] = x1px2 *  EV[3];
	  
	  x1px2 = ump_x1_1 * ump_x2_1;	       
	  x3[0] += x1px2 *  EV[4];
	  x3[1] += x1px2 *  EV[5];
	  x3[2] += x1px2 *  EV[6];
	  x3[3] += x1px2 *  EV[7];
	  
	  x1px2 = ump_x1_2 * ump_x2_2;	       
	  x3[0] += x1px2 *  EV[8];
	  x3[1] += x1px2 *  EV[9];
	  x3[2] += x1px2 *  EV[10];
	  x3[3] += x1px2 *  EV[11];
	  
	  x1px2 = ump_x1_3 * ump_x2_3;	      
	  x3[0] += x1px2 *  EV[12];
	  x3[1] += x1px2 *  EV[13];
	  x3[2] += x1px2 *  EV[14];
	  x3[3] += x1px2 *  EV[15];	      	 
	  
	  ex3[i] = 0;		  	    	   	    	      
	}          
      }
      break;
    case TIP_INNER:
      {		                                

	for (i = lower; i < n; i++) 
	  {		     		      	      
	    x1 = &(tipVector[4 * tipX1[i]]);  
	    x2 = &x2_start[4 * i];
	    x3 = &x3_start[4 * i];    
	    left = &left_start[cptr[i] * 24];	    
	    
	    ump_x1_0 = x1[0];
	    ump_x1_0 += x1[1] * *left++;
	    ump_x1_0 += x1[2] * *left++;
	    ump_x1_0 += x1[3] * *left++;	      		   	    
	    
	    ump_x1_1 =  x1[0];
	    ump_x1_1 += x1[1] * *left++;
	    ump_x1_1 += x1[2] * *left++;
	    ump_x1_1 += x1[3] * *left++;	      		  	 
	    
	    ump_x1_2 =  x1[0];
	    ump_x1_2 += x1[1] * *left++;
	    ump_x1_2 += x1[2] * *left++;
	    ump_x1_2 += x1[3] * *left++;	      		   	 
	    
	    ump_x1_3 =  x1[0];
	    ump_x1_3 += x1[1] * *left++;
	    ump_x1_3 += x1[2] * *left++;
	    ump_x1_3 += x1[3] * *left++;	      		   	  
	    
	    ump_x2_0 =  x2[0];
	    ump_x2_0 += x2[1] * *left++;
	    ump_x2_0 += x2[2] * *left++;
	    ump_x2_0 += x2[3] * *left++;		     	  	
	    
	    ump_x2_1 =  x2[0];
	    ump_x2_1 += x2[1] * *left++;
	    ump_x2_1 += x2[2] * *left++;
	    ump_x2_1 += x2[3] * *left++;		     	 	     
	    
	    ump_x2_2 =  x2[0];
	    ump_x2_2 += x2[1] * *left++;
	    ump_x2_2 += x2[2] * *left++;
	    ump_x2_2 += x2[3] * *left++;		     	   	      
	    
	    ump_x2_3 =  x2[0];
	    ump_x2_3 += x2[1] * *left++;
	    ump_x2_3 += x2[2] * *left++;
	    ump_x2_3 += x2[3] * *left++;
	    
	    x1px2 = ump_x1_0 * ump_x2_0;
	    x3[0] = x1px2 *  EV[0];
	    x3[1] = x1px2 *  EV[1];
	    x3[2] = x1px2 *  EV[2];
	    x3[3] = x1px2 *  EV[3];
	    
	    x1px2 = ump_x1_1 * ump_x2_1;
	    x3[0] += x1px2 *  EV[4];
	    x3[1] += x1px2 *  EV[5];
	    x3[2] += x1px2 *  EV[6];
	    x3[3] += x1px2 *  EV[7];
	    
	    x1px2 = ump_x1_2 * ump_x2_2;
	    x3[0] += x1px2 *  EV[8];
	    x3[1] += x1px2 *  EV[9];
	    x3[2] += x1px2 *  EV[10];
	    x3[3] += x1px2 *  EV[11];
	    
	    x1px2 = ump_x1_3 * ump_x2_3;
	    x3[0] += x1px2 *   EV[12];
	    x3[1] += x1px2 *   EV[13];
	    x3[2] += x1px2 *   EV[14];
	    x3[3] += x1px2 *   EV[15];
	    
	    ex3[i] = ex2[i];		  	      	       	       
	    
	    if (x3[0] < minlikelihood && x3[0] > minusminlikelihood &&
		x3[1] < minlikelihood && x3[1] > minusminlikelihood &&
		x3[2] < minlikelihood && x3[2] > minusminlikelihood &&
		x3[3] < minlikelihood && x3[3] > minusminlikelihood)	       
	      {	       	      
		x3[0]   *= twotothe256;
		x3[1]   *= twotothe256;
		x3[2]   *= twotothe256;		
		x3[3]   *= twotothe256;
		ex3[i]  += 1;			   
	      }	      	      	
	  }     	            
      }
      break;
    case INNER_INNER:     
      for (i = lower; i < n; i++) 
	{		     	    
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  x3 = &x3_start[4 * i];
	  left = &left_start[cptr[i] * 24];
	  
	  ump_x1_0 =  x1[0];
	  ump_x1_0 += x1[1] * *left++;
	  ump_x1_0 += x1[2] * *left++;
	  ump_x1_0 += x1[3] * *left++;	      		   	    
	  
	  ump_x1_1 =  x1[0];
	  ump_x1_1 += x1[1] * *left++;
	  ump_x1_1 += x1[2] * *left++;
	  ump_x1_1 += x1[3] * *left++;	      		  	 
	  
	  ump_x1_2 =  x1[0];
	  ump_x1_2 += x1[1] * *left++;
	  ump_x1_2 += x1[2] * *left++;
	  ump_x1_2 += x1[3] * *left++;	      		   	 
	  
	  ump_x1_3 =  x1[0];
	  ump_x1_3 += x1[1] * *left++;
	  ump_x1_3 += x1[2] * *left++;
	  ump_x1_3 += x1[3] * *left++;	      		   	  
	  
	  ump_x2_0 =  x2[0];
	  ump_x2_0 += x2[1] * *left++;
	  ump_x2_0 += x2[2] * *left++;
	  ump_x2_0 += x2[3] * *left++;		     	  	
	  
	  ump_x2_1 =  x2[0];
	  ump_x2_1 += x2[1] * *left++;
	  ump_x2_1 += x2[2] * *left++;
	  ump_x2_1 += x2[3] * *left++;		     	 	     
	  
	  ump_x2_2 =  x2[0];
	  ump_x2_2 += x2[1] * *left++;
	  ump_x2_2 += x2[2] * *left++;
	  ump_x2_2 += x2[3] * *left++;		     	   	      
	  
	  ump_x2_3 =  x2[0];
	  ump_x2_3 += x2[1] * *left++;
	  ump_x2_3 += x2[2] * *left++;
	  ump_x2_3 += x2[3] * *left++;
	  
	  x1px2 = ump_x1_0 * ump_x2_0;
	  x3[0] = x1px2 *  EV[0];
	  x3[1] = x1px2 *  EV[1];
	  x3[2] = x1px2 *  EV[2];
	  x3[3] = x1px2 *  EV[3];
	  
	  x1px2 = ump_x1_1 * ump_x2_1;
	  x3[0] += x1px2 *  EV[4];
	  x3[1] += x1px2 *  EV[5];
	  x3[2] += x1px2 *  EV[6];
	  x3[3] += x1px2 *  EV[7];
	  
	  x1px2 = ump_x1_2 * ump_x2_2;
	  x3[0] += x1px2 *  EV[8];
	  x3[1] += x1px2 *  EV[9];
	  x3[2] += x1px2 *  EV[10];
	  x3[3] += x1px2 *  EV[11];
	  
	  x1px2 = ump_x1_3 * ump_x2_3;
	  x3[0] += x1px2 *   EV[12];
	  x3[1] += x1px2 *   EV[13];
	  x3[2] += x1px2 *   EV[14];
	  x3[3] += x1px2 *   EV[15];
	  
	  ex3[i] = ex1[i] + ex2[i];	
	  
	  if (x3[0] < minlikelihood && x3[0] > minusminlikelihood &&
	      x3[1] < minlikelihood && x3[1] > minusminlikelihood &&
	      x3[2] < minlikelihood && x3[2] > minusminlikelihood &&
	      x3[3] < minlikelihood && x3[3] > minusminlikelihood)  
	    {	 	     
	      x3[0]   *= twotothe256;
	      x3[1]   *= twotothe256;
	      x3[2]   *= twotothe256;		
	      x3[3]   *= twotothe256;
	      ex3[i] += 1;
	    }	      	          	    
	}     	  
      break;
    default:
      assert(0);
    }

  free(left_start); 	   
}
  


static void newviewGTRCATMULT(traversalInfo *ti,  double *extEV,  double *extEI,  double *EIGN, 
			      double *rptr,  int *cptr, 
			      double *x1_start,  double *x2_start,  double *x3_start,  double *tipVector,
			      int    *ex1,  int *ex2,  int *ex3, char *tipX1, char *tipX2, int *modelptr,
			      int lower, int n,  int numberOfCategories, int numberOfModels, int multiBranch
			      )
{       
  double  
    *left, *left_start,
    *EI, *EV, *x1, *x2, *x3;
  double      
    ump_x1_1, ump_x1_2, ump_x1_3, ump_x1_0, 
    ump_x2_0, ump_x2_1, ump_x2_2, ump_x2_3, x1px2, z1 = 0.0, z2 = 0.0, ki,
    d1c, d1g, d1t, d2c, d2g, d2t,
    lz10, lz11, lz12, lz20, lz21, lz22;
  int 
    i, modelCounter, model;	             
                             
  if(!multiBranch)
    {
      z1  = ti->qz[0];      	  
      z2  = ti->rz[0];
    }
      
  left_start = left = (double *)malloc(24 * numberOfModels * numberOfCategories * sizeof(double));              
  
  for(modelCounter = 0; modelCounter < numberOfModels; modelCounter++)
    {
      if(multiBranch)
	{
	  z1  = ti->qz[modelCounter];	     	      
	  z2  = ti->rz[modelCounter];	     
	}
      
      lz10 = EIGN[modelCounter * 3] * z1;
      lz11 = EIGN[modelCounter * 3 + 1] * z1;
      lz12 = EIGN[modelCounter * 3 + 2] * z1;
      
      lz20 = EIGN[modelCounter * 3] * z2;
      lz21 = EIGN[modelCounter * 3 + 1] * z2;
      lz22 = EIGN[modelCounter * 3 + 2] * z2;        	  	 
      
      EI = &(extEI[modelCounter * 12]);
      
      for(i = 0; i < numberOfCategories; i++)
	{	
	  ki = rptr[i];
	  
	  d1c = exp (ki * lz10);
	  d1g = exp (ki * lz11);
	  d1t = exp (ki * lz12);	
	  
	  *left++ = d1c * EI[0];
	  *left++ = d1g * EI[1];
	  *left++ = d1t * EI[2];
	  
	  *left++ = d1c * EI[3];
	  *left++ = d1g * EI[4];
	  *left++ = d1t * EI[5];
	  
	  *left++ = d1c * EI[6];
	  *left++ = d1g * EI[7];
	  *left++ = d1t * EI[8];
	  
	  *left++ = d1c * EI[9];
	  *left++ = d1g * EI[10];
	  *left++ = d1t * EI[11];
	  
	  d2c = exp (ki * lz20);
	  d2g = exp (ki * lz21);
	  d2t = exp (ki * lz22);	
	  
	  *left++ = d2c * EI[0];
	  *left++ = d2g * EI[1];
	  *left++ = d2t * EI[2];
	  
	  *left++ = d2c * EI[3];
	  *left++ = d2g * EI[4];
	  *left++ = d2t * EI[5];
	  
	  *left++ = d2c * EI[6];
	  *left++ = d2g * EI[7];
	  *left++ = d2t * EI[8];
	  
	  *left++ = d2c * EI[9];
	  *left++ = d2g * EI[10];
	  *left++ = d2t * EI[11];
	}
    }
  
  switch(ti->tipCase)
    {
    case TIP_TIP:
      {	 	 	 

	for (i = lower; i < n; i++) 
	  {		   	    	   	
	    model = modelptr[i];	   
	    EV = &(extEV[model * 16]);
	    x1 = &(tipVector[model * 64 + 4 * tipX1[i]]);
	    x2 = &(tipVector[model * 64 + 4 * tipX2[i]]);
	    
	    x3 = &x3_start[4 * i];
	    left = &left_start[model * 24 * numberOfCategories + cptr[i] * 24];
	    
	    ump_x1_0 =  x1[0];
	    ump_x1_0 += x1[1] * *left++;
	    ump_x1_0 += x1[2] * *left++;
	    ump_x1_0 += x1[3] * *left++;	      		   	    
	    
	    ump_x1_1 =  x1[0];
	    ump_x1_1 += x1[1] * *left++;
	    ump_x1_1 += x1[2] * *left++;
	    ump_x1_1 += x1[3] * *left++;	      		  	 
	    
	    ump_x1_2 =  x1[0];
	    ump_x1_2 += x1[1] * *left++;
	    ump_x1_2 += x1[2] * *left++;
	    ump_x1_2 += x1[3] * *left++;	      		   	 
	    
	    ump_x1_3 =  x1[0];
	    ump_x1_3 += x1[1] * *left++;
	    ump_x1_3 += x1[2] * *left++;
	    ump_x1_3 += x1[3] * *left++;	      		   	  
	    
	    ump_x2_0 =  x2[0];
	    ump_x2_0 += x2[1] * *left++;
	    ump_x2_0 += x2[2] * *left++;
	    ump_x2_0 += x2[3] * *left++;		     	  	
	    
	    ump_x2_1 =  x2[0];
	    ump_x2_1 += x2[1] * *left++;
	    ump_x2_1 += x2[2] * *left++;
	    ump_x2_1 += x2[3] * *left++;		     	 	     
	    
	    ump_x2_2 =  x2[0];
	    ump_x2_2 += x2[1] * *left++;
	    ump_x2_2 += x2[2] * *left++;
	    ump_x2_2 += x2[3] * *left++;		     	   	      
	    
	    ump_x2_3 =  x2[0];
	    ump_x2_3 += x2[1] * *left++;
	    ump_x2_3 += x2[2] * *left++;
	    ump_x2_3 += x2[3] * *left++;		     	  	       	   
	    
	    x1px2 = ump_x1_0 * ump_x2_0;	      
	    x3[0] = x1px2 *  EV[0];
	    x3[1] = x1px2 *  EV[1];
	    x3[2] = x1px2 *  EV[2];
	    x3[3] = x1px2 *  EV[3];
	    
	    x1px2 = ump_x1_1 * ump_x2_1;	       
	    x3[0] += x1px2 *  EV[4];
	    x3[1] += x1px2 *  EV[5];
	    x3[2] += x1px2 *  EV[6];
	    x3[3] += x1px2 *  EV[7];
	    
	    x1px2 = ump_x1_2 * ump_x2_2;	       
	    x3[0] += x1px2 *  EV[8];
	    x3[1] += x1px2 *  EV[9];
	    x3[2] += x1px2 *  EV[10];
	    x3[3] += x1px2 *  EV[11];
	    
	    x1px2 = ump_x1_3 * ump_x2_3;	      
	    x3[0] += x1px2 *  EV[12];
	    x3[1] += x1px2 *  EV[13];
	    x3[2] += x1px2 *  EV[14];
	    x3[3] += x1px2 *  EV[15];	      
	    
	    ex3[i] = 0;	
	  }
      }
      break;
    case TIP_INNER:
      {		
  
	for (i = lower; i < n; i++) 
	  {		     	
	    model = modelptr[i];	   
	    EV = &(extEV[model * 16]);
	    x1 = &(tipVector[model * 64 + 4 * tipX1[i]]);     
	    x2 = &x2_start[4 * i];
	    x3 = &x3_start[4 * i];	    
	   
	    left = &left_start[model * 24 * numberOfCategories + cptr[i] * 24];	    
  	    	     
	    ump_x1_0 = x1[0];
	    ump_x1_0 += x1[1] * *left++;
	    ump_x1_0 += x1[2] * *left++;
	    ump_x1_0 += x1[3] * *left++;	      		   	    
	    
	    ump_x1_1 =  x1[0];
	    ump_x1_1 += x1[1] * *left++;
	    ump_x1_1 += x1[2] * *left++;
	    ump_x1_1 += x1[3] * *left++;	      		  	 
	    
	    ump_x1_2 =  x1[0];
	    ump_x1_2 += x1[1] * *left++;
	    ump_x1_2 += x1[2] * *left++;
	    ump_x1_2 += x1[3] * *left++;	      		   	 
	    
	    ump_x1_3 =  x1[0];
	    ump_x1_3 += x1[1] * *left++;
	    ump_x1_3 += x1[2] * *left++;
	    ump_x1_3 += x1[3] * *left++;	      		   	  
	    
	    ump_x2_0 =  x2[0];
	    ump_x2_0 += x2[1] * *left++;
	    ump_x2_0 += x2[2] * *left++;
	    ump_x2_0 += x2[3] * *left++;		     	  	
	    
	    ump_x2_1 =  x2[0];
	    ump_x2_1 += x2[1] * *left++;
	    ump_x2_1 += x2[2] * *left++;
	    ump_x2_1 += x2[3] * *left++;		     	 	     
	    
	    ump_x2_2 =  x2[0];
	    ump_x2_2 += x2[1] * *left++;
	    ump_x2_2 += x2[2] * *left++;
	    ump_x2_2 += x2[3] * *left++;		     	   	      
	    
	    ump_x2_3 =  x2[0];
	    ump_x2_3 += x2[1] * *left++;
	    ump_x2_3 += x2[2] * *left++;
	    ump_x2_3 += x2[3] * *left++;
	    
	    x1px2 = ump_x1_0 * ump_x2_0;
	    x3[0] = x1px2 *  EV[0];
	    x3[1] = x1px2 *  EV[1];
	    x3[2] = x1px2 *  EV[2];
	    x3[3] = x1px2 *  EV[3];
	    
	    x1px2 = ump_x1_1 * ump_x2_1;
	    x3[0] += x1px2 *  EV[4];
	    x3[1] += x1px2 *  EV[5];
	    x3[2] += x1px2 *  EV[6];
	    x3[3] += x1px2 *  EV[7];
	    
	    x1px2 = ump_x1_2 * ump_x2_2;
	    x3[0] += x1px2 *  EV[8];
	    x3[1] += x1px2 *  EV[9];
	    x3[2] += x1px2 *  EV[10];
	    x3[3] += x1px2 *  EV[11];
	    
	    x1px2 = ump_x1_3 * ump_x2_3;
	    x3[0] += x1px2 *   EV[12];
	    x3[1] += x1px2 *   EV[13];
	    x3[2] += x1px2 *   EV[14];
	    x3[3] += x1px2 *   EV[15];
	    
	    ex3[i] = ex2[i];		  	      	       	       
	    
	    if (x3[0] < minlikelihood && x3[0] > minusminlikelihood &&
		x3[1] < minlikelihood && x3[1] > minusminlikelihood &&
		x3[2] < minlikelihood && x3[2] > minusminlikelihood &&
		x3[3] < minlikelihood && x3[3] > minusminlikelihood)	       
	      {	       
		x3[0]   *= twotothe256;
		x3[1]   *= twotothe256;
		x3[2]   *= twotothe256;		
		x3[3]   *= twotothe256;
		ex3[i]  += 1;			   
	      }	
	  }     	 
      }
      break;
    case INNER_INNER:
   
      for (i = lower; i < n; i++) 
	{	
	  model = modelptr[i];	
	  EV = &(extEV[model * 16]);
	  left = &left_start[model * 24 * numberOfCategories + cptr[i] * 24];
	  x1 = &x1_start[4 * i];
	  x2 = &x2_start[4 * i];
	  x3 = &x3_start[4 * i];
	    
	  ump_x1_0 =  x1[0];
	  ump_x1_0 += x1[1] * *left++;
	  ump_x1_0 += x1[2] * *left++;
	  ump_x1_0 += x1[3] * *left++;	      		   	    
	  
	  ump_x1_1 =  x1[0];
	  ump_x1_1 += x1[1] * *left++;
	  ump_x1_1 += x1[2] * *left++;
	  ump_x1_1 += x1[3] * *left++;	      		  	 
	  
	  ump_x1_2 =  x1[0];
	  ump_x1_2 += x1[1] * *left++;
	  ump_x1_2 += x1[2] * *left++;
	  ump_x1_2 += x1[3] * *left++;	      		   	 
	  
	  ump_x1_3 =  x1[0];
	  ump_x1_3 += x1[1] * *left++;
	  ump_x1_3 += x1[2] * *left++;
	  ump_x1_3 += x1[3] * *left++;	      		   	  
	  
	  ump_x2_0 =  x2[0];
	  ump_x2_0 += x2[1] * *left++;
	  ump_x2_0 += x2[2] * *left++;
	  ump_x2_0 += x2[3] * *left++;		     	  	
	  
	  ump_x2_1 =  x2[0];
	  ump_x2_1 += x2[1] * *left++;
	  ump_x2_1 += x2[2] * *left++;
	  ump_x2_1 += x2[3] * *left++;		     	 	     
	  
	  ump_x2_2 =  x2[0];
	  ump_x2_2 += x2[1] * *left++;
	  ump_x2_2 += x2[2] * *left++;
	  ump_x2_2 += x2[3] * *left++;		     	   	      
	  
	  ump_x2_3 =  x2[0];
	  ump_x2_3 += x2[1] * *left++;
	  ump_x2_3 += x2[2] * *left++;
	  ump_x2_3 += x2[3] * *left++;
	  
	  x1px2 = ump_x1_0 * ump_x2_0;
	  x3[0] = x1px2 *  EV[0];
	  x3[1] = x1px2 *  EV[1];
	  x3[2] = x1px2 *  EV[2];
	  x3[3] = x1px2 *  EV[3];
	  
	  x1px2 = ump_x1_1 * ump_x2_1;
	  x3[0] += x1px2 *  EV[4];
	  x3[1] += x1px2 *  EV[5];
	  x3[2] += x1px2 *  EV[6];
	  x3[3] += x1px2 *  EV[7];
	  
	  x1px2 = ump_x1_2 * ump_x2_2;
	  x3[0] += x1px2 *  EV[8];
	  x3[1] += x1px2 *  EV[9];
	  x3[2] += x1px2 *  EV[10];
	  x3[3] += x1px2 *  EV[11];
	  
	  x1px2 = ump_x1_3 * ump_x2_3;
	  x3[0] += x1px2 *   EV[12];
	  x3[1] += x1px2 *   EV[13];
	  x3[2] += x1px2 *   EV[14];
	  x3[3] += x1px2 *   EV[15];
	  
	  ex3[i] = ex1[i] + ex2[i];	
	  
	  if (x3[0] < minlikelihood && x3[0] > minusminlikelihood &&
	      x3[1] < minlikelihood && x3[1] > minusminlikelihood &&
	      x3[2] < minlikelihood && x3[2] > minusminlikelihood &&
	      x3[3] < minlikelihood && x3[3] > minusminlikelihood)  
	    {	     	      	     
	      x3[0]   *= twotothe256;
	      x3[1]   *= twotothe256;
	      x3[2]   *= twotothe256;		
	      x3[3]   *= twotothe256;
	      ex3[i]  += 1;
	    }	      	             
	}     	  
      break;
    default:
      assert(0);
    }

  free(left_start); 	 
}


static void newviewGTRGAMMA(traversalInfo *ti,
			       double *x1_start, double *x2_start, double *x3_start,
			       double *EIGN, double *EV, double *EI, double *gammaRates, double *tipVector,
			       int    *ex1, int *ex2, int *ex3, char *tipX1, char *tipX2,
			       int lower, const int n, double z1, double z2
			       )
{      
  double  
    *left, *right, *left_start, *right_start,
    *x1, *x2, *x3;      
  double  
    ump_x1_1, ump_x1_2, ump_x1_3, ump_x1_0, 
    ump_x2_0, ump_x2_1, ump_x2_2, ump_x2_3, x1px2,
    ki, d1c, d1g, d1t, d2c, d2g, d2t;
  int i;   
  
  left_start  = left  = (double *)malloc(48 * sizeof(double));
  right_start = right = (double *)malloc(48 * sizeof(double));
      
  for(i = 0; i < 4; i++)
    {       
      ki = gammaRates[i];	  	
      
      d1c = exp (EIGN[0] * ki * z1);
      d1g = exp (EIGN[1] * ki * z1);
      d1t = exp (EIGN[2] * ki * z1);	
      
      *left++ = d1c * EI[0];
      *left++ = d1g * EI[1];
      *left++ = d1t * EI[2];
      
      *left++ = d1c * EI[3];
      *left++ = d1g * EI[4];
      *left++ = d1t * EI[5];
      
      *left++ = d1c * EI[6];
      *left++ = d1g * EI[7];
      *left++ = d1t * EI[8];
      
      *left++ = d1c * EI[9];
      *left++ = d1g * EI[10];
      *left++ = d1t * EI[11];
      
      d2c = exp (EIGN[0] * ki * z2);
      d2g = exp (EIGN[1] * ki * z2);
      d2t = exp (EIGN[2] * ki * z2);		
      
      *right++ = d2c * EI[0];
      *right++ = d2g * EI[1];
      *right++ = d2t * EI[2];
      
      *right++ = d2c * EI[3];
      *right++ = d2g * EI[4];
      *right++ = d2t * EI[5];
      
      *right++ = d2c * EI[6];
      *right++ = d2g * EI[7];
      *right++ = d2t * EI[8];
      
      *right++ = d2c * EI[9];
      *right++ = d2g * EI[10];
      *right++ = d2t * EI[11];
    }                 
        	 
  

  switch(ti->tipCase)
    {
    case TIP_TIP:
      {      
	double *uX1, umpX1[256], *uX2, umpX2[256];
	uX1 = &umpX1[16];
	uX2 = &umpX2[16];	 

	for(i = 1; i < 16; i++)
	  {	    	     
	    x1 = &(tipVector[i * 4]);	    
	    
	    left = left_start;
	    right = right_start;
	    
	    ump_x1_0 =  x1[1] * *left++;
	    ump_x1_0 += x1[2] * *left++;
	    ump_x1_0 += x1[3] * *left++;	      			
	    ump_x1_0 += x1[0];
	    
	    *uX1++ = ump_x1_0;
	    
	    ump_x1_1 =  x1[1] * *left++;
	    ump_x1_1 += x1[2] * *left++;
	    ump_x1_1 += x1[3] * *left++;	      			
	    ump_x1_1 += x1[0];
	    
	    *uX1++ = ump_x1_1;
	    
	    ump_x1_2 =  x1[1] * *left++;
	    ump_x1_2 += x1[2] * *left++;
	    ump_x1_2 += x1[3] * *left++;		
	    ump_x1_2 += x1[0];
	    
	    *uX1++ = ump_x1_2;
	    
	    ump_x1_3 =  x1[1] * *left++;
	    ump_x1_3 += x1[2] * *left++;
	    ump_x1_3 += x1[3] * *left++;		 
	    ump_x1_3 += x1[0];
	    
	    *uX1++ = ump_x1_3;		
	    
	    ump_x2_0 =  x1[1] * *right++;
	    ump_x2_0 += x1[2] * *right++;
	    ump_x2_0 += x1[3] * *right++;	
	    ump_x2_0 += x1[0];
	    
	    *uX2++ = ump_x2_0;
	    
	    ump_x2_1 =  x1[1] * *right++;
	    ump_x2_1 += x1[2] * *right++;
	    ump_x2_1 += x1[3] * *right++;	
	    ump_x2_1 += x1[0];	 
	    
	    *uX2++ = ump_x2_1;
	    
	    ump_x2_2 =  x1[1] * *right++;
	    ump_x2_2 += x1[2] * *right++;
	    ump_x2_2 += x1[3] * *right++;	       
	    ump_x2_2 += x1[0];	  
	    
	    *uX2++ = ump_x2_2;
	    
	    ump_x2_3 =  x1[1] * *right++;
	    ump_x2_3 += x1[2] * *right++;
	    ump_x2_3 += x1[3] * *right++;		     	       
	    ump_x2_3 += x1[0];	    
	    
	    *uX2++ = ump_x2_3;
	    
	    ump_x1_0 =  x1[1] * *left++;
	    ump_x1_0 += x1[2] * *left++;
	    ump_x1_0 += x1[3] * *left++;	      			
	    ump_x1_0 += x1[0];
	    
	    *uX1++ = ump_x1_0;
	    
	    ump_x1_1 =  x1[1] * *left++;
	    ump_x1_1 += x1[2] * *left++;
	    ump_x1_1 += x1[3] * *left++;	      			
	    ump_x1_1 += x1[0];
	    
	    *uX1++ = ump_x1_1;
	    
	    ump_x1_2 =  x1[1] * *left++;
	    ump_x1_2 += x1[2] * *left++;
	    ump_x1_2 += x1[3] * *left++;		
	    ump_x1_2 += x1[0];
	    
	    *uX1++ = ump_x1_2;
	    
	    ump_x1_3 =  x1[1] * *left++;
	    ump_x1_3 += x1[2] * *left++;
	    ump_x1_3 += x1[3] * *left++;		 
	    ump_x1_3 += x1[0];
	    
	    *uX1++ = ump_x1_3;		
	    
	    ump_x2_0 =  x1[1] * *right++;
	    ump_x2_0 += x1[2] * *right++;
	    ump_x2_0 += x1[3] * *right++;	
	    ump_x2_0 += x1[0];
	    
	    *uX2++ = ump_x2_0;
	    
	    ump_x2_1 =  x1[1] * *right++;
	    ump_x2_1 += x1[2] * *right++;
	    ump_x2_1 += x1[3] * *right++;	
	    ump_x2_1 += x1[0];	 
	    
	    *uX2++ = ump_x2_1;
	    
	    ump_x2_2 =  x1[1] * *right++;
	    ump_x2_2 += x1[2] * *right++;
	    ump_x2_2 += x1[3] * *right++;	       
	    ump_x2_2 += x1[0];	  
	    
	    *uX2++ = ump_x2_2;
	    
	    ump_x2_3 =  x1[1] * *right++;
	    ump_x2_3 += x1[2] * *right++;
	    ump_x2_3 += x1[3] * *right++;		     	       
	    ump_x2_3 += x1[0];	    
	    
	    *uX2++ = ump_x2_3;
	    
	    ump_x1_0 =  x1[1] * *left++;
	    ump_x1_0 += x1[2] * *left++;
	    ump_x1_0 += x1[3] * *left++;	      			
	    ump_x1_0 += x1[0];
	    
	    *uX1++ = ump_x1_0;
	    
	    ump_x1_1 =  x1[1] * *left++;
	    ump_x1_1 += x1[2] * *left++;
	    ump_x1_1 += x1[3] * *left++;	      			
	    ump_x1_1 += x1[0];
	    
	    *uX1++ = ump_x1_1;
	    
	    ump_x1_2 =  x1[1] * *left++;
	    ump_x1_2 += x1[2] * *left++;
	    ump_x1_2 += x1[3] * *left++;		
	    ump_x1_2 += x1[0];
	    
	    *uX1++ = ump_x1_2;
	    
	    ump_x1_3 =  x1[1] * *left++;
	    ump_x1_3 += x1[2] * *left++;
	    ump_x1_3 += x1[3] * *left++;		 
	    ump_x1_3 += x1[0];
	    
	    *uX1++ = ump_x1_3;		
	    
	    ump_x2_0 =  x1[1] * *right++;
	    ump_x2_0 += x1[2] * *right++;
	    ump_x2_0 += x1[3] * *right++;	
	    ump_x2_0 += x1[0];
	    
	    *uX2++ = ump_x2_0;
	    
	    ump_x2_1 =  x1[1] * *right++;
	    ump_x2_1 += x1[2] * *right++;
	    ump_x2_1 += x1[3] * *right++;	
	    ump_x2_1 += x1[0];	 
	    
	    *uX2++ = ump_x2_1;
	    
	    ump_x2_2 =  x1[1] * *right++;
	    ump_x2_2 += x1[2] * *right++;
	    ump_x2_2 += x1[3] * *right++;	       
	    ump_x2_2 += x1[0];	  
	    
	    *uX2++ = ump_x2_2;
	    
	    ump_x2_3 =  x1[1] * *right++;
	    ump_x2_3 += x1[2] * *right++;
	    ump_x2_3 += x1[3] * *right++;		     	       
	    ump_x2_3 += x1[0];	    
	    
	    *uX2++ = ump_x2_3;
	    
	    ump_x1_0 =  x1[1] * *left++;
	    ump_x1_0 += x1[2] * *left++;
	    ump_x1_0 += x1[3] * *left++;	      			
	    ump_x1_0 += x1[0];
	    
	    *uX1++ = ump_x1_0;
	    
	    ump_x1_1 =  x1[1] * *left++;
	    ump_x1_1 += x1[2] * *left++;
	    ump_x1_1 += x1[3] * *left++;	      			
	    ump_x1_1 += x1[0];
	    
	    *uX1++ = ump_x1_1;
	    
	    ump_x1_2 =  x1[1] * *left++;
	    ump_x1_2 += x1[2] * *left++;
	    ump_x1_2 += x1[3] * *left++;		
	    ump_x1_2 += x1[0];
	    
	    *uX1++ = ump_x1_2;
	    
	    ump_x1_3 =  x1[1] * *left++;
	    ump_x1_3 += x1[2] * *left++;
	    ump_x1_3 += x1[3] * *left++;		 
	    ump_x1_3 += x1[0];
	    
	    *uX1++ = ump_x1_3;		
	    
	    ump_x2_0 =  x1[1] * *right++;
	    ump_x2_0 += x1[2] * *right++;
	    ump_x2_0 += x1[3] * *right++;	
	    ump_x2_0 += x1[0];
	    
	    *uX2++ = ump_x2_0;
	    
	    ump_x2_1 =  x1[1] * *right++;
	    ump_x2_1 += x1[2] * *right++;
	    ump_x2_1 += x1[3] * *right++;	
	    ump_x2_1 += x1[0];	 
	    
	    *uX2++ = ump_x2_1;
	    
	    ump_x2_2 =  x1[1] * *right++;
	    ump_x2_2 += x1[2] * *right++;
	    ump_x2_2 += x1[3] * *right++;	       
	    ump_x2_2 += x1[0];	  
	    
	    *uX2++ = ump_x2_2;
	    
	    ump_x2_3 =  x1[1] * *right++;
	    ump_x2_3 += x1[2] * *right++;
	    ump_x2_3 += x1[3] * *right++;		     	       
	    ump_x2_3 += x1[0];	    
	    
	    *uX2++ = ump_x2_3;	       
	  }			
#ifdef _USE_OMP   	 	
#pragma omp parallel for private(x3, uX1, uX2, x1px2)
	for (i = 0; i < n; i++)
#else
	for (i = lower; i < n; i++) 
#endif
	  {		     
	     uX1 = &umpX1[16 * tipX1[i]];
	     uX2 = &umpX2[16 * tipX2[i]];
	     x3 = &x3_start[16 * i];
	  			
	     x1px2 = *uX1++ * *uX2++;
	     x3[0] = x1px2 * EV[0];
	     x3[1] = x1px2 * EV[1];
	     x3[2] = x1px2 * EV[2];
	     x3[3] = x1px2 * EV[3];
	     
	     x1px2 = *uX1++ * *uX2++;
	     x3[0] += x1px2 * EV[4];
	     x3[1] += x1px2 * EV[5];
	     x3[2] += x1px2 * EV[6];
	     x3[3] += x1px2 * EV[7];
	    
	     x1px2 = *uX1++ * *uX2++;
	     x3[0] += x1px2 * EV[8];
	     x3[1] += x1px2 * EV[9];
	     x3[2] += x1px2 * EV[10];
	     x3[3] += x1px2 * EV[11];
	     
	     x1px2 = *uX1++ * *uX2++;
	     x3[0] += x1px2 * EV[12];
	     x3[1] += x1px2 * EV[13];
	     x3[2] += x1px2 * EV[14];
	     x3[3] += x1px2 * EV[15];

	     /* rate 1 */
	     	   	
	     x1px2 = *uX1++ * *uX2++;/*ump_x1_0 * ump_x2_0;*/
	     x3[4] = x1px2 * EV[0];
	     x3[5] = x1px2 * EV[1];
	     x3[6] = x1px2 * EV[2];
	     x3[7] = x1px2 * EV[3];
	     
	     x1px2 = *uX1++  * *uX2++;/*ump_x1_1 * ump_x2_1;*/
	     x3[4] += x1px2 * EV[4];
	     x3[5] += x1px2 * EV[5];
	     x3[6] += x1px2 * EV[6];
	     x3[7] += x1px2 * EV[7];
	     
	     x1px2 = *uX1++ * *uX2++;/*ump_x1_2 * ump_x2_2;*/
	     x3[4] += x1px2 * EV[8];
	     x3[5] += x1px2 * EV[9];
	     x3[6] += x1px2 * EV[10];
	     x3[7] += x1px2 * EV[11];
	    
	     x1px2 = *uX1++ * *uX2++;/*ump_x1_3 * ump_x2_3;*/
	     x3[4] += x1px2 * EV[12];
	     x3[5] += x1px2 * EV[13];
	     x3[6] += x1px2 * EV[14];
	     x3[7] += x1px2 * EV[15];

	     /* rate 2 */
	     	     	    
	     x1px2 = *uX1++ * *uX2++;/*ump_x1_0 * ump_x2_0;*/
	     x3[8] = x1px2  * EV[0];
	     x3[9] = x1px2  * EV[1];
	     x3[10] = x1px2 * EV[2];
	     x3[11] = x1px2 * EV[3];
	
	     x1px2 = *uX1++ * *uX2++;/*ump_x1_1 * ump_x2_1;*/
	     x3[8] += x1px2  * EV[4];
	     x3[9] += x1px2  * EV[5];
	     x3[10] += x1px2 * EV[6];
	     x3[11] += x1px2 * EV[7];
	    
	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_2 * ump_x2_2;*/
	     x3[8] += x1px2  * EV[8];
	     x3[9] += x1px2  * EV[9];
	     x3[10] += x1px2 * EV[10];
	     x3[11] += x1px2 * EV[11];
	     
	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_3 * ump_x2_3;*/
	     x3[8] += x1px2  * EV[12];
	     x3[9] += x1px2  * EV[13];
	     x3[10] += x1px2 * EV[14];
	     x3[11] += x1px2 * EV[15];

	     /* rate 3 */	    
	     
	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_0 * ump_x2_0;*/
	     x3[12] = x1px2 * EV[0];
	     x3[13] = x1px2 * EV[1];
	     x3[14] = x1px2 * EV[2];
	     x3[15] = x1px2 * EV[3];
	     
	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_1 * ump_x2_1;*/
	     x3[12] += x1px2 * EV[4];
	     x3[13] += x1px2 * EV[5];
	     x3[14] += x1px2 * EV[6];
	     x3[15] += x1px2 * EV[7];
	    
	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_2 * ump_x2_2;*/
	     x3[12] += x1px2 * EV[8];
	     x3[13] += x1px2 * EV[9];
	     x3[14] += x1px2 * EV[10];
	     x3[15] += x1px2 * EV[11];
	     
	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_3 * ump_x2_3;*/
	     x3[12] += x1px2 * EV[12];
	     x3[13] += x1px2 * EV[13];
	     x3[14] += x1px2 * EV[14];
	     x3[15] += x1px2 * EV[15];

	     /********************************************************************************/

	     ex3[i] = 0;		  		    	     	       	       
	   }                
      }
      break;      
    case TIP_INNER:
      {	 	 
	double *uX1, umpX1[256];	
	 
	uX1 = &umpX1[16];
	 	
	 for(i = 1; i < 16; i++)
	   {	    	     	     
	     x1 = &(tipVector[i * 4]);	     
	     left = left_start;
	     
	     ump_x1_0 =  x1[1] * *left++;
	     ump_x1_0 += x1[2] * *left++;
	     ump_x1_0 += x1[3] * *left++;	       
	     ump_x1_0 += x1[0];
	     
	     *uX1++ = ump_x1_0;
	     
	     ump_x1_1 =  x1[1] * *left++;
	     ump_x1_1 += x1[2] * *left++;
	     ump_x1_1 += x1[3] * *left++;		 
	     ump_x1_1 += x1[0];	 
	     
	     *uX1++ = ump_x1_1;
	     
	     ump_x1_2 =  x1[1] * *left++;
	     ump_x1_2 += x1[2] * *left++;
	     ump_x1_2 += x1[3] * *left++;		
	     ump_x1_2 += x1[0];	  
	     
	     *uX1++ = ump_x1_2;
	     
	     ump_x1_3 =  x1[1] * *left++;
	     ump_x1_3 += x1[2] * *left++;
	     ump_x1_3 += x1[3] * *left++;		     	      
	     ump_x1_3 += x1[0];	    
	     
	     *uX1++ = ump_x1_3;
	     
	     ump_x1_0 =  x1[1] * *left++;
	     ump_x1_0 += x1[2] * *left++;
	     ump_x1_0 += x1[3] * *left++;	       
	     ump_x1_0 += x1[0];
	     
	     *uX1++ = ump_x1_0;
	     
	     ump_x1_1 =  x1[1] * *left++;
	     ump_x1_1 += x1[2] * *left++;
	     ump_x1_1 += x1[3] * *left++;		 
	     ump_x1_1 += x1[0];	 
	     
	     *uX1++ = ump_x1_1;
	     
	     ump_x1_2 =  x1[1] * *left++;
	     ump_x1_2 += x1[2] * *left++;
	     ump_x1_2 += x1[3] * *left++;		
	     ump_x1_2 += x1[0];	  
	     
	     *uX1++ = ump_x1_2;
	     
	     ump_x1_3 =  x1[1] * *left++;
	     ump_x1_3 += x1[2] * *left++;
	     ump_x1_3 += x1[3] * *left++;		     	      
	     ump_x1_3 += x1[0];	    
	     
	     *uX1++ = ump_x1_3;
	     
	     ump_x1_0 =  x1[1] * *left++;
	     ump_x1_0 += x1[2] * *left++;
	     ump_x1_0 += x1[3] * *left++;	       
	     ump_x1_0 += x1[0];
	     
	     *uX1++ = ump_x1_0;
	     
	     ump_x1_1 =  x1[1] * *left++;
	     ump_x1_1 += x1[2] * *left++;
	     ump_x1_1 += x1[3] * *left++;		 
	     ump_x1_1 += x1[0];	 
	     
	     *uX1++ = ump_x1_1;
	     
	     ump_x1_2 =  x1[1] * *left++;
	     ump_x1_2 += x1[2] * *left++;
	     ump_x1_2 += x1[3] * *left++;		
	     ump_x1_2 += x1[0];	  
	     
	     *uX1++ = ump_x1_2;
	     
	     ump_x1_3 =  x1[1] * *left++;
	     ump_x1_3 += x1[2] * *left++;
	     ump_x1_3 += x1[3] * *left++;		     	      
	     ump_x1_3 += x1[0];	    
	     
	     *uX1++ = ump_x1_3;
	     
	     ump_x1_0 =  x1[1] * *left++;
	     ump_x1_0 += x1[2] * *left++;
	     ump_x1_0 += x1[3] * *left++;	       
	     ump_x1_0 += x1[0];
	     
	     *uX1++ = ump_x1_0;
	     
	     ump_x1_1 =  x1[1] * *left++;
	     ump_x1_1 += x1[2] * *left++;
	     ump_x1_1 += x1[3] * *left++;		 
	     ump_x1_1 += x1[0];	 
	     
	     *uX1++ = ump_x1_1;
	     
	     ump_x1_2 =  x1[1] * *left++;
	     ump_x1_2 += x1[2] * *left++;
	     ump_x1_2 += x1[3] * *left++;		
	     ump_x1_2 += x1[0];	  
	     
	     *uX1++ = ump_x1_2;
	     
	     ump_x1_3 =  x1[1] * *left++;
	     ump_x1_3 += x1[2] * *left++;
	     ump_x1_3 += x1[3] * *left++;		     	      
	     ump_x1_3 += x1[0];	    
	     
	     *uX1++ = ump_x1_3;	    
	   }	            	    	       	
#ifdef _USE_OMP	 
#pragma omp parallel for private(x1px2, x2, x3, uX1, right, ump_x2_0, ump_x2_1, ump_x2_2, ump_x2_3)
	 for (i = 0; i < n; i++)
#else
	 for (i = lower; i < n; i++) 
#endif
	   {	
	     right = right_start;
	     uX1 = &umpX1[16 * tipX1[i]];			    
	     x2 = &x2_start[16 * i];
	     x3 = &x3_start[16 * i];
	     	  					 	 	    	    	     	     
	     ump_x2_0 = x2[1] * *right++;
	     ump_x2_0 += x2[2] * *right++;
	     ump_x2_0 += x2[3] * *right++;		     	   
	     ump_x2_0 += x2[0];
	     
	     ump_x2_1 = x2[1] * *right++;
	     ump_x2_1 += x2[2] * *right++;
	     ump_x2_1 +=  x2[3] * *right++;		     	    
	     ump_x2_1 += x2[0];	 
	     
	     ump_x2_2 =  x2[1] * *right++;
	     ump_x2_2 += x2[2] * *right++;
	     ump_x2_2 += x2[3] * *right++;		     	     
	     ump_x2_2 += x2[0];	  
	     
	     ump_x2_3 =  x2[1] * *right++;
	     ump_x2_3 += x2[2] * *right++;
	     ump_x2_3 += x2[3] * *right++;		     	     
	     ump_x2_3 += x2[0];	    	  		   	  	    
	     
	     x1px2 = *uX1++ * ump_x2_0;
	     x3[0] = x1px2 * EV[0];
	     x3[1] = x1px2 * EV[1];
	     x3[2] = x1px2 * EV[2];
	     x3[3] = x1px2 * EV[3];
	     
	     x1px2 = *uX1++ * ump_x2_1;
	     x3[0] += x1px2 * EV[4];
	     x3[1] += x1px2 * EV[5];
	     x3[2] += x1px2 * EV[6];
	     x3[3] += x1px2 * EV[7];
	     
	     x1px2 = *uX1++ * ump_x2_2;
	     x3[0] += x1px2 *  EV[8];
	     x3[1] += x1px2 *  EV[9];
	     x3[2] += x1px2 *  EV[10];
	     x3[3] += x1px2 *  EV[11];
	     
	     x1px2 = *uX1++ * ump_x2_3;
	     x3[0] += x1px2 *  EV[12];
	     x3[1] += x1px2 *  EV[13];
	     x3[2] += x1px2 *  EV[14];
	     x3[3] += x1px2 *  EV[15];
	     	      	    	    	     	     
	     ump_x2_0 = x2[5] * *right++;
	     ump_x2_0 += x2[6] * *right++;
	     ump_x2_0 += x2[7] * *right++;		     	  
	     ump_x2_0 += x2[4];
	     
	     ump_x2_1 = x2[5] * *right++;
	     ump_x2_1 += x2[6] * *right++;
	     ump_x2_1 +=  x2[7] * *right++;		     	   
	     ump_x2_1 += x2[4];	 
	     
	     ump_x2_2 = x2[5] * *right++;
	     ump_x2_2 += x2[6] * *right++;
	     ump_x2_2 +=  x2[7] * *right++;		     	  
	     ump_x2_2 += x2[4];	  
	     
	     ump_x2_3 = x2[5] * *right++;
	     ump_x2_3 += x2[6] * *right++;
	     ump_x2_3 += x2[7] * *right++;		     	 
	     ump_x2_3 += x2[4];	    	  		   	  	    
	     
	     x1px2 = *uX1++ * ump_x2_0;
	     x3[4] = x1px2 * EV[0];
	     x3[5] = x1px2 * EV[1];
	     x3[6] = x1px2 * EV[2];
	     x3[7] = x1px2 * EV[3];
	     
	     x1px2 = *uX1++ * ump_x2_1;
	     x3[4] += x1px2 * EV[4];
	     x3[5] += x1px2 * EV[5];
	     x3[6] += x1px2 * EV[6];
	     x3[7] += x1px2 * EV[7];
	     
	     x1px2 = *uX1++ * ump_x2_2;
	     x3[4] += x1px2 *  EV[8];
	     x3[5] += x1px2 *  EV[9];
	     x3[6] += x1px2 *  EV[10];
	     x3[7] += x1px2 *  EV[11];
	     
	     x1px2 = *uX1++ * ump_x2_3;
	     x3[4] += x1px2 *   EV[12];
	     x3[5] += x1px2 *   EV[13];
	     x3[6] += x1px2 *   EV[14];
	     x3[7] += x1px2 *   EV[15];	     	 
	     
	     ump_x2_0 = x2[9] * *right++;
	     ump_x2_0 += x2[10] * *right++;
	     ump_x2_0 += x2[11] * *right++;		     	   
	     ump_x2_0 += x2[8];
	     
	     ump_x2_1 = x2[9] * *right++;
	     ump_x2_1 += x2[10] * *right++;
	     ump_x2_1 +=  x2[11] * *right++;		     	    
	     ump_x2_1 += x2[8];	 
	     
	     ump_x2_2 = x2[9] * *right++;
	     ump_x2_2 += x2[10] * *right++;
	     ump_x2_2 +=  x2[11] * *right++;		     	     
	     ump_x2_2 += x2[8];	  
	     
	     ump_x2_3 = x2[9] * *right++;
	     ump_x2_3 += x2[10] * *right++;
	     ump_x2_3 += x2[11] * *right++;		     	   
	     ump_x2_3 += x2[8];	    	  		   	  	    
	     
	     x1px2 = *uX1++ * ump_x2_0;
	     x3[8] = x1px2 *  EV[0];
	     x3[9] = x1px2 *  EV[1];
	     x3[10] = x1px2 * EV[2];
	     x3[11] = x1px2 *  EV[3];
	     
	     x1px2 = *uX1++ * ump_x2_1;
	     x3[8] += x1px2  *  EV[4];
	     x3[9] += x1px2 *   EV[5];
	     x3[10] += x1px2 *  EV[6];
	     x3[11] += x1px2 *   EV[7];
	     
	     x1px2 = *uX1++ * ump_x2_2;
	     x3[8] += x1px2 *  EV[8];
	     x3[9] += x1px2*   EV[9];
	     x3[10] += x1px2 *  EV[10];
	     x3[11] += x1px2 *  EV[11];
	     
	     x1px2 = *uX1++ * ump_x2_3;
	     x3[8] += x1px2 *   EV[12];
	     x3[9] += x1px2 *   EV[13];
	     x3[10] += x1px2 *   EV[14];
	     x3[11] += x1px2 *   EV[15];
	     	     	 		 	    	    	     	     
	     ump_x2_0 = x2[13] * *right++;
	     ump_x2_0 += x2[14] * *right++;
	     ump_x2_0 += x2[15] * *right++;		     	   
	     ump_x2_0 += x2[12];
	     
	     ump_x2_1 = x2[13] * *right++;
	     ump_x2_1 += x2[14] * *right++;
	     ump_x2_1 +=  x2[15] * *right++;		     	    
	     ump_x2_1 += x2[12];	 
	     
	     ump_x2_2 = x2[13] * *right++;
	     ump_x2_2 += x2[14] * *right++;
	     ump_x2_2 +=  x2[15] * *right++;		     	     
	     ump_x2_2 += x2[12];	  
	     
	     ump_x2_3 = x2[13] * *right++;
	     ump_x2_3 += x2[14] * *right++;
	     ump_x2_3 += x2[15] * *right++;		     	    
	     ump_x2_3 += x2[12];	    	  		   	  	    
	     
	     x1px2 = *uX1++ * ump_x2_0;
	     x3[12] = x1px2 *  EV[0];
	     x3[13] = x1px2 *  EV[1];
	     x3[14] = x1px2 * EV[2];
	     x3[15] = x1px2 *  EV[3];
	     
	     x1px2 = *uX1++ * ump_x2_1;
	     x3[12] += x1px2  *  EV[4];
	     x3[13] += x1px2 *   EV[5];
	     x3[14] += x1px2 *  EV[6];
	     x3[15] += x1px2 *   EV[7];
	     
	     x1px2 = *uX1++ * ump_x2_2;
	     x3[12] += x1px2 *  EV[8];
	     x3[13] += x1px2*   EV[9];
	     x3[14] += x1px2 *  EV[10];
	     x3[15] += x1px2 *  EV[11];
	     
	     x1px2 = *uX1++ * ump_x2_3;
	     x3[12] += x1px2 *   EV[12];
	     x3[13] += x1px2 *   EV[13];
	     x3[14] += x1px2 *   EV[14];
	     x3[15] += x1px2 *   EV[15];
	     
	     ex3[i] = ex2[i];		  
	     
	     if (ABS(x3[0]) < minlikelihood && ABS(x3[2]) < minlikelihood && ABS(x3[1]) < minlikelihood && ABS(x3[3]) < minlikelihood &&
		 ABS(x3[4]) < minlikelihood && ABS(x3[6]) < minlikelihood && ABS(x3[5]) < minlikelihood && ABS(x3[7]) < minlikelihood &&
		 ABS(x3[8]) < minlikelihood && ABS(x3[10]) < minlikelihood && ABS(x3[9]) < minlikelihood && ABS(x3[11]) < minlikelihood &&
		 ABS(x3[12]) < minlikelihood && ABS(x3[14]) < minlikelihood && ABS(x3[13]) < minlikelihood && ABS(x3[15]) < minlikelihood)	         
	       {	     	    	 		    		 
		 x3[0]   *= twotothe256;
		 x3[1]   *= twotothe256;
		 x3[2]   *= twotothe256;
		 x3[3]   *= twotothe256;
		 
		 x3[4]   *= twotothe256;
		 x3[5]   *= twotothe256;
		 x3[6]   *= twotothe256;
		 x3[7]   *= twotothe256;
		 
		 x3[8]   *= twotothe256;
		 x3[9]   *= twotothe256;
		 x3[10]   *= twotothe256;
		 x3[11]   *= twotothe256;
		 
		 x3[12]   *= twotothe256;
		 x3[13]   *= twotothe256;
		 x3[14]   *= twotothe256;
		 x3[15]   *= twotothe256;
		 
		 ex3[i] += 1;
	       }	 	      	   
	   }
      }
      break;
    case INNER_INNER:

#ifdef _USE_OMP	
#pragma omp parallel for private(x1, x2, x3, left, right, ump_x2_0, ump_x2_1, ump_x2_2, ump_x2_3,ump_x1_0, ump_x1_1, ump_x1_2, ump_x1_3, x1px2)    
      for (i = 0; i < n; i++)
#else
     for (i = lower; i < n; i++) 
#endif
       {	
	 left = left_start;
	 right = right_start;	      
	 x1 = &x1_start[16 * i];
	 x2 = &x2_start[16 * i];
	 x3 = &x3_start[16 * i];
	 		 
	 /* Rate cat 0 */       
	 ump_x1_0 = x1[1] * *left++;
	 ump_x1_0 += x1[2] * *left++;
	 ump_x1_0 += x1[3]* *left++;	      	      
	 ump_x1_0 += x1[0];
	 
	 ump_x1_1 =  x1[1] * *left++;
	 ump_x1_1 += x1[2] * *left++;
	 ump_x1_1 += x1[3]* *left++;	      	       
	 ump_x1_1 += x1[0];
	 
	 ump_x1_2 = x1[1] * *left++;
	 ump_x1_2 += x1[2] * *left++;
	 ump_x1_2 += x1[3] * *left++;	      		
	 ump_x1_2 += x1[0];
	 
	 ump_x1_3 =  x1[1] * *left++;
	 ump_x1_3 += x1[2] * *left++;
	 ump_x1_3 += x1[3] * *left++;	      	      
	 ump_x1_3 += x1[0];
	 
	 ump_x2_0 = x2[1] * *right++;
	 ump_x2_0 += x2[2] * *right++;
	 ump_x2_0 += x2[3] * *right++;		     	
	 ump_x2_0 += x2[0];
	 
	 ump_x2_1 = x2[1] * *right++;
	 ump_x2_1 += x2[2] * *right++;
	 ump_x2_1 +=  x2[3] * *right++;		     	
	 ump_x2_1 += x2[0];	 
	 
	 ump_x2_2 = x2[1] * *right++;
	 ump_x2_2 += x2[2] * *right++;
	 ump_x2_2 +=  x2[3] * *right++;		     	
	 ump_x2_2 += x2[0];	  
	 
	 ump_x2_3 = x2[1] * *right++;
	 ump_x2_3 += x2[2] * *right++;
	 ump_x2_3 += x2[3] * *right++;		     	
	 ump_x2_3 += x2[0];	    	  		   	  	    
	 
	 x1px2 = ump_x1_0 * ump_x2_0;
	 x3[0] = x1px2 * EV[0];
	 x3[1] = x1px2 * EV[1];
	 x3[2] = x1px2 * EV[2];
	 x3[3] = x1px2 * EV[3];
	 
	 x1px2 = ump_x1_1 * ump_x2_1;
	 x3[0] += x1px2 * EV[4];
	 x3[1] += x1px2 * EV[5];
	 x3[2] += x1px2 * EV[6];
	 x3[3] += x1px2 * EV[7];
	 
	 x1px2 = ump_x1_2 * ump_x2_2;
	 x3[0] += x1px2 * EV[8];
	 x3[1] += x1px2 * EV[9];
	 x3[2] += x1px2 * EV[10];
	 x3[3] += x1px2 * EV[11];
	 
	 x1px2 = ump_x1_3 * ump_x2_3;
	 x3[0] += x1px2 * EV[12];
	 x3[1] += x1px2 * EV[13];
	 x3[2] += x1px2 * EV[14];
	 x3[3] += x1px2 * EV[15];
	 
	 /* rate 1 */
	        
	 ump_x1_0  = x1[5] * *left++;
	 ump_x1_0 += x1[6] * *left++;
	 ump_x1_0 += x1[7] * *left++;	      	     
	 ump_x1_0 += x1[4];
	 
	 ump_x1_1 =  x1[5] * *left++;
	 ump_x1_1 += x1[6] * *left++;
	 ump_x1_1 += x1[7] * *left++;	      	      
	 ump_x1_1 += x1[4];
	 
	 ump_x1_2  = x1[5] * *left++;
	 ump_x1_2 += x1[6] * *left++;
	 ump_x1_2 += x1[7] * *left++;	      	      
	 ump_x1_2 += x1[4];
	 
	 ump_x1_3 =  x1[5] * *left++;
	 ump_x1_3 += x1[6] * *left++;
	 ump_x1_3 += x1[7] * *left++;	      	
	 ump_x1_3 += x1[4];
	 
	 ump_x2_0  = x2[5] * *right++;
	 ump_x2_0 += x2[6] * *right++;
	 ump_x2_0 += x2[7] * *right++;		           
	 ump_x2_0 += x2[4];
	 
	 ump_x2_1  = x2[5] * *right++;
	 ump_x2_1 += x2[6] * *right++;
	 ump_x2_1 += x2[7] * *right++;		           
	 ump_x2_1 += x2[4];	 
	 
	 ump_x2_2  = x2[5] * *right++;
	 ump_x2_2 += x2[6] * *right++;
	 ump_x2_2 += x2[7] * *right++;		           
	 ump_x2_2 += x2[4];	  
	 
	 ump_x2_3  = x2[5] * *right++;
	 ump_x2_3 += x2[6] * *right++;
	 ump_x2_3 += x2[7] * *right++;		            
	 ump_x2_3 += x2[4];	    	  		   	  	    
	 
	 x1px2 = ump_x1_0 * ump_x2_0;
	 x3[4] = x1px2 * EV[0];
	 x3[5] = x1px2 * EV[1];
	 x3[6] = x1px2 * EV[2];
	 x3[7] = x1px2 * EV[3];
	 
	 x1px2 = ump_x1_1 * ump_x2_1;
	 x3[4] += x1px2 * EV[4];
	 x3[5] += x1px2 * EV[5];
	 x3[6] += x1px2 * EV[6];
	 x3[7] += x1px2 * EV[7];
	 
	 x1px2 = ump_x1_2 * ump_x2_2;
	 x3[4] += x1px2 *  EV[8];
	 x3[5] += x1px2 *  EV[9];
	 x3[6] += x1px2 *  EV[10];
	 x3[7] += x1px2 *  EV[11];
	 
	 x1px2 = ump_x1_3 * ump_x2_3;
	 x3[4] += x1px2 *  EV[12];
	 x3[5] += x1px2 *  EV[13];
	 x3[6] += x1px2 *  EV[14];
	 x3[7] += x1px2 *  EV[15];
	 
	 /* rate 2 */
	       
	 ump_x1_0  = x1[9]  * *left++;
	 ump_x1_0 += x1[10] * *left++;
	 ump_x1_0 += x1[11] * *left++;	      	     
	 ump_x1_0 += x1[8];
	 
	 ump_x1_1  =  x1[9] * *left++;
	 ump_x1_1 += x1[10] * *left++;
	 ump_x1_1 += x1[11] * *left++;	      	      
	 ump_x1_1 += x1[8];
	 
	 ump_x1_2  = x1[9]  * *left++;
	 ump_x1_2 += x1[10] * *left++;
	 ump_x1_2 += x1[11] * *left++;	      	     
	 ump_x1_2 += x1[8];
	 
	 ump_x1_3  = x1[9]  * *left++;
	 ump_x1_3 += x1[10] * *left++;
	 ump_x1_3 += x1[11] * *left++;	      	     
	 ump_x1_3 += x1[8];
	 	 
	 ump_x2_0  = x2[9]  * *right++;
	 ump_x2_0 += x2[10] * *right++;
	 ump_x2_0 += x2[11] * *right++;		     
	 ump_x2_0 += x2[8];
	 
	 ump_x2_1  = x2[9]  * *right++;
	 ump_x2_1 += x2[10] * *right++;
	 ump_x2_1 += x2[11] * *right++;		         
	 ump_x2_1 += x2[8];	 
	 
	 ump_x2_2  = x2[9]  * *right++;
	 ump_x2_2 += x2[10] * *right++;
	 ump_x2_2 += x2[11] * *right++;		            
	 ump_x2_2 += x2[8];	  
	 
	 ump_x2_3 = x2[9] * *right++;
	 ump_x2_3 += x2[10] * *right++;
	 ump_x2_3 += x2[11] * *right++;		          
	 ump_x2_3 += x2[8];	    	  		   	  	    
	 
	 x1px2 = ump_x1_0 * ump_x2_0;
	 x3[8]  = x1px2 * EV[0];
	 x3[9]  = x1px2 * EV[1];
	 x3[10] = x1px2 * EV[2];
	 x3[11] = x1px2 * EV[3];
	 
	 x1px2 = ump_x1_1 * ump_x2_1;
	 x3[8]  += x1px2 * EV[4];
	 x3[9]  += x1px2 * EV[5];
	 x3[10] += x1px2 * EV[6];
	 x3[11] += x1px2 * EV[7];
	 
	 x1px2 = ump_x1_2 * ump_x2_2;
	 x3[8]  += x1px2 * EV[8];
	 x3[9]  += x1px2 * EV[9];
	 x3[10] += x1px2 * EV[10];
	 x3[11] += x1px2 * EV[11];
	 
	 x1px2 = ump_x1_3 * ump_x2_3;
	 x3[8]  += x1px2 * EV[12];
	 x3[9]  += x1px2 * EV[13];
	 x3[10] += x1px2 * EV[14];
	 x3[11] += x1px2 * EV[15];
	 
	 /* rate 3 */
      
	 ump_x1_0  = x1[13] * *left++;
	 ump_x1_0 += x1[14] * *left++;
	 ump_x1_0 += x1[15] * *left++;	      	     
	 ump_x1_0 += x1[12];
	 
	 ump_x1_1 =  x1[13] * *left++;
	 ump_x1_1 += x1[14] * *left++;
	 ump_x1_1 += x1[15] * *left++;	      	 
	 ump_x1_1 += x1[12];
	 
	 ump_x1_2  = x1[13] * *left++;
	 ump_x1_2 += x1[14] * *left++;
	 ump_x1_2 += x1[15] * *left++;	      	      
	 ump_x1_2 += x1[12];
	 
	 ump_x1_3  = x1[13] * *left++;
	 ump_x1_3 += x1[14] * *left++;
	 ump_x1_3 += x1[15] * *left++;	      	     
	 ump_x1_3 += x1[12];
	 
	 ump_x2_0  = x2[13] * *right++;
	 ump_x2_0 += x2[14] * *right++;
	 ump_x2_0 += x2[15] * *right++;		         
	 ump_x2_0 += x2[12];
	 
	 ump_x2_1  = x2[13] * *right++;
	 ump_x2_1 += x2[14] * *right++;
	 ump_x2_1 += x2[15] * *right++;		           
	 ump_x2_1 += x2[12];	 
	 
	 ump_x2_2  = x2[13] * *right++;
	 ump_x2_2 += x2[14] * *right++;
	 ump_x2_2 += x2[15] * *right++;		         
	 ump_x2_2 += x2[12];	  
	 
	 ump_x2_3  = x2[13] * *right++;
	 ump_x2_3 += x2[14] * *right++;
	 ump_x2_3 += x2[15] * *right++;		            
	 ump_x2_3 += x2[12];	    	  		   	  	    
	 
	 x1px2 = ump_x1_0 * ump_x2_0;
	 x3[12] = x1px2 * EV[0];
	 x3[13] = x1px2 * EV[1];
	 x3[14] = x1px2 * EV[2];
	 x3[15] = x1px2 * EV[3];
	 
	 x1px2 = ump_x1_1 * ump_x2_1;
	 x3[12] += x1px2 * EV[4];
	 x3[13] += x1px2 * EV[5];
	 x3[14] += x1px2 * EV[6];
	 x3[15] += x1px2 * EV[7];
	 
	 x1px2 = ump_x1_2 * ump_x2_2;
	 x3[12] += x1px2 * EV[8];
	 x3[13] += x1px2 * EV[9];
	 x3[14] += x1px2 * EV[10];
	 x3[15] += x1px2 * EV[11];
	 
	 x1px2 = ump_x1_3 * ump_x2_3;
	 x3[12] += x1px2 * EV[12];
	 x3[13] += x1px2 * EV[13];
	 x3[14] += x1px2 * EV[14];
	 x3[15] += x1px2 * EV[15];
	 
	 /********************************************************************************/

	 ex3[i] = ex1[i] + ex2[i];		  
	 
	 if (ABS(x3[0]) < minlikelihood && ABS(x3[2]) < minlikelihood && ABS(x3[1]) < minlikelihood && ABS(x3[3]) < minlikelihood &&
	     ABS(x3[4]) < minlikelihood && ABS(x3[6]) < minlikelihood && ABS(x3[5]) < minlikelihood && ABS(x3[7]) < minlikelihood &&
	     ABS(x3[8]) < minlikelihood && ABS(x3[10]) < minlikelihood && ABS(x3[9]) < minlikelihood && ABS(x3[11]) < minlikelihood &&
	     ABS(x3[12]) < minlikelihood && ABS(x3[14]) < minlikelihood && ABS(x3[13]) < minlikelihood && ABS(x3[15]) < minlikelihood)	
	   {			    	      
	     x3[0]   *= twotothe256;
	     x3[1]   *= twotothe256;
	     x3[2]   *= twotothe256;
	     x3[3]   *= twotothe256;
	     
	     x3[4]   *= twotothe256;
	     x3[5]   *= twotothe256;
	     x3[6]   *= twotothe256;
	     x3[7]   *= twotothe256;
	     
	     x3[8]   *= twotothe256;
	     x3[9]   *= twotothe256;
	     x3[10]   *= twotothe256;
	     x3[11]   *= twotothe256;
	     
	     x3[12]   *= twotothe256;
	     x3[13]   *= twotothe256;
	     x3[14]   *= twotothe256;
	     x3[15]   *= twotothe256;
	     
	     ex3[i] += 1;
	   }		 	          
       }
     break;
    default:
      assert(0);
    }
  
  free(left_start); 
  free(right_start);         
}

static void newviewGTRGAMMAMULT(traversalInfo *ti,
				   double *x1_start, double *x2_start, double *x3_start,
				   double *extEIGN, double *extEV, double *extEI, double *gammaRates, double *tipVector,
				   int    *ex1, int *ex2, int *ex3, char *tipX1, char *tipX2, int *modelptr,
				int lower, int n, int numberOfModels, int multiBranch
				   )
{ 
  double  
    *left, *right, *left_start, *right_start,
    *x1, *x2, *x3, *EIGN, *EI, *EV;    
  double  
    ump_x1_1, ump_x1_2, ump_x1_3, ump_x1_0, 
    ump_x2_0, ump_x2_1, ump_x2_2, ump_x2_3, x1px2,
    z1 = 0.0, 
    z2 = 0.0, 
    ki,
    d1c, d1g, d1t, d2c, d2g, d2t;                
  int 
    model, i;

  if(!multiBranch)
    {
      z1  = ti->qz[0];     
      z2  = ti->rz[0];     
    }
  
  left_start =  left = (double *) malloc(48 * numberOfModels * sizeof(double));
  right_start = right = (double *)malloc(48 * numberOfModels * sizeof(double));  

  for(model = 0; model < numberOfModels; model++)
    {
      if(multiBranch)
	{
	  z1  = ti->qz[model];	
	  z2  = ti->rz[model];	  
	}
      
      EIGN = &(extEIGN[3 * model]);
      
      left  = &left_start[model * 48];
      right = &right_start[model * 48];      
      
      EI = &(extEI[model * 12]);

      for(i = 0; i < 4; i++)
	{       
	  ki = gammaRates[model * 4 + i];	  	 
	  
	  d1c = exp (EIGN[0] * ki * z1);
	  d1g = exp (EIGN[1] * ki * z1);
	  d1t = exp (EIGN[2] * ki * z1);	
	  
	  *left++ = d1c * EI[0];
	  *left++ = d1g * EI[1];
	  *left++ = d1t * EI[2];
	  
	  *left++ = d1c * EI[3];
	  *left++ = d1g * EI[4];
	  *left++ = d1t * EI[5];
	  
	  *left++ = d1c * EI[6];
	  *left++ = d1g * EI[7];
	  *left++ = d1t * EI[8];
	  
	  *left++ = d1c * EI[9];
	  *left++ = d1g * EI[10];
	  *left++ = d1t * EI[11];
	  	
	  d2c = exp (EIGN[0] * ki * z2);
	  d2g = exp (EIGN[1] * ki * z2);
	  d2t = exp (EIGN[2] * ki * z2);		
	  
	  *right++ = d2c * EI[0];
	  *right++ = d2g * EI[1];
	  *right++ = d2t * EI[2];
	  
	  *right++ = d2c * EI[3];
	  *right++ = d2g * EI[4];
	  *right++ = d2t * EI[5];
	  
	  *right++ = d2c * EI[6];
	  *right++ = d2g * EI[7];
	  *right++ = d2t * EI[8];
	  
	  *right++ = d2c * EI[9];
	  *right++ = d2g * EI[10];
	  *right++ = d2t * EI[11];
	}                 
    }
  
  
  switch(ti->tipCase)
    {
    case TIP_TIP:
      {	     
	double *uX1, *umpX1, *uX2, *umpX2;
      
	umpX1 = (double *)malloc(256 * numberOfModels * sizeof(double));
	umpX2 = (double *)malloc(256 * numberOfModels * sizeof(double));
           
	for(model = 0; model < numberOfModels; model++)
	{
	  uX1 = &umpX1[256 * model + 16];
	  uX2 = &umpX2[256 * model + 16];
	  
	  for(i = 1; i < 16; i++)
	    {	    	     
	      x1 = &(tipVector[model * 64 + i * 4]);
	      x2 = &(tipVector[model * 64 + i * 4]);
	      
	      left  = &left_start[model * 48];
	      right = &right_start[model * 48];
	      
	      ump_x1_0 =  x1[1] * *left++;
	      ump_x1_0 += x1[2] * *left++;
	      ump_x1_0 += x1[3] * *left++;	      			
	      ump_x1_0 += x1[0];
	      
	      *uX1++ = ump_x1_0;
	      
	      ump_x1_1 =  x1[1] * *left++;
	      ump_x1_1 += x1[2] * *left++;
	      ump_x1_1 += x1[3] * *left++;	      			
	      ump_x1_1 += x1[0];
	      
		 *uX1++ = ump_x1_1;
		 
		 ump_x1_2 =  x1[1] * *left++;
		 ump_x1_2 += x1[2] * *left++;
		 ump_x1_2 += x1[3] * *left++;		
		 ump_x1_2 += x1[0];
		 
		 *uX1++ = ump_x1_2;
		 
		 ump_x1_3 =  x1[1] * *left++;
		 ump_x1_3 += x1[2] * *left++;
		 ump_x1_3 += x1[3] * *left++;		 
		 ump_x1_3 += x1[0];
		 
		 *uX1++ = ump_x1_3;		
		 
		 ump_x2_0 =  x2[1] * *right++;
		 ump_x2_0 += x2[2] * *right++;
		 ump_x2_0 += x2[3] * *right++;	
		 ump_x2_0 += x2[0];
		 
		 *uX2++ = ump_x2_0;
		 
		 ump_x2_1 =  x2[1] * *right++;
		 ump_x2_1 += x2[2] * *right++;
		 ump_x2_1 += x2[3] * *right++;	
		 ump_x2_1 += x2[0];	 
		 
		 *uX2++ = ump_x2_1;
		 
		 ump_x2_2 =  x2[1] * *right++;
		 ump_x2_2 += x2[2] * *right++;
		 ump_x2_2 += x2[3] * *right++;	       
		 ump_x2_2 += x2[0];	  
		 
		 *uX2++ = ump_x2_2;
		 
		 ump_x2_3 =  x2[1] * *right++;
		 ump_x2_3 += x2[2] * *right++;
		 ump_x2_3 += x2[3] * *right++;		     	       
		 ump_x2_3 += x2[0];	    
		 
		 *uX2++ = ump_x2_3;
		 
		 ump_x1_0 =  x1[1] * *left++;
		 ump_x1_0 += x1[2] * *left++;
		 ump_x1_0 += x1[3] * *left++;	      			
		 ump_x1_0 += x1[0];
		 
		 *uX1++ = ump_x1_0;
		 
		 ump_x1_1 =  x1[1] * *left++;
		 ump_x1_1 += x1[2] * *left++;
		 ump_x1_1 += x1[3] * *left++;	      			
		 ump_x1_1 += x1[0];
		 
		 *uX1++ = ump_x1_1;
		 
		 ump_x1_2 =  x1[1] * *left++;
		 ump_x1_2 += x1[2] * *left++;
		 ump_x1_2 += x1[3] * *left++;		
		 ump_x1_2 += x1[0];
		 
		 *uX1++ = ump_x1_2;
		 
		 ump_x1_3 =  x1[1] * *left++;
		 ump_x1_3 += x1[2] * *left++;
		 ump_x1_3 += x1[3] * *left++;		 
		 ump_x1_3 += x1[0];
		 
		 *uX1++ = ump_x1_3;		
		 
		 ump_x2_0 =  x2[1] * *right++;
		 ump_x2_0 += x2[2] * *right++;
		 ump_x2_0 += x2[3] * *right++;	
		 ump_x2_0 += x2[0];
		 
		 *uX2++ = ump_x2_0;
		 
		 ump_x2_1 =  x2[1] * *right++;
		 ump_x2_1 += x2[2] * *right++;
		 ump_x2_1 += x2[3] * *right++;	
		 ump_x2_1 += x2[0];	 
		 
		 *uX2++ = ump_x2_1;
		 
		 ump_x2_2 =  x2[1] * *right++;
		 ump_x2_2 += x2[2] * *right++;
		 ump_x2_2 += x2[3] * *right++;	       
		 ump_x2_2 += x2[0];	  
		 
		 *uX2++ = ump_x2_2;
		 
		 ump_x2_3 =  x2[1] * *right++;
		 ump_x2_3 += x2[2] * *right++;
		 ump_x2_3 += x2[3] * *right++;		     	       
		 ump_x2_3 += x2[0];	    
		 
		 *uX2++ = ump_x2_3;
		 
		 ump_x1_0 =  x1[1] * *left++;
		 ump_x1_0 += x1[2] * *left++;
		 ump_x1_0 += x1[3] * *left++;	      			
		 ump_x1_0 += x1[0];
		 
		 *uX1++ = ump_x1_0;
		 
		 ump_x1_1 =  x1[1] * *left++;
		 ump_x1_1 += x1[2] * *left++;
		 ump_x1_1 += x1[3] * *left++;	      			
		 ump_x1_1 += x1[0];
		 
		 *uX1++ = ump_x1_1;
		 
		 ump_x1_2 =  x1[1] * *left++;
		 ump_x1_2 += x1[2] * *left++;
		 ump_x1_2 += x1[3] * *left++;		
		 ump_x1_2 += x1[0];
		 
		 *uX1++ = ump_x1_2;
		 
		 ump_x1_3 =  x1[1] * *left++;
		 ump_x1_3 += x1[2] * *left++;
		 ump_x1_3 += x1[3] * *left++;		 
		 ump_x1_3 += x1[0];
		 
		 *uX1++ = ump_x1_3;		
		 
		 ump_x2_0 =  x2[1] * *right++;
		 ump_x2_0 += x2[2] * *right++;
		 ump_x2_0 += x2[3] * *right++;	
		 ump_x2_0 += x2[0];
		 
		 *uX2++ = ump_x2_0;
		 
		 ump_x2_1 =  x2[1] * *right++;
		 ump_x2_1 += x2[2] * *right++;
		 ump_x2_1 += x2[3] * *right++;	
		 ump_x2_1 += x2[0];	 
		 
		 *uX2++ = ump_x2_1;
		 
		 ump_x2_2 =  x2[1] * *right++;
		 ump_x2_2 += x2[2] * *right++;
		 ump_x2_2 += x2[3] * *right++;	       
		 ump_x2_2 += x2[0];	  
		 
		 *uX2++ = ump_x2_2;
		 
		 ump_x2_3 =  x2[1] * *right++;
		 ump_x2_3 += x2[2] * *right++;
		 ump_x2_3 += x2[3] * *right++;		     	       
		 ump_x2_3 += x2[0];	    
		 
		 *uX2++ = ump_x2_3;
		 
		 ump_x1_0 =  x1[1] * *left++;
		 ump_x1_0 += x1[2] * *left++;
		 ump_x1_0 += x1[3] * *left++;	      			
		 ump_x1_0 += x1[0];
		 
		 *uX1++ = ump_x1_0;
		 
		 ump_x1_1 =  x1[1] * *left++;
		 ump_x1_1 += x1[2] * *left++;
		 ump_x1_1 += x1[3] * *left++;	      			
		 ump_x1_1 += x1[0];
		 
		 *uX1++ = ump_x1_1;
		 
		 ump_x1_2 =  x1[1] * *left++;
		 ump_x1_2 += x1[2] * *left++;
		 ump_x1_2 += x1[3] * *left++;		
		 ump_x1_2 += x1[0];
		 
		 *uX1++ = ump_x1_2;
		 
		 ump_x1_3 =  x1[1] * *left++;
		 ump_x1_3 += x1[2] * *left++;
		 ump_x1_3 += x1[3] * *left++;		 
		 ump_x1_3 += x1[0];
		 
		 *uX1++ = ump_x1_3;		
		 
		 ump_x2_0 =  x2[1] * *right++;
		 ump_x2_0 += x2[2] * *right++;
		 ump_x2_0 += x2[3] * *right++;	
		 ump_x2_0 += x2[0];
		 
		 *uX2++ = ump_x2_0;
		 
		 ump_x2_1 =  x2[1] * *right++;
		 ump_x2_1 += x2[2] * *right++;
		 ump_x2_1 += x2[3] * *right++;	
		 ump_x2_1 += x2[0];	 
		 
		 *uX2++ = ump_x2_1;
		 
		 ump_x2_2 =  x2[1] * *right++;
		 ump_x2_2 += x2[2] * *right++;
		 ump_x2_2 += x2[3] * *right++;	       
		 ump_x2_2 += x2[0];	  
		 
		 *uX2++ = ump_x2_2;
		 
		 ump_x2_3 =  x2[1] * *right++;
		 ump_x2_3 += x2[2] * *right++;
		 ump_x2_3 += x2[3] * *right++;		     	       
		 ump_x2_3 += x2[0];	    
		 
		 *uX2++ = ump_x2_3;	       
	       }			
	}
	     	 

	 for (i = lower; i < n; i++) 
	   {	
	     model = modelptr[i];
	     
	     uX1 = &umpX1[256 * model + 16 * tipX1[i]];
	     uX2 = &umpX2[256 * model + 16 * tipX2[i]];
	     x3 = &x3_start[16 * i];

	     EV = &(extEV[model * 16]);
		
	     x1px2 = *uX1++ * *uX2++;
	     x3[0] = x1px2 * EV[0];
	     x3[1] = x1px2 * EV[1];
	     x3[2] = x1px2 * EV[2];
	     x3[3] = x1px2 * EV[3];
	     
	     x1px2 = *uX1++ * *uX2++;
	     x3[0] += x1px2  *  EV[4];
	     x3[1] += x1px2 *   EV[5];
	     x3[2] += x1px2 *  EV[6];
	     x3[3] += x1px2 *   EV[7];
	    
	     x1px2 = *uX1++ * *uX2++;
	     x3[0] += x1px2 *  EV[8];
	     x3[1] += x1px2*   EV[9];
	     x3[2] += x1px2 *  EV[10];
	     x3[3] += x1px2 *  EV[11];
	     
	     x1px2 = *uX1++ * *uX2++;
	     x3[0] += x1px2 *   EV[12];
	     x3[1] += x1px2 *   EV[13];
	     x3[2] += x1px2 *   EV[14];
	     x3[3] += x1px2 *   EV[15];

	     /* rate 1 */
	     	    

	     x1px2 = *uX1++ * *uX2++;/*ump_x1_0 * ump_x2_0;*/
	     x3[4] = x1px2 *  EV[0];
	     x3[5] = x1px2 *  EV[1];
	     x3[6] = x1px2 *  EV[2];
	     x3[7] = x1px2 *  EV[3];
	     
	     x1px2 = *uX1++  * *uX2++;/*ump_x1_1 * ump_x2_1;*/
	     x3[4] += x1px2 * EV[4];
	     x3[5] += x1px2 * EV[5];
	     x3[6] += x1px2 * EV[6];
	     x3[7] += x1px2 * EV[7];
	     
	     x1px2 = *uX1++ * *uX2++;/*ump_x1_2 * ump_x2_2;*/
	     x3[4] += x1px2 * EV[8];
	     x3[5] += x1px2 * EV[9];
	     x3[6] += x1px2 * EV[10];
	     x3[7] += x1px2 * EV[11];
	    
	     x1px2 = *uX1++ * *uX2++;/*ump_x1_3 * ump_x2_3;*/
	     x3[4] += x1px2 *   EV[12];
	     x3[5] += x1px2 *   EV[13];
	     x3[6] += x1px2 *   EV[14];
	     x3[7] += x1px2 *   EV[15];

	     /* rate 2 */	     	     	     

	     x1px2 = *uX1++ * *uX2++;/*ump_x1_0 * ump_x2_0;*/
	     x3[8] = x1px2 *  EV[0];
	     x3[9] = x1px2 *  EV[1];
	     x3[10] = x1px2 * EV[2];
	     x3[11] = x1px2 *  EV[3];
	
	     x1px2 = *uX1++ * *uX2++;/*ump_x1_1 * ump_x2_1;*/
	     x3[8] += x1px2  *  EV[4];
	     x3[9] += x1px2 *   EV[5];
	     x3[10] += x1px2 *  EV[6];
	     x3[11] += x1px2 *   EV[7];
	    
	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_2 * ump_x2_2;*/
	     x3[8] += x1px2 *  EV[8];
	     x3[9] += x1px2*   EV[9];
	     x3[10] += x1px2 *  EV[10];
	     x3[11] += x1px2 *  EV[11];
	     
	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_3 * ump_x2_3;*/
	     x3[8] += x1px2 *   EV[12];
	     x3[9] += x1px2 *   EV[13];
	     x3[10] += x1px2 *   EV[14];
	     x3[11] += x1px2 *   EV[15];

	     /* rate 3 */	    	    

	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_0 * ump_x2_0;*/
	     x3[12] = x1px2 *  EV[0];
	     x3[13] = x1px2 *  EV[1];
	     x3[14] = x1px2 * EV[2];
	     x3[15] = x1px2 *  EV[3];
	     
	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_1 * ump_x2_1;*/
	     x3[12] += x1px2  *  EV[4];
	     x3[13] += x1px2 *   EV[5];
	     x3[14] += x1px2 *  EV[6];
	     x3[15] += x1px2 *   EV[7];
	    
	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_2 * ump_x2_2;*/
	     x3[12] += x1px2 *  EV[8];
	     x3[13] += x1px2*   EV[9];
	     x3[14] += x1px2 *  EV[10];
	     x3[15] += x1px2 *  EV[11];
	     
	     x1px2 =  *uX1++ * *uX2++;/*ump_x1_3 * ump_x2_3;*/
	     x3[12] += x1px2 *   EV[12];
	     x3[13] += x1px2 *   EV[13];
	     x3[14] += x1px2 *   EV[14];
	     x3[15] += x1px2 *   EV[15];

	     /********************************************************************************/

	     ex3[i] = 0;		  		    	     	       	       
	   }
         	
	  free(umpX1);
	  free(umpX2);	
      }
      break;    
    case TIP_INNER:
      {	 	 	   
	 double *uX1, *umpX1;	
	     
	 umpX1 = (double *)malloc(256 * numberOfModels * sizeof(double));
	 
	 for(model = 0; model < numberOfModels; model++)
	   {
	     uX1 = &umpX1[256 * model + 16];

	     for(i = 1; i < 16; i++)
	       {	    	     		 
		 x1 = &(tipVector[model * 64 + i * 4]);	     
		 left = &left_start[model * 48];
		 
		 ump_x1_0 =  x1[1] * *left++;
		 ump_x1_0 += x1[2] * *left++;
		 ump_x1_0 += x1[3] * *left++;	       
		 ump_x1_0 += x1[0];
		 
		 *uX1++ = ump_x1_0;
		 
		 ump_x1_1 =  x1[1] * *left++;
		 ump_x1_1 += x1[2] * *left++;
		 ump_x1_1 += x1[3] * *left++;		 
		 ump_x1_1 += x1[0];	 
	     
		 *uX1++ = ump_x1_1;
		 
		 ump_x1_2 =  x1[1] * *left++;
		 ump_x1_2 += x1[2] * *left++;
		 ump_x1_2 += x1[3] * *left++;		
		 ump_x1_2 += x1[0];	  
	     
		 *uX1++ = ump_x1_2;
		 
		 ump_x1_3 =  x1[1] * *left++;
		 ump_x1_3 += x1[2] * *left++;
		 ump_x1_3 += x1[3] * *left++;		     	      
		 ump_x1_3 += x1[0];	    
		 
		 *uX1++ = ump_x1_3;
		 
		 ump_x1_0 =  x1[1] * *left++;
		 ump_x1_0 += x1[2] * *left++;
		 ump_x1_0 += x1[3] * *left++;	       
		 ump_x1_0 += x1[0];
		 
		 *uX1++ = ump_x1_0;
		 
		 ump_x1_1 =  x1[1] * *left++;
		 ump_x1_1 += x1[2] * *left++;
		 ump_x1_1 += x1[3] * *left++;		 
		 ump_x1_1 += x1[0];	 
		 
		 *uX1++ = ump_x1_1;
		 
		 ump_x1_2 =  x1[1] * *left++;
		 ump_x1_2 += x1[2] * *left++;
		 ump_x1_2 += x1[3] * *left++;		
		 ump_x1_2 += x1[0];	  
		 
		 *uX1++ = ump_x1_2;
		 
		 ump_x1_3 =  x1[1] * *left++;
		 ump_x1_3 += x1[2] * *left++;
		 ump_x1_3 += x1[3] * *left++;		     	      
		 ump_x1_3 += x1[0];	    
		 
		 *uX1++ = ump_x1_3;

		 ump_x1_0 =  x1[1] * *left++;
		 ump_x1_0 += x1[2] * *left++;
		 ump_x1_0 += x1[3] * *left++;	       
		 ump_x1_0 += x1[0];
		 
		 *uX1++ = ump_x1_0;
		 
		 ump_x1_1 =  x1[1] * *left++;
		 ump_x1_1 += x1[2] * *left++;
		 ump_x1_1 += x1[3] * *left++;		 
		 ump_x1_1 += x1[0];	 
		 
		 *uX1++ = ump_x1_1;
		 
		 ump_x1_2 =  x1[1] * *left++;
		 ump_x1_2 += x1[2] * *left++;
		 ump_x1_2 += x1[3] * *left++;		
		 ump_x1_2 += x1[0];	  
		 
		 *uX1++ = ump_x1_2;
		 
		 ump_x1_3 =  x1[1] * *left++;
		 ump_x1_3 += x1[2] * *left++;
		 ump_x1_3 += x1[3] * *left++;		     	      
		 ump_x1_3 += x1[0];	    
		 
		 *uX1++ = ump_x1_3;
		 
		 ump_x1_0 =  x1[1] * *left++;
		 ump_x1_0 += x1[2] * *left++;
		 ump_x1_0 += x1[3] * *left++;	       
		 ump_x1_0 += x1[0];
		 
		 *uX1++ = ump_x1_0;
		 
		 ump_x1_1 =  x1[1] * *left++;
		 ump_x1_1 += x1[2] * *left++;
		 ump_x1_1 += x1[3] * *left++;		 
		 ump_x1_1 += x1[0];	 
		 
		 *uX1++ = ump_x1_1;
		 
		 ump_x1_2 =  x1[1] * *left++;
		 ump_x1_2 += x1[2] * *left++;
		 ump_x1_2 += x1[3] * *left++;		
		 ump_x1_2 += x1[0];	  
		 
		 *uX1++ = ump_x1_2;
		 
		 ump_x1_3 =  x1[1] * *left++;
		 ump_x1_3 += x1[2] * *left++;
		 ump_x1_3 += x1[3] * *left++;		     	      
		 ump_x1_3 += x1[0];	    
		 
		 *uX1++ = ump_x1_3;	    
	       }
	   }
	 
	 
	 for (i = lower; i < n; i++) 
	   {	
	     model = modelptr[i];
	     right = &right_start[model * 48];
	     uX1 = &umpX1[model * 256 + 16 * tipX1[i]];	       

	     x2 = &x2_start[16 * i];
	     x3 = &x3_start[16 * i];
	     EV = &(extEV[model * 16]);						 	 	    	    	     	    

	     ump_x2_0 =  x2[1] * *right++;
	     ump_x2_0 += x2[2] * *right++;
	     ump_x2_0 += x2[3] * *right++;		     	   
	     ump_x2_0 += x2[0];
	     
	     ump_x2_1 =  x2[1] * *right++;
	     ump_x2_1 += x2[2] * *right++;
	     ump_x2_1 += x2[3] * *right++;		     	    
	     ump_x2_1 += x2[0];	 
		 
	     ump_x2_2 = x2[1] * *right++;
	     ump_x2_2 += x2[2] * *right++;
	     ump_x2_2 +=  x2[3] * *right++;		     	     
	     ump_x2_2 += x2[0];	  
	     
	     ump_x2_3 = x2[1] * *right++;
	     ump_x2_3 += x2[2] * *right++;
	     ump_x2_3 += x2[3] * *right++;		     	     
	     ump_x2_3 += x2[0];	    	  		   	  	    
		 
	     x1px2 = *uX1++ * ump_x2_0;
	     x3[0] = x1px2 *  EV[0];
	     x3[1] = x1px2 *  EV[1];
	     x3[2] = x1px2 * EV[2];
	     x3[3] = x1px2 *  EV[3];
	     
	     x1px2 = *uX1++ * ump_x2_1;
	     x3[0] += x1px2  *  EV[4];
	     x3[1] += x1px2 *   EV[5];
	     x3[2] += x1px2 *  EV[6];
	     x3[3] += x1px2 *   EV[7];
		 
	     x1px2 = *uX1++ * ump_x2_2;
	     x3[0] += x1px2 *  EV[8];
	     x3[1] += x1px2*   EV[9];
	     x3[2] += x1px2 *  EV[10];
	     x3[3] += x1px2 *  EV[11];
	     
	     x1px2 = *uX1++ * ump_x2_3;
	     x3[0] += x1px2 *   EV[12];
	     x3[1] += x1px2 *   EV[13];
	     x3[2] += x1px2 *   EV[14];
	     x3[3] += x1px2 *   EV[15];		 		        

	     ump_x2_0 = x2[5] * *right++;
	     ump_x2_0 += x2[6] * *right++;
	     ump_x2_0 += x2[7] * *right++;		     	  
	     ump_x2_0 += x2[4];
	     
	     ump_x2_1 = x2[5] * *right++;
	     ump_x2_1 += x2[6] * *right++;
	     ump_x2_1 +=  x2[7] * *right++;		     	   
	     ump_x2_1 += x2[4];	 
	     
	     ump_x2_2 = x2[5] * *right++;
	     ump_x2_2 += x2[6] * *right++;
	     ump_x2_2 +=  x2[7] * *right++;		     	  
	     ump_x2_2 += x2[4];	  
	     
	     ump_x2_3 = x2[5] * *right++;
	     ump_x2_3 += x2[6] * *right++;
	     ump_x2_3 += x2[7] * *right++;		     	 
	     ump_x2_3 += x2[4];	    	  		   	  	    
	     
	     x1px2 = *uX1++ * ump_x2_0;
	     x3[4] = x1px2 *  EV[0];
	     x3[5] = x1px2 *  EV[1];
	     x3[6] = x1px2 * EV[2];
	     x3[7] = x1px2 *  EV[3];
	     
	     x1px2 = *uX1++ * ump_x2_1;
	     x3[4] += x1px2  *  EV[4];
	     x3[5] += x1px2 *   EV[5];
	     x3[6] += x1px2 *  EV[6];
	     x3[7] += x1px2 *   EV[7];
	     
	     x1px2 = *uX1++ * ump_x2_2;
	     x3[4] += x1px2 *  EV[8];
	     x3[5] += x1px2*   EV[9];
	     x3[6] += x1px2 *  EV[10];
	     x3[7] += x1px2 *  EV[11];
	     
	     x1px2 = *uX1++ * ump_x2_3;
	     x3[4] += x1px2 *   EV[12];
	     x3[5] += x1px2 *   EV[13];
	     x3[6] += x1px2 *   EV[14];
	     x3[7] += x1px2 *   EV[15];	     	   
		 
	     ump_x2_0 = x2[9] * *right++;
	     ump_x2_0 += x2[10] * *right++;
	     ump_x2_0 += x2[11] * *right++;		     	   
	     ump_x2_0 += x2[8];
	     
	     ump_x2_1 = x2[9] * *right++;
	     ump_x2_1 += x2[10] * *right++;
	     ump_x2_1 +=  x2[11] * *right++;		     	    
	     ump_x2_1 += x2[8];	 
	     
	     ump_x2_2 = x2[9] * *right++;
	     ump_x2_2 += x2[10] * *right++;
	     ump_x2_2 +=  x2[11] * *right++;		     	     
	     ump_x2_2 += x2[8];	  
	     
	     ump_x2_3 = x2[9] * *right++;
	     ump_x2_3 += x2[10] * *right++;
	     ump_x2_3 += x2[11] * *right++;		     	   
	     ump_x2_3 += x2[8];	    	  		   	  	    
	     
	     x1px2 = *uX1++ * ump_x2_0;
	     x3[8] = x1px2 *  EV[0];
	     x3[9] = x1px2 *  EV[1];
	     x3[10] = x1px2 * EV[2];
	     x3[11] = x1px2 *  EV[3];
	     
	     x1px2 = *uX1++ * ump_x2_1;
	     x3[8] += x1px2  *  EV[4];
	     x3[9] += x1px2 *   EV[5];
	     x3[10] += x1px2 *  EV[6];
	     x3[11] += x1px2 *   EV[7];
	     
	     x1px2 = *uX1++ * ump_x2_2;
	     x3[8] += x1px2 *  EV[8];
	     x3[9] += x1px2*   EV[9];
	     x3[10] += x1px2 *  EV[10];
	     x3[11] += x1px2 *  EV[11];
	     
	     x1px2 = *uX1++ * ump_x2_3;
	     x3[8] += x1px2 *   EV[12];
	     x3[9] += x1px2 *   EV[13];
	     x3[10] += x1px2 *   EV[14];
	     x3[11] += x1px2 *   EV[15];
	     
      	     ump_x2_0 = x2[13] * *right++;
	     ump_x2_0 += x2[14] * *right++;
	     ump_x2_0 += x2[15] * *right++;		     	   
	     ump_x2_0 += x2[12];
	     
	     ump_x2_1 = x2[13] * *right++;
	     ump_x2_1 += x2[14] * *right++;
	     ump_x2_1 +=  x2[15] * *right++;		     	    
	     ump_x2_1 += x2[12];	 
	     
	     ump_x2_2 = x2[13] * *right++;
	     ump_x2_2 += x2[14] * *right++;
	     ump_x2_2 +=  x2[15] * *right++;		     	     
	     ump_x2_2 += x2[12];	  
	     
	     ump_x2_3 = x2[13] * *right++;
	     ump_x2_3 += x2[14] * *right++;
	     ump_x2_3 += x2[15] * *right++;		     	    
	     ump_x2_3 += x2[12];	    	  		   	  	    
	     
	     x1px2 = *uX1++ * ump_x2_0;
	     x3[12] = x1px2 *  EV[0];
	     x3[13] = x1px2 *  EV[1];
	     x3[14] = x1px2 * EV[2];
	     x3[15] = x1px2 *  EV[3];
	     
	     x1px2 = *uX1++ * ump_x2_1;
	     x3[12] += x1px2  *  EV[4];
	     x3[13] += x1px2 *   EV[5];
	     x3[14] += x1px2 *  EV[6];
	     x3[15] += x1px2 *   EV[7];
	     
	     x1px2 = *uX1++ * ump_x2_2;
	     x3[12] += x1px2 *  EV[8];
	     x3[13] += x1px2*   EV[9];
	     x3[14] += x1px2 *  EV[10];
	     x3[15] += x1px2 *  EV[11];
	     
	     x1px2 = *uX1++ * ump_x2_3;
	     x3[12] += x1px2 *   EV[12];
	     x3[13] += x1px2 *   EV[13];
	     x3[14] += x1px2 *   EV[14];
	     x3[15] += x1px2 *   EV[15];
	     
	     ex3[i] = ex2[i];		  
	     
	     if (ABS(x3[0]) < minlikelihood && ABS(x3[2]) < minlikelihood && ABS(x3[1]) < minlikelihood && ABS(x3[3]) < minlikelihood &&
		 ABS(x3[4]) < minlikelihood && ABS(x3[6]) < minlikelihood && ABS(x3[5]) < minlikelihood && ABS(x3[7]) < minlikelihood &&
		 ABS(x3[8]) < minlikelihood && ABS(x3[10]) < minlikelihood && ABS(x3[9]) < minlikelihood && ABS(x3[11]) < minlikelihood &&
		 ABS(x3[12]) < minlikelihood && ABS(x3[14]) < minlikelihood && ABS(x3[13]) < minlikelihood && ABS(x3[15]) < minlikelihood)	         
	       {	     	    	 		    
		 
		 x3[0]   *= twotothe256;
		 x3[1]   *= twotothe256;
		 x3[2]   *= twotothe256;
		 x3[3]   *= twotothe256;
		 
		 x3[4]   *= twotothe256;
		 x3[5]   *= twotothe256;
		 x3[6]   *= twotothe256;
		 x3[7]   *= twotothe256;
		 
		 x3[8]   *= twotothe256;
		 x3[9]   *= twotothe256;
		 x3[10]   *= twotothe256;
		 x3[11]   *= twotothe256;
		 
		 x3[12]   *= twotothe256;
		 x3[13]   *= twotothe256;
		 x3[14]   *= twotothe256;
		 x3[15]   *= twotothe256;
		 
		 ex3[i] += 1;
	       }	 	      	   
	   }
	 	
	 free(umpX1);	 
      }
      break;
    case INNER_INNER:
    
     for (i = lower; i < n; i++) 
       {	
	 model = modelptr[i];
	
	 left  = &left_start[model * 48];
	 right = &right_start[model * 48];
		      
	 x1 = &x1_start[16 * i];
	 x2 = &x2_start[16 * i];
	 x3 = &x3_start[16 * i];
	 		        
	 EV = &(extEV[model * 16]);	

	 ump_x1_0  = x1[1] * *left++;
	 ump_x1_0 += x1[2] * *left++;
	 ump_x1_0 += x1[3] * *left++;	      	      
	 ump_x1_0 += x1[0];
		 
	 ump_x1_1  = x1[1] * *left++;
	 ump_x1_1 += x1[2] * *left++;
	 ump_x1_1 += x1[3] * *left++;	      	       
	 ump_x1_1 += x1[0];
	 
	 ump_x1_2  = x1[1] * *left++;
	 ump_x1_2 += x1[2] * *left++;
	 ump_x1_2 += x1[3] * *left++;	      		
	 ump_x1_2 += x1[0];
		 
	 ump_x1_3  = x1[1] * *left++;
	 ump_x1_3 += x1[2] * *left++;
	 ump_x1_3 += x1[3] * *left++;	      	      
	 ump_x1_3 += x1[0];
		 	 	 
	 ump_x2_0 = x2[1] * *right++;
	 ump_x2_0 += x2[2] * *right++;
	 ump_x2_0 += x2[3] * *right++;		     	
	 ump_x2_0 += x2[0];
	 
	 ump_x2_1 = x2[1] * *right++;
	 ump_x2_1 += x2[2] * *right++;
	 ump_x2_1 +=  x2[3] * *right++;		     	
	 ump_x2_1 += x2[0];	 
	 
	 ump_x2_2 = x2[1] * *right++;
	 ump_x2_2 += x2[2] * *right++;
	 ump_x2_2 +=  x2[3] * *right++;		     	
	 ump_x2_2 += x2[0];	  
	 
	 ump_x2_3 = x2[1] * *right++;
	 ump_x2_3 += x2[2] * *right++;
	 ump_x2_3 += x2[3] * *right++;		     	
	 ump_x2_3 += x2[0];	    	  		   	  	    
	 
	 x1px2 = ump_x1_0 * ump_x2_0;
	 x3[0] = x1px2 *  EV[0];
	 x3[1] = x1px2 *  EV[1];
	 x3[2] = x1px2 * EV[2];
	 x3[3] = x1px2 *  EV[3];
	 
	 x1px2 = ump_x1_1 * ump_x2_1;
	 x3[0] += x1px2  *  EV[4];
	 x3[1] += x1px2 *   EV[5];
	 x3[2] += x1px2 *  EV[6];
	 x3[3] += x1px2 *   EV[7];
	 
	 x1px2 = ump_x1_2 * ump_x2_2;
	 x3[0] += x1px2 *  EV[8];
	 x3[1] += x1px2*   EV[9];
	 x3[2] += x1px2 *  EV[10];
	 x3[3] += x1px2 *  EV[11];
	 
	 x1px2 = ump_x1_3 * ump_x2_3;
	 x3[0] += x1px2 *   EV[12];
	 x3[1] += x1px2 *   EV[13];
	 x3[2] += x1px2 *   EV[14];
	 x3[3] += x1px2 *   EV[15];

	 /* rate 1 */	

	 ump_x1_0 = x1[5] * *left++;
	 ump_x1_0 += x1[6] * *left++;
	 ump_x1_0 += x1[7]* *left++;	      	     
	 ump_x1_0 += x1[4];
	 
	 ump_x1_1 =  x1[5] * *left++;
	 ump_x1_1 += x1[6] * *left++;
	 ump_x1_1 += x1[7]* *left++;	      	      
	 ump_x1_1 += x1[4];
	 
	 ump_x1_2 = x1[5] * *left++;
	 ump_x1_2 += x1[6] * *left++;
	 ump_x1_2 += x1[7] * *left++;	      	      
	 ump_x1_2 += x1[4];
	 
	 ump_x1_3 =  x1[5] * *left++;
	 ump_x1_3 += x1[6] * *left++;
	 ump_x1_3 += x1[7] * *left++;	      	
	 ump_x1_3 += x1[4];
	 
	 ump_x2_0 = x2[5] * *right++;
	 ump_x2_0 += x2[6] * *right++;
	 ump_x2_0 += x2[7] * *right++;		           
	 ump_x2_0 += x2[4];
	 
	 ump_x2_1 = x2[5] * *right++;
	 ump_x2_1 += x2[6] * *right++;
	 ump_x2_1 +=  x2[7] * *right++;		           
	 ump_x2_1 += x2[4];	 
	 
	 ump_x2_2 = x2[5] * *right++;
	 ump_x2_2 += x2[6] * *right++;
	 ump_x2_2 +=  x2[7] * *right++;		           
	 ump_x2_2 += x2[4];	  
		 
	 ump_x2_3 = x2[5] * *right++;
	 ump_x2_3 += x2[6] * *right++;
	 ump_x2_3 += x2[7] * *right++;		            
	 ump_x2_3 += x2[4];	    	  		   	  	    
	 
	 x1px2 = ump_x1_0 * ump_x2_0;
	 x3[4] = x1px2 *  EV[0];
	 x3[5] = x1px2 *  EV[1];
	 x3[6] = x1px2 * EV[2];
	 x3[7] = x1px2 *  EV[3];
	 
	 x1px2 = ump_x1_1 * ump_x2_1;
	 x3[4] += x1px2  *  EV[4];
	 x3[5] += x1px2 *   EV[5];
	 x3[6] += x1px2 *  EV[6];
	 x3[7] += x1px2 *   EV[7];
	 
	 x1px2 = ump_x1_2 * ump_x2_2;
	 x3[4] += x1px2 *  EV[8];
	 x3[5] += x1px2*   EV[9];
	 x3[6] += x1px2 *  EV[10];
	 x3[7] += x1px2 *  EV[11];
	 
	 x1px2 = ump_x1_3 * ump_x2_3;
	 x3[4] += x1px2 *   EV[12];
	 x3[5] += x1px2 *   EV[13];
	 x3[6] += x1px2 *   EV[14];
	 x3[7] += x1px2 *   EV[15];
	 
	 /* rate 2 */	 	
	 
	 ump_x1_0 = x1[9] * *left++;
	 ump_x1_0 += x1[10] * *left++;
	 ump_x1_0 += x1[11] * *left++;	      	     
	 ump_x1_0 += x1[8];
	 
	 ump_x1_1 =  x1[9] * *left++;
	 ump_x1_1 += x1[10] * *left++;
	 ump_x1_1 += x1[11]* *left++;	      	      
	 ump_x1_1 += x1[8];
	 
	 ump_x1_2 = x1[9] * *left++;
	 ump_x1_2 += x1[10] * *left++;
	 ump_x1_2 += x1[11] * *left++;	      	     
	 ump_x1_2 += x1[8];
	 
	 ump_x1_3 =  x1[9] * *left++;
	 ump_x1_3 += x1[10] * *left++;
	 ump_x1_3 += x1[11] * *left++;	      	     
	 ump_x1_3 += x1[8];
	 	 
	 ump_x2_0 = x2[9] * *right++;
	 ump_x2_0 += x2[10] * *right++;
	 ump_x2_0 += x2[11] * *right++;		     
	 ump_x2_0 += x2[8];
	
	 ump_x2_1 = x2[9] * *right++;
	 ump_x2_1 += x2[10] * *right++;
	 ump_x2_1 +=  x2[11] * *right++;		         
	 ump_x2_1 += x2[8];	 
	 
	 ump_x2_2 = x2[9] * *right++;
	 ump_x2_2 += x2[10] * *right++;
	 ump_x2_2 +=  x2[11] * *right++;		            
	 ump_x2_2 += x2[8];	  
	 
	 ump_x2_3 = x2[9] * *right++;
	 ump_x2_3 += x2[10] * *right++;
	 ump_x2_3 += x2[11] * *right++;		          
	 ump_x2_3 += x2[8];	    	  		   	  	    
	 
	 x1px2 = ump_x1_0 * ump_x2_0;
	 x3[8] = x1px2 *  EV[0];
	 x3[9] = x1px2 *  EV[1];
	 x3[10] = x1px2 * EV[2];
	 x3[11] = x1px2 *  EV[3];
	 
	 x1px2 = ump_x1_1 * ump_x2_1;
	 x3[8] += x1px2  *  EV[4];
	 x3[9] += x1px2 *   EV[5];
	 x3[10] += x1px2 *  EV[6];
	 x3[11] += x1px2 *   EV[7];
	 
	 x1px2 = ump_x1_2 * ump_x2_2;
	 x3[8] += x1px2 *  EV[8];
	 x3[9] += x1px2*   EV[9];
	 x3[10] += x1px2 *  EV[10];
	 x3[11] += x1px2 *  EV[11];
	 
	 x1px2 = ump_x1_3 * ump_x2_3;
	 x3[8] += x1px2 *   EV[12];
	 x3[9] += x1px2 *   EV[13];
	 x3[10] += x1px2 *   EV[14];
	 x3[11] += x1px2 *   EV[15];
	 
	 /* rate 3 */	

	 ump_x1_0 = x1[13] * *left++;
	 ump_x1_0 += x1[14] * *left++;
	 ump_x1_0 += x1[15]* *left++;	      	     
	 ump_x1_0 += x1[12];
	 
	 ump_x1_1 =  x1[13] * *left++;
	 ump_x1_1 += x1[14] * *left++;
	 ump_x1_1 += x1[15]* *left++;	      	 
	 ump_x1_1 += x1[12];
	 
	 ump_x1_2 = x1[13] * *left++;
	 ump_x1_2 += x1[14] * *left++;
	 ump_x1_2 += x1[15] * *left++;	      	      
	 ump_x1_2 += x1[12];
	 
	 ump_x1_3 =  x1[13] * *left++;
	 ump_x1_3 += x1[14] * *left++;
	 ump_x1_3 += x1[15] * *left++;	      	     
	 ump_x1_3 += x1[12];
	 
	 ump_x2_0 = x2[13] * *right++;
	 ump_x2_0 += x2[14] * *right++;
	 ump_x2_0 += x2[15] * *right++;		         
	 ump_x2_0 += x2[12];
	 
	 ump_x2_1 = x2[13] * *right++;
	 ump_x2_1 += x2[14] * *right++;
	 ump_x2_1 +=  x2[15] * *right++;		           
	 ump_x2_1 += x2[12];	 
	 
	 ump_x2_2 = x2[13] * *right++;
	 ump_x2_2 += x2[14] * *right++;
	 ump_x2_2 +=  x2[15] * *right++;		         
	 ump_x2_2 += x2[12];	  
	 
	 ump_x2_3 = x2[13] * *right++;
	 ump_x2_3 += x2[14] * *right++;
	 ump_x2_3 += x2[15] * *right++;		            
	 ump_x2_3 += x2[12];	    	  		   	  	    
	 
	 x1px2 = ump_x1_0 * ump_x2_0;
	 x3[12] = x1px2 *  EV[0];
	 x3[13] = x1px2 *  EV[1];
	 x3[14] = x1px2 * EV[2];
	 x3[15] = x1px2 *  EV[3];
	 
	 x1px2 = ump_x1_1 * ump_x2_1;
	 x3[12] += x1px2  *  EV[4];
	 x3[13] += x1px2 *   EV[5];
	 x3[14] += x1px2 *  EV[6];
	 x3[15] += x1px2 *   EV[7];
	 
	 x1px2 = ump_x1_2 * ump_x2_2;
	 x3[12] += x1px2 *  EV[8];
	 x3[13] += x1px2*   EV[9];
	 x3[14] += x1px2 *  EV[10];
	 x3[15] += x1px2 *  EV[11];
	 
	 x1px2 = ump_x1_3 * ump_x2_3;
	 x3[12] += x1px2 *   EV[12];
	 x3[13] += x1px2 *   EV[13];
	 x3[14] += x1px2 *   EV[14];
	 x3[15] += x1px2 *   EV[15];
	 
	 /********************************************************************************/
	 
	 ex3[i] = ex1[i] + ex2[i];		  
	 
	 if (ABS(x3[0]) < minlikelihood && ABS(x3[2]) < minlikelihood && ABS(x3[1]) < minlikelihood && ABS(x3[3]) < minlikelihood &&
	     ABS(x3[4]) < minlikelihood && ABS(x3[6]) < minlikelihood && ABS(x3[5]) < minlikelihood && ABS(x3[7]) < minlikelihood &&
	     ABS(x3[8]) < minlikelihood && ABS(x3[10]) < minlikelihood && ABS(x3[9]) < minlikelihood && ABS(x3[11]) < minlikelihood &&
	     ABS(x3[12]) < minlikelihood && ABS(x3[14]) < minlikelihood && ABS(x3[13]) < minlikelihood && ABS(x3[15]) < minlikelihood)	
	   {			    
	     x3[0]   *= twotothe256;
	     x3[1]   *= twotothe256;
	     x3[2]   *= twotothe256;
	     x3[3]   *= twotothe256;
	     
	     x3[4]   *= twotothe256;
	     x3[5]   *= twotothe256;
	     x3[6]   *= twotothe256;
	     x3[7]   *= twotothe256;
	     
	     x3[8]   *= twotothe256;
	     x3[9]   *= twotothe256;
	     x3[10]   *= twotothe256;
	     x3[11]   *= twotothe256;
	     
	     x3[12]   *= twotothe256;
	     x3[13]   *= twotothe256;
	     x3[14]   *= twotothe256;
	     x3[15]   *= twotothe256;
	     
	     ex3[i] += 1;
	   }
       }
     break;
    default:
      assert(0);
    }
     
  free(left_start); 
  free(right_start); 
}




static void newviewGTRCATPROT(traversalInfo *ti, double *extEV,  double *extEI,  double *EIGN,
				 double *rptr,  int *cptr,
				 double *x1, double *x2, double *x3, double *tipVector,				
				 int    *ex1, int *ex2, int *ex3, char *tipX1, char *tipX2,
				 int lower, int n,  int numberOfCategories, double z1, double z2)
{
  double  
    *left, *left_start, *right, *right_start, *EV, *EI, *v;
  double      
    ump_x1[20], ump_x2[20],  
    d1[19], d2[19], lz1[19], lz2[19];   
  double
    x1px2, ki;
  int i, l, scale;	                                           
       
  left_start  = left  = (double *)malloc(380 * numberOfCategories * sizeof(double));
  right_start = right = (double *)malloc(380 * numberOfCategories * sizeof(double));
  
  for(i = 0; i < 19; i++)
    lz1[i] = EIGN[i] * z1;
          
  for(i = 0; i < 19; i++)
    lz2[i] = EIGN[i] * z2;      
       
  for(i = 0; i < numberOfCategories; i++)
    {	
      ki = rptr[i];
	   	   	   
      d1[0] = exp (ki * lz1[0]);
      d1[1] = exp (ki * lz1[1]);
      d1[2] = exp (ki * lz1[2]);
      d1[3] = exp (ki * lz1[3]);
      d1[4] = exp (ki * lz1[4]);
      d1[5] = exp (ki * lz1[5]);
      d1[6] = exp (ki * lz1[6]);
      d1[7] = exp (ki * lz1[7]);
      d1[8] = exp (ki * lz1[8]);
      d1[9] = exp (ki * lz1[9]);
      d1[10] = exp (ki * lz1[10]);
      d1[11] = exp (ki * lz1[11]);
      d1[12] = exp (ki * lz1[12]);
      d1[13] = exp (ki * lz1[13]);
      d1[14] = exp (ki * lz1[14]);
      d1[15] = exp (ki * lz1[15]);
      d1[16] = exp (ki * lz1[16]);
      d1[17] = exp (ki * lz1[17]);
      d1[18] = exp (ki * lz1[18]);
      
      d2[0] = exp (ki * lz2[0]);
      d2[1] = exp (ki * lz2[1]);
      d2[2] = exp (ki * lz2[2]);
      d2[3] = exp (ki * lz2[3]);
      d2[4] = exp (ki * lz2[4]);
      d2[5] = exp (ki * lz2[5]);
      d2[6] = exp (ki * lz2[6]);
      d2[7] = exp (ki * lz2[7]);
      d2[8] = exp (ki * lz2[8]);
      d2[9] = exp (ki * lz2[9]);
      d2[10] = exp (ki * lz2[10]);
      d2[11] = exp (ki * lz2[11]);
      d2[12] = exp (ki * lz2[12]);
      d2[13] = exp (ki * lz2[13]);
      d2[14] = exp (ki * lz2[14]);
      d2[15] = exp (ki * lz2[15]);
      d2[16] = exp (ki * lz2[16]);
      d2[17] = exp (ki * lz2[17]);
      d2[18] = exp (ki * lz2[18]);	   	      	 
	   	   	 
      EI = extEI;
  	         
      for(l = 0; l < 20; l++)
	{	       
	 *left++ = d1[0] * *EI++;	     	     		   	      
	 *left++ = d1[1] * *EI++;	 
	 *left++ = d1[2] * *EI++;	     	     		   	      
	 *left++ = d1[3] * *EI++;	 
	 *left++ = d1[4] * *EI++;	     	     		   	      
	 *left++ = d1[5] * *EI++;	 
	 *left++ = d1[6] * *EI++;	     	     		   	      
	 *left++ = d1[7] * *EI++;
	 *left++ = d1[8] * *EI++;	     	     		   	      
	 *left++ = d1[9] * *EI++;	 
	 *left++ = d1[10] * *EI++;	     	     		   	      
	 *left++ = d1[11] * *EI++;	 
	 *left++ = d1[12] * *EI++;	     	     		   	      
	 *left++ = d1[13] * *EI++;	 
	 *left++ = d1[14] * *EI++;	     	     		   	      
	 *left++ = d1[15] * *EI++;
	 *left++ = d1[16] * *EI++;	     	     		   	      
	 *left++ = d1[17] * *EI++;	 
	 *left++ = d1[18] * *EI++;	     	     		   	      	      
	}

      EI = extEI;
      
      for(l = 0; l < 20; l++)
	{
	  *right++ = d2[0] * *EI++;	     	     		   	      
	  *right++ = d2[1] * *EI++;	 
	  *right++ = d2[2] * *EI++;	     	     		   	      
	  *right++ = d2[3] * *EI++;	 
	  *right++ = d2[4] * *EI++;	     	     		   	      
	  *right++ = d2[5] * *EI++;	 
	  *right++ = d2[6] * *EI++;	     	     		   	      
	  *right++ = d2[7] * *EI++;
	  *right++ = d2[8] * *EI++;	     	     		   	      
	  *right++ = d2[9] * *EI++;	 
	  *right++ = d2[10] * *EI++;	     	     		   	      
	  *right++ = d2[11] * *EI++;	 
	  *right++ = d2[12] * *EI++;	     	     		   	      
	  *right++ = d2[13] * *EI++;	 
	  *right++ = d2[14] * *EI++;	     	     		   	      
	  *right++ = d2[15] * *EI++;
	  *right++ = d2[16] * *EI++;	     	     		   	      
	  *right++ = d2[17] * *EI++;	 
	  *right++ = d2[18] * *EI++;	     		   	      
	}	  
    }
  
  switch(ti->tipCase)	 
    {
    case TIP_TIP:       
      {	  	 	 	   	   	        

	for (i = lower; i < n; i++) 
	  {		   	    	   	   	     	  
	    left = &left_start[cptr[i] * 380];	      
	    
	    v = &(tipVector[20 * tipX1[i]]);	     	       	       
	    
	    for(l = 0; l < 20; l++)
	      {
		ump_x1[l] = v[0];		 		  		     
		ump_x1[l] += v[1] * *left++;
		ump_x1[l] += v[2] * *left++;
		ump_x1[l] += v[3] * *left++;
		ump_x1[l] += v[4] * *left++;
		ump_x1[l] += v[5] * *left++;
		ump_x1[l] += v[6] * *left++;
		ump_x1[l] += v[7] * *left++;
		ump_x1[l] += v[8] * *left++;
		ump_x1[l] += v[9] * *left++;
		ump_x1[l] += v[10] * *left++;
		ump_x1[l] += v[11] * *left++;
		ump_x1[l] += v[12] * *left++;
		ump_x1[l] += v[13] * *left++;
		ump_x1[l] += v[14] * *left++;
		ump_x1[l] += v[15] * *left++;
		ump_x1[l] += v[16] * *left++;
		ump_x1[l] += v[17] * *left++;
		ump_x1[l] += v[18] * *left++;
		ump_x1[l] += v[19] * *left++;		       		     
	      }
	    
	    v = &(tipVector[20 * tipX2[i]]);
	    
	    left = &right_start[cptr[i] * 380];
	    
	    for(l = 0; l < 20; l++)
	      {
		ump_x2[l] = v[0];		   		  
		ump_x2[l] += v[1] * *left++;
		ump_x2[l] += v[2] * *left++;
		ump_x2[l] += v[3] * *left++;
		ump_x2[l] += v[4] * *left++;
		ump_x2[l] += v[5] * *left++;
		ump_x2[l] += v[6] * *left++;
		ump_x2[l] += v[7] * *left++;
		ump_x2[l] += v[8] * *left++;
		ump_x2[l] += v[9] * *left++;
		ump_x2[l] += v[10] * *left++;
		ump_x2[l] += v[11] * *left++;
		ump_x2[l] += v[12] * *left++;
		ump_x2[l] += v[13] * *left++;
		ump_x2[l] += v[14] * *left++;
		ump_x2[l] += v[15] * *left++;
		ump_x2[l] += v[16] * *left++;
		ump_x2[l] += v[17] * *left++;
		ump_x2[l] += v[18] * *left++;
		ump_x2[l] += v[19] * *left++;		       		     
	      }	      
	    
	    v = &x3[20 * i];
	    EV = extEV;

	    x1px2 = ump_x1[0] * ump_x2[0];
	    v[0] = x1px2 *  *EV++;	
	    v[1] = x1px2 *  *EV++;
	    v[2] = x1px2 *  *EV++;	
	    v[3] = x1px2 *  *EV++;
	    v[4] = x1px2 *  *EV++;	
	    v[5] = x1px2 *  *EV++;
	    v[6] = x1px2 *  *EV++;	
	    v[7] = x1px2 *  *EV++;
	    v[8] = x1px2 *  *EV++;	
	    v[9] = x1px2 *  *EV++;
	    v[10] = x1px2 *  *EV++;	
	    v[11] = x1px2 *  *EV++;
	    v[12] = x1px2 *  *EV++;	
	    v[13] = x1px2 *  *EV++;
	    v[14] = x1px2 *  *EV++;	
	    v[15] = x1px2 *  *EV++;
	    v[16] = x1px2 *  *EV++;	
	    v[17] = x1px2 *  *EV++;
	    v[18] = x1px2 *  *EV++;	
	    v[19] = x1px2 *  *EV++;
	    
	    for(l = 1; l < 20; l++)
	      {
		x1px2 = ump_x1[l] * ump_x2[l];
		
		v[0] += x1px2 *  *EV++;	
		v[1] += x1px2 *  *EV++;
		v[2] += x1px2 *  *EV++;	
		v[3] += x1px2 *  *EV++;
		v[4] += x1px2 *  *EV++;	
		v[5] += x1px2 *  *EV++;
		v[6] += x1px2 *  *EV++;	
		v[7] += x1px2 *  *EV++;
		v[8] += x1px2 *  *EV++;	
		v[9] += x1px2 *  *EV++;
		v[10] += x1px2 *  *EV++;	
		v[11] += x1px2 *  *EV++;
		v[12] += x1px2 *  *EV++;	
		v[13] += x1px2 *  *EV++;
		v[14] += x1px2 *  *EV++;	
		v[15] += x1px2 *  *EV++;
		v[16] += x1px2 *  *EV++;	
		v[17] += x1px2 *  *EV++;
		v[18] += x1px2 *  *EV++;	
		v[19] += x1px2 *  *EV++;			 		     
	      }
	    
	    /* NO SCALING AT TIPS */
	    	    
	    ex3[i] = 0;		  	    	   	    	 
	  }       
      }
      break;
    case TIP_INNER:    
	 {		   

	   for (i = lower; i < n; i++) 
	     {			     	     	     
	       v = &(tipVector[20 * tipX1[i]]);     
	       left = &left_start[cptr[i] * 380];	    	       	   
	     	       
	       for(l = 0; l < 20; l++)
		 {
		   ump_x1[l] = v[0];		 		  		     
		   ump_x1[l] += v[1] * *left++;
		   ump_x1[l] += v[2] * *left++;
		   ump_x1[l] += v[3] * *left++;
		   ump_x1[l] += v[4] * *left++;
		   ump_x1[l] += v[5] * *left++;
		   ump_x1[l] += v[6] * *left++;
		   ump_x1[l] += v[7] * *left++;
		   ump_x1[l] += v[8] * *left++;
		   ump_x1[l] += v[9] * *left++;
		   ump_x1[l] += v[10] * *left++;
		   ump_x1[l] += v[11] * *left++;
		   ump_x1[l] += v[12] * *left++;
		   ump_x1[l] += v[13] * *left++;
		   ump_x1[l] += v[14] * *left++;
		   ump_x1[l] += v[15] * *left++;
		   ump_x1[l] += v[16] * *left++;
		   ump_x1[l] += v[17] * *left++;
		   ump_x1[l] += v[18] * *left++;
		   ump_x1[l] += v[19] * *left++;		       		     
		 }
	       
	       v = &x2[20 * i];
	       left = &right_start[cptr[i] * 380];

	       for(l = 0; l < 20; l++)
		 {
		   ump_x2[l] = v[0];		   		  
		   ump_x2[l] += v[1] * *left++;
		   ump_x2[l] += v[2] * *left++;
		   ump_x2[l] += v[3] * *left++;
		   ump_x2[l] += v[4] * *left++;
		   ump_x2[l] += v[5] * *left++;
		   ump_x2[l] += v[6] * *left++;
		   ump_x2[l] += v[7] * *left++;
		   ump_x2[l] += v[8] * *left++;
		   ump_x2[l] += v[9] * *left++;
		   ump_x2[l] += v[10] * *left++;
		   ump_x2[l] += v[11] * *left++;
		   ump_x2[l] += v[12] * *left++;
		   ump_x2[l] += v[13] * *left++;
		   ump_x2[l] += v[14] * *left++;
		   ump_x2[l] += v[15] * *left++;
		   ump_x2[l] += v[16] * *left++;
		   ump_x2[l] += v[17] * *left++;
		   ump_x2[l] += v[18] * *left++;
		   ump_x2[l] += v[19] * *left++;		       		     
		 }
	      
	       v = &x3[20 * i];
	       EV = extEV;

	       x1px2 = ump_x1[0] * ump_x2[0];
	       v[0] = x1px2 *  *EV++;	
	       v[1] = x1px2 *  *EV++;
	       v[2] = x1px2 *  *EV++;	
	       v[3] = x1px2 *  *EV++;
	       v[4] = x1px2 *  *EV++;	
	       v[5] = x1px2 *  *EV++;
	       v[6] = x1px2 *  *EV++;	
	       v[7] = x1px2 *  *EV++;
	       v[8] = x1px2 *  *EV++;	
	       v[9] = x1px2 *  *EV++;
	       v[10] = x1px2 *  *EV++;	
	       v[11] = x1px2 *  *EV++;
	       v[12] = x1px2 *  *EV++;	
	       v[13] = x1px2 *  *EV++;
	       v[14] = x1px2 *  *EV++;	
	       v[15] = x1px2 *  *EV++;
	       v[16] = x1px2 *  *EV++;	
	       v[17] = x1px2 *  *EV++;
	       v[18] = x1px2 *  *EV++;	
	       v[19] = x1px2 *  *EV++;
	     
	       for(l = 1; l < 20; l++)
		 {
		   x1px2 = ump_x1[l] * ump_x2[l];
		   
		   v[0] += x1px2 *  *EV++;	
		   v[1] += x1px2 *  *EV++;
		   v[2] += x1px2 *  *EV++;	
		   v[3] += x1px2 *  *EV++;
		   v[4] += x1px2 *  *EV++;	
		   v[5] += x1px2 *  *EV++;
		   v[6] += x1px2 *  *EV++;	
		   v[7] += x1px2 *  *EV++;
		   v[8] += x1px2 *  *EV++;	
		   v[9] += x1px2 *  *EV++;
		   v[10] += x1px2 *  *EV++;	
		   v[11] += x1px2 *  *EV++;
		   v[12] += x1px2 *  *EV++;	
		   v[13] += x1px2 *  *EV++;
		   v[14] += x1px2 *  *EV++;	
		   v[15] += x1px2 *  *EV++;
		   v[16] += x1px2 *  *EV++;	
		   v[17] += x1px2 *  *EV++;
		   v[18] += x1px2 *  *EV++;	
		   v[19] += x1px2 *  *EV++;			 		     
		 }
		   		   	       
	       scale = 1;
	       for(l = 0; scale && (l < 20); l++)
		 scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));
	       
	       ex3[i] = ex2[i];		  	      	       	       
	       
	       if(scale)
		 {			  
		   for(l = 0; l < 20; l++)
		     v[l] *= twotothe256;		   
		   ex3[i] += 1;
		 }	            	       
	     }     	 	  
	 }
	 break;
    case INNER_INNER:

      for (i = lower; i < n; i++) 
	{		 	  	 
	  v = &x1[20 * i];   	   
	  left = &left_start[cptr[i] * 380];
	   
	  for(l = 0; l < 20; l++)
	    {
	      ump_x1[l] = v[0];		 		  		     
	      ump_x1[l] += v[1] * *left++;
	      ump_x1[l] += v[2] * *left++;
	      ump_x1[l] += v[3] * *left++;
	      ump_x1[l] += v[4] * *left++;
	      ump_x1[l] += v[5] * *left++;
	      ump_x1[l] += v[6] * *left++;
	      ump_x1[l] += v[7] * *left++;
	      ump_x1[l] += v[8] * *left++;
	      ump_x1[l] += v[9] * *left++;
	      ump_x1[l] += v[10] * *left++;
	      ump_x1[l] += v[11] * *left++;
	      ump_x1[l] += v[12] * *left++;
	      ump_x1[l] += v[13] * *left++;
	      ump_x1[l] += v[14] * *left++;
	      ump_x1[l] += v[15] * *left++;
	      ump_x1[l] += v[16] * *left++;
	      ump_x1[l] += v[17] * *left++;
	      ump_x1[l] += v[18] * *left++;
	      ump_x1[l] += v[19] * *left++;		       		     
	    }
	  
	  v = &x2[20 * i];	   	    
	  left = &right_start[cptr[i] * 380];
	  
	  for(l = 0; l < 20; l++)
	    {
	      ump_x2[l] = v[0];		   		  
	      ump_x2[l] += v[1] * *left++;
	      ump_x2[l] += v[2] * *left++;
	      ump_x2[l] += v[3] * *left++;
	      ump_x2[l] += v[4] * *left++;
	      ump_x2[l] += v[5] * *left++;
	      ump_x2[l] += v[6] * *left++;
	      ump_x2[l] += v[7] * *left++;
	      ump_x2[l] += v[8] * *left++;
	      ump_x2[l] += v[9] * *left++;
	      ump_x2[l] += v[10] * *left++;
	      ump_x2[l] += v[11] * *left++;
	      ump_x2[l] += v[12] * *left++;
	      ump_x2[l] += v[13] * *left++;
	      ump_x2[l] += v[14] * *left++;
	      ump_x2[l] += v[15] * *left++;
	      ump_x2[l] += v[16] * *left++;
	      ump_x2[l] += v[17] * *left++;
	      ump_x2[l] += v[18] * *left++;
	      ump_x2[l] += v[19] * *left++;		       		     
	    }
	  	  
	   v = &x3[20 * i];
	   EV = extEV;

	   x1px2 = ump_x1[0] * ump_x2[0];
	   v[0] = x1px2 *  *EV++;	
	   v[1] = x1px2 *  *EV++;
	   v[2] = x1px2 *  *EV++;	
	   v[3] = x1px2 *  *EV++;
	   v[4] = x1px2 *  *EV++;	
	   v[5] = x1px2 *  *EV++;
	   v[6] = x1px2 *  *EV++;	
	   v[7] = x1px2 *  *EV++;
	   v[8] = x1px2 *  *EV++;	
	   v[9] = x1px2 *  *EV++;
	   v[10] = x1px2 *  *EV++;	
	   v[11] = x1px2 *  *EV++;
	   v[12] = x1px2 *  *EV++;	
	   v[13] = x1px2 *  *EV++;
	   v[14] = x1px2 *  *EV++;	
	   v[15] = x1px2 *  *EV++;
	   v[16] = x1px2 *  *EV++;	
	   v[17] = x1px2 *  *EV++;
	   v[18] = x1px2 *  *EV++;	
	   v[19] = x1px2 *  *EV++;
	   
	   for(l = 1; l < 20; l++)
	     {
	       x1px2 = ump_x1[l] * ump_x2[l];
	       
	       v[0] += x1px2 *  *EV++;	
	       v[1] += x1px2 *  *EV++;
	       v[2] += x1px2 *  *EV++;	
	       v[3] += x1px2 *  *EV++;
	       v[4] += x1px2 *  *EV++;	
	       v[5] += x1px2 *  *EV++;
	       v[6] += x1px2 *  *EV++;	
	       v[7] += x1px2 *  *EV++;
	       v[8] += x1px2 *  *EV++;	
	       v[9] += x1px2 *  *EV++;
	       v[10] += x1px2 *  *EV++;	
	       v[11] += x1px2 *  *EV++;
	       v[12] += x1px2 *  *EV++;	
	       v[13] += x1px2 *  *EV++;
	       v[14] += x1px2 *  *EV++;	
	       v[15] += x1px2 *  *EV++;
	       v[16] += x1px2 *  *EV++;	
	       v[17] += x1px2 *  *EV++;
	       v[18] += x1px2 *  *EV++;	
	       v[19] += x1px2 *  *EV++;			 		     
	     }
	   
	   scale = 1;
	   for(l = 0; scale && (l < 20); l++)
	     scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));

	   ex3[i] = ex1[i] + ex2[i];		  
	   
	   if(scale)	
	     {	 	       
	       for(l = 0; l < 20; l++)
		 v[l] *= twotothe256;
	       ex3[i] += 1;
	     }	            
	} 
      break;
    default:
      assert(0);
    }
  free(left_start); 	
  free(right_start);  
}


static void newviewGTRCATPROTMULT(traversalInfo *ti, double *extEV,  double *extEI,  double *EIGN,
				     double *rptr,  int *cptr,
				     double *x1, double *x2, double *x3, double *tipVector,				
				     int    *ex1, int *ex2, int *ex3, char *tipX1, char *tipX2, int *modelptr,
				     int lower, int n,  int numberOfCategories, int numberOfModels, int multiBranch)
{   
  double  *left, *left_start, *right, *right_start,  *EI, *EV, *v;
  double      
    ump_x1[20], ump_x2[20], lz1[19], lz2[19],d1[19], d2[19]; 
  double 
    x1px2, z1 = 0.0, z2 = 0.0, ki;
  int i, l, scale, modelCounter, model;              

  if(!multiBranch)
    {
      z1  = ti->qz[0];	  	   
      z2  = ti->rz[0];	 
    }
  
  left_start  = left  = (double *)malloc(380 * numberOfModels * numberOfCategories * sizeof(double));
  right_start = right = (double *)malloc(380 * numberOfModels * numberOfCategories * sizeof(double));
                    
  for(modelCounter = 0; modelCounter < numberOfModels; modelCounter++)
    {       
      if(multiBranch)
	{
	  z1  = ti->qz[modelCounter];	      	       
	  z2  = ti->rz[modelCounter];	        
	}
      
      for(i = 0; i < 19; i++)
	lz1[i] = EIGN[modelCounter * 19 + i] * z1;

      for(i = 0; i < 19; i++)
	lz2[i] = EIGN[modelCounter * 19 + i] * z2;                    
     
      for(i = 0; i < numberOfCategories; i++)
	{	
	  ki = rptr[i];
	   	   	   
	  d1[0] = exp (ki * lz1[0]);
	  d1[1] = exp (ki * lz1[1]);
	  d1[2] = exp (ki * lz1[2]);
	  d1[3] = exp (ki * lz1[3]);
	  d1[4] = exp (ki * lz1[4]);
	  d1[5] = exp (ki * lz1[5]);
	  d1[6] = exp (ki * lz1[6]);
	  d1[7] = exp (ki * lz1[7]);
	  d1[8] = exp (ki * lz1[8]);
	  d1[9] = exp (ki * lz1[9]);
	  d1[10] = exp (ki * lz1[10]);
	  d1[11] = exp (ki * lz1[11]);
	  d1[12] = exp (ki * lz1[12]);
	  d1[13] = exp (ki * lz1[13]);
	  d1[14] = exp (ki * lz1[14]);
	  d1[15] = exp (ki * lz1[15]);
	  d1[16] = exp (ki * lz1[16]);
	  d1[17] = exp (ki * lz1[17]);
	  d1[18] = exp (ki * lz1[18]);
	       
	  d2[0] = exp (ki * lz2[0]);
	  d2[1] = exp (ki * lz2[1]);
	  d2[2] = exp (ki * lz2[2]);
	  d2[3] = exp (ki * lz2[3]);
	  d2[4] = exp (ki * lz2[4]);
	  d2[5] = exp (ki * lz2[5]);
	  d2[6] = exp (ki * lz2[6]);
	  d2[7] = exp (ki * lz2[7]);
	  d2[8] = exp (ki * lz2[8]);
	  d2[9] = exp (ki * lz2[9]);
	  d2[10] = exp (ki * lz2[10]);
	  d2[11] = exp (ki * lz2[11]);
	  d2[12] = exp (ki * lz2[12]);
	  d2[13] = exp (ki * lz2[13]);
	  d2[14] = exp (ki * lz2[14]);
	  d2[15] = exp (ki * lz2[15]);
	  d2[16] = exp (ki * lz2[16]);
	  d2[17] = exp (ki * lz2[17]);
	  d2[18] = exp (ki * lz2[18]);	   	      	 
	  
	  EI = &(extEI[modelCounter * 380]);

	  for(l = 0; l < 20; l++)
	    {	       
	      *left++ = d1[0] * *EI++;	     	     		   	      
	      *left++ = d1[1] * *EI++;	 
	      *left++ = d1[2] * *EI++;	     	     		   	      
	      *left++ = d1[3] * *EI++;	 
	      *left++ = d1[4] * *EI++;	     	     		   	      
	      *left++ = d1[5] * *EI++;	 
	      *left++ = d1[6] * *EI++;	     	     		   	      
	      *left++ = d1[7] * *EI++;
	      *left++ = d1[8] * *EI++;	     	     		   	      
	      *left++ = d1[9] * *EI++;	 
	      *left++ = d1[10] * *EI++;	     	     		   	      
	      *left++ = d1[11] * *EI++;	 
	      *left++ = d1[12] * *EI++;	     	     		   	      
	      *left++ = d1[13] * *EI++;	 
	      *left++ = d1[14] * *EI++;	     	     		   	      
	      *left++ = d1[15] * *EI++;
	      *left++ = d1[16] * *EI++;	     	     		   	      
	      *left++ = d1[17] * *EI++;	 
	      *left++ = d1[18] * *EI++;	     	     		   	      	      
	    }

	  EI = &(extEI[modelCounter * 380]);
	   
	  for(l = 0; l < 20; l++)
	    {
	      *right++ = d2[0] * *EI++;	     	     		   	      
	      *right++ = d2[1] * *EI++;	 
	      *right++ = d2[2] * *EI++;	     	     		   	      
	      *right++ = d2[3] * *EI++;	 
	      *right++ = d2[4] * *EI++;	     	     		   	      
	      *right++ = d2[5] * *EI++;	 
	      *right++ = d2[6] * *EI++;	     	     		   	      
	      *right++ = d2[7] * *EI++;
	      *right++ = d2[8] * *EI++;	     	     		   	      
	      *right++ = d2[9] * *EI++;	 
	      *right++ = d2[10] * *EI++;	     	     		   	      
	      *right++ = d2[11] * *EI++;	 
	      *right++ = d2[12] * *EI++;	     	     		   	      
	      *right++ = d2[13] * *EI++;	 
	      *right++ = d2[14] * *EI++;	     	     		   	      
	      *right++ = d2[15] * *EI++;
	      *right++ = d2[16] * *EI++;	     	     		   	      
	      *right++ = d2[17] * *EI++;	 
	      *right++ = d2[18] * *EI++;    	     		   	      
	    }	  
	}
    }
  
  switch(ti->tipCase)
    {
    case TIP_TIP:
      {	  	 	 	   	   	        	 
	   
	for (i = lower; i < n; i++) 
	  {	 
	    model = modelptr[i];	   	    	   	   	     
	    EV = &(extEV[model * 400]);
	    left = &left_start[model * 380 * numberOfCategories + cptr[i] * 380];	      	       
	    v = &(tipVector[460 * model + 20 * tipX1[i]]);	     	       	       

	    for(l = 0; l < 20; l++)
	      {
		ump_x1[l] = v[0];		 		  		     
		ump_x1[l] += v[1] * *left++;
		ump_x1[l] += v[2] * *left++;
		ump_x1[l] += v[3] * *left++;
		ump_x1[l] += v[4] * *left++;
		ump_x1[l] += v[5] * *left++;
		ump_x1[l] += v[6] * *left++;
		ump_x1[l] += v[7] * *left++;
		ump_x1[l] += v[8] * *left++;
		ump_x1[l] += v[9] * *left++;
		ump_x1[l] += v[10] * *left++;
		ump_x1[l] += v[11] * *left++;
		ump_x1[l] += v[12] * *left++;
		ump_x1[l] += v[13] * *left++;
		ump_x1[l] += v[14] * *left++;
		ump_x1[l] += v[15] * *left++;
		ump_x1[l] += v[16] * *left++;
		ump_x1[l] += v[17] * *left++;
		ump_x1[l] += v[18] * *left++;
		ump_x1[l] += v[19] * *left++;		       		     
	      }
	    
	    v = &(tipVector[460 * model + 20 * tipX2[i]]);

	    left = &right_start[model * 380 * numberOfCategories + cptr[i] * 380];

	    for(l = 0; l < 20; l++)
	      {
		ump_x2[l] = v[0];		   		  
		ump_x2[l] += v[1] * *left++;
		ump_x2[l] += v[2] * *left++;
		ump_x2[l] += v[3] * *left++;
		ump_x2[l] += v[4] * *left++;
		ump_x2[l] += v[5] * *left++;
		ump_x2[l] += v[6] * *left++;
		ump_x2[l] += v[7] * *left++;
		ump_x2[l] += v[8] * *left++;
		ump_x2[l] += v[9] * *left++;
		ump_x2[l] += v[10] * *left++;
		ump_x2[l] += v[11] * *left++;
		ump_x2[l] += v[12] * *left++;
		ump_x2[l] += v[13] * *left++;
		ump_x2[l] += v[14] * *left++;
		ump_x2[l] += v[15] * *left++;
		ump_x2[l] += v[16] * *left++;
		ump_x2[l] += v[17] * *left++;
		ump_x2[l] += v[18] * *left++;
		ump_x2[l] += v[19] * *left++;		       		     
	      }	      
	    
	    v = &x3[20 * i];
	    
	    x1px2 = ump_x1[0] * ump_x2[0];
	    v[0] = x1px2 *  *EV++;	
	    v[1] = x1px2 *  *EV++;
	    v[2] = x1px2 *  *EV++;	
	    v[3] = x1px2 *  *EV++;
	    v[4] = x1px2 *  *EV++;	
	    v[5] = x1px2 *  *EV++;
	    v[6] = x1px2 *  *EV++;	
	    v[7] = x1px2 *  *EV++;
	    v[8] = x1px2 *  *EV++;	
	    v[9] = x1px2 *  *EV++;
	    v[10] = x1px2 *  *EV++;	
	    v[11] = x1px2 *  *EV++;
	    v[12] = x1px2 *  *EV++;	
	    v[13] = x1px2 *  *EV++;
	    v[14] = x1px2 *  *EV++;	
	    v[15] = x1px2 *  *EV++;
	    v[16] = x1px2 *  *EV++;	
	    v[17] = x1px2 *  *EV++;
	    v[18] = x1px2 *  *EV++;	
	    v[19] = x1px2 *  *EV++;
	    	    
	    for(l = 1; l < 20; l++)
	      {
		x1px2 = ump_x1[l] * ump_x2[l];
		
		v[0] += x1px2 *  *EV++;	
		v[1] += x1px2 *  *EV++;
		v[2] += x1px2 *  *EV++;	
		v[3] += x1px2 *  *EV++;
		v[4] += x1px2 *  *EV++;	
		v[5] += x1px2 *  *EV++;
		v[6] += x1px2 *  *EV++;	
		v[7] += x1px2 *  *EV++;
		v[8] += x1px2 *  *EV++;	
		v[9] += x1px2 *  *EV++;
		v[10] += x1px2 *  *EV++;	
		v[11] += x1px2 *  *EV++;
		v[12] += x1px2 *  *EV++;	
		v[13] += x1px2 *  *EV++;
		v[14] += x1px2 *  *EV++;	
		v[15] += x1px2 *  *EV++;
		v[16] += x1px2 *  *EV++;	
		v[17] += x1px2 *  *EV++;
		v[18] += x1px2 *  *EV++;	
		v[19] += x1px2 *  *EV++;			 		     
	      }
		   
	    /* NO SCALING AT TIPS */
	      
	    ex3[i] = 0;		  	    	   	    	 
	  }
      }
      break;	      
    case TIP_INNER:
      {	

	for (i = lower; i < n; i++) 
	  {			     	     	    
	    model = modelptr[i];
	    EV = &(extEV[model * 400]);	      
	    v = &(tipVector[model * 460 + 20 * tipX1[i]]);     
	    left = &left_start[model * 380 * numberOfCategories + cptr[i] * 380];	    	       	   
	      
	    for(l = 0; l < 20; l++)
	      {
		ump_x1[l] = v[0];		 		  		     
		ump_x1[l] += v[1] * *left++;
		ump_x1[l] += v[2] * *left++;
		ump_x1[l] += v[3] * *left++;
		ump_x1[l] += v[4] * *left++;
		ump_x1[l] += v[5] * *left++;
		ump_x1[l] += v[6] * *left++;
		ump_x1[l] += v[7] * *left++;
		ump_x1[l] += v[8] * *left++;
		ump_x1[l] += v[9] * *left++;
		ump_x1[l] += v[10] * *left++;
		ump_x1[l] += v[11] * *left++;
		ump_x1[l] += v[12] * *left++;
		ump_x1[l] += v[13] * *left++;
		ump_x1[l] += v[14] * *left++;
		ump_x1[l] += v[15] * *left++;
		ump_x1[l] += v[16] * *left++;
		ump_x1[l] += v[17] * *left++;
		ump_x1[l] += v[18] * *left++;
		ump_x1[l] += v[19] * *left++;		       		     
	      }
	    
	    v = &x2[20 * i];
	    
	    left = &right_start[model * 380 * numberOfCategories + cptr[i] * 380];	     
	    
	    for(l = 0; l < 20; l++)
	      {
		ump_x2[l] = v[0];		   		  
		ump_x2[l] += v[1] * *left++;
		ump_x2[l] += v[2] * *left++;
		ump_x2[l] += v[3] * *left++;
		ump_x2[l] += v[4] * *left++;
		ump_x2[l] += v[5] * *left++;
		ump_x2[l] += v[6] * *left++;
		ump_x2[l] += v[7] * *left++;
		ump_x2[l] += v[8] * *left++;
		ump_x2[l] += v[9] * *left++;
		ump_x2[l] += v[10] * *left++;
		ump_x2[l] += v[11] * *left++;
		ump_x2[l] += v[12] * *left++;
		ump_x2[l] += v[13] * *left++;
		ump_x2[l] += v[14] * *left++;
		ump_x2[l] += v[15] * *left++;
		ump_x2[l] += v[16] * *left++;
		ump_x2[l] += v[17] * *left++;
		ump_x2[l] += v[18] * *left++;
		ump_x2[l] += v[19] * *left++;		       		     
	      }
	    
	    v = &x3[20 * i];	    
	    
	    x1px2 = ump_x1[0] * ump_x2[0];    
	    v[0] = x1px2 *  *EV++;		    
	    v[1] = x1px2 *  *EV++;
	    v[2] = x1px2 *  *EV++;	
	    v[3] = x1px2 *  *EV++;
	    v[4] = x1px2 *  *EV++;	
	    v[5] = x1px2 *  *EV++;
	    v[6] = x1px2 *  *EV++;	
	    v[7] = x1px2 *  *EV++;
	    v[8] = x1px2 *  *EV++;	
	    v[9] = x1px2 *  *EV++;
	    v[10] = x1px2 *  *EV++;	
	    v[11] = x1px2 *  *EV++;
	    v[12] = x1px2 *  *EV++;	
	    v[13] = x1px2 *  *EV++;
	    v[14] = x1px2 *  *EV++;	
	    v[15] = x1px2 *  *EV++;
	    v[16] = x1px2 *  *EV++;	
	    v[17] = x1px2 *  *EV++;
	    v[18] = x1px2 *  *EV++;	
	    v[19] = x1px2 *  *EV++;	     	 
	    
	    for(l = 1; l < 20; l++)
	      {
		x1px2 = ump_x1[l] * ump_x2[l];
		
		v[0] += x1px2 *  *EV++;	
		v[1] += x1px2 *  *EV++;
		v[2] += x1px2 *  *EV++;	
		v[3] += x1px2 *  *EV++;
		v[4] += x1px2 *  *EV++;	
		v[5] += x1px2 *  *EV++;
		v[6] += x1px2 *  *EV++;	
		v[7] += x1px2 *  *EV++;
		v[8] += x1px2 *  *EV++;	
		v[9] += x1px2 *  *EV++;
		v[10] += x1px2 *  *EV++;	
		v[11] += x1px2 *  *EV++;
		v[12] += x1px2 *  *EV++;	
		v[13] += x1px2 *  *EV++;
		v[14] += x1px2 *  *EV++;	
		v[15] += x1px2 *  *EV++;
		v[16] += x1px2 *  *EV++;	
		v[17] += x1px2 *  *EV++;
		v[18] += x1px2 *  *EV++;	
		v[19] += x1px2 *  *EV++;			 		     
	      }
	    	    
	    scale = 1;
	    for(l = 0; scale && (l < 20); l++)
	      scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));	       	   

	    ex3[i] = ex2[i];		  	      	       	       
	    
	    if(scale)
	      {			  
		for(l = 0; l < 20; l++)
		  v[l] *= twotothe256;		   
		ex3[i] += 1;
	      }	            	       
	  }     	 	
      }
      break;
    case INNER_INNER:
                
      for (i = lower; i < n; i++) 
	{		 	  	 
	  model = modelptr[i];
	  v = &x1[20 * i];   	   
	  left = &left_start[model * 380 * numberOfCategories + cptr[i] * 380];
	   
	  for(l = 0; l < 20; l++)
	    {
	      ump_x1[l] = v[0];		 		  		     
	      ump_x1[l] += v[1] * *left++;
	      ump_x1[l] += v[2] * *left++;
	      ump_x1[l] += v[3] * *left++;
	      ump_x1[l] += v[4] * *left++;
	      ump_x1[l] += v[5] * *left++;
	      ump_x1[l] += v[6] * *left++;
	      ump_x1[l] += v[7] * *left++;
	      ump_x1[l] += v[8] * *left++;
	      ump_x1[l] += v[9] * *left++;
	      ump_x1[l] += v[10] * *left++;
	      ump_x1[l] += v[11] * *left++;
	      ump_x1[l] += v[12] * *left++;
	      ump_x1[l] += v[13] * *left++;
	      ump_x1[l] += v[14] * *left++;
	      ump_x1[l] += v[15] * *left++;
	      ump_x1[l] += v[16] * *left++;
	      ump_x1[l] += v[17] * *left++;
	      ump_x1[l] += v[18] * *left++;
	      ump_x1[l] += v[19] * *left++;		       		     
	    }
	  
	  v = &x2[20 * i];	   	    
	  left = &right_start[model * 380 * numberOfCategories + cptr[i] * 380];

	  for(l = 0; l < 20; l++)
	    {
	      ump_x2[l] = v[0];		   		  
	      ump_x2[l] += v[1] * *left++;
	      ump_x2[l] += v[2] * *left++;
	      ump_x2[l] += v[3] * *left++;
	      ump_x2[l] += v[4] * *left++;
	      ump_x2[l] += v[5] * *left++;
	      ump_x2[l] += v[6] * *left++;
	      ump_x2[l] += v[7] * *left++;
	      ump_x2[l] += v[8] * *left++;
	      ump_x2[l] += v[9] * *left++;
	      ump_x2[l] += v[10] * *left++;
	      ump_x2[l] += v[11] * *left++;
	      ump_x2[l] += v[12] * *left++;
	      ump_x2[l] += v[13] * *left++;
	      ump_x2[l] += v[14] * *left++;
	      ump_x2[l] += v[15] * *left++;
	      ump_x2[l] += v[16] * *left++;
	      ump_x2[l] += v[17] * *left++;
	      ump_x2[l] += v[18] * *left++;
	      ump_x2[l] += v[19] * *left++;		       		     
	    }
	   
	  EV = &(extEV[model * 400]);
	  v = &x3[20 * i];
	   
	  x1px2 = ump_x1[0] * ump_x2[0];
	  v[0] = x1px2 *  *EV++;	
	  v[1] = x1px2 *  *EV++;
	  v[2] = x1px2 *  *EV++;	
	  v[3] = x1px2 *  *EV++;
	  v[4] = x1px2 *  *EV++;	
	  v[5] = x1px2 *  *EV++;
	  v[6] = x1px2 *  *EV++;	
	  v[7] = x1px2 *  *EV++;
	  v[8] = x1px2 *  *EV++;	
	  v[9] = x1px2 *  *EV++;
	  v[10] = x1px2 *  *EV++;	
	  v[11] = x1px2 *  *EV++;
	  v[12] = x1px2 *  *EV++;	
	  v[13] = x1px2 *  *EV++;
	  v[14] = x1px2 *  *EV++;	
	  v[15] = x1px2 *  *EV++;
	  v[16] = x1px2 *  *EV++;	
	  v[17] = x1px2 *  *EV++;
	  v[18] = x1px2 *  *EV++;	
	  v[19] = x1px2 *  *EV++;
	  	  
	  for(l = 1; l < 20; l++)
	    {
	      x1px2 = ump_x1[l] * ump_x2[l];
	      
	      v[0] += x1px2 *  *EV++;	
	      v[1] += x1px2 *  *EV++;
	      v[2] += x1px2 *  *EV++;	
	      v[3] += x1px2 *  *EV++;
	      v[4] += x1px2 *  *EV++;	
	      v[5] += x1px2 *  *EV++;
	      v[6] += x1px2 *  *EV++;	
	      v[7] += x1px2 *  *EV++;
	      v[8] += x1px2 *  *EV++;	
	      v[9] += x1px2 *  *EV++;
	      v[10] += x1px2 *  *EV++;	
	      v[11] += x1px2 *  *EV++;
	      v[12] += x1px2 *  *EV++;	
	      v[13] += x1px2 *  *EV++;
	      v[14] += x1px2 *  *EV++;	
	      v[15] += x1px2 *  *EV++;
	      v[16] += x1px2 *  *EV++;	
	      v[17] += x1px2 *  *EV++;
	      v[18] += x1px2 *  *EV++;	
	      v[19] += x1px2 *  *EV++;			 		     
	    }	 
	  
	  scale = 1;
	  for(l = 0; scale && (l < 20); l++)
	    scale = ((v[l] < minlikelihood) && (v[l] > minusminlikelihood));
	  
	  ex3[i] = ex1[i] + ex2[i];		  
	  
	  if(scale)	
	    {	 	       
	      for(l = 0; l < 20; l++)
		v[l] *= twotothe256;
	      ex3[i] += 1;
	    }	            
	} 
      break;
    default:
      assert(0);
    }
  free(left_start); 	
  free(right_start);  
}


static void newviewGTRGAMMAPROT(traversalInfo *ti,
				   double *x1, double *x2, double *x3,
				   double *EIGN, double *extEV, double *extEI, double *gammaRates, double *tipVector,
				   int    *ex1, int *ex2, int *ex3, char *tipX1, char *tipX2,
				   int lower, int n, double z1, double z2)
{  
  double  *left, *right, *left_start, *right_start, *uX1, *uX2, *v, *EV,*EI;
  double x1px2;   
  int  i, l, k, scale;   
  double  ki, ump_x1_0;       
  double  d1[19], d2[19];  
  double *vl, *vr, al, ar;         
      
  left_start  = left  = (double *)malloc(1520 * sizeof(double));
  right_start = right = (double *)malloc(1520 * sizeof(double));
      
  for(i = 0; i < 4;i++)
    {       
      ki = gammaRates[i];	    
	  
      for(k = 0; k < 19; k++)
	{
	  d1[k] = exp (EIGN[k] * ki * z1);
	  d2[k] = exp (EIGN[k] * ki * z2);
	}
	
      EI = extEI;
      for(k = 0; k < 20; k++)
	{
	  *left++ = d1[0] * *EI++;
	  *left++ = d1[1] * *EI++;
	  *left++ = d1[2] * *EI++;
	  *left++ = d1[3] * *EI++;
	  *left++ = d1[4] * *EI++;
	  *left++ = d1[5] * *EI++;
	  *left++ = d1[6] * *EI++;
	  *left++ = d1[7] * *EI++;
	  *left++ = d1[8] * *EI++;
	  *left++ = d1[9] * *EI++;
	  *left++ = d1[10] * *EI++;
	  *left++ = d1[11] * *EI++;
	  *left++ = d1[12] * *EI++;
	  *left++ = d1[13] * *EI++;
	  *left++ = d1[14] * *EI++;
	  *left++ = d1[15] * *EI++;
	  *left++ = d1[16] * *EI++;
	  *left++ = d1[17] * *EI++;
	  *left++ = d1[18] * *EI++;	      
	}
            
      EI = extEI;
      
      for(k = 0; k < 20; k++)
	{
	  *right++ = d2[0] * *EI++;
	  *right++ = d2[1] * *EI++;
	  *right++ = d2[2] * *EI++;
	  *right++ = d2[3] * *EI++;
	  *right++ = d2[4] * *EI++;
	  *right++ = d2[5] * *EI++;
	  *right++ = d2[6] * *EI++;
	  *right++ = d2[7] * *EI++;
	  *right++ = d2[8] * *EI++;
	  *right++ = d2[9] * *EI++;
	  *right++ = d2[10] * *EI++;
	  *right++ = d2[11] * *EI++;
	  *right++ = d2[12] * *EI++;
	  *right++ = d2[13] * *EI++;
	  *right++ = d2[14] * *EI++;
	  *right++ = d2[15] * *EI++;
	  *right++ = d2[16] * *EI++;
	  *right++ = d2[17] * *EI++;
	  *right++ = d2[18] * *EI++;	      
	}	
    }          
  
  switch(ti->tipCase)
    {      
    case TIP_TIP:
      {      
	double umpX1[1840], umpX2[1840];

	uX1 = umpX1;
	uX2 = umpX2;
	
	for(i = 0; i < 23; i++)
	  {	    	     
	    v = &(tipVector[20 * i]);	     
	    left = left_start;
	    right = right_start;
	    
	    for(k = 0; k < 80; k++)
	      {
		ump_x1_0 = v[0];
		for(l = 1; l < 20; l++)
		  ump_x1_0 +=  v[l] * *left++;
		*uX1++ = ump_x1_0;
	      }
	    
	    for(k = 0; k < 80; k++)
	      {
		ump_x1_0 = v[0];
		for(l = 1; l < 20; l++)
		  ump_x1_0 +=  v[l] * *right++;
		*uX2++ = ump_x1_0;
	      }
	  }	
	

	for (i = lower; i < n; i++) 
	  {		     
	    uX1 = &umpX1[80 * tipX1[i]];
	    uX2 = &umpX2[80 * tipX2[i]];
	     
	    for(k = 0; k < 4; k++)
	      {
		EV = extEV;
		x1px2 = *uX1++ * *uX2++;
		v = &(x3[80 * i + k * 20]);
		for(l = 0; l < 20; l++)
		  v[l] = x1px2 *  *EV++;
		for(l = 0; l < 19; l++)
		  {
		    x1px2 = *uX1++ * *uX2++;
		    v[0] += x1px2 *   *EV++;
		    v[1] += x1px2 *   *EV++;
		    v[2] += x1px2 *   *EV++;
		    v[3] += x1px2 *   *EV++;
		    v[4] += x1px2 *   *EV++;
		    v[5] += x1px2 *   *EV++;
		    v[6] += x1px2 *   *EV++;
		    v[7] += x1px2 *   *EV++;
		    v[8] += x1px2 *   *EV++;
		    v[9] += x1px2 *   *EV++;
		    v[10] += x1px2 *  *EV++;
		    v[11] += x1px2 *  *EV++;
		    v[12] += x1px2 *  *EV++;
		    v[13] += x1px2 *  *EV++;
		    v[14] += x1px2 *  *EV++;
		    v[15] += x1px2 *  *EV++;
		    v[16] += x1px2 *  *EV++;
		    v[17] += x1px2 *  *EV++;
		    v[18] += x1px2 *  *EV++;
		    v[19] += x1px2 *  *EV++;
		  }
	      }

	    ex3[i] = 0;		  		    	     	            
	  }         	
      }
      break;
    case TIP_INNER:
      {	    
	double umpX1[1840], ump_x2[20];	
	uX1 = umpX1;
	 
	for(i = 0; i < 23; i++)
	  {	    	     		 
	    v = &(tipVector[20 * i]);	     
	    left = left_start;
		       
	    for(k = 0; k < 80; k++)
	      {		 		 	
		ump_x1_0 = v[0];
		for(l = 1; l < 20; l++)		 
		  ump_x1_0 +=  v[l] * *left++;		 
		*uX1++ = ump_x1_0;
	      }
	  }
	       	    	       	


	for (i = lower; i < n; i++) 
	  {		   
	    right = right_start;
	    uX1 = &umpX1[80 * tipX1[i]];			    	    
	    
	    for(k = 0; k < 4; k++)
	      {
		v = &(x2[80 * i + k * 20]);	    
		for(l = 0; l < 20; l++)
		  {		     
		    ump_x2[l] =  v[0];		 
		    ump_x2[l] += v[1] * *right++;
		    ump_x2[l] += v[2] * *right++;
		    ump_x2[l] += v[3] * *right++;
		    ump_x2[l] += v[4] * *right++;
		    ump_x2[l] += v[5] * *right++;
		    ump_x2[l] += v[6] * *right++;
		    ump_x2[l] += v[7] * *right++;
		    ump_x2[l] += v[8] * *right++;
		    ump_x2[l] += v[9] * *right++;
		    ump_x2[l] += v[10] * *right++;
		    ump_x2[l] += v[11] * *right++;
		    ump_x2[l] += v[12] * *right++;
		    ump_x2[l] += v[13] * *right++;
		    ump_x2[l] += v[14] * *right++;
		    ump_x2[l] += v[15] * *right++;
		    ump_x2[l] += v[16] * *right++;
		    ump_x2[l] += v[17] * *right++;
		    ump_x2[l] += v[18] * *right++;
		    ump_x2[l] += v[19] * *right++;
		  }
		 	     
		uX2 = ump_x2;
		EV = extEV;
		x1px2 = *uX1++ * *uX2++;
		v = &(x3[80 * i + 20 * k]);
		for(l = 0; l < 20; l++)
		  v[l] = x1px2 *  *EV++;
		
		for(l = 0; l < 19; l++)
		  {
		    x1px2 = *uX1++ * *uX2++;
		    v[0] += x1px2 *   *EV++;
		    v[1] += x1px2 *   *EV++;
		    v[2] += x1px2 *   *EV++;
		    v[3] += x1px2 *   *EV++;
		    v[4] += x1px2 *   *EV++;
		    v[5] += x1px2 *   *EV++;
		    v[6] += x1px2 *   *EV++;
		    v[7] += x1px2 *   *EV++;
		    v[8] += x1px2 *   *EV++;
		    v[9] += x1px2 *   *EV++;
		    v[10] += x1px2 *  *EV++;
		    v[11] += x1px2 *  *EV++;
		    v[12] += x1px2 *  *EV++;
		    v[13] += x1px2 *  *EV++;
		    v[14] += x1px2 *  *EV++;
		    v[15] += x1px2 *  *EV++;
		    v[16] += x1px2 *  *EV++;
		    v[17] += x1px2 *  *EV++;
		    v[18] += x1px2 *  *EV++;
		    v[19] += x1px2 *  *EV++;
		  }
	      }	
	    
	    ex3[i] = ex2[i];		
	    v = &x3[80 * i];
	    scale = 1;
	    for(l = 0; scale && (l < 80); l++)
	      scale = (ABS(v[l]) <  minlikelihood);
	    
	    if (scale)	         
	      {	     	    	 
		for(l = 0; l < 80; l++)
		  v[l] *= twotothe256;		   
		ex3[i] += 1;  
	      }	 	           	 
	  }
      }
      break;
    case INNER_INNER:	      

      for (i = lower; i < n; i++) 
       {		 
	 left = left_start;
	 right = right_start;

	 for(k = 0; k < 4; k++)
	   {
	     EV = extEV;
	     vl = &(x1[80 * i + 20 * k]);
	     vr = &(x2[80 * i + 20 * k]);
	     v =  &(x3[80 * i + 20 * k]);
	     	    
	     for(l = 0; l < 20; l++)	       
	       v[l] = 0;		 
	       
	     for(l = 0; l < 20; l++)
	       {
		 al =  vl[0];		 
		 al += vl[1] * *left++;
		 al += vl[2] * *left++;
		 al += vl[3] * *left++;
		 al += vl[4] * *left++;
		 al += vl[5] * *left++;
		 al += vl[6] * *left++;
		 al += vl[7] * *left++;
		 al += vl[8] * *left++;
		 al += vl[9] * *left++;
		 al += vl[10] * *left++;
		 al += vl[11] * *left++;
		 al += vl[12] * *left++;
		 al += vl[13] * *left++;
		 al += vl[14] * *left++;
		 al += vl[15] * *left++;
		 al += vl[16] * *left++;
		 al += vl[17] * *left++;
		 al += vl[18] * *left++;
		 al += vl[19] * *left++;

		 ar =  vr[0];		 
		 ar += vr[1] * *right++;
		 ar += vr[2] * *right++;
		 ar += vr[3] * *right++;
		 ar += vr[4] * *right++;
		 ar += vr[5] * *right++;
		 ar += vr[6] * *right++;
		 ar += vr[7] * *right++;
		 ar += vr[8] * *right++;
		 ar += vr[9] * *right++;
		 ar += vr[10] * *right++;
		 ar += vr[11] * *right++;
		 ar += vr[12] * *right++;
		 ar += vr[13] * *right++;
		 ar += vr[14] * *right++;
		 ar += vr[15] * *right++;
		 ar += vr[16] * *right++;
		 ar += vr[17] * *right++;
		 ar += vr[18] * *right++;
		 ar += vr[19] * *right++;
		 
		 x1px2 = al * ar;
		 v[0] += x1px2 *   *EV++;
		 v[1] += x1px2 *   *EV++;
		 v[2] += x1px2 *   *EV++;
		 v[3] += x1px2 *   *EV++;
		 v[4] += x1px2 *   *EV++;
		 v[5] += x1px2 *   *EV++;
		 v[6] += x1px2 *   *EV++;
		 v[7] += x1px2 *   *EV++;
		 v[8] += x1px2 *   *EV++;
		 v[9] += x1px2 *   *EV++;
		 v[10] += x1px2 *  *EV++;
		 v[11] += x1px2 *  *EV++;
		 v[12] += x1px2 *  *EV++;
		 v[13] += x1px2 *  *EV++;
		 v[14] += x1px2 *  *EV++;
		 v[15] += x1px2 *  *EV++;
		 v[16] += x1px2 *  *EV++;
		 v[17] += x1px2 *  *EV++;
		 v[18] += x1px2 *  *EV++;
		 v[19] += x1px2 *  *EV++;
	       }
	   }
	 
	 ex3[i] = ex1[i] + ex2[i];
	 v = &(x3[80 * i]);
	 scale = 1;
	 for(l = 0; scale && (l < 80); l++)
	   scale = ((ABS(v[l]) <  minlikelihood));
	 
	 if (scale)	         
	   {	     	    	 
	     for(l = 0; l < 80; l++)
	       v[l] *= twotothe256;		   
	     ex3[i] += 1;  
	   }		 	 	 	   
       }
      break;
    default:
      assert(0);
    }
     
  free(left_start); 
  free(right_start); 
}

static void newviewGTRGAMMAPROTMULT(traversalInfo *ti,
				       double *x1, double *x2, double *x3,
				       double *extEIGN, double *extEV, double *extEI, double *gammaRates, double *tipVector,
				       int    *ex1, int *ex2, int *ex3, char *tipX1, char *tipX2, int *modelptr,
				       int lower, int n, int numberOfModels, int multiBranch)
{ 
  double  *left, *right, *left_start, *right_start, *uX1, *uX2, *v, *EV, *EI, *EIGN;
  double x1px2;   
  int  i, l, k, scale;   
  double   
    z1 = 0.0, 
    z2 = 0.0, 
    ki, ump_x1_0;       
  double  d1[19], d2[19]; 
  double *vl, *vr, al, ar; 
  int model;    
     
  if(!multiBranch)
    {
      z1  = ti->qz[0];            
      z2  = ti->rz[0];      
    }
  
  left_start  = left  = (double *)malloc(numberOfModels * 1520 * sizeof(double));
  right_start = right = (double *)malloc(numberOfModels * 1520 * sizeof(double));
  
  for(model = 0; model < numberOfModels; model++)
    {
      if(multiBranch)
	{
	  z1  = ti->qz[model];	 	  
	  z2  = ti->rz[model];	  
	}
      
      EIGN = &(extEIGN[19 * model]);
      
      left  = &left_start[model * 1520];
      right = &right_start[model * 1520];
      
      for(i = 0; i < 4;i++)
	{       	   
	  ki = gammaRates[model * 4 + i];	    
	  
	  for(k = 0; k < 19; k++)
	    {
	      d1[k] = exp (EIGN[k] * ki * z1);
	      d2[k] = exp (EIGN[k] * ki * z2);
	    }
	  
	  EI = &(extEI[model * 380]);
	  
	  for(k = 0; k < 20; k++)
	    {
	      *left++ = d1[0] * *EI++;
	      *left++ = d1[1] * *EI++;
	      *left++ = d1[2] * *EI++;
	      *left++ = d1[3] * *EI++;
	      *left++ = d1[4] * *EI++;
	      *left++ = d1[5] * *EI++;
	      *left++ = d1[6] * *EI++;
	      *left++ = d1[7] * *EI++;
	      *left++ = d1[8] * *EI++;
	      *left++ = d1[9] * *EI++;
	      *left++ = d1[10] * *EI++;
	      *left++ = d1[11] * *EI++;
	      *left++ = d1[12] * *EI++;
	      *left++ = d1[13] * *EI++;
	      *left++ = d1[14] * *EI++;
	      *left++ = d1[15] * *EI++;
	      *left++ = d1[16] * *EI++;
	      *left++ = d1[17] * *EI++;
	      *left++ = d1[18] * *EI++;	      
	    }
	  
	  EI =  &(extEI[model * 380]);
	  
	  for(k = 0; k < 20; k++)
	    {
	      *right++ = d2[0] * *EI++;
	      *right++ = d2[1] * *EI++;
	      *right++ = d2[2] * *EI++;
	      *right++ = d2[3] * *EI++;
	      *right++ = d2[4] * *EI++;
	      *right++ = d2[5] * *EI++;
	      *right++ = d2[6] * *EI++;
	      *right++ = d2[7] * *EI++;
	      *right++ = d2[8] * *EI++;
	      *right++ = d2[9] * *EI++;
	      *right++ = d2[10] * *EI++;
	      *right++ = d2[11] * *EI++;
	      *right++ = d2[12] * *EI++;
	      *right++ = d2[13] * *EI++;
	      *right++ = d2[14] * *EI++;
	      *right++ = d2[15] * *EI++;
	      *right++ = d2[16] * *EI++;
	      *right++ = d2[17] * *EI++;
	      *right++ = d2[18] * *EI++;	      
	    }	
	}               
    }    
  
  switch(ti->tipCase)
    {
    case TIP_TIP:
      {	      
	double *umpX1, *umpX2;
      
	umpX1 = (double *)malloc(1840 * numberOfModels * sizeof(double));
	umpX2 = (double *)malloc(1840 * numberOfModels * sizeof(double));
      
	for(model = 0; model < numberOfModels; model++)
	  {
	    uX1 = &umpX1[1840 * model];
	    uX2 = &umpX2[1840 * model];
	  
	    for(i = 0; i < 23; i++)
	      {	    	     
		v = &(tipVector[model * 460 + i * 20]);	     
		left = &left_start[model * 1520];
		right = &right_start[model * 1520];
	      
		for(k = 0; k < 80; k++)
		  {
		    ump_x1_0 = v[0];
		    for(l = 1; l < 20; l++)
		      ump_x1_0 +=  v[l] * *left++;
		    *uX1++ = ump_x1_0;
		  }
		
		for(k = 0; k < 80; k++)
		  {
		    ump_x1_0 = v[0];
		    for(l = 1; l < 20; l++)
		      ump_x1_0 +=  v[l] * *right++;
		    *uX2++ = ump_x1_0;
		  }
	      }	
	  }
	
	 
	for(i = lower; i < n; i++) 
	  {
	    model = modelptr[i];
	  
	    uX1 = &umpX1[1840 * model + 80 * tipX1[i]];
	    uX2 = &umpX2[1840 * model + 80 * tipX2[i]];
	  	  
	    for(k = 0; k < 4; k++)
	      {
		EV = &(extEV[model * 400]);
		x1px2 = *uX1++ * *uX2++;
		v = &(x3[80 * i + k * 20]);
		for(l = 0; l < 20; l++)
		  v[l] = x1px2 *  *EV++;
		for(l = 0; l < 19; l++)
		  {
		    x1px2 = *uX1++ * *uX2++;
		    v[0] += x1px2 *   *EV++;
		    v[1] += x1px2 *   *EV++;
		    v[2] += x1px2 *   *EV++;
		    v[3] += x1px2 *   *EV++;
		    v[4] += x1px2 *   *EV++;
		    v[5] += x1px2 *   *EV++;
		    v[6] += x1px2 *   *EV++;
		    v[7] += x1px2 *   *EV++;
		    v[8] += x1px2 *   *EV++;
		    v[9] += x1px2 *   *EV++;
		    v[10] += x1px2 *  *EV++;
		    v[11] += x1px2 *  *EV++;
		    v[12] += x1px2 *  *EV++;
		    v[13] += x1px2 *  *EV++;
		    v[14] += x1px2 *  *EV++;
		    v[15] += x1px2 *  *EV++;
		    v[16] += x1px2 *  *EV++;
		    v[17] += x1px2 *  *EV++;
		    v[18] += x1px2 *  *EV++;
		    v[19] += x1px2 *  *EV++;
		  }
	      }
	  
	    ex3[i] = 0;		  		    	     	            
	  }

	free(umpX1);
	free(umpX2);
      }
      break;      
    case TIP_INNER:      
      {	 	
	double *umpX1, ump_x2[20];	
      
	umpX1 = (double *)malloc(1840 * numberOfModels * sizeof(double));
      
	for(model = 0; model < numberOfModels; model++)
	  {
	    uX1 = &umpX1[1840 * model];
	    
	    for(i = 0; i < 23; i++)
	      {	    	     		 
		v = &(tipVector[model * 460 + 20 * i]);	     
		left = &left_start[model * 1520];
		
		for(k = 0; k < 80; k++)
		  {		 		 	
		    ump_x1_0 = v[0];
		    for(l = 1; l < 20; l++)		 
		      ump_x1_0 +=  v[l] * *left++;		 
		    *uX1++ = ump_x1_0;
		  }
	      }
	  }
	     
     
	for (i = lower; i < n; i++) 
	  {		   
	    model = modelptr[i];
	    
	    right = &right_start[model * 1520];
	    uX1 = &umpX1[1840 * model + 80 * tipX1[i]];			    	    
	    
	    for(k = 0; k < 4; k++)
	      {
		v = &(x2[80 * i + k * 20]);	    
		for(l = 0; l < 20; l++)
		  {		     
		    ump_x2[l] =  v[0];		 
		    ump_x2[l] += v[1] * *right++;
		    ump_x2[l] += v[2] * *right++;
		    ump_x2[l] += v[3] * *right++;
		    ump_x2[l] += v[4] * *right++;
		    ump_x2[l] += v[5] * *right++;
		    ump_x2[l] += v[6] * *right++;
		    ump_x2[l] += v[7] * *right++;
		    ump_x2[l] += v[8] * *right++;
		    ump_x2[l] += v[9] * *right++;
		    ump_x2[l] += v[10] * *right++;
		    ump_x2[l] += v[11] * *right++;
		    ump_x2[l] += v[12] * *right++;
		    ump_x2[l] += v[13] * *right++;
		    ump_x2[l] += v[14] * *right++;
		    ump_x2[l] += v[15] * *right++;
		    ump_x2[l] += v[16] * *right++;
		    ump_x2[l] += v[17] * *right++;
		    ump_x2[l] += v[18] * *right++;
		    ump_x2[l] += v[19] * *right++;
		  }
		 	     
		uX2 = ump_x2;
		EV = &(extEV[model * 400]);
		x1px2 = *uX1++ * *uX2++;
		v = &(x3[80 * i + 20 * k]);
		for(l = 0; l < 20; l++)
		  v[l] = x1px2 *  *EV++;
		
		for(l = 0; l < 19; l++)
		  {
		    x1px2 = *uX1++ * *uX2++;
		    v[0] += x1px2 *   *EV++;
		    v[1] += x1px2 *   *EV++;
		    v[2] += x1px2 *   *EV++;
		    v[3] += x1px2 *   *EV++;
		    v[4] += x1px2 *   *EV++;
		    v[5] += x1px2 *   *EV++;
		    v[6] += x1px2 *   *EV++;
		    v[7] += x1px2 *   *EV++;
		    v[8] += x1px2 *   *EV++;
		    v[9] += x1px2 *   *EV++;
		    v[10] += x1px2 *  *EV++;
		    v[11] += x1px2 *  *EV++;
		    v[12] += x1px2 *  *EV++;
		    v[13] += x1px2 *  *EV++;
		    v[14] += x1px2 *  *EV++;
		    v[15] += x1px2 *  *EV++;
		    v[16] += x1px2 *  *EV++;
		    v[17] += x1px2 *  *EV++;
		    v[18] += x1px2 *  *EV++;
		    v[19] += x1px2 *  *EV++;
		  }
	      }	
	    
	    ex3[i] = ex2[i];		
	    v = &x3[80 * i];
	    scale = 1;
	    for(l = 0; scale && (l < 80); l++)
	      scale = (ABS(v[l]) <  minlikelihood);
	  
	    if (scale)	         
	      {	     	    	 
		for(l = 0; l < 80; l++)
		  v[l] *= twotothe256;		   
		ex3[i] += 1;  
	      }	 	           	 
	  }           
	free(umpX1);     
      }
      break;
    case INNER_INNER:
	     
      for (i = lower; i < n; i++) 
	{		
	  model = modelptr[i];
      
	  left =  &left_start[model * 1520];
	  right = &right_start[model * 1520];
            
	  for(k = 0; k < 4; k++)
	    {
	      EV =  &(extEV[model * 400]);
	      vl = &(x1[80 * i + 20 * k]);
	      vr = &(x2[80 * i + 20 * k]);
	      v =  &(x3[80 * i + 20 * k]);
	      
	      for(l = 0; l < 20; l++)	       
		v[l] = 0;		 
	      
	      for(l = 0; l < 20; l++)
		{
		  al =  vl[0];		 
		  al += vl[1] * *left++;
		  al += vl[2] * *left++;
		  al += vl[3] * *left++;
		  al += vl[4] * *left++;
		  al += vl[5] * *left++;
		  al += vl[6] * *left++;
		  al += vl[7] * *left++;
		  al += vl[8] * *left++;
		  al += vl[9] * *left++;
		  al += vl[10] * *left++;
		  al += vl[11] * *left++;
		  al += vl[12] * *left++;
		  al += vl[13] * *left++;
		  al += vl[14] * *left++;
		  al += vl[15] * *left++;
		  al += vl[16] * *left++;
		  al += vl[17] * *left++;
		  al += vl[18] * *left++;
		  al += vl[19] * *left++;
		  
		  ar =  vr[0];		 
		  ar += vr[1] * *right++;
		  ar += vr[2] * *right++;
		  ar += vr[3] * *right++;
		  ar += vr[4] * *right++;
		  ar += vr[5] * *right++;
		  ar += vr[6] * *right++;
		  ar += vr[7] * *right++;
		  ar += vr[8] * *right++;
		  ar += vr[9] * *right++;
		  ar += vr[10] * *right++;
		  ar += vr[11] * *right++;
		  ar += vr[12] * *right++;
		  ar += vr[13] * *right++;
		  ar += vr[14] * *right++;
		  ar += vr[15] * *right++;
		  ar += vr[16] * *right++;
		  ar += vr[17] * *right++;
		  ar += vr[18] * *right++;
		  ar += vr[19] * *right++;
		  
		  x1px2 = al * ar;
		  v[0] += x1px2 *   *EV++;
		  v[1] += x1px2 *   *EV++;
		  v[2] += x1px2 *   *EV++;
		  v[3] += x1px2 *   *EV++;
		  v[4] += x1px2 *   *EV++;
		  v[5] += x1px2 *   *EV++;
		  v[6] += x1px2 *   *EV++;
		  v[7] += x1px2 *   *EV++;
		  v[8] += x1px2 *   *EV++;
		  v[9] += x1px2 *   *EV++;
		  v[10] += x1px2 *  *EV++;
		  v[11] += x1px2 *  *EV++;
		  v[12] += x1px2 *  *EV++;
		  v[13] += x1px2 *  *EV++;
		  v[14] += x1px2 *  *EV++;
		  v[15] += x1px2 *  *EV++;
		  v[16] += x1px2 *  *EV++;
		  v[17] += x1px2 *  *EV++;
		  v[18] += x1px2 *  *EV++;
		  v[19] += x1px2 *  *EV++;
		}
	    }
	 
	  ex3[i] = ex1[i] + ex2[i];
	  v = &(x3[80 * i]);
	  scale = 1;
	  for(l = 0; scale && (l < 80); l++)
	    scale = ((ABS(v[l]) <  minlikelihood));
	  
	  if(scale)	         
	    {	     	    	 
	      for(l = 0; l < 80; l++)
		v[l] *= twotothe256;		   
	      ex3[i] += 1;  
	    }		 	 	 	   
	}
      break;
    default:
      assert(0);
    }
     
  free(left_start); 
  free(right_start);  
}


#ifdef _MULTI_GENE
void computeMultiTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int model)
{
  if(isTip(p->number, maxTips)) 
    {
      assert(p->backs[model]);
      return;   
    }

  if(!p->backs[model])
    {
      assert(0);     
    }
  else
    {       
      nodeptr q = p->next->backs[model];
      nodeptr r = p->next->next->backs[model];

      assert(p == p->next->next->next);
      
      assert(q && r);
      
      if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
	{	  
	  while (! p->xs[model])
	    {
	      if (! p->xs[model])
		getxsnode(p, model); 	   
	    }

	  assert(p->xs[model]);
	  
	  ti[*counter].tipCase = TIP_TIP; 
	  ti[*counter].pNumber = p->number;
	  ti[*counter].qNumber = q->number;
	  ti[*counter].rNumber = r->number;
	  
	  {
	    double z;
	    z = q->z[model];
	    z = (z > zmin) ? log(z) : log(zmin);
	    ti[*counter].qz[model] = z;
	    
	    z = r->z[model];
	    z = (z > zmin) ? log(z) : log(zmin);
	    ti[*counter].rz[model] = z;	    
	  }     
	  
	  *counter = *counter + 1;	  
	}  
      else
	{
	  if(isTip(r->number, maxTips) || isTip(q->number, maxTips))
	    {		
	      nodeptr tmp;
	      
	      if(isTip(r->number, maxTips))
		{
		  tmp = r;
		  r = q;
		  q = tmp;
		}
	      
	      while ((! p->xs[model]) || (! r->xs[model])) 
		{	 	    
		  if (! r->xs[model]) 
		    computeMultiTraversalInfo(r, ti, counter, maxTips, model);
		  if (! p->xs[model]) 
		    getxsnode(p, model);	
		}
	      
	      assert(p->xs[model] && r->xs[model]);

	      ti[*counter].tipCase = TIP_INNER; 
	      ti[*counter].pNumber = p->number;
	      ti[*counter].qNumber = q->number;
	      ti[*counter].rNumber = r->number;
	      
	      {
		double z;
		z = q->z[model];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].qz[model] = z;
		
		z = r->z[model];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].rz[model] = z;		
	      }   
	      
	      *counter = *counter + 1;	      
	    }
	  else
	    {	 
	      
	      while ((! p->xs[model]) || (! q->xs[model]) || (! r->xs[model])) 
		{
		  if (! q->xs[model]) 
		    computeMultiTraversalInfo(q, ti, counter, maxTips, model);
		  if (! r->xs[model]) 
		    computeMultiTraversalInfo(r, ti, counter, maxTips, model);
		  if (! p->xs[model]) 
		    getxsnode(p, model);	
		}
	      assert(p->xs[model] && r->xs[model] && q->xs[model]);
	      ti[*counter].tipCase = INNER_INNER; 
	      ti[*counter].pNumber = p->number;
	      ti[*counter].qNumber = q->number;
	      ti[*counter].rNumber = r->number;
	      
	      {
		double z;
		z = q->z[model];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].qz[model] = z;
		
		z = r->z[model];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].rz[model] = z;		
	      }   
	      
	      *counter = *counter + 1;	      
	    }
	}
    }
}


#endif

void computeTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches)
{
  if(isTip(p->number, maxTips))
    return;

  {     
    int i;
    nodeptr q = p->next->back;
    nodeptr r = p->next->next->back;
    
    if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
      {	  
	while (! p->x)
	 {
	   if (! p->x)
	     getxnode(p); 	   
	 }

	ti[*counter].tipCase = TIP_TIP; 
	ti[*counter].pNumber = p->number;
	ti[*counter].qNumber = q->number;
	ti[*counter].rNumber = r->number;
	for(i = 0; i < numBranches; i++)
	  {
	    double z;
	    z = q->z[i];
	    z = (z > zmin) ? log(z) : log(zmin);
	    ti[*counter].qz[i] = z;

	    z = r->z[i];
	    z = (z > zmin) ? log(z) : log(zmin);
	    ti[*counter].rz[i] = z;	    
	  }     
	*counter = *counter + 1;
      }  
    else
      {
	if(isTip(r->number, maxTips) || isTip(q->number, maxTips))
	  {		
	    nodeptr tmp;

	    if(isTip(r->number, maxTips))
	      {
		tmp = r;
		r = q;
		q = tmp;
	      }

	    while ((! p->x) || (! r->x)) 
	      {	 	    
		if (! r->x) 
		  computeTraversalInfo(r, ti, counter, maxTips, numBranches);
		if (! p->x) 
		  getxnode(p);	
	      }
	    	   
	    ti[*counter].tipCase = TIP_INNER; 
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;
	    for(i = 0; i < numBranches; i++)
	      {
		double z;
		z = q->z[i];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].qz[i] = z;
		
		z = r->z[i];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].rz[i] = z;		
	      }   
	    
	    *counter = *counter + 1;
	  }
	else
	  {	 

	    while ((! p->x) || (! q->x) || (! r->x)) 
	      {
		if (! q->x) 
		  computeTraversalInfo(q, ti, counter, maxTips, numBranches);
		if (! r->x) 
		  computeTraversalInfo(r, ti, counter, maxTips, numBranches);
		if (! p->x) 
		  getxnode(p);	
	      }
   
	    ti[*counter].tipCase = INNER_INNER; 
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;
	    for(i = 0; i < numBranches; i++)
	      {
		double z;
		z = q->z[i];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].qz[i] = z;
		
		z = r->z[i];
		z = (z > zmin) ? log(z) : log(zmin);
		ti[*counter].rz[i] = z;		
	      }   
	    
	    *counter = *counter + 1;
	  }
      }    
  }

}


static void newviewMixedData(int model, tree *tr, traversalInfo *tInfo)
{
  double 
    *x1_start = (double*)NULL,
    *x2_start = (double*)NULL, 
    *x3_start = (double*)NULL;
  int 
    *ex1 = (int*)NULL, 
    *ex2 = (int*)NULL, 
    *ex3 = (int*)NULL;  
  char 
    *tipX1 = (char *)NULL,
    *tipX2 = (char *)NULL;

  int branchIndex;
  /*int l = tr->modelIndices[model][0];
    int u = tr->modelIndices[model][1];	  	*/
  int l = tr->partitionData[model].lower;
  int u = tr->partitionData[model].upper;


  int width  = u - l;
  /*int offset = tr->modelOffsets[model]; */
  int offset = tr->partitionData[model].modelOffset;
  
  if(tr->multiBranch)
    branchIndex = model;
  else
    branchIndex = 0;
  
  switch(tInfo->tipCase)
    {
    case TIP_TIP:
      tipX1    = tr->yVector[tInfo->qNumber];		 
      tipX1    = &tipX1[l];
      
      tipX2    = tr->yVector[tInfo->rNumber];
      tipX2    = &tipX2[l];
      
      x3_start = getLikelihoodArray(tInfo->pNumber,  tr->mxtips, tr->xVector);
      ex3      = getScalingArray(tInfo->pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
      x3_start = &x3_start[offset];
      ex3      = &ex3[l];
      
      break;
    case TIP_INNER:
      tipX1    = tr->yVector[tInfo->qNumber];	 
      tipX1    = &tipX1[l];
      
      x2_start = getLikelihoodArray(tInfo->rNumber,  tr->mxtips, tr->xVector);	  
      ex2      = getScalingArray(tInfo->rNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
      x2_start = &x2_start[offset];
      ex2      = &ex2[l];
      
      x3_start = getLikelihoodArray(tInfo->pNumber,  tr->mxtips, tr->xVector);
      ex3      = getScalingArray(tInfo->pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
      x3_start = &x3_start[offset];
      ex3      = &ex3[l];
      break;
    case INNER_INNER:	 
      x1_start = getLikelihoodArray(tInfo->qNumber,  tr->mxtips, tr->xVector);
      ex1      = getScalingArray(tInfo->qNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
      x1_start = &x1_start[offset];
      ex1      = &ex1[l];
      
      x2_start = getLikelihoodArray(tInfo->rNumber,  tr->mxtips, tr->xVector);
      ex2      = getScalingArray(tInfo->rNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);	  
      x2_start = &x2_start[offset];
      ex2      = &ex2[l];
      
      x3_start = getLikelihoodArray(tInfo->pNumber,  tr->mxtips, tr->xVector);
      ex3      = getScalingArray(tInfo->pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
      x3_start = &x3_start[offset];
      ex3      = &ex3[l];
      break;
    default:
      assert(0);
    }
	      
  switch(/*tr->dataType[model]*/ tr->partitionData[model].dataType)
    { 
    case DNA_DATA:
      switch(tr->rateHetModel)
	{
	case CAT:
	  newviewGTRCAT(tInfo,  &(tr->EV_DNA[model * 16]), &(tr->EI_DNA[model * 12]), 
			&(tr->EIGN_DNA[model * 3]), &(tr->cdta->patrat[0]), &(tr->cdta->rateCategory[l]), 
			x1_start, x2_start, x3_start, 
			&(tr->tipVectorDNA[model * 64]),
			ex1, ex2, ex3, 
			tipX1, tipX2,
			0, width, tr->NumberOfCategories, tInfo->qz[branchIndex], tInfo->rz[branchIndex]);
	  break;
	case GAMMA:
	case GAMMA_I:
	   newviewGTRGAMMA(tInfo,
			   x1_start, x2_start, x3_start,
			   &(tr->EIGN_DNA[model * 3]), &(tr->EV_DNA[model * 16]), 
			   &(tr->EI_DNA[model * 12]), &(tr->gammaRates[model * 4]), 
			   &(tr->tipVectorDNA[model * 64]),
			   ex1, ex2, ex3, tipX1, tipX2,
			   0, width, tInfo->qz[branchIndex], tInfo->rz[branchIndex]);
	  break;	
	default:
	  assert(0);
	}
      break;
    case AA_DATA: 
      switch(tr->rateHetModel)
	{
	case CAT:
	  newviewGTRCATPROT(tInfo,  &(tr->EV_AA[400 * model]), &(tr->EI_AA[380 * model]), 
			    &(tr->EIGN_AA[19 * model]), &(tr->cdta->patrat[0]), 
			    &(tr->cdta->rateCategory[l]), 
			    x1_start, x2_start, x3_start, 
			    &(tr->tipVectorAA[model * 460]),
			    ex1, ex2, ex3, 
			    tipX1, tipX2,
			    0, width, tr->NumberOfCategories, tInfo->qz[branchIndex], tInfo->rz[branchIndex]);
	  break;
	case GAMMA:	
	case GAMMA_I:
	  newviewGTRGAMMAPROT(tInfo,
			      x1_start, x2_start, x3_start,
			      &(tr->EIGN_AA[model * 19]), &(tr->EV_AA[model * 400]), 
			      &(tr->EI_AA[model * 380]), &(tr->gammaRates[model * 4]), 
			      &(tr->tipVectorAA[model * 460]),
			      ex1, ex2, ex3, tipX1, tipX2,
			      0, width, tInfo->qz[branchIndex], tInfo->rz[branchIndex]);
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

void newviewIterative (tree *localTree, int startIndex, int endIndex)
{	   
  /* LTD */

  traversalInfo *ti   = localTree->td[0].ti;
  int i;  
   
  for(i = 1; i < localTree->td[0].count; i++)
    {    
      traversalInfo *tInfo = &ti[i];
      
      if(localTree->mixedData)
	{                                            	  
	  int model;
	  
	  assert(0);

	  for(model = 0; model < localTree->NumberOfModels; model++)
	    newviewMixedData(model, localTree, tInfo);	     
	}
      else
	{
	  double 
	    *x1_start = (double*)NULL,
	    *x2_start = (double*)NULL, 
	    *x3_start = (double*)NULL;
	  int 
	    *ex1 = (int*)NULL, 
	    *ex2 = (int*)NULL, 
	    *ex3 = (int*)NULL;  
	  char 
	    *tipX1 = (char *)NULL,
	    *tipX2 = (char *)NULL;
	  
	  switch(tInfo->tipCase)
	    {
	    case TIP_TIP:
	      tipX1    = &localTree->strided_yVector[tInfo->qNumber][startIndex];	      
	      tipX2    = &localTree->strided_yVector[tInfo->rNumber][startIndex];
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber,  localTree->mxtips, localTree->xVector);
	      ex3      = getScalingArray(tInfo->pNumber,     localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      break;
	    case TIP_INNER:
	      tipX1    = &localTree->strided_yVector[tInfo->qNumber][startIndex];
	      
	      x2_start = getLikelihoodArray(tInfo->rNumber,  localTree->mxtips, localTree->xVector);	  
	      ex2      = getScalingArray(tInfo->rNumber,     localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber,   localTree->mxtips, localTree->xVector);
	      ex3      = getScalingArray(tInfo->pNumber,      localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      break;
	    case INNER_INNER:	 
	      x1_start = getLikelihoodArray(tInfo->qNumber,  localTree->mxtips, localTree->xVector);
	      ex1      = getScalingArray(tInfo->qNumber,     localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      x2_start = getLikelihoodArray(tInfo->rNumber,  localTree->mxtips, localTree->xVector);
	      ex2      = getScalingArray(tInfo->rNumber,     localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber, localTree->mxtips, localTree->xVector); 
	      ex3      = getScalingArray(tInfo->pNumber,    localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      break;
	    default:
	      assert(0);
	    }

	  switch(localTree->likelihoodFunction)
	    {
	    case GTRCAT:	 
	      newviewGTRCAT(tInfo,  localTree->EV_DNA, localTree->EI_DNA, localTree->EIGN_DNA, 
			    localTree->strided_patrat, &(localTree->strided_rateCategory[startIndex]), 
			    x1_start, x2_start, x3_start, localTree->tipVectorDNA,
			    ex1, ex2, ex3, tipX1, tipX2,
			    0, (endIndex - startIndex), localTree->NumberOfCategories, tInfo->qz[0], tInfo->rz[0]
			    );   
	      break;
	    case GTRCATMULT:		     	     
	      newviewGTRCATMULT(tInfo,  localTree->EV_DNA, localTree->EI_DNA, localTree->EIGN_DNA, 
				localTree->strided_patrat, 
				&(localTree->strided_rateCategory[startIndex]), 
				x1_start, x2_start, x3_start, localTree->tipVectorDNA,
				ex1, ex2, ex3, tipX1, tipX2, &(localTree->strided_model[startIndex]),
				0, (endIndex - startIndex), localTree->NumberOfCategories, localTree->NumberOfModels, 
				localTree->multiBranch);
	      break;
	    case PROTCAT:	    
	      newviewGTRCATPROT(tInfo,  localTree->EV_AA, localTree->EI_AA, localTree->EIGN_AA, localTree->strided_patrat, 
				&(localTree->strided_rateCategory[startIndex]), 
				x1_start, x2_start, x3_start, localTree->tipVectorAA,
				ex1, ex2, ex3, tipX1, tipX2,
				0, (endIndex - startIndex), localTree->NumberOfCategories, tInfo->qz[0], tInfo->rz[0]);
	      break;
	    case PROTCATMULT:
	      newviewGTRCATPROTMULT(tInfo,  localTree->EV_AA, localTree->EI_AA, localTree->EIGN_AA, localTree->strided_patrat, 
				    &(localTree->strided_rateCategory[startIndex]), 
				    x1_start, x2_start, x3_start, localTree->tipVectorAA,
				    ex1, ex2, ex3, tipX1, tipX2, &(localTree->strided_model[startIndex]),
				    0, (endIndex - startIndex), localTree->NumberOfCategories, localTree->NumberOfModels, 
				    localTree->multiBranch); 
	      break;
	    case GTRGAMMA:
	    case GTRGAMMAI:
	      newviewGTRGAMMA(tInfo,
			      x1_start, x2_start, x3_start,
			      localTree->EIGN_DNA, localTree->EV_DNA, localTree->EI_DNA, localTree->gammaRates, 
			      localTree->tipVectorDNA,
			      ex1, ex2, ex3, tipX1, tipX2,
			      0, (endIndex - startIndex), tInfo->qz[0], tInfo->rz[0]);
	      break;
	    case GTRGAMMAMULT:
	    case GTRGAMMAMULTI:
	      newviewGTRGAMMAMULT(tInfo,
				  x1_start, x2_start, x3_start,
				  localTree->EIGN_DNA, localTree->EV_DNA, localTree->EI_DNA, localTree->gammaRates, 
				  localTree->tipVectorDNA,
				  ex1, ex2, ex3, tipX1, tipX2, &(localTree->strided_model[startIndex]),
				  0, (endIndex - startIndex), localTree->NumberOfModels, localTree->multiBranch);      
	      break;
	    case PROTGAMMA:
	    case PROTGAMMAI:
	      newviewGTRGAMMAPROT(tInfo,
				  x1_start, x2_start, x3_start,
				  localTree->EIGN_AA, localTree->EV_AA, localTree->EI_AA, localTree->gammaRates, 
				  localTree->tipVectorAA,
				  ex1, ex2, ex3, tipX1, tipX2,
				  0, (endIndex - startIndex), tInfo->qz[0], tInfo->rz[0]);     
	      break;
	    case PROTGAMMAMULT:
	    case PROTGAMMAMULTI:
	      newviewGTRGAMMAPROTMULT(tInfo,
				      x1_start, x2_start, x3_start,
				      localTree->EIGN_AA, localTree->EV_AA, localTree->EI_AA, localTree->gammaRates, 
				      localTree->tipVectorAA,
				      ex1, ex2, ex3, tipX1, tipX2, &(localTree->strided_model[startIndex]),
				      0, (endIndex - startIndex), localTree->NumberOfModels, localTree->multiBranch);       
	      break;
	    default:
	      assert(0);
	    }
	}
    }
}


#else

void newviewIterative (tree *tr, int startIndex, int endIndex)
{	   
  traversalInfo *ti   = tr->td[0].ti;
  int i;  
   
  for(i = 1; i < tr->td[0].count; i++)
    {    
      traversalInfo *tInfo = &ti[i];
      
      if(tr->mixedData)
	{                                            	  
	  int model;
	  
	  for(model = 0; model < tr->NumberOfModels; model++)
	    newviewMixedData(model, tr, tInfo);
	     
	}
      else
	{
	  double 
	    *x1_start = (double*)NULL,
	    *x2_start = (double*)NULL, 
	    *x3_start = (double*)NULL;
	  int 
	    *ex1 = (int*)NULL, 
	    *ex2 = (int*)NULL, 
	    *ex3 = (int*)NULL;  
	  char 
	    *tipX1 = (char *)NULL,
	    *tipX2 = (char *)NULL;
	  
	  switch(tInfo->tipCase)
	    {
	    case TIP_TIP:
	      tipX1    = tr->yVector[tInfo->qNumber];
	      
	      tipX2    = tr->yVector[tInfo->rNumber];
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber,  tr->mxtips, tr->xVector);
	      ex3      = getScalingArray(tInfo->pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
	      
	      break;
	    case TIP_INNER:
	      tipX1    = tr->yVector[tInfo->qNumber];
	      
	      x2_start = getLikelihoodArray(tInfo->rNumber,  tr->mxtips, tr->xVector);	  
	      ex2      = getScalingArray(tInfo->rNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber,  tr->mxtips, tr->xVector);
	      ex3      = getScalingArray(tInfo->pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
	      
	      break;
	    case INNER_INNER:	 
	      x1_start = getLikelihoodArray(tInfo->qNumber,  tr->mxtips, tr->xVector);
	      ex1      = getScalingArray(tInfo->qNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);	 
	      
	      x2_start = getLikelihoodArray(tInfo->rNumber,  tr->mxtips, tr->xVector);
	      ex2      = getScalingArray(tInfo->rNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);	  
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber,  tr->mxtips, tr->xVector);
	      ex3      = getScalingArray(tInfo->pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
	      
	      break;
	    default:
	      assert(0);
	    }

	  switch(tr->likelihoodFunction)
	    {
	    case GTRCAT:	 
	      newviewGTRCAT(tInfo,  tr->EV_DNA, tr->EI_DNA, tr->EIGN_DNA, tr->cdta->patrat, tr->cdta->rateCategory, 
			    x1_start, x2_start, x3_start, tr->tipVectorDNA,
			    ex1, ex2, ex3, tipX1, tipX2,
			    startIndex, endIndex, tr->NumberOfCategories, tInfo->qz[0], tInfo->rz[0]
			    );   
	      break;
	    case GTRCATMULT:	   
	      newviewGTRCATMULT(tInfo,  tr->EV_DNA, tr->EI_DNA, tr->EIGN_DNA, tr->cdta->patrat, tr->cdta->rateCategory, 
				x1_start, x2_start, x3_start, tr->tipVectorDNA,
				ex1, ex2, ex3, tipX1, tipX2, tr->model,
				startIndex, endIndex, tr->NumberOfCategories, tr->NumberOfModels, tr->multiBranch);
	      break;
	    case PROTCAT:	    
	      newviewGTRCATPROT(tInfo,  tr->EV_AA, tr->EI_AA, tr->EIGN_AA, tr->cdta->patrat, tr->cdta->rateCategory, 
				x1_start, x2_start, x3_start, tr->tipVectorAA,
				ex1, ex2, ex3, tipX1, tipX2,
				startIndex, endIndex, tr->NumberOfCategories, tInfo->qz[0], tInfo->rz[0]);
	      break;
	    case PROTCATMULT:
	      newviewGTRCATPROTMULT(tInfo,  tr->EV_AA, tr->EI_AA, tr->EIGN_AA, tr->cdta->patrat, tr->cdta->rateCategory, 
				    x1_start, x2_start, x3_start, tr->tipVectorAA,
				    ex1, ex2, ex3, tipX1, tipX2, tr->model,
				    startIndex, endIndex, tr->NumberOfCategories, tr->NumberOfModels, tr->multiBranch); 
	      break;
	    case GTRGAMMA:
	    case GTRGAMMAI:
	      newviewGTRGAMMA(tInfo,
			      x1_start, x2_start, x3_start,
			      tr->EIGN_DNA, tr->EV_DNA, tr->EI_DNA, tr->gammaRates, tr->tipVectorDNA,
			      ex1, ex2, ex3, tipX1, tipX2,
			      startIndex, endIndex, tInfo->qz[0], tInfo->rz[0]);
	      break;
	    case GTRGAMMAMULT:
	    case GTRGAMMAMULTI:
	      newviewGTRGAMMAMULT(tInfo,
				  x1_start, x2_start, x3_start,
				  tr->EIGN_DNA, tr->EV_DNA, tr->EI_DNA, tr->gammaRates, tr->tipVectorDNA,
				  ex1, ex2, ex3, tipX1, tipX2, tr->model,
				  startIndex, endIndex, tr->NumberOfModels, tr->multiBranch);      
	      break;
	    case PROTGAMMA:
	    case PROTGAMMAI:
	      newviewGTRGAMMAPROT(tInfo,
				  x1_start, x2_start, x3_start,
				  tr->EIGN_AA, tr->EV_AA, tr->EI_AA, tr->gammaRates, tr->tipVectorAA,
				  ex1, ex2, ex3, tipX1, tipX2,
				  startIndex, endIndex, tInfo->qz[0], tInfo->rz[0]);     
	      break;
	    case PROTGAMMAMULT:
	    case PROTGAMMAMULTI:
	      newviewGTRGAMMAPROTMULT(tInfo,
				      x1_start, x2_start, x3_start,
				      tr->EIGN_AA, tr->EV_AA, tr->EI_AA, tr->gammaRates, tr->tipVectorAA,
				      ex1, ex2, ex3, tipX1, tipX2, tr->model,
				      startIndex, endIndex, tr->NumberOfModels, tr->multiBranch);       
	      break;
	    default:
	      assert(0);
	    }
	}
    }
}
#endif

void newviewGeneric (tree *tr, nodeptr p)
{
  if(isTip(p->number, tr->rdta->numsp)) 
    return;

#ifdef _MULTI_GENE
  if(tr->doMulti)
    {
      int model;

      for(model = 0; model < tr->numBranches; model++)
	{
	  tr->td[model].count = 1;
	  computeMultiTraversalInfo(p, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->rdta->numsp, model);   
	  if(tr->td[model].count > 1)
	    newviewIterativePartition(tr, tr->partitionData[model].lower, tr->partitionData[model].lower, model);
	}
      
    }
  else
#endif
  { 	           
    tr->td[0].count = 1;
    computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp, tr->numBranches);   

    if(tr->td[0].count > 1)
      {
#ifdef _USE_PTHREADS  
	masterBarrier(THREAD_NEWVIEW, tr);         
#else
	newviewIterative(tr, 0, tr->cdta->endsite);   
#endif
      }   
  }
}




#ifdef _LOCAL_DATA

void newviewIterativePartition(tree *localTree, int lower, int upper, int model)
{	   
  traversalInfo *ti   = localTree->td[0].ti;
  int i, branchIndex;  
       
  if(localTree->multiBranch)
    branchIndex = model;
  else
    branchIndex = 0;
 

  for(i = 1; i < localTree->td[0].count; i++)
    {    
      traversalInfo *tInfo = &ti[i];

      if(localTree->mixedData)
	{
	  assert(0);
	  newviewMixedData(model, localTree, tInfo);
	}
      else
	{
	  double 
	    *x1_start = (double*)NULL,
	    *x2_start = (double*)NULL, 
	    *x3_start = (double*)NULL;
	  int 
	    *ex1 = (int*)NULL, 
	    *ex2 = (int*)NULL, 
	    *ex3 = (int*)NULL;  
	  char 
	    *tipX1 = (char *)NULL,
	    *tipX2 = (char *)NULL;
	  
	  switch(tInfo->tipCase)
	    {
	    case TIP_TIP:
	      tipX1    = localTree->strided_yVector[tInfo->qNumber];	      
	      tipX2    = localTree->strided_yVector[tInfo->rNumber];
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber,  localTree->mxtips, localTree->xVector);
	      ex3      = getScalingArray(tInfo->pNumber,     localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      break;
	    case TIP_INNER:
	      tipX1    = localTree->strided_yVector[tInfo->qNumber];
	      
	      x2_start = getLikelihoodArray(tInfo->rNumber,  localTree->mxtips, localTree->xVector);	  
	      ex2      = getScalingArray(tInfo->rNumber,     localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber,   localTree->mxtips, localTree->xVector);
	      ex3      = getScalingArray(tInfo->pNumber,      localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      break;
	    case INNER_INNER:	 
	      x1_start = getLikelihoodArray(tInfo->qNumber,  localTree->mxtips, localTree->xVector);
	      ex1      = getScalingArray(tInfo->qNumber,     localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      x2_start = getLikelihoodArray(tInfo->rNumber,  localTree->mxtips, localTree->xVector);
	      ex2      = getScalingArray(tInfo->rNumber,     localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber, localTree->mxtips, localTree->xVector); 
	      ex3      = getScalingArray(tInfo->pNumber,    localTree->mySpan, localTree->mxtips, localTree->expArray);
	      
	      break;
	    default:
	      assert(0);
	    }
      
	  switch(localTree->likelihoodFunction)
	    {
	    case GTRCAT:
	    case GTRCATMULT:
	      newviewGTRCAT(tInfo,  &(localTree->EV_DNA[model * 16]), &(localTree->EI_DNA[model * 12]), 
			    &(localTree->EIGN_DNA[model * 3]), 
			    localTree->strided_patrat, localTree->strided_rateCategory, 
			    x1_start, x2_start, x3_start, &(localTree->tipVectorDNA[model * 64]),
			    ex1, ex2, ex3, tipX1, tipX2,
			    lower, upper, localTree->NumberOfCategories, tInfo->qz[branchIndex], tInfo->rz[branchIndex]
			    );   
	      break;
	    case GTRGAMMA:  /* needed for rate opt*/
	    case GTRGAMMAI: /* needed for rate opt*/
	    case GTRGAMMAMULT:
	    case GTRGAMMAMULTI:
	      newviewGTRGAMMA(tInfo,
			      x1_start, x2_start, x3_start,
			      &(localTree->EIGN_DNA[model * 3]), &(localTree->EV_DNA[model * 16]), 
			      &(localTree->EI_DNA[model * 12]), &(localTree->gammaRates[model * 4]), 
			      &(localTree->tipVectorDNA[model * 64]),
			      ex1, ex2, ex3, tipX1, tipX2,
			      lower, upper, tInfo->qz[branchIndex], tInfo->rz[branchIndex]);
	      break;
	    case PROTGAMMA:  /* needed for rate opt*/
	    case PROTGAMMAI: /* needed for rate opt*/
	    case PROTGAMMAMULT:
	    case PROTGAMMAMULTI:
	      newviewGTRGAMMAPROT(tInfo,
				  x1_start, x2_start, x3_start,
				  &(localTree->EIGN_AA[model * 19]), &(localTree->EV_AA[model * 400]), 
				  &(localTree->EI_AA[model * 380]), &(localTree->gammaRates[model * 4]), 
				  &(localTree->tipVectorAA[model * 460]),
				  ex1, ex2, ex3, tipX1, tipX2,
				  lower, upper, tInfo->qz[branchIndex], tInfo->rz[branchIndex]);     
	      break;	 
	    default:
	      assert(0);
	    }
	}
    }
}


#else

void newviewIterativePartition(tree *tr, int lower, int upper, int model)
{	   
  traversalInfo *ti;
  int i, branchIndex, count;  
       
#ifdef _MULTI_GENE
  if(tr->doMulti)
    {
      assert(tr->multiBranch);
      count = tr->td[model].count;
      ti  = &(tr->td[model].ti[0]);
    }
  else    
#endif
    {
      count = tr->td[0].count;    
      ti    = tr->td[0].ti;
    }

      

  if(tr->multiBranch)
    branchIndex = model;
  else
    branchIndex = 0;

  for(i = 1; i < count; i++)
    {    
      traversalInfo *tInfo = &ti[i];

      if(tr->mixedData)
	newviewMixedData(model, tr, tInfo);
      else
	{
	  double 
	    *x1_start = (double*)NULL,
	    *x2_start = (double*)NULL, 
	    *x3_start = (double*)NULL;
	  int 
	    *ex1 = (int*)NULL, 
	    *ex2 = (int*)NULL, 
	    *ex3 = (int*)NULL;  
	  char 
	    *tipX1 = (char *)NULL,
	    *tipX2 = (char *)NULL;
	  
#ifdef _MULTI_GENE
	  /*if(tr->doMulti)
	    printf("Doing %d %d into %d\n", tInfo->qNumber, tInfo->rNumber, tInfo->pNumber);*/
#endif

	  switch(tInfo->tipCase)
	    {
	    case TIP_TIP:
#ifdef _MULTI_GENE
	      /*printf("TIP_TIP %d %d\n", tInfo->qNumber, tInfo->rNumber);*/
#endif
	      tipX1    = tr->yVector[tInfo->qNumber];
	      
	      tipX2    = tr->yVector[tInfo->rNumber];
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber,  tr->mxtips, tr->xVector);
	      ex3      = getScalingArray(tInfo->pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
	      
	      break;
	    case TIP_INNER:
#ifdef _MULTI_GENE
	      /*printf("TIP_INNER %d %d\n", tInfo->qNumber, tInfo->rNumber);*/
#endif
	      tipX1    = tr->yVector[tInfo->qNumber];	 
	      
	      x2_start = getLikelihoodArray(tInfo->rNumber,  tr->mxtips, tr->xVector);	  
	      ex2      = getScalingArray(tInfo->rNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber,  tr->mxtips, tr->xVector);
	      ex3      = getScalingArray(tInfo->pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
	      
	      break;
	    case INNER_INNER:
#ifdef _MULTI_GENE
	      /*printf("INNER_INNER %d %d\n", tInfo->qNumber, tInfo->rNumber);*/
#endif	 
	      x1_start = getLikelihoodArray(tInfo->qNumber,  tr->mxtips, tr->xVector);
	      ex1      = getScalingArray(tInfo->qNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);	 
	      
	      x2_start = getLikelihoodArray(tInfo->rNumber,  tr->mxtips, tr->xVector);
	      ex2      = getScalingArray(tInfo->rNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);	  
	      
	      x3_start = getLikelihoodArray(tInfo->pNumber,  tr->mxtips, tr->xVector);
	      ex3      = getScalingArray(tInfo->pNumber, tr->cdta->endsite, tr->mxtips, tr->expArray);
	      
	      break;
	    default:
	      assert(0);
	    }
      
      
	  switch(tr->likelihoodFunction)
	    {
	    case GTRCAT:
	    case GTRCATMULT:
	      newviewGTRCAT(tInfo,  &(tr->EV_DNA[model * 16]), &(tr->EI_DNA[model * 12]), &(tr->EIGN_DNA[model * 3]), 
			    tr->cdta->patrat, tr->cdta->rateCategory, 
			    x1_start, x2_start, x3_start, &(tr->tipVectorDNA[model * 64]),
			    ex1, ex2, ex3, tipX1, tipX2,
			    lower, upper, tr->NumberOfCategories, tInfo->qz[branchIndex], tInfo->rz[branchIndex]
			    );   
	      break;
	    case GTRGAMMA:  /* needed for rate opt*/
	    case GTRGAMMAI: /* needed for rate opt*/
	    case GTRGAMMAMULT:
	    case GTRGAMMAMULTI:
	      newviewGTRGAMMA(tInfo,
			      x1_start, x2_start, x3_start,
			      &(tr->EIGN_DNA[model * 3]), &(tr->EV_DNA[model * 16]), 
			      &(tr->EI_DNA[model * 12]), &(tr->gammaRates[model * 4]), 
			      &(tr->tipVectorDNA[model * 64]),
			      ex1, ex2, ex3, tipX1, tipX2,
			      lower, upper, tInfo->qz[branchIndex], tInfo->rz[branchIndex]);
	      break;
	    case PROTGAMMA:  /* needed for rate opt*/
	    case PROTGAMMAI: /* needed for rate opt*/
	    case PROTGAMMAMULT:
	    case PROTGAMMAMULTI:
	      newviewGTRGAMMAPROT(tInfo,
				  x1_start, x2_start, x3_start,
				  &(tr->EIGN_AA[model * 19]), &(tr->EV_AA[model * 400]), 
				  &(tr->EI_AA[model * 380]), &(tr->gammaRates[model * 4]), 
				  &(tr->tipVectorAA[model * 460]),
				  ex1, ex2, ex3, tipX1, tipX2,
				  lower, upper, tInfo->qz[branchIndex], tInfo->rz[branchIndex]);     
	      break;	 
	    default:
	      assert(0);
	    }
	}
    }
}

#endif



void newviewPartitionGeneric (tree *tr, nodeptr p, int model)
{
  if(isTip(p->number, tr->rdta->numsp)) 
    return;
  
  { 	          
   
#ifndef _USE_PTHREADS
    int lower = tr->partitionData[model].lower;
    int upper = tr->partitionData[model].upper;
#endif
   
    tr->td[0].count = 1;
    tr->modelNumber    = model;
    computeTraversalInfo(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp, tr->numBranches);   

    if(tr->td[0].count > 1)
      {
#ifdef _USE_PTHREADS  
	masterBarrier(THREAD_NEWVIEW_PARTITION, tr);  	 
#else
	newviewIterativePartition(tr, lower, upper, model);   
#endif
      }  
  }
}



