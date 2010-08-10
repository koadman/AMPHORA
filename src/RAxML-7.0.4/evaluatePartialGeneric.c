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





/********************** GTRCAT ***************************************/




#ifdef WIN32
static void computeVectorGTRCATPROT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
				    traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
				    char **yVector, int mxtips)
#else
static inline void computeVectorGTRCATPROT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
					   traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
					   char **yVector, int mxtips)
#endif
{       
  double   *x1, *x2, *x3;  
  int ex3,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &(lVector[20 * (pNumber  - mxtips)]);  
  ex3 = pNumber - mxtips;   

  switch(ti->tipCase)
    {
    case TIP_TIP:    
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(tipVector[20 * yVector[rNumber][i]]);
      eVector[ex3] = 0;
      break;
    case TIP_INNER:     
      x1 = &(tipVector[20 * yVector[qNumber][i]]);
      x2 = &(  lVector[20 * (rNumber - mxtips)]);     
      eVector[ex3] = eVector[rNumber - mxtips];           
      break;
    case INNER_INNER:            
      x1 = &(lVector[20 * (qNumber - mxtips)]);
      x2 = &(lVector[20 * (rNumber - mxtips)]);
      eVector[ex3] = eVector[qNumber - mxtips] + eVector[rNumber - mxtips];            
      break;    
    default:
      assert(0);
    }
     
  {
    double  d1[19], d2[19], ump_x1[20], ump_x2[20], x1px2, lz1, lz2;
    double *left, *right, *eptr, *eptr2;
    int l, scale;
     
    lz1 = qz * ki;            
    lz2 = rz * ki;
    
    left  = x1;
    right = x2;
    
    d1[0] = left[1] * exp(EIGN[0] * lz1);
    d1[1] = left[2] * exp(EIGN[1] * lz1);
    d1[2] = left[3] * exp(EIGN[2] * lz1);
    d1[3] = left[4] * exp(EIGN[3] * lz1);
    d1[4] = left[5] * exp(EIGN[4] * lz1);
    d1[5] = left[6] * exp(EIGN[5] * lz1);
    d1[6] = left[7] * exp(EIGN[6] * lz1);
    d1[7] = left[8] * exp(EIGN[7] * lz1);
    d1[8] = left[9] * exp(EIGN[8] * lz1);
    d1[9] = left[10] * exp(EIGN[9] * lz1);
    d1[10] = left[11] * exp(EIGN[10] * lz1);
    d1[11] = left[12] * exp(EIGN[11] * lz1);
    d1[12] = left[13] * exp(EIGN[12] * lz1);
    d1[13] = left[14] * exp(EIGN[13] * lz1);
    d1[14] = left[15] * exp(EIGN[14] * lz1);
    d1[15] = left[16] * exp(EIGN[15] * lz1);
    d1[16] = left[17] * exp(EIGN[16] * lz1);
    d1[17] = left[18] * exp(EIGN[17] * lz1);
    d1[18] = left[19] * exp(EIGN[18] * lz1);
    
    d2[0] = right[1] * exp(EIGN[0] * lz2);
    d2[1] = right[2] * exp(EIGN[1] * lz2);
    d2[2] = right[3] * exp(EIGN[2] * lz2);
    d2[3] = right[4] * exp(EIGN[3] * lz2);
    d2[4] = right[5] * exp(EIGN[4] * lz2);
    d2[5] = right[6] * exp(EIGN[5] * lz2);
    d2[6] = right[7] * exp(EIGN[6] * lz2);
    d2[7] = right[8] * exp(EIGN[7] * lz2);
    d2[8] = right[9] * exp(EIGN[8] * lz2);
    d2[9] = right[10] * exp(EIGN[9] * lz2);
    d2[10] = right[11] * exp(EIGN[10] * lz2);
    d2[11] = right[12] * exp(EIGN[11] * lz2);
    d2[12] = right[13] * exp(EIGN[12] * lz2);
    d2[13] = right[14] * exp(EIGN[13] * lz2);
    d2[14] = right[15] * exp(EIGN[14] * lz2);
    d2[15] = right[16] * exp(EIGN[15] * lz2);
    d2[16] = right[17] * exp(EIGN[16] * lz2);
    d2[17] = right[18] * exp(EIGN[17] * lz2);
    d2[18] = right[19] * exp(EIGN[18] * lz2);
      
    eptr = EI;
    for(l = 0; l < 20; l++)
      {
	ump_x1[l] = left[0];    
	ump_x1[l] += d1[0] * *eptr++;
	ump_x1[l] += d1[1] * *eptr++;
	ump_x1[l] += d1[2] * *eptr++;
	ump_x1[l] += d1[3] * *eptr++;
	ump_x1[l] += d1[4] * *eptr++;
	ump_x1[l] += d1[5] * *eptr++;
	ump_x1[l] += d1[6] * *eptr++;
	ump_x1[l] += d1[7] * *eptr++;
	ump_x1[l] += d1[8] * *eptr++;
	ump_x1[l] += d1[9] * *eptr++;
	ump_x1[l] += d1[10] * *eptr++;
	ump_x1[l] += d1[11] * *eptr++;
	ump_x1[l] += d1[12] * *eptr++;
	ump_x1[l] += d1[13] * *eptr++;
	ump_x1[l] += d1[14] * *eptr++;
	ump_x1[l] += d1[15] * *eptr++;
	ump_x1[l] += d1[16] * *eptr++;
	ump_x1[l] += d1[17] * *eptr++;
	ump_x1[l] += d1[18] * *eptr++;   
      }
      
    eptr = EI;
    for(l = 0; l < 20; l++)
      {
	ump_x2[l] = right[0];
	ump_x2[l] += d2[0] * *eptr++;
	ump_x2[l] += d2[1] * *eptr++;
	ump_x2[l] += d2[2] * *eptr++;
	ump_x2[l] += d2[3] * *eptr++;
	ump_x2[l] += d2[4] * *eptr++;
	ump_x2[l] += d2[5] * *eptr++;
	ump_x2[l] += d2[6] * *eptr++;
	ump_x2[l] += d2[7] * *eptr++;
	ump_x2[l] += d2[8] * *eptr++;
	ump_x2[l] += d2[9] * *eptr++;
	ump_x2[l] += d2[10] * *eptr++;
	ump_x2[l] += d2[11] * *eptr++;
	ump_x2[l] += d2[12] * *eptr++;
	ump_x2[l] += d2[13] * *eptr++;
	ump_x2[l] += d2[14] * *eptr++;
	ump_x2[l] += d2[15] * *eptr++;
	ump_x2[l] += d2[16] * *eptr++;
	ump_x2[l] += d2[17] * *eptr++;
	ump_x2[l] += d2[18] * *eptr++;   
      }
    
    left = x3;
    eptr2 = EV;

    x1px2 = ump_x1[0] * ump_x2[0];  
    left[0] = x1px2 * *eptr2++;
    left[1] = x1px2 * *eptr2++;
    left[2] = x1px2 * *eptr2++;
    left[3] = x1px2 * *eptr2++;
    left[4] = x1px2 * *eptr2++;
    left[5] = x1px2 * *eptr2++;
    left[6] = x1px2 * *eptr2++;
    left[7] = x1px2 * *eptr2++;
    left[8] = x1px2 * *eptr2++;
    left[9] = x1px2 * *eptr2++;
    left[10] = x1px2 * *eptr2++;
    left[11] = x1px2 * *eptr2++;
    left[12] = x1px2 * *eptr2++;
    left[13] = x1px2 * *eptr2++;
    left[14] = x1px2 * *eptr2++;
    left[15] = x1px2 * *eptr2++;
    left[16] = x1px2 * *eptr2++;
    left[17] = x1px2 * *eptr2++;
    left[18] = x1px2 * *eptr2++;
    left[19] = x1px2 * *eptr2++;
    
    for(l = 1; l < 20; l++)
      {
	x1px2 = ump_x1[l] * ump_x2[l];
	left[0] += x1px2 * *eptr2++;
	left[1] += x1px2 * *eptr2++;
	left[2] += x1px2 * *eptr2++;
	left[3] += x1px2 * *eptr2++;
	left[4] += x1px2 * *eptr2++;
	left[5] += x1px2 * *eptr2++;
	left[6] += x1px2 * *eptr2++;
	left[7] += x1px2 * *eptr2++;
	left[8] += x1px2 * *eptr2++;
	left[9] += x1px2 * *eptr2++;
	left[10] += x1px2 * *eptr2++;
	left[11] += x1px2 * *eptr2++;
	left[12] += x1px2 * *eptr2++;
	left[13] += x1px2 * *eptr2++;
	left[14] += x1px2 * *eptr2++;
	left[15] += x1px2 * *eptr2++;
	left[16] += x1px2 * *eptr2++;
	left[17] += x1px2 * *eptr2++;
	left[18] += x1px2 * *eptr2++;
	left[19] += x1px2 * *eptr2++;
      }
    
    scale = 1;
    for(l = 0; scale && (l < 20); l++)
      scale = ((left[l] < minlikelihood) && (left[l] > minusminlikelihood));	       	      	      	       	       
    
    if(scale)
      {	
	for(l = 0; l < 20; l++)
	  left[l] *= twotothe256;		   
	eVector[ex3] += 1;
      }
    
    return;      
  }
}

#ifdef WIN32
static void computeVectorGTRCAT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
				traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
				char **yVector, int mxtips)
#else
static inline void computeVectorGTRCAT(double *lVector, int *eVector, double ki, int i, double qz, double rz,
				       traversalInfo *ti, double *EIGN, double *EI, double *EV, double *tipVector, 
				       char **yVector, int mxtips)
#endif
{       
  double  d1[3], d2[3],  ump_x1_1, ump_x1_2, ump_x1_3, ump_x1_0, 
    ump_x2_0, ump_x2_1, ump_x2_2, ump_x2_3, x1px2, lz1, lz2; 
  double *x1, *x2, *x3;
  int ex3,
    pNumber = ti->pNumber,
    rNumber = ti->rNumber,
    qNumber = ti->qNumber;
 
  x3  = &lVector[4 * (pNumber  - mxtips)];  
  ex3 = pNumber - mxtips;   

  switch(ti->tipCase)
    {
    case TIP_TIP:     
      x1 = &(tipVector[4 * yVector[qNumber][i]]);
      x2 = &(tipVector[4 * yVector[rNumber][i]]);
      eVector[ex3] = 0;
      break;
    case TIP_INNER:     
      x1 = &(tipVector[4 * yVector[qNumber][i]]);
      x2 = &lVector[4 * (rNumber - mxtips)];     
      eVector[ex3] = eVector[rNumber - mxtips];           
      break;
    case INNER_INNER:            
      x1 = &lVector[4 * (qNumber - mxtips)];
      x2 = &lVector[4 * (rNumber - mxtips)];
      eVector[ex3] = eVector[qNumber - mxtips] + eVector[rNumber - mxtips];            
      break;
    default:
      assert(0);
    }
     
  lz1 = qz * ki;  
  lz2 = rz * ki;
  
  d1[0] = x1[1] * exp(EIGN[0] * lz1);
  d2[0] = x2[1] * exp(EIGN[0] * lz2);	    
  d1[1] = x1[2] * exp(EIGN[1] * lz1);
  d2[1] = x2[2] * exp(EIGN[1] * lz2);
  d1[2] = x1[3] * exp(EIGN[2] * lz1);
  d2[2] = x2[3] * exp(EIGN[2] * lz2); 
  
  ump_x1_0  = d1[0] * EI[0];
  ump_x1_0 += d1[1] * EI[1];
  ump_x1_0 += d1[2] * EI[2];	      	
  ump_x1_0 += x1[0];
  
  ump_x1_1  = d1[0] * EI[3];
  ump_x1_1 += d1[1] * EI[4];
  ump_x1_1 += d1[2] * EI[5];	      	
  ump_x1_1 += x1[0];
  
  ump_x1_2  = d1[0] * EI[6];
  ump_x1_2 += d1[1] * EI[7];
  ump_x1_2 += d1[2] * EI[8];	      	
  ump_x1_2 += x1[0];
  
  ump_x1_3  = d1[0] * EI[9];
  ump_x1_3 += d1[1] * EI[10];
  ump_x1_3 += d1[2] * EI[11];	      	
  ump_x1_3 += x1[0]; 
  	 	    	    	     	     
  ump_x2_0  = d2[0] * EI[0];
  ump_x2_0 += d2[1] * EI[1];
  ump_x2_0 += d2[2] * EI[2];		     
  ump_x2_0 += x2[0];
  
  ump_x2_1  = d2[0] * EI[3];
  ump_x2_1 += d2[1] * EI[4];
  ump_x2_1 += d2[2] * EI[5];		     
  ump_x2_1 += x2[0];	 
  
  ump_x2_2  = d2[0] * EI[6];
  ump_x2_2 += d2[1] * EI[7];
  ump_x2_2 += d2[2] * EI[8];		     
  ump_x2_2 += x2[0];	  
		   
  ump_x2_3  = d2[0] * EI[9];
  ump_x2_3 += d2[1] * EI[10];
  ump_x2_3 += d2[2] * EI[11];		     
  ump_x2_3 += x2[0];	    	  		   	  
    
  x1px2 = ump_x1_0 * ump_x2_0;
  x3[0] = x1px2 *  EV[0];
  x3[1] = x1px2 *  EV[1];
  x3[2] = x1px2 *  EV[2];
  x3[3] = x1px2 *  EV[3];
  
  x1px2 = ump_x1_1 * ump_x2_1;
  x3[0] += x1px2  *  EV[4];
  x3[1] += x1px2 *   EV[5];
  x3[2] += x1px2 *   EV[6];
  x3[3] += x1px2 *   EV[7];
  
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
       
  if (x3[0] < minlikelihood && x3[0] > minusminlikelihood &&
      x3[1] < minlikelihood && x3[1] > minusminlikelihood &&
      x3[2] < minlikelihood && x3[2] > minusminlikelihood &&
      x3[3] < minlikelihood && x3[3] > minusminlikelihood)
    {	     
      x3[0]   *= twotothe256;
      x3[1]   *= twotothe256;
      x3[2]   *= twotothe256;     
      x3[3]   *= twotothe256;
      eVector[ex3] += 1;
    }	              

  return;
}


static double evaluatePartialGTRCAT(int i, double ki, int counter,  traversalInfo *ti, double qz,
				    int w, double *EIGN, double *EI, double *EV,
				    double *tipVector, char **yVector, 
				    int branchReference, int mxtips)
{
  double lz, term;       
  double  d[3];
  double   *x1, *x2; 
  int scale, k;
  double *lVector = (double *)malloc(sizeof(double) * 4 * mxtips); 
  int *eVector    = (int *)malloc(sizeof(int) * mxtips);
  traversalInfo *trav = &ti[0];
 
  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[4 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    computeVectorGTRCAT(lVector, eVector, ki, i, ti[k].qz[branchReference], ti[k].rz[branchReference], &ti[k], 
			EIGN, EI, EV, 
			tipVector, yVector, mxtips);       
   
  x2 = &lVector[4 * (trav->qNumber - mxtips)];

  scale = eVector[trav->qNumber - mxtips];      

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
       
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;  
  
  d[0] = exp (EIGN[0] * lz);
  d[1] = exp (EIGN[1] * lz);
  d[2] = exp (EIGN[2] * lz);       	   
  
  term =  x1[0] * x2[0];
  term += x1[1] * x2[1] * d[0];
  term += x1[2] * x2[2] * d[1];
  term += x1[3] * x2[3] * d[2];     

  term = log(term) + (scale * log(minlikelihood));   

  term = term * w;

  free(lVector);
  free(eVector);

  return  term;
}



static double evaluatePartialGTRCATPROT(int i, double ki, int counter,  traversalInfo *ti, double qz,
					int w, double *EIGN, double *EI, double *EV,
					double *tipVector, char **yVector, 
					int branchReference, int mxtips)
{
  double lz, term, *left, *right;       
  double  *ds, d[19];
  double   *x1, *x2; 
  int scale, k;
  double *lVector = (double *)malloc(sizeof(double) * 20 * mxtips);
  int *eVector    = (int *)malloc(sizeof(int) * mxtips); 
  traversalInfo *trav = &ti[0];

  ds = d;

  assert(isTip(trav->pNumber, mxtips));
     
  x1 = &(tipVector[20 *  yVector[trav->pNumber][i]]);   

  for(k = 1; k < counter; k++)                
    computeVectorGTRCATPROT(lVector, eVector, ki, i, ti[k].qz[branchReference], ti[k].rz[branchReference], 
			    &ti[k], EIGN, EI, EV, 
			    tipVector, yVector, mxtips);       
   
  x2 = &lVector[20 * (trav->qNumber - mxtips)];

  scale = eVector[trav->qNumber - mxtips];      

  assert(0 <=  (trav->qNumber - mxtips) && (trav->qNumber - mxtips) < mxtips);  
  
  if(qz < zmin) 
    lz = zmin;
  lz  = log(qz); 
  lz *= ki;
  
  d[0] = exp (EIGN[0] * lz);
  d[1] = exp (EIGN[1] * lz);
  d[2] = exp (EIGN[2] * lz);
  d[3] = exp (EIGN[3] * lz);
  d[4] = exp (EIGN[4] * lz);
  d[5] = exp (EIGN[5] * lz);
  d[6] = exp (EIGN[6] * lz);
  d[7] = exp (EIGN[7] * lz);
  d[8] = exp (EIGN[8] * lz);
  d[9] = exp (EIGN[9] * lz);
  d[10] = exp (EIGN[10] * lz);
  d[11] = exp (EIGN[11] * lz);
  d[12] = exp (EIGN[12] * lz);
  d[13] = exp (EIGN[13] * lz);
  d[14] = exp (EIGN[14] * lz);
  d[15] = exp (EIGN[15] * lz);
  d[16] = exp (EIGN[16] * lz);
  d[17] = exp (EIGN[17] * lz);
  d[18] = exp (EIGN[18] * lz);
  
  
  left  = x1;
  right = x2;
  
  term =  *left++ * *right++;  
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++;
  
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++;
  term += *left++ * *right++ * *ds++; 
  
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++;
  term += *left++ * *right++ * *ds++; 
  
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++; 
  term += *left++ * *right++ * *ds++;
  term += *left++ * *right++ * *ds++;   

  term = log(term) + (scale * log(minlikelihood));   

  term = term * w;

  free(lVector);
  free(eVector);

  return  term;
}






/*********************************************************************************************/



void computeFullTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches)
{
  if(isTip(p->number, maxTips))
    return; 

  {     
    int i;
    nodeptr q = p->next->back;
    nodeptr r = p->next->next->back;

    /* set xnode info at this point */

    p->x = 1;
    p->next->x = 0;
    p->next->next->x = 0;     

    if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
      {	  
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
	    
	    computeFullTraversalInfo(r, ti, counter, maxTips, numBranches);	
	    	   
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
	    computeFullTraversalInfo(q, ti, counter, maxTips, numBranches);	       
	    computeFullTraversalInfo(r, ti, counter, maxTips, numBranches);
	   
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

void determineFullTraversal(nodeptr p, tree *tr)
{
  nodeptr q = p->back;
  int k;

  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  for(k = 0; k < tr->numBranches; k++)        
    tr->td[0].ti[0].qz[k] = q->z[k];    

  assert(isTip(p->number, tr->rdta->numsp));

  tr->td[0].count = 1; 
  computeFullTraversalInfo(q, &(tr->td[0].ti[0]),  &(tr->td[0].count), tr->rdta->numsp, tr->numBranches); 
  computeFullTraversalInfo(p, &(tr->td[0].ti[0]),  &(tr->td[0].count), tr->rdta->numsp, tr->numBranches);


}


#ifdef _MULTI_GENE

void computeFullMultiTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int model)
{
  if(isTip(p->number, maxTips))
    return; 

  {        
    /* set xnode info at this point */

    /*p->x = 1;
    p->next->x = 0;
    p->next->next->x = 0;       */

    if(p->backs[model])
      {
	nodeptr q = p->next->backs[model];
	nodeptr r = p->next->next->backs[model];
	assert(p == p->next->next->next);
	p->xs[model] = 1;
	p->next->xs[model] = 0;
	p->next->next->xs[model] = 0;
	
	if(isTip(r->number, maxTips) && isTip(q->number, maxTips))
	  {	  
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
		
		computeFullMultiTraversalInfo(r, ti, counter, maxTips, model);	
		
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
		computeFullMultiTraversalInfo(q, ti, counter, maxTips, model);	       
		computeFullMultiTraversalInfo(r, ti, counter, maxTips, model);
		
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
    else
      {	
	p->xs[model] = 0;
	p->next->xs[model] = 0;
	p->next->next->xs[model] = 0;
	assert(p == p->next->next->next);

	computeFullMultiTraversalInfo(p->next->back, ti, counter, maxTips, model);
	computeFullMultiTraversalInfo(p->next->next->back, ti, counter, maxTips, model);
      }
  }
}


void determineFullMultiTraversal(tree *tr)
{
  int model;

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      nodeptr p = tr->startVector[model];

      nodeptr q = p->backs[model];     
      
      tr->td[model].ti[0].pNumber = p->number;
      tr->td[model].ti[0].qNumber = q->number;
           
      tr->td[model].ti[0].qz[model] = q->z[model];    
      
      assert(isTip(p->number, tr->rdta->numsp));           

      tr->td[model].count = 1;

      computeFullMultiTraversalInfo(q, &(tr->td[model].ti[0]),  &(tr->td[model].count), tr->rdta->numsp, model); 
      computeFullMultiTraversalInfo(p, &(tr->td[model].ti[0]),  &(tr->td[model].count), tr->rdta->numsp, model);

      /*printf("%d %d\n", model,  tr->td[model].count);*/
    }
}
#endif

#ifdef _LOCAL_DATA

double evaluatePartialGeneric (tree *localTree, int i, double ki)
{
  double result;
  int model, branchReference;   
 
  model = localTree->strided_model[i];
  
  if(localTree->multiBranch)
    branchReference = model;
  else
    branchReference = 0;

  if(localTree->mixedData)
    {
      assert(0);

      /*switch(tr->partitionData[model].dataType)
	{
	case DNA_DATA:
	result = evaluatePartialGTRCAT(i, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
	tr->cdta->aliaswgt[i],
	&(tr->EIGN_DNA[model * 3]), &(tr->EI_DNA[model * 12]), &(tr->EV_DNA[model * 16]),
	&(tr->tipVectorDNA[model * 64]),
	tr->yVector, branchReference, tr->mxtips);
	break;
	case AA_DATA:
	result = evaluatePartialGTRCATPROT(i, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
	tr->cdta->aliaswgt[i],
	&(tr->EIGN_AA[model * 19]), &(tr->EI_AA[model * 380]), 
	&(tr->EV_AA[model * 400]),
	&(tr->tipVectorAA[model * 460]), 
	tr->yVector, branchReference, tr->mxtips);
	break;
	default:
	assert(0);
	}
      */
    }
  else
    {
    
      switch(localTree->likelihoodFunction)
	{
	case GTRCAT:
	  result = evaluatePartialGTRCAT(i, ki, localTree->td[0].count, localTree->td[0].ti, localTree->td[0].ti[0].qz[0], 
					 localTree->strided_aliaswgt[i],
					 &(localTree->EIGN_DNA[0]), &(localTree->EI_DNA[0]), &(localTree->EV_DNA[0]),
					 &(localTree->tipVectorDNA[0]), 
					 localTree->strided_yVector, 0, localTree->mxtips);     
	  break;
	case GTRCATMULT:
	  result = evaluatePartialGTRCAT(i, ki, localTree->td[0].count, localTree->td[0].ti, localTree->td[0].ti[0].qz[branchReference], 
					 localTree->strided_aliaswgt[i],
					 &(localTree->EIGN_DNA[model * 3]), &(localTree->EI_DNA[model * 12]), 
					 &(localTree->EV_DNA[model * 16]),
					 &(localTree->tipVectorDNA[model * 64]), 
					 localTree->strided_yVector, branchReference, localTree->mxtips);
	  break;     
	case PROTCAT:
	  result = evaluatePartialGTRCATPROT(i, ki, localTree->td[0].count, localTree->td[0].ti, 
					     localTree->td[0].ti[0].qz[0], localTree->strided_aliaswgt[i],
					     &(localTree->EIGN_AA[0]), &(localTree->EI_AA[0]), &(localTree->EV_AA[0]),
					     &(localTree->tipVectorAA[0]), 
					     localTree->strided_yVector, 0, localTree->mxtips);     
	  break;
	case PROTCATMULT:
	  result = evaluatePartialGTRCATPROT(i, ki, localTree->td[0].count, localTree->td[0].ti, localTree->td[0].ti[0].qz[branchReference], 
					     localTree->strided_aliaswgt[i],
					     &(localTree->EIGN_AA[model * 19]), &(localTree->EI_AA[model * 380]), 
					     &(localTree->EV_AA[model * 400]),
					     &(localTree->tipVectorAA[model * 460]), 
					     localTree->strided_yVector, branchReference, localTree->mxtips);  
	  break;
	default:
	  assert(0);
	}
    }

  return result;
}

#else

double evaluatePartialGeneric (tree *tr, int i, double ki)
{
  double result;
  int model, branchReference;   
 
  model = tr->model[i];
  
  if(tr->multiBranch)
    branchReference = model;
  else
    branchReference = 0;

  if(tr->mixedData)
    {
    

      switch(tr->partitionData[model].dataType)
	{
	case DNA_DATA:
	  result = evaluatePartialGTRCAT(i, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					 tr->cdta->aliaswgt[i],
					 &(tr->EIGN_DNA[model * 3]), &(tr->EI_DNA[model * 12]), &(tr->EV_DNA[model * 16]),
					 &(tr->tipVectorDNA[model * 64]),
					 tr->yVector, branchReference, tr->mxtips);
	  break;
	case AA_DATA:
	  result = evaluatePartialGTRCATPROT(i, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					     tr->cdta->aliaswgt[i],
					     &(tr->EIGN_AA[model * 19]), &(tr->EI_AA[model * 380]), 
					     &(tr->EV_AA[model * 400]),
					     &(tr->tipVectorAA[model * 460]), 
					     tr->yVector, branchReference, tr->mxtips);
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
	  result = evaluatePartialGTRCAT(i, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[0], tr->cdta->aliaswgt[i],
					 &(tr->EIGN_DNA[0]), &(tr->EI_DNA[0]), &(tr->EV_DNA[0]),
					 &(tr->tipVectorDNA[0]), 
					 tr->yVector, 0, tr->mxtips);     
	  break;
	case GTRCATMULT:
	  result = evaluatePartialGTRCAT(i, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					 tr->cdta->aliaswgt[i],
					 &(tr->EIGN_DNA[model * 3]), &(tr->EI_DNA[model * 12]), &(tr->EV_DNA[model * 16]),
					 &(tr->tipVectorDNA[model * 64]), 
					 tr->yVector, branchReference, tr->mxtips);
	  break;     
	case PROTCAT:
	  result = evaluatePartialGTRCATPROT(i, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[0], tr->cdta->aliaswgt[i],
					     &(tr->EIGN_AA[0]), &(tr->EI_AA[0]), &(tr->EV_AA[0]),
					     &(tr->tipVectorAA[0]), 
					     tr->yVector, 0, tr->mxtips);     
	  break;
	case PROTCATMULT:
	  result = evaluatePartialGTRCATPROT(i, ki, tr->td[0].count, tr->td[0].ti, tr->td[0].ti[0].qz[branchReference], 
					     tr->cdta->aliaswgt[i],
					     &(tr->EIGN_AA[model * 19]), &(tr->EI_AA[model * 380]), &(tr->EV_AA[model * 400]),
					     &(tr->tipVectorAA[model * 460]), 
					     tr->yVector, branchReference, tr->mxtips);  
	  break;
	default:
	  assert(0);
	}
    }

  return result;
}

#endif
