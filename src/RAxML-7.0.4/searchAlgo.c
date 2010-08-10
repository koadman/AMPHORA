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

#ifndef WIN32
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h> 
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>



#include "axml.h"

extern int Thorough;
extern infoList iList;
extern char inverseMeaningDNA[16];
extern char seq_file[1024];
extern char resultFileName[1024];
extern char tree_file[1024];
extern FILE *INFILE;

boolean initrav (tree *tr, nodeptr p)
{ 
  nodeptr  q;
  
  if (!isTip(p->number, tr->rdta->numsp)) 
    {      
      q = p->next;
      
      do 
	{	   
	  if (! initrav(tr, q->back))  return FALSE;		   
	  q = q->next;	
	} 
      while (q != p);  
      
      newviewGeneric(tr, p);
    }
  
  return TRUE;
} 





boolean initravDIST (tree *tr, nodeptr p, int distance)
  {
    nodeptr  q;

    if (/*! p->tip*/ !isTip(p->number, tr->rdta->numsp) && distance > 0) 
      {
      
	q = p->next;
      
	do 
	  {
	    if (! initravDIST(tr, q->back, --distance))  return FALSE;	
	    q = q->next;	
	  } 
	while (q != p);
      
      
	newviewGeneric(tr, p);
      }

    return TRUE;
  } /* initrav */

void initravPartition (tree *tr, nodeptr p, int model)
{
  nodeptr  q;
  
  if (/*!p->tip*/ !isTip(p->number, tr->rdta->numsp)) 
    {      
      q = p->next;      
      do 
	{
	  initravPartition(tr, q->back, model);
	  q = q->next;	
	} 
      while (q != p);
    
      newviewPartitionGeneric(tr, p, model);       
    }
} 





boolean update(tree *tr, nodeptr p)
{       
  nodeptr  q;
  boolean smoothed;
  int i;
  double   z[NUM_BRANCHES], z0[NUM_BRANCHES];
  
  q = p->back;   

  for(i = 0; i < tr->numBranches; i++)
    z0[i] = q->z[i];

  makenewzGeneric(tr, p, q, z0, newzpercycle, z);  
  
  smoothed = tr->smoothed;

  for(i = 0; i < tr->numBranches; i++)
    {      
      if(ABS(z[i] - z0[i]) > deltaz)  
	smoothed = FALSE;
      p->z[i] = q->z[i] = z[i];      
    }
  
  tr->smoothed = smoothed;
  
  return TRUE;
}




boolean smooth (tree *tr, nodeptr p)
{
  nodeptr  q;
  
  if (! update(tr, p))               return FALSE; /*  Adjust branch */
  if (! isTip(p->number, tr->rdta->numsp)) 
    {                                  /*  Adjust descendants */
      q = p->next;
      while (q != p) 
	{
	  if (! smooth(tr, q->back))   return FALSE;
	  q = q->next;
	}	
      

      newviewGeneric(tr, p);
    }
  
  return TRUE;
} 


boolean smoothTree (tree *tr, int maxtimes)
{
  nodeptr  p, q;   
   
  p = tr->start;
  
  while (--maxtimes >= 0) 
    {
      tr->smoothed = TRUE;
      if (! smooth(tr, p->back))       return FALSE;
      if (!isTip(p->number, tr->rdta->numsp)) 
	{
	  q = p->next;
	  while (q != p) 
	    {
	      if (! smooth(tr, q->back))   return FALSE;
	      q = q->next;
	    }
	}
      if (tr->smoothed)  break;
    }
  
  return TRUE;
} 



boolean localSmooth (tree *tr, nodeptr p, int maxtimes)
{ 
  nodeptr  q;
  
  if (/*p->tip*/ isTip(p->number, tr->rdta->numsp)) return FALSE;            /* Should be an error */
  
  while (--maxtimes >= 0) 
    {
      tr->smoothed = TRUE;
      q = p;
      do 
	{
	  if (! update(tr, q)) return FALSE;
	  q = q->next;
        } 
      while (q != p);
      if (tr->smoothed)  break;
    }

  tr->smoothed = FALSE;             /* Only smooth locally */
  return TRUE;
}






static void resetInfoList(void)
{
  int i;

  iList.valid = 0;

  for(i = 0; i < iList.n; i++)    
    {
      iList.list[i].node = (nodeptr)NULL;
      iList.list[i].likelihood = unlikely;
    }    
}

void initInfoList(int n)
{
  int i;

  iList.n = n;
  iList.valid = 0;
  iList.list = (bestInfo *)malloc(sizeof(bestInfo) * n);

  for(i = 0; i < n; i++)
    {
      iList.list[i].node = (nodeptr)NULL;
      iList.list[i].likelihood = unlikely;
    }
}

void freeInfoList(void)
{ 
  free(iList.list);   
}


void insertInfoList(nodeptr node, double likelihood)
{
  int i;
  int min = 0;
  double min_l =  iList.list[0].likelihood;

  for(i = 1; i < iList.n; i++)
    {
      if(iList.list[i].likelihood < min_l)
	{
	  min = i;
	  min_l = iList.list[i].likelihood;
	}
    }

  if(likelihood > min_l)
    {
      iList.list[min].likelihood = likelihood;
      iList.list[min].node = node;
      iList.valid += 1;
    }

  if(iList.valid > iList.n)
    iList.valid = iList.n;
}


boolean smoothRegion (tree *tr, nodeptr p, int region)
{ 
  nodeptr  q;
  
  if (! update(tr, p))               return FALSE; /*  Adjust branch */

  if(region > 0)
    {
      if (/*! p->tip*/  !isTip(p->number, tr->rdta->numsp)) 
	{                                 
	  q = p->next;
	  while (q != p) 
	    {
	      if (! smoothRegion(tr, q->back, --region))   return FALSE;
	      q = q->next;
	    }	
	  
	    newviewGeneric(tr, p);
	}
    }
  
  return TRUE;
}

boolean regionalSmooth (tree *tr, nodeptr p, int maxtimes, int region)
  {
    nodeptr  q;

    if (isTip(p->number, tr->rdta->numsp)) return FALSE;            /* Should be an error */

    while (--maxtimes >= 0) 
      {
	tr->smoothed = TRUE;
	q = p;
	do 
	  {
	    if (! smoothRegion(tr, q, region)) return FALSE;
	    q = q->next;
	  } 
	while (q != p);
	if (tr->smoothed)  break;
      }

    tr->smoothed = FALSE;             /* Only smooth locally */
    return TRUE;
  } /* localSmooth */





nodeptr  removeNodeBIG (tree *tr, nodeptr p, int numBranches)
{  
  double   zqr[NUM_BRANCHES], result[NUM_BRANCHES];
  nodeptr  q, r;
  int i;
        
  q = p->next->back;
  r = p->next->next->back;
  
  for(i = 0; i < numBranches; i++)
    zqr[i] = q->z[i] * r->z[i];        
   
  makenewzGeneric(tr, q, r, zqr, iterations, result);   

  for(i = 0; i < numBranches; i++)        
    tr->zqr[i] = result[i];

  hookup(q, r, result, numBranches); 
      
  p->next->next->back = p->next->back = (node *) NULL;

  return  q; 
}

nodeptr  removeNodeRestoreBIG (tree *tr, nodeptr p)
{
  nodeptr  q, r;
        
  q = p->next->back;
  r = p->next->next->back;  

  newviewGeneric(tr, q);
  newviewGeneric(tr, r);
  
  hookup(q, r, tr->currentZQR, tr->numBranches);

  p->next->next->back = p->next->back = (node *) NULL;
     
  return  q;
}


boolean insertBIG (tree *tr, nodeptr p, nodeptr q, int numBranches)
{
  nodeptr  r, s;
  int i;
  
  r = q->back;
  s = p->back;
      
  for(i = 0; i < numBranches; i++)
    tr->lzi[i] = q->z[i];
  
  if(Thorough)
    { 
      double  zqr[NUM_BRANCHES], zqs[NUM_BRANCHES], zrs[NUM_BRANCHES], lzqr, lzqs, lzrs, lzsum, lzq, lzr, lzs, lzmax;      
      double defaultArray[NUM_BRANCHES];	
      double e1[NUM_BRANCHES], e2[NUM_BRANCHES], e3[NUM_BRANCHES];
      double *qz;
      
      qz = q->z;
      
      for(i = 0; i < numBranches; i++)
	defaultArray[i] = defaultz;
      
      makenewzGeneric(tr, q, r, qz, iterations, zqr);           
      makenewzGeneric(tr, q, s, defaultArray, iterations, zqs);                  
      makenewzGeneric(tr, r, s, defaultArray, iterations, zrs);
      
      
      for(i = 0; i < numBranches; i++)
	{
	  lzqr = (zqr[i] > zmin) ? log(zqr[i]) : log(zmin); 
	  lzqs = (zqs[i] > zmin) ? log(zqs[i]) : log(zmin);
	  lzrs = (zrs[i] > zmin) ? log(zrs[i]) : log(zmin);
	  lzsum = 0.5 * (lzqr + lzqs + lzrs);
	  
	  lzq = lzsum - lzrs;
	  lzr = lzsum - lzqs;
	  lzs = lzsum - lzqr;
	  lzmax = log(zmax);
	  
	  if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
	  else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
	  else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          
	  
	  e1[i] = exp(lzq);
	  e2[i] = exp(lzr);
	  e3[i] = exp(lzs);
	}
      hookup(p->next,       q, e1, numBranches);
      hookup(p->next->next, r, e2, numBranches);
      hookup(p,             s, e3, numBranches);      		  
    }
  else
    {       
      double  z[NUM_BRANCHES]; 
      
      for(i = 0; i < numBranches; i++)
	{
	  z[i] = sqrt(q->z[i]);      
	  
	  if(z[i] < zmin) 
	    z[i] = zmin;
	  if(z[i] > zmax)
	    z[i] = zmax;
	}
      
      hookup(p->next,       q, z, tr->numBranches);
      hookup(p->next->next, r, z, tr->numBranches);	                         
    }
  
  newviewGeneric(tr, p);
  
  if(Thorough)
    {     
      localSmooth(tr, p, smoothings);   
      for(i = 0; i < numBranches; i++)
	{
	  tr->lzq[i] = p->next->z[i];
	  tr->lzr[i] = p->next->next->z[i];
	  tr->lzs[i] = p->z[i];            
	}
    }           
  
  return  TRUE;
}

boolean insertRestoreBIG (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r, s;
  
  r = q->back;
  s = p->back;

  if(Thorough)
    {                        
      hookup(p->next,       q, tr->currentLZQ, tr->numBranches);
      hookup(p->next->next, r, tr->currentLZR, tr->numBranches);
      hookup(p,             s, tr->currentLZS, tr->numBranches);      		  
    }
  else
    {       
      double  z[NUM_BRANCHES];
      int i;
      
      for(i = 0; i < tr->numBranches; i++)
	{
	  double zz;
	  zz = sqrt(q->z[i]);     
	  if(zz < zmin) 
	    zz = zmin;
	  if(zz > zmax)
	    zz = zmax;
  	  z[i] = zz;
	}

      hookup(p->next,       q, z, tr->numBranches);
      hookup(p->next->next, r, z, tr->numBranches);
    }   
    
  newviewGeneric(tr, p);
       
  return  TRUE;
}


void restoreTopologyOnly(tree *tr, bestlist *bt)
{ 
  nodeptr p = tr->removeNode;
  nodeptr q = tr->insertNode;
  double qz[NUM_BRANCHES], pz[NUM_BRANCHES], p1z[NUM_BRANCHES], p2z[NUM_BRANCHES];
  nodeptr p1, p2, r, s;
  double currentLH = tr->likelihood;
  int i;
      
  p1 = p->next->back;
  p2 = p->next->next->back;
  
  for(i = 0; i < tr->numBranches; i++)
    {
      p1z[i] = p1->z[i];
      p2z[i] = p2->z[i];
    }
  
  hookup(p1, p2, tr->currentZQR, tr->numBranches);
  
  p->next->next->back = p->next->back = (node *) NULL;             
  for(i = 0; i < tr->numBranches; i++)
    {
      qz[i] = q->z[i];
      pz[i] = p->z[i];           
    }
  
  r = q->back;
  s = p->back;
  
  if(Thorough)
    {                        
      hookup(p->next,       q, tr->currentLZQ, tr->numBranches);
      hookup(p->next->next, r, tr->currentLZR, tr->numBranches);
      hookup(p,             s, tr->currentLZS, tr->numBranches);      		  
    }
  else
    { 	
      double  z[NUM_BRANCHES];	
      for(i = 0; i < tr->numBranches; i++)
	{
	  z[i] = sqrt(q->z[i]);      
	  if(z[i] < zmin)
	    z[i] = zmin;
	  if(z[i] > zmax)
	    z[i] = zmax;
	}
      hookup(p->next,       q, z, tr->numBranches);
      hookup(p->next->next, r, z, tr->numBranches);
    }     
  
  tr->likelihood = tr->bestOfNode;
  
  saveBestTree(bt, tr);
  
  tr->likelihood = currentLH;
  
  hookup(q, r, qz, tr->numBranches);
  
  p->next->next->back = p->next->back = (nodeptr) NULL;
  
  if(Thorough)    
    hookup(p, s, pz, tr->numBranches);          
      
  hookup(p->next,       p1, p1z, tr->numBranches); 
  hookup(p->next->next, p2, p2z, tr->numBranches);      
}


boolean testInsertBIG (tree *tr, nodeptr p, nodeptr q)
{
  double  qz[NUM_BRANCHES], pz[NUM_BRANCHES];
  nodeptr  r;
  boolean doIt = TRUE;
  double startLH = tr->endLH;
  int i;
  
  r = q->back; 
  for(i = 0; i < tr->numBranches; i++)
    {
      qz[i] = q->z[i];
      pz[i] = p->z[i];
    }
  
  if(tr->grouped)
    {
      int rNumber, qNumber, pNumber;
      
      doIt = FALSE;
      
      rNumber = tr->constraintVector[r->number];
      qNumber = tr->constraintVector[q->number];
      pNumber = tr->constraintVector[p->number];
      
      if(pNumber == -9)
	pNumber = checker(tr, p->back);
      if(pNumber == -9)
	doIt = TRUE;
      else
	{
	  if(qNumber == -9)
	    qNumber = checker(tr, q);
	  
	  if(rNumber == -9)
	    rNumber = checker(tr, r);
	  
	  if(pNumber == rNumber || pNumber == qNumber)
	    doIt = TRUE;       
	}
    }
  
  if(doIt)
    {     
      if (! insertBIG(tr, p, q, tr->numBranches))       return FALSE;         
      
      evaluateGeneric(tr, p->next->next);       
      		   
      if(tr->likelihood > tr->bestOfNode)
	{
	  tr->bestOfNode = tr->likelihood;
	  tr->insertNode = q;
	  tr->removeNode = p;   
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      tr->currentZQR[i] = tr->zqr[i];           
	      tr->currentLZR[i] = tr->lzr[i];
	      tr->currentLZQ[i] = tr->lzq[i];
	      tr->currentLZS[i] = tr->lzs[i];      
	    }
	}
      
      if(tr->likelihood > tr->endLH)
	{			  
	  tr->insertNode = q;
	  tr->removeNode = p;   
	  for(i = 0; i < tr->numBranches; i++)
	    tr->currentZQR[i] = tr->zqr[i];      
	  tr->endLH = tr->likelihood;                      
	}        
      
      hookup(q, r, qz, tr->numBranches);
      
      p->next->next->back = p->next->back = (nodeptr) NULL;
      
      if(Thorough)
	{
	  nodeptr s = p->back;
	  hookup(p, s, pz, tr->numBranches);      
	} 
      
      if((tr->doCutoff) && (tr->likelihood < startLH))
	{
	  tr->lhAVG += (startLH - tr->likelihood);
	  tr->lhDEC++;
	  if((startLH - tr->likelihood) >= tr->lhCutoff)
	    return FALSE;	    
	  else
	    return TRUE;
	}
      else
	return TRUE;
    }
  else
    return TRUE;  
}





 
void addTraverseBIG(tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav)
{  
  if (--mintrav <= 0) 
    {              
      if (! testInsertBIG(tr, p, q))  return;        
    }
  
  if ((!isTip(q->number, tr->rdta->numsp)) && (--maxtrav > 0)) 
    {    
      addTraverseBIG(tr, p, q->next->back, mintrav, maxtrav);
      addTraverseBIG(tr, p, q->next->next->back, mintrav, maxtrav);    
    }
} 





int rearrangeBIG(tree *tr, nodeptr p, int mintrav, int maxtrav)   
{  
  double   p1z[NUM_BRANCHES], p2z[NUM_BRANCHES], q1z[NUM_BRANCHES], q2z[NUM_BRANCHES];
  nodeptr  p1, p2, q, q1, q2;
  int      mintrav2, i;  
  boolean doP = TRUE, doQ = TRUE;
  
  if (maxtrav < 1 || mintrav > maxtrav)  return 0;
  q = p->back;
  
  if(tr->constrained)
    {
      if(! tipHomogeneityChecker(tr, p->back, 0))
	doP = FALSE;
      
      if(! tipHomogeneityChecker(tr, q->back, 0))
	doQ = FALSE;
      
      if(doQ == FALSE && doP == FALSE)
	return 0;
    }
  
  if (/*!p->tip*/ !isTip(p->number, tr->rdta->numsp) && doP) 
    {     
      p1 = p->next->back;
      p2 = p->next->next->back;
      
      /*if (! p1->tip || ! p2->tip) */
      if(!isTip(p1->number, tr->rdta->numsp) || !isTip(p2->number, tr->rdta->numsp))
	{
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      p1z[i] = p1->z[i];
	      p2z[i] = p2->z[i];	   	   
	    }
	  
	  if (! removeNodeBIG(tr, p,  tr->numBranches)) return badRear;
	  
	  if (/*! p1->tip*/ !isTip(p1->number, tr->rdta->numsp)) 
	    {
	      addTraverseBIG(tr, p, p1->next->back,
			     mintrav, maxtrav);         
	      addTraverseBIG(tr, p, p1->next->next->back,
			     mintrav, maxtrav);          
	    }
	  
	  if (/*! p2->tip*/ !isTip(p2->number, tr->rdta->numsp)) 
	    {
	      addTraverseBIG(tr, p, p2->next->back,
			     mintrav, maxtrav);
	      addTraverseBIG(tr, p, p2->next->next->back,
			     mintrav, maxtrav);          
	    }
	  	  
	  hookup(p->next,       p1, p1z, tr->numBranches); 
	  hookup(p->next->next, p2, p2z, tr->numBranches);	   	    	    
	  initravDIST(tr, p, 1);	   	    
	}
    }  
  
  if (/*! q->tip*/ !isTip(q->number, tr->rdta->numsp) && maxtrav > 0 && doQ) 
    {
      q1 = q->next->back;
      q2 = q->next->next->back;
      
      /*if (((!q1->tip) && (!q1->next->back->tip || !q1->next->next->back->tip)) ||
	((!q2->tip) && (!q2->next->back->tip || !q2->next->next->back->tip))) */
      if (
	  (
	   ! isTip(q1->number, tr->rdta->numsp) && 
	   (! isTip(q1->next->back->number, tr->rdta->numsp) || ! isTip(q1->next->next->back->number, tr->rdta->numsp))
	   )
	  ||
	  (
	   ! isTip(q2->number, tr->rdta->numsp) && 
	   (! isTip(q2->next->back->number, tr->rdta->numsp) || ! isTip(q2->next->next->back->number, tr->rdta->numsp))
	   )
	  )
	{
	  
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      q1z[i] = q1->z[i];
	      q2z[i] = q2->z[i];
	    }
	  
	  if (! removeNodeBIG(tr, q, tr->numBranches)) return badRear;
	  
	  mintrav2 = mintrav > 2 ? mintrav : 2;
	  
	  if (/*! q1->tip*/ !isTip(q1->number, tr->rdta->numsp)) 
	    {
	      addTraverseBIG(tr, q, q1->next->back,
			     mintrav2 , maxtrav);
	      addTraverseBIG(tr, q, q1->next->next->back,
			     mintrav2 , maxtrav);         
	    }
	  
	  if (/*! q2->tip*/ ! isTip(q2->number, tr->rdta->numsp)) 
	    {
	      addTraverseBIG(tr, q, q2->next->back,
			     mintrav2 , maxtrav);
	      addTraverseBIG(tr, q, q2->next->next->back,
			     mintrav2 , maxtrav);          
	    }	   
	  
	  hookup(q->next,       q1, q1z, tr->numBranches); 
	  hookup(q->next->next, q2, q2z, tr->numBranches);
	  
	  initravDIST(tr, q, 1); 	   
	}
    } 
  
  return  1;
} 




double treeOptimizeRapid(tree *tr, int mintrav, int maxtrav, analdef *adef, bestlist *bt)
{
  int i, index,
    *perm = (int*)NULL;   

  nodeRectifier(tr);

  if (maxtrav > tr->ntips - 3)  
    maxtrav = tr->ntips - 3;  
    
  resetInfoList();
  
  resetBestTree(bt);
 
  tr->startLH = tr->endLH = tr->likelihood;
 
  if(tr->doCutoff)
    {
      if(tr->bigCutoff)
	{	  
	  if(tr->itCount == 0)    
	    tr->lhCutoff = 0.5 * (tr->likelihood / -1000.0);    
	  else    		 
	    tr->lhCutoff = 0.5 * ((tr->lhAVG) / ((double)(tr->lhDEC))); 	  
	}
      else
	{
	  if(tr->itCount == 0)    
	    tr->lhCutoff = tr->likelihood / -1000.0;    
	  else    		 
	    tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));   
	}    

      tr->itCount = tr->itCount + 1;
      tr->lhAVG = 0;
      tr->lhDEC = 0;
    }
  
  if(adef->permuteTreeoptimize)
    {
      int n = tr->mxtips + tr->mxtips - 2;   
      perm = (int *)malloc(sizeof(int) * (n + 1));
      makePermutation(perm, n, adef);
    }

  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
    {           
      tr->bestOfNode = unlikely;          

      if(adef->permuteTreeoptimize)
	index = perm[i];
      else
	index = i;
      
      if(rearrangeBIG(tr, tr->nodep[index], mintrav, maxtrav))
	{    
	  if(Thorough)
	    {
	      if(tr->endLH > tr->startLH)                 	
		{			   	     
		  restoreTreeFast(tr);	 	 
		  tr->startLH = tr->endLH = tr->likelihood;	 
		  saveBestTree(bt, tr);
		}
	      else
		{ 		  
		  if(tr->bestOfNode != unlikely)		    	     
		    restoreTopologyOnly(tr, bt);		    
		}	   
	    }
	  else
	    {
	      insertInfoList(tr->nodep[index], tr->bestOfNode);	    
	      if(tr->endLH > tr->startLH)                 	
		{		      
		  restoreTreeFast(tr);	  	      
		  tr->startLH = tr->endLH = tr->likelihood;	  	 	  	  	  	  	  	  
		}	    	  
	    }
	}     
    }
   


  if(!Thorough)
    {           
      Thorough = 1;  
      
      for(i = 0; i < iList.valid; i++)
	{      
	 
	  tr->bestOfNode = unlikely;
	  
	  if(rearrangeBIG(tr, iList.list[i].node, mintrav, maxtrav))
	    {	  
	      if(tr->endLH > tr->startLH)                 	
		{	 	     
		  restoreTreeFast(tr);	 	 
		  tr->startLH = tr->endLH = tr->likelihood;	 
		  saveBestTree(bt, tr);
		}
	      else
		{ 
	      
		  if(tr->bestOfNode != unlikely)
		    {	     
		      restoreTopologyOnly(tr, bt);
		    }	
		}      
	    }
	}       
          
      Thorough = 0;
    }

  if(adef->permuteTreeoptimize)
    free(perm);

  return tr->startLH;     
}


boolean testInsertRestoreBIG (tree *tr, nodeptr p, nodeptr q)
{    
  if(Thorough)
    {
      if (! insertBIG(tr, p, q, tr->numBranches))       return FALSE;    
      
      evaluateGeneric(tr, p->next->next);          


      /*	
	if (! insertRestoreBIG(tr, p, q))       return FALSE;
	
	{
	nodeptr x, y;
	
	x = p->next->next;
	y = p->back;
	while ((! x->x) || (! y->x)) 
	{
	if(! (x->x))
	newviewGeneric(tr, x);
	if (! (y->x)) 
	newviewGeneric(tr, y);
	}
	}
	
	tr->likelihood = tr->endLH;
      */

    }
  else
    {
      if (! insertRestoreBIG(tr, p, q))       return FALSE;
      
      {
	nodeptr x, y;
	x = p->next->next;
	y = p->back;
			
	if(! isTip(x->number, tr->rdta->numsp) && isTip(y->number, tr->rdta->numsp))
	  {
	    while ((! x->x)) 
	      {
		if (! (x->x))
		  newviewGeneric(tr, x);		     
	      }
	  }
	
	if(isTip(x->number, tr->rdta->numsp) && !isTip(y->number, tr->rdta->numsp))
	  {
	    while ((! y->x)) 
	      {		  
		if (! (y->x))
		  newviewGeneric(tr, y);
	      }
	  }
	
	if(!isTip(x->number, tr->rdta->numsp) && !isTip(y->number, tr->rdta->numsp))
	  {
	    while ((! x->x) || (! y->x)) 
	      {
		if (! (x->x))
		  newviewGeneric(tr, x);
		if (! (y->x))
		  newviewGeneric(tr, y);
	      }
	  }				      	
	
      }
	
      tr->likelihood = tr->endLH;
    }
     
  return TRUE;
} 

void restoreTreeFast(tree *tr)
{
  removeNodeRestoreBIG(tr, tr->removeNode);    
  testInsertRestoreBIG(tr, tr->removeNode, tr->insertNode);
}


int determineRearrangementSetting(tree *tr,  analdef *adef, bestlist *bestT, bestlist *bt)
{
  int i, mintrav, maxtrav, bestTrav, impr, index, MaxFast,
    *perm = (int*)NULL;
  double startLH; 
  boolean cutoff;  

  MaxFast = 26;

  startLH = tr->likelihood;

  cutoff = tr->doCutoff;
  tr->doCutoff = FALSE;
 
    
  mintrav = 1;
  maxtrav = 5;

  bestTrav = maxtrav = 5;

  impr = 1;

  resetBestTree(bt);

  if(adef->permuteTreeoptimize)
    {
      int n = tr->mxtips + tr->mxtips - 2;   
      perm = (int *)malloc(sizeof(int) * (n + 1));
      makePermutation(perm, n, adef);
    }
  

  while(impr && maxtrav < MaxFast)
    {	
      recallBestTree(bestT, 1, tr);     
      
      /* TODO, why are nodes not rectified here ? */
      
      if (maxtrav > tr->ntips - 3)  
	maxtrav = tr->ntips - 3;    
 
      tr->startLH = tr->endLH = tr->likelihood;
          
      for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
	{                

	  if(adef->permuteTreeoptimize)
	    index = perm[i];
	  else
	    index = i;	 	 

	  tr->bestOfNode = unlikely;
	  if(rearrangeBIG(tr, tr->nodep[index], mintrav, maxtrav))
	    {	     
	      if(tr->endLH > tr->startLH)                 	
		{		 	 	      
		  restoreTreeFast(tr);	        	  	 	  	      
		  tr->startLH = tr->endLH = tr->likelihood;	  	 	  	  	  	  	  	  	      
		}	         	       	
	    }
	}
      
      treeEvaluate(tr, 0.25);
      saveBestTree(bt, tr);                                    

      /*printf("DETERMINE_BEST: %d %f\n", maxtrav, tr->likelihood);*/

      if(tr->likelihood > startLH)
	{	 
	  startLH = tr->likelihood; 	  	  	  
	  printLog(tr, adef, FALSE);	  
	  bestTrav = maxtrav;	 
	  impr = 1;
	}
      else
	{
	  impr = 0;
	}
      maxtrav += 5;
      
      if(tr->doCutoff)
	{
	  tr->lhCutoff = (tr->lhAVG) / ((double)(tr->lhDEC));       
  
	  tr->itCount =  tr->itCount + 1;
	  tr->lhAVG = 0;
	  tr->lhDEC = 0;
	}
    }

  recallBestTree(bt, 1, tr);   
  tr->doCutoff = cutoff;

  if(adef->permuteTreeoptimize)
    free(perm);

  
  return bestTrav;     
}


#ifdef _MULTI_GENE
static void analyzeMultiGene(tree *tr)
{
  int model, i, j;
  boolean complete = FALSE;
  int modelCount[NUM_BRANCHES];
  int totalweight   = 0;
  int missingweight = 0;
  
  for(i = 0; i < tr->cdta->endsite; i++)    
    totalweight += tr->cdta->aliaswgt[i];

  totalweight *= tr->mxtips;
  

  for(i = 0; i < NUM_BRANCHES; i++)
    modelCount[i] = 0;

  assert(tr->NumberOfModels > 1 && tr->multiBranch);
  for(i = 1; i <= tr->mxtips; i++)
    {
      char *tip = tr->yVector[i];
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  int lower = tr->partitionData[model].lower;
	  int upper = tr->partitionData[model].upper;
	  char missing = 1;
	  int mw = 0;

	  for(j = lower; j < upper && missing; j++)
	    {
	      mw += tr->cdta->aliaswgt[j];
	      if(tip[j] != 15)		
		missing = 0;		
	    }

	  tr->tipMissing[i][model] = missing;

	  if(missing)
	    {
	      missingweight += mw;
	      /*printf("Seq %d part %d completely missing\n", i, model);	      */
	    }
	  else
	    modelCount[model] = modelCount[model] + 1;
	}
    }

  printf("\nSwitching to multi-gene model\n");
  printf("Sampling-induced gapyness: %f\n\n", (double)missingweight / (double)totalweight);

  /* now let's see if we can find at least one complete sequence */ 

  /*for(i = 0; i < tr->NumberOfModels; i++)
    printf("Partition %d: %d sequences\n", i, modelCount[i]);*/

  for(i = 1; i <= tr->mxtips; i++)
    {
      complete = TRUE;

      for(model = 0; model < tr->NumberOfModels && complete; model++)	
	{
	  if(tr->tipMissing[i][model])
	    complete = FALSE;	
	}

      if(complete)
	{
	  int j;
	  tr->start = tr->nodep[i];
	  /*printf("%d: ", i);*/
	  for(j = 0; j < tr->NumberOfModels; j++)
	    {
	      /* printf("%d ", tr->tipMissing[i][j]); */
	      tr->startVector[j] = tr->nodep[i];
	    }
	  /* printf("\n"); */
	  break;
	}
      else
	{
	  /*printf("Tip %d missing part %d\n", i, model);*/
	}
    }

  /*if(complete)
    return;
    else*/
  tr->start = (nodeptr)NULL;
    {
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  for(i = 1; i <= tr->mxtips; i++)
	    {
	      if(!tr->tipMissing[i][model])
		{
		  /*printf("START %d\n", i);*/
		  tr->startVector[model] = tr->nodep[i];
		  break;
		}
	    }
	}
      /*printf("TODO\n");
	assert(0);*/
    }
}



static int containsModel(nodeptr p, int model, tree *tr)
{
  if(isTip(p->number, tr->mxtips))
    {
      if(tr->tipMissing[p->number][model])
	return 0;
      else
	return 1;
    }
  else
    {
      assert(p = p->next->next->next);
      return (containsModel(p->next->back, model, tr) || containsModel(p->next->next->back, model, tr));
    }
}


static void reduceTreeModelREC(nodeptr p, nodeptr reference, int model, tree *tr)
{

  if(isTip(p->number, tr->mxtips))
    {
      assert(!tr->tipMissing[p->number][model]);

      p->backs[model]         = reference;
      reference->backs[model] = p;
      p->z[model] = reference->z[model] = defaultz;
    }
  else
    {
      nodeptr q = p->next;
      nodeptr r = p->next->next;
      
      int left  = containsModel(q->back, model, tr);      
      int right = containsModel(r->back, model, tr);

      assert(p = p->next->next->next);

      if(left && right)
	{
	  p->backs[model]         = reference;
	  reference->backs[model] = p;
	  p->z[model] = reference->z[model] = defaultz;

	  reduceTreeModelREC(q->back, q,  model, tr);   	 	
	  reduceTreeModelREC(r->back, r,  model, tr);
	} 
      else
	{
	  if(left || right)
	    {
	      if(left)
		{
		  reduceTreeModelREC(q->back, reference,  model, tr);
		  /* contained in q schtrawutsni */
		}
	      else
		{
		  reduceTreeModelREC(r->back, reference,  model, tr);
		  /* contained in r schtrawutsni */
		}
	    }
	  else
	    {
	      assert(0);
	      /* should not get here */
	    }
	}     
    }
}


static void reduceTreeModel(tree *tr, int model)
{
  assert(isTip(tr->startVector[model]->number, tr->mxtips));

  reduceTreeModelREC(tr->startVector[model]->back, tr->startVector[model], model, tr);
}


static void printTreeModelRec(nodeptr p, int model, tree *tr, int *tips)
{
  if(isTip(p->number, tr->mxtips))
    {
      assert(!tr->tipMissing[p->number][model]);
      printf("%s,", tr->nameList[p->number]);
      *tips = *tips + 1;
    }
  else
    {
      if(p->next->backs[model])
	printTreeModelRec(p->next->backs[model], model, tr, tips);
      if(p->next->next->backs[model])
	printTreeModelRec(p->next->next->backs[model], model, tr, tips);
    }
}

static void printTreeModel(tree *tr, int model)
{
  int tips = 0;

  printf("Tree %d:", model);
  
  printTreeModelRec(tr->startVector[model], model, tr, &tips);
  printTreeModelRec(tr->startVector[model]->backs[model], model, tr, &tips);
  printf(" %d tips \n", tips);
}

static void treeReduction(tree *tr)
{
  int model;

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      reduceTreeModel(tr, model);
      /*printTreeModel(tr, model);*/
    }

}

#endif


void computeBIGRAPID (tree *tr, analdef *adef) 
{ 
  int i,  impr, bestTrav,
    rearrangementsMax = 0, 
    rearrangementsMin = 0;
   
  double lh, previousLh, difference, epsilon;              
  bestlist *bestT, *bt;  
  
  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);
      
  bt = (bestlist *) malloc(sizeof(bestlist));      
  bt->ninit = 0;
  initBestTree(bt, 20, tr->mxtips); 

  initInfoList(50);
 
  difference = 10.0;
  epsilon = 0.01;    
    
  Thorough = 0;  
 
#ifndef _MULTI_GENE
  optimizeModel(tr, adef);       
  treeEvaluate(tr, 2);   
 
  printLog(tr, adef, FALSE);  
#endif

  saveBestTree(bestT, tr);
  
#ifdef _MULTI_GENE
  {
    double t;
    int k;

    printf("Standard Method:\n");

    resetBranches(tr);
    t = gettime();
    for(k = 0; k < 50; k++)
      evaluateGenericInitrav(tr, tr->start); 
    t = gettime() - t;
    printf("50 tree traversals: \t %f %f secs\n", tr->likelihood, t);
     
    t = gettime();
    treeEvaluate(tr, 2.0);
    t = gettime() - t;
    printf("Br-Len Opt: \t\t %f %f secs\n", tr->likelihood, t);

    resetBranches(tr);
    analyzeMultiGene(tr);
    treeReduction(tr);

    t = gettime();

    tr->doMulti = 1;
    for(k = 0; k < 50; k++)
      {
	/*
	  result = 0.0;
	  determineFullMultiTraversal(tr);   
	  for(model = 0; model < tr->NumberOfModels; model++)          
	  result += evaluateIterativePartition(tr, tr->partitionData[model].lower, tr->partitionData[model].upper, model);       
	*/
	evaluateGenericInitrav(tr, tr->start); 

      }
    t = gettime() - t;
    printf("50 tree traversals: \t %f %f secs\n", tr->likelihood, t);
    /*for(k = 0; k < tr->NumberOfModels; k++)
      printf("Part[%d]: %f\n", k,  tr->perPartitionLH[k]);*/
    
    t = gettime();
    /*printf("D1\n");*/
    treeEvaluateMulti(tr, 2.0);
    /*printf("D2\n");*/
    t = gettime() - t;
    printf("Br-Len Opt: \t\t %f %f secs\n", tr->likelihood, t);
    


    /*

    for(k = 1; k < 100; k++)
      {
	double erg = evaluateGeneric(tr, tr->nodep[k]);
	printf("%d: %f\n", k, erg);
	}*/

    exit(1);
  }
#endif

  if(!adef->initialSet)   
    bestTrav = adef->bestTrav = determineRearrangementSetting(tr, adef, bestT, bt);                   
  else
    bestTrav = adef->bestTrav = adef->initial;



  saveBestTree(bestT, tr); 
  impr = 1;
  if(tr->doCutoff)
    tr->itCount = 0;
 
  while(impr)
    {              
      recallBestTree(bestT, 1, tr);     
      optimizeModel(tr, adef);            
      treeEvaluate(tr, 2);	 	      
      saveBestTree(bestT, tr);           
      printLog(tr, adef, FALSE);            
      printResult(tr, adef, FALSE);    
      lh = previousLh = tr->likelihood;
   
     
      treeOptimizeRapid(tr, 1, bestTrav, adef, bt);   
      
      impr = 0;
	  
      for(i = 1; i <= bt->nvalid; i++)
	{	    		  	   
	  recallBestTree(bt, i, tr);	 
	  treeEvaluate(tr, 0.25);	    	 		      	 

	  difference = ((tr->likelihood > previousLh)? 
			tr->likelihood - previousLh: 
			previousLh - tr->likelihood); 	    
	  if(tr->likelihood > lh && difference > epsilon)
	    {
	      impr = 1;	       
	      lh = tr->likelihood;	       	     
	      saveBestTree(bestT, tr);
	    }	   	   
	}	
    }

  Thorough = 1;
  impr = 1;
  
  while(1)
    {		
      recallBestTree(bestT, 1, tr);    
      if(impr)
	{	    
	  printResult(tr, adef, FALSE);
	  rearrangementsMin = 1;
	  rearrangementsMax = adef->stepwidth;	    
	}			  			
      else
	{		       	   
	  rearrangementsMax += adef->stepwidth;
	  rearrangementsMin += adef->stepwidth; 	        	      
	  if(rearrangementsMax > adef->max_rearrange)	     	     	 
	    goto cleanup; 	   
	}
      
      optimizeModel(tr, adef);       
      treeEvaluate(tr, 2.0);	      
      previousLh = lh = tr->likelihood;	      
      saveBestTree(bestT, tr);     
      printLog(tr, adef, FALSE);

      treeOptimizeRapid(tr, rearrangementsMin, rearrangementsMax, adef, bt);
	
      impr = 0;			      	
		
      for(i = 1; i <= bt->nvalid; i++)
	{	

	  recallBestTree(bt, i, tr);	    
	  
	  treeEvaluate(tr, 0.25);	    	 
	  
	  difference = ((tr->likelihood > previousLh)? 
			tr->likelihood - previousLh: 
			previousLh - tr->likelihood); 	    
	  if(tr->likelihood > lh && difference > epsilon)
	    {
	      impr = 1;	       
	      lh = tr->likelihood;	  	     
	      saveBestTree(bestT, tr);
	    }	   	   
	}	
    }

 cleanup:   
  freeBestTree(bestT);
  free(bestT);
  freeBestTree(bt);
  free(bt);
  freeInfoList();
  printLog(tr, adef, FALSE);
  printResult(tr, adef, FALSE);
}


void computeBIGRAPIDMULTIBOOT (tree *tr, analdef *adef) 
{ 
  int i,  impr, bestTrav,
    rearrangementsMax = 0, 
    rearrangementsMin = 0;
   
  double lh, previousLh, difference, epsilon;              
  bestlist *bestT, *bt;    
  
  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);
      
  bt = (bestlist *) malloc(sizeof(bestlist));      
  bt->ninit = 0;
  initBestTree(bt, 20, tr->mxtips); 

  initInfoList(50);
 
  difference = 10.0;
  epsilon = 0.01;    
    
  Thorough = 0; 
 
  treeEvaluate(tr, 2);   
   
  printLog(tr, adef, FALSE);  

  saveBestTree(bestT, tr);

  if(!adef->initialSet)   
    bestTrav = adef->bestTrav = determineRearrangementSetting(tr, adef, bestT, bt);                   
  else
    bestTrav = adef->bestTrav = adef->initial;

  saveBestTree(bestT, tr); 
  impr = 1;
  if(tr->doCutoff)
    tr->itCount = 0;

  while(impr)
    {              
      recallBestTree(bestT, 1, tr);      
      treeEvaluate(tr, 2);	
      /*printf("%f\n", tr->likelihood);*/
      saveBestTree(bestT, tr);     
      printLog(tr, adef, FALSE);     
      printResult(tr, adef, FALSE);    
      lh = previousLh = tr->likelihood;
         
      treeOptimizeRapid(tr, 1, bestTrav, adef, bt);   
      
      impr = 0;
	  
      for(i = 1; i <= bt->nvalid; i++)
	{	    		  	   
	  recallBestTree(bt, i, tr);	    
	  treeEvaluate(tr, 0.25);	    	 	
	      
	  difference = ((tr->likelihood > previousLh)? 
			tr->likelihood - previousLh: 
			previousLh - tr->likelihood); 	    
	  if(tr->likelihood > lh && difference > epsilon)
	    {
	      impr = 1;	       
	      lh = tr->likelihood;	       	     
	      saveBestTree(bestT, tr);
	    }	   	   
	}	
    }

  Thorough = 1;
  impr = 1;
  
  while(1)
    {		
      recallBestTree(bestT, 1, tr);    
      if(impr)
	{	    
	  printResult(tr, adef, FALSE);
	  rearrangementsMin = 1;
	  rearrangementsMax = adef->stepwidth;	    
	}			  			
      else
	{		       	   
	  rearrangementsMax += adef->stepwidth;
	  rearrangementsMin += adef->stepwidth; 	        	      
	  if(rearrangementsMax > adef->max_rearrange)	     	     	 
	    goto cleanup; 	   
	}
               
      treeEvaluate(tr, 2.0);
      /*printf("%f\n", tr->likelihood);*/
      previousLh = lh = tr->likelihood;	      
      saveBestTree(bestT, tr);     
      printLog(tr, adef, FALSE);

      treeOptimizeRapid(tr, rearrangementsMin, rearrangementsMax, adef, bt);
	
      impr = 0;			      	
		
      for(i = 1; i <= bt->nvalid; i++)
	{	

	  recallBestTree(bt, i, tr);	    
	  
	  treeEvaluate(tr, 0.25);	    	 
	  
	  difference = ((tr->likelihood > previousLh)? 
			tr->likelihood - previousLh: 
			previousLh - tr->likelihood); 	    
	  if(tr->likelihood > lh && difference > epsilon)
	    {
	      impr = 1;	       
	      lh = tr->likelihood;	  	     
	      saveBestTree(bestT, tr);
	    }	   	   
	}	
    }

 cleanup:   
  freeBestTree(bestT);
  free(bestT);
  freeBestTree(bt);
  free(bt);
  freeInfoList();
  printLog(tr, adef, FALSE);
  printResult(tr, adef, FALSE);
}


boolean treeEvaluate (tree *tr, double smoothFactor)       /* Evaluate a user tree */
  { /* treeEvaluate */
    
    /*double inLH = tr->likelihood;*/

    if (! smoothTree(tr, (int)((double)smoothings * smoothFactor))) 
      {
	return FALSE;      
      }
      
    evaluateGeneric(tr, tr->start);   
   

    /*    if(inLH > tr->likelihood)
      {
	printf("FATAL error in treeEvaluate %.20f <-> %.20f factor %d\n", inLH, tr->likelihood, (int)((double)smoothings * smoothFactor));
	}*/

    return TRUE;
  } /* treeEvaluate */


/************* per partition branch length optimization ****************************/

static boolean updatePartition(tree *tr, nodeptr p, int model)
{
  nodeptr  q;
  double   z0, z;
	
  q = p->back;
  z0 = q->z[0];
        
  z = makenewzPartitionGeneric(tr, p, q, z0, newzpercycle, model);
    
  p->z[0] = q->z[0] = z;
  if (ABS(z - z0) > deltaz)  tr->smoothed = FALSE;
        
  return TRUE; 
}

static boolean smoothPartition(tree *tr, nodeptr p, int model)
{
  nodeptr  q;
  
  if(! updatePartition(tr, p, model))               
    return FALSE; 

  if(/*! p->tip*/ !isTip(p->number, tr->rdta->numsp)) 
    {                    
      q = p->next;
      while (q != p) 
	{
	  if (! smoothPartition(tr, q->back, model))   
	    return FALSE;
	  q = q->next;
	}	

      newviewPartitionGeneric(tr, p, model);
    }

  return TRUE;
} 


static boolean smoothTreePartition(tree *tr, int maxtimes, int model)
{
  nodeptr  p, q;    

  p = tr->start;

  while(--maxtimes >= 0) 
    {
      tr->smoothed = TRUE;

      if(! smoothPartition(tr, p->back, model))      
	return FALSE;

      if(/*! p->tip*/ !isTip(p->number, tr->rdta->numsp)) 
	  {
	    q = p->next;
	    while (q != p) 
	      {
		if(! smoothPartition(tr, q->back, model))   
		  return FALSE;
		q = q->next;
	      }
	  }
      if (tr->smoothed)  break;
    }

  return TRUE;
}


boolean treeEvaluatePartition(tree *tr, double smoothFactor, int model)
{    
  if(! smoothTreePartition(tr, (int)((double)smoothings * smoothFactor), model))     
    return FALSE;          
      
  evaluatePartitionGeneric(tr, tr->start, model);    
   
  return TRUE;
}

/******** Mehring Algo DEVEL ********************************************/


static boolean insertLight (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r, s;
  int i;
  
  r = q->back;
  s = p->back;
      
  for(i = 0; i < tr->numBranches; i++)
    tr->lzi[i] = q->z[i];
  
  if(Thorough)
    { 
      double  zqr[NUM_BRANCHES], zqs[NUM_BRANCHES], zrs[NUM_BRANCHES], lzqr, lzqs, lzrs, lzsum, lzq, lzr, lzs, lzmax;      
      double defaultArray[NUM_BRANCHES];	
      double e1[NUM_BRANCHES], e2[NUM_BRANCHES], e3[NUM_BRANCHES];
      double *qz;
      
      qz = q->z;
      
      for(i = 0; i < tr->numBranches; i++)
	defaultArray[i] = defaultz;
      
      makenewzGeneric(tr, q, r, qz, iterations, zqr);           
      makenewzGeneric(tr, q, s, defaultArray, iterations, zqs);                  
      makenewzGeneric(tr, r, s, defaultArray, iterations, zrs);
      
      
      for(i = 0; i < tr->numBranches; i++)
	{
	  lzqr = (zqr[i] > zmin) ? log(zqr[i]) : log(zmin); 
	  lzqs = (zqs[i] > zmin) ? log(zqs[i]) : log(zmin);
	  lzrs = (zrs[i] > zmin) ? log(zrs[i]) : log(zmin);
	  lzsum = 0.5 * (lzqr + lzqs + lzrs);
	  
	  lzq = lzsum - lzrs;
	  lzr = lzsum - lzqs;
	  lzs = lzsum - lzqr;
	  lzmax = log(zmax);
	  
	  if      (lzq > lzmax) {lzq = lzmax; lzr = lzqr; lzs = lzqs;} 
	  else if (lzr > lzmax) {lzr = lzmax; lzq = lzqr; lzs = lzrs;}
	  else if (lzs > lzmax) {lzs = lzmax; lzq = lzqs; lzr = lzrs;}          
	  
	  e1[i] = exp(lzq);
	  e2[i] = exp(lzr);
	  e3[i] = exp(lzs);
	}
      hookup(p->next,       q, e1, tr->numBranches);
      hookup(p->next->next, r, e2, tr->numBranches);
      hookup(p,             s, e3, tr->numBranches);      		  
    }
  else
    {       
      double  z[NUM_BRANCHES]; 
      
      for(i = 0; i < tr->numBranches; i++)
	{
	  z[i] = sqrt(q->z[i]);      
	  
	  if(z[i] < zmin) 
	    z[i] = zmin;
	  if(z[i] > zmax)
	    z[i] = zmax;
	}
      
      hookup(p->next,       q, z, tr->numBranches);
      hookup(p->next->next, r, z, tr->numBranches);	                         
    }
  
  newviewGeneric(tr, p);
  
  if(Thorough)
    {     
      localSmooth(tr, p, smoothings);
      /*regionalSmooth (tr, p, int maxtimes, 10);*/
      

      for(i = 0; i < tr->numBranches; i++)
	{
	  tr->lzq[i] = p->next->z[i];
	  tr->lzr[i] = p->next->next->z[i];
	  tr->lzs[i] = p->z[i];            
	}
    }           
  
  return  TRUE;
}


static boolean testInsertLight(tree *tr, nodeptr p, nodeptr q, int *count)
{
  double  qz[NUM_BRANCHES], pz[NUM_BRANCHES];
  nodeptr  r;
  boolean doIt = TRUE; 
  int i;
  
  r = q->back; 
  for(i = 0; i < tr->numBranches; i++)
    {
      qz[i] = q->z[i];
      pz[i] = p->z[i];
    }
  
  if(tr->grouped)
    {
      int rNumber, qNumber, pNumber;
      
      doIt = FALSE;
      
      rNumber = tr->constraintVector[r->number];
      qNumber = tr->constraintVector[q->number];
      pNumber = tr->constraintVector[p->number];
      
      if(pNumber == -9)
	pNumber = checker(tr, p->back);
      if(pNumber == -9)
	doIt = TRUE;
      else
	{
	  if(qNumber == -9)
	    qNumber = checker(tr, q);
	  
	  if(rNumber == -9)
	    rNumber = checker(tr, r);
	  
	  if(pNumber == rNumber || pNumber == qNumber)
	    doIt = TRUE;       
	}
    }
  
  if(doIt)
    {     
      if (! insertLight(tr, p, q))       return FALSE;         
      
      evaluateGeneric(tr, p->next->next);       
      
      *count = *count + 1;
      		   
      if(tr->likelihood > tr->bestOfNode)
	{
	  tr->bestOfNode = tr->likelihood;
	  tr->insertNode = q;
	  tr->removeNode = p;   
	  for(i = 0; i < tr->numBranches; i++)
	    {
	      tr->currentZQR[i] = tr->zqr[i];           
	      tr->currentLZR[i] = tr->lzr[i];
	      tr->currentLZQ[i] = tr->lzq[i];
	      tr->currentLZS[i] = tr->lzs[i];      
	    }
	}
      
         
      
      hookup(q, r, qz, tr->numBranches);
      
      p->next->next->back = p->next->back = (nodeptr) NULL;
      
      if(Thorough)
	{
	  nodeptr s = p->back;
	  hookup(p, s, pz, tr->numBranches);      
	} 
           
      return TRUE;
    }
  else
    return TRUE;  
}





static void addTraverseLight(tree *tr, nodeptr p, nodeptr q, int *count)
{                
  if (! testInsertLight(tr, p, q, count))         
    assert(0);

  if(!isTip(q->number, tr->rdta->numsp))
    {   
      addTraverseLight(tr, p, q->next->back, count);
      addTraverseLight(tr, p, q->next->next->back, count);    
    }
} 




static void findBestPosition(tree *tr, insertionBranch *ib, 
			     int numberOfInsertionBranches, nodeptr testSequence, nodeptr q, nodeptr r,
			     nodeptr s)
{  
  int count = 0;
  int i;
  boolean found = FALSE;     

  /* modOpt(tr, adef); */
 
  assert(testSequence != tr->start);
 
  /*printf("ENTER: %f\n", tr->likelihood);*/

  tr->bestOfNode = unlikely;

  assert(r->back == s && s->back == r);

  newviewGeneric(tr, r);
  newviewGeneric(tr, s);
  
  hookupDefault(q, testSequence, tr->numBranches);

  testInsertLight(tr, q, r, &count);   

  if(!isTip(r->number, tr->mxtips))
    {
      addTraverseLight(tr, q, r->next->back, &count);
      addTraverseLight(tr, q, r->next->next->back, &count);   
    }
  
  if(!isTip(s->number, tr->mxtips))
    {
      addTraverseLight(tr, q, s->next->back, &count);
      addTraverseLight(tr, q, s->next->next->back, &count);   
    }

  
  /*printf("Trees analyzed: %d Best %f insertion at node %d (back %d)\n", count, tr->bestOfNode, 
    tr->insertNode->number, tr->insertNode->back->number);   */
  
         

  for(i = 0; i < numberOfInsertionBranches && !found; i++)
    {
      if(ib[i].p == tr->insertNode)
	{
	  ib[i].freq = ib[i].freq + 1;
	  found = TRUE;
	}
      else
	{
	  if(ib[i].p == (nodeptr)NULL)
	    break;
	}
    }

  if(!found)
    {
      assert(i < numberOfInsertionBranches);
      ib[i].freq = 1;
      ib[i].p = tr->insertNode;
    } 
}




void determineSequencePosition(tree *tr, analdef *adef)
{
  int i;
  int numberOfInsertionBranches;  
  FILE *f;
  nodeptr testSequence, q, r, s;
  int *originalRateCategories = (int*)malloc(tr->cdta->endsite * sizeof(int));      
  int *originalInvariant      = (int*)malloc(tr->cdta->endsite * sizeof(int));

  assert(adef->restart && adef->outgroup);

  adef->outgroup = FALSE;
  tr->doCutoff   = FALSE;
  Thorough = 1;

  numberOfInsertionBranches = 2 * tr->mxtips - 5; 
  tr->ib = (insertionBranch *)malloc(sizeof(insertionBranch) * numberOfInsertionBranches);
  
  for(i = 0; i < numberOfInsertionBranches; i++)
    {
      tr->ib[i].p    = (nodeptr)NULL;
      tr->ib[i].freq = 0;
    }


  if(tr->outgroupNums[0] == 1)
    tr->start = tr->nodep[2];

  testSequence = tr->nodep[tr->outgroupNums[0]];
  
  assert(testSequence != tr->start);

  q = testSequence->back;
 
  r = q->next->back;
  s = q->next->next->back;
  assert(q->next->next->next->back == testSequence);

  q->next->back       = (nodeptr)NULL;
  q->next->next->back = (nodeptr)NULL;

  hookupDefault(r, s, tr->numBranches);
  
  modOpt(tr, adef);

  memcpy(originalRateCategories, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);
  memcpy(originalInvariant,      tr->invariant,          sizeof(int) * tr->cdta->endsite);

  if(adef->rapidBoot)
    {            
      for(i = 0; i < adef->multipleRuns; i++)
	{
	  double t = gettime();	  	
	  	  
	  computeNextReplicate(tr, adef, originalRateCategories, originalInvariant);
	  resetBranches(tr);
	  evaluateGenericInitrav(tr, tr->start);    
	  treeEvaluate(tr, 1);

	  findBestPosition(tr, tr->ib, numberOfInsertionBranches, testSequence, q, r, s);	 
	  t = gettime() - t;
	  printf("Replicate[%d]: %f seconds\n", i, t);	 
	}     
    }
  else
    findBestPosition(tr, tr->ib, numberOfInsertionBranches, testSequence, q, r, s);

  i = 0;
  while(i < numberOfInsertionBranches && tr->ib[i].p)
    {
      double support = ((double)(tr->ib[i].freq)) / ((double) (adef->multipleRuns));
#ifdef WIN32
      tr->ib[i].freq = (int)floor(0.5 + support * 100.0);
#else
      tr->ib[i].freq = (int)(0.5 + support * 100.0);
#endif           
      i++;
    }

  f = fopen(resultFileName, "w");

  Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, 
	      FALSE, FALSE, FALSE, adef, SUMMARIZE_LH);  

  fprintf(f, "%s", tr->tree_string);

  fclose(f);
  printf("\nResult written to file: %s\n", resultFileName);

  exit(0);
}


static void markTips(nodeptr p, int *perm, int maxTips)
{
  if(isTip(p->number, maxTips))
    {
      perm[p->number] = 1;
      return;
    }
  else
    {
      nodeptr q = p->next;
      while(q != p)
	{
	  markTips(q->back, perm, maxTips);
	  q = q->next;
	}
      
    }
}

static double testInsertRob (tree *tr, nodeptr p, nodeptr q)
{
  double  qz[NUM_BRANCHES], pz[NUM_BRANCHES];
  nodeptr  r; 
  double result;
  int i;
    
  r = q->back; 

  for(i = 0; i < tr->numBranches; i++)
    {
      qz[i] = q->z[i];
      pz[i] = p->z[i];
    }
  
  assert(!tr->grouped);

  insertBIG(tr, p, q, tr->numBranches);       
  
  result = evaluateGeneric(tr, p->next->next);             		              
  
  hookup(q, r, qz, tr->numBranches);
  
  p->next->next->back = p->next->back = (nodeptr) NULL;
  
  if(Thorough)
    {
      nodeptr s = p->back;
      hookup(p, s, pz, tr->numBranches);      
    } 
  
  return result;
}


 


static void addTraverseRob(tree *tr, nodeptr r, nodeptr q, int *inserts,  int numberOfTipsForInsertion, int *count,
			   double *insertLikelihoods, nodeptr *insertNodes)
{                
  int i;

  for(i = 0; i < numberOfTipsForInsertion; i++)
    {     
      double result;
      nodeptr p = tr->nodep[inserts[i]];     

      hookupDefault(p, r, tr->numBranches);      
      result = testInsertRob(tr, r, q);
      
      if(result > insertLikelihoods[i])
	{
	  insertLikelihoods[i] = result;
	  insertNodes[i]       = q;
	}
    }

  /*printf("Sierra %d\n", *count);*/
  
  *count = *count + 1;

  if(!isTip(q->number, tr->rdta->numsp))
    {   
      addTraverseRob(tr, r, q->next->back,       inserts,  numberOfTipsForInsertion, count, insertLikelihoods, insertNodes);
      addTraverseRob(tr, r, q->next->next->back, inserts,  numberOfTipsForInsertion, count, insertLikelihoods, insertNodes);    
    }
} 

static void insertInsertionPoint(insertionPoints *ip, nodeptr p)
{
  int pos = 0;

  for(pos = 0; pos <  ip->max; pos++)
    {
      if(ip->ib[pos].p == (nodeptr)NULL)
	break;     	
      if(ip->ib[pos].p == p)
	break;
    }  

  if(pos < ip->max)
    {
      if(ip->ib[pos].p == (nodeptr)NULL)
	{	
	  ip->ib[pos].p = p;
	  ip->ib[pos].freq = 1;	 
	}
      else
	{
	  /*printf("adding at %d\n", pos, ip->max);*/
	  assert( ip->ib[pos].p == p);
	  ip->ib[pos].freq = ip->ib[pos].freq + 1;
	}
    }
  else
    {
      int i;
      insertionBranch *new;

      printf("extending\n");

      new = (insertionBranch *)malloc(sizeof(insertionBranch) * 2 * ip->max);
      memcpy(new, ip->ib, sizeof(insertionBranch) * ip->max);
      for(i = ip->max; i < 2 * ip->max; i++)
	{
	  new[i].p    = (nodeptr)NULL;
	  new[i].freq = 0;
	}

      new[pos].p    = p;
      new[pos].freq = 1;
      
      ip->max   = 2 * ip->max;

      free(ip->ib);
      ip->ib = new;
    }

}



static void getIP(nodeptr p, tree *tr, int *count, int *inserts, int *frequencies)
{
  int i, j;

  *count = 0;

  for(j = 0; j < tr->numberOfTipsForInsertion; j++)
    {
      i = 0;
      while(i < tr->ip[j].max && tr->ip[j].ib[i].p)
	{
	  assert(tr->ip[j].ib[i].freq > 0 && tr->ip[j].ib[i].p);
	  if(p == tr->ip[j].ib[i].p || p->back == tr->ip[j].ib[i].p)
	    {
	      inserts[*count]     = tr->ip[j].reference;
	      frequencies[*count] = tr->ip[j].ib[i].freq;
	      *count = *count + 1;
	      break;
	    }	    
	  i++;
	}
    }

  /*printf("Count %d %d\n", *count, p->number);*/



  return;  
}

static char *Tree2StringML(char *treestr, tree *tr, nodeptr p, int *inserts, int *frequencies)
{  
  char  *nameptr;  
  int   count, i;
  
  assert(tr->ip);

  getIP(p, tr, &count, inserts, frequencies);
  
      
  if(isTip(p->number, tr->rdta->numsp)) 
    {	 
      if(count > 0)
	{
	  if(count > 1)
	    printf("TIP > 1\n");

	  nameptr = tr->nameList[p->number];   
	  getIP(p, tr, &count, inserts, frequencies);	 

	  if(count == 1)	    
	    sprintf(treestr, "(%s:%d,%s)", tr->nameList[inserts[0]], frequencies[0], nameptr);
	  else
	    {	     
	      sprintf(treestr, "(%s:%d", tr->nameList[inserts[0]], frequencies[0]);
	      while (*treestr) treestr++;

	      for(i = 1; i < count; i++)
		{
		  sprintf(treestr, ",%s:%d", tr->nameList[inserts[i]], frequencies[i]);
		  while (*treestr) treestr++;
		}			    
  
	      sprintf(treestr, "),%s)", nameptr);
	      while (*treestr) treestr++;
	    }	
	}
      else
	{
	  nameptr = tr->nameList[p->number];     
	  sprintf(treestr, "%s", nameptr);
	}

      while (*treestr) treestr++;
    }
  else 
    {           
      if(count > 0)
	{	
	  if(count > 1)
	    printf("INNER > 1\n");

	  *treestr++ = '(';
	  
	  if(count == 1)
	    {
	      sprintf(treestr, "%s:%d", tr->nameList[inserts[0]], frequencies[0]);
	      while (*treestr) treestr++;
	    }
	  else
	    {
	      *treestr++ = '(';
	      sprintf(treestr, "%s:%d", tr->nameList[inserts[0]], frequencies[0]);
	      while (*treestr) treestr++;
	      for(i = 1;  i < count; i++)
		{
		  sprintf(treestr, ",%s:%d", tr->nameList[inserts[i]], frequencies[i]);	 
		  while (*treestr) treestr++;
		}
	      *treestr++ = ')';
	    }
	  

	  *treestr++ = ',';
	  *treestr++ = '(';
	  treestr = Tree2StringML(treestr, tr, p->next->back, inserts, frequencies);
	  *treestr++ = ',';
	  treestr = Tree2StringML(treestr, tr, p->next->next->back, inserts, frequencies);
	  if(p == tr->start->back) 
	    {
	      *treestr++ = ',';
	      treestr = Tree2StringML(treestr, tr, p->back, inserts, frequencies);
	    }
	  *treestr++ = ')'; 
	  *treestr++ = ')'; 
	}
      else
	{	 
	  *treestr++ = '(';
	  treestr = Tree2StringML(treestr, tr, p->next->back, inserts, frequencies);
	  *treestr++ = ',';
	  treestr = Tree2StringML(treestr, tr, p->next->next->back, inserts, frequencies);
	  if(p == tr->start->back) 
	    {
	      *treestr++ = ',';
	      treestr = Tree2StringML(treestr, tr, p->back, inserts, frequencies);
	    }
	  *treestr++ = ')';                
	}    
    }

  if(p == tr->start->back) 
    {	
      if(count > 0)
	{
	  /*sprintf(treestr, "%d;\n", insertionFrequency);*/
	  sprintf(treestr, ";\n");
	}
      else
	{	 	 
	  sprintf(treestr, ";\n");	 	  
	}
    }
  else 
    {      
      if(count > 0)
	{
	  /*sprintf(treestr, "%d", insertionFrequency);*/
	  sprintf(treestr, "%s", "\0");
	}
      else
	{	 	  
	  sprintf(treestr, "%s", "\0");	    
	}    
    }

  while (*treestr) treestr++;
  return  treestr;
}


void rapidML_Addition(tree *tr, analdef *adef)
{
  int i, j, *perm, *inserts;
  int numberOfTipsForInsertion = 0;
  int count;  
  double  *insertLikelihoods;
  insertionPoints *ip;
  nodeptr *insertNodes;  
  nodeptr r, q;

  printf("Current Tips: %d\n", tr->ntips);
  
  assert(adef->restart);

  adef->outgroup = FALSE;
  tr->doCutoff   = FALSE;

  evaluateGenericInitrav(tr, tr->start);
  printf("Init: %f\n", tr->likelihood);
  
  
  if(adef->model == M_PROTGAMMA || adef->model == M_GTRGAMMA)    
    modOpt(tr, adef);    
  else    
    quickOpt(tr, adef);    
 
  printf("Opt: %f\n", tr->likelihood);
 
  Thorough = 0;


  perm    = (int *)calloc(tr->mxtips + 1, sizeof(int));
  inserts = (int *)calloc(tr->mxtips, sizeof(int));

  markTips(tr->start,       perm, tr->mxtips);
  markTips(tr->start->back, perm ,tr->mxtips);

  numberOfTipsForInsertion = 0;

  for(i = 1; i <= tr->mxtips; i++)
    {
      if(perm[i] == 0)
	{
	  inserts[numberOfTipsForInsertion] = i;
	  numberOfTipsForInsertion++;
	}
    }
  
  /*numberOfTipsForInsertion = 10;*/

  printf("RAxML will insert %d Sequences\n",  numberOfTipsForInsertion);

  insertLikelihoods = (double*)malloc(sizeof(double) *  numberOfTipsForInsertion);
  insertNodes       = (nodeptr*)malloc(sizeof(nodeptr) *  numberOfTipsForInsertion);
  ip                = (insertionPoints *)malloc(sizeof(insertionPoints) *  numberOfTipsForInsertion);
    
  for(i = 0; i < numberOfTipsForInsertion; i++) 
    {            
      ip[i].max       = 8;
      ip[i].reference = inserts[i];
      ip[i].ib        = (insertionBranch *)malloc(sizeof(insertionBranch) * 8);
      
      for(j = 0; j < 8; j++)
	{
	  ip[i].ib[j].p    = (nodeptr)NULL;
	  ip[i].ib[j].freq = 0;
	}  
    }
  
 
  r = tr->nodep[(tr->nextnode)++];     
  q = findAnyTip(tr->start, tr->rdta->numsp);
  q = q->back;       
   

  if(adef->boot)
    {            
      for(i = 0; i < adef->multipleRuns; i++)
	{	  
	  if(i > 0)
	    {
	      makeboot(adef, tr);  	      
	      /*initModel(tr, tr->rdta, tr->cdta, adef);*/
	      evaluateGenericInitrav(tr, q->back);
	      printf("Init: %f\n", tr->likelihood);
	      /*if(adef->model == M_PROTGAMMA || adef->model == M_GTRGAMMA)    
		modOpt(tr, adef);    
	      else    
		quickOpt(tr, adef);  
	      */
	    }
	  
	  printf("Iteration %d: %f\n", i, tr->likelihood);      
      
	  count = 0;
      
	  for(j = 0; j < numberOfTipsForInsertion; j++) 
	    {
	      insertLikelihoods[j] = unlikely;
	      insertNodes[j]     = (nodeptr)NULL;
	    }

	  addTraverseRob(tr, r, q, inserts,  numberOfTipsForInsertion, &count, insertLikelihoods, insertNodes);
  
	  for(j = 0; j < numberOfTipsForInsertion; j++) 
	    {
	      insertInsertionPoint(&ip[j], insertNodes[j]);
	      /*printf("%s %f %d\n", tr->nameList[inserts[j]], insertLikelihoods[j], insertNodes[j]->number);*/
	    }	 
	}
    }
  else    
    { 
      evaluateGenericInitrav(tr, q->back);          
      
      count = 0;
      
      for(j = 0; j < numberOfTipsForInsertion; j++) 
	{
	  insertLikelihoods[j] = unlikely;
	  insertNodes[j]     = (nodeptr)NULL;
	}

      addTraverseRob(tr, r, q, inserts,  numberOfTipsForInsertion, &count, insertLikelihoods, insertNodes);
  
      for(j = 0; j < numberOfTipsForInsertion; j++) 
	{
	  insertInsertionPoint(&ip[j], insertNodes[j]);
	  printf("%s %f %d\n", tr->nameList[inserts[j]], insertLikelihoods[j], insertNodes[j]->number);
	}
    }

 
 

  for(j = 0; j <  numberOfTipsForInsertion; j++) 
    {
      i = 0;
      while(i < ip[j].max && ip[j].ib[i].p)
	{
	  double support = ((double)(ip[j].ib[i].freq)) / ((double) (adef->multipleRuns));
#ifdef WIN32
	  ip[j].ib[i].freq = (int)floor(0.5 + support * 100.0);
#else
	  ip[j].ib[i].freq = (int)(0.5 + support * 100.0);
#endif           
	  i++;
	}
    }
  
  tr->ip = ip;
  tr->numberOfTipsForInsertion = numberOfTipsForInsertion;

  {
    int *insertRefs, *frequencies;

    insertRefs  = (int*)malloc(sizeof(int) * numberOfTipsForInsertion);
    frequencies = (int*)malloc(sizeof(int) * numberOfTipsForInsertion);

    Tree2StringML(tr->tree_string, tr, q, insertRefs, frequencies);
    printf("%s \n", tr->tree_string);
  }


  exit(0);
}


#ifdef _MULTI_GENE

static void sumGAMMA(int tipCase, double *sumtable, double *x1_start, double *x2_start, double *tipVector,
		     char *tipX1, char *tipX2, int lower, int n)
{
  double *x1, *x2, *sum;
  int i;

  switch(tipCase)
    {
    case TIP_TIP:         
      for (i = lower; i < n; i++) 
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
      for (i = lower; i < n; i++) 
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


static void coreGTRGAMMA(int lower, int upper, double *sumtable, 
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

      
      
      inv_Li = 1.0 / inv_Li;
      
      dlnLidlz   *= inv_Li;
      d2lnLidlz2 *= inv_Li;
      
      dlnLdlz  += wrptr[i] * dlnLidlz;
      d2lnLdlz2 += wrptr[i] * (d2lnLidlz2 - dlnLidlz * dlnLidlz);
    }
  
  

  *d1 = dlnLdlz;
  *d2 = d2lnLdlz2;

  free(diagptable_start);
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
	    	 	 	    	  	  
	  coreGTRGAMMA(lower, upper, tr->sumBuffer, 
		       &dlnLdlz, &d2lnLdlz2, &(tr->EIGN_DNA[model * 3]), &(tr->gammaRates[model * 4]), lz, 
		       tr->cdta->aliaswgt);

	    	
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
		

static double evaluateGTRGAMMA(int *ex1, int *ex2, int *wptr,
			       double *x1_start, double *x2_start, double *EIGN, double *gammaRates, 
			       double *tipVector, double pz, 
			       char *tipX1, int lower, int n)
{
  double   sum = 0.0, z, lz, term, ki;    
  int     i;
  double  *diagptable, *x1, *x2; 
           
  assert(x1_start != x2_start);
  assert(ex1 != ex2);

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

	  

	  term = log(0.25 * term) + ex2[i] * log(minlikelihood);	 

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
	  
	  term = log(0.25 * term) + (ex1[i] + ex2[i]) * log(minlikelihood);
	  
	  sum += wptr[i] * term;
	}
            
      free(diagptable); 
    	
      return  sum;
    }
} 




static boolean updateMulti(tree *tr, nodeptr p, int model)
{       
  nodeptr  q;
  boolean smoothed;  
  double   z, z0;
  
  q = p->backs[model];   

  assert(p->backs[model]);
  assert(q->backs[model] == p);
  assert(q->z[model] == p->z[model]);

  z0 = q->z[model];
 
  /*printf("%d %d\n", p->number, q->number);*/

  tr->td[model].ti[0].pNumber = p->number;
  tr->td[model].ti[0].qNumber = q->number;	  
  tr->td[model].ti[0].qz[model] =  q->z[model];	  
  tr->td[model].count = 1;	       

  if(!p->xs[model])
    computeMultiTraversalInfo(p, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->rdta->numsp, model);  
  if(!q->xs[model])
    computeMultiTraversalInfo(q, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->rdta->numsp, model);

  newviewIterativePartition(tr, tr->partitionData[model].lower, tr->partitionData[model].upper, model);

  /*{
    double 
      result = 0.0,
      *x1_start = (double*)NULL, 
      *x2_start = (double*)NULL;
    int    
      *ex1 = (int*)NULL, 
      *ex2 = (int*)NULL;
    char *tip = (char*)NULL;   
    int pNumber, qNumber;
    double pz;

    pNumber = p->number;
    qNumber = q->number;
    pz      = q->z[model];

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
 	
    result = evaluateGTRGAMMA(ex1, ex2, tr->cdta->aliaswgt,
			      x1_start, x2_start, 
			      &(tr->EIGN_DNA[model * 3]), 
			      &(tr->gammaRates[model * 4]), 
			      &(tr->tipVectorDNA[model * 64]), 
			      pz, 
			      tip, tr->partitionData[model].lower, tr->partitionData[model].upper); 
    printf("RRRRRR %f\n", result);
  }
  z = z0;*/


  {
    double   
      *x1_start = (double*)NULL, 
      *x2_start = (double*)NULL;
    char 
      *tipX1 = (char*)NULL,
      *tipX2 = (char*)NULL;      
    int tipCase;
    int pNumber, qNumber; 


    pNumber = p->number;
    qNumber = q->number;    

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

    
    sumGAMMAMULT(tipCase, tr->sumBuffer, x1_start, x2_start, &(tr->tipVectorDNA[0]), tipX1, tipX2,  tr->model,
		 tr->partitionData[model].lower, tr->partitionData[model].upper);     
  }


  /*z = makenewzPartitionGeneric(tr, p, q, z0, newzpercycle, model);  */
  
  {
    double result = 0.0;

    topLevelMakenewzPartition(tr, tr->partitionData[model].lower, tr->partitionData[model].upper, model, 
			      z0, newzpercycle, &result);

    z = result;
  }
  
  smoothed = tr->smoothed;

  if (ABS(z - z0) > deltaz)  
    smoothed = FALSE;
  p->z[model] = q->z[model] = z;          
  
  tr->smoothed = smoothed;
  
  return TRUE;
}




static boolean smoothMulti (tree *tr, nodeptr p, int model, int *count)
{ 
  if(isTip(p->number, tr->rdta->numsp))
    assert(p->backs[model]);
  else
    assert(p->backs[model] && p->next->backs[model] && p->next->next->backs[model]);  

  if (! updateMulti(tr, p, model))
    {
      assert(0);
      return FALSE;
    }

  *count = *count + 1;

  /*if(*count == 22)
    exit(1);*/

  if(!isTip(p->number, tr->rdta->numsp)) 
   {              
     if(!smoothMulti(tr, p->next->backs[model], model, count))
       {
	 assert(0);
	 return FALSE;
       }        
       
     if (!smoothMulti(tr, p->next->next->backs[model], model, count))   
       {
	 assert(0);
	 return FALSE;      	
       }
    
     
     tr->td[model].count = 1;
     
     assert(!p->xs[model]);
     if(!p->xs[model])
       {
	 computeMultiTraversalInfo(p, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->rdta->numsp, model);
	 newviewIterativePartition(tr, tr->partitionData[model].lower, tr->partitionData[model].upper, model); 
       }
      
     /*updateMulti(tr, p, model);*/
     /*printf("HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH %d\n", *count);*/
   }
  
  return TRUE;
} 


static boolean smoothTreeMulti(tree *tr, int maxtimes, int model)
{ 
  nodeptr  p;    
  int count;

  p = tr->startVector[model];

  assert(p->backs[model]);
  
  assert(isTip(p->number, tr->rdta->numsp));
 
  while (--maxtimes >= 0) 
    {
      tr->smoothed = TRUE;

      count = 0;
      if (! smoothMulti(tr, p->backs[model], model, &count))
	{
	  exit(1);
	  return FALSE;       
	}
      if(!isTip(p->number, tr->rdta->numsp))
	assert(0);
      if(tr->smoothed)  
	break;      
    }

  return TRUE;
}


boolean treeEvaluateMulti(tree *tr, double smoothFactor)
{
  int model;

  double result = 0.0;

  for(model = 0; model < tr->NumberOfModels; model++)
    {
      nodeptr p, q;
      double erg;
      tr->modelNumber    = model;
      p = tr->startVector[model];
      q = p->backs[model];
      assert(p->backs[model]);
      assert(q->backs[model] == p);
     

      if (! smoothTreeMulti(tr, (int)((double)smoothings * smoothFactor), model))
	{
	  return FALSE;      
	}               
      
      assert(p->backs[model] == q && q->backs[model] == p);     
      
      tr->td[model].ti[0].pNumber = p->number;
      tr->td[model].ti[0].qNumber = q->number;	  
      tr->td[model].ti[0].qz[model] =  q->z[model];     
      tr->td[model].count = 1;

      assert(q->z[model] == p->z[model]);
      assert(isTip(p->number, tr->mxtips));

      /*if(!q->xs[model])*/
      computeFullMultiTraversalInfo(q, &(tr->td[model].ti[0]), &(tr->td[model].count), tr->rdta->numsp, model);

      result += (erg = evaluateIterativePartition(tr, 
						  tr->partitionData[model].lower, tr->partitionData[model].upper, model));
      /*printf("%d %f\n", model, erg);*/
    }
   
  tr->likelihood = result;
  /*evaluateGenericInitrav(tr, tr->start); 
    printf("%f\n", tr->likelihood);*/

    return TRUE;
  } /* treeEvaluate */


#endif
