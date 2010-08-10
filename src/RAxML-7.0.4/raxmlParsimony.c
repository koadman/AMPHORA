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
#include <sys/times.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>  
#endif

#include <limits.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "axml.h"

extern double masterTime;
extern const int protTipParsimonyValue[23];

#ifdef _USE_PTHREADS
extern int NumberOfThreads;
extern int *reductionBufferParsimony;
#endif


/********************************DNA FUNCTIONS *****************************************************************/

static void computeTraversalInfoParsimony(nodeptr p, traversalInfo *ti, int *counter, int maxTips)
{
  if(isTip(p->number, maxTips))
    return;

  {         
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
		  computeTraversalInfoParsimony(r, ti, counter, maxTips);
		if (! p->x) 
		  getxnode(p);	
	      }
	    	   
	    ti[*counter].tipCase = TIP_INNER; 
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;	   	    
	    *counter = *counter + 1;
	  }
	else
	  {	 

	    while ((! p->x) || (! q->x) || (! r->x)) 
	      {
		if (! q->x) 
		  computeTraversalInfoParsimony(q, ti, counter, maxTips);
		if (! r->x) 
		  computeTraversalInfoParsimony(r, ti, counter, maxTips);
		if (! p->x) 
		  getxnode(p);	
	      }
   
	    ti[*counter].tipCase = INNER_INNER; 
	    ti[*counter].pNumber = p->number;
	    ti[*counter].qNumber = q->number;
	    ti[*counter].rNumber = r->number;	   
	    *counter = *counter + 1;
	  }
      }    
  }
}



static void newviewParsimonyDNA(traversalInfo *tInfo, char *right, char *left, 
				parsimonyVector *rightVector, parsimonyVector *leftVector, 
				parsimonyVector *thisVector, int lower, int upper)
{
  int i, ts; 
  int le, ri, t; 

  switch(tInfo->tipCase)
    {
    case TIP_TIP:     
      for(i = lower; i < upper; i++)
	{	 
	  le = left[i];
	  ri = right[i];
	  
	  t = le & ri;
	  
	  ts = 0;
	  
	  if(!t)
	    {
	      t = le | ri;
	      ts = 1;
	    }
	  
	  thisVector[i].parsimonyScore = ts;	       	    
	  thisVector[i].parsimonyState = t;	     
	}
      break;
    case TIP_INNER:      
      for(i = lower; i < upper; i++)
	{
	  le = left[i];
	  ri = rightVector[i].parsimonyState;
	  
	  t = le & ri;
	  
	  ts = rightVector[i].parsimonyScore;
	  
	  if(!t)
	    {
	      t = le | ri;
	      ts++;
	    }
	  
	  thisVector[i].parsimonyScore = ts;	       	    
	  thisVector[i].parsimonyState = t;	    	     
	}
      break;      
    case INNER_INNER:    
      for(i = lower; i < upper; i++)
	{
	  le = leftVector[i].parsimonyState;
	  ri = rightVector[i].parsimonyState;
	  
	  t = le & ri;
	  
	  ts = rightVector[i].parsimonyScore + leftVector[i].parsimonyScore;
	  
	  if(!t)
	    {
	      t = le | ri;
	      ts++;
	    }
	  
	  thisVector[i].parsimonyScore = ts;	       	    
	  thisVector[i].parsimonyState = t;	 
	  
	}
      break;
    default:
      assert(0);
    }
}

static void newviewParsimonyPROT(traversalInfo *tInfo, char *right, char *left, 
				 parsimonyVector *rightVector, parsimonyVector *leftVector, 
				 parsimonyVector *thisVector, int lower, int upper)
{
  int i, ts; 
  int le, ri, t;

  switch(tInfo->tipCase)
    {
    case TIP_TIP:
      for(i = lower; i < upper; i++)
	{	 
	  le = protTipParsimonyValue[left[i]];
	  ri = protTipParsimonyValue[right[i]];	  
	  
	  t = le & ri;
	  
	  ts = 0;
	  
	  if(!t)
	    {
	      t = le | ri;
	      ts = 1;
	    }
	  
	  thisVector[i].parsimonyScore = ts;	       	    
	  thisVector[i].parsimonyState = t;	     
	}
      break;
    case TIP_INNER:
      for(i = lower; i < upper; i++)
	{	  
	  le = protTipParsimonyValue[left[i]];       

	  ri = rightVector[i].parsimonyState;
	  
	  t = le & ri;
	  
	  ts = rightVector[i].parsimonyScore;
	  
	  if(!t)
	    {
	      t = le | ri;
	      ts++;
	    }
	  
	  thisVector[i].parsimonyScore = ts;	       	    
	  thisVector[i].parsimonyState = t;	    	     
	}
      break;      
    case INNER_INNER:
      for(i = lower; i < upper; i++)
	{
	  le = leftVector[i].parsimonyState;
	  ri = rightVector[i].parsimonyState;
	  
	  t = le & ri;
	  
	  ts = rightVector[i].parsimonyScore + leftVector[i].parsimonyScore;
	  
	  if(!t)
	    {
	      t = le | ri;
	      ts++;
	    }
	  
	  thisVector[i].parsimonyScore = ts;	       	    
	  thisVector[i].parsimonyState = t;	 
	  
	}
      break;
    default:
      assert(0);
    }
}




static int evalDNA(char *right, parsimonyVector *rightVector,parsimonyVector *leftVector, int lower, int upper, int *wptr)
{
  int i, sum, acc = 0;
  int le, ri;

  if(right)
    {
      for(i = lower; i < upper; i++)
	{
	  le = leftVector[i].parsimonyState;
	  ri = right[i];
	  
	  sum = leftVector[i].parsimonyScore;
	  
	  if(!(le & ri))
	    sum++;
	  
	  acc += wptr[i] * sum;	  		 	       
	} 
    }
  else
    {
      for(i = lower; i < upper; i++)
	{
	  le = leftVector[i].parsimonyState;
	  ri = rightVector[i].parsimonyState;	     	
	  
	  sum = rightVector[i].parsimonyScore + leftVector[i].parsimonyScore;
	     
	  if(!(le & ri))
	    sum++;
	  
	  acc += wptr[i] * sum;     
	}     
    }

  return acc;
}

static int evalPROT(char *right, parsimonyVector *rightVector, parsimonyVector *leftVector, int lower, int upper, int *wptr)
{
  int i, sum, acc = 0;
  int le, ri;

  if(right)
    {
      for(i = lower; i < upper; i++)
	{
	  le = leftVector[i].parsimonyState;
	  ri = protTipParsimonyValue[right[i]];
	  
	  
	  
	  sum = leftVector[i].parsimonyScore;
	  
	  if(!(le & ri))
	    sum++;
	  
	  acc += wptr[i] * sum;	  		 	       
	} 
    }
  else
    {
      for(i = lower; i < upper; i++)
	{
	  le = leftVector[i].parsimonyState;
	  ri = rightVector[i].parsimonyState;	     	
	  
	  sum = rightVector[i].parsimonyScore + leftVector[i].parsimonyScore;
	     
	  if(!(le & ri))
	    sum++;
	  
	  acc += wptr[i] * sum;     
	}     
    }


  return acc;
}

#ifdef _LOCAL_DATA

void newviewParsimonyIterative(tree *localTree, int startIndex, int endIndex)
{  
  traversalInfo *ti   = localTree->td[0].ti;
  int i;  

  for(i = 1; i < localTree->td[0].count; i++)
    {    
      traversalInfo *tInfo = &ti[i];
      char 
	*right = (char*)NULL, 
	*left  = (char*)NULL;      
      parsimonyVector 
	*rightVector = (parsimonyVector *)NULL, 
	*leftVector  = (parsimonyVector *)NULL, 
	*thisVector  = (parsimonyVector *)NULL;

      switch(tInfo->tipCase)
	{
	case TIP_TIP:
	  left       = &localTree->strided_yVector[tInfo->qNumber][startIndex];
	  right      = &localTree->strided_yVector[tInfo->rNumber][startIndex];
	  thisVector = &(localTree->parsimonyData[localTree->mySpan * (tInfo->pNumber - localTree->mxtips - 1)]);
	  break;
	case TIP_INNER:
	  left        = &localTree->strided_yVector[tInfo->qNumber][startIndex];
	  rightVector = &(localTree->parsimonyData[localTree->mySpan * (tInfo->rNumber - localTree->mxtips - 1)]);
	  thisVector  = &(localTree->parsimonyData[localTree->mySpan * (tInfo->pNumber - localTree->mxtips - 1)]); 
	  break;
	case INNER_INNER:
	  leftVector  = &(localTree->parsimonyData[localTree->mySpan * (tInfo->qNumber - localTree->mxtips - 1)]);
	  rightVector = &(localTree->parsimonyData[localTree->mySpan * (tInfo->rNumber - localTree->mxtips - 1)]);
	  thisVector  = &(localTree->parsimonyData[localTree->mySpan * (tInfo->pNumber - localTree->mxtips - 1)]);
	  break;
	default:
	  assert(0);
	}
      
      if(localTree->mixedData)
	{	  
	  int i;
	  
	  for(i = 0; i < localTree->NumberOfModels; i++)
	    {	     
	      int l = localTree->partitionData[i].lower;
	      int u = localTree->partitionData[i].upper;
	      
	      switch(localTree->partitionData[i].dataType)
		{       	
		case AA_DATA: 
		  newviewParsimonyPROT(tInfo, right, left, rightVector, leftVector, thisVector, l, u);
		  break;
		case DNA_DATA:
		  newviewParsimonyDNA(tInfo, right, left, rightVector, leftVector, thisVector, l, u);
		  break;
		default:
		  assert(0);
		}   	
	    }    
	}
      else
	{	 
	  switch(localTree->partitionData[0].dataType)
	    {
	    case AA_DATA:            
	      newviewParsimonyPROT(tInfo, right, left, rightVector, leftVector, thisVector, 0, (endIndex- startIndex));	 
	      break;
	    case DNA_DATA:	     
	      newviewParsimonyDNA(tInfo, right, left, rightVector, leftVector, thisVector, 0, (endIndex - startIndex));	 
	      break;
	    default:
	      assert(0);
	    }           
	}           
    }
}

int evaluateParsimonyIterative(tree *localTree, int lower, int upper)
{
  int pNumber, qNumber, result;   
  char *right = (char *)NULL; 
  parsimonyVector 
    *rightVector = (parsimonyVector *)NULL, 
    *leftVector  = (parsimonyVector *)NULL;      

  pNumber = localTree->td[0].ti[0].pNumber;
  qNumber = localTree->td[0].ti[0].qNumber;

 
  newviewParsimonyIterative(localTree, lower, upper);

  if(isTip(pNumber, localTree->mxtips) || isTip(qNumber, localTree->mxtips))
    {	        	    
      if(isTip(qNumber, localTree->mxtips))
	{	
	  leftVector = &(localTree->parsimonyData[localTree->mySpan * (pNumber - localTree->mxtips - 1)]);
	  right = &localTree->strided_yVector[qNumber][lower];
	}           
      else
	{
	  leftVector = &(localTree->parsimonyData[localTree->mySpan * (qNumber - localTree->mxtips - 1)]);
	  right = &localTree->strided_yVector[pNumber][lower];
	}
    }
  else
    {           
      leftVector  = &(localTree->parsimonyData[localTree->mySpan * (pNumber - localTree->mxtips - 1)]);
      rightVector = &(localTree->parsimonyData[localTree->mySpan * (qNumber - localTree->mxtips - 1)]);         
    }
  
  if(localTree->mixedData)
    {     
      int i, partialResult = 0;

      for(i = 0; i < localTree->NumberOfModels; i++)
	{	  
	  int l = localTree->partitionData[i].lower;
	  int u = localTree->partitionData[i].upper;

	  switch(localTree->partitionData[i].dataType)
	    {       	
	    case AA_DATA: 
	      partialResult += evalPROT(right, rightVector, leftVector, l, u, localTree->strided_aliaswgt);
	      break;
	    case DNA_DATA:	      
	      partialResult += evalDNA(right, rightVector, leftVector, l, u, localTree->strided_aliaswgt);	      	    
	      break;
	    default:
	      assert(0);
	    } 	 
	}
      result = partialResult;
    }
  else
    {
      switch(localTree->partitionData[0].dataType)
	{
	case AA_DATA: 
	  result = evalPROT(right, rightVector, leftVector, 0, (upper - lower), &(localTree->strided_aliaswgt[lower]));
	  break;
	case DNA_DATA:
	  result = evalDNA(right, rightVector, leftVector, 0, (upper - lower), &(localTree->strided_aliaswgt[lower]));
	  break;
	default:
	  assert(0);
	}    
    }
  
  return result;
}

#else

void newviewParsimonyIterative(tree *tr, int startIndex, int endIndex)
{  
  traversalInfo *ti   = tr->td[0].ti;
  int i;

  for(i = 1; i < tr->td[0].count; i++)
    {    
      traversalInfo *tInfo = &ti[i];
      char 
	*right = (char*)NULL, 
	*left  = (char*)NULL;      
      parsimonyVector 
	*rightVector = (parsimonyVector *)NULL, 
	*leftVector  = (parsimonyVector *)NULL, 
	*thisVector  = (parsimonyVector *)NULL;

      switch(tInfo->tipCase)
	{
	case TIP_TIP:
	  left       = tr->yVector[tInfo->qNumber];
	  right      = tr->yVector[tInfo->rNumber];
	  thisVector = &(tr->parsimonyData[tr->parsimonyLength * (tInfo->pNumber - tr->mxtips - 1)]);
	  break;
	case TIP_INNER:
	  left        = tr->yVector[tInfo->qNumber];
	  rightVector = &(tr->parsimonyData[tr->parsimonyLength * (tInfo->rNumber - tr->mxtips - 1)]);
	  thisVector  = &(tr->parsimonyData[tr->parsimonyLength * (tInfo->pNumber - tr->mxtips - 1)]); 
	  break;
	case INNER_INNER:
	  leftVector  = &(tr->parsimonyData[tr->parsimonyLength * (tInfo->qNumber - tr->mxtips - 1)]);
	  rightVector = &(tr->parsimonyData[tr->parsimonyLength * (tInfo->rNumber - tr->mxtips - 1)]);
	  thisVector  = &(tr->parsimonyData[tr->parsimonyLength * (tInfo->pNumber - tr->mxtips - 1)]);
	  break;
	default:
	  assert(0);
	}
      
      if(tr->mixedData)
	{	  
	  int i;
	  
	  for(i = 0; i < tr->NumberOfModels; i++)
	    {	     
	      int l = tr->partitionData[i].lower;
	      int u = tr->partitionData[i].upper;
	      
	      switch(tr->partitionData[i].dataType)
		{       	
		case AA_DATA: 
		  newviewParsimonyPROT(tInfo, right, left, rightVector, leftVector, thisVector, l, u);
		  break;
		case DNA_DATA:
		  newviewParsimonyDNA(tInfo, right, left, rightVector, leftVector, thisVector, l, u);
		  break;
		default:
		  assert(0);
		}   	
	    }    
	}
      else
	{
	  switch(tr->partitionData[0].dataType)
	    {
	    case AA_DATA:            
	      newviewParsimonyPROT(tInfo, right, left, rightVector, leftVector, thisVector, startIndex, endIndex);	 
	      break;
	    case DNA_DATA:
	      newviewParsimonyDNA(tInfo, right, left, rightVector, leftVector, thisVector, startIndex, endIndex);	 
	      break;
	    default:
	      assert(0);
	    }           
	}           
    }
}

int evaluateParsimonyIterative(tree *tr, int lower, int upper)
{
  int pNumber, qNumber, result;   
  char *right = (char *)NULL; 
  parsimonyVector 
    *rightVector = (parsimonyVector *)NULL, 
    *leftVector  = (parsimonyVector *)NULL;   
  
  pNumber = tr->td[0].ti[0].pNumber;
  qNumber = tr->td[0].ti[0].qNumber;

  newviewParsimonyIterative(tr, lower, upper);

  if(isTip(pNumber, tr->rdta->numsp) || isTip(qNumber, tr->rdta->numsp))
    {	        	    
      if(isTip(qNumber, tr->rdta->numsp))
	{	
	  leftVector = &(tr->parsimonyData[tr->parsimonyLength * (pNumber - tr->mxtips - 1)]);
	  right = tr->yVector[qNumber];
	}           
      else
	{
	  leftVector = &(tr->parsimonyData[tr->parsimonyLength * (qNumber - tr->mxtips - 1)]);
	  right = tr->yVector[pNumber];
	}
    }
  else
    {           
      leftVector  = &(tr->parsimonyData[tr->parsimonyLength * (pNumber - tr->mxtips - 1)]);
      rightVector = &(tr->parsimonyData[tr->parsimonyLength * (qNumber - tr->mxtips - 1)]);         
    }
  
  if(tr->mixedData)
    {     
      int i, partialResult = 0;

      for(i = 0; i < tr->NumberOfModels; i++)
	{	  
	  int l = tr->partitionData[i].lower;
	  int u = tr->partitionData[i].upper;

	  switch(tr->partitionData[i].dataType)
	    {       	
	    case AA_DATA: 
	      partialResult += evalPROT(right, rightVector, leftVector, l, u, tr->cdta->aliaswgt);
	      break;
	    case DNA_DATA:
	      partialResult += evalDNA(right, rightVector, leftVector, l, u, tr->cdta->aliaswgt);
	      break;
	    default:
	      assert(0);
	    } 	 
	}
      result = partialResult;
    }
  else
    {
      switch(tr->partitionData[0].dataType)
	{
	case AA_DATA: 
	  result = evalPROT(right, rightVector, leftVector, lower, upper, tr->cdta->aliaswgt);
	  break;
	case DNA_DATA:
	  result = evalDNA(right, rightVector, leftVector, lower, upper, tr->cdta->aliaswgt);
	  break;
	default:
	  assert(0);
	}    
    }

  return result;
}


#endif

static int evaluateParsimony(tree *tr, nodeptr p)
{
  int result;
  nodeptr q = p->back;
  tr->td[0].ti[0].pNumber = p->number;
  tr->td[0].ti[0].qNumber = q->number;

  tr->td[0].count = 1;

  if(!p->x)
    computeTraversalInfoParsimony(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp);
  if(!q->x)
    computeTraversalInfoParsimony(q, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp); 

#ifdef _USE_PTHREADS
  {
    int i;    

    masterBarrier(THREAD_EVALUATE_PARSIMONY, tr);    
    
    for(i = 0, result = 0; i < NumberOfThreads; i++)
      result += reductionBufferParsimony[i]; 
  }
#else
  result = evaluateParsimonyIterative(tr, 0, tr->parsimonyLength);
#endif

  return result;
}


static void newviewParsimony(tree *tr, nodeptr  p)
{     
  if(isTip(p->number, tr->rdta->numsp))
    return;
  
  tr->td[0].count = 1;
  computeTraversalInfoParsimony(p, &(tr->td[0].ti[0]), &(tr->td[0].count), tr->rdta->numsp);              

  if(tr->td[0].count > 1)
    {
#ifdef _USE_PTHREADS
      masterBarrier(THREAD_NEWVIEW_PARSIMONY, tr); 
#else
      newviewParsimonyIterative(tr, 0, tr->parsimonyLength);    
#endif
    }
}





/****************************************************************************************************************************************/

static void initravParsimonyNormal(tree *tr, nodeptr p)
{
  nodeptr  q;
  
  if (! isTip(p->number, tr->rdta->numsp)) 
    {
      q = p->next;
      
      do 
	{
	  initravParsimonyNormal(tr, q->back);
	  q = q->next;	
	} 
      while (q != p);
      
      newviewParsimony(tr, p);	      
    }
}


static void initravParsimony(tree *tr, nodeptr p, int *constraintVector)
{
  nodeptr  q;
  
  if (! isTip(p->number, tr->rdta->numsp)) 
    {    
      q = p->next;
      
      do 
	{
	  initravParsimony(tr, q->back, constraintVector);
	  q = q->next;	
	} 
      while (q != p);
      
      newviewParsimony(tr, p);	      
    }
  else
    constraintVector[p->number] = 1;
}

static void insertParsimony (tree *tr, nodeptr p, nodeptr q)
{
  nodeptr  r;
  
  r = q->back;
  
  hookupDefault(p->next,       q, tr->numBranches);
  hookupDefault(p->next->next, r, tr->numBranches); 
   
  newviewParsimony(tr, p);     
} 

static void insertRandom (nodeptr p, nodeptr q, int numBranches)
{
  nodeptr  r;
  
  r = q->back;
  
  hookupDefault(p->next,       q, numBranches);
  hookupDefault(p->next->next, r, numBranches); 
} 


static nodeptr buildNewTip (tree *tr, nodeptr p)
{ 
  nodeptr  q;

  q = tr->nodep[(tr->nextnode)++];
  hookupDefault(p, q, tr->numBranches);
  return  q;
} 

static void buildSimpleTree (tree *tr, int ip, int iq, int ir)
{    
  nodeptr  p, s;
  int  i;
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  hookupDefault(p, tr->nodep[iq], tr->numBranches);
  s = buildNewTip(tr, tr->nodep[ir]);
  insertParsimony(tr, s, p);
}


static void buildSimpleTreeRandom (tree *tr, int ip, int iq, int ir)
{    
  nodeptr  p, s;
  int  i;
  
  i = MIN(ip, iq);
  if (ir < i)  i = ir; 
  tr->start = tr->nodep[i];
  tr->ntips = 3;
  p = tr->nodep[ip];
  hookupDefault(p, tr->nodep[iq], tr->numBranches);
  s = buildNewTip(tr, tr->nodep[ir]);
  insertRandom( s, p, tr->numBranches);
}

int checker(tree *tr, nodeptr p)
{
  int group = tr->constraintVector[p->number];

  if(isTip(p->number, tr->rdta->numsp))
    {
      group = tr->constraintVector[p->number];
      return group;
    }
  else
    {
      if(group != -9) 
	return group;

      group = checker(tr, p->next->back);
      if(group != -9) 
	return group;

      group = checker(tr, p->next->next->back);
      if(group != -9) 
	return group;

      return -9;
    }
}


static void testInsertParsimony (tree *tr, nodeptr p, nodeptr q)
{ 
  int mp;
  boolean doIt = TRUE;
  nodeptr  r = q->back;
 
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
      insertParsimony(tr, p, q);   
      mp = evaluateParsimony(tr, p->next->next);          
      
      if(mp < tr->bestParsimony)
	{
	  tr->bestParsimony = mp;
	  tr->insertNode = q;
	  tr->removeNode = p;
	}
      
      hookupDefault(q, r, tr->numBranches);
      p->next->next->back = p->next->back = (nodeptr) NULL;

    }      

  return;
} 


static void restoreTreeParsimony(tree *tr, nodeptr p, nodeptr q)
{ 
  insertParsimony(tr, p, q);  

  if(! isTip(p->number, tr->rdta->numsp) && isTip(q->number, tr->rdta->numsp))
    {
      while ((! p->x)) 
	{
	  if (! (p->x))
	    newviewParsimony(tr, p);		     
	}
    }
  if(isTip(p->number, tr->rdta->numsp) && ! isTip(q->number, tr->rdta->numsp))
    {
      while ((! q->x)) 
	{		  
	  if (! (q->x)) 
	    newviewParsimony(tr, q);
	}
    }
  if(! isTip(p->number, tr->rdta->numsp) && ! isTip(q->number, tr->rdta->numsp))
    {
      while ((! p->x) || (! q->x)) 
	{
	  if (! (p->x))
	    newviewParsimony(tr, p);
	  if (! (q->x))
	    newviewParsimony(tr, q);
	}
    }	
}


static int markBranches(nodeptr *branches, nodeptr p, int *counter, int numsp)
{
  if(isTip(p->number, numsp))
    return 0;
  else
    {
      branches[*counter] = p->next;
      branches[*counter + 1] = p->next->next;
      
      *counter = *counter + 2;
      
      return ((2 + markBranches(branches, p->next->back, counter, numsp) + 
	       markBranches(branches, p->next->next->back, counter, numsp)));
    }
}

static void addTraverseParsimony (tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav)
{       
  if (--mintrav <= 0)           
    testInsertParsimony(tr, p, q);	            

  if ((! isTip(q->number, tr->rdta->numsp)) && (--maxtrav > 0)) 
    {	
      addTraverseParsimony(tr, p, q->next->back, mintrav, maxtrav);	
      addTraverseParsimony(tr, p, q->next->next->back, mintrav, maxtrav);              	
    }
}


nodeptr findAnyTip(nodeptr p, int numsp)
{ 
  return  isTip(p->number, numsp) ? p : findAnyTip(p->next->back, numsp);
} 


int randomInt(int n)
{
  return rand() %n;
}

void makePermutation(int *perm, int n, analdef *adef)
{    
  int  i, j, k;

#ifdef PARALLEL
   srand((unsigned int) gettimeSrand()); 
#else
  if(adef->parsimonySeed == 0)   
    srand((unsigned int) gettimeSrand());          
#endif

  for (i = 1; i <= n; i++)    
    perm[i] = i;               

  for (i = 1; i <= n; i++) 
    {
#ifdef PARALLEL      
      k        = randomInt(n + 1 - i);
#else
     if(adef->parsimonySeed == 0) 
       k        = randomInt(n + 1 - i);
     else
       k =  (int)((double)(n + 1 - i) * randum(&adef->parsimonySeed));
#endif
     assert(i + k <= n);

     j        = perm[i];
     perm[i]     = perm[i + k];
     perm[i + k] = j; 
    }

  /*  for(i = 1; i <= n; i++)
    printf("%d ", perm[i]);
    printf("\n");*/

}

static void initravDISTParsimony (tree *tr, nodeptr p, int distance)
{
  nodeptr  q;

  if (! isTip(p->number, tr->rdta->numsp) && distance > 0) 
    {      
      q = p->next;      
      do 
	{
	  initravDISTParsimony(tr, q->back, --distance);	
	  q = q->next;	
	} 
      while (q != p);
      
      
      newviewParsimony(tr, p);	      
    } 
}


static nodeptr  removeNodeParsimony (nodeptr p, int numBranches)
{ 
  nodeptr  q, r;         

  q = p->next->back;
  r = p->next->next->back;   
    
  hookupDefault(q, r, numBranches);

  p->next->next->back = p->next->back = (node *) NULL;
  
  return  q;
}



boolean tipHomogeneityChecker(tree *tr, nodeptr p, int grouping)
{
  if(isTip(p->number, tr->rdta->numsp))
    {
      if(tr->constraintVector[p->number] != grouping) 
	return FALSE;
      else 
	return TRUE;
    }
  else
    {   
      return  (tipHomogeneityChecker(tr, p->next->back, grouping) && tipHomogeneityChecker(tr, p->next->next->back,grouping));      
    }
}

static int rearrangeParsimony(tree *tr, nodeptr p, int mintrav, int maxtrav)  
{   
  nodeptr  p1, p2, q, q1, q2;
  int      mintrav2;
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

  if (!isTip(p->number, tr->rdta->numsp) && doP) 
    {     
      p1 = p->next->back;
      p2 = p->next->next->back;
      
      if (! isTip(p1->number, tr->rdta->numsp) || ! isTip(p2->number, tr->rdta->numsp)) 
	{	  	  
	  removeNodeParsimony(p, tr->numBranches);
	  
	  if (! isTip(p1->number, tr->rdta->numsp)) 
	    {
	      addTraverseParsimony(tr, p, p1->next->back, mintrav, maxtrav);         
	      addTraverseParsimony(tr, p, p1->next->next->back, mintrav, maxtrav);          
	    }

	  if (! isTip(p2->number, tr->rdta->numsp)) 
	    {
	      addTraverseParsimony(tr, p, p2->next->back, mintrav, maxtrav);
	      addTraverseParsimony(tr, p, p2->next->next->back, mintrav, maxtrav);          
	    }
	    
	   
	  hookupDefault(p->next,       p1, tr->numBranches); 
	  hookupDefault(p->next->next, p2, tr->numBranches);	   	    	    
	  initravDISTParsimony(tr, p, 1);   
	}
    }  
       
  if (! isTip(q->number, tr->rdta->numsp) && maxtrav > 0 && doQ) 
    {
      q1 = q->next->back;
      q2 = q->next->next->back;

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

	  removeNodeParsimony(q, tr->numBranches);
	  
	  mintrav2 = mintrav > 2 ? mintrav : 2;

	  if (! isTip(q1->number, tr->rdta->numsp)) 
	    {
	      addTraverseParsimony(tr, q, q1->next->back, mintrav2 , maxtrav);
	      addTraverseParsimony(tr, q, q1->next->next->back, mintrav2 , maxtrav);         
	    }

	  if (! isTip(q2->number, tr->rdta->numsp)) 
	    {
	      addTraverseParsimony(tr, q, q2->next->back, mintrav2 , maxtrav);
	      addTraverseParsimony(tr, q, q2->next->next->back, mintrav2 , maxtrav);          
	    }	   
	   
	  hookupDefault(q->next,       q1, tr->numBranches); 
	  hookupDefault(q->next->next, q2, tr->numBranches);
	   
	  initravDISTParsimony(tr, q, 1); 	   
	}
    }

  return 1;
} 


static void restoreTreeRearrangeParsimony(tree *tr)
{    
  removeNodeParsimony(tr->removeNode, tr->numBranches);  
  restoreTreeParsimony(tr, tr->removeNode, tr->insertNode);  
}

static void allocNodexParsimony(tree *tr)
{
  nodeptr  p;
  int  i;  

#ifdef _LOCAL_DATA
  masterBarrier(THREAD_PREPARE_PARSIMONY, tr);  
#else
  tr->parsimonyData = (parsimonyVector *)malloc(sizeof(parsimonyVector) * tr->mxtips * tr->parsimonyLength); 
#endif
  
  for (i = tr->mxtips + 1; (i <= 2*(tr->mxtips) - 2); i++) 
    {       
      p = tr->nodep[i];    

      p->x             = 1;
      p->next->x       = 0;
      p->next->next->x = 0;
    }
}


static void freeNodexParsimony (tree *tr)
{
  nodeptr  p;
  int  i;  
  
#ifdef _LOCAL_DATA
  masterBarrier(THREAD_FINISH_PARSIMONY, tr);   
#else
  free(tr->parsimonyData); 
#endif
  
  for (i = tr->mxtips + 1; (i <= 2*(tr->mxtips) - 2); i++) 
    {
      p = tr->nodep[i];
      while(!p->x)
	p = p->next;    
      p->x = 0;
    }
}

#ifdef WIN32
static void switchTipEntries(int number, int position1, int position2, 
			     char *y0, int originalCrunchedLength, int numsp)
#else
static inline void switchTipEntries(int number, int position1, int position2, 
				    char *y0, int originalCrunchedLength, int numsp)
#endif
{
  char buf;
  char *ref = &y0[originalCrunchedLength * (number - 1)];

  assert(number <= numsp && number > 0);
  assert(position1 <  originalCrunchedLength && position2 < originalCrunchedLength);
  assert(position1 >= 0 && position2 >= 0);

  buf = ref[position1];
  ref[position1] = ref[position2];
  ref[position2] = buf;
}





static void sortInformativeSites(tree *tr, int *informative)
{
  int i, l, j;

  for(i = 0; i < tr->rdta->numsp; i++)
    {
      char *yPos    = &(tr->rdta->y0[tr->originalCrunchedLength * i]);        
      
      for(j = 0, l = 0; j < tr->cdta->endsite; j++)
	{	
	  if(informative[j])	  
	    {	     
	      yPos[l++] = yPos[j];
	    }
	}               
    }
    
  for(j = 0, l = 0; j < tr->cdta->endsite; j++)
    {     
      if(informative[j])	
	{
	  tr->cdta->aliaswgt[l]     = tr->cdta->aliaswgt[j];	
	  tr->model[l]              = tr->model[j];	 
	  tr->dataVector[l]           = tr->dataVector[j];
	  l++;
	}
    }    
}

/* TODO should re-visit this one day, not sure that I am getting all uninformative sites */

static void determineUninformativeSites(tree *tr, int *informative)
{
  int i, j;
  int check[23];
  int nucleotide;
  int informativeCounter;       
  int number = 0;

  /* 
     Not all characters are useful in constructing a parsimony tree. 
     Invariant characters, those that have the same state in all taxa, 
     are obviously useless and are ignored by the method. Characters in 
     which a state occurs in only one taxon are also ignored. 
     All these characters are called parsimony uninformative.
  */

  for(i = 0; i < tr->cdta->endsite; i++)
    {
      switch(tr->dataVector[i])
	{
	case AA_DATA: 
	  for(j = 0; j < 23; j++)
	    check[j] = 0;
	  
	  for(j = 1; j <= tr->mxtips; j++)
	    {	     
	      nucleotide = tr->yVector[tr->nodep[j]->number][i];	      
	      check[nucleotide] = check[nucleotide] + 1;
	      assert(nucleotide < 23 && nucleotide >= 0);
	    }
	  
	  informativeCounter = 0;
	  
	  for(j = 0; j < 22; j++)
	    {
	      if(check[j] > 0)
		informativeCounter++;    
	    } 
	  
	  if(informativeCounter <= 1)
	    {
	      informative[i] = 0;
	      number++;
	    }
	  else
	    {
	      boolean isInformative = FALSE;
	      for(j = 0; j < 22 && !(isInformative); j++)
		{
		  if(check[j] > 1)
		    isInformative = TRUE;
		} 

	      if(isInformative)
		informative[i] = 1; 
	      else
		{
		  informative[i] = 0; 
		  number++;
		}	  
	    }           
	  break;
	case DNA_DATA:     
	  for(j = 1; j < 16; j++)
	    check[j] = 0;
	  
	  for(j = 1; j <= tr->mxtips; j++)
	    {	   
	      nucleotide = tr->yVector[tr->nodep[j]->number][i];	    
	      check[nucleotide] =  check[nucleotide] + 1;
	      assert(nucleotide < 16 && nucleotide >= 0);
	    }
	  
	  informativeCounter = 0;
	  
	  for(j = 1; j < 15; j++)
	    {
	      if(check[j] > 0)
		informativeCounter++;    
	    } 
	  
	  if(informativeCounter <= 1)
	    {
	      informative[i] = 0;
	      number++;
	    }
	  else
	    {
	      boolean isInformative = FALSE;
	      for(j = 1; j < 15 && !(isInformative); j++)
		{
		  if(check[j] > 1)
		    isInformative = TRUE;
		} 

	      if(isInformative)
		informative[i] = 1; 
	      else
		{
		  informative[i] = 0; 
		  number++;
		}	  
	    }           	 
	  break;
	default:
	  assert(0);
	}
    }
 
  sortInformativeSites(tr, informative);

  /*printf("Uninformative Patterns: %d\n", number);*/
  
  tr->parsimonyLength = tr->cdta->endsite - number;
}



void makeRandomTree(tree *tr, analdef *adef)
{  
  nodeptr p, f, randomBranch;    
  int nextsp;
  int *perm, branchCounter;
  nodeptr *branches;
  
  branches = (nodeptr *)malloc(sizeof(nodeptr) * (2 * tr->mxtips));
  perm = (int *)malloc((tr->mxtips + 1) * sizeof(int));                         
  
  makePermutation(perm, tr->mxtips, adef);              
  
  tr->ntips = 0;       	       
  tr->nextnode = tr->mxtips + 1;    
  
  buildSimpleTreeRandom(tr, perm[1], perm[2], perm[3]);
  
  while (tr->ntips < tr->mxtips) 
    {	       
      tr->bestParsimony = INT_MAX;
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];
      
      /*printf("ADDING SPECIES %d\n", nextsp);*/
      
      buildNewTip(tr, p);  	
      
      f = findAnyTip(tr->start, tr->rdta->numsp);
      f = f->back;
      
      branchCounter = 1;
      branches[0] = f;
      markBranches(branches, f, &branchCounter, tr->rdta->numsp);

      assert(branchCounter == ((2 * (tr->ntips - 1)) - 3));
      
      randomBranch = branches[randomInt(branchCounter)];
      
      insertRandom(p->back, randomBranch, tr->numBranches);
      
    }
  free(perm);            
  free(branches);
}


static void reorderNodes(tree *tr, nodeptr *np, nodeptr p, int *count)
{
  int i, found = 0;

  if(isTip(p->number, tr->rdta->numsp))    
    return;
  else
    {              
      for(i = tr->mxtips + 1; (i <= (tr->mxtips + tr->mxtips - 1)) && (found == 0); i++)
	{
	  if (p == np[i] || p == np[i]->next || p == np[i]->next->next)
	    {
	      if(p == np[i])			       
		tr->nodep[*count + tr->mxtips + 1] = np[i];		 		
	      else
		{
		  if(p == np[i]->next)		  
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next;		     	   
		  else		   
		    tr->nodep[*count + tr->mxtips + 1] = np[i]->next->next;		    		    
		}

	      found = 1;	      	     
	      *count = *count + 1;
	    }
	} 
      
      assert(found != 0);
     
      reorderNodes(tr, np, p->next->back, count);     
      reorderNodes(tr, np, p->next->next->back, count);                
    }
}

void nodeRectifier(tree *tr)
{
  nodeptr *np = (nodeptr *)malloc(2 * tr->mxtips * sizeof(nodeptr));
  int i;
  int count = 0;
  
  tr->start       = tr->nodep[1];
  tr->rooted      = FALSE;

  /* TODO why is tr->rooted set to FALSE here ?*/
  
  for(i = tr->mxtips + 1; i <= (tr->mxtips + tr->mxtips - 1); i++)
    np[i] = tr->nodep[i];           
  
  reorderNodes(tr, np, tr->start->back, &count);  

 
  free(np);
}


void makeParsimonyTree(tree *tr, analdef *adef)
{   
  nodeptr  p, f;    
  int  i, nextsp, mintrav, maxtrav, randomMP, startMP;
  int *perm, *informative, *aliaswgt, *model, *dataVector;
  char *parsimonyBuffer;     

  /* stuff for informative sites */

  informative = (int *)malloc(sizeof(int) * tr->cdta->endsite);
  
  aliaswgt    = (int *)malloc(sizeof(int) * tr->cdta->endsite);
  memcpy(aliaswgt, tr->cdta->aliaswgt, sizeof(int) * tr->cdta->endsite);

  model       = (int *)malloc(sizeof(int) * tr->cdta->endsite);
  memcpy(model, tr->model, sizeof(int) * tr->cdta->endsite);

  dataVector    = (int *)malloc(sizeof(int) * tr->cdta->endsite);
  memcpy(dataVector, tr->dataVector, sizeof(int) * tr->cdta->endsite);
  
  parsimonyBuffer = (char *)malloc(tr->originalCrunchedLength * tr->rdta->numsp * sizeof(char));
  memcpy(parsimonyBuffer, tr->rdta->y0, tr->originalCrunchedLength * tr->rdta->numsp * sizeof(char)); 
  
  /* end */

  perm = (int *)malloc((tr->mxtips + 1) * sizeof(int));  
  
  determineUninformativeSites(tr, informative);    

  fixModelIndices(tr, adef, tr->parsimonyLength); 

  makePermutation(perm, tr->mxtips, adef);
  
  allocNodexParsimony(tr);
  
  tr->ntips = 0;    
  
  tr->nextnode = tr->mxtips + 1;       
  buildSimpleTree(tr, perm[1], perm[2], perm[3]);      

  while (tr->ntips < tr->mxtips) 
    {	
      tr->bestParsimony = INT_MAX;
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];
      
      /* printf("ADDING SPECIES %d\n", p->number); */
      
      buildNewTip(tr, p);
      
      f = findAnyTip(tr->start, tr->rdta->numsp);
      f = f->back;
      addTraverseParsimony(tr, p->back, f, 1, tr->ntips - 2);               	 	
      restoreTreeParsimony(tr, p->back, tr->insertNode);		      
      
      /* printf("MP %d\n", tr->bestParsimony); */
      
      assert(INT_MAX - tr->bestParsimony >= 1000);	  
    }    
 
  free(perm);    
  
  nodeRectifier(tr);   
  initravParsimonyNormal(tr, tr->start);
  initravParsimonyNormal(tr, tr->start->back);               
  
 
  mintrav = 1;
  maxtrav = 20;
  randomMP = tr->bestParsimony;        
  
  do
    {
      startMP = randomMP;
      
      for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
	{
	  rearrangeParsimony(tr, tr->nodep[i], mintrav, maxtrav);
	  if(tr->bestParsimony < randomMP)
	    {		
	      restoreTreeRearrangeParsimony(tr);
	      randomMP = tr->bestParsimony;
	    }
	}      		  	   
    }
  while(randomMP < startMP);
      
  /* printf("REARRANGEMENT MP Score %d Time %f\n", tr->bestParsimony, gettime() - masterTime); */
   
  nodeRectifier(tr);       
    
  /* repair */
  
  memcpy(tr->cdta->aliaswgt, aliaswgt, sizeof(int) * tr->cdta->endsite);

  memcpy(tr->model,    model   , sizeof(int) * tr->cdta->endsite);
 
  memcpy(tr->dataVector, dataVector, sizeof(int) * tr->cdta->endsite);
  
  memcpy(tr->rdta->y0, parsimonyBuffer, tr->originalCrunchedLength * tr->rdta->numsp * sizeof(char));

  fixModelIndices(tr, adef, tr->cdta->endsite);

  /* 
     for(i = 0; i < tr->NumberOfModels; i++)
     printf("%d: %d %d\n", i, tr->modelIndices[i][0],  tr->modelIndices[i][1]);
     printf("%d: %d %d\n", i, tr->partitionData[i].lower,  tr->partitionData[i].upper); 
  */

  /* repair end */  

  free(informative);    
  free(parsimonyBuffer);
  free(model);
  free(dataVector);
  free(aliaswgt);

  freeNodexParsimony(tr);
} 



void makeParsimonyTreeIncomplete(tree *tr, analdef *adef)
{   
  nodeptr  p, f;    
  int  i, j, k, nextsp, mintrav, maxtrav, randomMP, startMP;
  int *perm, *informative, *aliaswgt, *model, *dataVector;
  char *parsimonyBuffer;     

  /* stuff for informative sites */  

  informative = (int *)malloc(sizeof(int) * tr->cdta->endsite);
  
  aliaswgt    = (int *)malloc(sizeof(int) * tr->cdta->endsite);
  memcpy(aliaswgt, tr->cdta->aliaswgt, sizeof(int) * tr->cdta->endsite);

  model       = (int *)malloc(sizeof(int) * tr->cdta->endsite);
  memcpy(model, tr->model, sizeof(int) * tr->cdta->endsite);

  dataVector    = (int *)malloc(sizeof(int) * tr->cdta->endsite);
  memcpy(dataVector, tr->dataVector, sizeof(int) * tr->cdta->endsite);
  
  parsimonyBuffer = (char *)malloc(tr->originalCrunchedLength * tr->rdta->numsp * sizeof(char));
  memcpy(parsimonyBuffer, tr->rdta->y0, tr->originalCrunchedLength * tr->rdta->numsp * sizeof(char)); 
  
  /* end */

  perm = (int *)malloc((tr->mxtips + 1) * sizeof(int));

  if(!tr->grouped)
    {     
      for(i = 1; i <= tr->mxtips; i++)
	tr->constraintVector[i] = 0;
    }  
  
  determineUninformativeSites(tr, informative);

  fixModelIndices(tr, adef, tr->parsimonyLength);

  allocNodexParsimony(tr);
  
  if(!tr->grouped)
    {
      initravParsimony(tr, tr->start,       tr->constraintVector);
      initravParsimony(tr, tr->start->back, tr->constraintVector);
    }
  else
    {
      initravParsimonyNormal(tr, tr->start);
      initravParsimonyNormal(tr, tr->start->back);      
    }
    
  /* printf("Incomplete Parsimony score %d\n", evaluateParsimony(tr, tr->start)); */
  
  j = tr->ntips + 1;
  if(!tr->grouped)
    {
      for(i = 1; i <= tr->mxtips; i++)      
	if(tr->constraintVector[i] == 0) perm[j++] = i;	    	  
    }
  else
    {
      for(i = 1; i <= tr->mxtips; i++)      
	{
	  if(tr->constraintVector[i] == -1) 
	    {
	      perm[j++] = i;		
	      tr->constraintVector[i] = -9;
	    }
	}
    }
     
#ifdef PARALLEL
  srand((unsigned int) gettimeSrand()); 
#else
  if(adef->parsimonySeed == 0)   
    srand((unsigned int) gettimeSrand());         
#endif

  for (i = tr->ntips + 1; i <= tr->mxtips; i++) 
    {
#ifdef PARALLEL         
      k        = randomInt(tr->mxtips + 1 - i);
#else
      if(adef->parsimonySeed == 0) 
	k        = randomInt(tr->mxtips + 1 - i);
      else
	k =  (int)((double)(tr->mxtips + 1 - i) * randum(&adef->parsimonySeed));
#endif
      assert(i + k <= tr->mxtips);
      j        = perm[i];
      perm[i]     = perm[i + k];
      perm[i + k] = j; 
    }             
 
#ifdef DEBUG_CONSTRAINTS        
  for(i = 1; i <= tr->mxtips; i++)     
    printf("TIP %s %d\n", tr->nameList[i], tr->constraintVector[i]);              
#endif

  while (tr->ntips < tr->mxtips) 
    {	
      tr->bestParsimony = INT_MAX;
      nextsp = ++(tr->ntips);             
      p = tr->nodep[perm[nextsp]];
      
      /*printf("ADDING SPECIES %d %s\n", perm[nextsp], tr->nameList[perm[nextsp]]);*/
      
      buildNewTip(tr, p);      
      
      if(tr->grouped)
	{
	  int number = p->back->number;
	  tr->constraintVector[number] = -9;
	}
      
      f = findAnyTip(tr->start, tr->rdta->numsp);
      f = f->back;      
      
      if(tr->grouped)
	{
	  tr->grouped = FALSE;
	  addTraverseParsimony(tr, p->back, f, 1, tr->ntips - 2);  
	  tr->grouped = TRUE;
	}
      else
	addTraverseParsimony(tr, p->back, f, 1, tr->ntips - 2);
      
      restoreTreeParsimony(tr, p->back, tr->insertNode);		      
      
      /* printf("MP %d\n", tr->bestParsimony); */
      
      assert(INT_MAX - tr->bestParsimony >= 1000);	 
    }               
  
  free(perm);
  
  /* printf("RANDOM ADDITION MP Score %d Time %f\n", tr->bestParsimony, gettime() - masterTime); */   

  nodeRectifier(tr);   
  initravParsimonyNormal(tr, tr->start);
  initravParsimonyNormal(tr, tr->start->back); 
  
  if(adef->mode != PARSIMONY_ADDITION)
    {
      mintrav = 1;
      maxtrav = 20;
      randomMP = tr->bestParsimony;        
      
      do
	{
	  startMP = randomMP;
	  
	  for(i = 1; i <= tr->mxtips + tr->mxtips - 2; i++)
	    {		
	      if(rearrangeParsimony(tr, tr->nodep[i], mintrav, maxtrav))
		{
		  if(tr->bestParsimony < randomMP)
		    {		
		      restoreTreeRearrangeParsimony(tr);
		      randomMP = tr->bestParsimony;
		    }				
		}		    		     		 
	    }   
	}
      while(randomMP < startMP);
      
      /*printf("REARRANGEMENT MP Score %d Time %f\n", tr->bestParsimony, gettime() - masterTime);*/
    }
  else
    {               
      return;
    }
  
  nodeRectifier(tr);  
  
   /* repair */
  
  memcpy(tr->cdta->aliaswgt, aliaswgt, sizeof(int) * tr->cdta->endsite);

  memcpy(tr->model,    model   , sizeof(int) * tr->cdta->endsite);
 
  memcpy(tr->dataVector, dataVector, sizeof(int) * tr->cdta->endsite);
  
  memcpy(tr->rdta->y0, parsimonyBuffer, tr->originalCrunchedLength * tr->rdta->numsp * sizeof(char));

  fixModelIndices(tr, adef, tr->cdta->endsite);

  /* repair end */  

  free(informative);    
  free(parsimonyBuffer);
  free(model);
  free(dataVector);
  free(aliaswgt); 
  
  

  freeNodexParsimony(tr);              
} 


