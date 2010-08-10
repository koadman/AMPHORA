/*  RAxML-HPC, a program for sequential and parallel estimation of phylogenetic trees 
 *  Copyright March 2006 by Alexandros Stamatakis
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
 *  stamatak@ics.forth.gr
 *
 *  When publishing work that is based on the results from RAxML-VI-HPC please cite:
 *  
 *  Alexandros Stamatakis: "An Efficient Program for phylogenetic Inference Using Simulated Annealing". 
 *  Proceedings of IPDPS2005,  Denver, Colorado, April 2005.
 *  
 *  AND
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

#include <limits.h>
#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "axml.h"

extern FILE *INFILE;
extern char run_id[128];
extern char workdir[1024];
extern char bootStrapFile[1024];
extern char tree_file[1024];
extern char infoFileName[1024];
extern char resultFileName[1024];


int countTips(nodeptr p, int numsp)
{
  if(isTip(p->number, numsp))
    return 1;
  {
    nodeptr q;
    int tips = 0;

    q = p->next;
    while(q != p)
      { 
	tips += countTips(q->back, numsp);
	q = q->next;
      } 
    
    return tips;
  }
}

static void getTips(nodeptr p, int *c, int *entries, int numsp)
{
  nodeptr q;   

  if(isTip(p->number, numsp))
    {
      entries[*c] = p->number;
      *c = *c + 1;
      return;
    } 
    
  q = p->next;
  while(q != p)
    { 
      getTips(q->back, c, entries, numsp);
      q = q->next;
    }   

  return;
}

static int intCompare(const void *p1, const void *p2)
{
 int *rc1 = (int *)p1;
 int *rc2 = (int *)p2;

 int i = *rc1;
 int j = *rc2;
  
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}



static void makeBipartitionsRec(nodeptr p, bList *blThis, int *bCountThis, int numsp)
{  
  if(/*p->tip*/ isTip(p->number, numsp))
    return;
  {
    nodeptr q;
    int l, r;
    int c;   
    
    if(/*!p->back->tip*/ ! isTip(p->back->number, numsp))
      {       
	l = countTips(p, numsp);
	r = countTips(p->back, numsp);
	c = 0;
	     
	if(l < r)
	  {
	    blThis[*bCountThis].entries = (int *)malloc(l * sizeof(int));	   
	    getTips(p, &c, blThis[*bCountThis].entries, numsp);
	  }
	else
	  {
	    blThis[*bCountThis].entries = (int *)malloc(r * sizeof(int));
	    getTips(p->back, &c, blThis[*bCountThis].entries, numsp);
	  }
	
	blThis[*bCountThis].length = c;      

	qsort((blThis[*bCountThis].entries), c, sizeof(int), intCompare);
	blThis[*bCountThis].p = p;
	blThis[*bCountThis].pNum = p->number;
	blThis[*bCountThis].qNum = p->back->number;
	*bCountThis = *bCountThis + 1;
      }
  
    q = p->next;
    while(q != p)
      {
	makeBipartitionsRec(q->back, blThis, bCountThis, numsp);
	q = q->next;
      } 
    return;
  }
}


static int bListCompare(const void *p1, const void *p2)
{
 bList *rc1 = (bList *)p1;
 bList *rc2 = (bList *)p2;

 int i = rc1->length;
 int j = rc2->length;
  
  if (i > j)
    return (1);
  if (i < j)
    return (-1);
  return (0);
}



static bList *bipartitionList(tree *tr,  int *bCountThis)
{
  int i, n = tr->mxtips - 3;
  bList *blThis;

  blThis = (bList *)malloc(sizeof(bList) * n);       
  *bCountThis = 0;

  for(i = 0; i < n; i++)
    {        
      blThis[i].p = (nodeptr) NULL;
      blThis[i].support = 0;
      blThis[i].length = 0;
      blThis[i].entries = (int *)NULL;
      blThis[i].pNum = 0;
      blThis[i].qNum = 0;
    }
      
  makeBipartitionsRec(tr->nodep[1]->back, blThis, bCountThis, tr->rdta->numsp);
            
  qsort(&(blThis[0]), *bCountThis, sizeof(bList), bListCompare);
        
  return blThis;            
}


/*void printBlist(bList *blThis, int n)
{
  int i, j;

  for(i = 0; i < n; i++)	    
    {
      printf("%d %d %d: (", i, blThis[i].length, blThis[i].support);
      for(j = 0; j < blThis[i].length; j++)
	{
	  if(j == (blThis[i].length - 1))
	    printf("%d)\n", blThis[i].entries[j]);
	  else
	     printf("%d, ", blThis[i].entries[j]);
	}
    }    
    }*/

static void freeBList(bList *blThis, int n)
{
  int i;
 
  for(i = 0; i < n; i++)	             
    {    
      free(blThis[i].entries);	      
    }  
}

static void updateReferenceList(bList *referenceList, int referenceListLength, bList *currentList, int currentListLength)
{
  int i, j, k, length;
  boolean found;
  int f = 0;

  /*  printf("%d %d\n", currentListLength, referenceListLength);*/

  for(i = 0; i < currentListLength; i++)
    {
      j = 0;
      length = currentList[i].length;
      while(length > referenceList[j].length)
	j++;

      /*printf("%d at %d\n", length, j);*/

      while(j < referenceListLength && length == referenceList[j].length)
	{	
	  k = 0;
	  found = TRUE;
	  while(k < length && found)
	    {
	      if(currentList[i].entries[k] != referenceList[j].entries[k])
		found = FALSE;
	      k++;
	    }

	  if(found)
	    {
	      referenceList[j].support = referenceList[j].support + 1;	      	      
	      f++;
	      /*goto foundThisOne;*/
	      break;
	    }
	  j++;
	}
      /*foundThisOne:	         */
    }
  /*printf("FOUND %d\n", f);*/
}


void calcBipartitions(tree *tr, analdef *adef, char *bestTreeFileName, char *bootStrapFileName)
{
  bList *ML_Tree = (bList *)NULL, *BOOT_Tree = (bList *)NULL;
  int countML_Tree = 0, countBOOT_Tree = 0, numberOfTrees = 0, i;
  char ch;
  FILE *infoFile;  
 
  INFILE = fopen(bestTreeFileName, "r");
  treeReadTopologyOnly(INFILE, tr, adef, FALSE);
  fclose(INFILE); 
  
  ML_Tree = bipartitionList(tr, &countML_Tree);
  
  INFILE = fopen(bootStrapFileName, "r");
       

  while((ch = getc(INFILE)) != EOF)
    {
      if(ch == ';')
	numberOfTrees++;
    }	 
  rewind(INFILE);
 
/*  if(!adef->allInOne)
    {
      infoFile = fopen(infoFileName, "a");
      printf("\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
      fprintf(infoFile, "\n\nAA Found %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
      fclose(infoFile);
    }
*/
  for(i = 0; i < numberOfTrees; i++)
    {
      /*printf("Tree %d\n", i);*/
      treeReadTopologyOnly(INFILE, tr, adef, FALSE);
      BOOT_Tree = bipartitionList(tr, &countBOOT_Tree);
      updateReferenceList(ML_Tree, countML_Tree, BOOT_Tree, countBOOT_Tree);  
      freeBList(BOOT_Tree, countBOOT_Tree);
      free(BOOT_Tree);
    }

  fclose(INFILE);

  /*printBlist(ML_Tree, countML_Tree);*/
 
  INFILE = fopen(bestTreeFileName, "r");
  treeReadTopologyOnly(INFILE, tr, adef, TRUE);
  fclose(INFILE);

  /*for(i = 0; i < countML_Tree; i++)
    {
      p = ML_Tree[i].p;
      p->z = p->back->z =  ((double)ML_Tree[i].support) / ((double) numberOfTrees);
      }            */

  tr->ML_Tree = ML_Tree;
  tr->countML_Tree = countML_Tree;
  tr->numberOfTrees = numberOfTrees; 

  printBipartitionResult(tr, adef, TRUE);  
  freeBList(ML_Tree, countML_Tree);  
  free(ML_Tree);
}

/*
void drawBipartitionsOnTree(tree *tr, analdef *adef, topolRELL_LIST *rl, int numberOfTrees)
{
  bList *ML_Tree = (bList *)NULL, *BOOT_Tree = (bList *)NULL;
  nodeptr p, q;
  int countML_Tree = 0, countBOOT_Tree = 0, i;
  char ch;
  FILE *infoFile;     
  topolRELL_LIST *referenceTree;     

  referenceTree = (topolRELL_LIST *)malloc(sizeof(topolRELL_LIST));
  initTL(referenceTree, tr, 1);  
  saveTL(referenceTree, tr, 0);  
  restoreTL(referenceTree, tr, 0);  

  ML_Tree = bipartitionList(tr, &countML_Tree); 

  for(i = 0; i < numberOfTrees; i++)
    {      
      restoreTL(rl, tr, i);
      BOOT_Tree = bipartitionList(tr, &countBOOT_Tree);
      updateReferenceList(ML_Tree, countML_Tree, BOOT_Tree, countBOOT_Tree);  
      freeBList(BOOT_Tree, countBOOT_Tree);      
      free(BOOT_Tree);
    } 

  tr->ML_Tree = ML_Tree;
  tr->countML_Tree = countML_Tree;
  tr->numberOfTrees = numberOfTrees;

  restoreTL(referenceTree, tr, 0);
  for(i = tr->mxtips + 1; i <= tr->mxtips + tr->mxtips - 2; i++)
    {
      int res = -1, res2 = -1;
      int countTips = 0;
      nodeptr tips[2];

      p = tr->nodep[i];
      if(!isTip(p->back->number, numsp))
	{
	  res = lookupBipartition(tr, p);	 
	  p->support = p->back->support = res;
	}
      else		
	tips[countTips++] = p;
	  	
      q = p->next;
      while(q != p)
	{
	  if(! isTip(q->back->number, numsp))
	    {
	      res2 = lookupBipartition(tr, q);
	      q->support = q->back->support = res2;	    
	    }
	  else
	    tips[countTips++] = q;
	      	  
	  q = q->next;	  
	}

      assert(countTips <= 2);
      
      if(countTips == 1)
	tips[0]->support = (res + res2) / 2;
      else
	{
	  if(res == -1)
	    res = res2;
	  assert(res != -1);
	  tips[0] = res;
	  tips[1] = res;
	}       
    }


freeTL(referenceTree, tr);  
  free(referenceTree); 

  freeBList(ML_Tree, countML_Tree);  
  free(ML_Tree); 
  }*/


/*************************************************************************************************************/


static const unsigned char bitmask[8] = {1, 2, 4, 8, 16, 32, 64, 128};

#ifdef WIN32
static void set_bit(unsigned char *vector, int pos)
#else
static inline void set_bit(unsigned char *vector, int pos)
#endif
{
   vector[pos / BITS_BYTE] |= bitmask[pos % BITS_BYTE];
}

#ifdef WIN32
static int get_bit(unsigned char *vector, int pos)
#else
static inline int get_bit(unsigned char *vector, int pos)
#endif
{
  if(vector[pos / BITS_BYTE] == 0)
    return 0;
  return ((vector[pos / BITS_BYTE] & bitmask[pos % BITS_BYTE]) == bitmask[pos % BITS_BYTE]);
} 

static void permute(int *perm, int n, long *seed)
{
  int  i, j, k;
 
  for (i = 0; i < n; i++) 
    {
      k =  (int)((double)(n - i) * randum(seed));
      j        = perm[i];    
      perm[i]     = perm[i + k];
      perm[i + k] = j; 
      /*assert(i + k < n);*/
    }
}

static double testFreq(double *vect1, double *vect2, int n)
{
  int i;
  double
    avg1 = 0.0, 
    avg2 = 0.0,
    sum_xy = 0.0, 
    sum_x  = 0.0, 
    sum_y  = 0.0,
    corr   = 0.0;

  for(i = 0; i < n; i++)
    {	     
      avg1 += vect1[i];
      avg2 += vect2[i];
    }
      
  avg1 /= ((double)n);
  avg2 /= ((double)n);

  printf("Average %f %f\n", avg1, avg2);
    
  for(i = 0; i < n; i++)
    {
      sum_xy += ((vect1[i] - avg1) * (vect2[i] - avg2));
      sum_x  += ((vect1[i] - avg1) * (vect1[i] - avg1));
      sum_y  += ((vect2[i] - avg2) * (vect2[i] - avg2));	 
    }

  corr = sum_xy / (sqrt(sum_x) * sqrt(sum_y));
   
#ifndef WIN32
  if(isnan(corr))
    {
      printf("Numerical Error pearson correlation is not a number\n");
      assert(0);
    }
#endif

  return corr;
}

static double testFreq2(double *vect1, double *vect2, int n, int reps, double bootstopCutoff)
{
  int i;
  double  
    sum_xy = 0.0, 
    sum_x  = 0.0, 
    sum_y  = 0.0,
    corr   = 0.0;  
  double 
    av1 = 0.0,
    av2 = 0.0;
 
 
  if(bootstopCutoff > 0.0)
    {
      const double scale  = 1.0 / ((double)reps * 0.5);
      const double upper_bound = 0.5 + bootstopCutoff;
      const double lower_bound = 0.5 - bootstopCutoff;
      int 
	count = 0,
	verify = 0;
      char *isSet = (char *)calloc(n, sizeof(char));
      
      assert(bootstopCutoff > 0.0 && bootstopCutoff <= 0.5);
      assert(lower_bound < upper_bound);
      assert(lower_bound >= 0.0 && upper_bound <= 1.0);     

      for(i = 0; i < n; i++)
	{
	  vect1[i] *= scale;
	  vect2[i] *= scale;
	}
      
      for(i = 0; i < n; i++)
	{
	  if((vect1[i] >= lower_bound && vect1[i] <= upper_bound) || (vect2[i] >= lower_bound && vect2[i] <= upper_bound))
	    {
	      av1 += vect1[i];
	      av2 += vect2[i];
	      isSet[i] = 1;
	      count++;
	    }
	} 
      
      if(count == 0)
	{	  
	  return 0.0;
	}
      
      av1 /= ((double)count);
      av2 /= ((double)count);
      
      for(i = 0; i < n; i++)
	{
	  if(isSet[i])
	    {
	      sum_xy += ((vect1[i] - av1) * (vect2[i] - av2));
	      sum_x  += ((vect1[i] - av1) * (vect1[i] - av1));
	      sum_y  += ((vect2[i] - av2) * (vect2[i] - av2));	 
	      /*printf("%f %f ", vect1[i], vect2[i]);*/
	      verify++;
	    }
	}
      /*printf("\n");*/
      free(isSet);

      assert(verify == count);

     

      corr = sum_xy / (sqrt(sum_x) * sqrt(sum_y));            
#ifndef WIN32
      if(isnan(corr))
	{
	  printf("Numerical Error pearson correlation is not a number\n");
	  assert(0);
	}
#endif
       
      return corr;
    }
  else
    {     
      for(i = 0; i < n; i++)
	{	  
	  av1 += vect1[i];
	  av2 += vect2[i];       
	} 
      
      av1 /= ((double)n);
      av2 /= ((double)n);
      
      for(i = 0; i < n; i++)
	{	  
	  sum_xy += ((vect1[i] - av1) * (vect2[i] - av2));
	  sum_x  += ((vect1[i] - av1) * (vect1[i] - av1));
	  sum_y  += ((vect2[i] - av2) * (vect2[i] - av2));	 	    
	}

     

      corr = sum_xy / (sqrt(sum_x) * sqrt(sum_y));
#ifndef WIN32
      if(isnan(corr))
	{
	  printf("Numerical Error pearson correlation is not a number\n");
	  assert(0);
	}
#endif
      return corr;
    }
}

/*static double testFreq3(double *vect1, double *vect2, double avg1, double avg2, int n)
{
  int i;
  double  
    sum_xy = 0.0, 
    sum_x  = 0.0, 
    sum_y  = 0.0,
    corr   = 0.0;  
  printf("Average %f %f\n", avg1, avg2);
  for(i = 0; i < n; i++)
    {
      sum_xy += ((vect1[i] - avg1) * (vect2[i] - avg2));
      sum_x  += ((vect1[i] - avg1) * (vect1[i] - avg1));
      sum_y  += ((vect2[i] - avg2) * (vect2[i] - avg2));	 
    }

  corr = sum_xy / (sqrt(sum_x) * sqrt(sum_y));
   
  return corr;
}
*/

static double probks(double alam)
{
  int i;
  double 
    EPS1=0.001, 
    EPS2=1.E-8, 
    FAC = 2.0,
    A2 = -2.0 *  (alam * alam),
    PROBKS=0.0,
    TERMBF=0.0,
    TERM = 0.0;
  
  for(i = 1; i <= 100; i++)
    {
      TERM = FAC * exp(A2 * i * i);
      PROBKS = PROBKS + TERM;
      if((fabs(TERM) < EPS1 * TERMBF) || (fabs(TERM) < EPS2 * PROBKS))
	return PROBKS;      
      FAC =-FAC;
      TERMBF= fabs(TERM);
    }

  PROBKS = 1.0;
  
  return PROBKS;
}



static double ksTwo(int *vect1, int *vect2, int n)
{
  int     
    j1 = 1, 
    j2 = 1;
  
  double 
    en1 = (double)n,
    en2 = (double)n,
    f01 = 0.0, 
    f02 = 0.0, 
    fn1 = 0.0,
    fn2 = 0.0,
    d1  = 0.0,
    d2  = 0.0,
    d   = 0.0,
    dt  = 0.0,
    result;

  qsort(vect1, n, sizeof(int), intCompare);
  qsort(vect2, n, sizeof(int), intCompare);
  
  while(j1 <= n && j2 <= n)
    {
      if(vect1[j1 - 1] < vect2[j2 - 1])
	{
	  fn1 = j1 / en1;
	  d1 = fabs(fn1 - f02);
	  d2 = fabs(f01 - f02);
	  if(d1 > d2)
	    dt = d1;
	  else
	    dt = d2;
	  if(dt > d)
	    d = dt;
	  f01 = fn1;
	  j1++;
	}
      else
	{
	  fn2 = j2/en2;
	  d1 = fabs(fn2 - f01);
	  d2 = fabs(f02 - f01);
	  if(d1 > d2) 
	    dt = d1;
	  else
	    dt = d2;
	  if(dt > d)
	    d = dt;
	  f02 = fn2;
	  j2++;
	}
    } 
	  
  result = probks(sqrt(en1 * en2 / (en1 + en2)) * d);
 
  return result;
}

static void computeAllLHs(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int 
    numberOfTrees = 0,   
    i;
  char ch; 
  double 
    bestLH = unlikely;    
  bestlist *bestT;
  FILE *infoFile, *result;
  

  infoFile = fopen(infoFileName, "a");
  result   = fopen(resultFileName, "w");

  bestT = (bestlist *) malloc(sizeof(bestlist));
  bestT->ninit = 0;
  initBestTree(bestT, 1, tr->mxtips);

  allocNodex(tr, adef); 

  INFILE = fopen(bootStrapFileName, "r");       
  while((ch = getc(INFILE)) != EOF)
    {
      if(ch == ';')
	numberOfTrees++;
    }	 
  rewind(INFILE);
 
  printf("\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
  fprintf(infoFile, "\n\nBB Found %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
 
  for(i = 0; i < numberOfTrees; i++)
    {              
      treeReadLen(INFILE, tr, adef);      
      
      if(i == 0)
	{
	  modOpt(tr, adef);
	  printf("Model optimization, first Tree: %f\n", tr->likelihood);
	  fprintf(infoFile, "Model optimization, first Tree: %f\n", tr->likelihood);
	  bestLH = tr->likelihood;
	  resetBranches(tr);
	}
      
      treeEvaluate(tr, 2);
      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, 
		  TRUE, adef, SUMMARIZE_LH);
                 
      fprintf(result, "%s", tr->tree_string);
      
      saveBestTree(bestT, tr);

      if(tr->likelihood > bestLH)		
	bestLH   = tr->likelihood;	
      printf("Tree %d Likelihood %f\n", i, tr->likelihood);
      fprintf(infoFile, "Tree %d Likelihood %f\n", i, tr->likelihood);
    }        
    
  recallBestTree(bestT, 1, tr);  
  evaluateGeneric(tr, tr->start);
  printf("Model optimization, %f <-> %f\n", bestLH, tr->likelihood); 
  fprintf(infoFile, "Model optimization, %f <-> %f\n", bestLH, tr->likelihood); 
  modOpt(tr, adef);
  treeEvaluate(tr, 2);
  printf("Model optimization, %f <-> %f\n", bestLH, tr->likelihood);
  fprintf(infoFile, "Model optimization, %f <-> %f\n", bestLH, tr->likelihood); 

  printf("\nAll evaluated trees with branch lengths written to File: %s\n", resultFileName);
  fprintf(infoFile, "\nAll evaluated trees with branch lengths written to File: %s\n", resultFileName);

  fclose(INFILE); 
  fclose(infoFile);
  fclose(result);
  exit(0);
}



static void addBipartitionsFaster(BL *b, int number, int *candidate, int candidateLength, int modulo)
{
  int i, *entries, k;   
  const int bCount = b->count;
  bipList *bibPtr = b->b;
 
  for(i = 0; i < bCount; i++)
    {                     
      if(candidateLength == bibPtr[i].length)
	{	
	  entries = bibPtr[i].entries;
	  
	  for(k = 0; k < candidateLength; k++)	    
	    if(candidate[k] != entries[k])
	      break;	    
	  
	  if(k == candidateLength)
	    {
	      if(modulo == 0)		
		set_bit(bibPtr[i].isSet, number);	     	     		
	      else
		{
		  if(modulo == 1)
		    bibPtr[i].freq1 =  bibPtr[i].freq1 + 1;
		  else
		    bibPtr[i].freq2 =  bibPtr[i].freq2 + 1;
		}
	      return;	     
	    }	        
	}
    }
           
  if(b->count == b->n)
    {       
      bipList *bip = (bipList *)malloc(sizeof(bipList) *  b->n * 2);      
      
      memcpy(bip, b->b, b->n * sizeof(bipList));
      
      for(k = b->n; k  < (b->n * 2); k++)
	{		
	  bip[k].length  = 0;	     
	  bip[k].entries = (int *)NULL;
	  if(modulo == 0)		
	    bip[k].isSet   = (unsigned char*)calloc(b->treeVectorLength, sizeof(unsigned char));		
	  else
	    {
	      bip[k].freq1 = 0;
	      bip[k].freq2 = 0;
	    }
	}
      free(b->b);
      b->b = bip;
      b->n *= 2;	      
    }   
  
  b->b[b->count].entries = (int *)malloc(sizeof(int) * candidateLength);
  b->b[b->count].length  = candidateLength;
  
  if(modulo == 0)   
    set_bit(b->b[b->count].isSet, number);    
  else
    {
      if(modulo == 1)
	b->b[b->count].freq1 =  b->b[b->count].freq1 + 1;
      else
	b->b[b->count].freq2 =  b->b[b->count].freq2 + 1;
    }
  
  memcpy(b->b[b->count].entries, candidate, candidateLength * sizeof(int));
  
  b->count =  b->count + 1;
}



double addTime = 0.0;

static void addRecBL(nodeptr p, BL *b, int number, int modulo, int numberOfTips, int numsp)
{  
  if(/*p->tip*/ isTip(p->number, numsp))
    return;
  {
    nodeptr q;    
    
    if(/*!p->back->tip*/ ! isTip(p->back->number, numsp))
      {       
	int *entries;
	int length = 0; 
	int l, r;  

	l = countTips(p, numsp);
	if(l < ((numberOfTips/2) + 1))
	  {
	    entries = (int *)malloc(l * sizeof(int));	   
	    getTips(p, &length, entries, numsp);
	  }
	else
	  {
	    r = numberOfTips - l;
	    entries = (int *)malloc(r * sizeof(int));
	    getTips(p->back, &length, entries, numsp);	  
	  }	    

	qsort(entries, length, sizeof(int), intCompare);

	{
	  double t = gettime();
	  addBipartitionsFaster(b, number, entries, length, modulo);
	  addTime += (gettime() - t);
	}

	free(entries);	
      }
  
    q = p->next;
    while(q != p)
      {
        addRecBL(q->back, b, number, modulo, numberOfTips, numsp);
	q = q->next;
      } 
    return;
  }
}


static void compareBips(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int 
    numberOfTreesAll = 0, 
    numberOfTreesStop = 0,
    i;
  char ch; 
  double *vect1, *vect2, p;
  int 
    bipAll = 0,
    bipStop = 0;
  char bipFileName[1024];
  FILE *outf; 
  BL *b = (BL *)malloc(sizeof(BL));

  b->n = 100 * tr->mxtips;
  b->count = 0;
 
  b->treeVectorLength = 0;      

  b->b = (bipList *)malloc(sizeof(bipList) *  b->n);

  for(i = 0; i < b->n; i++)
    {     
      b->b[i].length = 0;
      b->b[i].freq1  = 0;
      b->b[i].freq2  = 0;
      b->b[i].entries = (int *)NULL;
    }            
               
  INFILE = fopen(bootStrapFileName, "r");       
  while((ch = getc(INFILE)) != EOF)
    {
      if(ch == ';')
	numberOfTreesAll++;
    }	 
  rewind(INFILE);
      
  printf("\n\nFound %d trees in File %s\n\n", numberOfTreesAll, bootStrapFileName);
              
  for(i = 0; i < numberOfTreesAll; i++)
    {              
      treeReadTopologyOnly(INFILE, tr, adef, FALSE);
      addRecBL(tr->nodep[1]->back, b, i, 1, tr->mxtips, tr->rdta->numsp); 
    }
	  
  fclose(INFILE); 

  /* do BOOTSTOP ********************************************************************************************************/  

  INFILE = fopen(tree_file, "r");       
  while((ch = getc(INFILE)) != EOF)
    {
      if(ch == ';')
	numberOfTreesStop++;
    }	 
  rewind(INFILE);
      
  printf("\n\nFound %d trees in File %s\n\n", numberOfTreesStop, tree_file);
       
  for(i = 0; i < numberOfTreesStop; i++)
    {              
      treeReadTopologyOnly(INFILE, tr, adef, FALSE);
      addRecBL(tr->nodep[1]->back, b, i, 2, tr->mxtips, tr->rdta->numsp); 
    }
	  
  fclose(INFILE); 

  /***************************************************************************************************/
   

  vect1 = (double *)malloc(sizeof(double) * b->count);
  vect2 = (double *)malloc(sizeof(double) * b->count);

  strcpy(bipFileName,         workdir);  
  strcat(bipFileName,         "RAxML_bipartitionFrequencies.");
  strcat(bipFileName,         run_id);

  outf = fopen(bipFileName, "w");

  for(i = 0; i < b->count; i++)	     
    {
      vect1[i] = ((double)b->b[i].freq1) / ((double)numberOfTreesAll);
      if(b->b[i].freq1 > 0)
	bipAll++;
      vect2[i] = ((double)b->b[i].freq2)/ ((double)numberOfTreesStop);
      if(b->b[i].freq2 > 0)
	bipStop++;
      fprintf(outf, "%f %f\n", vect1[i], vect2[i]);
    }
  
  fclose(outf);

  p = testFreq(vect1, vect2, b->count);

  printf("Pearson: %f Bipartitions-All: %d Bipartitions-Stop: %d\n", p, bipAll, bipStop);

  /*p = pairedT(vect1, vect2, b->count);
     
  printf("pairedT: %f Bipartitions-All: %d Bipartitions-Stop: %d\n", p, bipAll, bipStop);*/

  exit(0);
}








void computeBootStopOnly(tree *tr, analdef *adef, char *bootStrapFileName)
{
  if(adef->bootStopOnly >= 3)
    {
      if(adef->bootStopOnly == 3)
	compareBips(tr, adef, bootStrapFileName);
      if(adef->bootStopOnly == 4)
	computeAllLHs(tr, adef, bootStrapFileName);	
    }
  else
    {
      int numberOfTrees = 0, i;
      boolean stop = FALSE;
      char ch;     
      BL *b = (BL *)malloc(sizeof(BL));
      double t, addT = 0.0, calcT = 0.0, readT = 0.0;
      int checkEvery;
      int treesAdded = 0;
     

      INFILE = fopen(bootStrapFileName, "r");       
      while((ch = getc(INFILE)) != EOF)
	{
	  if(ch == ';')
	    numberOfTrees++;
	}	 
      rewind(INFILE);
      
      printf("\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
      
      assert(sizeof(unsigned char) == 1);

      switch(adef->bootStopOnly)
	{
	case 1:
	  checkEvery = FC_SPACING;
	  b->n = FC_INIT * tr->mxtips;
	  break;
	case 2:         
	  checkEvery = BC_SPACING;
	  b->n = BC_INIT * tr->mxtips;
	  break;
	default:
	  assert(0);
	}
     
      b->count = 0;    
      b->treeVectorLength = (numberOfTrees / BITS_BYTE) + 1;         
      b->b = (bipList *)malloc(sizeof(bipList) *  b->n);

      for(i = 0; i < b->n; i++)
	{     
	  b->b[i].length = 0;
	  b->b[i].isSet = (unsigned char*)calloc(b->treeVectorLength, sizeof(unsigned char));
	  b->b[i].entries = (int *)NULL;
	}                       
      
      printf("# Trees \t Average Pearson Coefficient \t # Permutations: pearson >= %f\n", 
	     FC_LOWER);

      for(i = 1; i <= numberOfTrees && !stop; i++)
	{              
	  t = gettime();
	  treeReadTopologyOnly(INFILE, tr, adef, FALSE);	  	
	  readT += (gettime() - t);  
	  t = gettime();
	  addRecBL(tr->nodep[1]->back, b, (i - 1), 0, tr->mxtips, tr->rdta->numsp); 
	  addT += (gettime() - t);      
	  treesAdded++;	
	  
	  if(i > START_BSTOP_TEST && i % checkEvery == 0)
	    {		     
	      t = gettime();
	      if(adef->bootStopOnly == 1)
		{
		  int k, j, l;
		  int *perm =  (int *)malloc(sizeof(int) * i);
		  long seed = 12345;
		  double result;
		  double avg = 0;
		  double *vect1, *vect2;
		  int countBetter = 0;				  

		  for(j = 0; j < i; j++)
		    perm[j] = j;
		  
		  for(k = 0; k < BOOTSTOP_PERMUTATIONS; k++)
		    {   
		      unsigned char *set;		      
		     
		      permute(perm, i, &seed);
		      
		      vect1 = (double *)calloc(b->count, sizeof(double));
		      vect2 = (double *)calloc(b->count, sizeof(double));
		     		      
		      for(j = 0; j < b->count; j++)
			{		       
			  set = b->b[j].isSet;
			  
			  for(l = 0; l < i; l++)
			    {			     
			      if(get_bit(set,l))
				{
				  if(perm[l] % 2 == 0)
				    vect1[j] = vect1[j] + 1.0;				   
				  else
				    vect2[j] = vect2[j] + 1.0;				  
				}			     
			    }		       
			}		    
				  		      
		      result = testFreq2(vect1, vect2, b->count, i, adef->bootstopCutoff);
		     
		      if(result >= FC_LOWER)
			countBetter++;
		      
		      avg += result;
		      
		      free(vect1);		  
		      free(vect2);		 
		    }
		 
	      
		  avg /= BOOTSTOP_PERMUTATIONS;
		 
		  printf("%d \t\t\t %f \t\t\t\t %d\n", i, avg, countBetter);
		  free(perm);
	       
		  stop = (countBetter >= FC_THRESHOLD && avg >= FC_LOWER);
		}
	      if(adef->bootStopOnly == 2)
		{
		  int k, j, l;
		  int *perm =  (int *)malloc(sizeof(int) * i);
		  long seed = 12345;
		  double result;
		  double avg = 0;
		  int *vect1, *vect2;		 	      		
		  
		  for(j = 0; j < i; j++)
		    perm[j] = j;		 
		  

		  for(k = 0; k < BOOTSTOP_PERMUTATIONS; k++)
		    {     
		      unsigned char *set;		     
		      vect1 = (int *)calloc(i/2, sizeof(int));
		      vect2 = (int *)calloc(i/2, sizeof(int));	
	    		     
		      permute(perm, i, &seed);		     
		      		      		      		     
		      for(l = 0; l < b->count; l++)
			{			
			  boolean v1 = TRUE, v2 = TRUE;
			  set =  b->b[l].isSet;

			  for(j = 0; (j < i) && (v1 || v2); j++)
			    {			     
			      if(get_bit(set, j))
				{
				  if(perm[j] % 2 == 0)
				    {			
				      if(v1)
					{
					  vect1[perm[j]/2] += 1;
					  v1 = FALSE;
					}
				    }
				  else
				    {
				      if(v2)
					{
					  vect2[perm[j]/2] += 1;
					  v2 = FALSE;
					}
				    }
				}
			    }
			}		     		      	      		      
		     
		      result = ksTwo(vect1, vect2, i/2);
		      		     
		      avg += result;
		      free(vect1);		  
		      free(vect2);
		    }
		  		  		  
		  free(perm);
		  avg /= BOOTSTOP_PERMUTATIONS;
		 
		  printf("%d %f\n", i, avg);
		  stop = (avg <= BC_THRESHOLD);
		}

	      calcT += (gettime() - t);	     	     
	    }	 	   
	}
      printf("\n\nExecution time Analysis for development purposes:\n");
      printf("Tree-Reading: %f secs \nIdentify Bipartitions: %f secs \nAdd Bipartitions: %f secs \nCompute Statistics: %f secs\n\n", readT, addT - addTime, addTime, calcT);
      if(stop)
	{
	  char bootstopFileName[1028]; 
	  FILE *outf;
	  int bootStopNumber = 0;
	  
	  printf("Stopping after %d trees\n", treesAdded);
	  
	  strcpy(bootstopFileName,         workdir);  
	  strcat(bootstopFileName,         "RAxML_bootStop.");
	  strcat(bootstopFileName,         run_id);
	  
	  outf = fopen(bootstopFileName, "w");
	  
	  rewind(INFILE);
	  while(((ch = getc(INFILE)) != EOF) && (bootStopNumber < treesAdded))
	    {
	      fprintf(outf, "%c", ch);
	      if(ch == ';')
		bootStopNumber++;		   
	    }
	  fprintf(outf, "\n");
	  fclose(outf);
	  fclose(INFILE);
	  exit(0);
	}	  
      else
	{
	  printf("Bootstopping test did not converge after %d trees\n", treesAdded);
	  fclose(INFILE); 
	  exit(0);
	}
    }
}

boolean bootStop(tree *tr, BL *b, int numberOfTrees, double *pearsonAverage, double bootstopCutoff)
{
  int n = numberOfTrees + 1;
  addRecBL(tr->nodep[1]->back, b, numberOfTrees, 0, tr->mxtips, tr->rdta->numsp);  

  if((n > START_BSTOP_TEST) && (n % FC_SPACING == 0))
    {     
      int k, j, l;
      int *perm =  (int *)malloc(sizeof(int) * n);
      long seed = 12345;
      double result;
      double avg = 0;
      double *vect1, *vect2;
      int countBetter = 0;	         
      
      for(j = 0; j < n; j++)
	perm[j] = j;
      
      for(k = 0; k < BOOTSTOP_PERMUTATIONS; k++)
	{    
	  unsigned char *set;	 
	    
	  permute(perm, n, &seed);
	  
	  vect1 = (double *)calloc(b->count, sizeof(double));
	  vect2 = (double *)calloc(b->count, sizeof(double));
	 
	  for(j = 0; j < b->count; j++)
	    {
	      set = b->b[j].isSet;
	      
	      for(l = 0; l < n; l++)
		{	       		  
		  if(get_bit(set,l))
		    {
		      if(perm[l] % 2 == 0)
			vect1[j] = vect1[j] + 1;
		      else
			vect2[j] = vect2[j] + 1;
		    }
		}	      
	    }
	  
	  result = testFreq2(vect1, vect2, b->count, n, bootstopCutoff);
	 	 
	  if(result >= FC_LOWER)
	    countBetter++;
	  avg += result;
	  
	  free(vect1);		  
	  free(vect2);		 
	}
      
      avg /= BOOTSTOP_PERMUTATIONS;
           
      free(perm);

      *pearsonAverage = avg;      

      if(countBetter >= FC_THRESHOLD)
	return TRUE;
      else
	return FALSE;

    }
  else
    return FALSE;
}





