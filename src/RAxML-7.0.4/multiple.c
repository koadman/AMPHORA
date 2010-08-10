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
#include <unistd.h>
#endif

#include <math.h>
#include <time.h> 
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include "axml.h"


extern int  optimizeRatesInvocations;
extern int  optimizeRateCategoryInvocations;
extern int  optimizeAlphaInvocations;
extern int  optimizeInvarInvocations;
extern int  checkPointCounter;
extern int  Thorough;
extern int  partCount;
extern char tree_file[1024];

extern double masterTime;

extern FILE   *INFILE, *permutationFile, *logFile, *infoFile;

extern char seq_file[1024];
extern char permFileName[1024], resultFileName[1024], 
  logFileName[1024], checkpointFileName[1024], infoFileName[1024], run_id[128], workdir[1024], bootStrapFile[1024], bootstrapFileName[1024], bipartitionsFileName[1024]; 




void catToGamma(tree *tr, analdef *adef)
{
  freeNodex(tr); 
  
  assert(tr->rateHetModel == CAT);
  assert(adef->model == M_GTRCAT || adef->model == M_PROTCAT);
  
  if(adef->useInvariant)
    tr->rateHetModel = GAMMA_I;
  else
    tr->rateHetModel = GAMMA;

  switch(adef->model)
    {
    case M_GTRCAT:
      adef->model = M_GTRGAMMA;	     
      break;
    case M_PROTCAT:
      adef->model = M_PROTGAMMA;
      break;
    default:
      assert(0);
    }
#ifdef _LOCAL_DATA
  tr->currentModel = adef->model;
#endif

  assignLikelihoodFunctions(tr, adef);
  allocNodex(tr, adef);  
}

void gammaToCat(tree *tr, analdef *adef)
{
  freeNodex(tr);   	  

  assert(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I);
  assert(adef->model == M_GTRGAMMA || adef->model == M_PROTGAMMA);
   
  tr->rateHetModel = CAT;

  switch(adef->model)
    {
    case M_GTRGAMMA:
      adef->model = M_GTRCAT;	     
      break;
    case M_PROTGAMMA:
      adef->model = M_PROTCAT;
      break;
    default:
      assert(0);
    }
#ifdef _LOCAL_DATA
  tr->currentModel = adef->model;
#endif
  assignLikelihoodFunctions(tr, adef);
  allocNodex(tr, adef);  
}

static void gammaToParsimony(tree *tr, analdef *adef)
{
  /* TODO-PTHREADS check */

  freeNodex(tr);   	  

  assert(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I);
  assert(adef->model == M_GTRGAMMA || adef->model == M_PROTGAMMA);
   
  tr->rateHetModel = CAT;

  switch(adef->model)
    {
    case M_GTRGAMMA:
      adef->model = M_GTRCAT;	     
      break;
    case M_PROTGAMMA:
      adef->model = M_PROTCAT;
      break;
    default:
      assert(0);
    }
#ifdef _LOCAL_DATA
  tr->currentModel = adef->model;
#endif
}

static void gammaToParsimonyNoFree(tree *tr, analdef *adef)
{ 
  /* TODO PTHREADS check */

  assert(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I);
  assert(adef->model == M_GTRGAMMA || adef->model == M_PROTGAMMA);
   
  tr->rateHetModel = CAT;

  switch(adef->model)
    {
    case M_GTRGAMMA:
      adef->model = M_GTRCAT;	     
      break;
    case M_PROTGAMMA:
      adef->model = M_PROTCAT;
      break;
    default:
      assert(0);
    }
#ifdef _LOCAL_DATA
  tr->currentModel = adef->model;
#endif
}

static void singleBootstrap(tree *tr, int i, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int k;

  tr->treeID = i;
  tr->checkPointCounter = 0;
     
  if(i > 0)
    {
      makeboot(adef, tr);  	      
      initModel(tr, rdta, cdta, adef);                                   
    }
      
  getStartingTree(tr, adef);
  computeBIGRAPID(tr, adef);         

  if(adef->bootstrapBranchLengths)
    {
      if(adef->model == M_GTRGAMMA || adef->model == M_PROTGAMMA)     	    
	modOpt(tr, adef);		      	    	    
      else
	{	  	      
	  if(adef->useMixedModel)
	    {		  
	      tr->likelihood = unlikely;
	       
	      catToGamma(tr, adef);

	      initModel(tr, rdta, cdta, adef);	  	  	  
	      modOpt(tr, adef);			

	      gammaToParsimonyNoFree(tr, adef);
 		  	      
	      /*
		if(adef->model == M_GTRGAMMA)
		adef->model = M_GTRCAT;
		else
		adef->model = M_PROTCAT;
	      */
	    }	     	
	}
    }  
#ifndef PARALLEL
   printBootstrapResult(tr, adef, TRUE);
#endif
   freeNodex(tr);                         
   
   for(k = 1; k <= rdta->sites; k++)
     {	  
       rdta->wgt[k] = rdta->wgt2[k] = 1;	 	  
       cdta->aliaswgt[k] = 0;
     }
}

static void multipleBootstrap(tree *tr, int i, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int k, l, j, max;
  topolRELL_LIST *rl; 
  double bestLH;

  tr->treeID = i;
  tr->checkPointCounter = 0;
     
  if(i > 0)
    {
      makeboot(adef, tr);  	      
      initModel(tr, rdta, cdta, adef);                                   
    }
  
  rl = (topolRELL_LIST *)malloc(sizeof(topolRELL_LIST));
  initTL(rl, tr, adef->multiBoot);
  
  for(l = 0; l < adef->multiBoot; l++)
    {         
      getStartingTree(tr, adef);

      if(l == 0)
	{	
	  if(tr->rateHetModel == CAT && (adef->model == M_GTRCAT || adef->model == M_PROTCAT))
	    catToGamma(tr, adef); 	  
	  for(j = 0; j < tr->cdta->endsite; j++)
	    tr->cdta->wr[j] = tr->cdta->aliaswgt[j];
	  modOpt(tr, adef);
	  /* printf("ModOpt, replicate %d: %f\n", i, tr->likelihood); */
	  if(tr->rateHetModel == CAT && (adef->model == M_GTRCAT || adef->model == M_PROTCAT))
	    categorizeGeneric(tr, tr->start);       
	}

      if(tr->rateHetModel == CAT && (adef->model == M_GTRCAT || adef->model == M_PROTCAT))
	gammaToCat(tr, adef);	 


      evaluateGenericInitrav(tr, tr->start);
          

      computeBIGRAPIDMULTIBOOT(tr, adef);
      saveTL(rl, tr, l);       
      if(l < adef->multiBoot - 1)
	freeNodex(tr); 
    }  
  if(tr->rateHetModel == CAT && (adef->model == M_GTRCAT || adef->model == M_PROTCAT))
    catToGamma(tr, adef); 

  for(j = 0; j < tr->cdta->endsite; j++)
    tr->cdta->wr[j] = tr->cdta->aliaswgt[j];

  bestLH = unlikely;
  max = -1;

  for(l = 0; l < adef->multiBoot; l++)
    {      
      restoreTL(rl, tr, l);


      evaluateGenericInitrav(tr, tr->start);

      treeEvaluate(tr, 0.5);
      rl->t[l]->likelihood = unlikely;
      /* printf("TREE %d GAMMA %f\n", l, tr->likelihood); */
      if(tr->likelihood > bestLH)
	{
	  bestLH = tr->likelihood;
	  max = l;
	}
     
      saveTL(rl, tr, l);     
    }

  
  restoreTL(rl, tr, max);  


  evaluateGenericInitrav(tr, tr->start);
  
  /* printf("Best under GAMMA: %f\n", tr->likelihood); */

  if(adef->reallyThoroughBoot)
    {
      treeOptimizeThorough(tr, 1, 10);


      evaluateGenericInitrav(tr, tr->start);
    
      /* printf("AFTER THOROUGH: %f\n", tr->likelihood); */
    }
    
#ifndef PARALLEL
  printBootstrapResult(tr, adef, TRUE);
#endif
              
  for(k = 1; k <= rdta->sites; k++)
    {	  
      rdta->wgt[k] = rdta->wgt2[k] = 1;	 	  
      cdta->aliaswgt[k] = 0;
    }
  
  freeNodex(tr);   
  freeTL(rl, tr);   
  free(rl);  
}





/***************************** EXPERIMENTAL FUNCTIONS ********************************************************/




static int compareTopolRell(const void *p1, const void *p2)
{
  topolRELL **rc1 = (topolRELL **)p1;
  topolRELL **rc2 = (topolRELL **)p2;
 
  double i = (*rc1)->likelihood;
  double j = (*rc2)->likelihood;
  
  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}


void fixModelIndices(tree *tr, analdef *adef, int endsite)
{
  if(adef->useMultipleModel)
    {
      int i, model;     

      tr->partitionData[0].lower = 0;
     
      model = tr->model[0];
      i = 1;

      while(i < endsite)
	{
	  if(tr->model[i] != model)
	    {	      
	      tr->partitionData[model].upper = i;
	      tr->partitionData[model + 1].lower = i;
	      model = tr->model[i];
	    }
	  i++;
	}

      
      tr->partitionData[tr->NumberOfModels - 1].upper = endsite;
    }
  else
    {     
      tr->partitionData[0].lower = 0;
      tr->partitionData[0].upper = endsite;     
    }

  if(tr->mixedData)    
    calculateModelOffsets(tr);   
}

void reductionCleanup(tree *tr, analdef *adef, int *originalRateCategories, int *originalInvariant)
{
  int j;

  tr->cdta->endsite = tr->originalCrunchedLength;

  memcpy(tr->cdta->aliaswgt, tr->originalWeights, sizeof(int) * tr->cdta->endsite);
  memcpy(tr->model, tr->originalModel, sizeof(int) * tr->cdta->endsite);
  memcpy(tr->dataVector, tr->originalDataVector,  sizeof(int) * tr->cdta->endsite);

  memcpy(tr->cdta->rateCategory, originalRateCategories, sizeof(int) * tr->cdta->endsite);
  memcpy(tr->invariant,          originalInvariant,      sizeof(int) * tr->cdta->endsite); 
  
  for (j = 0; j < tr->originalCrunchedLength; j++) 
    {	
      double temp, wtemp;     
      temp = tr->cdta->patrat[originalRateCategories[j]];
      tr->cdta->wr[j]  = wtemp = temp * tr->cdta->aliaswgt[j];
      tr->cdta->wr2[j] = temp * wtemp;
    }                     
      
  if(!adef->computeELW)
    assert(adef->model == M_GTRCAT || adef->model == M_PROTCAT);
      
  memcpy(tr->rdta->y0, tr->rdta->yBUF, tr->rdta->numsp * tr->cdta->endsite * sizeof(char));  
      
  tr->cdta->endsite = tr->originalCrunchedLength;
  fixModelIndices(tr, adef, tr->originalCrunchedLength);    
#ifdef _LOCAL_DATA
  masterBarrier(THREAD_NEXT_REPLICATE, tr);
#endif
}


void makeboot (analdef *adef, tree *tr)
{
  int  i, nonzero, j, model, *weightBuffer, pos, w, endsite, l;
#ifdef PARALLEL
  long seed;
#endif       

  for(j = 0; j < tr->originalCrunchedLength; j++)
    tr->cdta->aliaswgt[j] = 0;

  if(tr->NumberOfModels > 1)
    {      	  
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  nonzero = 0;        	 

	  for (j = 0; j < tr->originalCrunchedLength; j++)  
	    {
	      if(tr->originalModel[j] == model)
		nonzero += tr->originalWeights[j];
	    }				          
	     
	  weightBuffer = (int *)calloc(nonzero, sizeof(int));	 
#ifdef PARALLEL 
	  seed = (long) gettimeSrand();
	  for (j = 0; j < nonzero; j++)
	    weightBuffer[(int) (nonzero*randum(& seed))]++;
#else  
	  for (j = 0; j < nonzero; j++)
	    weightBuffer[(int) (nonzero*randum(& adef->boot))]++;                  
#endif
	      
	  pos = 0;	      
		
	  for(j = 0; j < tr->originalCrunchedLength; j++) 
	    {
	      if(model == tr->originalModel[j])
		{
		  for(w = 0; w < tr->originalWeights[j]; w++)	  	  	 
		    {
		      tr->cdta->aliaswgt[j] += weightBuffer[pos];
		      pos++;		      
		    }				   
		}
	    }  
	  free(weightBuffer);	  
	}
    }
  else
    {
      weightBuffer = (int*)calloc(tr->fullSites, sizeof(int));
#ifdef PARALLEL 
      seed = (long) gettimeSrand();
      for (j = 0; j < tr->fullSites; j++)
	weightBuffer[(int) (tr->fullSites * randum(& seed))]++;
#else        
      for(j = 0; j < tr->fullSites; j++)
	weightBuffer[(int)(tr->fullSites * randum(& adef->boot))]++;      
#endif      
      
      pos = 0;
      for(j = 0; j < tr->originalCrunchedLength; j++)	  
	for(w = 0; w < tr->originalWeights[j]; w++)	  	  	 
	  {
	    tr->cdta->aliaswgt[j] += weightBuffer[pos];
	    pos++;
	  }
      
      free(weightBuffer);
    }            


  endsite = 0;
  for (j = 0; j < tr->originalCrunchedLength; j++) 
    {	      
      if(tr->cdta->aliaswgt[j] > 0)
	endsite++;      
    }               

 
                  
  for(i = 0; i < tr->rdta->numsp; i++)
    {
      char *yPos    = &(tr->rdta->y0[tr->originalCrunchedLength * i]);
      char *origSeq = &(tr->rdta->yBUF[tr->originalCrunchedLength * i]);
      int l, j;
      
      for(j = 0, l = 0; j < tr->originalCrunchedLength; j++)
	{	
	  if(tr->cdta->aliaswgt[j] > 0)	  
	    {	     
	      yPos[l++] = origSeq[j];
	    }
	}               
    }
      
  for(j = 0, l = 0; j < tr->originalCrunchedLength; j++)
    {     
      if(tr->cdta->aliaswgt[j] > 0)	
	{
	  tr->cdta->aliaswgt[l]     = tr->cdta->aliaswgt[j];	
	  tr->model[l]              = tr->originalModel[j];
	  tr->dataVector[l]         = tr->originalDataVector[j];
	  l++;
	}
    }

  tr->cdta->endsite = endsite; 
  fixModelIndices(tr, adef, endsite);   
#ifdef _LOCAL_DATA
  masterBarrier(THREAD_NEXT_REPLICATE, tr);
#endif  
}




void computeNextReplicate(tree *tr, analdef *adef, int *originalRateCategories, int *originalInvariant)
{ 
  int pos, nonzero, j, model, w;   
  int *weightBuffer, endsite;                
  int *weights, i, l;  

  assert(adef->rapidBoot != 0);

  for(j = 0; j < tr->originalCrunchedLength; j++)
    tr->cdta->aliaswgt[j] = 0;

  if(tr->NumberOfModels > 1)
    {      	  
      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  nonzero = 0;        	 

	  for (j = 0; j < tr->originalCrunchedLength; j++)  
	    {
	      if(tr->originalModel[j] == model)
		nonzero += tr->originalWeights[j];
	    }				          
	     
	  weightBuffer = (int *)calloc(nonzero, sizeof(int));	 
	      
	  for (j = 0; j < nonzero; j++)
	    weightBuffer[(int) (nonzero*randum(& adef->rapidBoot))]++;                  
	      
	  pos = 0;	      
		
	  for(j = 0; j < tr->originalCrunchedLength; j++) 
	    {
	      if(model == tr->originalModel[j])
		{
		  for(w = 0; w < tr->originalWeights[j]; w++)	  	  	 
		    {
		      tr->cdta->aliaswgt[j] += weightBuffer[pos];
		      pos++;		      
		    }				   
		}
	    }  
	  free(weightBuffer);	  
	}
    }
  else
    {
      weightBuffer = (int*)calloc(tr->fullSites, sizeof(int));
      
      for(j = 0; j < tr->fullSites; j++)
	weightBuffer[(int)(tr->fullSites * randum(& adef->rapidBoot))]++;      
      
      pos = 0;
      for(j = 0; j < tr->originalCrunchedLength; j++)	  
	for(w = 0; w < tr->originalWeights[j]; w++)	  	  	 
	  {
	    tr->cdta->aliaswgt[j] += weightBuffer[pos];
	    pos++;
	  }
      
      free(weightBuffer);
    }            

  

  endsite = 0;
  
  for (j = 0; j < tr->originalCrunchedLength; j++) 
    {	
      double temp, wtemp;
      if(tr->cdta->aliaswgt[j] > 0)
	endsite++;
      temp = tr->cdta->patrat[originalRateCategories[j]];
      tr->cdta->wr[j]  = wtemp = temp * tr->cdta->aliaswgt[j];
      tr->cdta->wr2[j] = temp * wtemp;
    }          
  
  weights = tr->cdta->aliaswgt;

#ifndef _VINCENT  
  if(!adef->computeELW && !(adef->mode == MEHRING_ALGO))
    assert(adef->model == M_GTRCAT || adef->model == M_PROTCAT);
#endif      

  for(i = 0; i < tr->rdta->numsp; i++)
    {     
      char *yPos    = &(tr->rdta->y0[tr->originalCrunchedLength * i]);
      char *origSeq = &(tr->rdta->yBUF[tr->originalCrunchedLength * i]);
      int l, j;
      
      for(j = 0, l = 0; j < tr->originalCrunchedLength; j++)      
	if(tr->cdta->aliaswgt[j] > 0)	  	    
	  yPos[l++] = origSeq[j];	                   
    }

  for(j = 0, l = 0; j < tr->originalCrunchedLength; j++)
    {     
      if(weights[j])	
	{
	  tr->cdta->aliaswgt[l]     = tr->cdta->aliaswgt[j];
	  tr->cdta->wr[l]           = tr->cdta->wr[j];
	  tr->cdta->wr2[l]          = tr->cdta->wr2[j];
	  tr->model[l]              = tr->originalModel[j];
	  tr->dataVector[l]         = tr->originalDataVector[j];
	  tr->cdta->rateCategory[l] = originalRateCategories[j];           
	  tr->invariant[l]          = originalInvariant[j];
	  l++;
	}
    }

  tr->cdta->endsite = endsite;
  fixModelIndices(tr, adef, endsite);
#ifdef _LOCAL_DATA
  masterBarrier(THREAD_NEXT_REPLICATE, tr);
#endif              
}

void quickOpt(tree *tr, analdef *adef)
{	   	  	    
  int l;	  
  double sl, ol;	  	 
    	  
  /* printf("START %f\n", tr->likelihood); */
	   
  l = 0;
	  
  do
    {
      sl = tr->likelihood;
      quickAndDirtyOptimization(tr, adef);
      treeEvaluate(tr, 0.25);
      /* printf("%d %f\n", l, tr->likelihood); */
      ol = tr->likelihood;
      l++;
    }
  while((ol > sl) && (ol - sl > 1.0));  
} 
  


typedef struct {
  double *tipVectorDNA;
  double *tipVectorAA;    
  double *EV_DNA;         
  double *EV_AA;           	  
  double *EI_DNA;
  double *EI_AA;           	  
  double *EIGN_DNA;       
  double *EIGN_AA;        	  
  double *frequencies_DNA;
  double *frequencies_AA;  	  
  double *initialRates_DNA;
  double *initialRates_AA;    
} modelParams;


static void allocParams(modelParams *params, tree *tr)
{
  params->tipVectorDNA    = (double *)malloc(tr->NumberOfModels * 64 * sizeof(double));
  params->tipVectorAA     = (double *)malloc(tr->NumberOfModels * 460 * sizeof(double));
  params->EV_DNA          = (double *)malloc(tr->NumberOfModels * 16 * sizeof(double));
  params->EV_AA           = (double *)malloc(tr->NumberOfModels * 400 * sizeof(double));	  
  params->EI_DNA          = (double *)malloc(tr->NumberOfModels * 12 * sizeof(double));
  params->EI_AA           = (double *)malloc(tr->NumberOfModels * 380 * sizeof(double));  	  
  params->EIGN_DNA        = (double *)malloc(tr->NumberOfModels * 3 * sizeof(double));
  params->EIGN_AA         = (double *)malloc(tr->NumberOfModels * 19  * sizeof(double));
  params->frequencies_DNA = (double *)malloc(tr->NumberOfModels * 4 * sizeof(double));
  params->frequencies_AA  = (double *)malloc(tr->NumberOfModels * 20  * sizeof(double));
  params->initialRates_DNA = (double *)malloc(tr->NumberOfModels * 5 * sizeof(double));
  params->initialRates_AA  = (double *)malloc(tr->NumberOfModels * 190 * sizeof(double));	    
}

static void freeParams(modelParams *params)
{
  free(params->tipVectorDNA);
  free(params->tipVectorAA);
  free(params->EV_DNA);
  free(params->EV_AA);
  free(params->EI_DNA);
  free(params->EI_AA);
  free(params->EIGN_DNA);
  free(params->EIGN_AA);
  free(params->frequencies_DNA);
  free(params->frequencies_AA);
  free(params->initialRates_DNA);
  free(params->initialRates_AA);
}

static void storeParams(modelParams *params, tree *tr)
{
   memcpy(params->tipVectorDNA,     tr->tipVectorDNA,     tr->NumberOfModels * 64 * sizeof(double));
   memcpy(params->tipVectorAA,      tr->tipVectorAA,      tr->NumberOfModels * 460 * sizeof(double));
   memcpy(params->EV_DNA,           tr->EV_DNA,           tr->NumberOfModels * 16 * sizeof(double));
   memcpy(params->EV_AA,            tr->EV_AA,            tr->NumberOfModels * 400 * sizeof(double));
   memcpy(params->EI_DNA,           tr->EI_DNA,           tr->NumberOfModels * 12 * sizeof(double));
   memcpy(params->EI_AA,            tr->EI_AA,            tr->NumberOfModels * 380 * sizeof(double));
   memcpy(params->EIGN_DNA,         tr->EIGN_DNA,         tr->NumberOfModels * 3 * sizeof(double));
   memcpy(params->EIGN_AA,          tr->EIGN_AA,          tr->NumberOfModels * 19  * sizeof(double));
   memcpy(params->frequencies_DNA,  tr->frequencies_DNA,  tr->NumberOfModels * 4 * sizeof(double));
   memcpy(params->frequencies_AA,   tr->frequencies_AA,   tr->NumberOfModels * 20  * sizeof(double));	  
   memcpy(params->initialRates_DNA, tr->initialRates_DNA, tr->NumberOfModels * 5 * sizeof(double));
   memcpy(params->initialRates_AA,  tr->initialRates_AA,  tr->NumberOfModels * 190 * sizeof(double));
}

static void loadParams(modelParams *params, tree *tr)
{
   memcpy(tr->tipVectorDNA,     params->tipVectorDNA,     tr->NumberOfModels * 64 * sizeof(double));
   memcpy(tr->tipVectorAA,      params->tipVectorAA,      tr->NumberOfModels * 460 * sizeof(double));
   memcpy(tr->EV_DNA,           params->EV_DNA,           tr->NumberOfModels * 16 * sizeof(double));
   memcpy(tr->EV_AA,            params->EV_AA,            tr->NumberOfModels * 400 * sizeof(double));
   memcpy(tr->EI_DNA,           params->EI_DNA,           tr->NumberOfModels * 12 * sizeof(double));
   memcpy(tr->EI_AA,            params->EI_AA,            tr->NumberOfModels * 380 * sizeof(double));
   memcpy(tr->EIGN_DNA,         params->EIGN_DNA,         tr->NumberOfModels * 3 * sizeof(double));
   memcpy(tr->EIGN_AA,          params->EIGN_AA,          tr->NumberOfModels * 19  * sizeof(double));
   memcpy(tr->frequencies_DNA,  params->frequencies_DNA,  tr->NumberOfModels * 4 * sizeof(double));
   memcpy(tr->frequencies_AA,   params->frequencies_AA,   tr->NumberOfModels * 20  * sizeof(double));	  
   memcpy(tr->initialRates_DNA, params->initialRates_DNA, tr->NumberOfModels * 5 * sizeof(double));
   memcpy(tr->initialRates_AA,  params->initialRates_AA,  tr->NumberOfModels * 190 * sizeof(double));
}


#ifndef PARALLEL




void doAllInOne(tree *tr, analdef *adef)
{
  int i, j, n, sites, bestIndex, model, bootstrapsPerformed;
  double loopTime; 
  int      *originalRateCategories;
  int      *originalInvariant;
  int      slowSearches, fastEvery = 5;
  topolRELL_LIST *rl;  
  double bestLH, mlTime, overallTime;  
  long radiusSeed = adef->rapidBoot;
  FILE *infoFile, *f;
  char bestTreeFileName[1024];  
  BL *b = (BL *)NULL;
  boolean bootStopIt = FALSE;
  double pearsonAverage = 0.0;
  modelParams *catParams   = (modelParams *)malloc(sizeof(modelParams));
  modelParams *gammaParams = (modelParams *)malloc(sizeof(modelParams));

  allocParams(catParams,   tr);
  allocParams(gammaParams, tr);

  n = adef->multipleRuns;

  if(adef->bootStopping)
    {
      b = (BL *)malloc(sizeof(BL));
      
      b->count = 0;
      b->n = FC_INIT * tr->mxtips;   
      b->treeVectorLength = (n / BITS_BYTE) + 1;      
      b->b = (bipList *)malloc(sizeof(bipList) *  b->n);

      for(i = 0; i < b->n; i++)
	{     
	  b->b[i].length = 0;
	  b->b[i].isSet = (unsigned char*)calloc(b->treeVectorLength, sizeof(unsigned char));
	  b->b[i].entries = (int *)NULL;
	}  
    }

  rl = (topolRELL_LIST *)malloc(sizeof(topolRELL_LIST));
  initTL(rl, tr, n);
     
  originalRateCategories = (int*)malloc(tr->cdta->endsite * sizeof(int));      
  originalInvariant      = (int*)malloc(tr->cdta->endsite * sizeof(int));

  sites = tr->cdta->endsite;             

  if(adef->model == M_PROTGAMMA || adef->model == M_GTRGAMMA || 
     tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
    {
      assert(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I);
      assert(adef->model == M_GTRGAMMA || adef->model == M_PROTGAMMA);

      tr->rateHetModel = CAT;

      switch(adef->model)
	{
	case M_GTRGAMMA:
	  adef->model = M_GTRCAT;	     
	  break;
	case M_PROTGAMMA:
	  adef->model = M_PROTCAT;
	  break;
	default:
	  assert(0);
	}

      initModel(tr, tr->rdta, tr->cdta, adef);

      printf("\nSwitching from GAMMA to CAT for rapid Bootstrap, final ML search will be conducted under the %s model you specified\n", 
	     (adef->useInvariant)?"GAMMA+P-Invar":"GAMMA");
      infoFile = fopen(infoFileName, "a");
      fprintf(infoFile, "\nSwitching from GAMMA to CAT for rapid Bootstrap, final ML search will be conducted under the %s model you specified\n", 
	      (adef->useInvariant)?"GAMMA+P-Invar":"GAMMA");
      fclose(infoFile);           
    }
 
  for(i = 0; i < n && !bootStopIt; i++)
    {                
      tr->treeID = i;
      tr->checkPointCounter = 0;
         
      loopTime = gettime();      
     
      if(i % 10 == 0)
	{
	  if(i > 0)
	    {
	      freeNodex(tr);
	      reductionCleanup(tr, adef, originalRateCategories, originalInvariant);
	    }

	  

	  if(adef->grouping || adef->constraint)
	    {
	      FILE *f = fopen(tree_file, "r");	

	      assert(adef->restart);
	      partCount = 0;
	      if (! treeReadLenMULT(f, tr, adef))
		exit(-1);
	      /*printf("Constraint read: %d\n", partCount);*/
	      fclose(f);
	    }
	  else
	    makeParsimonyTree(tr, adef);
	  allocNodex(tr, adef); 
	  tr->likelihood = unlikely;
	  if(i == 0)
	    {
	      double t;
	          

	      onlyInitrav(tr, tr->start);

	      /* debug */

	      /*evaluateGeneric(tr, tr->start);
		printf("LH %f\n", tr->likelihood);*/
	      
	      /* debug-end */

	      treeEvaluate(tr, 1);	     	
	      /*printf("LH2 %f\n", tr->likelihood);*/
	      

	      t = gettime();    

	      quickOpt(tr, adef);	    
 
	      printf("\nTime for BS model parameter optimization %f\n", gettime() - t);
	      infoFile = fopen(infoFileName, "a");
	      fprintf(infoFile, "\nTime for BS model parameter optimization %f\n", gettime() - t);
	      fclose(infoFile);
	      
	      memcpy(originalRateCategories, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);
	      memcpy(originalInvariant,      tr->invariant,          sizeof(int) * tr->cdta->endsite);

	      if(adef->bootstrapBranchLengths)
		{
		  storeParams(catParams, tr);
		  assert(tr->cdta->endsite == tr->originalCrunchedLength);
		  catToGamma(tr, adef);
		  modOpt(tr, adef);
		  storeParams(gammaParams, tr);
		  gammaToCat(tr, adef);
		  loadParams(catParams, tr);		  
		}
	    }	  	  
	}

      computeNextReplicate(tr, adef, originalRateCategories, originalInvariant); 
      resetBranches(tr);

      evaluateGenericInitrav(tr, tr->start);
    
      treeEvaluate(tr, 1);    	             

      computeBOOTRAPID(tr, adef, &radiusSeed);                        	  
      saveTL(rl, tr, i);

      if(adef->bootstrapBranchLengths)
	{
	  double lh = tr->likelihood;
	  int    endsite;
	 
	  loadParams(gammaParams, tr);
     
	  
	  endsite = tr->cdta->endsite;
	  tr->cdta->endsite = tr->originalCrunchedLength;
	  catToGamma(tr, adef);
	  tr->cdta->endsite = endsite;

	  resetBranches(tr);
	  treeEvaluate(tr, 2.0);
	  
	  endsite = tr->cdta->endsite;
	  tr->cdta->endsite = tr->originalCrunchedLength;
	  gammaToCat(tr, adef);
	  tr->cdta->endsite = endsite;	 	    
	
	  loadParams(catParams, tr);

	 
	  tr->likelihood = lh;
	}
      
      printBootstrapResult(tr, adef, TRUE); 

      loopTime = gettime() - loopTime; 
      writeInfoFile(adef, tr, loopTime); 
     
      if(adef->bootStopping)	
	bootStopIt = bootStop(tr, b, i, &pearsonAverage, adef->bootstopCutoff);
    }   

  bootstrapsPerformed = i;

  freeParams(catParams);
  free(catParams);

  freeParams(gammaParams);
  free(gammaParams);

  if(adef->bootStopping)
    {
      for(j = 0; j < b->n; j++)
	{
	  if(b->b[j].entries != (int *)NULL)
	    free(b->b[j].entries);
	  free(b->b[j].isSet);	 
	}
      free(b->b);
      free(b);
    }

  if(!adef->allInOne)
    {      
      double t = gettime() - masterTime;
      infoFile = fopen(infoFileName, "a");

      printf("\n\n");
      fprintf(infoFile, "\n\n");

      if(adef->bootStopping)
	{
	  printf("Stopped Rapid BS search after %d replicates with Bootstopping criterion\n", bootstrapsPerformed);
	  printf("Pearson Average of 100 random splits: %f\n", pearsonAverage);
	  fprintf(infoFile, "Stopped Rapid BS search after %d replicates with Bootstopping criterion\n", 
		  bootstrapsPerformed);	  
	  fprintf(infoFile, "Pearson Average of 100 random splits: %f\n", pearsonAverage);
	}

      printf("Overall Time for %d Rapid Bootstraps %f\n", bootstrapsPerformed, t);
      fprintf(infoFile, "Overall Time for %d Rapid Bootstraps %f\n", bootstrapsPerformed, t);

      printf("Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bootstrapsPerformed)));
      fprintf(infoFile, "Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bootstrapsPerformed)));

      printf("All %d bootstrapped trees written to: %s\n", bootstrapsPerformed, bootstrapFileName);
      fprintf(infoFile, "All %d bootstrapped trees written to: %s\n", bootstrapsPerformed, bootstrapFileName);           
     
      fclose(infoFile);
      exit(0);
    }
  else
    {
      double t = gettime() - masterTime;
      infoFile = fopen(infoFileName, "a");

      printf("\n\n");
      fprintf(infoFile, "\n\n");

      if(adef->bootStopping)
	{
	  printf("Stopped Rapid BS search after %d replicates with Bootstopping criterion\n", bootstrapsPerformed);
	  printf("Pearson Average of 100 random splits: %f\n", pearsonAverage);
	  fprintf(infoFile, "Stopped Rapid BS search after %d replicates with Bootstopping criterion\n", 
		  bootstrapsPerformed);	 
	  fprintf(infoFile, "Pearson Average of 100 random splits: %f\n", pearsonAverage);
	}

      printf("Overall Time for %d Rapid Bootstraps %f\n", bootstrapsPerformed, t);
      fprintf(infoFile, "Overall Time for %d Rapid Bootstraps %f\n", bootstrapsPerformed, t);

      printf("Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bootstrapsPerformed)));                
      fprintf(infoFile, "Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bootstrapsPerformed)));	           
     
      fclose(infoFile);
    }
  
  /* ML-search */ 

  mlTime = gettime();

  infoFile = fopen(infoFileName, "a");
  printf("\nStarting ML Search ...\n\n");
  fprintf(infoFile, "\nStarting ML Search ...\n\n");
  fclose(infoFile);

  /***CLEAN UP reduction stuff */  

  reductionCleanup(tr, adef, originalRateCategories, originalInvariant);  

  /****/
 
  catToGamma(tr, adef);     	   
  
  restoreTL(rl, tr, 0);

  resetBranches(tr);

  evaluateGenericInitrav(tr, tr->start);   
  
  modOpt(tr, adef);
 
  evaluateGenericInitrav(tr, tr->start);   

  categorizeGeneric(tr, tr->start);        

  gammaToCat(tr, adef);      

#ifdef _DEBUG_AA
  {
    double *vector = (double *)malloc(sizeof(double) * tr->cdta->endsite);

    for(i = 0; i < tr->NumberOfCategories; i++)
      printf("%d: %f\n", i, tr->cdta->patrat[i]);

   
    evaluateGenericInitrav(tr, tr->start);
    
    printf("INITRAV %f\n", tr->likelihood);
    
    evaluateGenericVector (tr, tr->start, vector);
    for(i = 0; i < tr->cdta->endsite; i++)
      printf("%d %f %f %f %d\n", tr->cdta->rateCategory[i], vector[i], tr->cdta->wr[i], tr->cdta->wr2[i], 
	     tr->cdta->aliaswgt[i]);
   
  }
#endif


  if(adef->bootStopping)
    {
      if(bootstrapsPerformed <= 100)
	fastEvery = 5;
      else
	fastEvery = bootstrapsPerformed / 20;     
    }
  else    
    fastEvery = 5;    



#ifdef _DEBUG_AA
  for(i = 0; i < bootstrapsPerformed; i++)
    {
      restoreTL(rl, tr, i);
      resetBranches(tr);	 
      evaluateGenericInitrav(tr, tr->start);
      printf("INI %d %f\n", i, tr->likelihood);
      treeEvaluate(tr, 1); 	
      assert(!isnan(tr->likelihood));
      printf("EVAL %d %f\n", i, tr->likelihood);
    }
#endif

  for(i = 0; i < bootstrapsPerformed; i++)
    {            
      rl->t[i]->likelihood = unlikely;
    
      if(i % fastEvery == 0)
	{	  	 
	  restoreTL(rl, tr, i); 	 	    	   	
	  
	  /* DEBUG: does that help? */
	  resetBranches(tr);	 

	  evaluateGenericInitrav(tr, tr->start);
	  
	  treeEvaluate(tr, 1); 		 
	  
	  optimizeRAPID(tr, adef);	  			         

	  saveTL(rl, tr, i);   	 
	}    
    }     

  infoFile = fopen(infoFileName, "a");
  printf("Fast ML optimization finished\n\n");
  fprintf(infoFile, "Fast ML optimization finished\n\n");
  fclose(infoFile);  
  qsort(&(rl->t[0]), bootstrapsPerformed, sizeof(topolRELL*), compareTopolRell);
     
  catToGamma(tr, adef);  
  
  restoreTL(rl, tr, 0);

  resetBranches(tr);

  evaluateGenericInitrav(tr, tr->start);

#ifdef _DEBUG_AA  
  printf("I1: %f\n", tr->likelihood);
#endif

  modOpt(tr, adef);

#ifdef _DEBUG_AA
  printf("I2: %f\n", tr->likelihood);
#endif

  evaluateGenericInitrav(tr, tr->start);
  
  categorizeGeneric(tr, tr->start);   
       	  	 	      
  gammaToCat(tr, adef);    
  
  slowSearches = bootstrapsPerformed / 5;
  if(bootstrapsPerformed % 5 != 0)
    slowSearches++;

  slowSearches  = MIN(slowSearches, 10); 


  for(i = 0; i < slowSearches; i++)
    {           
      restoreTL(rl, tr, i);     
      rl->t[i]->likelihood = unlikely;  

      evaluateGenericInitrav(tr, tr->start);

#ifdef _DEBUG_AA
      printf("S1: %f\n", tr->likelihood);
#endif

      treeEvaluate(tr, 1.0);   

#ifdef _DEBUG_AA
      printf("S2: %f\n", tr->likelihood);
#endif

      thoroughOptimization(tr, adef, rl, i);    
         
#ifdef _DEBUG_AA
      printf("S3: %f\n", tr->likelihood);
#endif
   }

  infoFile = fopen(infoFileName, "a");
  printf("Slow ML optimization finished\n\n");
  fprintf(infoFile, "Slow ML optimization finished\n\n");
  fclose(infoFile);

  /*************************************************************************************************************/  
  
  catToGamma(tr, adef);    

  bestIndex = -1;
  bestLH = unlikely;
    
  for(i = 0; i < slowSearches; i++)
    { 
     
#ifdef _DEBUG_AA
      printf("TL: rlMax %d rlmembers %d\n", rl->max, rl->members);
#endif

      restoreTL(rl, tr, i);
      resetBranches(tr);

      evaluateGenericInitrav(tr, tr->start);

#ifdef _DEBUG_AA
      printf("R1: %f\n", tr->likelihood);
#endif

      treeEvaluate(tr, 2);

#ifdef _DEBUG_AA
      printf("R2: %f\n", tr->likelihood);
#endif

      infoFile = fopen(infoFileName, "a"); 

#ifdef _DEBUG_AA
      printf("R3: %f\n", tr->likelihood);
#endif

      printf("Slow ML Search %d Likelihood: %f\n", i, tr->likelihood);
      fprintf(infoFile, "Slow ML Search %d Likelihood: %f\n", i, tr->likelihood);     
      fclose(infoFile);
      if(tr->likelihood > bestLH)
	{
	  bestLH = tr->likelihood;
	  bestIndex = i;
	}
    }
  
  restoreTL(rl, tr, bestIndex);
  resetBranches(tr);

  evaluateGenericInitrav(tr, tr->start);
 

  treeEvaluate(tr, 2); 
         
  Thorough = 1;
  tr->doCutoff = FALSE;  
  
  treeOptimizeThorough(tr, 1, 10);
  modOpt(tr, adef);
 

  infoFile = fopen(infoFileName, "a"); 
  printf("\nFinal ML Optimization Likelihood: %f\n", tr->likelihood);
  fprintf(infoFile, "\nFinal ML Optimization Likelihood: %f\n", tr->likelihood);
  
  printf("\nModel Information:\n\n");
  fprintf(infoFile, "\nModel Information:\n\n");

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

      fprintf(infoFile, "Model Parameters of Partition %d, Name: %s, Type of Data: %s\n", 
	      model, tr->partitionData[model].partitionName, typeOfData);
      fprintf(infoFile, "alpha: %f\n", tr->alphas[model]);

      printf("Model Parameters of Partition %d, Name: %s, Type of Data: %s\n", 
	     model, tr->partitionData[model].partitionName, typeOfData);
      printf("alpha: %f\n", tr->alphas[model]);
      
      if(adef->useInvariant)
	{
	  fprintf(infoFile, "invar: %f\n", tr->invariants[model]);    
	  printf("invar: %f\n", tr->invariants[model]);    
	}
                 
      if(adef->perGeneBranchLengths)
	tl = treeLength(tr, model);
      else
	tl = treeLength(tr, 0);

      fprintf(infoFile, "Tree-Length: %f\n", tl);    
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
	      {
		fprintf(infoFile, "rate %s: %f\n", names[k], tr->initialRates_DNA[model * DNA_RATES + k]);			    
		printf("rate %s: %f\n", names[k], tr->initialRates_DNA[model * DNA_RATES + k]);
	      }
	
	    fprintf(infoFile, "rate %s: %f\n", names[5], 1.0);
	    printf("rate %s: %f\n", names[5], 1.0);
	  }      
	  break;
	default:
	  assert(0);
	}

      fprintf(infoFile, "\n");
      printf("\n");
    }		    		  
  

  fclose(infoFile); 
  
  
  strcpy(bestTreeFileName, workdir); 
  strcat(bestTreeFileName, "RAxML_bestTree.");
  strcat(bestTreeFileName,         run_id);
   
  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH);
  f = fopen(bestTreeFileName, "w");
  fprintf(f, "%s", tr->tree_string);
  fclose(f);

  if(adef->perGeneBranchLengths)
    printTreePerGene(tr, adef, bestTreeFileName, "w");

  infoFile = fopen(infoFileName, "a"); 
  overallTime = gettime() - masterTime;
  mlTime    = gettime() - mlTime;
  printf("\nML search took %f secs or %f hours\n", mlTime, mlTime / 3600.0);
  fprintf(infoFile, "\nML search took %f secs or %f hours\n", mlTime, mlTime / 3600.0);
  printf("\nCombined Bootstrap and ML search took %f secs or %f hours\n", overallTime, overallTime / 3600.0);
  fprintf(infoFile, "\nCombined Bootstrap and ML search took %f secs or %f hours\n", overallTime, overallTime / 3600.0);  

  printf("\nDrawing Bootstrap Support Values on best-scoring ML tree ...\n\n");
  fprintf(infoFile, "Drawing Bootstrap Support Values on best-scoring ML tree ...\n\n");  

  fclose(infoFile);
    
  freeTL(rl, tr);   
  free(rl);       
  
  calcBipartitions(tr, adef, bestTreeFileName, bootstrapFileName);  

  overallTime = gettime() - masterTime;
  printf("Program execution info written to %s\n", infoFileName);
  printf("All %d bootstrapped trees written to: %s\n\n", adef->multipleRuns, bootstrapFileName);
  printf("Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
  if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
    printf("Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n", bestTreeFileName,  bestTreeFileName, 
	   tr->NumberOfModels - 1);    
  printf("Best-scoring ML tree with support values written to: %s\n\n", bipartitionsFileName);
  printf("Overall execution time for full ML analysis: %f secs or %f hours or %f days\n\n", overallTime, overallTime/3600.0, overallTime/86400.0);
  
  infoFile = fopen(infoFileName, "a");
  fprintf(infoFile, "All %d bootstrapped trees written to: %s\n\n", adef->multipleRuns, bootstrapFileName);
  fprintf(infoFile, "Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
  if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
    fprintf(infoFile, "Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n", bestTreeFileName,  bestTreeFileName, 
	   tr->NumberOfModels - 1);    
  fprintf(infoFile, "Best-scoring ML tree with support values written to: %s\n\n", bipartitionsFileName);
  fprintf(infoFile, "Overall execution time for full ML analysis: %f secs or %f hours or %f days\n\n", overallTime, overallTime/3600.0, overallTime/86400.0);
  fclose(infoFile);
   
  exit(0); 
}

#ifdef _VINCENT

void doAllInOneVincent(tree *tr, analdef *adef)
{
  int i, j, n, sites, bestIndex, model, bootstrapsPerformed;
  double loopTime; 
  int      *originalRateCategories;
  int      slowSearches, fastEvery = 5;
  topolRELL_LIST *rl;  
  double bestLH, mlTime, overallTime;  
  long radiusSeed = adef->rapidBoot;
  FILE *infoFile, *f;
  char bestTreeFileName[1024];  
  BL *b = (BL *)NULL;
  boolean bootStopIt = FALSE;
  double pearsonAverage = 0.0;
  boolean optimizeModel = adef->optimizeBSmodel;

  n = adef->multipleRuns;   
     
  originalRateCategories = (int*)calloc(tr->cdta->endsite, sizeof(int));      

  sites = tr->cdta->endsite;             

  assert((adef->model == M_PROTGAMMA || adef->model == M_GTRGAMMA) && 
	 (tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I));         
  initModel(tr, tr->rdta, tr->cdta, adef);   
  
  makeParsimonyTree(tr, adef);

  allocNodex(tr, adef);

  if(!optimizeModel)
    {      
      modOpt(tr, adef);
      printf("Initial Model Optimization: %f\n", tr->likelihood);
    }

  for(i = 0; i < n; i++)
    {                
      tr->treeID = i;
      tr->checkPointCounter = 0;
      
      loopTime = gettime();                     

      computeNextReplicate(tr, adef, originalRateCategories);           

      
      if(i % 10 == 0)
	{	 
	  freeNodex(tr);	       
	  makeParsimonyTree(tr, adef);
	  allocNodex(tr, adef);	     
	}
      else
	{
	  freeNodex(tr);	       	  
	  allocNodex(tr, adef);
	}
           
      if(optimizeModel)
	{	  	 
	  initModel(tr, tr->rdta, tr->cdta, adef);
	  modOpt(tr, adef);       
	  /*printf("Bootstrap[%d]: Intermediate Model Optimization: %f\n", i, tr->likelihood);*/
	}
      else
	{
	  resetBranches(tr);
	  tr->likelihood = unlikely;
	  onlyInitrav(tr, tr->start);
	  treeEvaluate(tr, 1);
	  /*printf("Bootstrap[%d]: Without Model Optimization: %f\n", i, tr->likelihood);*/
	}

      computeBOOTRAPID(tr, adef, &radiusSeed);                        	        
      printBootstrapResult(tr, adef, TRUE);
 
      loopTime = gettime() - loopTime; 

      writeInfoFile(adef, tr, loopTime);            
    }   

  bootstrapsPerformed = i;

  {      
    double t = gettime() - masterTime;
    infoFile = fopen(infoFileName, "a");
    
    printf("\n\n");
    fprintf(infoFile, "\n\n");       

    printf("Overall Time for %d Rapid Bootstraps %f\n", bootstrapsPerformed, t);
    fprintf(infoFile, "Overall Time for %d Rapid Bootstraps %f\n", bootstrapsPerformed, t);
    
    printf("Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bootstrapsPerformed)));
    fprintf(infoFile, "Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bootstrapsPerformed)));
    
    printf("All %d bootstrapped trees written to: %s\n", bootstrapsPerformed, bootstrapFileName);
    fprintf(infoFile, "All %d bootstrapped trees written to: %s\n", bootstrapsPerformed, bootstrapFileName);           
    
    fclose(infoFile);
    exit(0);
  }
}
 


#endif

/*******************************************EXPERIMENTAL FUNCTIONS END *****************************************************/





void doBootstrap(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int i, n;
  double loopTime;
  n = adef->multipleRuns;
          
  for(i = 0; i < n; i++)
    {    
      loopTime = gettime();
        
      if(adef->multiBoot < 2)
	singleBootstrap(tr, i, adef, rdta, cdta);     
      else
	multipleBootstrap(tr, i, adef, rdta, cdta);
           	                         
      loopTime = gettime() - loopTime;
      
      writeInfoFile(adef, tr, loopTime);           
    }      
}


void doInference(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int i, n;
  double loopTime;

  n = adef->multipleRuns;
       
  for(i = 0; i < n; i++)
    {                
      tr->treeID = i;
      tr->checkPointCounter = 0;
         
      loopTime = gettime();
               
      if(i > 0)	       
	initModel(tr, rdta, cdta, adef);             
             
      getStartingTree(tr, adef);                      
      
      computeBIGRAPID(tr, adef);                     

      if(adef->model == M_GTRGAMMA || adef->model == M_PROTGAMMA || 
	 tr->rateHetModel == GAMMA ||  tr->rateHetModel == GAMMA_I)
	{	 
	  modOpt(tr, adef);	
	  printLog(tr, adef, TRUE);
	  printResult(tr, adef, TRUE);
	  loopTime = gettime() - loopTime;	 
	  tr->likelihoods[i] = tr->likelihood;
	  freeNodex(tr); 
	}
      else
	{	  
	  if(adef->useMixedModel)
	    {
	      tr->likelihood = unlikely;	      	     

	      catToGamma(tr, adef);

	      initModel(tr, rdta, cdta, adef);	  	  	  
	      modOpt(tr, adef);	
	      printLog(tr, adef, TRUE);
	      printResult(tr, adef, TRUE);
	      loopTime = gettime() - loopTime;	
	      tr->likelihoods[i] = tr->likelihood;
	 
	      gammaToParsimony(tr, adef);
	    }
	  else
	    {
	      loopTime = gettime() - loopTime;       
	      tr->likelihoods[i] = tr->likelihood;
	      freeNodex(tr);
	    }	
	}

      writeInfoFile(adef, tr, loopTime);                                     
    } 
}
#else



#include <mpi.h>

extern int processID;
extern int numOfWorkers;

static void sendTree(tree *tr, analdef *adef, double t, boolean finalPrint, int tag)
{
  int bufferSize, i, bufCount;
  double *buffer;
  char *tree_ptr;

  bufferSize = tr->treeStringLength + 4 + tr->NumberOfModels + tr->NumberOfModels;

  buffer = (double *)malloc(sizeof(double) * bufferSize);
  
  bufCount = 0;
  
  buffer[bufCount++] = (double) adef->bestTrav;
  buffer[bufCount++] = (double) tr->treeID;
  buffer[bufCount++] = tr->likelihood;
  buffer[bufCount++] = t;

  for(i = 0; i < tr->NumberOfModels; i++)        
    buffer[bufCount++] = tr->alphas[i];

  for(i = 0; i < tr->NumberOfModels; i++)        
    buffer[bufCount++] = tr->invariants[i];
    
    
  if(adef->boot || adef->rapidBoot)
    {
     if(adef->bootstrapBranchLengths)
       Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH);
     else
       Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES);
    }
  else
    {
      if((adef->model == M_GTRCAT || adef->model == M_PROTCAT) && (adef->useMixedModel == 0))
	Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES);
      else
	Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH);
    }

  tree_ptr = tr->tree_string;

  while(*tree_ptr != ';')    
    buffer[bufCount++] = (double)*tree_ptr++;        
 
  buffer[bufCount++] = (double)(';');
  buffer[bufCount++] = (double)('\n');

  MPI_Send(buffer, bufferSize, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
  
  free(buffer);
}

static void receiveTree(tree *tr, analdef *adef, int workerID, double *t, int tag)
{
  int bufferSize, i, bufCount;
  double *buffer, *buf_ptr;
  char *tree_ptr, content;
  MPI_Status msgStatus; 

  bufferSize = tr->treeStringLength + 4 + tr->NumberOfModels + tr->NumberOfModels;

  buffer = (double *)malloc(sizeof(double) * bufferSize);

  MPI_Recv(buffer, bufferSize, MPI_DOUBLE, workerID, tag, MPI_COMM_WORLD, &msgStatus);
  
  bufCount = 0;
  
  adef->bestTrav = (int)buffer[bufCount++]; 
  tr->treeID     = (int) buffer[bufCount++];
  tr->likelihood = buffer[bufCount++];
  *t = buffer[bufCount++];

  tr->likelihoods[tr->treeID] = tr->likelihood;

  for(i = 0; i < tr->NumberOfModels; i++)    
    tr->alphas[i] = buffer[bufCount++];

  for(i = 0; i < tr->NumberOfModels; i++)    
    tr->invariants[i] = buffer[bufCount++];

  buf_ptr = &buffer[bufCount];
  tree_ptr = tr->tree_string;

  while((content = (char)(buffer[bufCount++])) != ';')
    {      
      *tree_ptr++ = content;
    }
  
  *tree_ptr++ = ';';
  *tree_ptr++ = '\n';
#ifdef DEBUG
  printf("Received tree %s\n", tr->tree_string);
#endif 
  free(buffer);
}




void doBootstrap(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int i, n, dummy;
  double loopTime;
  MPI_Status msgStatus; 

  n = adef->multipleRuns;
          
  if(processID == 0)
    {
      int jobsSent = 0;
      int jobsReceived = n;

      while(jobsReceived > 0)
	{
	  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus);
	  switch(msgStatus.MPI_TAG)
	    {
	    case JOB_REQUEST:
#ifdef DEBUG
	      printf("Master receiving work request from worker %d\n",  msgStatus.MPI_SOURCE);
#endif	      
	      MPI_Recv(&dummy, 1, MPI_INT, msgStatus.MPI_SOURCE, JOB_REQUEST, MPI_COMM_WORLD, &msgStatus);
	       if(jobsSent < n)
		 {
		   MPI_Send(&jobsSent, 1, MPI_INT, msgStatus.MPI_SOURCE, COMPUTE_TREE, MPI_COMM_WORLD);
#ifdef DEBUG
		   printf("Master sending job %d to worker %d\n",  jobsSent, msgStatus.MPI_SOURCE);
#endif
		   jobsSent++;
		 }
	       break;
	    case TREE:
#ifdef DEBUG
	      printf("--------> Master receiving tree from worker %d\n",  msgStatus.MPI_SOURCE);	
#endif
	      receiveTree(tr, adef, msgStatus.MPI_SOURCE, &loopTime, TREE);	     	   	      
	      printBootstrapResult(tr, adef, TRUE);
	      printf("Bootstrap[%d] completed\n", tr->treeID);	 
	      writeInfoFile(adef, tr, loopTime);
	      jobsReceived--;
	      if(jobsSent < n)
		{
		  MPI_Send(&jobsSent, 1, MPI_INT, msgStatus.MPI_SOURCE, COMPUTE_TREE, MPI_COMM_WORLD);
#ifdef DEBUG
		  printf("Master sending job %d to worker %d\n",  jobsSent, msgStatus.MPI_SOURCE);
#endif
		  jobsSent++;
		}
	      break;
	    }
	}
      
       for(i = 1; i < numOfWorkers; i++)
	{
	  MPI_Send(&dummy, 1, MPI_INT, i, FINALIZE, MPI_COMM_WORLD);
#ifdef DEBUG
	  printf("Master sending FINALIZE to worker %d\n",  i);
#endif
	}
       return;
    }
  else
    {
      int treeCounter = 0;

      MPI_Send(&dummy, 1, MPI_INT, 0, JOB_REQUEST, MPI_COMM_WORLD);
#ifdef DEBUG
      printf("Worker %d sending job request to master\n",  processID);
#endif      
       while(1)
	{	
	  MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus); 
	  	 
	  switch(msgStatus.MPI_TAG)
	    {
	    case COMPUTE_TREE: 
	      MPI_Recv(&dummy, 1, MPI_INT, 0, COMPUTE_TREE, MPI_COMM_WORLD, &msgStatus);	      
#ifdef DEBUG
	      printf("Worker %d receiving job %d from master\n",  processID, dummy);
#endif	
	      loopTime = masterTime = gettime();
	      
	      if(adef->multiBoot < 2)
		singleBootstrap(tr, dummy, adef, rdta, cdta);     
	      else
		multipleBootstrap(tr, dummy, adef, rdta, cdta);

	      treeCounter++;
	      loopTime = gettime() - loopTime;
	      sendTree(tr, adef, loopTime, TRUE, TREE);
	      break;
	    case FINALIZE:
	      MPI_Recv(&dummy, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);
#ifdef DEBUG
	      printf("Worker %d receiving FINALIZE %d\n",  processID);
#endif
	      return;
	    }
	}
    }  
}

void doInference(tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta)
{
  int i, n, dummy;
  double loopTime;
  MPI_Status msgStatus; 

  n = adef->multipleRuns;
          
  if(processID == 0)
    {
      int jobsSent = 0;
      int jobsReceived = n;

      while(jobsReceived > 0)
	{
	  MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus);
	  switch(msgStatus.MPI_TAG)
	    {
	    case JOB_REQUEST:
#ifdef DEBUG
	      printf("Master receiving work request from worker %d\n",  msgStatus.MPI_SOURCE);
#endif	      
	      MPI_Recv(&dummy, 1, MPI_INT, msgStatus.MPI_SOURCE, JOB_REQUEST, MPI_COMM_WORLD, &msgStatus);
	      if(jobsSent < n)
		{
		  MPI_Send(&jobsSent, 1, MPI_INT, msgStatus.MPI_SOURCE, COMPUTE_TREE, MPI_COMM_WORLD);
#ifdef DEBUG
		  printf("Master snding job %d to worker %d\n",  jobsSent, msgStatus.MPI_SOURCE);
#endif
		  jobsSent++;
		}
	       break;
	    case TREE:
#ifdef DEBUG
	      printf("--------> Master receiving tree from worker %d\n",  msgStatus.MPI_SOURCE);	
#endif
	      receiveTree(tr, adef, msgStatus.MPI_SOURCE, &loopTime, TREE);	     	   	      	      
	      printf("Inference[%d] completed\n", tr->treeID);	 
	      writeInfoFile(adef, tr, loopTime);
	      jobsReceived--;
	      if(jobsSent < n)
		{
		  MPI_Send(&jobsSent, 1, MPI_INT, msgStatus.MPI_SOURCE, COMPUTE_TREE, MPI_COMM_WORLD);
#ifdef DEBUG
		  printf("Master sending job %d to worker %d\n",  jobsSent, msgStatus.MPI_SOURCE);
#endif
		  jobsSent++;
		}
	      break;
	    }
	}
      
       for(i = 1; i < numOfWorkers; i++)
	{
	  MPI_Send(&dummy, 1, MPI_INT, i, FINALIZE, MPI_COMM_WORLD);
#ifdef DEBUG
	  printf("Master sending FINALIZE to worker %d\n",  i);
#endif
	}
       return;
    }
  else
    {
      int treeCounter = 0;

      MPI_Send(&dummy, 1, MPI_INT, 0, JOB_REQUEST, MPI_COMM_WORLD);
#ifdef DEBUG
      printf("Worker %d sending job request to master\n",  processID);
#endif      
       while(1)
	{	
	  MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus); 
	  	 
	  switch(msgStatus.MPI_TAG)
	    {
	    case COMPUTE_TREE: 
	      MPI_Recv(&dummy, 1, MPI_INT, 0, COMPUTE_TREE, MPI_COMM_WORLD, &msgStatus);	      
#ifdef DEBUG
	      printf("Worker %d receiving job %d from master\n",  processID, dummy);
#endif
	      loopTime =  masterTime = gettime();

	      tr->treeID = dummy;
	      tr->checkPointCounter = 0;
	                    
	      if(treeCounter > 0)	       
		initModel(tr, rdta, cdta, adef); 

	      treeCounter++;

	      getStartingTree(tr, adef);  
     
	      computeBIGRAPID(tr, adef);                     

	      if(adef->model == M_GTRGAMMA || adef->model == M_PROTGAMMA)
		{
		  modOpt(tr, adef);			
		  printLog(tr, adef, TRUE);
		  printResult(tr, adef, TRUE);
		  loopTime = gettime() - loopTime;	 		 
		  freeNodex(tr); 
		}
	      else
		{	  
		  if(adef->useMixedModel)
		    {
		      tr->likelihood = unlikely;

		      catToGamma(tr, adef);

		      initModel(tr, rdta, cdta, adef);	  	  	  
		      modOpt(tr, adef);	
		      printLog(tr, adef, TRUE);
		      printResult(tr, adef, TRUE);
		      loopTime = gettime() - loopTime;			    

		      gammaToParsimony(tr, adef);		      
		    }
		  else
		    {
		      loopTime = gettime() - loopTime;       		      
		      freeNodex(tr);
		    }	
		}
	     
	      sendTree(tr, adef, loopTime, TRUE, TREE);
	      break;
	    case FINALIZE:
	      MPI_Recv(&dummy, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);
#ifdef DEBUG
	      printf("Worker %d receiving FINALIZE %d\n",  processID);
#endif
	      return;
	    }
	}
    }  
}


static void allInOneMaster(tree *tr, analdef *adef)
{  
  MPI_Status msgStatus;
  double loopTime; 
  int workers = numOfWorkers - 1;
  int i, width, dummy, whoHasBestTree = -1;
  int bsCount, bsTotal;
  int mlCount, mlTotal;
  int finished = FALSE;  
  FILE *infoFile;
  double bestLikelihood = unlikely;  

  if(adef->multipleRuns % workers == 0)
    width = adef->multipleRuns / workers;
  else
    width = (adef->multipleRuns / workers) + 1;    		      

  bsCount = 0;
  mlCount = 0;
  bsTotal = width * workers;
  mlTotal = workers;

  free(tr->likelihoods);
  tr->likelihoods = (double *)malloc((width * workers + workers) * sizeof(double));

  if(adef->model == M_PROTGAMMA || adef->model == M_GTRGAMMA || 
     tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
    {     
      printf("\nSwitching from GAMMA to CAT for rapid Bootstrap, final ML search will be conducted under the %s model you specified\n", 
	     (adef->useInvariant)?"GAMMA+P-Invar":"GAMMA");
      infoFile = fopen(infoFileName, "a");
      fprintf(infoFile, 
	      "\nSwitching from GAMMA to CAT for rapid Bootstrap, final ML search will be conducted under the %s model you specified\n", 
	      (adef->useInvariant)?"GAMMA+P-Invar":"GAMMA");
      fclose(infoFile);           
    }

  while(! finished)
    {
      MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus);
      switch(msgStatus.MPI_TAG)
	{
	case BS_TREE:
	  receiveTree(tr, adef, msgStatus.MPI_SOURCE, &loopTime, BS_TREE);
#ifdef DEBUG	 
	  printf("Received BS TREE %d %f\n", bsCount, tr->likelihood);
#endif
	  printBootstrapResult(tr, adef, TRUE);
	  writeInfoFile(adef, tr, loopTime);
	  bsCount++;
	  break;
	case ML_TREE: 
	  receiveTree(tr, adef, msgStatus.MPI_SOURCE, &loopTime, ML_TREE);
	  if(tr->likelihood > bestLikelihood)
	    {
	      bestLikelihood = tr->likelihood;
	      whoHasBestTree = msgStatus.MPI_SOURCE;
	    }
#ifdef DEBUG
	  printf("Received ML TREE %d %f ID %d\n", mlCount, tr->likelihood, tr->treeID);
#endif
	  mlCount++;
	  break;
	default:
	  assert(0);
	}
      if(adef->allInOne)	
	finished = (bsCount == bsTotal && mlCount == mlTotal);
      else
	finished = (bsCount == bsTotal);

      if(adef->allInOne && bsCount == bsTotal && mlCount == 0)
	{
	  double t = gettime() - masterTime;
	  infoFile = fopen(infoFileName, "a");

	  printf("\n\n");
	  fprintf(infoFile, "\n\n");

	  if(adef->bootStopping)
	    {
	      assert(0);	     
	    }

	  printf("Overall Time for %d Rapid Bootstraps %f\n", bsCount, t);
	  fprintf(infoFile, "Overall Time for %d Rapid Bootstraps %f\n", bsCount, t);

	  printf("Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bsCount)));                
	  fprintf(infoFile, "Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bsCount)));	           
     
	  fclose(infoFile);    
	}

      if(finished)
	{
	  if(adef->allInOne)
	    {
	      double overallTime;
	      char bestTreeFileName[1024];

	      strcpy(bestTreeFileName, workdir);
	      strcat(bestTreeFileName, "RAxML_bestTree.");
	      strcat(bestTreeFileName,         run_id);

	      assert(whoHasBestTree > 0 && whoHasBestTree < numOfWorkers);
#ifdef DEBUG
	      printf("worker %d xas mpest tri with %f \n", whoHasBestTree, bestLikelihood);
#endif	      
	      MPI_Send(&dummy, 1, MPI_INT, whoHasBestTree, PRINT_TREE, MPI_COMM_WORLD);
	      MPI_Recv(&dummy, 1, MPI_INT, whoHasBestTree, I_PRINTED_IT,  MPI_COMM_WORLD, &msgStatus);

	      overallTime = gettime() - masterTime;
	      printf("Program execution info written to %s\n", infoFileName);
	      printf("All %d bootstrapped trees written to: %s\n\n", bsCount, bootstrapFileName);
	      printf("Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
	      if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
		printf("Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n", 
		       bestTreeFileName,  bestTreeFileName, tr->NumberOfModels - 1);    
	      printf("Best-scoring ML tree with support values written to: %s\n\n", bipartitionsFileName);
	      printf("Overall execution time for full ML analysis: %f secs or %f hours or %f days\n\n", 
		     overallTime, overallTime/3600.0, overallTime/86400.0);
  
	      infoFile = fopen(infoFileName, "a");
	      fprintf(infoFile, "All %d bootstrapped trees written to: %s\n\n", bsCount, bootstrapFileName);
	      fprintf(infoFile, "Best-scoring ML tree written to: %s\n\n", bestTreeFileName);
	      if(adef->perGeneBranchLengths && tr->NumberOfModels > 1)    
		fprintf(infoFile, "Per-Partition branch lengths of best-scoring ML tree written to %s.PARTITION.0 to  %s.PARTITION.%d\n\n"
			, bestTreeFileName,  bestTreeFileName, tr->NumberOfModels - 1);    
	      fprintf(infoFile, "Best-scoring ML tree with support values written to: %s\n\n", 
		      bipartitionsFileName);
	      fprintf(infoFile, "Overall execution time for full ML analysis: %f secs or %f hours or %f days\n\n", 
		      overallTime, overallTime/3600.0, overallTime/86400.0);
	      fclose(infoFile);   
	    }
	  else
	    {
	      double t = gettime() - masterTime;
	      infoFile = fopen(infoFileName, "a");

	      printf("\n\n");
	      fprintf(infoFile, "\n\n");

	      if(adef->bootStopping)
		{
		  assert(0);
		}

	      printf("Overall Time for %d Rapid Bootstraps %f\n", bsCount, t);
	      fprintf(infoFile, "Overall Time for %d Rapid Bootstraps %f\n", bsCount, t);

	      printf("Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bsCount)));
	      fprintf(infoFile, "Average Time per Rapid Bootstrap %f\n", (double)(t/((double)bsCount)));

	      printf("All %d bootstrapped trees written to: %s\n", bsCount, bootstrapFileName);
	      fprintf(infoFile, "All %d bootstrapped trees written to: %s\n", bsCount, bootstrapFileName);           
     
	      fclose(infoFile);	      
	    }

	  for(i = 1; i < numOfWorkers; i++)	    
	    MPI_Send(&dummy, 1, MPI_INT, i, FINALIZE, MPI_COMM_WORLD);	 
	}
    }
  
  MPI_Finalize();
  exit(0);
}

static void allInOneWorker(tree *tr, analdef *adef)
{
  int dummy, NumberOfLocalTrees;
  MPI_Status msgStatus; 

  if(adef->multipleRuns % (numOfWorkers - 1)  == 0)
    NumberOfLocalTrees = adef->multipleRuns / (numOfWorkers - 1);
  else
    NumberOfLocalTrees = (adef->multipleRuns / (numOfWorkers - 1)) + 1;

#ifdef DEBUG
  printf("Worker %d %d\n", processID, NumberOfLocalTrees);
#endif
  /* re-initialize adfe->rapidBoot, otherwise the workers will be doing the exact same replicates */

  adef->rapidBoot = (long)gettimeSrand();

  /* the one below is kind of an ugly fix, but who cares */

  tr->treeID = NumberOfLocalTrees * (processID - 1);

  	  
  {
    int i, n, sites, bestIndex, model, bootstrapsPerformed;
    double loopTime = 0.0; 
    int      *originalRateCategories;
    int      *originalInvariant;
    int      slowSearches, fastEvery = 5;
    topolRELL_LIST *rl;  
    double bestLH, mlTime;  
    long radiusSeed = adef->rapidBoot;
    FILE *infoFile, *f;
    char bestTreeFileName[1024];  
    modelParams *catParams   = (modelParams *)malloc(sizeof(modelParams));
    modelParams *gammaParams = (modelParams *)malloc(sizeof(modelParams));

    allocParams(catParams,   tr);
    allocParams(gammaParams, tr);
   
    n = NumberOfLocalTrees;   

    rl = (topolRELL_LIST *)malloc(sizeof(topolRELL_LIST));
    initTL(rl, tr, n);
     
    originalRateCategories = (int*)malloc(tr->cdta->endsite * sizeof(int));
    originalInvariant      = (int*)malloc(tr->cdta->endsite * sizeof(int));

    sites = tr->cdta->endsite;             

    if(adef->model == M_PROTGAMMA || adef->model == M_GTRGAMMA || 
       tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
      {
	assert(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I);
	assert(adef->model == M_GTRGAMMA || adef->model == M_PROTGAMMA);

	tr->rateHetModel = CAT;
	
	switch(adef->model)
	  {
	  case M_GTRGAMMA:
	    adef->model = M_GTRCAT;	     
	    break;
	  case M_PROTGAMMA:
	    adef->model = M_PROTCAT;
	    break;
	  default:
	    assert(0);
	  }
	
	initModel(tr, tr->rdta, tr->cdta, adef);                
      }
 
    for(i = 0; i < n; i++)
      {                
	tr->checkPointCounter = 0;
	
	loopTime = gettime();      
	
	if(i % 10 == 0)
	  {
	    if(i > 0)
	      {
		freeNodex(tr);
		reductionCleanup(tr, adef, originalRateCategories, originalInvariant);
	      }
	    
	    makeParsimonyTree(tr, adef);
	    allocNodex(tr, adef); 
	    tr->likelihood = unlikely;
	    if(i == 0)
	      {
		double t;
		
		onlyInitrav(tr, tr->start);
		
		treeEvaluate(tr, 1);	     	     
		
		t = gettime();    
		
		quickOpt(tr, adef);	     		
		memcpy(originalRateCategories, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);
		memcpy(originalInvariant,      tr->invariant,          sizeof(int) * tr->cdta->endsite);

		if(adef->bootstrapBranchLengths)
		  {
		    storeParams(catParams, tr);
		    assert(tr->cdta->endsite == tr->originalCrunchedLength);
		    catToGamma(tr, adef);
		    modOpt(tr, adef);
		    storeParams(gammaParams, tr);
		    gammaToCat(tr, adef);
		    loadParams(catParams, tr);		  
		  }
	      }	  	  
	  }
	
	computeNextReplicate(tr, adef, originalRateCategories, originalInvariant); 
	resetBranches(tr);
	
	evaluateGenericInitrav(tr, tr->start);    
	
	treeEvaluate(tr, 1);    	             
	
	computeBOOTRAPID(tr, adef, &radiusSeed);                        	  
	saveTL(rl, tr, i);
      
	if(adef->bootstrapBranchLengths)
	  {
	    double lh = tr->likelihood;
	    int    endsite;
	    
	    loadParams(gammaParams, tr);
	    
	    
	    endsite = tr->cdta->endsite;
	    tr->cdta->endsite = tr->originalCrunchedLength;
	    catToGamma(tr, adef);
	    tr->cdta->endsite = endsite;
	    
	    resetBranches(tr);
	    treeEvaluate(tr, 2.0);
	    
	    endsite = tr->cdta->endsite;
	    tr->cdta->endsite = tr->originalCrunchedLength;
	    gammaToCat(tr, adef);
	    tr->cdta->endsite = endsite;	 	    
	    
	    loadParams(catParams, tr);
	    	    
	    tr->likelihood = lh;
	  }
      

	loopTime = gettime() - loopTime;
	sendTree(tr, adef, loopTime, TRUE, BS_TREE);	
     
	if(adef->bootStopping)
	  {
	    assert(0);
	    /*bootStopIt = bootStop(tr, b, i, &pearsonAverage);*/
	  }
	tr->treeID = tr->treeID + 1;
      }  
 
    bootstrapsPerformed = i;
        
    if(!adef->allInOne)
      {       	
	MPI_Recv(&dummy, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);
	MPI_Finalize();	 
	exit(0);
      }
    
    freeParams(catParams);
    free(catParams);

    freeParams(gammaParams);
    free(gammaParams);

    tr->treeID = NumberOfLocalTrees * (numOfWorkers - 1) + (processID - 1);
   
    mlTime = gettime();
  
#ifdef DEBUG
    printf("\nWorker %d Starting ML Search ...\n\n", processID);   
#endif

    reductionCleanup(tr, adef, originalRateCategories, originalInvariant);    
 
    catToGamma(tr, adef);     	   
    
    restoreTL(rl, tr, 0);

    resetBranches(tr);
    
    evaluateGenericInitrav(tr, tr->start);   
    
    modOpt(tr, adef);
    
    evaluateGenericInitrav(tr, tr->start);
    
    categorizeGeneric(tr, tr->start);        
    
    gammaToCat(tr, adef);      
         
    fastEvery = 5;    

    for(i = 0; i < bootstrapsPerformed; i++)
      {            
	rl->t[i]->likelihood = unlikely;
	
	if(i % fastEvery == 0)
	  {	 
	    restoreTL(rl, tr, i); 	 	    	   
	    
	    resetBranches(tr);

	    evaluateGenericInitrav(tr, tr->start);
	    
	    treeEvaluate(tr, 1); 		 
	    
	    optimizeRAPID(tr, adef);	  		
	    saveTL(rl, tr, i);      
	  }    
      }     

#ifdef DEBUG  
    printf("Worker %d Fast ML optimization finished\n\n", processID);
#endif
   
    qsort(&(rl->t[0]), bootstrapsPerformed, sizeof(topolRELL*), compareTopolRell);
     
    catToGamma(tr, adef);  
    
    restoreTL(rl, tr, 0);

    resetBranches(tr);
    
    evaluateGenericInitrav(tr, tr->start); 
    
    modOpt(tr, adef);
    
    evaluateGenericInitrav(tr, tr->start);
    
    categorizeGeneric(tr, tr->start);   
    
    gammaToCat(tr, adef);    
  
    slowSearches = bootstrapsPerformed / 5;
    if(bootstrapsPerformed % 5 != 0)
      slowSearches++;

    slowSearches  = MIN(slowSearches, 10); 

    for(i = 0; i < slowSearches; i++)
      {           
	restoreTL(rl, tr, i);     
	rl->t[i]->likelihood = unlikely;
	
	evaluateGenericInitrav(tr, tr->start);
	
	treeEvaluate(tr, 1.0);   
	thoroughOptimization(tr, adef, rl, i);             
      }
 
#ifdef DEBUG
    printf("Worker %d Slow ML optimization finished\n\n", processID);
#endif     
    catToGamma(tr, adef);    
    
    bestIndex = -1;
    bestLH = unlikely;
    
    for(i = 0; i < slowSearches; i++)
      {      
	restoreTL(rl, tr, i);
	resetBranches(tr);

	evaluateGenericInitrav(tr, tr->start);
	
	treeEvaluate(tr, 2);
#ifdef DEBUG	
	printf("Worker %d Slow ML Search %d Likelihood: %f\n", processID, i, tr->likelihood);
#endif	
	if(tr->likelihood > bestLH)
	  {
	    bestLH = tr->likelihood;
	    bestIndex = i;
	  }
      }
    
    restoreTL(rl, tr, bestIndex);
    resetBranches(tr);

    evaluateGenericInitrav(tr, tr->start);
    
    treeEvaluate(tr, 2); 
    
    Thorough = 1;
    tr->doCutoff = FALSE;  
    
    treeOptimizeThorough(tr, 1, 10);
    modOpt(tr, adef);

    sendTree(tr, adef, loopTime, TRUE, ML_TREE);
   
    while(1)
      {
	MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &msgStatus);
	switch(msgStatus.MPI_TAG)
	  {
	  case FINALIZE:	
	    MPI_Recv(&dummy, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);
	    MPI_Finalize();  
	    exit(0);
	    break;
	  case PRINT_TREE:
	    MPI_Recv(&dummy, 1, MPI_INT, 0, PRINT_TREE, MPI_COMM_WORLD, &msgStatus);
#ifdef DEBUG
	    printf("Ox malaka eimai o %d kai prepi na printaro to dentro %f re pousti \n", processID, tr->likelihood);
#endif
	    infoFile = fopen(infoFileName, "a"); 
	    
	    printf("\nFinal ML Optimization Likelihood: %f\n", tr->likelihood);
	    fprintf(infoFile, "\nFinal ML Optimization Likelihood: %f\n", tr->likelihood);
	    
	    printf("\nModel Information:\n\n");
	    fprintf(infoFile, "\nModel Information:\n\n");
	    
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
		
		fprintf(infoFile, "Model Parameters of Partition %d, Name: %s, Type of Data: %s\n", 
			model, tr->partitionData[model].partitionName, typeOfData);
		fprintf(infoFile, "alpha: %f\n", tr->alphas[model]);
		
		printf("Model Parameters of Partition %d, Name: %s, Type of Data: %s\n", 
		       model, tr->partitionData[model].partitionName, typeOfData);
		printf("alpha: %f\n", tr->alphas[model]);
		
		if(adef->useInvariant)
		  {
		    fprintf(infoFile, "invar: %f\n", tr->invariants[model]);    
		    printf("invar: %f\n", tr->invariants[model]);    
		  }
		
		if(adef->perGeneBranchLengths)
		  tl = treeLength(tr, model);
		else
		  tl = treeLength(tr, 0);
		
		fprintf(infoFile, "Tree-Length: %f\n", tl);    
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
			{
			  fprintf(infoFile, "rate %s: %f\n", names[k], tr->initialRates_DNA[model * DNA_RATES + k]);			    
			  printf("rate %s: %f\n", names[k], tr->initialRates_DNA[model * DNA_RATES + k]);
			}
		      
		      fprintf(infoFile, "rate %s: %f\n", names[5], 1.0);
		      printf("rate %s: %f\n", names[5], 1.0);
		    }      
		    break;
		  default:
		    assert(0);
		  }
		
		fprintf(infoFile, "\n");
		printf("\n");
	      }		    		  
	    
	    fclose(infoFile); 
	    
	    strcpy(bestTreeFileName, workdir); 
	    strcat(bestTreeFileName, "RAxML_bestTree.");
	    strcat(bestTreeFileName,         run_id);
	    
	    Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, TRUE, adef, SUMMARIZE_LH);
	    f = fopen(bestTreeFileName, "w");
	    fprintf(f, "%s", tr->tree_string);
	    fclose(f);

	    if(adef->perGeneBranchLengths)
	      printTreePerGene(tr, adef, bestTreeFileName, "w");
	    
	    infoFile = fopen(infoFileName, "a"); 
	    
	    
	    printf("\nDrawing Bootstrap Support Values on best-scoring ML tree ...\n\n");
	    fprintf(infoFile, "Drawing Bootstrap Support Values on best-scoring ML tree ...\n\n");  
	    
	    fclose(infoFile);
	    
	    freeTL(rl, tr);   
	    free(rl);       
	    
	    calcBipartitions(tr, adef, bestTreeFileName, bootstrapFileName);  
	     
	    

	    MPI_Send(&dummy, 1, MPI_INT, 0, I_PRINTED_IT, MPI_COMM_WORLD);
	    break;
	  default:
	    assert(0);
	  }
      }
  }
}


void doAllInOne(tree *tr, analdef *adef)
{
  if(processID == 0)
    allInOneMaster(tr, adef);
  else
    allInOneWorker(tr, adef);
}


#endif
