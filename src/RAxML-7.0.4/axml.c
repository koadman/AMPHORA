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

#ifdef WIN32
#include <direct.h>
#endif

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

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef _USE_OMP
#include <omp.h>
#endif

#ifdef _USE_PTHREADS
#include <pthread.h>
#endif


#include "axml.h"
#include "globalVariables.h"


/***************** UTILITY FUNCTIONS **************************/





double gettime(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 60.0*localtm.tm_min + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec * 0.000001;
#endif
}

int gettimeSrand(void)
{
#ifdef WIN32
  time_t tp;
  struct tm localtm;
  tp = time(NULL);
  localtm = *localtime(&tp);
  return 24*60*60*localtm.tm_yday + 60*60*localtm.tm_hour + 60*localtm.tm_min  + localtm.tm_sec;
#else
  struct timeval ttime;
  gettimeofday(&ttime , NULL);
  return ttime.tv_sec + ttime.tv_usec;
#endif
}

double randum (long  *seed)
    /* random number generator, modified to use 12 bit chunks */
  { /* randum */
    long  sum, mult0, mult1, seed0, seed1, seed2, newseed0, newseed1, newseed2;
    double res;

    mult0 = 1549;
    seed0 = *seed & 4095;
    sum  = mult0 * seed0;
    newseed0 = sum & 4095;
    sum >>= 12;
    seed1 = (*seed >> 12) & 4095;
    mult1 =  406;
    sum += mult0 * seed1 + mult1 * seed0;
    newseed1 = sum & 4095;
    sum >>= 12;
    seed2 = (*seed >> 24) & 255;
    sum += mult0 * seed2 + mult1 * seed1;
    newseed2 = sum & 255;

    *seed = newseed2 << 24 | newseed1 << 12 | newseed0;
    res = 0.00390625 * (newseed2 + 0.000244140625 * (newseed1 + 0.000244140625 * newseed0));     

    return res;   
  } /* randum */

static int filexists(char *filename)
{
  FILE *fp;
  int res;
  fp = fopen(filename,"r");
  
  if(fp) 
    {
      res = 1;
      fclose(fp);
    }
  else 
    res = 0;
       
  return res;
} 

/********************* END UTILITY FUNCTIONS ********************/


/******************************some functions for the likelihood computation ****************************/

#ifdef WIN32
boolean isTip(int number, int maxTips)
#else
inline boolean isTip(int number, int maxTips)
#endif
{  
  assert(number > 0);

  if(number <= maxTips)
    return TRUE;
  else
    return FALSE;
}


#ifdef WIN32
double *getLikelihoodArray(int number, int mxtips, double **xVector)
#else
inline double *getLikelihoodArray(int number, int mxtips, double **xVector)
#endif
{
  return (xVector[number - mxtips - 1]);
}

#ifdef WIN32
int *getScalingArray(int number, int endsite, int mxtips, int *scalingArray)
#else
inline int *getScalingArray(int number, int endsite, int mxtips, int *scalingArray)
#endif
{  
  return &(scalingArray[endsite * (number - mxtips - 1)]);
}


#ifdef _MULTI_GENE
void getxsnode (nodeptr p, int model)  
{  
  assert(p->xs[model] || p->next->xs[model] || p->next->next->xs[model]);
  assert(p->xs[model] + p->next->xs[model] + p->next->next->xs[model] == 1);
  
  assert(p == p->next->next->next);

  p->xs[model] = 1;
  
  if(p->next->xs[model])
    {      
      p->next->xs[model] = 0;
      return;
    }
  else
    {
      p->next->next->xs[model] = 0;
      return;
    }  

  assert(0);
}

#endif

void getxnode (nodeptr p)  
{ 
  nodeptr  s;
 
  if ((s = p->next)->x || (s = s->next)->x) 
    {
      p->x = s->x;
      s->x = 0;
    }     
  
  assert(p->x);
}





void hookup (nodeptr p, nodeptr q, double *z, int numBranches)
{
  int i;

  p->back = q;
  q->back = p;
    
  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = z[i];
}

void hookupDefault (nodeptr p, nodeptr q, int numBranches)
{
  int i;
  
  p->back = q;
  q->back = p;
    
  for(i = 0; i < numBranches; i++)
    p->z[i] = q->z[i] = defaultz;  
}


/******************************some functions for the likelihood computation ****************************/




/***********************reading and initializing input ******************/

static void getnums (rawdata *rdta)
{    
  if (fscanf(INFILE, "%d %d", & rdta->numsp, & rdta->sites) != 2) 
    {
      if(processID == 0)
	printf("ERROR: Problem reading number of species and sites\n");
      errorExit(-1);
    }
    
  if (rdta->numsp < 4) 
    {
      if(processID == 0)
	printf("TOO FEW SPECIES\n");
      errorExit(-1);
    }

  if (rdta->sites < 1) 
    {
      if(processID == 0)
	printf("TOO FEW SITES\n");
      errorExit(-1);
    }
  
  return;
}


static boolean digitchar (int ch) 
{
  return (ch >= '0' && ch <= '9'); 
}


boolean whitechar (int ch)
{ 
  return (ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r'); /* PC-LINEBREAK*/
}


static void uppercase (int *chptr)
{
  int  ch;
  
  ch = *chptr;
  if ((ch >= 'a' && ch <= 'i') || (ch >= 'j' && ch <= 'r')
      || (ch >= 's' && ch <= 'z'))
    *chptr = ch + 'A' - 'a';
} 




static void getyspace (rawdata *rdta)
{
  long   size;
  int    i;
  char *y0;

  if (! (rdta->y = (char **) malloc((rdta->numsp + 1) * sizeof(char *)))) 
    {
      printf("ERROR: Unable to obtain space for data array pointers\n");
      exit(-1);
    }

  size = 4 * (rdta->sites / 4 + 1);

  if (! (y0 = (char *) malloc((rdta->numsp + 1) * size * sizeof(char)))) 
    {
      printf("ERROR: Unable to obtain space for data array\n");
      exit(-1);
    }

  rdta->y0 = y0;

  for (i = 0; i <= rdta->numsp; i++) 
    {
      rdta->y[i] = y0;
      y0 += size;
    }

  return;
}

static boolean setupTree (tree *tr, analdef *adef)
{
  nodeptr  p0, p, q;
  int      
    i, 
    j, 
    tips, 
    inter;

  tr->expArray        = (int*)NULL;
  tr->likelihoodArray = (double*)NULL;
  tr->sumBuffer       = (double*)NULL;

  tr->bigCutoff = FALSE;
  tr->mixedData = FALSE;

  tr->ib = (insertionBranch *)NULL;
  tr->ip = (insertionPoints *)NULL;
  tr->numberOfTipsForInsertion = 0;

  for(i = 0; i < NUM_BRANCHES; i++)
    tr->partitionContributions[i] = -1.0;    
  

  if(!adef->useMultipleModel)
    tr->NumberOfModels = 1;

  if(adef->grouping)
    tr->grouped = TRUE;
  else
    tr->grouped = FALSE;
  
  if(adef->constraint)
    tr->constrained = TRUE;
  else
    tr->constrained = FALSE;
  
  tr->treeID = 0;
  
  tips  = tr->mxtips;
  inter = tr->mxtips - 1;
  
  tr->tipVectorDNA    = (double *)malloc(tr->NumberOfModels * 64 * sizeof(double));  
  tr->tipVectorAA     = (double *)malloc(tr->NumberOfModels * 460 * sizeof(double)); 

  tr->EV_DNA          = (double *)malloc(tr->NumberOfModels * 16 * sizeof(double));
  tr->EV_AA           = (double *)malloc(tr->NumberOfModels * 400 * sizeof(double));

  tr->EI_DNA          = (double *)malloc(tr->NumberOfModels * 12 * sizeof(double));
  tr->EI_AA           = (double *)malloc(tr->NumberOfModels * 380 * sizeof(double));  

  tr->EIGN_DNA        = (double *)malloc(tr->NumberOfModels * 3 * sizeof(double));
  tr->EIGN_AA         = (double *)malloc(tr->NumberOfModels * 19  * sizeof(double));
  
  tr->frequencies_DNA = (double *)malloc(tr->NumberOfModels * 4 * sizeof(double));
  tr->frequencies_AA  = (double *)malloc(tr->NumberOfModels * 20  * sizeof(double));
  

  tr->initialRates_DNA = (double *)malloc(tr->NumberOfModels * 5 * sizeof(double));
  tr->initialRates_AA = (double *)malloc(tr->NumberOfModels * 190 * sizeof(double));   
  
  tr->xVector      = (double **)malloc(tr->mxtips * sizeof(double *));
  tr->yVector      = (char **)  malloc((tr->mxtips + 1) * sizeof(char *));

  tr->gammaRates   = (double *)malloc(tr->NumberOfModels * 4 * sizeof(double));
  tr->alphas       = (double *)malloc(tr->NumberOfModels * sizeof(double));
  tr->invariants   = (double *)malloc(tr->NumberOfModels * sizeof(double));
  tr->fracchanges  = (double *)malloc(tr->NumberOfModels * sizeof(double));
  tr->likelihoods  = (double *)malloc(adef->multipleRuns * sizeof(double)); 
 
  /* 
     TODO: br-lens values, support values, what happens when phylo classification algo is used
     -> worst case analysis?

     TODO fix for two algos below 
     assert(adef->mode != MEHRING_ALGO && adef->rapidML_Addition != TRUE);
  */

  tr->treeStringLength = tr->mxtips * (nmlngth+128) + 256 + tr->mxtips * 2;
  
  tr->tree_string  = (char*)malloc(tr->treeStringLength * sizeof(char));

  /*TODO, must that be so long ?*/  
 

  tr->td[0].count = 0;
  tr->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * tr->mxtips);

  for(i = 0; i < tr->NumberOfModels; i++)
    tr->fracchanges[i] = -1.0;
  tr->fracchange = -1.0;


#ifdef _MULTI_GENE
  tr->doMulti = 0;
  {
    int k;
    
    for(k = 0; k < tr->numBranches; k++)
      {
	tr->td[k].count = 0;
	tr->td[k].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * tr->mxtips);
      }
  }
#endif


  tr->constraintVector = (int *)malloc((2 * tr->mxtips) * sizeof(int));

  tr->nameList = (char **)malloc(sizeof(char *) * (tips + 1));    
             
  if (!(p0 = (nodeptr) malloc((tips + 3*inter) * sizeof(node)))) 
    {
      printf("ERROR: Unable to obtain sufficient tree memory\n");
      return  FALSE;
    }

  if (!(tr->nodep = (nodeptr *) malloc((2*tr->mxtips) * sizeof(nodeptr)))) 
    {
      printf("ERROR: Unable to obtain sufficient tree memory, too\n");
      return  FALSE;
    }
    
  tr->nodep[0] = (node *) NULL;    /* Use as 1-based array */

  for (i = 1; i <= tips; i++) 
    {
      p = p0++;
      p->x      =  0;    
      p->number =  i;
      p->next   =  p;
      p->back   = (node *) NULL;          
#ifdef  _MULTI_GENE
      {
	int k;
	
	for(k = 0; k < tr->numBranches; k++)
	  {
	    p->xs[k]    = 0;
	    p->backs[k] = (nodeptr)NULL;
	  }
      }
#endif
      tr->nodep[i] = p;
    }

  for (i = tips + 1; i <= tips + inter; i++) 
    {
      q = (node *) NULL;
      for (j = 1; j <= 3; j++) 
	{
	  p = p0++;
	  p->x      =  0;	  
	  p->number = i;
	  p->next   = q;

	  p->back   = (node *) NULL;
#ifdef  _MULTI_GENE
	  {
	    int k;
		  
	    for(k = 0; k < tr->numBranches; k++)
	      {
		p->xs[k]    = 0;
		p->backs[k] = (nodeptr)NULL;
	      }
	  }
#endif
	  q = p;
	}
      p->next->next->next = p;
      tr->nodep[i] = p;
    }

  tr->likelihood  = unlikely;
  tr->start       = (node *) NULL;
  tr->ntips       = 0;
  tr->nextnode    = 0;  
  tr->smoothed    = FALSE;
  
  return TRUE;
} 


static void checkTaxonName(char *buffer, int len)
{
  int i;

  for(i = 0; i < len - 1; i++)
    {
      boolean valid;

      switch(buffer[i])
	{
	case '\0':  
	case '\t':  
	case '\n':  
	case '\r': 
	case ' ':
	case ':':  
	case ',':   
	case '(':   
	case ')':  
	case ';':
	case '[':
	case ']':
	  valid = FALSE;
	  break;
	default:
	  valid = TRUE;
	}
     
      if(!valid)
	{
	  printf("ERROR: Taxon Name \"%s\" is invalid at position %d, it contains illegal character %c\n", buffer, i, buffer[i]);
	  printf("Illegal characters in taxon-names are: tabulators, carriage returns, spaces, \":\", \",\", \")\", \"(\", \";\", \"]\", \"[\"\n");
	  printf("Exiting\n");
	  exit(-1);
	}

    }
  assert(buffer[len - 1] == '\0');
}

static boolean getdata(analdef *adef, rawdata *rdta, tree *tr)
{
  int   i, j, basesread, basesnew, ch, my_i, meaning;
  int   meaningAA[256], meaningDNA[256];
  boolean  allread, firstpass;
  char buffer[300]; 
  int len;
  unsigned long total = 0;
  unsigned long gaps  = 0;
  int gapValueAA, gapValueDNA;
     
  for (i = 0; i <= 255; i++) 
    {
      meaningAA[i] = -1;
      meaningDNA[i] = -1;
    }

  /* AA data */

  meaningAA['A'] =  0;  /* alanine */
  meaningAA['R'] =  1;  /* arginine */
  meaningAA['N'] =  2;  /*  asparagine*/
  meaningAA['D'] =  3;  /* aspartic */
  meaningAA['C'] =  4;  /* cysteine */
  meaningAA['Q'] =  5;  /* glutamine */
  meaningAA['E'] =  6;  /* glutamic */
  meaningAA['G'] =  7;  /* glycine */
  meaningAA['H'] =  8;  /* histidine */
  meaningAA['I'] =  9;  /* isoleucine */
  meaningAA['L'] =  10; /* leucine */
  meaningAA['K'] =  11; /* lysine */
  meaningAA['M'] =  12; /* methionine */
  meaningAA['F'] =  13; /* phenylalanine */
  meaningAA['P'] =  14; /* proline */
  meaningAA['S'] =  15; /* serine */
  meaningAA['T'] =  16; /* threonine */
  meaningAA['W'] =  17; /* tryptophan */
  meaningAA['Y'] =  18; /* tyrosine */
  meaningAA['V'] =  19; /* valine */
  meaningAA['B'] =  20;/* asparagine, aspartic 2 and 3*/
  meaningAA['Z'] =  21;/*21 glutamine glutamic 5 and 6*/
  meaningAA['X'] =  meaningAA['?'] = meaningAA['*'] = meaningAA['-'] = 22; /* all = 1.0 */
  gapValueAA = 22;

  /* DNA data */

  meaningDNA['A'] =  1;
  meaningDNA['B'] = 14;
  meaningDNA['C'] =  2;
  meaningDNA['D'] = 13;
  meaningDNA['G'] =  4;
  meaningDNA['H'] = 11;
  meaningDNA['K'] = 12;
  meaningDNA['M'] =  3;
  meaningDNA['N'] = 15;
  meaningDNA['O'] = 15;
  meaningDNA['R'] =  5;
  meaningDNA['S'] =  6;
  meaningDNA['T'] =  8;
  meaningDNA['U'] =  8;
  meaningDNA['V'] =  7;
  meaningDNA['W'] =  9;
  meaningDNA['X'] = 15;
  meaningDNA['Y'] = 10;     
  meaningDNA['-'] = 15;	
  meaningDNA['?'] = 15;
  gapValueDNA = 15;  

  /*******************************************************************/
  
  basesread = basesnew = 0;

  allread = FALSE;
  firstpass = TRUE;
  ch = ' ';

  while (! allread) 
    {
      for (i = 1; i <= tr->mxtips; i++) 
	{   	  
	  if (firstpass) 
	    {                      	       
	      ch = getc(INFILE);
	      while(ch == ' ' || ch == '\n' || ch == '\t' || ch == '\r') /* PC-LINEBREAK*/
		{
		  ch = getc(INFILE);		  
		}	      
	      my_i = 0;	      

	      do 
		{
		  buffer[my_i] = ch;		  
		  ch = getc(INFILE);		   
		  my_i++;
		  if(my_i >= nmlngth)
		    {
		      if(processID == 0)
			{
			  printf("Taxon Name to long at taxon %d, adapt constant nmlngth in\n", i);
			  printf("axml.h, current setting %d\n", nmlngth);
			}
		      errorExit(-1);
		    }		 
		}
	      while(ch !=  ' ' && ch != '\n' && ch != '\t' && ch != '\r');
			      
	      buffer[my_i] = '\0';	      
	      len = strlen(buffer) + 1;	
	      checkTaxonName(buffer, len);
	      tr->nameList[i] = (char *)malloc(sizeof(char) * len);	      
	      strcpy(tr->nameList[i], buffer);				
	    }

	  j = basesread;
	  while ((j < rdta->sites) && ((ch = getc(INFILE)) != EOF) && (ch != '\n') && (ch != '\r')) /* PC-LINEBREAK*/
	    {
	      uppercase(& ch);

	      assert(tr->dataVector[j + 1] != -1);

	      switch(tr->dataVector[j + 1])
		{
		case DNA_DATA:
		  meaning = meaningDNA[ch];
		  break;
		case AA_DATA:
		  meaning = meaningAA[ch];
		  break;
		default:
		  assert(0);
		}

	      if ((meaning != -1) || ch == '.') 
		{
		  j++;
		  if (ch == '.') 
		    {
		      if (i != 1) 
			ch = rdta->y[1][j];
		      else 
			{
			  printf("ERROR: Dot (.) found at site %d of sequence 1\n", j + 1);
			  return  FALSE;
			}
		    }
		  rdta->y[i][j] = ch;
		}
	      else 
		{
		  if(whitechar(ch) || digitchar(ch)) ;
		  else 
		    {
		      printf("ERROR: Bad base (%c) at site %d of sequence %d\n",
			     ch, j + 1, i);
		      return  FALSE;
		    }
		}
	    }

	  if (ch == EOF) 
	    {
	      printf("ERROR: End-of-file at site %d of sequence %d\n", j + 1, i);
	      return  FALSE;
	    }

	  if (! firstpass && (j == basesread)) 
	    i--; 
	  else 
	    {
	      if (i == 1) 
		basesnew = j;
	      else 
		if (j != basesnew) 
		  {
		    printf("ERROR: Sequences out of alignment\n");		    
		    printf("%d (instead of %d) residues read in sequence %d %s\n",
			   j - basesread, basesnew - basesread, i, tr->nameList[i]);
		    return  FALSE;
		  }
	    }
	  while (ch != '\n' && ch != EOF && ch != '\r') ch = getc(INFILE);  /* flush line *//* PC-LINEBREAK*/
	}

      firstpass = FALSE;
      basesread = basesnew;
      allread = (basesread >= rdta->sites);
    }
     
  for(j = 1; j <= tr->mxtips; j++)    
    for(i = 1; i <= rdta->sites; i++) 	
      {
	assert(tr->dataVector[i] != -1);
	
	switch(tr->dataVector[i])
	  {
	  case DNA_DATA:
	    meaning = meaningDNA[rdta->y[j][i]];
	    if(meaning == gapValueDNA)
	      gaps++;
	    break;
	  case AA_DATA:
	    meaning = meaningAA[rdta->y[j][i]];
	    if(meaning == gapValueAA)
	      gaps++;
	    break;
	  default:
	    assert(0);
	  }     
	
	total++;
	rdta->y[j][i] = meaning;	    
      }
  
  adef->gapyness = (double)gaps / (double)total;

  return  TRUE;
}



static void inputweights (rawdata *rdta)    
{
  int i, w, fres;
  FILE *weightFile;
  int *wv = (int *)malloc(sizeof(int) *  rdta->sites + 1);
    
  weightFile = fopen(weightFileName, "r");
  if (!weightFile)
    {
      if(processID == 0)
	printf( "Could not open weight file: %s\n", weightFileName);
      errorExit(-1);
    }
  
  i = 1;
  
  while((fres = fscanf(weightFile,"%d", &w)) != EOF)
    {
      if(!fres)
	{
	  if(processID == 0)
	    printf("error reading weight file probably encountered a non-integer weight value\n");
	  errorExit(-1);
	}
      wv[i] = w;
      i++;	
    }
    
  if(i != (rdta->sites + 1))
    {
      if(processID == 0)
	printf("number %d of weights not equal to number %d of alignment columns\n", i, rdta->sites);
      errorExit(-1);
    }
 
  for(i = 1; i <= rdta->sites; i++)     
    rdta->wgt[i] = wv[i];         
  
  fclose(weightFile);
  free(wv);
}



static void getinput(analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{ 
  int i;  

  getnums(rdta);

  tr->mxtips         = rdta->numsp;
  rdta->wgt          = (int *)    malloc((rdta->sites + 1) * sizeof(int));
  rdta->wgt2         = (int *)    malloc((rdta->sites + 1) * sizeof(int));
  cdta->alias        = (int *)    malloc((rdta->sites + 1) * sizeof(int)); 
  cdta->aliaswgt     = (int *)    malloc((rdta->sites + 1) * sizeof(int));      
  cdta->rateCategory = (int *)    malloc((rdta->sites + 1) * sizeof(int)); 
  tr->model          = (int *)    calloc((rdta->sites + 1), sizeof(int));
  tr->dataVector     = (int *)    malloc((rdta->sites + 1) * sizeof(int));
  cdta->wr           = (double *) malloc((rdta->sites + 1) * sizeof(double));  
  cdta->wr2          = (double *) malloc((rdta->sites + 1) * sizeof(double));  	
  cdta->patrat       = (double *) malloc((rdta->sites + 1) * sizeof(double));
  cdta->patratStored = (double *) malloc((rdta->sites + 1) * sizeof(double));             

      
  if(!adef->useWeightFile)
    {
      for (i = 1; i <= rdta->sites; i++) 
	rdta->wgt[i] = 1;         
    }
  else    
    inputweights(rdta);   
  
  tr->multiBranch = 0;
  tr->numBranches = 1;

  if(adef->useMultipleModel)  
    {
      int ref;
     
      parsePartitions(adef, rdta, tr);     
     
      for(i = 1; i <= rdta->sites; i++)
	{
	  ref = tr->model[i];
	  tr->dataVector[i] = tr->partitionData[ref].dataType;
	}                     
    }
  else
    {     
      int dataType;

      tr->partitionData  = (pInfo*)malloc(sizeof(pInfo));
      tr->partitionData[0].partitionName = (char*)malloc(128 * sizeof(char));    
      strcpy(tr->partitionData[0].partitionName, "No Name Provided");

      tr->partitionData[0].protModels = adef->proteinMatrix;
      tr->partitionData[0].protFreqs  = adef->protEmpiricalFreqs;


      tr->NumberOfModels = 1;     

      if(adef->model == M_PROTCAT || adef->model == M_PROTGAMMA)
	dataType = AA_DATA;       
      else
	dataType = DNA_DATA;
      
      /* INIT data-type, model, dataVector for good */
      /* those values will be constant throughout the */
      /* inference process */
      
      tr->partitionData[0].dataType = dataType;     
      
      for(i = 0; i <= rdta->sites; i++)
	{
	  tr->dataVector[i] = dataType;      
	  tr->model[i]      = 0;
	}
    }  

#ifdef _MULTI_GENE
  {
    int i;

    tr->startVector = (nodeptr *)malloc(sizeof(nodeptr) * tr->NumberOfModels);
    tr->tipMissing = (char **)malloc(sizeof(char *) * (tr->mxtips + 1));
    
    for(i = 0; i <= tr->mxtips; i++)
      tr->tipMissing[i] = (char *)malloc(sizeof(char) * (tr->NumberOfModels));
  }				     
#endif
 

  getyspace(rdta);
  setupTree(tr, adef);
    
  if(!getdata(adef, rdta, tr))
    {
      printf("Problem reading alignment file \n");
      errorExit(1);
    }
        
  return;
} 


  

static void sitesort(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)   
{ 
  int  gap, i, j, jj, jg, k, n, nsp;
  int  *index, *category, *superCategory;
  boolean  flip, tied;
  char  **data;
    
  if(adef->useMultipleModel)
    {
      superCategory = tr->dataVector;
      category      = tr->model;
    }
  else
    {
      category      = (int*)NULL;
      superCategory = (int*)NULL;
    }

  index    = cdta->alias;
  data     = rdta->y;
  n        = rdta->sites;
  nsp      = rdta->numsp;
  index[0] = -1;

  
  if((adef->mode != OPTIMIZE_RATES))
    {
      for (gap = n / 2; gap > 0; gap /= 2) 
	{
	  for (i = gap + 1; i <= n; i++) 
	    {
	      j = i - gap;
	      
	      do 
		{
		  jj = index[j];
		  jg = index[j+gap];
		  if(adef->useMultipleModel)
		    {	
		      assert(superCategory[jj] != -1 &&
			     superCategory[jg] != -1 &&
			     category[jj] != -1      &&
			     category[jg] != -1);
		      
		      if(superCategory[jj] > superCategory[jg])
			{
			  flip = TRUE;
			  tied = FALSE;
			}
		      else
			{
			  flip = (category[jj] >  category[jg]);
			  tied = (category[jj] == category[jg]);		    
			}
		    }
		  else
		    {		    
		      flip = 0;
		      tied = 1;
		    }
		  
		  for (k = 1; (k <= nsp) && tied; k++) 
		    {
		      flip = (data[k][jj] >  data[k][jg]);
		      tied = (data[k][jj] == data[k][jg]);
		    }
		  
		  if (flip) 
		    {
		      index[j]     = jg;
		      index[j+gap] = jj;
		      j -= gap;
		    }
		} 
	      while (flip && (j > 0));	      
	    }  
	}
    }  
}


static void sitecombcrunch (rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)    
{ 
  int  i, sitei, j, sitej, k;   
  boolean  tied;
  int *aliasModel; 
  int *aliasSuperModel;
    


  if(adef->useMultipleModel)
    {
      aliasSuperModel = (int*)malloc(sizeof(int) * (rdta->sites + 1));
      aliasModel      = (int*)malloc(sizeof(int) * (rdta->sites + 1));
    }
  else
    {
      aliasModel      = (int*)NULL;
      aliasSuperModel = (int*)NULL;
    }
  
  i = 0;    
  cdta->alias[0]    = cdta->alias[1];
  cdta->aliaswgt[0] = 0;   
    
  for (j = 1; j <= rdta->sites; j++) 
    {
      sitei = cdta->alias[i];
      sitej = cdta->alias[j];
      if(adef->mode == OPTIMIZE_RATES)
	tied = 0;
      else
	{	 
	  if(adef->useMultipleModel)
	    {
	      tied = (tr->model[sitei] == tr->model[sitej]);	    
	      if(tied)
		assert(tr->dataVector[sitei] == tr->dataVector[sitej]);  
	    }
	  else	      
	    tied = 1;	      
	}
      
      for (k = 1; tied && (k <= rdta->numsp); k++)
	tied = (rdta->y[k][sitei] == rdta->y[k][sitej]);
      
      if (tied) 
	{
	  cdta->aliaswgt[i] += rdta->wgt2[sitej];	   
	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];	  	     	  
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
      else 
	{
	  if (cdta->aliaswgt[i] > 0) i++;
	  cdta->aliaswgt[i] = rdta->wgt2[sitej];
	  cdta->alias[i] = sitej;
	  if(adef->useMultipleModel)
	    {
	      aliasModel[i]      = tr->model[sitej];		   
	      aliasSuperModel[i] = tr->dataVector[sitej];
	    }
	}
    }

  cdta->endsite = i;
  if (cdta->aliaswgt[i] > 0) cdta->endsite++;       
    
  if(adef->useMultipleModel)
    {       
      for(i = 0; i <= rdta->sites; i++)
	{
	  tr->model[i]      = aliasModel[i];	  
	  tr->dataVector[i] = aliasSuperModel[i];
	}
    }
  
  if(adef->useMultipleModel)
    {
      free(aliasModel);	
      free(aliasSuperModel);
    }

  /*for(i = 0; i < cdta->endsite; i++)
    {
      printf("%d %d %d\n",  i, tr->model[i], tr->dataVector[i]);
      }*/
} 


static boolean makeweights (analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)    
{
  int  i;

 
  for (i = 1; i <= rdta->sites; i++)  
    rdta->wgt2[i] = rdta->wgt[i];     

  for (i = 1; i <= rdta->sites; i++)  
    cdta->alias[i] = i;
        
  sitesort(rdta, cdta, tr, adef);
  sitecombcrunch(rdta, cdta, tr, adef);
      
  return TRUE;
} 




static boolean makevalues(rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef)   
{   
  int  i, j, model, fullSites = 0, modelCounter;   

  char *y    = (char *)malloc(rdta->numsp * cdta->endsite * sizeof(char));
  char *yBUF = (char *)malloc(rdta->numsp * cdta->endsite * sizeof(char));

  for (i = 1; i <= rdta->numsp; i++)       
    for (j = 0; j < cdta->endsite; j++) 	  	          
      y[((i - 1) * cdta->endsite) + j] = rdta->y[i][cdta->alias[j]];	                

  free(rdta->y0);
  free(rdta->y);
    
  rdta->y0 = y;
  memcpy(yBUF, y, rdta->numsp * cdta->endsite * sizeof(char));
  rdta->yBUF = yBUF;
                      
  if(!adef->useMultipleModel)
    tr->NumberOfModels = 1;

  if(adef->useMultipleModel)
    {      
      int *perm          = (int*)malloc(sizeof(int) * tr->NumberOfModels);     
      pInfo *partitionData = (pInfo*)malloc(sizeof(pInfo) * tr->NumberOfModels);
     
      tr->partitionData[0].lower = 0;
           
      model        = tr->model[0];
      modelCounter = 0;
      perm[modelCounter] = model;
      i            = 1;                 

      while(i <  cdta->endsite)
	{
	  if(tr->model[i] != model)
	    {	     
	      tr->partitionData[modelCounter].upper     = i;
	      tr->partitionData[modelCounter + 1].lower = i;

	      model = tr->model[i];
	      perm[modelCounter + 1] = model;
	      modelCounter++;
	    }
	  i++;
	}
      
     
      tr->partitionData[tr->NumberOfModels - 1].upper = cdta->endsite;
     
      memcpy(partitionData, tr->partitionData, tr->NumberOfModels * sizeof(pInfo));
      
      for(i = 0; i < tr->NumberOfModels; i++)
	{	 	  
	  tr->partitionData[i].dataType   = partitionData[perm[i]].dataType;
	  tr->partitionData[i].protModels = partitionData[perm[i]].protModels;
	  tr->partitionData[i].protFreqs  = partitionData[perm[i]].protFreqs;
	}            

      model        = tr->model[0];
      modelCounter = 0;
      tr->model[0] = modelCounter;    
      i            = 1;                 

      while(i < cdta->endsite)
	{
	  if(tr->model[i] != model)
	    {	     
	      model = tr->model[i];	      
	      modelCounter++;
	      tr->model[i] = modelCounter;
	    }
	  else
	    tr->model[i] = modelCounter;
	  i++;
	}      

      for(i = 0; i < (cdta->endsite - 1); i++)
	{
	  if(tr->dataVector[i] != tr->dataVector[i + 1])
	    {
	      tr->mixedData = TRUE;	      
	      break;
	    }
	}

      free(perm);       
      free(partitionData);
    }
  else
    {                
      tr->partitionData[0].lower = 0;
      tr->partitionData[0].upper = cdta->endsite;     
    }

#ifdef _USE_PTHREADS
  /* 
     TODO-MIX Have to seriously think about whether 
     to implement this or not 
  */
     
  if(tr->mixedData)
    {
      printf("No Pthreads implementation for mixed data available right now \n");
      assert(0);
    }
#endif


  tr->rdta       = rdta;
  tr->cdta       = cdta;   
  
  tr->invariant          = (int *)malloc(cdta->endsite * sizeof(int));
  tr->originalDataVector = (int *)malloc(cdta->endsite * sizeof(int));
  tr->originalModel      = (int *)malloc(cdta->endsite * sizeof(int));
  tr->originalWeights    = (int *)malloc(cdta->endsite * sizeof(int));

  memcpy(tr->originalModel, tr->model,            cdta->endsite * sizeof(int)); 
  memcpy(tr->originalDataVector, tr->dataVector,  cdta->endsite * sizeof(int));
  memcpy(tr->originalWeights, tr->cdta->aliaswgt, cdta->endsite * sizeof(int));
  
  tr->originalCrunchedLength = tr->cdta->endsite;
  for(i = 0; i < tr->cdta->endsite; i++)
    fullSites += tr->cdta->aliaswgt[i];

  tr->fullSites = fullSites;
  
  for(i = 0; i < rdta->numsp; i++)
    tr->yVector[i + 1] = &(rdta->y0[tr->originalCrunchedLength * i]);

  return TRUE;
} 







static int sequenceSimilarity(char *tipJ, char *tipK, int n)
{
  int i;
  
  for(i = 0; i < n; i++)    
    if(*tipJ++ != *tipK++)	
      return 0;	
      
  return 1;
}

static void checkSequences(tree *tr, rawdata *rdta, analdef *adef)
{
  int n = tr->mxtips + 1; 
  int i, j;
  int *omissionList     = (int *)malloc(n * sizeof(int));
  int *undeterminedList = (int *)malloc((rdta->sites + 1)* sizeof(int));
  int *modelList        = (int *)malloc((rdta->sites + 1)* sizeof(int)); 
  int count = 0;
  int countNameDuplicates = 0;
  int countUndeterminedColumns = 0;
  int countOnlyGaps = 0;
  int modelCounter = 1;
  char undetermined_AA  = 22;
  char undetermined_DNA = 15;
  char *tipI, *tipJ;
  FILE *f;  


  if(processID == 0)	      
    f = fopen(infoFileName, "a");
  else
    f = (FILE *)NULL; 

  for(i = 1; i < n; i++)       
    omissionList[i] = 0;              

  for(i = 0; i < rdta->sites + 1; i++)
    undeterminedList[i] = 0;
      
  for(i = 1; i < n; i++)
    {
      for(j = i + 1; j < n; j++)
	if(strcmp(tr->nameList[i], tr->nameList[j]) == 0)
	  {
	    countNameDuplicates++;
	    if(processID == 0)
	      {
		printf("Sequence names of taxon %d and %d are identical, they are both called %s\n", i, j, tr->nameList[i]);
		fprintf(f, "Sequence names of taxon %d and %d are identical, they are both called %s\n", i, j, tr->nameList[i]);
	      }
	  }
    }
	  
  if(countNameDuplicates > 0)
    {
      if(processID == 0)
	{
	  printf("ERROR: Found %d taxa that had equal names in the alignment, exiting...\n", countNameDuplicates);
	  fprintf(f, "ERROR: Found %d taxa that had equal names in the alignment, exiting...\n", countNameDuplicates);
	  fclose(f);
	}
      errorExit(-1);
    }

  for(i = 1; i < n; i++)
    {
      j = 1;
      
      while(j <= rdta->sites)
	{
	  if(tr->dataVector[j] == DNA_DATA && rdta->y[i][j] != undetermined_DNA)
	    break;
	  if(tr->dataVector[j] == AA_DATA && rdta->y[i][j] != undetermined_AA)
	    break;
	  j++;
	}

      if(j == (rdta->sites + 1))
	{       
	  if(processID == 0)
	    {
	      printf("ERROR: Sequence %s consists entirely of undetermined values which will be treated as missing data\n",      
		     tr->nameList[i]);
	      fprintf(f, "ERROR: Sequence %s consists entirely of undetermined values which will be treated as missing data\n",      
		      tr->nameList[i]);	      
	    }
	  countOnlyGaps++;
	}
      
    }
  
  if(countOnlyGaps > 0)
    {
      if(processID == 0)
	{
	  printf("ERROR: Found %d sequences that consist entirely of undetermined values, exiting...\n", countOnlyGaps);
	  fprintf(f, "ERROR: Found %d sequences that consist entirely of undetermined values, exiting...\n", countOnlyGaps);
	  fclose(f);
	}
      errorExit(-1);
    }

  for(i = 0; i <= rdta->sites; i++)
    modelList[i] = -1;

  for(i = 1; i <= rdta->sites; i++)
    {    
      j = 1;
     
      while(j < n)
	{
	  if(tr->dataVector[i] == DNA_DATA && rdta->y[j][i] != undetermined_DNA)
	    break;
	  if(tr->dataVector[i] == AA_DATA && rdta->y[j][i] != undetermined_AA)
	    break;	
	  j++;
	}
      
      if(j == n)
	{
	  undeterminedList[i] = 1;
	  if(processID == 0)
	    {
	      printf("IMPORTANT WARNING: Alignment column %d contains only undetermined values which will be treated as missing data\n", 
		     i);
	      fprintf(f, "IMPORTANT WARNING: Alignment column %d contains only undetermined values which will be treated as missing data\n", i);
	    }
	  countUndeterminedColumns++;	  
	}
      else
	{
	  if(adef->useMultipleModel)
	    {
	      modelList[modelCounter] = tr->model[i];
	      modelCounter++;
	    }
	}
    }
  

  for(i = 1; i < n; i++)
    {
      if(omissionList[i] == 0)
	{
	  tipI = &(rdta->y[i][1]);

	  for(j = i + 1; j < n; j++)
	    {
	      if(omissionList[j] == 0)
		{
		  tipJ = &(rdta->y[j][1]);
		  if(sequenceSimilarity(tipI, tipJ, rdta->sites))
		    {
		      if(processID == 0)
			{
			  printf("\n\nIMPORTANT WARNING: Sequences %s and %s are exactly identical\n", tr->nameList[i], tr->nameList[j]);
			  fprintf(f, "\n\nIMPORTANT WARNING: Sequences %s and %s are exactly identical\n", tr->nameList[i], tr->nameList[j]);
			}
		      omissionList[j] = 1;
		      count++;
		    }
		}
	    }
	}
    }

  if(count > 0 || countUndeterminedColumns > 0)
    {
      char noDupFile[2048];
      char noDupModels[2048];
              
      if(count > 0)
	{
	  if(processID == 0)
	    {
	      printf("\n");
	      
	      printf("IMPORTANT WARNING\n");
	      
	      printf("Found %d %s that %s exactly identical to other sequences in the alignment.\n", count, (count == 1)?"sequence":"sequences", (count == 1)?"is":"are");
	      printf("Normally they should be excluded from the analysis.\n\n");
	      
	      fprintf(f, "\n");
	      
	      fprintf(f, "IMPORTANT WARNING\n");
	      
	      fprintf(f, "Found %d %s that %s exactly identical to other sequences in the alignment.\n", count, (count == 1)?"sequence":"sequences", (count == 1)?"is":"are");
	      fprintf(f, "Normally they should be excluded from the analysis.\n\n");
	    }
	}
      
      if(countUndeterminedColumns > 0)
	{
	  if(processID == 0)
	    {
	      printf("\n");
	      
	      printf("IMPORTANT WARNING\n");
	      
	      printf("Found %d %s that %s only undetermined values which will be treated as missing data.\n", 
		     countUndeterminedColumns, (countUndeterminedColumns == 1)?"column":"columns", (countUndeterminedColumns == 1)?"contains":"contain");
	      printf("Normally these columns should be excluded from the analysis.\n\n");
	      
	      fprintf(f, "\n");
	      
	      fprintf(f, "IMPORTANT WARNING\n");
	      
	      fprintf(f, "Found %d %s that %s only undetermined values which will be treated as missing data.\n", 
		      countUndeterminedColumns, (countUndeterminedColumns == 1)?"column":"columns", (countUndeterminedColumns == 1)?"contains":"contain");
	      fprintf(f, "Normally these columns should be excluded from the analysis.\n\n");      	  
	    }
	}

      strcpy(noDupFile, seq_file);
      strcat(noDupFile, ".reduced");

      strcpy(noDupModels, modelFileName);
      strcat(noDupModels, ".reduced");

      if(processID == 0)
	{

	  if(adef->useMultipleModel && !filexists(noDupModels) && countUndeterminedColumns)
	    {      
	      FILE *newFile = fopen(noDupModels, "w");

	      printf("\nJust in case you might need it, a mixed model file with \n");
	      printf("model assignments for undetermined columns removed is printed to file %s\n",noDupModels);

	      fprintf(f, "\nJust in case you might need it, a mixed model file with \n");
	      fprintf(f, "model assignments for undetermined columns removed is printed to file %s\n",noDupModels);
	      
 
	      for(i = 0; i < tr->NumberOfModels; i++)
		{
		  boolean modelStillExists = FALSE;
		  
		  for(j = 1; (j <= rdta->sites) && (!modelStillExists); j++)
		    {
		      if(modelList[j] == i)
			modelStillExists = TRUE;
		    }

		  if(modelStillExists)
		    {	  
		      char *protModels[10] = {"DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM"};
		      int k = 1;
		      int lower, upper;
		      int parts = 0;


		      switch(tr->partitionData[i].dataType)
			{
			case AA_DATA:		      		     
			  {
			    char AAmodel[1024];
			    
			    strcpy(AAmodel, protModels[tr->partitionData[i].protModels]);
			    if(tr->partitionData[i].protFreqs)
			      strcat(AAmodel, "F");		  
			  
			    fprintf(newFile, "%s, ", AAmodel);
			  }
			  break;
			case DNA_DATA:
			  fprintf(newFile, "DNA, ");
			  break;
			default:
			  assert(0);
			}

		      fprintf(newFile, "%s = ", tr->partitionData[i].partitionName);
		      
		      while(k <= rdta->sites)
			{
			  if(modelList[k] == i)
			    {
			      lower = k;
			      while((modelList[k + 1] == i) && (k <= rdta->sites))		      			
				k++;
			      upper = k;
			      
			      if(lower == upper)		  
				{
				  if(parts == 0)
				    fprintf(newFile, "%d", lower);
				  else
				    fprintf(newFile, ",%d", lower);
				}
			      else
				{
				  if(parts == 0)
				    fprintf(newFile, "%d-%d", lower, upper);
				  else
				    fprintf(newFile, ",%d-%d", lower, upper);
				}		  
			      parts++;
			    }
			  k++;
			}
		      fprintf(newFile, "\n");
		    }		  
		}	
	      fclose(newFile);
	    }
	  else
	    {
	      if(adef->useMultipleModel)
		{
		  printf("\n A mixed model file with model assignments for undetermined\n");
		  printf("columns removed has already been printed to  file %s\n",noDupModels);

		  fprintf(f, "\n A mixed model file with model assignments for undetermined\n");
		  fprintf(f, "columns removed has already been printed to  file %s\n",noDupModels);
		}	      
	    }
	     

	  if(!filexists(noDupFile))
	    {
	      FILE *newFile;
	      
	      printf("Just in case you might need it, an alignment file with \n");
	      if(count && !countUndeterminedColumns)
		printf("sequence duplicates removed is printed to file %s\n", noDupFile);
	      if(!count && countUndeterminedColumns)
		printf("undetermined columns removed is printed to file %s\n", noDupFile);
	      if(count && countUndeterminedColumns)
		printf("sequence duplicates and undetermined columns removed is printed to file %s\n", noDupFile);
	      
	      fprintf(f, "Just in case you might need it, an alignment file with \n");
	      if(count && !countUndeterminedColumns)
		fprintf(f, "sequence duplicates removed is printed to file %s\n", noDupFile);
	      if(!count && countUndeterminedColumns)
		fprintf(f, "undetermined columns removed is printed to file %s\n", noDupFile);
	      if(count && countUndeterminedColumns)
		fprintf(f, "sequence duplicates and undetermined columns removed is printed to file %s\n", noDupFile);
	      
	      newFile = fopen(noDupFile, "w");
	      
	      fprintf(newFile, "%d %d\n", tr->mxtips - count, rdta->sites - countUndeterminedColumns);
	      
	      for(i = 1; i < n; i++)
		{
		  if(!omissionList[i])
		    {
		      fprintf(newFile, "%s ", tr->nameList[i]);
		      tipI =  &(rdta->y[i][1]);

		      for(j = 0; j < rdta->sites; j++)
			{
			  if(undeterminedList[j + 1] == 0)
			    {
			      switch(tr->dataVector[j + 1])
				{
				case AA_DATA:
				  fprintf(newFile, "%c", inverseMeaningPROT[tipI[j]]);
				  break;
				case DNA_DATA:
				  fprintf(newFile, "%c", inverseMeaningDNA[tipI[j]]);
				  break;
				default:
				  assert(0);
				}
			    }
			}
		      		      
		      fprintf(newFile, "\n");
		    }
		}
	      
	      fclose(newFile);	    
	    }
	  else
	    {
	      if(count && !countUndeterminedColumns)
		printf("An alignment file with sequence duplicates removed has already\n");
	      if(!count && countUndeterminedColumns)
		printf("An alignment file with undetermined columns removed has already\n");
	      if(count && countUndeterminedColumns)
		printf("An alignment file with undetermined columns and sequence duplicates removed has already\n");
	      
	      printf("been printed to file %s\n",  noDupFile);
	      
	      if(count && !countUndeterminedColumns)
		fprintf(f, "An alignment file with sequence duplicates removed has already\n");
	      if(!count && countUndeterminedColumns)
		fprintf(f, "An alignment file with undetermined columns removed has already\n");
	      if(count && countUndeterminedColumns)
		fprintf(f, "An alignment file with undetermined columns and sequence duplicates removed has already\n");
	      
	      fprintf(f, "been printed to file %s\n",  noDupFile);
	    }
	}
    }


  free(undeterminedList);
  free(omissionList);
  free(modelList);
  if(processID == 0)	      
    fclose(f);
}



static float dist(int i, int j, const int sites, const float nDouble, char **y)
{
  int k, count;  
  char *tipI = &(y[i + 1][1]);
  char *tipJ = &(y[j + 1][1]); 

  for(k = 0, count = 0; k < sites; k ++)
    if(tipI[k] == tipJ[k])
      count++;
  
  return (((float)count) * nDouble);
}

static void distArray(int i, const int sites, const float nDouble, int n, float *ref, int *omitted, char **y)
{
  int k, l;  
  char *tipI = &(y[i + 1][1]); 
  
  for(l = 0; l < n; l++)
    {
      if((!omitted[l]) && (l != i))
	{
	  char *tipJ = &(y[l + 1][1]); 
	  int count = 0;
	  for(k = 0, count = 0; k < sites; k ++)
	    if(tipI[k] == tipJ[k])
	      count++;
	  ref[l] = (((float)count) * nDouble);
	}
    }
}



static int qtCompare(const void *p1, const void *p2)
{
 qtData *rc1 = (qtData *)p1;
 qtData *rc2 = (qtData *)p2;

  float i = rc1->val;
  float j = rc2->val;
  
  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}


static qtList * clusterQT_LARGE(int n, float thres, int *ccc, rawdata *rdta)
{  
  int clusterCount;
  int i, j;
  int *omitted, *current, *best;  
  qtList *clusters = (qtList*)NULL;
  const float nDouble = 1.0 / (float)(rdta->sites);
  double t = gettime();  
  const int sites = rdta->sites; 
  char **y = rdta->y; 
  float *ref;
  qtData *candidates;
  
  candidates = (qtData *)malloc(sizeof(qtData) * n);
  clusters = (qtList *)malloc(sizeof(qtList) * n);
  omitted = (int*)calloc(n, sizeof(int));
  current = (int*)malloc(sizeof(int) * n);
  best    = (int*)malloc(sizeof(int) * n);            
  ref     = (float*)malloc(sizeof(float) * n);
  clusterCount = 0;


  for(i = 0; i < n; i++)
    {
      if(!omitted[i])
	{
	  int entCount = 0;	  	     
	  int countCandidates = 0;	  

	  current[entCount++] = i;	      
	  omitted[i] = 1;

	  distArray(i, sites, nDouble, n, ref, omitted, y);       
	 		 
	  for(j = 0; j < n; j++)		
	    {
	      if(!omitted[j] && i != j)
		{
		  float temp;			
		  
		  if((temp = ref[j]) >= thres)
		    {
		      candidates[countCandidates].val    = temp;
		      candidates[countCandidates].number = j;			
		      countCandidates++;
		    }
		}
	    }
		  
	  if(countCandidates > 0)
	    {
	      qsort(candidates, countCandidates, sizeof(qtData), qtCompare);	 	      
	      
	      for(j = 0; j < countCandidates; j++)
		{
		  int k;	       
		  
		  for(k = 0; k < entCount; k++)					     
		    if(dist(current[k], candidates[j].number, sites, nDouble, y) < thres)
		      break;
		  
		  if(k == entCount)
		    {
		      current[entCount++] = candidates[j].number;		     		    		     
		      omitted[candidates[j].number] = 1;
		    }	      	      	
		}
	    }
	  
	  clusters[clusterCount].entries = (int *)malloc(sizeof(int) * entCount);  	      
	  memcpy(clusters[clusterCount].entries, current, entCount * sizeof(int));
	  clusters[clusterCount++].count = entCount;
	}
    }
	    
  printf("Time %f\n", gettime() - t);
  printf("FOUND %d Clusters\n", clusterCount);  
    

   if(1)
    {
      int ver = 0;
      int check = 0;
      int total = 0;
      for(i = 0; i < n; i++)
	ver += i;
      
      for(i = 0; i < clusterCount; i++)
	{	 	  
	  {
	    int k;
	    for(j = 0; j < clusters[i].count; j++)
	      for(k = 0; k < clusters[i].count; k++)	   	    
		assert(dist(clusters[i].entries[j], clusters[i].entries[k],sites, nDouble, y)  >=  thres);		
	  }
	  
	  for(j = 0; j < clusters[i].count; j++)
	    {
	      check += clusters[i].entries[j];	     
	      total++;
	    }	
	}
      assert(ver == check);          
      printf("Total: %d\n", total);
    }



  for(i = 0; i < clusterCount; i++)
    {
      float max  = 0.0;
      int length = clusters[i].count;
      int pos    = -1;      
      int *c     = clusters[i].entries;
      int buf;

      if(length > 2)
	{
	  for(j = 0; j < length; j++)
	    {
	      int k;
	      float avg = 0.0;

	      for(k = 0; k < length; k++)
		{
		  if(j != k)		    		     
		    avg += dist(c[j], c[k], sites, nDouble, y);	
		}	      	     
	       
	      if(avg > max)
		{
		  max = avg;
		  pos = j;
		}
	    }	  
	  if(pos > 0)
	    {	      
	      buf    = c[0];
	      c[0]   = c[pos];
	      c[pos] = buf;	 
	    }
	}
      for(j = 0; j < length; j++)
	c[j] = c[j] + 1;	
    }

  free(candidates);
  free(omitted);
  free(current);
  free(best);
  free(ref);
  *ccc = clusterCount;
  return clusters;
}



static qtList * clusterQT(float **d, int n, float thres, int *ccc)
{  
  int clusterCount;
  int i, j;
  int *omitted, *current, *best;  
  int total = 0;  
  qtList *clusters = (qtList*)NULL;
  double t = gettime();
  
  clusters = (qtList *)malloc(sizeof(qtList) * n);
  omitted = (int*)calloc(n, sizeof(int));
  current = (int*)malloc(sizeof(int) * n);
  best    = (int*)malloc(sizeof(int) * n);            

  clusterCount = 0;

  while(1)
    {
      int max = -1;
      int maxPos = -1;

      for(i = 0; i < n; i++)
	{
	  if(!omitted[i])
	    {
	      int entCount = 0;	  	     
	      int *inSet = (int *)calloc(n, sizeof(int));	      
	      boolean aboveThres = TRUE;

	      current[entCount++] = i;
	      inSet[i] = 1;

	      while(aboveThres)
		{	
		  float dm = -1.0;
		  int dmPos = -1;
	 		 
		  for(j = 0; j < n; j++)		
		    if(i != j && (!omitted[j]) && (!inSet[j]) && d[i][j] > dm)
		      {
			dm    = d[i][j];
			dmPos = j;
		      }
		  
		  if(dmPos == -1)
		    aboveThres = FALSE;
		  else
		    {
		      for(j = 0; j < entCount && aboveThres; j++)			
			if(d[current[j]][dmPos] < thres)
			  aboveThres = FALSE;
		      
		      if(aboveThres)
			{
			  current[entCount++] = dmPos;
			  inSet[dmPos] = 1;
			}
		    }
		}	      	      	      
	      
	      if(entCount > max)
		{
		  max    = entCount;
		  maxPos = i;
		  memcpy(best, current, entCount * sizeof(int));
		}
	      free(inSet);
	    }
	}

      if(maxPos == -1)
	break;

      clusters[clusterCount].entries = (int *)malloc(sizeof(int) * max);
      memcpy(clusters[clusterCount].entries, best, max * sizeof(int));

      for(i = 0; i < max; i++)			 
	omitted[best[i]]    = 1;
	
      clusters[clusterCount++].count = max;               
    }
  
  printf("Time %f\n", gettime() - t);
  printf("FOUND %d Clusters\n", clusterCount);
  
  if(1)
    {
      int ver = 0;
      int check = 0;
      for(i = 0; i < n; i++)
	ver += i;
      
      for(i = 0; i < clusterCount; i++)
	{
	  /*printf("Cluster %d:", i);*/
	  
	  {
	    int k;
	    for(j = 0; j < clusters[i].count; j++)
	      for(k = 0; k < clusters[i].count; k++)	   	    
		assert(d[clusters[i].entries[j]][clusters[i].entries[k]] >=  thres);		
	  }
	  
	  for(j = 0; j < clusters[i].count; j++)
	    {
	      check += clusters[i].entries[j];
	      /*printf("%d ", clusters[i].entries[j]);*/
	      total++;
	    }
	  /*printf("\n");*/
	}
      assert(ver == check);
      
      /*printf("TOTAL: %d\n", total);*/
    }

  for(i = 0; i < clusterCount; i++)
    {
      float max  = 0.0;
      int length = clusters[i].count;
      int pos    = -1;    
      int *c     = clusters[i].entries;
      int buf;

      if(length > 2)
	{
	  for(j = 0; j < length; j++)
	    {
	      int k;
	      float avg = 0.0;
	      for(k = 0; k < length; k++)
		{		  
		  if(j != k)		    		      
		    avg += d[c[j]][c[k]];		   
		}	      	      	      

	      if(avg > max)
		{
		  max = avg;
		  pos = j;
		}
	    }	  
	  /*printf("Cluster %d length %d avg %f\n", i, length, max);*/
	  
	  if(pos > 0)
	    {
	      /*printf("Cluster %d siwtching %d <-> %d\n", i, 0, pos);*/
	      buf    = c[0];
	      c[0]   = c[pos];
	      c[pos] = buf;
	    }
	}
      for(j = 0; j < length; j++)
	c[j] = c[j] + 1;	
    }

  free(omitted);
  free(current);
  free(best);
  *ccc = clusterCount;
  return clusters;
}

static void reduceBySequenceSimilarity(tree *tr, rawdata *rdta, analdef *adef)
{
  int n = tr->mxtips + 1; 
  int i, j;
  int *omissionList     = (int *)malloc(n * sizeof(int));
  int *undeterminedList = (int *)malloc((rdta->sites + 1)* sizeof(int));
  int *modelList        = (int *)malloc((rdta->sites + 1)* sizeof(int));  
  int countNameDuplicates = 0;
  int countUndeterminedColumns = 0;
  int countOnlyGaps = 0;
  int modelCounter = 1;
  char buf[16], outName[1024];
  char undetermined_AA  = 22;
  char undetermined_DNA = 15;
  char *tipI;
  qtList *clusters = (qtList*)NULL;
  FILE *f, *assoc;  
  int numberOfClusters = 0;
  int nonTrivial = 0;
  
  strcpy(outName,         workdir);  
  strcat(outName,         "RAxML_reducedList.");
  strcat(outName,         run_id);

  if(processID == 0)	      
    f = fopen(infoFileName, "a");
  else
    f = (FILE *)NULL;

 

  for(i = 1; i < n; i++)       
    omissionList[i] = 0;              

  for(i = 0; i < rdta->sites + 1; i++)
    undeterminedList[i] = 0;
      
  for(i = 1; i < n; i++)
    {
      for(j = i + 1; j < n; j++)
	if(strcmp(tr->nameList[i], tr->nameList[j]) == 0)
	  {
	    countNameDuplicates++;
	    if(processID == 0)
	      {
		printf("Sequence names of taxon %d and %d are identical, they are both called %s\n", i, j, tr->nameList[i]);
		fprintf(f, "Sequence names of taxon %d and %d are identical, they are both called %s\n", i, j, tr->nameList[i]);
	      }
	  }
    }
	  
  if(countNameDuplicates > 0)
    {
      if(processID == 0)
	{
	  printf("ERROR: Found %d taxa that had equal names in the alignment, exiting...\n", countNameDuplicates);
	  fprintf(f, "ERROR: Found %d taxa that had equal names in the alignment, exiting...\n", countNameDuplicates);
	  fclose(f);
	}
      errorExit(-1);
    }

  for(i = 1; i < n; i++)
    {
      j = 1;

      while(j <= rdta->sites)
	{
	  if(tr->dataVector[j] == DNA_DATA && rdta->y[i][j] != undetermined_DNA)
	    break;
	  if(tr->dataVector[j] == AA_DATA && rdta->y[i][j] != undetermined_AA)
	    break;
	  j++;
	}     

      if(j == (rdta->sites + 1))
	{       
	  if(processID == 0)
	    {
	      printf("ERROR: Sequence %s consists entirely of undetermined values which will be treated as missing data\n",      tr->nameList[i]);
	      fprintf(f, "ERROR: Sequence %s consists entirely of undetermined values which will be treated as missing data\n",      tr->nameList[i]);	      
	    }
	  countOnlyGaps++;
	}
      
    }
  
  if(countOnlyGaps > 0)
    {
      if(processID == 0)
	{
	  printf("ERROR: Found %d sequences that consist entirely of undetermined values, exiting...\n", countOnlyGaps);
	  fprintf(f, "ERROR: Found %d sequences that consist entirely of undetermined values, exiting...\n", countOnlyGaps);
	  fclose(f);
	}
      errorExit(-1);
    }

  for(i = 0; i <= rdta->sites; i++)
    modelList[i] = -1;

  for(i = 1; i <= rdta->sites; i++)
    {    
      j = 1;
     
      while(j < n)
	{
	  if(tr->dataVector[i] == DNA_DATA && rdta->y[j][i] != undetermined_DNA)
	    break;
	  if(tr->dataVector[i] == AA_DATA && rdta->y[j][i] != undetermined_AA)
	    break;	
	  j++;
	}
      
      if(j == n)
	{
	  undeterminedList[i] = 1;
	  if(processID == 0)
	    {
	      printf("IMPORTANT WARNING: Alignment column %d contains only undetermined values which will be treated as missing data\n", i);
	      fprintf(f, "IMPORTANT WARNING: Alignment column %d contains only undetermined values which will be treated as missing data\n", i);
	    }
	  countUndeterminedColumns++;	  
	}
      else
	{
	  if(adef->useMultipleModel)
	    {
	      modelList[modelCounter] = tr->model[i];
	      modelCounter++;
	    }
	}
    }  

  switch(adef->similarityFilterMode)
    {
    case SMALL_DATA:
      {
	float **d;   
	int n = tr->mxtips;
	int i, j;
	double t = gettime();
	float nDouble = 1.0 / (float)(rdta->sites);    
	int sites = rdta->sites;
	char *tipI, *tipJ;
	
	
	d = (float **)malloc(sizeof(float *) * n);
	for(i = 0; i < n; i++)
	  d[i] = (float *)malloc(sizeof(float) * n);
	
	for(i = 0; i < n; i++)
	  {
	    d[i][i] = 1.0;	
	    tipI = &(rdta->y[i + 1][1]);
	    for(j = i + 1; j < n; j++)
	      {
		int k;
		int count = 0;
		tipJ = &(rdta->y[j + 1][1]);	    
		for(k = 0; k < sites; k++)
		  if(tipJ[k] == tipI[k])
		    count++;	   	   	   
		
		d[i][j] = ((float)count * nDouble);
		d[j][i] = d[i][j];
	      }
	  }
	
	printf("DistMat %f\n", gettime() - t);
		   
	t = gettime();
	clusters = clusterQT(d, n, (float)(adef->sequenceSimilarity), &numberOfClusters);
	printf("QT %f %d\n", gettime() - t, numberOfClusters);       
      }
      break;
    case LARGE_DATA:
      {	  
	double t;

	t = gettime();
	clusters = clusterQT_LARGE(tr->mxtips, (float)(adef->sequenceSimilarity), &numberOfClusters, rdta);
	printf("QT %f %d\n", gettime() - t, numberOfClusters);             
      }
      break;
    default:
      assert(0);
    }

  assoc = fopen(outName, "w");

  for(i = 0; i < numberOfClusters; i++)
    {
      int length = clusters[i].count;
      int *c     = clusters[i].entries;
      int j;
      
      if(length > 1)
	{	 
	  fprintf(assoc, "%s:%s", tr->nameList[c[0]], tr->nameList[c[1]]);
	  for(j = 2; j < length; j++)
	    fprintf(assoc, ",%s", tr->nameList[c[j]]);
	  fprintf(assoc, "\n");
	 
	  nonTrivial++;
	}		      		     
    }

  fclose(assoc);
  

  if(nonTrivial > 0 || countUndeterminedColumns > 0)
    {
      char noDupFile[2048];
      char noDupModels[2048];
              
      if(nonTrivial > 0)
	{
	  if(processID == 0)
	    {
	      printf("\n");	      	   
	      
	      printf("Found %d non-trival clusters, reduction to %d sequences\n", nonTrivial, numberOfClusters);
	      
	      fprintf(f, "\n");
	      
	      fprintf(f, "Found %d non-trival clusters, reduction to %d sequences\n", nonTrivial, numberOfClusters);
	    }
	}
      
      if(countUndeterminedColumns > 0)
	{
	  if(processID == 0)
	    {
	      printf("\n");
	      
	      printf("IMPORTANT WARNING\n");
	      
	      printf("Found %d %s that %s only undetermined values which will be treated as missing data.\n", 
		     countUndeterminedColumns, (countUndeterminedColumns == 1)?"column":"columns", (countUndeterminedColumns == 1)?"contains":"contain");
	      printf("Normally these columns should be excluded from the analysis.\n\n");
	      
	      fprintf(f, "\n");
	      
	      fprintf(f, "IMPORTANT WARNING\n");
	      
	      fprintf(f, "Found %d %s that %s only undetermined values which will be treated as missing data.\n", 
		      countUndeterminedColumns, (countUndeterminedColumns == 1)?"column":"columns", (countUndeterminedColumns == 1)?"contains":"contain");
	      fprintf(f, "Normally these columns should be excluded from the analysis.\n\n");      	  
	    }
	}

      sprintf(buf, "%f", adef->sequenceSimilarity);

      strcpy(noDupFile, seq_file);
      strcat(noDupFile, ".reducedBy.");
      strcat(noDupFile, buf);
      

      strcpy(noDupModels, modelFileName);
      strcat(noDupModels, ".reducedBy.");
      strcat(noDupModels, buf);
	      

      if(processID == 0)
	{
	  if(adef->useMultipleModel && !filexists(noDupModels) && countUndeterminedColumns)
	    {      
	      FILE *newFile = fopen(noDupModels, "w");

	      printf("\nJust in case you might need it, a mixed model file with \n");
	      printf("model assignments for undetermined columns removed is printed to file %s\n",noDupModels);

	      fprintf(f, "\nJust in case you might need it, a mixed model file with \n");
	      fprintf(f, "model assignments for undetermined columns removed is printed to file %s\n",noDupModels);
	      
 
	      for(i = 0; i < tr->NumberOfModels; i++)
		{
		  boolean modelStillExists = FALSE;
		  
		  for(j = 1; (j <= rdta->sites) && (!modelStillExists); j++)
		    {
		      if(modelList[j] == i)
			modelStillExists = TRUE;
		    }

		  if(modelStillExists)
		    {	  
		      char *protModels[10] = {"DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM"};
		      int k = 1;
		      int lower, upper;
		      int parts = 0;		      		     		      

		      switch(tr->partitionData[i].dataType)
			{
			case AA_DATA:		      		     
			  {
			    char AAmodel[1024];
			    
			    strcpy(AAmodel, protModels[tr->partitionData[i].protModels]);
			    if(tr->partitionData[i].protFreqs)
			      strcat(AAmodel, "F");		  
			  
			    fprintf(newFile, "%s, ", AAmodel);
			  }
			  break;
			case DNA_DATA:
			  fprintf(newFile, "DNA, ");
			  break;
			default:
			  assert(0);
			}

		      fprintf(newFile, "%s = ", tr->partitionData[i].partitionName);


		      while(k <= rdta->sites)
			{
			  if(modelList[k] == i)
			    {
			      lower = k;
			      while((modelList[k + 1] == i) && (k <= rdta->sites))		      			
				k++;
			      upper = k;
			      
			      if(lower == upper)		  
				{
				  if(parts == 0)
				    fprintf(newFile, "%d", lower);
				  else
				    fprintf(newFile, ",%d", lower);
				}
			      else
				{
				  if(parts == 0)
				    fprintf(newFile, "%d-%d", lower, upper);
				  else
				    fprintf(newFile, ",%d-%d", lower, upper);
				}		  
			      parts++;
			    }
			  k++;
			}
		      fprintf(newFile, "\n");
		    }		  
		}	
	      fclose(newFile);
	    }
	  else
	    {
	      if(adef->useMultipleModel)
		{
		  printf("\n A mixed model file with model assignments for undetermined\n");
		  printf("columns removed has already been printed to  file %s\n",noDupModels);

		  fprintf(f, "\n A mixed model file with model assignments for undetermined\n");
		  fprintf(f, "columns removed has already been printed to  file %s\n",noDupModels);
		}	      
	    }
	     

	  if(!filexists(noDupFile))
	    {
	      FILE *newFile;
	      
	      printf("Just in case you might need it, an alignment file with \n");
	      if(nonTrivial && !countUndeterminedColumns)
		printf("similar sequences removed is printed to file %s\n", noDupFile);
	      if(!nonTrivial && countUndeterminedColumns)
		printf("undetermined columns removed is printed to file %s\n", noDupFile);
	      if(nonTrivial && countUndeterminedColumns)
		printf("similar sequences and undetermined columns removed is printed to file %s\n", noDupFile);
	      
	      fprintf(f, "Just in case you might need it, an alignment file with \n");
	      if(nonTrivial && !countUndeterminedColumns)
		fprintf(f, "similar sequences removed is printed to file %s\n", noDupFile);
	      if(!nonTrivial && countUndeterminedColumns)
		fprintf(f, "undetermined columns removed is printed to file %s\n", noDupFile);
	      if(nonTrivial && countUndeterminedColumns)
		fprintf(f, "similar sequences and undetermined columns removed is printed to file %s\n", noDupFile);
	      
	      newFile = fopen(noDupFile, "w");
	      
	      fprintf(newFile, "%d %d\n", numberOfClusters, rdta->sites - countUndeterminedColumns);
	      
	      for(i = 0; i < numberOfClusters; i++)
		{
		  
		  fprintf(newFile, "%s ", tr->nameList[clusters[i].entries[0]]);
		  tipI =  &(rdta->y[clusters[i].entries[0]][1]);

		   for(j = 0; j < rdta->sites; j++)
		     {
		       if(undeterminedList[j + 1] == 0)
			 {
			   switch(tr->dataVector[j + 1])
			     {
			     case AA_DATA:
			       fprintf(newFile, "%c", inverseMeaningPROT[tipI[j]]);
			       break;
			     case DNA_DATA:
			       fprintf(newFile, "%c", inverseMeaningDNA[tipI[j]]);
			       break;
			     default:
			       assert(0);
			     }
			 }
		     }
		 
		  fprintf(newFile, "\n");		
		}	      
	      fclose(newFile);	    
	    }
	  else
	    {
	      if(nonTrivial && !countUndeterminedColumns)
		printf("An alignment file with similar sequences removed has already\n");
	      if(!nonTrivial && countUndeterminedColumns)
		printf("An alignment file with undetermined columns removed has already\n");
	      if(nonTrivial && countUndeterminedColumns)
		printf("An alignment file with undetermined columns and similar sequences removed has already\n");
	      
	      printf("been printed to file %s\n",  noDupFile);
	      
	      if(nonTrivial && !countUndeterminedColumns)
		fprintf(f, "An alignment file with similar sequences removed has already\n");
	      if(!nonTrivial && countUndeterminedColumns)
		fprintf(f, "An alignment file with undetermined columns removed has already\n");
	      if(nonTrivial && countUndeterminedColumns)
		fprintf(f, "An alignment file with undetermined columns and similar sequences removed has already\n");
	      
	      fprintf(f, "been printed to file %s\n",  noDupFile);
	    }
	}
    }


  free(undeterminedList);
  free(omissionList);
  free(modelList);
  if(processID == 0)	      
    fclose(f);
}




static void generateBS(tree *tr, analdef *adef)
{
  int i, j, k, w;
  int count;
  char outName[1024], buf[16];
  FILE *of;

  assert(adef->boot != 0);

  for(i = 0; i < adef->multipleRuns; i++)
    {     
      makeboot(adef, tr);

      count = 0;
      for(j = 0; j < tr->cdta->endsite; j++)	
	count += tr->cdta->aliaswgt[j];    	         
      
      assert(count == tr->rdta->sites);
     
      strcpy(outName, workdir);
      strcat(outName, seq_file);
      strcat(outName, ".BS");
      sprintf(buf, "%d", i);
      strcat(outName, buf);
      printf("Printing replicate %d to %s\n", i, outName);

      of = fopen(outName, "w");

      fprintf(of, "%d %d\n", tr->mxtips, count);  
      
      for(j = 1; j <= tr->mxtips; j++)
	{	
	  char *tip   =  tr->yVector[tr->nodep[j]->number];	   
	  fprintf(of, "%s ", tr->nameList[j]);	
	  
	  for(k = 0; k < tr->cdta->endsite; k++)
	    {
	      switch(tr->dataVector[k])
		{
		case DNA_DATA:
		  for(w = 0; w < tr->cdta->aliaswgt[k]; w++)		  
		    fprintf(of, "%c", inverseMeaningDNA[tip[k]]);	
		  break;
		case AA_DATA:
		  for(w = 0; w < tr->cdta->aliaswgt[k]; w++)
		    fprintf(of, "%c", inverseMeaningPROT[tip[k]]);
		  break;
		default:
		  assert(0);
		}
	    }
	  
	  fprintf(of, "\n");	
	}
      fclose(of);
    }  
}



     

static void splitMultiGene(tree *tr, rawdata *rdta)
{
  int i, l;  
  int n = rdta->sites + 1;
  int *modelFilter = (int *)malloc(sizeof(int) * n);
  int length, k;
  char *tip;
  FILE *outf;
  char outFileName[2048];
  char buf[16];
  
  for(i = 0; i < tr->NumberOfModels; i++)
    {      
      strcpy(outFileName, seq_file);
      sprintf(buf, "%d", i);
      strcat(outFileName, ".GENE.");
      strcat(outFileName, buf);
      outf = fopen(outFileName, "w");
      length = 0;
      for(k = 1; k < n; k++)
	{
	  if(tr->model[k] == i)
	    {
	      modelFilter[k] = 1;
	      length++;
	    }
	  else
	    modelFilter[k] = -1;
	}

      fprintf(outf, "%d %d\n", rdta->numsp, length);
      
      for(l = 1; l <= rdta->numsp; l++)
	{
	  fprintf(outf, "%s ", tr->nameList[l]);

	  tip = &(rdta->y[l][0]);	    

	  for(k = 1; k < n; k++)
	    {
	      if(modelFilter[k] == 1)
		{
		  switch(tr->dataVector[k])
		    {
		    case AA_DATA:
		      fprintf(outf, "%c", inverseMeaningPROT[tip[k]]);
		      break;
		    case DNA_DATA:
		      fprintf(outf, "%c", inverseMeaningDNA[tip[k]]);
		      break;
		    default:
		      assert(0);
		    }		 
		}
	    }
	  fprintf(outf, "\n");

	}
      
      fclose(outf);

      printf("Wrote individual gene/partition alignment to file %s\n", outFileName);
    }  
            
  free(modelFilter);
  printf("Wrote all %d individual gene/partition alignments\n", tr->NumberOfModels);
  printf("Exiting normally\n");
}



void calculateModelOffsets(tree *tr)
{
  int 
    i, 
    patterns,
    currentOffset = 0,
    dnaSpan,
    aaSpan; 
  
  switch(tr->rateHetModel)
    {
    case CAT:
      dnaSpan = 4;
      aaSpan  = 20;
      break;
    case GAMMA:
    case GAMMA_I:
      dnaSpan = 16;
      aaSpan = 80;
      break;
    default:
      assert(0);
    } 

  tr->partitionData[0].modelOffset = currentOffset;
  
  for(i = 1; i < tr->NumberOfModels; i++)
    {     
      patterns =  tr->partitionData[i - 1].upper - tr->partitionData[i - 1].lower;
      switch(tr->partitionData[i - 1].dataType)
	{
	case AA_DATA:
	  currentOffset += aaSpan * patterns;
	  break;
	case DNA_DATA:
	  currentOffset += dnaSpan * patterns;
	  break;
	default:
	  assert(0);
	}

      tr->partitionData[i].modelOffset = currentOffset;      
    } 
}


void allocNodex (tree *tr, analdef *adef)
{
  nodeptr  p;
  int  i;       

  assert(tr->expArray == (int*)NULL);
  assert(tr->likelihoodArray == (double*)NULL);
  assert(tr->sumBuffer == (double *)NULL);

#ifdef _LOCAL_DATA
  tr->currentModel = adef->model;
  masterBarrier(THREAD_ALLOC_LIKELIHOOD, tr);  
#else
  {
    int span;

    tr->expArray = (int *)malloc(tr->cdta->endsite * tr->mxtips * sizeof(int));

    if(tr->mixedData)
      {
	tr->numberOfProteinPositions    = 0;
	tr->numberOfNucleotidePositions = 0;
	for(i = 0; i < tr->cdta->endsite; i++)
	  {
	    switch(tr->dataVector[i])
	      {
	      case AA_DATA:
		tr->numberOfProteinPositions++;
		break;
	      case DNA_DATA:
		tr->numberOfNucleotidePositions++;
		break;
	      default:
		assert(0);
	      }
	  }      
	
	switch(adef->model)
	  {	 
	  case M_PROTCAT:	
	  case M_GTRCAT:
	    span = tr->numberOfNucleotidePositions * 4 + tr->numberOfProteinPositions * 20;
	    tr->likelihoodArray = (double *)malloc(tr->mxtips * span * sizeof(double));	 
	    tr->sumBuffer  = (double *)malloc(span * sizeof(double));
	    break;  	        
	  case M_PROTGAMMA:	 	
	  case M_GTRGAMMA:
	    span = tr->numberOfNucleotidePositions * 16 + tr->numberOfProteinPositions * 80;
	    tr->likelihoodArray = (double *)malloc(tr->mxtips * span * sizeof(double));	
	    tr->sumBuffer  = (double *)malloc(span * sizeof(double));
	    break;	 
	  default:	
	    assert(0);
	  } 
	
	calculateModelOffsets(tr);
	/*printf("DNA %d AA %d\n",  tr->numberOfNucleotidePositions, tr->numberOfProteinPositions);     */
      }
    else
      {
	switch(adef->model)
	  {	 
	  case M_PROTCAT:	 
	    span = 20 * tr->cdta->endsite;
	    tr->likelihoodArray = (double *)malloc(span * tr->mxtips * sizeof(double));
	    tr->sumBuffer  = (double *)malloc(span * sizeof(double));
	    break;  	        
	  case  M_PROTGAMMA:
	    span = 80 * tr->cdta->endsite;
	    tr->likelihoodArray = (double *)malloc(span * tr->mxtips * sizeof(double));
	    tr->sumBuffer  = (double *)malloc(span * sizeof(double));
	    break;	
	  case M_GTRGAMMA:
	    span = 16 * tr->cdta->endsite;
	    tr->likelihoodArray = (double *)malloc(span * tr->mxtips * sizeof(double));
	    tr->sumBuffer  = (double *)malloc(span * sizeof(double));
	    break;
	  case M_GTRCAT:
	    span = 4 * tr->cdta->endsite;
	    tr->likelihoodArray = (double *)malloc(span * tr->mxtips * sizeof(double));
	    tr->sumBuffer  = (double *)malloc(span * sizeof(double));	 
	    break;	 
	  default:	
	    assert(0);
	  }      
      }
    
    for(i = 0; i < tr->mxtips; i++)
      tr->xVector[i] = &(tr->likelihoodArray[i * span]);
  }
#endif

  for (i = tr->mxtips + 1; (i <= 2*(tr->mxtips) - 2); i++) 
    {    
      p = tr->nodep[i];                
      p->x = 1;	          
    }      
}


void freeNodex(tree *tr)
{
  nodeptr  p; 
  int  i;   
       
#ifdef _LOCAL_DATA
  masterBarrier(THREAD_FREE_LIKELIHOOD, tr);
#else
  free(tr->expArray); 
  free(tr->likelihoodArray);
  free(tr->sumBuffer);
  tr->expArray        = (int*)NULL;
  tr->likelihoodArray = (double*)NULL;
  tr->sumBuffer       = (double*)NULL;
#endif
  
  for (i = tr->mxtips + 1; (i <= 2*(tr->mxtips) - 2); i++) 
    {
      p = tr->nodep[i];
      while(!p->x)	
	p = p->next;   
      p->x             = 0;
      p->next->x       = 0;
      p->next->next->x = 0;
    }
}

static void initAdef(analdef *adef)
{
 
  adef->bootstrapBranchLengths = FALSE;
  adef->model                  = M_GTRCAT;  
  adef->max_rearrange          = 21;  
  adef->stepwidth              = 5;
  adef->initial                = adef->bestTrav = 10;
  adef->initialSet             = FALSE;
  adef->restart                = FALSE;
  adef->mode                   = BIG_RAPID_MODE;
  adef->categories             = 25; 
  adef->boot                   = 0;
  adef->rapidBoot              = 0;
  adef->useWeightFile          = FALSE;
  adef->checkpoints            = 0;
  adef->startingTreeOnly       = 0;
  adef->useMixedModel          = 0;
  adef->multipleRuns           = 1;
  adef->useMultipleModel       = FALSE;
  adef->likelihoodEpsilon      = 0.1;
  adef->constraint             = FALSE;
  adef->grouping               = FALSE; 
  adef->randomStartingTree     = FALSE;
  adef->categorizeGamma        = FALSE;
  adef->parsimonySeed          = 0;
  adef->proteinMatrix          = JTT;
  adef->protEmpiricalFreqs     = 0;    
  adef->outgroup               = FALSE;
  adef->useInvariant           = FALSE;
  adef->sequenceSimilarity     = 1.0;
  adef->permuteTreeoptimize    = FALSE;
  adef->useInvariant           = FALSE; 
  adef->allInOne               = FALSE; 
  adef->multiBoot              = 0;
  adef->likelihoodTest         = FALSE;
  adef->reallyThoroughBoot     = FALSE;
  adef->perGeneBranchLengths   = FALSE;
  adef->treeLength             = FALSE;
  adef->computePerSiteLLs      = FALSE;
  adef->generateBS             = FALSE;
  adef->bootStopOnly           = 0;
  adef->bootStopping           = FALSE;
  adef->gapyness               = 0.0;
  adef->similarityFilterMode   = 0;
  adef->bootstopCutoff         = 0.0;
  adef->useExcludeFile         = FALSE;
  adef->userProteinModel       = FALSE;
  adef->externalAAMatrix       = (double*)NULL;
  adef->rapidML_Addition       = FALSE;
  adef->computeELW             = FALSE;
  adef->computeDistance        = FALSE;
#ifdef _VINCENT
  adef->optimizeBSmodel        = TRUE;
#endif
}




static int modelExists(char *model, analdef *adef)
{
  int i;
  char *protModels[10] = {"DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM"};
  char thisModel[1024];
  

  /*********** DNA **********************/

  if(strcmp(model, "GTRGAMMAI\0") == 0)
    {
      adef->model = M_GTRGAMMA;
      adef->useInvariant = TRUE;
      return 1;
    }  

  if(strcmp(model, "GTRGAMMA\0") == 0)
    {
      adef->model = M_GTRGAMMA;
      return 1;
    }
  
  if(strcmp(model, "GTRCAT\0") == 0)
    {
      adef->model = M_GTRCAT;      
      return 1;
    }     
  
  if(strcmp(model, "GTRMIX\0") == 0)
    {
      adef->model = M_GTRCAT;
      adef->useMixedModel = 1;
      return 1;
    }

   if(strcmp(model, "GTRMIXI\0") == 0)
    {
      adef->model = M_GTRCAT;
      adef->useMixedModel = 1;
      adef->useInvariant = TRUE;
      return 1;
    }


  if(strcmp(model, "GTRCAT_GAMMA\0") == 0)
    {
      adef->model = M_GTRCAT;
      adef->useMixedModel = 1;
      adef->categorizeGamma = TRUE;
      return 1;
    }      

  if(strcmp(model, "GTRCAT_GAMMAI\0") == 0)
    {
      adef->model = M_GTRCAT;
      adef->useMixedModel = 1;
      adef->categorizeGamma = TRUE;
      adef->useInvariant = TRUE;
      return 1;
    } 

  /*************** AA GTR ********************/
	
  if(strcmp(model, "PROTCATGTR\0") == 0)
    {
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR;
      return 1;
    }
  if(strcmp(model, "PROTMIXGTR\0") == 0)
    {
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR;
      adef->useMixedModel = 1;
      return 1;
    }
  if(strcmp(model, "PROTGAMMAGTR\0") == 0)
    {
      adef->model = M_PROTGAMMA;
      adef->proteinMatrix = GTR;
      return 1;
    }

   if(strcmp(model, "PROTCAT_GAMMAGTR\0") == 0)
    {
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR;
      adef->useMixedModel = 1;
      adef->categorizeGamma = TRUE;
      return 1;
    }

   if(strcmp(model, "PROTCAT_GAMMAIGTR\0") == 0)
    {
      adef->model = M_PROTCAT;
      adef->proteinMatrix = GTR;
      adef->useMixedModel = 1;
      adef->categorizeGamma = TRUE;
      adef->useInvariant = TRUE;
      return 1;
    }
  
  /****************** AA ************************/
  
  for(i = 0; i < 10; i++)
    {
      /* check CAT */

      strcpy(thisModel, "PROTCAT");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  return 1;
	}

      /* check CATF */

      strcpy(thisModel, "PROTCAT");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  return 1;
	}

      /****************check MIX ************************/

      strcpy(thisModel, "PROTMIX");
      strcat(thisModel, protModels[i]);      
  
      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->useMixedModel = 1;
	  return 1;
	}

      /*check MIXI */

      strcpy(thisModel, "PROTMIXI");
      strcat(thisModel, protModels[i]);      
  
      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->useMixedModel = 1;
	  adef->useInvariant = TRUE;
	  return 1;
	}


      /* check MIXmodelF */

      strcpy(thisModel, "PROTMIX");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->useMixedModel = 1;
	  adef->protEmpiricalFreqs = 1;
	  return 1;
	}
      
      /* check MIXImodelF */

      strcpy(thisModel, "PROTMIXI");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->proteinMatrix = i;
	  adef->useMixedModel = 1;
	  adef->protEmpiricalFreqs = 1;
	   adef->useInvariant = TRUE;
	  return 1;
	}

      /****************check GAMMA ************************/

      strcpy(thisModel, "PROTGAMMA");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  return 1;
	}	

      /*check GAMMAI*/

      strcpy(thisModel, "PROTGAMMAI");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->useInvariant = TRUE;
	  return 1;
	}


      /* check GAMMAmodelF */

      strcpy(thisModel, "PROTGAMMA");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");
     
      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  return 1;
	}	

      /* check GAMMAImodelF */

      strcpy(thisModel, "PROTGAMMAI");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");
      
      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTGAMMA;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  adef->useInvariant = TRUE;
	  return 1;
	}	
      
      /****************check CAT_GAMMA ************************/
      

      strcpy(thisModel, "PROTCAT_GAMMA");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->useMixedModel = 1;
	  adef->categorizeGamma = TRUE;
	  adef->proteinMatrix = i;
	  return 1;
	}	

      /* check CAT_GAMMAI */

      strcpy(thisModel, "PROTCAT_GAMMAI");
      strcat(thisModel, protModels[i]);

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->useMixedModel = 1;
	  adef->categorizeGamma = TRUE;
	  adef->proteinMatrix = i;
	  adef->useInvariant = TRUE;
	  return 1;
	}	

      /* check CAT_GAMMAmodelF */
      
      strcpy(thisModel, "PROTCAT_GAMMA");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->useMixedModel = 1;
	  adef->categorizeGamma = TRUE;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  return 1;
	}	

      /*check CAT_GAMMAImodelF */

      strcpy(thisModel, "PROTCAT_GAMMAI");
      strcat(thisModel, protModels[i]);
      strcat(thisModel, "F");

      if(strcmp(model, thisModel) == 0)
	{
	  adef->model = M_PROTCAT;
	  adef->useMixedModel = 1;
	  adef->categorizeGamma = TRUE;
	  adef->proteinMatrix = i;
	  adef->protEmpiricalFreqs = 1;
	  adef->useInvariant = TRUE;
	  return 1;
	}	

      

    }

  /*********************************************************************************/

  
  
  return 0;
}



static int mygetopt(int argc, char **argv, char *opts, int *optind, char **optarg)
{
  static int sp = 1;
  register int c;
  register char *cp;

  if(sp == 1)
    {
      if(*optind >= argc || argv[*optind][0] != '-' || argv[*optind][1] == '\0')
	return -1;
    }
  else 
    {
      if(strcmp(argv[*optind], "--") == 0) 
	{
	  *optind =  *optind + 1;
	  return -1;
	}
    }

  c = argv[*optind][sp];
  if(c == ':' || (cp=strchr(opts, c)) == 0) 
    {
      printf(": illegal option -- %c \n", c);
      if(argv[*optind][++sp] == '\0') 
	{
	  *optind =  *optind + 1;
	  sp = 1;
	}
      return('?');
    }
  if(*++cp == ':') 
    {
      if(argv[*optind][sp+1] != '\0')
	{
	  *optarg = &argv[*optind][sp+1];
	  *optind =  *optind + 1;
	}
      else 
	{
	  *optind =  *optind + 1;
	  if(*optind >= argc) 
	    {
	      printf(": option requires an argument -- %c\n", c);
	      sp = 1;
	      return('?');
	    } 
	  else
	    {
	      *optarg = argv[*optind];
	      *optind =  *optind + 1;
	    }
	}
      sp = 1;
    } 
  else 
    {
      if(argv[*optind][++sp] == '\0') 
	{
	  sp = 1;
	  *optind =  *optind + 1;
	}
      *optarg = 0;
    }
  return(c);
  }

static void checkOutgroups(tree *tr, analdef *adef)
{
  if(adef->outgroup)
    {
      boolean found;
      int i, j;
      
      if(tr->numberOfOutgroups != 1 && adef->mode ==  MEHRING_ALGO)
	{
	  printf("Error, you must specify exactly one sequence via \"-o\" to \n");
	  printf("to run the sequence position determination algorithm\n");
	  exit(-1);
	}

      for(j = 0; j < tr->numberOfOutgroups; j++)
	{
	  found = FALSE;
	  for(i = 1; (i <= tr->mxtips) && !found; i++)
	    {
	      if(strcmp(tr->nameList[i], tr->outgroups[j]) == 0)
		{
		  tr->outgroupNums[j] = i;
		  found = TRUE;
		}
	    }
	  if(!found)
	    {
	      printf("Error, the outgroup name \"%s\" you specified can not be found in the alignment, exiting ....\n", tr->outgroups[j]);
	      errorExit(-1);
	    }
	}
    }
  
}

static void parseOutgroups(char *outgr, tree *tr)
{
  int count = 1, i, k;
  char name[nmlngth];

  i = 0;
  while(outgr[i] != '\0')
    {
      if(outgr[i] == ',')
	count++;
      i++;
    }

  tr->numberOfOutgroups = count;

  tr->outgroups = (char **)malloc(sizeof(char *) * count);

  for(i = 0; i < tr->numberOfOutgroups; i++)   
    tr->outgroups[i] = (char *)malloc(sizeof(char) * nmlngth);    

  tr->outgroupNums = (int *)malloc(sizeof(int) * count);
    
  i = 0;
  k = 0;
  count = 0;
  while(outgr[i] != '\0')
    {
      if(outgr[i] == ',')
	{	
	  name[k] = '\0';
	  strcpy(tr->outgroups[count], name);
	  count++;
	  k = 0;	 
	}
      else
	{
	  name[k] = outgr[i];
	  k++;
	}
      i++;
    }

  name[k] = '\0';
  strcpy(tr->outgroups[count], name);

  /*for(i = 0; i < tr->numberOfOutgroups; i++)
    printf("%d %s \n", i, tr->outgroups[i]);*/


  /*printf("%s \n", name);*/
}


/*********************************** OUTGROUP STUFF END *********************************************************/
static void printVersionInfo(void)
{
  printf("\nThis is %s version %s released by Alexandros Stamatakis in %s\n\n",  programName, programVersion, programDate);
}

static void printMinusFUsage(void)
{
  printf("\n");
  printf("              \"-f a\": rapid Bootstrap analysis and search for best-scoring ML tree in one program run\n");
  printf("              \"-f b\": draw bipartition information on a tree provided with \"-t\" based on multiple trees\n");
  printf("                      (e.g. form a bootstrap) in a file specifed by \"-z\"\n");
  printf("              \"-f c\": check if the alignment can be properly read by RAxML\n");
  printf("              \"-f d\": new rapid hill-climbing \n");
  printf("              \"-f e\": optimize model+branch lengths for given input tree under GAMMA/GAMMAI only\n"); 

  /*printf("              \"-f f\": optimize individual per-site evolutionary rates 
    on a fixed input tree and compute sliding window tree lengths\n");*/

  printf("              \"-f g\": compute per site log Likelihoods for one ore more trees passed via\n");
  printf("                      \"-z\" and write them to a file that can be read by CONSEL\n");
  printf("              \"-f h\": compute log likelihood test (SH-test) between best tree passed via \"-t\"\n");
  printf("                      and a bunch of other trees passed via \"-z\" \n");
  printf("              \"-f i\": perform a really thorough bootstrap, refinement of final BS tree under GAMMA and a\n");
  printf("                      more exhaustive algorithm\n");
  printf("              \"-f j\": generate a bunch of bootstrapped alignment files from an original alignemnt file\n");

  /* printf("              \"-f k\": \n"); */ 

  printf("              \"-f m\": Compare bipartitions between two bunches of trees passed via \"-t\" and \"-z\" \n");
  printf("                      respectively. This will return the Pearson correlation between all bipartitions found\n");
  printf("                      in the two tree files. A file called RAxML_bipartitionFrequencies.outpuFileName\n");
  printf("                      will be printed that contains the pair-wise bipartition frequencies of the two sets\n");
  printf("              \"-f n\": Compute the log likelihood score of all trees contained in a tree file provided by\n");
  printf("                      \"-z\" under GAMMA or GAMMA+P-Invar\n");
  printf("              \"-f o\": old and slower rapid hill-climbing \n");
  printf("              \"-f p\": perform pure stepwise MP addition of new sequences to an incomplete starting tree\n");
   
  /* printf("              \"-f r\": optimize individual per-site evolutionary rates on a fixed input tree\n");*/

  printf("              \"-f s\": split up a multi-gene partitioned alignment into the respective subalignments \n");
  printf("              \"-f t\": do randomized tree searches on one fixed starting tree\n");
  printf("              \"-f w\": compute ELW test on a bunch of trees passed via \"-z\" \n");
  printf("              \"-f x\": compute pair-wise ML distances, ML model parameters will be estimated on an MP \n");
  printf("                      starting tree or a user-defined tree passed via \"-t\", only allowed for GAMMA-based\n");
  printf("                      models of rate heterogeneity\n");
  printf("\n"); 
  printf("              DEFAULT: new rapid hill climbing\n");  
  printf("\n");
}


static void printREADME(void)
{
  printVersionInfo();  
  printf("\n");
  printf("Please also consult the RAxML-manual\n");
  printf("To report bugs send an email to Alexandros.Stamatakis@epfl.ch\n\n\n");
  
  printf("raxmlHPC[-MPI|-PTHREADS] -s sequenceFileName -n outputFileName -m substitutionModel\n");
  printf("                         [-a weightFileName] [-b bootstrapRandomNumberSeed] [-c numberOfCategories]\n");
  printf("                         [-d] [-e likelihoodEpsilon] [-E excludeFileName] [-f a|b|c|d|e|g|h|i|j|m|n|o|p|s|t|w|x]\n");
  printf("                         [-g groupingFileName] [-h] [-i initialRearrangementSetting] [-j] [-k] \n");
  printf("                         [-l sequenceSimilarityThreshold] [-L sequenceSimilarityThreshold] [-M]\n"); 
  printf("                         [-o outGroupName1[,outGroupName2[,...]]] [-p parsimonyRandomSeed] [-P proteinModel]\n");
  printf("                         [-q multipleModelFileName] [-r binaryConstraintTree] [-t userStartingTree]\n"); 
  printf("                         [-T numberOfThreads] [-u multiBootstrapSearches] [-v][-w workingDirectory]\n");
  printf("                         [-x rapidBootstrapRandomNumberSeed][-y][-z multipleTreesFile] [-#|-N numberOfRuns]\n");
  printf("\n");
  printf("      -a      Specify a column weight file name to assign individual weights to each column of \n");
  printf("              the alignment. Those weights must be integers separated by any type and number \n");
  printf("              of whitespaces whithin a separate file, see file \"example_weights\" for an example.\n");
  printf("\n");
  printf("      -b      Specify an integer number (random seed) and turn on bootstrapping\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -c      Specify number of distinct rate catgories for RAxML when modelOfEvolution\n");
  printf("              is set to GTRCAT or GTRMIX\n");
  printf("              Individual per-site rates are categorized into numberOfCategories rate \n");
  printf("              categories to accelerate computations. \n");
  printf("\n");
  printf("              DEFAULT: 25\n");
  printf("\n");
  printf("      -d      start ML optimization from random starting tree \n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -e      set model optimization precision in log likelihood units for final\n"); 
  printf("              optimization of tree topology under MIX/MIXI or GAMMA/GAMMAI\n");
  printf("\n");
  printf("              DEFAULT: 0.1   for models not using proportion of invariant sites estimate\n");
  printf("                       0.001 for models using proportion of invariant sites estimate\n");
  printf("\n");
  printf("      -E      specify an exclude file name, that contains a specification of alignment positions you wish to exclude.\n");
  printf("              Format is similar to Nexus, the file shall contain entries like \"100-200 300-400\", to exclude a\n");
  printf("              single column write, e.g., \"100-100\", if you use a mixed model, an appropriatly adapted model file\n");
  printf("              will be written.\n"); 
  printf("\n");
  printf("      -f      select algorithm:\n");

  printMinusFUsage();
 
  printf("\n");
  printf("      -g      specify the file name of a multifurcating constraint tree\n");
  printf("              this tree does not need to be comprehensive, i.e. must not contain all taxa\n");
  printf("\n");
  printf("      -h      Display this help message.\n");  
  printf("\n");
  printf("      -i      Initial rearrangement setting for the subsequent application of topological \n");
  printf("              changes phase\n");
  printf("\n");
  printf("              DEFAULT: determined by program\n");
  printf("\n");
  printf("      -j      Specifies if checkpoints will be written by the program. If checkpoints \n");
  printf("              (intermediate tree topologies) shall be written by the program specify \"-j\"\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");
  printf("      -k      Specifies that bootstrapped trees should be printed with branch lengths.\n");
  printf("              The bootstraps will run a bit longer, because model parameters will be optimized\n");
  printf("              at the end of each run. Use with CATMIX/PROTMIX or GAMMA/GAMMAI.\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");   
  printf("\n");
  printf("      -l      Specify a threshold for sequence similarity clustering. RAxML will then print out an alignment\n");
  printf("              to a file called sequenceFileName.reducedBy.threshold that only contains sequences <= the\n");
  printf("              specified thresold that must be between  0.0 and 1.0. RAxML uses the QT-clustering algorithm \n");
  printf("              to perform this task. In addition, a file called RAxML_reducedList.outputFileName will be written\n");
  printf("              that contains clustering information.\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");   
  printf("\n");
  printf("      -L      Same functionality as \"-l\" above, but uses a less exhasutive and thus faster clustering algorithm\n");
  printf("              This is intended for very large datasets with more than 20,000-30,000 sequences\n");
  printf("\n");
  printf("              DEFAULT: OFF\n"); 
  printf("\n");
  printf("      -m      Model of Nucleotide or Amino Acid Substitution: \n");  
  printf("\n");
  printf("              NUCLEOTIDES:\n\n");
  printf("                \"-m GTRCAT\"        : GTR + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                     evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                     rate categories for greater computational efficiency\n");
  printf("                                     if you do a multiple analysis with  \"-#\" or \"-N\" but without bootstrapping the program\n");
  printf("                                     will use GTRMIX instead\n");
  printf("                \"-m GTRGAMMA\"      : GTR + Optimization of substitution rates + GAMMA model of rate \n");
  printf("                                     heterogeneity (alpha parameter will be estimated)\n");
  printf("                \"-m GTRMIX\"        : Inference of the tree under GTRCAT\n");
  printf("                                     and thereafter evaluation of the final tree topology under GTRGAMMA\n");
  printf("                \"-m GTRCAT_GAMMA\"  : Inference of the tree with site-specific evolutionary rates.\n");
  printf("                                     However, here rates are categorized using the 4 discrete GAMMA rates.\n");
  printf("                                     Evaluation of the final tree topology under GTRGAMMA\n");
  printf("                \"-m GTRGAMMAI\"     : Same as GTRGAMMA, but with estimate of proportion of invariable sites \n");
  printf("                \"-m GTRMIXI\"       : Same as GTRMIX, but with estimate of proportion of invariable sites \n");
  printf("                \"-m GTRCAT_GAMMAI\" : Same as GTRCAT_GAMMA, but with estimate of proportion of invariable sites \n");
  printf("\n");
  printf("              AMINO ACIDS:\n\n");          
  printf("                \"-m PROTCATmatrixName[F]\"        : specified AA matrix + Optimization of substitution rates + Optimization of site-specific\n");
  printf("                                                   evolutionary rates which are categorized into numberOfCategories distinct \n");
  printf("                                                   rate categories for greater computational efficiency\n");
  printf("                                                   if you do a multiple analysis with  \"-#\" or \"-N\" but without bootstrapping the program\n");
  printf("                                                   will use PROTMIX... instead\n");
  printf("                \"-m PROTGAMMAmatrixName[F]\"      : specified AA matrix + Optimization of substitution rates + GAMMA model of rate \n");
  printf("                                                   heterogeneity (alpha parameter will be estimated)\n");
  printf("                \"-m PROTMIXmatrixName[F]\"        : Inference of the tree under specified AA matrix + CAT\n");
  printf("                                                   and thereafter evaluation of the final tree topology under specified AA matrix + GAMMA\n");
  printf("                \"-m PROTCAT_GAMMAmatrixName[F]\"  : Inference of the tree under specified AA matrix and site-specific evolutionary rates.\n");
  printf("                                                   However, here rates are categorized using the 4 discrete GAMMA rates.\n"); 
  printf("                                                   Evaluation of the final tree topology under specified AA matrix + GAMMA\n");
  printf("                \"-m PROTGAMMAImatrixName[F]\"     : Same as PROTGAMMAmatrixName[F], but with estimate of proportion of invariable sites \n");
  printf("                \"-m PROTMIXImatrixName[F]\"       : Same as PROTMIXmatrixName[F], but with estimate of proportion of invariable sites \n");
  printf("                \"-m PROTCAT_GAMMAImatrixName[F]\" : Same as PROTCAT_GAMMAmatrixName[F], but with estimate of proportion of invariable sites \n");
  printf("\n");
  printf("                Available AA substitution models: DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, GTR\n");
  printf("                With the optional \"F\" appendix you can specify if you want to use empirical base frequencies\n");
  printf("                Please not that for mixed models you can in addition specify the per-gene AA model in\n");
  printf("                the mixed model file (see manual for details)\n");
  printf("\n");
  printf("      -M      Switch on estimation of individual per-partition branch lengths. Only has effect when used in combination with \"-q\"\n");  
  printf("              Branch lengths for individual partitions will be printed to separate files\n");
  printf("              A weighted average of the branch lengths is computed by using the respective partition lengths\n");
  printf("\n"),
  printf("              DEFAULT: OFF\n");  
  printf("\n");
  printf("      -n      Specifies the name of the output file.\n"); 
  printf("\n"); 
  printf("      -o      Specify the name of a single outgrpoup or a comma-separated list of outgroups, eg \"-o Rat\" \n");
  printf("              or \"-o Rat,Mouse\", in case that multiple outgroups are not monophyletic the first name \n");
  printf("              in the list will be selected as outgroup, don't leave spaces between taxon names!\n");
  printf("\n");
  printf("      -q      Specify the file name which contains the assignment of models to alignment\n"); 
  printf("              partitions for multiple models of substitution. For the syntax of this file\n");
  printf("              please consult the manual.\n"); 
  printf("\n"); 
  printf("      -p      Specify a random number seed for the parsimony inferences. This allows you to reproduce your results\n");
  printf("              and will help me debug the program. This option HAS NO EFFECT in the parallel MPI version\n");
  printf("\n");
  printf("      -P      Specify the file name of a user-defined AA (Protein) substitution model. This file must contain\n");
  printf("              420 entries, the first 400 being the AA substitution rates (this must be a symmetric matrix) and the\n");
  printf("              last 20 are the empirical base frequencies\n");
  printf("\n");
  printf("      -r      Specify the file name of a binary constraint tree.\n");
  printf("              this tree does not need to be comprehensive, i.e. must not contain all taxa\n");
  printf("\n");  
  printf("      -s      Specify the name of the alignment data file in PHYLIP format\n");
  printf("\n");
  printf("      -t      Specify a user starting tree file name in Newick format\n");
  printf("\n");
  printf("      -T      PTHREADS VERSION ONLY! Specify the number of threads you want to run.\n");
  printf("              Make sure to set \"-T\" to at most the number of CPUs you have on your machine,\n");
  printf("              otherwise, there will be a huge performance decrease!\n");
  printf("\n");
  printf("      -u      Specify the number of multiple BS searches per replicate\n");
  printf("              to obtain better ML trees for each replicate\n");
  printf("\n");
  printf("              DEFAULT: One ML search per BS replicate\n");
  printf("\n");
  printf("      -v      Display version information\n");  
  printf("\n");
  printf("      -w      Name of the working directory where RAxML will write its output files\n");
  printf("\n");
  printf("              DEFAULT: current directory\n"); 
  printf("\n");
  printf("      -x      Specify an integer number (random seed) and turn on rapid bootstrapping\n");  
  printf("\n");
  printf("      -y      If you want to only compute a parsimony starting tree with RAxML specify \"-y\",\n");
  printf("              the program will exit after computation of the starting tree\n");
  printf("\n");
  printf("              DEFAULT: OFF\n");
  printf("\n");  
  printf("      -z      Specify the file name of a file containing multiple trees e.g. from a bootstrap\n");
  printf("              that shall be used to draw bipartition values onto a tree provided with \"-t\",\n");
  printf("              It can also be used to compute per site log likelihoods in combination with \"-f g\"\n");
  printf("              and to read a bunch of trees for a couple of other options (\"-f h\", \"-f m\", \"-f n\").\n");
  printf("\n"); 
  printf("      -#|-N   Specify the number of alternative runs on distinct starting trees\n");
  printf("              In combination with the \"-b\" option, this will invoke a multiple boostrap analysis\n");
  printf("              Note that \"-N\" has been added as an alternative since \"-#\" sometimes caused problems\n");
  printf("              with certain MPI job submission systems, since \"-#\" is often used to start comments\n");
  printf("\n");
  printf("              DEFAULT: 1 single analysis\n");
  printf("\n\n\n\n");

}


static void get_args(int argc, char *argv[], analdef *adef, tree *tr)
{
  int	optind = 1;
  int        c;
  boolean    bad_opt=FALSE;
  char       aut[256];
  char       buf[2048];
  char       *optarg;
  char       model[2048] = ""; 
  char       modelChar;
  double likelihoodEpsilon, sequenceSimilarity;
  int nameSet = 0, 
    alignmentSet = 0,    
    multipleRuns = 0,    
    constraintSet = 0,
    treeSet = 0,
    groupSet = 0,
    modelSet = 0,
    treesSet  = 0,
    multipleBoots = 0;
  long parsimonySeed = 0;
  run_id[0] = 0;
  workdir[0] = 0;
  seq_file[0] = 0;
  tree_file[0] = 0;
  model[0] = 0;
  weightFileName[0] = 0;
  modelFileName[0] = 0;

  /*********** tr inits **************/

#ifdef _USE_PTHREADS
  NumberOfThreads = 0;
#endif
  
  tr->doCutoff = TRUE;
 
  /********* tr inits end*************/

#ifdef _VINCENT
   while(!bad_opt && 
	 ((c = mygetopt(argc,argv,"T:E:N:u:l:x:X:z:g:r:e:a:b:c:f:i:m:t:w:s:n:o:L:B:q:#:p:vdyjhkM", &optind, &optarg))!=-1))
#else
  while(!bad_opt && 
	((c = mygetopt(argc,argv,"T:E:N:u:l:x:z:g:r:e:a:b:c:f:i:m:t:w:s:n:o:L:B:P:q:#:p:vdyjhkM", &optind, &optarg))!=-1))
#endif
    {
    switch(c) 
      { 
      case 'T':		
#ifdef _USE_PTHREADS
	sscanf(optarg,"%d", &NumberOfThreads);	
#else
	if(processID == 0)
	  {
	    printf("Option -T does not have any effect with the sequential or parallel MPI version.\n");
	    printf("It is used to specify the number of threads for the Pthreads-based parallelization\n");
	  }	 
#endif        
	break;
      case 'P':   
	strcpy(proteinModelFileName, optarg);
	adef->userProteinModel = TRUE;
	parseProteinModel(adef);
	break;
      case 'E':
	strcpy(excludeFileName, optarg);
	adef->useExcludeFile = TRUE;
	break;
      case 'M':
	adef->perGeneBranchLengths = TRUE;
	break;
      case 'u':
	sscanf(optarg,"%d", &multipleBoots);
	adef->multiBoot = multipleBoots;
	break;
      case 'o':
	{
	  char *outgroups;
	  outgroups = (char*)malloc(sizeof(char) * (strlen(optarg) + 1));
	  strcpy(outgroups, optarg);
	  parseOutgroups(outgroups, tr);
	  free(outgroups);
	  adef->outgroup = TRUE;
	}
	break;	
      case 'k':
	adef->bootstrapBranchLengths = TRUE;
	break;
      case 'z':	
	strcpy(bootStrapFile, optarg);
	treesSet = 1;
	break;               
      case 'd':
	adef->randomStartingTree = TRUE;
	break;
      case 'g':
	strcpy(tree_file, optarg);
	adef->grouping = TRUE;
	adef->restart  = TRUE;      
	groupSet = 1;
	break;
      case 'r':	
	strcpy(tree_file, optarg);
	adef->restart = TRUE;
	adef->constraint = TRUE;
	constraintSet = 1;
	break;
      case 'e':      
	sscanf(optarg,"%lf", &likelihoodEpsilon);
	adef->likelihoodEpsilon = likelihoodEpsilon;      
	break;
      case 'q':
	strcpy(modelFileName,optarg);
	adef->useMultipleModel = TRUE;      
        break;       
      case 'p':
	sscanf(optarg,"%ld", &parsimonySeed);
	adef->parsimonySeed = parsimonySeed;
	break;
      case 'N':
      case '#':
	/* TODO include auto in readme */
	if(sscanf(optarg,"%d", &multipleRuns) > 0)
	  {
	    adef->multipleRuns = multipleRuns;
	  }
	else
	  {
	    if((sscanf(optarg,"%s", aut) > 0) && 
	       (
		(strcmp(aut, "auto") == 0) ||
		(strcmp(aut, "Auto") == 0) ||
		(strcmp(aut, "AUTO") == 0) ||
		(strcmp(aut, "automatic") == 0) ||
		(strcmp(aut, "Automatic") == 0) ||
		(strcmp(aut, "AUTOMATIC") == 0)
		)
	       )
	      {
		adef->bootStopping = TRUE;
		adef->multipleRuns = 1000;	       	      
	      }
	    else
	      {
		if(processID == 0)
		  {
		    printf("Use -# or -N option either with an integer, e.g., -# 100 or with -# auto\n");
		    printf("or -N 100 or  -N auto respectively, note that auto will not work for the\n");
		    printf("MPI-based parallel version\n");
		  }		
		errorExit(0);
	      }
	  }
	break;
      case 'v':
	printVersionInfo();
	errorExit(0);
      case 'y':
	adef->startingTreeOnly = 1;
	break;
      case 'h':
	printREADME();
	errorExit(0);
      case 'j':	
	adef->checkpoints = 1;
	break;
      case 'a':
	strcpy(weightFileName,optarg);
	adef->useWeightFile = TRUE;
        break;            
      case 'b':
	sscanf(optarg,"%ld", &adef->boot);             
	break;     
      case 'x':       
	sscanf(optarg,"%ld", &adef->rapidBoot);   
#ifdef _VINCENT
	adef->optimizeBSmodel        = FALSE;
#endif   
	break;
#ifdef _VINCENT
      case 'X':       
	sscanf(optarg,"%ld", &adef->rapidBoot);   
	adef->optimizeBSmodel        = TRUE;
	break;
#endif
      case 'c':
	sscanf(optarg, "%d", &adef->categories);       
	break;	  
      case 'l':
	sscanf(optarg,"%lf", &sequenceSimilarity);
	adef->sequenceSimilarity = sequenceSimilarity; 
	adef->mode = SEQUENCE_SIMILARITY_FILTER;
	adef->similarityFilterMode = SMALL_DATA;	
	break;
      case 'L':
	sscanf(optarg,"%lf", &sequenceSimilarity);
	adef->sequenceSimilarity = sequenceSimilarity; 
	adef->mode = SEQUENCE_SIMILARITY_FILTER;
	adef->similarityFilterMode = LARGE_DATA;		
	break;
      case 'B':
	/* TODO include in readme */
	sscanf(optarg,"%lf", &(adef->bootstopCutoff));
	if(adef->bootstopCutoff <= 0.0)
	  {
	    printf("ERROR BootstopCutoff was set to %f, but must be greater than 0.0\n", 
		   adef->bootstopCutoff);
	    errorExit(-1);
	  }
	if(adef->bootstopCutoff == 0.5)
	  {
	    printf("\n\nWARNING: BootstopCutoff was set to %f, this is equivalent to default\n", adef->bootstopCutoff);
	    printf("Bootstopping without the \"-B\" option. Are you sure that this is \n");
	    printf("what you want to do?\n\n");
	  }
	if(adef->bootstopCutoff > 0.5)
	  {
	    printf("ERROR BootstopCutoff was set to %f, but must be smaller or equal to 0.5\n", 
		   adef->bootstopCutoff);
	    errorExit(-1);
	  }
	break;
      case 'f': 
	sscanf(optarg, "%c", &modelChar);
	switch(modelChar)
	  {
	  case 'a':
	    adef->allInOne = TRUE;
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    break;
	  case 'b': 
	    adef->mode = CALC_BIPARTITIONS; 
	    break;
	  case 'c':
	    adef->mode = CHECK_ALIGNMENT;
	    break;
	  case 'd': 
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    break;
	  case 'e': 
	    adef->mode = TREE_EVALUATION; 
	    break;
	  case 'f':
	    /* TODO include in readme */
	    adef->mode       = OPTIMIZE_RATES;
	    adef->treeLength = TRUE;
	    break;
	  case 'g':
	    adef->mode              = OPTIMIZE_RATES;
	    adef->computePerSiteLLs = TRUE;
	    break;
	  case 'h':
	    adef->mode = TREE_EVALUATION;
	    adef->likelihoodTest = TRUE;
	    break;
	  case 'i':
	    adef->reallyThoroughBoot = TRUE;
	    break;
	  case 'j':
	    adef->generateBS = TRUE;
	    break;
	  case 'k':
	    /* TODO include in readme */
	    adef->bootStopOnly = 1;
	    break;
	  case 'l':
	    /* TODO include in readme */
	    adef->bootStopOnly = 2;
	    break;
	  case 'm':
	    adef->bootStopOnly = 3;
	    break;
	  case 'n':
	    adef->bootStopOnly = 4;
	    break;
	  case 'o': 
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = FALSE;
	    break;
	  case 'p':
	    adef->mode =  PARSIMONY_ADDITION;
	    break;	 
	  case 'q':
	    /* TODO include in README */
	    adef->mode = MEHRING_ALGO;
	    break;
	  case 'r':
	    adef->mode =  OPTIMIZE_RATES;
	    break; 
	  case 's':
	    adef->mode = SPLIT_MULTI_GENE;
	    break;       	  	 
	  case 't':
	    adef->mode = BIG_RAPID_MODE;
	    tr->doCutoff = TRUE;
	    adef->permuteTreeoptimize = TRUE;	    
	    break;
	  case 'u':
	    /* TODO readme */
	    adef->mode = ARNDT_MODE;	 
	    break;
	  case 'v':
	    /* TODO README */
	    adef->rapidML_Addition = TRUE;
	    break;
	  case 'w':	  
	    adef->computeELW = TRUE;
	    break;
	  case 'x':
	    adef->mode = DISTANCE_MODE;
	    adef->computeDistance = TRUE;	    
	    break;
	  default: 
	    {
	      if(processID == 0)
		{
		  printf("Error select one of the following algorithms via -f :\n");
		  printMinusFUsage();		 		 
		}	     
	      errorExit(-1);
	    }
	  }
	break;      
      case 'i':
	sscanf(optarg, "%d", &adef->initial);
	adef->initialSet = TRUE;
	break;     
      case 'n':
        strcpy(run_id,optarg);
	nameSet = 1;
        break;
      case 'w':
        strcpy(workdir,optarg);
        break;                 
      case 't':
	strcpy(tree_file, optarg);
	adef->restart = TRUE;
	treeSet = 1;
	break;
      case 's':
	strcpy(seq_file, optarg);
	alignmentSet = 1;
	break;
      case 'm':
	strcpy(model,optarg);
	if(modelExists(model, adef) == 0)
	  {
	    if(processID == 0)
	      {
		printf("Model %s does not exist\n\n", model);
		printf("For DNA data use: GTRCAT        or GTRGAMMA      or\n");
		printf("                  GTRMIX        or GTRMIXI       or\n");
		printf("                  GTRGAMMAI     or GTRCAT_GAMMAI or\n");
		printf("                  GTRCAT_GAMMA\n\n");
		printf("For AA data use:  PROTCATmatrixName[F]        or PROTGAMMAmatrixName[F]      or\n");
		printf("                  PROTMIXmatrixName[F]        or PROTMIXImatrixName[F]       or\n");
		printf("                  PROTGAMMAImatrixName[F]     or PROTCAT_GAMMAImatrixName[F] or\n");
		printf("                  PROTCAT_GAMMAImatrixName[F]\n\n");                
		printf("The AA substitution matrix can be one of the following: \n");
		printf("DAYHOFF, DCMUT, JTT, MTREV, WAG, RTREV, CPREV, VT, BLOSUM62, MTMAM, GTR\n\n");
		printf("With the optional \"F\" appendix you can specify if you want to use empirical base frequencies\n");
		printf("Please note that for mixed models you can in addition specify the per-gene model in\n");
		printf("the mixed model file (see manual for details)\n");	    
	      }
	    errorExit(-1);
	  }      
	else	  
	  modelSet = 1;
	break;      
      default:	   
	errorExit(-1);              
    }    
  }     

#ifdef _USE_PTHREADS  
  if(NumberOfThreads < 2)
    {
      printf("\nThe number of threads is currently set to %d\n", NumberOfThreads);
      printf("Specify the number of threads to run via -T numberOfThreads\n");
      printf("NumberOfThreads must be set to an integer value greater than 1\n\n");
      errorExit(-1);
    }
#endif 




  if(adef->computeELW)
    {
      if(processID == 0)
	{
	  if(adef->boot == 0)
	    {
	      printf("Error, you must specify a bootstrap seed via \"-b\" to compute ELW statistics\n");
	      errorExit(-1);
	    }

	  if(adef->multipleRuns < 2)
	    {
	      printf("Error, you must specify the number of BS replicates via \"-#\" or \"-N\" to compute ELW statistics\n");
	      printf("it should be larger than 1, recommended setting is 100\n");
	      errorExit(-1);
	    }
	  
	  if(!treesSet)
	    {
	      printf("Error, you must specify an input file containing several candidate trees\n");
	      printf("via \"-z\" to compute ELW statistics.\n");
	      errorExit(-1);
	    }

	  if(!(adef->model == M_PROTGAMMA || adef->model == M_GTRGAMMA))
	    {
	      printf("Error ELW test can only be conducted undetr GAMMA or GAMMA+P-Invar models\n");
	      errorExit(-1);
	    }
	}
    }

  if(adef->mode == MEHRING_ALGO && !(adef->restart && adef->outgroup))
    {
      if(processID == 0)
	{
	  printf("\nTo use the sequence position determination algorithm you have to specify a starting tree with \"-t\" \n");
	  printf("and a taxon to be re-inserted with \"-o\" \n");
	  errorExit(-1);
	}      
    }
      


  if(((!adef->boot) && (!adef->rapidBoot)) && adef->bootStopping)
    {
      if(processID == 0)
	{
	  printf("Can't use automatic bootstopping without actually doing a Bootstrap\n");
	  printf("Specify either -x randomNumberSeed (rapid) or -b randomNumberSeed (standard)\n");
	  errorExit(-1);
	}      
    }

  if(adef->boot && adef->rapidBoot)
    {
      if(processID == 0)
	{
	  printf("Can't use standard and rapid BOOTSTRAP simultaneously\n");
	  errorExit(-1);
	}
    }

  if(adef->rapidBoot && !(adef->mode == MEHRING_ALGO))
    {
      if(processID == 0 && (adef->restart || treesSet) && !(groupSet || constraintSet))
	{
	  printf("Error, starting tree(s) will be ignored by rapid Bootstrapping\n");
	  errorExit(-1);
	}

      /*if(processID == 0 && (groupSet || constraintSet))
	{
	  printf("Error, constraint tree will be ignored by rapid Bootstrapping\n");
	  errorExit(-1);
	}            
      */
    }

  if(adef->allInOne && (adef->rapidBoot == 0))
    {
      if(processID == 0)
	{
	  printf("Error, to carry out an ML search after a rapid BS inference you must specify a random number seed with -x\n");
	  errorExit(-1);
	}
    }

  if(adef->mode == SEQUENCE_SIMILARITY_FILTER)
    {
      if(processID == 0)
	{
	  if(adef->sequenceSimilarity <= 0.0 || adef->sequenceSimilarity >= 1.0)
	    {
	      printf("\n ERROR: sequence similarity must be > 0.0 and < 1.0, exiting ...\n");
	      errorExit(-1);
	    }
	}
    }

  if(adef->mode == OPTIMIZE_RATES)
    {     
      if(adef->treeLength && !(adef->model == M_GTRGAMMA || adef->model == M_PROTGAMMA))
	{
	  if(processID == 0)	    
	    printf("\n ERROR: Tree-Length-based sliding window approach only allowed under GAMMA model of rate heterogeneity!\n");
	  errorExit(-1);	    
	}
      
      if(adef->computePerSiteLLs)
	{
	  if(!(adef->model == M_GTRGAMMA || adef->model == M_PROTGAMMA))
	    {
	      if(processID == 0)		
		printf("\n ERROR: Computation of per-site log LHs is only allowed under GAMMA model of rate heterogeneity!\n");
	      errorExit(-1);	    
	    }
	  if(!treesSet)
	    {
	      if(processID == 0)		
		printf("\n ERROR: For Computation of per-site log LHs you need to specify several input trees with \"-z\"\n");
	      errorExit(-1);	  
	    }
	}
      
      if(!adef->restart)
	{
	  if(!adef->computePerSiteLLs && processID == 0)
	    {	      
	      if(!adef->treeLength)	       
		printf("\n You need to specify an input tree with \"-t\" to optimize rates using \"-f r\"\n");
	      else
		printf("\n You need to specify an input tree with \"-t\" to optimize rates and compute the sliding window tree length using \"-f f\"\n");		      	    	      errorExit(-1);
	    }
	}
    } 

  if(adef->bootstrapBranchLengths && (adef->model == M_GTRCAT || adef->model == M_PROTCAT) && (!adef->useMixedModel))
    {
      if(processID == 0)
	{
	  printf("\nWARNING: you want to print out the branch lengths of your bootstrapped trees\n");	
	  printf("WARNING: However you have currently chosen one of the CAT models where the branch lengths\n");
	  printf("WARNING: are essentially meaningless, you should better use CATMIX/PROTMIX instead\n");
	}    
    }  
    
  if(adef->mode == SPLIT_MULTI_GENE && (!adef->useMultipleModel))
    {
      if(processID == 0)
	{
	  printf("\n  Error, you are trying to split a multi-gene alignment into individual genes with the \"-f s\" option\n");	
	  printf("Without specifying a multiple model file with \"-q modelFileName\" \n");
	}
      errorExit(-1);
    }

  if(adef->mode == CALC_BIPARTITIONS && !treesSet)
    {
      if(processID == 0)
	printf("\n  Error, in bipartition computation mode you must specify a file containing multiple trees with the \"-z\" option\n");
      errorExit(-1);
    }

  if(adef->mode == CALC_BIPARTITIONS && !adef->restart)
    {
      if(processID == 0)
	printf("\n  Error, in bipartition computation mode you must specify a tree on which bipartition information will be drawn with the \"-t\" option\n");
      errorExit(-1);
    }

  if(!modelSet)
    { 
      if(processID == 0)
	printf("\n Error, you must specify a model of substitution with the \"-m\" option\n");
      errorExit(-1);
    }
      
  if(adef->computeDistance)
    {
      if(adef->model == M_PROTCAT || adef->model == M_GTRCAT)
	{
	  if(processID == 0)
	    printf("\n Error pairwise distance computation only allowed for GAMMA-based models of rate heterogeneity\n");
	  errorExit(-1);
	}

      if(adef->restart)
	{
	  if(adef->randomStartingTree)
	    {
	      if(processID == 0)
		printf("\n Error pairwise distance computation not allowed for random starting trees\n");
	      errorExit(-1);
	    }
	  	       
	  if(adef->constraint)
	    {
	      if(processID == 0)
		printf("\n Error pairwise distance computation not allowed for binary backbone  constraint tree\n");
	      errorExit(-1);
	    }
	 
	  if(adef->grouping)
	    {
	      if(processID == 0)
		printf("\n Error pairwise distance computation not allowed for constraint tree\n");
	      errorExit(-1);
	    }
	      
	}

      if(adef->boot || adef->rapidBoot)
	{
	  if(processID == 0)
	    printf("\n Bootstrapping not implemented for pairwise distance computation\n");
	  errorExit(-1);
	}
    }

  


  if(adef->useMultipleModel && (adef->model == M_PROTGAMMA || adef->model == M_PROTCAT) && (adef->proteinMatrix == GTR))
    {
      if(processID == 0)
	printf("\n Error GTR model of AA substiution in combination with mixed models is currently not implemented\n");
      errorExit(-1);
    }  

 

  if(!adef->restart && adef->mode == PARSIMONY_ADDITION)
    {
       if(processID == 0)
	 {
	   printf("\n You need to specify an incomplete binary input tree with \"-t\" to execute \n");
	   printf(" RAxML MP stepwise addition with \"-f p\"\n");
	 }
      errorExit(-1);
    }

  

  if(adef->restart && adef->randomStartingTree)
    {
      if(processID == 0)
	{
	  if(adef->constraint)
	    {
	      printf("\n Error you specified a binary constraint tree with -r AND the computation\n");
	      printf("of a random starting tree with -d for the same run\n");
	    }
	  else
	    {
	      if(adef->grouping)
		{
		  printf("\n Error you specified a multifurcating constraint tree with -g AND the computation\n");
		  printf("of a random starting tree with -d for the same run\n");
		}
	      else
		{
		  printf("\n Error you specified a starting tree with -t AND the computation\n");
		  printf("of a random starting tree with -d for the same run\n");
		}
	    }
	}
      errorExit(-1);
    }

  if(treeSet && constraintSet)
    {
      if(processID == 0)
	printf("\n Error you specified a binary constraint tree AND a starting tree for the same run\n");
      errorExit(-1);
    }


  if(treeSet && groupSet)
    {
      if(processID == 0)
	printf("\n Error you specified a multifurcating constraint tree AND a starting tree for the same run\n");
      errorExit(-1);
    }


  if(groupSet && constraintSet)
    {
      if(processID == 0)
	printf("\n Error you specified a bifurcating constraint tree AND a multifurcating constraint tree for the same run\n");
      errorExit(-1);
    }

  if(adef->restart && adef->startingTreeOnly)
    {
      if(processID == 0)
	{
	  printf("\n Error conflicting options: you want to compute only a parsimony starting tree with -y\n");
	  printf(" while you actually specified a starting tree with -t %s\n", tree_file);
	}
      errorExit(-1);
    }
  
  if(adef->mode == TREE_EVALUATION && (!adef->restart))
    {
      if(processID == 0)
	printf("\n Error: please specify a treefile for the tree you want to evaluate with -t\n");
      errorExit(-1);
    }

#ifdef PARALLEL

  if(adef->mode == SPLIT_MULTI_GENE)
    {
      if(processID == 0)
	printf("Multi gene alignment splitting (-f s) not implemented for the MPI-Version\n");
      errorExit(-1);
    }

  if(adef->mode == TREE_EVALUATION)
    {
      if(processID == 0)
	printf("Tree Evaluation mode (-f e) noot implemented for the MPI-Version\n");
      errorExit(-1);
    } 

   if(adef->mode == CALC_BIPARTITIONS)
     {
       if(processID == 0)
	 printf("Computation of bipartitions (-f b) not implemented for the MPI-Version\n");
       errorExit(-1);
     }
   
   if(adef->multipleRuns == 1)
     {
       if(processID == 0)
	 {
	   printf("Error: you are running the parallel MPI program but only want to compute one tree\n");
	   printf("For the MPI version you must specify a number of trees greater than 1 with the -# or -N option\n");
	 }
       errorExit(-1);
     }
#endif

   
   
   if(adef->mode == TREE_EVALUATION && (adef->model == M_GTRCAT || adef->model == M_PROTCAT))
     {
       if(processID == 0)
	 {
	   printf("\n Error: No tree evaluation with GTRCAT/PROTCAT possible\n");
	   printf("the GTRCAT likelihood values are instable at present and should not\n");
	   printf("be used to compare trees based on ML values\n");
	 }
       errorExit(-1);
     }
  
  if(!nameSet)
    {
      if(processID == 0)
	printf("\n Error: please specify a name for this run with -n\n");
      errorExit(-1);
    }
    
  if(! alignmentSet)
    {
      if(processID == 0)
	printf("\n Error: please specify an alignment for this run with -s\n");
      errorExit(-1);
    }


#ifdef WIN32
  if(workdir[0]==0 || workdir[0] != '\\') 
    {
      getcwd(buf,sizeof(buf));
      if( buf[strlen(buf)-1] != '\\') strcat(buf,"\\");
      strcat(buf,workdir);
      if( buf[strlen(buf)-1] != '\\') strcat(buf,"\\");
      strcpy(workdir,buf);
    }
#else
  if(workdir[0]==0 || workdir[0] != '/') 
    {
      getcwd(buf,sizeof(buf));
      if( buf[strlen(buf)-1] != '/') strcat(buf,"/");
      strcat(buf,workdir);
      if( buf[strlen(buf)-1] != '/') strcat(buf,"/");
      strcpy(workdir,buf);
    }
#endif
 
  return;
}




void errorExit(int e)
{
#ifdef PARALLEL
  MPI_Status msgStatus; 
  int i, dummy;

  if(processID == 0)
    {      
      for(i = 1; i < numOfWorkers; i++)	
	MPI_Send(&dummy, 1, MPI_INT, i, FINALIZE, MPI_COMM_WORLD);
  
      MPI_Finalize();
      exit(e);
    }     
  else
    {	
      MPI_Recv(&dummy, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);     
      MPI_Finalize();
      exit(e);
    }
#else 
  exit(e);
#endif
}



static void makeFileNames(void)
{
  int infoFileExists = 0;
#ifdef PARALLEL
  MPI_Status msgStatus; 
#endif

  strcpy(permFileName,         workdir);    
  strcpy(resultFileName,       workdir);
  strcpy(logFileName,          workdir);
  strcpy(checkpointFileName,   workdir);
  strcpy(infoFileName,         workdir);
  strcpy(randomFileName,       workdir);  
  strcpy(bootstrapFileName,    workdir);
  strcpy(bipartitionsFileName, workdir);
  strcpy(ratesFileName,        workdir);
  strcpy(lengthFileName,       workdir);
  strcpy(lengthFileNameModel,  workdir);
  strcpy( perSiteLLsFileName,  workdir);

  strcat(permFileName,         "RAxML_parsimonyTree.");
  strcat(resultFileName,       "RAxML_result.");
  strcat(logFileName,          "RAxML_log.");
  strcat(checkpointFileName,   "RAxML_checkpoint.");
  strcat(infoFileName,         "RAxML_info.");
  strcat(randomFileName,       "RAxML_randomTree."); 
  strcat(bootstrapFileName,    "RAxML_bootstrap.");  
  strcat(bipartitionsFileName, "RAxML_bipartitions.");
  strcat(ratesFileName,        "RAxML_perSiteRates.");
  strcat(lengthFileName,       "RAxML_treeLength.");
  strcat(lengthFileNameModel,  "RAxML_treeLengthModel.");
  strcat( perSiteLLsFileName, "RAxML_perSiteLLs.");

  strcat(permFileName,         run_id);
  strcat(resultFileName,       run_id);
  strcat(logFileName,          run_id);
  strcat(checkpointFileName,   run_id);
  strcat(infoFileName,         run_id);    
  strcat(randomFileName,       run_id);   
  strcat(bootstrapFileName,    run_id); 
  strcat(bipartitionsFileName, run_id);
  strcat(ratesFileName,        run_id);
  strcat(lengthFileName,       run_id);
  strcat(lengthFileNameModel,  run_id);
  strcat(perSiteLLsFileName,   run_id);

  if(processID == 0)
    {
/*      infoFileExists = filexists(infoFileName); */

#ifdef PARALLEL
      {
	int i;

	for(i = 1; i < numOfWorkers; i++)	
	  MPI_Send(&infoFileExists, 1, MPI_INT, i, FINALIZE, MPI_COMM_WORLD);
      }
#endif

      if(infoFileExists)
	{
	  printf("RAxML output files with the run ID <%s> already exist \n", run_id);
	  printf("in directory %s ...... exiting\n", workdir);
#ifdef PARALLEL 	  
	  MPI_Finalize();
	  exit(-1);
#else
	  exit(-1);
#endif
	}     
    }
#ifdef PARALLEL
  else
    {	
      MPI_Recv(&infoFileExists, 1, MPI_INT, 0, FINALIZE, MPI_COMM_WORLD, &msgStatus);
      if(infoFileExists)
	{	 	  
	  MPI_Finalize();
	  exit(-1);
	}    
    }
#endif
}


static void readData(analdef *adef, rawdata *rdta, cruncheddata *cdta, tree *tr)
{
  INFILE = fopen(seq_file, "r");
  
  if (!INFILE)
    {
      if(processID == 0)
	printf( "Could not open sequence file: %s\n", seq_file);
      errorExit(-1);
    }   
  getinput(adef, rdta, cdta, tr); 
  
  fclose(INFILE);   
}



/***********************reading and initializing input ******************/


/********************PRINTING various INFO **************************************/


static void printModelAndProgramInfo(tree *tr, analdef *adef, int argc, char *argv[])
{ 
  if(processID == 0)
    {      
      int i, model;     
      FILE *infoFile = fopen(infoFileName, "a");  
      char modelType[128];

      if(adef->useInvariant)
	strcpy(modelType, "GAMMA+P-Invar");
      else
	strcpy(modelType, "GAMMA");     
      
      printf("\n\nYou are using %s version %s released by Alexandros Stamatakis in %s\n",  programName, programVersion, programDate);
      fprintf(infoFile, "\n\nYou are using %s version %s released by Alexandros Stamatakis in %s\n",  programName, programVersion, programDate);

      if(adef->mode == OPTIMIZE_RATES)
	{
	  printf("\nAlignment has %d columns\n\n",  tr->cdta->endsite);
	  fprintf(infoFile, "\nAlignment has %d columns\n\n",  tr->cdta->endsite);
	}
      else
	{
	  printf("\nAlignment has %d distinct alignment patterns\n\n",  tr->cdta->endsite);
	  fprintf(infoFile, "\nAlignment has %d distinct alignment patterns\n\n",  tr->cdta->endsite);
	}

      if(adef->useInvariant)
	{
	  printf("Found %d invariant alignment patterns that correspond to %d columns \n", tr->numberOfInvariableColumns, tr->weightOfInvariableColumns);
	  fprintf(infoFile, "Found %d invariant alignment patterns that correspond to %d columns \n", tr->numberOfInvariableColumns, tr->weightOfInvariableColumns);
	}

      printf("Proportion of gaps and completely undetermined characters in this alignment: %f\n", adef->gapyness);
      fprintf(infoFile, "Proportion of gaps and completely undetermined characters in this alignment: %f\n", 
	      adef->gapyness);

      switch(adef->mode)
	{       
	case DISTANCE_MODE:
	  printf("Computing pairwise distances\n");
	  fprintf(infoFile, "Computing pairwise distances\n");
	  break;
	case ARNDT_MODE:
	  printf("Arndt-Mode\n");
	  fprintf(infoFile, "Arndt-Mode\n");
	  break;
	case TREE_EVALUATION : 
	  printf("\nRAxML Model Optimization up to an accuracy of %f log likelihood units\n\n", 
		 adef->likelihoodEpsilon);
	  fprintf(infoFile, "\nRAxML Model Optimization up to an accuracy of %f log likelihood units\n\n", 
		  adef->likelihoodEpsilon);
	  break;     
	case  BIG_RAPID_MODE:
	  if(adef->rapidBoot)
	    {
	      if(adef->allInOne)
		{
		  printf("\nRAxML rapid bootstrapping and subsequent ML search\n\n"); 
		  fprintf(infoFile, "\nRAxML rapid bootstrapping and subsequent ML search\n\n");
		}
	      else
		{
		  printf("\nRAxML rapid bootstrapping algorithm\n\n"); 
		  fprintf(infoFile, "\nRAxML rapid bootstrapping algorithm\n\n");
		}
	    }
	  else
	    {
	      printf("\nRAxML rapid hill-climbing mode\n\n"); 
	      fprintf(infoFile, "\nRAxML rapid hill-climbing mode\n\n");
	    }
	  break;	 
	case CALC_BIPARTITIONS:
	  printf("\nRAxML Bipartition Computation: Drawing support values from trees in file %s onto tree in file %s\n\n", 
		 bootStrapFile, tree_file);
	  fprintf(infoFile, "\nRAxML Bipartition Computation: Drawing support values from trees in file %s onto tree in file %s\n\n", 
		  bootStrapFile, tree_file);
	  fclose(infoFile);
	  return;	
	case OPTIMIZE_RATES:        
	  if(!(adef->treeLength || adef->computePerSiteLLs))
	    {
	      printf("\nRAxML optimization of per-site evolutionary rates\n\n");
	      fprintf(infoFile,"\nRAxML optimization of per-site evolutionary rates\n\n");
	    }
	  if(adef->treeLength)
	    {
	      printf("\nRAxML optimization of per-site evolutionary rates and tree-length sliding window\n\n");
	      fprintf(infoFile,"\nRAxML optimization of per-site evolutionary rates and tree-length sliding window\n\n");
	    }
	  if(adef->computePerSiteLLs)
	    {
	      printf("\nRAxML computation of per-site log likelihoods\n\n");
	      fprintf(infoFile,"\nRAxML computation of per-site log likelihoods\n\n");
	    }
	  fclose(infoFile);	    
	  return;	
	case PARSIMONY_ADDITION:
	  printf("\nRAxML stepwise MP addition to incomplete starting tree\n\n");
	  fprintf(infoFile,"\nRAxML stepwise MP addition to incomplete starting tree\n\n");
	  fclose(infoFile);	    
	  return;
	case MEHRING_ALGO:
	  printf("\nRAxML single-sequence position determination algorithm\n\n");
	  fprintf(infoFile,"\nRAxML single-sequence position determination algorithm\n\n");
	  break;
	default:
	  printf("Oups, forgot to implement mode description %d exiting\n", adef->mode);
	  exit(-1);
	}
      
      if(tr->NumberOfModels > 1)
	{       	  
	  if(adef->perGeneBranchLengths)
	    {
	      printf("Partitioned Data Mode: Using %d distinct models/partitions with individual per partition branch length optimization\n", tr->NumberOfModels);
	      fprintf(infoFile, "Partitioned Data Mode: Using %d distinct models/partitions with individual per partition branch length optimization\n", tr->NumberOfModels);
	      printf(           "\n\n");
	      fprintf(infoFile, "\n\n");
	    }
	  else
	    {
	      printf("Partitioned Data Mode: Using %d distinct models/partitions with joint branch length optimization\n", 
		     tr->NumberOfModels);
	      fprintf(infoFile, "Partitioned Data Mode: Using %d distinct models/partitions with joint branch length optimization\n", 
		      tr->NumberOfModels);	  
	      printf(           "\n\n");
	      fprintf(infoFile, "\n\n");
	    }
	}
      
      if(adef->rapidBoot)
	{
	  if(adef->allInOne)
	    {
	      printf("\nExecuting %d rapid bootstrap inferences and thereafter a thorough ML search \n\n", adef->multipleRuns);
	      fprintf(infoFile, "\nExecuting %d rapid bootstrap inferences and thereafter a thorough ML search \n\n", adef->multipleRuns);
	    }
	  else
	    {
	      printf("\nExecuting %d rapid bootstrap inferences\n\n", adef->multipleRuns); 
	      fprintf(infoFile, "\nExecuting %d rapid bootstrap inferences\n\n", adef->multipleRuns);
	    }
	}
      else
	{
	  if(adef->boot)
	    {	 
	      if(adef->multipleRuns > 1)
		{
		  printf("Executing %d non-parametric bootstrap inferences\n\n", adef->multipleRuns);
		  fprintf(infoFile, "Executing %d non-parametric bootstrap inferences\n\n", adef->multipleRuns);
		}
	      else
		{
		  printf("Executing %d non-parametric bootstrap inference\n\n", adef->multipleRuns);
		  fprintf(infoFile, "Executing %d non-parametric bootstrap inference\n\n", adef->multipleRuns);
		}
	    }
	  else
	    {
	      char treeType[1024];
	      
	      if(adef->restart)
		strcpy(treeType, "user-specifed");
	      else
		{
		  if(adef->randomStartingTree)		
		    strcpy(treeType, "distinct complete random");
		  else
		    strcpy(treeType, "distinct randomized MP");
		}


	      if(adef->multipleRuns > 1)
		{
		  printf("Executing %d inferences on the original alignment using %d %s trees\n\n", 
			 adef->multipleRuns, adef->multipleRuns, treeType);
		  fprintf(infoFile, "Executing %d inferences on the original alignment using %d %s trees\n\n", 
			  adef->multipleRuns, adef->multipleRuns, treeType);
		}
	      else
		{
		  printf("Executing %d inference on the original alignment using a %s tree\n\n", 
			 adef->multipleRuns, treeType);
		  fprintf(infoFile, "Executing %d inference on the original alignment using a %s tree\n\n", 
			  adef->multipleRuns, treeType);
		}
	    }
	}         

      if(tr->rateHetModel == GAMMA || tr->rateHetModel == GAMMA_I)
	{
	   printf( "All free model parameters will be estimated by RAxML\n");  
	   printf( "%s model of rate heteorgeneity, ML estimate of alpha-parameter\n", modelType);
	   printf( "%s Model parameters will be estimated up to an accuracy of %2.10f Log Likelihood units\n\n", 
		   modelType, adef->likelihoodEpsilon);	
	   fprintf(infoFile, "All free model parameters will be estimated by RAxML\n");  
	   fprintf(infoFile, "%s model of rate heterogeneity, ML estimate of alpha-parameter\n", modelType);
	   fprintf(infoFile, "%s Model parameters will be estimated up to an accuracy of %2.10f Log Likelihood units\n\n", 
		   modelType, adef->likelihoodEpsilon);
	}
      else
	{
	  if(adef->useMixedModel)
	    { 		 
	      printf( "All free model parameters will be estimated by RAxML\n");
	      printf( "ML estimate of %d per site rate categories\n", adef->categories);      
	      printf( "Likelihood of final tree will be evaluated and optimized under %s\n", modelType);
	      printf( "Final model parameters will be estimated up to an accuracy of %2.10f Log Likelihood units\n\n", 
		      adef->likelihoodEpsilon);

	      fprintf(infoFile, "All free model parameters will be estimated by RAxML\n");
	      fprintf(infoFile, "ML estimate of %d per site rate categories\n", adef->categories);      
	      fprintf(infoFile, "Likelihood of final tree will be evaluated and optimized under %s\n", modelType);
	      fprintf(infoFile, "Final model parameters will be estimated up to an accuracy of %2.10f Log Likelihood units\n\n", 
		      adef->likelihoodEpsilon);		 
	    }
	  else
	    {
	      printf( "Approximation of rate heterogeneity only!\n"); 
	      printf( "All free model parameters will be estimated by RAxML\n");
	      printf( "ML estimate of %d per site rate categories\n", adef->categories);   
	      printf( "WARNING: CAT-based likelihood values should NEVER be used to COMPARE trees!\n\n");
	       
	      fprintf(infoFile, "Approximation of rate heterogeneity only!\n"); 
	      fprintf(infoFile, "All free model parameters will be estimated by RAxML\n");
	      fprintf(infoFile, "ML estimate of %d per site rate categories\n", adef->categories);   
	      fprintf(infoFile, "WARNING: CAT-based likelihood values should NEVER be used to COMPARE trees!\n\n");	      
	    }	
	}


      for(model = 0; model < tr->NumberOfModels; model++)
	{
	  printf("Partition: %d\n", model);
	  printf("Name: %s\n", tr->partitionData[model].partitionName);

	  fprintf(infoFile, "Partition: %d\n", model);
	  fprintf(infoFile, "Name: %s\n", tr->partitionData[model].partitionName);

	  switch(tr->partitionData[model].dataType)
	    {
	    case DNA_DATA:
	      printf("DataType: DNA\n");	      
	      printf("Substitution Matrix: GTR\n");
	     
	      if(adef->boot == 0)
		{
		  printf("Empirical Base Frequencies:\n");
		  printf("pi(A): %f pi(C): %f pi(G): %f pi(T): %f",  
			 tr->frequencies_DNA[model * 4 + 0], tr->frequencies_DNA[model * 4 + 1], 
			 tr->frequencies_DNA[model * 4 + 2], tr->frequencies_DNA[model * 4 + 3]);
		}
	      else
		{
		  printf("Empirical Base Frequencies will not be printed for Bootstrapping\n");
		}
	      
	      fprintf(infoFile, "DataType: DNA\n");	      
	      fprintf(infoFile, "Substitution Matrix: GTR\n");
	      
	      if(adef->boot == 0)
		{
		  fprintf(infoFile, "Empirical Base Frequencies:\n");
		  fprintf(infoFile, "pi(A): %f pi(C): %f pi(G): %f pi(T): %f",  
			  tr->frequencies_DNA[model * 4 + 0], tr->frequencies_DNA[model * 4 + 1], 
			  tr->frequencies_DNA[model * 4 + 2], tr->frequencies_DNA[model * 4 + 3]);
		}	
	      else
		{
		  fprintf(infoFile, "Empirical Base Frequencies will not be printed for Bootstrapping\n");
		}
	      break;
	    case AA_DATA:
	      {
		char *protStrings[10] = {"DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM"};
		char basesPROT[20] = {'A', 'R', 'N' , 'D', 'C', 'Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V'};		

		assert(tr->partitionData[model].protModels >= 0 && tr->partitionData[model].protModels < 10);		

		printf("DataType: AA\n");
		printf("Substitution Matrix: %s\n", protStrings[tr->partitionData[model].protModels]);
		printf("%s Base Frequencies:\n", (tr->partitionData[model].protFreqs == 1)?"Empirical":"Fixed");
	      
		if(adef->boot == 0)
		  {
		    int k;

		    for(k = 0; k < 20; k++)
		      {
			if(k % 4 == 0 && k > 0)			  
			  printf("\n");		      
			 			
			printf("pi(%c): %f ",  basesPROT[k], tr->frequencies_AA[model * 20 + k]);		  
		      }		
		  }
		else
		  {
		    printf("Base Frequencies will not be printed for Bootstrapping\n");
		  }
		
		fprintf(infoFile, "DataType: AA\n");
		fprintf(infoFile, "Substitution Matrix: %s\n", protStrings[tr->partitionData[model].protModels]);
		fprintf(infoFile, "%s Base Frequencies:\n", (tr->partitionData[model].protFreqs == 1)?"Empirical":"Fixed");
	      
		if(adef->boot == 0)
		  {
		    int k;

		    for(k = 0; k < 20; k++)
		      {
			if(k % 4 == 0)			  
			  fprintf(infoFile, "\n");		      
			 			
			fprintf(infoFile, "pi(%c): %f ",  basesPROT[k], tr->frequencies_AA[model * 20 + k]);		  
		      }		
		  }
		else
		  {
		    fprintf(infoFile, "Base Frequencies will not be printed for Bootstrapping\n");
		  }
		

	      }
	      break;
	    default:
	      assert(0);
	    }
	  
	  printf("\n\n\n");
	  fprintf(infoFile,"\n\n\n"); 
	}              
      
      printf("\n");
      fprintf(infoFile, "\n");
      
      fprintf(infoFile,"RAxML was called as follows:\n\n");
      for(i = 0; i < argc; i++)
	fprintf(infoFile,"%s ", argv[i]);
      fprintf(infoFile,"\n\n\n");  
      
      fclose(infoFile);
    }
}

void printResult(tree *tr, analdef *adef, boolean finalPrint)
{
  FILE *logFile;
  char temporaryFileName[1024] = "", treeID[64] = "";

  strcpy(temporaryFileName, resultFileName);

  switch(adef->mode)
    {         
    case TREE_EVALUATION:     
      

      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH);
            
      logFile = fopen(temporaryFileName, "w");
      fprintf(logFile, "%s", tr->tree_string);
      fclose(logFile);          

      if(adef->perGeneBranchLengths)
	printTreePerGene(tr, adef, temporaryFileName, "w");


      break;     
    case BIG_RAPID_MODE:
      if(!adef->boot)
	{
	  if(adef->multipleRuns > 1)
	    {	  	 	  
	      sprintf(treeID, "%d", tr->treeID);	  	  
	      strcat(temporaryFileName, ".RUN.");
	      strcat(temporaryFileName, treeID);	  	 	      	
	    }
	  
	  if((adef->model == M_GTRCAT || adef->model == M_PROTCAT) && (adef->useMixedModel == 0))	    
	    {
	      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES);	    	    
	      logFile = fopen(temporaryFileName, "w");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile);
	    }
	  else	    
	    {
	      if(finalPrint)
		{		  
		  Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, 
			      SUMMARIZE_LH);

		  logFile = fopen(temporaryFileName, "w");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);

		  if(adef->perGeneBranchLengths)
		    printTreePerGene(tr, adef, temporaryFileName, "w");
		}
	      else
		{
		  Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, 
			      NO_BRANCHES);
		  logFile = fopen(temporaryFileName, "w");
		  fprintf(logFile, "%s", tr->tree_string);
		  fclose(logFile);
		}
	    }	 
	}
      break;       
    default:
      printf("FATAL ERROR call to printResult from undefined STATE %d\n", adef->mode);
      exit(-1);
      break;
    }
}

void printBootstrapResult(tree *tr, analdef *adef, boolean finalPrint)
{
  if(processID == 0)
    {
      FILE *logFile;
      
      if(adef->mode == BIG_RAPID_MODE && (adef->boot || adef->rapidBoot))
	{
#ifndef PARALLEL
	  if(adef->bootstrapBranchLengths)	    
	    {	      
	      Tree2String(tr->tree_string, tr, tr->start->back, TRUE, TRUE, FALSE, FALSE, finalPrint, adef, SUMMARIZE_LH);	    
	      logFile = fopen(bootstrapFileName, "a");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile); 
	      if(adef->perGeneBranchLengths)
		printTreePerGene(tr, adef, bootstrapFileName, "a");
	    }
	  else
	    {	      	     
	      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES);
	      logFile = fopen(bootstrapFileName, "a");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile); 
	    }
#else	  
	  logFile = fopen(bootstrapFileName, "a");
	  fprintf(logFile, "%s", tr->tree_string);
	  fclose(logFile);     
#endif	  
	}
      else
	{
	  printf("FATAL ERROR in  printBootstrapResult\n");
	  exit(-1);	 
	}
    }
}



void printBipartitionResult(tree *tr, analdef *adef, boolean finalPrint)
{
  if(processID == 0 || adef->allInOne)
    {
      FILE *logFile;
      
      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, TRUE, finalPrint, adef, NO_BRANCHES);
      logFile = fopen(bipartitionsFileName, "a");
      fprintf(logFile, "%s", tr->tree_string);
      fclose(logFile);   
    }
}



void printLog(tree *tr, analdef *adef, boolean finalPrint)
{
  FILE *logFile;
  char temporaryFileName[1024] = "", checkPoints[1024] = "", treeID[64] = "";
  double lh, t;
  
  lh = tr->likelihood;
  t = gettime() - masterTime;

  strcpy(temporaryFileName, logFileName);
  strcpy(checkPoints,       checkpointFileName);

  switch(adef->mode)
    {    
    case TREE_EVALUATION:	 
      logFile = fopen(temporaryFileName, "a");

      printf("%f %f\n", t, lh);
      fprintf(logFile, "%f %f\n", t, lh);

      fclose(logFile);     
      break;      
    case BIG_RAPID_MODE:
      if(adef->boot || adef->rapidBoot)
	{
	  /* testing only printf("%f %f\n", t, lh);*/
	  /* NOTHING PRINTED so far */
	}
      else
	{
	  if(adef->multipleRuns > 1)
	    {	  	 	  
	      sprintf(treeID, "%d", tr->treeID);	  	  
	      strcat(temporaryFileName, ".RUN.");
	      strcat(temporaryFileName, treeID);	  	 
	      
	      strcat(checkPoints, ".RUN.");
	      strcat(checkPoints, treeID);	      	      
	    }


	  if(!adef->checkpoints)
	    {
	      logFile = fopen(temporaryFileName, "a");
#ifndef PARALLEL	      
	      printf("%f %f\n", t, lh);
#endif
	      fprintf(logFile, "%f %f\n", t, lh);
	      
	      fclose(logFile);
	    }
	  else
	    {
	      logFile = fopen(temporaryFileName, "a");
#ifndef PARALLEL	      
	      printf("%f %f %d\n", t, lh, tr->checkPointCounter);
#endif
	      fprintf(logFile, "%f %f %d\n", t, lh, tr->checkPointCounter);
	      
	      fclose(logFile);
	      
	      strcat(checkPoints, ".");

	      sprintf(treeID, "%d", tr->checkPointCounter);
	      strcat(checkPoints, treeID);
	     
	      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES);

	      logFile = fopen(checkPoints, "a");
	      fprintf(logFile, "%s", tr->tree_string);
	      fclose(logFile);

	      tr->checkPointCounter++;
	    }
	}
      break;       
    default:
      printf("FATAL ERROR call to printLog from undefined STATE %d\n", adef->mode);
      exit(-1);
      break;
    }
}



void printStartingTree(tree *tr, analdef *adef, boolean finalPrint)
{  
  if(adef->boot)
    {          
      /* not printing starting trees for bootstrap */
    }
  else
    {
      FILE *treeFile;
      char temporaryFileName[1024] = "", treeID[64] = "";
   
      Tree2String(tr->tree_string, tr, tr->start->back, FALSE, TRUE, FALSE, FALSE, finalPrint, adef, NO_BRANCHES);
          
      if(adef->randomStartingTree)	    
	strcpy(temporaryFileName, randomFileName);	    
      else
	strcpy(temporaryFileName, permFileName);

      if(adef->multipleRuns > 1)
	{	  	 	  
	  sprintf(treeID, "%d", tr->treeID);	  	  
	  strcat(temporaryFileName, ".RUN.");
	  strcat(temporaryFileName, treeID);	  	 
	}
     	  	 
      treeFile = fopen(temporaryFileName, "a");	
      fprintf(treeFile, "%s", tr->tree_string);
      fclose(treeFile);	    	
    }
}

void writeInfoFile(analdef *adef, tree *tr, double t)
{  
  if(processID == 0)
    {
      FILE *infoFile = fopen(infoFileName, "a");           

      switch(adef->mode)
	{
	case TREE_EVALUATION:	 	  
	  break;      
	case BIG_RAPID_MODE:
	  if(adef->boot || adef->rapidBoot)
	    {
	      if(!adef->initialSet)	       
		{
		  fprintf(infoFile, "Bootstrap[%d]: Time %f bootstrap likelihood %f, best rearrangement setting %d\n", tr->treeID, t, tr->likelihood,  adef->bestTrav);		
		  printf("Bootstrap[%d]: Time %f bootstrap likelihood %f, best rearrangement setting %d\n", tr->treeID, t, tr->likelihood,  adef->bestTrav);
		}
	      else		
		{
		  fprintf(infoFile, "Bootstrap[%d]: Time %f bootstrap likelihood %f\n", tr->treeID, t, tr->likelihood);	
		  printf("Bootstrap[%d]: Time %f bootstrap likelihood %f\n", tr->treeID, t, tr->likelihood);	
		}
	    }
	  else
	    {
	      if((adef->model == M_GTRCAT || adef->model == M_PROTCAT) && !adef->useMixedModel)
		{		  		 
		  if(!adef->initialSet)		   
		    fprintf(infoFile, "Inference[%d]: Time %f CAT-likelihood %f, best rearrangement setting %d\n", tr->treeID, t, tr->likelihood,  adef->bestTrav);
		  else		  
		    fprintf(infoFile, "Inference[%d]: Time %f CAT-likelihood %f\n", tr->treeID, t, tr->likelihood);		    				 
		}
	      else
		{	
		  int model;
		  char modelType[128];

		  if(adef->useInvariant)
		    strcpy(modelType, "GAMMA+P-Invar");
		  else
		    strcpy(modelType, "GAMMA");		  	       

		  if(!adef->initialSet)		    		     	  
		    fprintf(infoFile, "Inference[%d]: Time %f %s-likelihood %f, best rearrangement setting %d, ", 
			    tr->treeID, t, modelType, tr->likelihood,  adef->bestTrav);
		  else		  
		    fprintf(infoFile, "Inference[%d]: Time %f %s-likelihood %f, ", 
			    tr->treeID, t, modelType, tr->likelihood);		    		  
		   
		  for(model = 0; model < tr->NumberOfModels; model++)		    		    
		    {
		      fprintf(infoFile, "alpha[%d]: %f ", model, tr->alphas[model]);
		      if(adef->useInvariant)
			fprintf(infoFile, "invar[%d]: %f ", model, tr->invariants[model]);
#ifndef PARALLEL
		      if(adef->model == M_GTRCAT || adef->model == M_GTRGAMMA)
			{
			  int k;

			  fprintf(infoFile, "rates[%d] ac ag at cg ct gt: ",model);
			  for(k = 0; k < DNA_RATES; k++)			    
			    fprintf(infoFile, "%f ", tr->initialRates_DNA[model * DNA_RATES + k]);			    
			}
		      fprintf(infoFile, "1.0 ");
#endif
		    }
		    		  
		  fprintf(infoFile, "\n");
		}
	    }
	  break;
	default:
	  assert(0);
	}

      fclose(infoFile);
    }
}

static void finalizeInfoFile(tree *tr, analdef *adef)
{
  if(processID == 0)
    {
      FILE *infoFile = fopen(infoFileName, "a");
      double t;
      int model;

      t = gettime() - masterTime;

      switch(adef->mode)
	{       
	case TREE_EVALUATION :
	  printf("\n\nOverall Time for Tree Evaluation %f\n", t);	 
	  printf("Final GAMMA  likelihood: %f\n", tr->likelihood);

	  fprintf(infoFile, "\n\nOverall Time for Tree Evaluation %f\n", t);       
	  fprintf(infoFile, "Final GAMMA  likelihood: %f\n", tr->likelihood);

	  {
	    int 
	      params,
	      paramsBrLen;

	    if(tr->NumberOfModels == 1)
	      {
		if(adef->useInvariant)
		  {
		    params      = 1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */;
		    paramsBrLen = 1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */ + 
		      (2 * tr->mxtips - 3);
		  }
		else
		  {
		    params      = 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */;
		    paramsBrLen = 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */ + 
		      (2 * tr->mxtips - 3);
		  }
	      }
	    else
	      {
		if(tr->multiBranch)
		  {
		    if(adef->useInvariant)
		      {
			params      = tr->NumberOfModels * (1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */);
			paramsBrLen = tr->NumberOfModels * (1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */ + 
							    (2 * tr->mxtips - 3));
		      }
		    else
		      {
			params      = tr->NumberOfModels * (5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */);
			paramsBrLen = tr->NumberOfModels * (5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */ + 
							    (2 * tr->mxtips - 3));
		      }
		  }
		else
		  {
		    if(adef->useInvariant)
		      {
			params      = tr->NumberOfModels * (1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */);
			paramsBrLen = tr->NumberOfModels * (1 /* INVAR */ + 5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */) 
			  + (2 * tr->mxtips - 3);
		      }
		    else
		      {
			params      = tr->NumberOfModels * (5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */);
			paramsBrLen = tr->NumberOfModels * (5 /* RATES */ + 3 /* freqs */ + 1 /* alpha */) 
			  + (2 * tr->mxtips - 3);
		      }

		  }
	      }
		
	    if(!tr->mixedData && tr->partitionData[0].dataType == DNA_DATA)
	      {
		printf("Number of free parameters for AIC-TEST(BR-LEN): %d\n",    paramsBrLen);
		printf("Number of free parameters for AIC-TEST(NO-BR-LEN): %d\n", params);
		fprintf(infoFile, "Number of free parameters for AIC-TEST(BR-LEN): %d\n",    paramsBrLen);
		fprintf(infoFile, "Number of free parameters for AIC-TEST(NO-BR-LEN): %d\n", params);
	      }
	    
	  }
	
	  printf("\n\n");
	  fprintf(infoFile, "\n\n");

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
	  	
	  printf("Final tree written to:                 %s\n", resultFileName);  
	  printf("Execution Log File written to:         %s\n", logFileName);
	   
	  fprintf(infoFile, "Final tree written to:                 %s\n", resultFileName);  
	  fprintf(infoFile, "Execution Log File written to:         %s\n", logFileName);	  
	
	  break;  
	case  BIG_RAPID_MODE:
	  if(adef->boot)
	    {
	      printf("\n\nOverall Time for %d Bootstraps %f\n", adef->multipleRuns, t);
	      printf("\n\nAverage Time per Bootstrap %f\n", (double)(t/((double)adef->multipleRuns)));
	      printf("All %d bootstrapped trees written to: %s\n", adef->multipleRuns, bootstrapFileName);
	      
	      fprintf(infoFile, "\n\nOverall Time for %d Bootstraps %f\n", adef->multipleRuns, t);
	      fprintf(infoFile, "Average Time per Bootstrap %f\n", (double)(t/((double)adef->multipleRuns)));	     
	      fprintf(infoFile, "\n\nAll %d bootstrapped trees written to: %s\n", adef->multipleRuns, bootstrapFileName);	     
	    }
	  else
	    {
	      if(adef->multipleRuns > 1)
		{
		  double avgLH = 0;
		  double bestLH = unlikely;
		  int i, bestI  = 0;
		  
		  for(i = 0; i < adef->multipleRuns; i++)
		    {     
		      avgLH   += tr->likelihoods[i];
		      if(tr->likelihoods[i] > bestLH)
			{
			  bestLH = tr->likelihoods[i];
			  bestI  = i;
			}
		    }
		  avgLH /= ((double)adef->multipleRuns);

		  printf("\n\nOverall Time for %d Inferences %f\n", adef->multipleRuns, t);
		  printf("Average Time per Inference %f\n", (double)(t/((double)adef->multipleRuns)));		 
		  printf("Average Likelihood   : %f\n", avgLH);
		  printf("\n");
		  printf("Best Likelihood in run number %d: likelihood %f\n\n", bestI, bestLH);

		  if(adef->checkpoints)   
		    printf("Checkpoints written to:                 %s.RUN.%d.* to %d.*\n", checkpointFileName, 0, adef->multipleRuns - 1);  
		  if(!adef->restart)
		    {
		      if(adef->randomStartingTree)
			printf("Random starting trees written to:       %s.RUN.%d to %d\n", randomFileName, 0, adef->multipleRuns - 1);
		      else
			printf("Parsimony starting trees written to:    %s.RUN.%d to %d\n", permFileName, 0, adef->multipleRuns - 1);   	            
		    }					  
		  printf("Final trees written to:                 %s.RUN.%d to %d\n", resultFileName,  0, adef->multipleRuns - 1);  
		  printf("Execution Log Files written to:         %s.RUN.%d to %d\n", logFileName, 0, adef->multipleRuns - 1);   
		  printf("Execution information file written to:  %s\n", infoFileName);
		  

		  fprintf(infoFile, "\n\nOverall Time for %d Inferences %f\n", adef->multipleRuns, t);
		  fprintf(infoFile, "Average Time per Inference %f\n", (double)(t/((double)adef->multipleRuns)));		 
		  fprintf(infoFile, "Average Likelihood   : %f\n", avgLH);
		  fprintf(infoFile, "\n");
		  fprintf(infoFile, "Best Likelihood in run number %d: likelihood %f\n\n", bestI, bestLH); 
		  if(adef->checkpoints)   
		    fprintf(infoFile, "Checkpoints written to:                %s.RUN.%d.* to %d.*\n", checkpointFileName, 0, adef->multipleRuns - 1);  
		  if(!adef->restart)
		    {
		      if(adef->randomStartingTree)
			fprintf(infoFile, "Random starting trees written to:      %s.RUN.%d to %d\n", randomFileName, 0, adef->multipleRuns - 1);
		      else
			fprintf(infoFile, "Parsimony starting trees written to:   %s.RUN.%d to %d\n", permFileName, 0, adef->multipleRuns - 1);   	            
		    }					  
		  fprintf(infoFile, "Final trees written to:                %s.RUN.%d to %d\n", resultFileName,  0, adef->multipleRuns - 1);  
		  fprintf(infoFile, "Execution Log Files written to:        %s.RUN.%d to %d\n", logFileName, 0, adef->multipleRuns - 1);   
		  fprintf(infoFile, "Execution information file written to: %s\n", infoFileName);
		  		  		 		   
		}
	      else
		{
		  printf("\n\nOverall Time for 1 Inference %f\n", t);		  
		  printf("Likelihood   : %f\n", tr->likelihood);
		  printf("\n\n");	     

		  if(adef->checkpoints)   
		  printf("Checkpoints written to:                %s.*\n", checkpointFileName);  
		  if(!adef->restart)
		    {
		      if(adef->randomStartingTree)
			printf("Random starting tree written to:       %s\n", randomFileName);
		      else
			printf("Parsimony starting tree written to:    %s\n", permFileName);   	            
		    }					  
		  printf("Final tree written to:                 %s\n", resultFileName);  
		  printf("Execution Log File written to:         %s\n", logFileName);   
		  printf("Execution information file written to: %s\n",infoFileName);

		  
		  
		  fprintf(infoFile, "\n\nOverall Time for 1 Inference %f\n", t);		  
		  fprintf(infoFile, "Likelihood   : %f\n", tr->likelihood);
		  fprintf(infoFile, "\n\n");

		  if(adef->checkpoints)   
		    fprintf(infoFile, "Checkpoints written to:                %s.*\n", checkpointFileName);  
		  if(!adef->restart)
		    {
		      if(adef->randomStartingTree)
			fprintf(infoFile, "Random starting tree written to:       %s\n", randomFileName);
		      else
			fprintf(infoFile, "Parsimony starting tree written to:    %s\n", permFileName);   	            
		    }					  
		  fprintf(infoFile, "Final tree written to:                 %s\n", resultFileName);  
		  fprintf(infoFile, "Execution Log File written to:         %s\n", logFileName);   
		  fprintf(infoFile, "Execution information file written to: %s\n",infoFileName);
		  		 		  
		}
	    }
	    
	  break;	 
	case CALC_BIPARTITIONS:
	  printf("\n\nTime for Computation of Bipartitions %f\n", t);
	  printf("Tree with bipartitions written to file:  %s\n", bipartitionsFileName);
	  printf("Execution information file written to :  %s\n",infoFileName);
	

	  fprintf(infoFile, "\n\nTime for Computation of Bipartitions %f\n", t);
	  fprintf(infoFile, "Tree with bipartitions written to file:  %s\n", bipartitionsFileName);

	 
	  break;
	case OPTIMIZE_RATES:
	  if(! (adef->computePerSiteLLs || adef->treeLength))
	    {
	      printf("\n\nTime for Optimization of per-site rates %f\n", t);
	      printf("Optimized rates written to file:  %s\n", ratesFileName);
	      printf("Execution information file written to :  %s\n",infoFileName);
	     
	      
	      fprintf(infoFile, "\n\nTime for Optimization of per-site rates %f\n", t);
	      fprintf(infoFile, "Optimized rates written to file:  %s\n", ratesFileName);
	    }
	  
	   if(adef->treeLength)
	    {
	      printf("\n\nTime for Optimization of per-site rates and sliding window tree length %f\n", t);
	      printf("Optimized rates written to file:  %s, \n sliding window data written to %s and %s\n", ratesFileName, lengthFileName, lengthFileNameModel);
	      printf("Execution information file written to :  %s\n",infoFileName);
	     
	      
	      fprintf(infoFile, "\n\nTime for Optimization of per-site rates and sliding window tree length %f\n", t);
	      fprintf(infoFile, "Optimized rates written to file:  %s, \n sliding window data written to %s and %s\n", ratesFileName, lengthFileName, lengthFileNameModel);
	    }
    
	   if(adef->computePerSiteLLs)
	     {
	       printf("\n\nTime for Optimization of per-site log likelihoods %f\n", t);
	       printf("Per-site Log Likelihoods written to File %s in Tree-Puzzle format\n",  perSiteLLsFileName);
	       printf("Execution information file written to :  %s\n",infoFileName);
	     
	       
	       fprintf(infoFile, "\n\nTime for Optimization of per-site log likelihoods %f\n", t);
	       fprintf(infoFile, "Per-site Log Likelihoods written to File %s in Tree-Puzzle format\n",  perSiteLLsFileName);
	     }
	 
	  break;
	case PARSIMONY_ADDITION:
	  printf("\n\nTime for MP stepwise addition %f\n", t);	 
	  printf("Execution information file written to :  %s\n",infoFileName);
	  printf("Complete parsimony tree written to:      %s\n", permFileName); 

	 

	  fprintf(infoFile, "\n\nTime for MP stepwise addition %f\n", t);	 
	  fprintf(infoFile, "Complete parsimony tree written to:      %s\n", permFileName); 

	 
	  break;
	default:
	  assert(0);
	}
      fclose(infoFile);
    }

}



/********************PRINTING various INFO **************************************/

/************************************************************************************/

static void computeLHTest(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int numberOfTrees = 0, i;
  char ch; 
  double bestLH, currentLH, weightSum = 0.0;
  double *bestVector, *otherVector;


  bestVector = (double*)malloc(sizeof(double) * tr->cdta->endsite);
  otherVector = (double*)malloc(sizeof(double) * tr->cdta->endsite);

  for(i = 0; i < tr->cdta->endsite; i++)
    weightSum += (double)(tr->cdta->aliaswgt[i]);

  modOpt(tr, adef);
  printf("Model optimization, best Tree: %f\n", tr->likelihood);
  bestLH = tr->likelihood;


  evaluateGenericInitrav(tr, tr->start);   

  evaluateGenericVector(tr, tr->start, bestVector);
      
  INFILE = fopen(bootStrapFileName, "r");       
  while((ch = getc(INFILE)) != EOF)
    {
      if(ch == ';')
	numberOfTrees++;
    }	 
  rewind(INFILE);
 
  printf("Found %d trees in File %s\n", numberOfTrees, bootStrapFileName);
 
  for(i = 0; i < numberOfTrees; i++)
    {              
      treeReadLen(INFILE, tr, adef);      
      treeEvaluate(tr, 2);
      tr->start = tr->nodep[1];

      evaluateGenericInitrav(tr, tr->start);         

      currentLH = tr->likelihood;
      if(currentLH > bestLH)
	{
	  printf("Better tree found %d at %f\n", i, currentLH);
	  /*exit(1);*/
	}
      /*printf("Tree %d %f\n",i, tr->likelihood);*/
     
      evaluateGenericVector(tr, tr->start, otherVector);
     
      {
	 int j;
	 double temp, wtemp, sum, sum2, sd;

	 sum = 0.0;
	 sum2 = 0.0;

	 for (j = 0; j < tr->cdta->endsite; j++) 
	   {
	     temp  = bestVector[j] - otherVector[j];
	     wtemp = tr->cdta->aliaswgt[j] * temp;
	     sum  += wtemp;
	     sum2 += wtemp * temp;
	   }
	 
	 sd = sqrt( weightSum * (sum2 - sum*sum / weightSum)
		    / (weightSum - 1) );	       
       
	 printf("Tree: %d Likelihood: %f D(LH): %f SD: %f Significantly Worse: %s\n", i, currentLH, currentLH - bestLH, sd, (sum > 1.95996 * sd) ? "Yes" : " No");	 
      }
    }
      
  fclose(INFILE); 
  
  free(bestVector);
  free(otherVector);
  exit(0);
}

static void computePerSiteLLs(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int numberOfTrees = 0, i, j;
  char ch;
  double *otherVector;
  FILE *tlf;	   
	  
  tlf = fopen( perSiteLLsFileName, "w");
 
  otherVector = (double*)malloc(sizeof(double) * tr->cdta->endsite);  

  allocNodex(tr, adef); 
      
  INFILE = fopen(bootStrapFileName, "r");       
  while((ch = getc(INFILE)) != EOF)
    {
      if(ch == ';')
	numberOfTrees++;
    }	 
  rewind(INFILE);
 
  printf("Found %d trees in File %s\n", numberOfTrees, bootStrapFileName);
 
  fprintf(tlf, "  %d  %d\n", numberOfTrees, tr->cdta->endsite);

  for(i = 0; i < numberOfTrees; i++)
    {              
      treeReadLen(INFILE, tr, adef);      
      if(i == 0)	
	modOpt(tr, adef);	  	
      else
	treeEvaluate(tr, 2);

      printf("Tree %d: %f\n", i, tr->likelihood);

      tr->start = tr->nodep[1];

      evaluateGenericInitrav(tr, tr->start);      
     
      evaluateGenericVector(tr, tr->start, otherVector);           
   
      fprintf(tlf, "tr%d\t", i + 1);
      for(j = 0; j < tr->cdta->endsite; j++)	
	{
	  fprintf(tlf, "%f ", otherVector[j]);	
	}
      fprintf(tlf, "\n");     
    }
      
  fclose(INFILE); 
  fclose(tlf);  
  free(otherVector); 
}

#ifdef _USE_PTHREADS

#ifndef _MAC
#include <sched.h>

static void pinThread2Cpu(int myTid)
{
  char *coreSteppingStr;
  int myCore, len;
  cpu_set_t cpuMask;    

  coreSteppingStr = getenv("SCHEDULE");
  len = coreSteppingStr ? strlen(coreSteppingStr) : 0;
  
  myCore = myTid;

  if (myTid < len)
    {
      if ((coreSteppingStr[myTid] >= '0') && (coreSteppingStr[myTid] <= '9'))
	myCore = coreSteppingStr[myTid] - '0';

      if ((coreSteppingStr[myTid] >= 'a') && (coreSteppingStr[myTid] <= 'f'))
	myCore  = coreSteppingStr[myTid] - 'a' + 10;

      if ((coreSteppingStr[myTid] >= 'A') && (coreSteppingStr[myTid] <= 'F'))
	myCore = coreSteppingStr[myTid] - 'A' + 10;
    }
  
  CPU_ZERO(&cpuMask);
  CPU_SET(myCore, &cpuMask);

  if (sched_setaffinity(0, sizeof(cpuMask), &cpuMask))
    {
      printf("Error while scheduling Thread #%d to CPU %d\n",
	     myTid, myCore);
      exit(1);
    }

  /* printf("Scheduled Thread #%d to logical CPU %d\n", myTid, myCore); */
  return;
} 

#endif

typedef struct {
  tree *tr;
  int threadNumber;
} threadData;


static void calcBounds(int tid, const int n, int start, int end, int *l, int *u)
{      
  int span = end - start;

  /* LTD */
  /* assert(span % n == 0); */

  if(span % n == 0)    
    span = span / n;	
  else
    span = 1 + (span / n);
	
  *l = start + tid * span;
  if(tid == n - 1)
    *u = end;
  else
    *u = *l + span;   
}

#ifdef _LOCAL_DATA

static void strided_Bounds(int tid, int endsite, int n, int *startIndex, int *endIndex)
{
  int endsiteL = endsite - tid;

  if(endsiteL % n == 0)
    endsiteL = endsiteL / n;
  else
    endsiteL = 1 + (endsiteL / n);
  
  *startIndex = 0;
  *endIndex   = endsiteL;
}

static void collectDouble(double *dest, double *source, const int totallength, const int stride, const int offset)
{ 
  int     
    i = 0,
    k = offset;

  for(; k < totallength; i++, k += stride)    
    dest[k] = source[i];    
}




static void strideTips(char **dest, char **source, const int totallength, const int stride, 
		       const int offset, const int mxtips, int strideLength)
{ 
  int i, j, k;

  assert(offset < stride);

  for(i = 0; i < mxtips; i++)
    {      
      char *d   = &dest[i + 1][0];
      char *s   = &source[i + 1][0];      

      for(k = 0, j = offset; j < totallength; j += stride, k++)
	{
	  assert(k < strideLength);
	  d[k] = s[j];	
	}
    }
 
}

static void strideInt(int *dest, int *source, const int totallength, const int stride, const int offset)
{ 
  int i, k,
    *d = &dest[0];

  for(i = offset, k = 0; i < totallength; i += stride, k++) 
    d[k] = source[i];        
}

static void strideDouble(double *dest, double *source, const int totallength, const int stride, const int offset)
{ 
  int     i = offset;
  double *d = dest;

  for(; i < totallength; i += stride, d++)    
    *d = source[i];    
}


static void stridePartitionData(tree *localTree, int tid, int n, int length)
{
  int 
    i,
    endsite,
    dummy;

  strided_Bounds(tid, length,  n, &dummy, &endsite);

  /* printf("%d %d\n", endsite, tid); */

  for(i = 0; i < localTree->NumberOfModels; i++)
    {
      localTree->strided_partitionData[i].dataType   = localTree->partitionData[i].dataType;
      localTree->strided_partitionData[i].protModels = localTree->strided_partitionData[i].protModels;
      localTree->strided_partitionData[i].protFreqs  = localTree->strided_partitionData[i].protFreqs;     
    }

  if(localTree->NumberOfModels > 1)
    {
      int i, model;     

      localTree->strided_partitionData[0].lower = 0;
     
      model = localTree->strided_model[0];
      i = 1;

      while(i < endsite)
	{
	  if(localTree->strided_model[i] != model)
	    {	      
	      localTree->strided_partitionData[model].upper = i;
	      localTree->strided_partitionData[model + 1].lower = i;
	      model = localTree->strided_model[i];
	    }
	  i++;
	}

      
      localTree->strided_partitionData[localTree->NumberOfModels - 1].upper = endsite;

      /* 
	 for(i = 0; i < localTree->NumberOfModels; i++)
	 printf("%d %d %d\n", tid,  localTree->strided_partitionData[i].lower, localTree->strided_partitionData[i].upper);
      */
    }
  else
    {     
      localTree->strided_partitionData[0].lower = 0;
      localTree->strided_partitionData[0].upper = endsite;     
    }

}


inline static void sendTraversalInfo(tree *localTree, tree *tr)
{
  /* the one below is a hack we are re-assigning the local pointer to the global one
     the memcpy version below is just for testing and preparing the 
     fine-grained MPI BlueGene version */

  if(1)
    {
      localTree->td[0] = tr->td[0];
    }
  else
    {     
      localTree->td[0].count = tr->td[0].count;     
      memcpy(localTree->td[0].ti, tr->td[0].ti, localTree->td[0].count * sizeof(traversalInfo));
    }
}

#endif




static void execFunction(tree *tr, tree *localTree, const int startIndex, const int endIndex, 
			 const int parsimonyStartIndex, const int parsimonyEndIndex, int tid, int n)
{
  double result, dlnLdlz, d2lnLdlz2;
  int parsimonyResult;   

  /* new */
  int currentJob;
  
  currentJob = threadJob >> 16;

  /* new */
  switch(currentJob)      
    /*switch(threadJob) */ 
    {
    case THREAD_NEWVIEW:	
#ifdef _LOCAL_DATA  
      /* send */
      sendTraversalInfo(localTree, tr);     
     
      newviewIterative(localTree, startIndex,  endIndex);
#else  
      newviewIterative(tr,        startIndex,  endIndex);
#endif
      break;

      /*****************************************************/

    case THREAD_EVALUATE:
#ifdef _LOCAL_DATA     
      /* send */
      sendTraversalInfo(localTree, tr);     
      result = evaluateIterative(localTree, startIndex,  endIndex);
#else
      result = evaluateIterative(tr,        startIndex,  endIndex);
#endif
      
      /* receive */
      reductionBuffer[tid] = result;
      break;

      /*****************************************************/

    case THREAD_SUM_MAKENEWZ:                       
#ifdef _LOCAL_DATA    
      /* send */
      sendTraversalInfo(localTree, tr);
      makenewzIterative(localTree, startIndex,  endIndex);
#else           
      makenewzIterative(tr,        startIndex,  endIndex);
#endif
     
      break;

      /*****************************************************/

    case THREAD_MAKENEWZ: 
#ifdef _LOCAL_DATA
      
      /* send */

      localTree->modelNumber = tr->modelNumber;
      localTree->coreLZ      = tr->coreLZ;

      if(localTree->multiBranch)	  	  	 		  
	execCore(localTree, &dlnLdlz, &d2lnLdlz2, localTree->strided_partitionData[localTree->modelNumber].lower, 
		 localTree->strided_partitionData[localTree->modelNumber].upper, localTree->modelNumber); 	 	
      else	
	execCore(localTree, &dlnLdlz, &d2lnLdlz2, startIndex, endIndex, localTree->modelNumber);       
#else
      if(tr->multiBranch)
	{
	  int u, l;
	  int start = tr->partitionData[tr->modelNumber].lower; 
	  int end   = tr->partitionData[tr->modelNumber].upper;
	  
	  calcBounds(tid, n, start, end, &l, &u);	  	 	   
	  
	  
	  execCore(tr, &dlnLdlz, &d2lnLdlz2, l, u, tr->modelNumber);       
	}
      else
	{
	  execCore(tr, &dlnLdlz, &d2lnLdlz2, startIndex, endIndex, tr->modelNumber); 
	}
#endif

      
      /* receive */
      reductionBuffer[tid]    = dlnLdlz;
      reductionBufferTwo[tid] = d2lnLdlz2;
      break;

   /*****************************************************/

    case THREAD_SUM_MAKENEWZ_PARTITION:
       {
	 int u, l, start, end;
	 
#ifdef _LOCAL_DATA
	 /* TODO */
	 assert(0);
	 /* only required for a rarely used, undocumented function */
#endif
	
	start = tr->partitionData[tr->modelNumber].lower;
	end   = tr->partitionData[tr->modelNumber].upper;
	calcBounds(tid, n, start, end, &l, &u);
	 	
	makenewzIterativePartition(tr, l, u, tr->modelNumber);
       }
      break;

      /*****************************************************/

    case THREAD_MAKENEWZ_PARTITION:
      {	
	int u, l, start, end;

#ifdef _LOCAL_DATA
	/* TODO */
	assert(0);
	/* only required for a rarely used, undocumented function */
#endif

	
	start = tr->partitionData[tr->modelNumber].lower; 
	end   = tr->partitionData[tr->modelNumber].upper;	 
	calcBounds(tid, n, start, end, &l, &u);	   
	
	execCorePartition(tr, &dlnLdlz, &d2lnLdlz2, l, u, tr->modelNumber); 
	reductionBuffer[tid]    = dlnLdlz;
	reductionBufferTwo[tid] = d2lnLdlz2;      
      }
      break;   

      /*****************************************************/

    case THREAD_NEWVIEW_PARTITION:     	
#ifdef _LOCAL_DATA
      /* send */
      localTree->modelNumber = tr->modelNumber;
      sendTraversalInfo(localTree, tr);
           	
      newviewIterativePartition(localTree, localTree->strided_partitionData[localTree->modelNumber].lower, 
				localTree->strided_partitionData[localTree->modelNumber].upper, localTree->modelNumber);     
#else
      {
	int u, l;
	int start = tr->partitionData[tr->modelNumber].lower; 
	int end   = tr->partitionData[tr->modelNumber].upper; 
	calcBounds(tid, n, start, end, &l, &u);
	
	newviewIterativePartition(tr, l, u, tr->modelNumber);	      
      }
#endif	     
      break;

      /*****************************************************/

    case THREAD_EVALUATE_PARTITION:      	 
#ifdef _LOCAL_DATA
      /* send */

      localTree->modelNumber =  tr->modelNumber;	     	
      sendTraversalInfo(localTree, tr);
      

     

      result = evaluateIterativePartition(localTree, localTree->strided_partitionData[localTree->modelNumber].lower, 
					  localTree->strided_partitionData[localTree->modelNumber].upper, 
					  localTree->modelNumber);      
#else
      {
	int u, l;
	int start = tr->partitionData[tr->modelNumber].lower;
	int end   = tr->partitionData[tr->modelNumber].upper;
	calcBounds(tid, n, start, end, &l, &u);	
	result = evaluateIterativePartition(tr, l, u, tr->modelNumber);
      }
#endif	

      /* receive */
      reductionBuffer[tid] = result;      
      break;
      
      /*****************************************************/

    case THREAD_RATE_CATS:		
#ifdef _LOCAL_DATA
      
      /* send */

      sendTraversalInfo(localTree, tr);           
      localTree->lower_spacing = tr->lower_spacing;
      localTree->upper_spacing = tr->upper_spacing;           

      optRateCat_LOCAL(localTree, startIndex, endIndex, 
		       localTree->lower_spacing, localTree->upper_spacing, localTree->strided_lhs);       
      
      /* receive */

      collectDouble(tr->cdta->patrat,       localTree->strided_patrat,       tr->cdta->endsite, n, tid);    
      collectDouble(tr->cdta->patratStored, localTree->strided_patratStored, tr->cdta->endsite, n, tid);      
      collectDouble(tr->lhs,                localTree->strided_lhs,          tr->cdta->endsite, n, tid);        
#else
      {
	int i;
	for(i = 0; i < tr->cdta->endsite; i++)
	  if(i % n == tid)
	    optRateCat(tr, i, tr->lower_spacing, tr->upper_spacing, tr->lhs); 
      }
#endif      
      break;

      /*****************************************************/

    case THREAD_NEWVIEW_PARSIMONY:     

#ifdef _LOCAL_DATA    
      /* send */
      sendTraversalInfo(localTree, tr);  
      newviewParsimonyIterative(localTree, parsimonyStartIndex, parsimonyEndIndex); 
#else      
      newviewParsimonyIterative(tr,        parsimonyStartIndex, parsimonyEndIndex); 
#endif
     
      break;

      /*****************************************************/

    case THREAD_EVALUATE_PARSIMONY:      	

#ifdef _LOCAL_DATA 
      /* send */
      sendTraversalInfo(localTree, tr);   
      parsimonyResult = evaluateParsimonyIterative(localTree, parsimonyStartIndex, parsimonyEndIndex);
#else          
      parsimonyResult = evaluateParsimonyIterative(tr,        parsimonyStartIndex, parsimonyEndIndex);
#endif
      
      /* receive */
      reductionBufferParsimony[tid] = parsimonyResult;
      
      break;

      /*****************************************************/

    case THREAD_EVALUATE_VECTOR:         
#ifdef _LOCAL_DATA
      sendTraversalInfo(localTree, tr);
      evaluateGenericVectorIterative(localTree, startIndex, endIndex);

      collectDouble(tr->siteLL_Vector, localTree->strided_siteLL_Vector,       tr->cdta->endsite, n, tid);
      /* TODO */
      /* assert(0);*/
      /* rarely used function */
#else

      evaluateGenericVectorIterative(tr, startIndex, endIndex);
#endif
      break;

      /*****************************************************/

    case THREAD_CATEGORIZE:          
#ifdef _LOCAL_DATA
      {    	
	int i;

	/* send */
	sendTraversalInfo(localTree, tr);	

	for(i = 0; i < localTree->NumberOfModels; i++)
	  {      
	    localTree->strided_patrat[i * 4]     = localTree->gammaRates[i * 4];
	    localTree->strided_patrat[i * 4 + 1] = localTree->gammaRates[i * 4 + 1];
	    localTree->strided_patrat[i * 4 + 2] = localTree->gammaRates[i * 4 + 2];
	    localTree->strided_patrat[i * 4 + 3] = localTree->gammaRates[i * 4 + 3];
	    assert(i * 4 + 3 < localTree->originalCrunchedLength);
	  }
	
	localTree->NumberOfCategories = 4 * localTree->NumberOfModels;
	categorizeIterative(localTree, startIndex, endIndex);

	 for(i = startIndex; i < endIndex; i++)
	   {
	     double temp, wtemp;
	     temp = localTree->gammaRates[localTree->strided_rateCategory[i]];     
	     localTree->strided_wr[i]  = wtemp = temp * localTree->strided_aliaswgt[i];
	     localTree->strided_wr2[i] = temp * wtemp;
	   }
      }
#else
      categorizeIterative(tr, startIndex, endIndex);      
#endif
      break;

      /*****************************************************/

      
#ifdef _LOCAL_DATA                  
    case THREAD_PREPARE_PARSIMONY: 
      /*printf("THREAD_PREPARE_PARSIMONY\n"); */
      if(tid > 0)
	{      
	  localTree->parsimonyLength = tr->parsimonyLength;         
	  memcpy(localTree->partitionData ,  tr->partitionData, sizeof(pInfo) * localTree->NumberOfModels);	         	
	}
     
      strideInt(localTree->strided_model, tr->model,  
		localTree->originalCrunchedLength, n, tid);
      strideInt(localTree->strided_dataVector, tr->dataVector,  
		localTree->originalCrunchedLength, n, tid);
      stridePartitionData(localTree, tid, n, localTree->parsimonyLength);
      
      strideTips(localTree->strided_yVector, tr->yVector, localTree->originalCrunchedLength, n, tid, 
		 localTree->mxtips, localTree->strideLength);
      strideInt(localTree->strided_aliaswgt, tr->cdta->aliaswgt,  
		localTree->originalCrunchedLength, n, tid);

      localTree->mySpan          = 1 + (localTree->parsimonyLength / n);	
      localTree->parsimonyData   = (parsimonyVector *)malloc(sizeof(parsimonyVector) * 
							     localTree->mxtips * localTree->mySpan);         
      
      break;

      /*****************************************************/

    case  THREAD_FINISH_PARSIMONY:      
      free(localTree->parsimonyData);
                       
      if(tid > 0)
	{
	  localTree->cdta->endsite = tr->cdta->endsite;
	  memcpy(localTree->partitionData ,  tr->partitionData, sizeof(pInfo) * localTree->NumberOfModels); 
	}
      strideInt(localTree->strided_model, tr->model,  
		localTree->originalCrunchedLength, n, tid);
      strideInt(localTree->strided_dataVector, tr->dataVector,  
		localTree->originalCrunchedLength, n, tid);
      stridePartitionData(localTree, tid, n, localTree->cdta->endsite);
     
      strideTips(localTree->strided_yVector, tr->yVector, localTree->originalCrunchedLength, n, tid, 
		 localTree->mxtips, localTree->strideLength);
      strideInt(localTree->strided_aliaswgt, tr->cdta->aliaswgt,  
		localTree->originalCrunchedLength, n, tid);

      break;
      
      /*****************************************************/
      
    case THREAD_ALLOC_LIKELIHOOD:     
      /*printf("THREAD_ALLOC_LIKELIHOOD\n");*/
      {
	int span, i;	

	localTree->likelihoodFunction = tr->likelihoodFunction;
	localTree->currentModel       = tr->currentModel;
	localTree->cdta->endsite      = tr->cdta->endsite;
	localTree->mySpan             = 1 + (localTree->cdta->endsite / n);
		
	localTree->expArray = (int *)malloc(localTree->mySpan * localTree->mxtips * sizeof(int));	

	if(localTree->mixedData)
	  {
	    assert(0);	    	   
	  }
	else
	  {
	    switch(localTree->currentModel)
	      {	 	   
	      case M_PROTCAT:	 
		span = 20 * localTree->mySpan;	
		localTree->sumBuffer  = (double *)malloc(20 * localTree->strideLength * sizeof(double));
		break;  	        
	      case  M_PROTGAMMA:
		span = 80 * localTree->mySpan; 
		localTree->sumBuffer  = (double *)malloc(80 * localTree->strideLength * sizeof(double));		
		break;	
	      case M_GTRGAMMA:
		span = 16 * localTree->mySpan;		
		localTree->sumBuffer  = (double *)malloc(16 * localTree->strideLength * sizeof(double));
		break;
	      case M_GTRCAT:
		span = 4 * localTree->mySpan;		
		localTree->sumBuffer  = (double *)malloc(4 * localTree->strideLength * sizeof(double));
		break;	 
	      default:	
		assert(0);
	      }
	    localTree->likelihoodArray = (double *)malloc(span * localTree->mxtips * sizeof(double));
	  }
	
	
	for(i = 0; i < localTree->mxtips; i++)
	  localTree->xVector[i] = &(localTree->likelihoodArray[i * span]);
      }
      break;

      /*****************************************************/

    case THREAD_FREE_LIKELIHOOD:    
      free(localTree->expArray); 
      free(localTree->likelihoodArray);
      free(localTree->sumBuffer);
      localTree->expArray        = (int*)NULL;
      localTree->likelihoodArray = (double*)NULL;
      localTree->sumBuffer       = (double*)NULL;
      break;

      /*****************************************************/

    case THREAD_COPY_REVERSIBLE:     
      if(tid > 0)
	{
	  memcpy(localTree->tipVectorDNA, tr->tipVectorDNA, localTree->NumberOfModels * 64 * sizeof(double));  
	  memcpy(localTree->tipVectorAA,  tr->tipVectorAA,  localTree->NumberOfModels * 460 * sizeof(double)); 
	  
	  memcpy(localTree->EV_DNA, tr->EV_DNA, localTree->NumberOfModels * 16 * sizeof(double));
	  memcpy(localTree->EV_AA,  tr->EV_AA,  localTree->NumberOfModels * 400 * sizeof(double));
	  
	  memcpy(localTree->EI_DNA, tr->EI_DNA, localTree->NumberOfModels * 12 * sizeof(double));
	  memcpy(localTree->EI_AA,  tr->EI_AA,  localTree->NumberOfModels * 380 * sizeof(double));  
	  
	  memcpy(localTree->EIGN_DNA, tr->EIGN_DNA, localTree->NumberOfModels * 3 * sizeof(double));
	  memcpy(localTree->EIGN_AA,  tr->EIGN_AA,  localTree->NumberOfModels * 19  * sizeof(double));
	}   
      break;

      /*****************************************************/

    case THREAD_COPY_RATE_CATS:          
      if(tid > 0)
	localTree->NumberOfCategories = tr->NumberOfCategories;

      strideInt(localTree->strided_rateCategory, tr->cdta->rateCategory,  
		localTree->originalCrunchedLength, n, tid);

      memcpy(localTree->strided_patrat, tr->cdta->patrat, localTree->originalCrunchedLength * sizeof(double));      

      strideDouble(localTree->strided_patratStored, tr->cdta->patratStored,  
		   localTree->originalCrunchedLength, n, tid);
      strideDouble(localTree->strided_wr, tr->cdta->wr,  
		   localTree->originalCrunchedLength, n, tid);
      strideDouble(localTree->strided_wr2, tr->cdta->wr2,  
		   localTree->originalCrunchedLength, n, tid);	 	   
      
      break;
      
      /*****************************************************/

    case THREAD_COPY_GAMMA_RATES:       
      if(tid > 0)
	memcpy(localTree->gammaRates, tr->gammaRates, localTree->NumberOfModels * 4 * sizeof(double));
      break;

      /*****************************************************/

    case THREAD_COPY_INVARIANTS:     
      if(tid > 0)
	memcpy(localTree->invariants, tr->invariants, localTree->NumberOfModels * sizeof(double));
      break;

      /*****************************************************/

    case THREAD_COPY_INIT_MODEL: 
      /* printf("THREAD_COPY_INIT_MODEL\n"); */
      if(tid > 0)
	{	
	  localTree->NumberOfCategories = tr->NumberOfCategories;
	  localTree->likelihoodFunction = tr->likelihoodFunction;
	  localTree->cdta->endsite      = tr->cdta->endsite;
	  
	  memcpy(localTree->tipVectorDNA, tr->tipVectorDNA, localTree->NumberOfModels * 64 * sizeof(double));  
	  memcpy(localTree->tipVectorAA,  tr->tipVectorAA,  localTree->NumberOfModels * 460 * sizeof(double)); 
	  
	  memcpy(localTree->EV_DNA, tr->EV_DNA, localTree->NumberOfModels * 16 * sizeof(double));
	  memcpy(localTree->EV_AA,  tr->EV_AA,  localTree->NumberOfModels * 400 * sizeof(double));
	  
	  memcpy(localTree->EI_DNA, tr->EI_DNA, localTree->NumberOfModels * 12 * sizeof(double));
	  memcpy(localTree->EI_AA,  tr->EI_AA,  localTree->NumberOfModels * 380 * sizeof(double));  
	  
	  memcpy(localTree->EIGN_DNA, tr->EIGN_DNA, localTree->NumberOfModels * 3 * sizeof(double));
	  memcpy(localTree->EIGN_AA,  tr->EIGN_AA,  localTree->NumberOfModels * 19  * sizeof(double));
	  	 	  
	  memcpy(localTree->frequencies_DNA, tr->frequencies_DNA, localTree->NumberOfModels * 4 * sizeof(double));
	  memcpy(localTree->frequencies_AA,  tr->frequencies_AA,  localTree->NumberOfModels * 20  * sizeof(double));
	  
	  memcpy(localTree->invariants, tr->invariants, localTree->NumberOfModels * sizeof(double));
	  
	  memcpy(localTree->gammaRates, tr->gammaRates, localTree->NumberOfModels * 4 * sizeof(double));

	  memcpy(localTree->partitionData ,  tr->partitionData, sizeof(pInfo) * localTree->NumberOfModels);	
	}
           
      strideInt(localTree->strided_model, tr->model,  
		localTree->originalCrunchedLength, n, tid);  
      strideInt(localTree->strided_dataVector, tr->dataVector,  
		localTree->originalCrunchedLength, n, tid);
      stridePartitionData(localTree, tid, n, localTree->cdta->endsite);

      strideInt(localTree->strided_rateCategory, tr->cdta->rateCategory,  
		localTree->originalCrunchedLength, n, tid);
     
      strideInt(localTree->strided_aliaswgt, tr->cdta->aliaswgt,  
		localTree->originalCrunchedLength, n, tid);
     
      
      strideInt(localTree->strided_invariant, tr->invariant,  
		localTree->originalCrunchedLength, n, tid);

      memcpy(localTree->strided_patrat, tr->cdta->patrat, localTree->originalCrunchedLength * sizeof(double));    

      strideDouble(localTree->strided_patratStored, tr->cdta->patratStored,  
		   localTree->originalCrunchedLength, n, tid);
      strideDouble(localTree->strided_wr, tr->cdta->wr,  
		   localTree->originalCrunchedLength, n, tid);
      strideDouble(localTree->strided_wr2, tr->cdta->wr2,  
		   localTree->originalCrunchedLength, n, tid);
      
      strideTips(localTree->strided_yVector, tr->yVector, localTree->originalCrunchedLength, n, tid, 
		 localTree->mxtips, localTree->strideLength);       
      break;

      /*****************************************************/
          
    case THREAD_NEXT_REPLICATE: 
      /* printf("THREAD_NEXT_REPLICATE\n"); */
      if(tid > 0)
	{
	  localTree->cdta->endsite = tr->cdta->endsite;	       
	  memcpy(localTree->partitionData ,  tr->partitionData, sizeof(pInfo) * localTree->NumberOfModels);             
	}

      strideInt(localTree->strided_model, tr->model,  
		localTree->originalCrunchedLength, n, tid);
      strideInt(localTree->strided_dataVector, tr->dataVector,  
		localTree->originalCrunchedLength, n, tid);
      stridePartitionData(localTree, tid, n, localTree->cdta->endsite);

      strideInt(localTree->strided_aliaswgt, tr->cdta->aliaswgt,  
		localTree->originalCrunchedLength, n, tid);
      strideInt(localTree->strided_rateCategory, tr->cdta->rateCategory,  
		localTree->originalCrunchedLength, n, tid);

      strideDouble(localTree->strided_wr, tr->cdta->wr,  
		   localTree->originalCrunchedLength, n, tid);
      strideDouble(localTree->strided_wr2, tr->cdta->wr2,  
		   localTree->originalCrunchedLength, n, tid);
      
      memcpy(localTree->strided_patrat, tr->cdta->patrat, localTree->originalCrunchedLength * sizeof(double));
     
      strideTips(localTree->strided_yVector, tr->yVector, localTree->originalCrunchedLength, n, tid, 
		 localTree->mxtips, localTree->strideLength);       
          
      break;
#endif
    default:
      assert(0);
    }
}



void masterBarrier(int jobType, tree *tr) 
{
  const int n = NumberOfThreads;
  int startIndex, endIndex, i, sum,
    parsimonyStartIndex, parsimonyEndIndex; 
 
  /* new */
  jobCycle = !jobCycle;   
  threadJob = (jobType << 16) + jobCycle;

  /*
    old
    threadJob = jobType;    
    jobCycle = !jobCycle;   
  */

#ifdef _LOCAL_DATA
  strided_Bounds(0, tr->cdta->endsite,   n, &startIndex, &endIndex);
  strided_Bounds(0, tr->parsimonyLength, n, &parsimonyStartIndex, &parsimonyEndIndex); 
#else
  calcBounds(0, n, 0, tr->parsimonyLength, &parsimonyStartIndex, &parsimonyEndIndex);
  calcBounds(0, n, 0, tr->cdta->endsite, &startIndex, &endIndex);  
#endif
  
  execFunction(tr, tr, startIndex, endIndex, parsimonyStartIndex, parsimonyEndIndex, 0, n);      

  do
    {     
      for(i = 1, sum = 1; i < n; i++)
	sum += barrierBuffer[i];
    }
  while(sum < n);

  for(i = 1; i < n; i++)
    barrierBuffer[i] = 0;
  /*threadJob = -1;   */
}


#ifdef _LOCAL_DATA

static void allocStrides(tree *tr)
{
  int i; 

  if(tr->numBranches < NUM_BRANCHES)
    {
      printf("PERFORMANCE WARNING: for optimal efficiency on this dataset\n");
      printf("set NUM_BRANCHES to %d in file axml.h an re-compile\n", tr->numBranches);     
    }

  tr->strideLength =  1 + (tr->originalCrunchedLength / NumberOfThreads);  

  tr->strided_y0         = (char *)malloc(tr->strideLength * tr->mxtips * sizeof(char));
  tr->strided_yVector    = (char **)malloc((tr->mxtips + 1) * sizeof(char *));
  
  for(i = 0; i <  tr->mxtips; i++)
    tr->strided_yVector[i + 1] = &(tr->strided_y0[tr->strideLength * i]);

  tr->strided_aliaswgt      = (int *)malloc(sizeof(int) *  tr->strideLength);
  tr->strided_invariant     = (int *)malloc(sizeof(int) *  tr->strideLength);
  tr->strided_model         = (int *)malloc(sizeof(int) *  tr->strideLength);
  tr->strided_rateCategory  = (int *)malloc(sizeof(int) *  tr->strideLength);
  tr->strided_dataVector    = (int *)malloc(sizeof(int) *  tr->strideLength);

  tr->strided_wr            = (double *)malloc(sizeof(double) *  tr->strideLength);
  tr->strided_wr2           = (double *)malloc(sizeof(double) *  tr->strideLength);   
  tr->strided_siteLL_Vector = (double *)malloc(sizeof(double) *  tr->strideLength);

  /* this is a bit ugly here */

  tr->strided_patrat       = (double *)malloc(sizeof(double) *  tr->originalCrunchedLength);



  tr->strided_patratStored = (double *)malloc(sizeof(double) *  tr->strideLength);
  tr->strided_lhs          = (double *)malloc(sizeof(double) *  tr->strideLength);

  tr->strided_partitionData = (pInfo*)malloc(sizeof(pInfo) * tr->NumberOfModels); 
}

#endif

static void *likelihoodThread(void *tData)
{ 
  threadData *td = (threadData*)tData; 
  tree *tr = td->tr; 
  tree *localTree = (tree *)malloc(sizeof(tree));

  int 
    parsimonyStartIndex, 
    parsimonyEndIndex,
    startIndex, 
    endIndex, 
    myCycle = 0;

  const int n = NumberOfThreads;
  const int tid             = td->threadNumber;

#ifdef _LOCAL_DATA
  cruncheddata *cdta = (cruncheddata *)malloc(sizeof(cruncheddata));
#endif 	  

#ifndef _MAC
  pinThread2Cpu(tid);
#endif

#ifdef _LOCAL_DATA 
  localTree->expArray        = (int*)NULL;
  localTree->likelihoodArray = (double*)NULL;
  localTree->sumBuffer       = (double*)NULL;
  localTree->cdta = cdta;  
  localTree->mixedData               = tr->mixedData;
  localTree->NumberOfModels          = tr->NumberOfModels;
  localTree->mxtips                  = tr->mxtips;    
  localTree->originalCrunchedLength  = tr->originalCrunchedLength;  
  localTree->multiBranch             = tr->multiBranch;
  localTree->numBranches             = tr->numBranches;

  localTree->tipVectorDNA    = (double *)malloc(localTree->NumberOfModels * 64 * sizeof(double));  
  localTree->tipVectorAA     = (double *)malloc(localTree->NumberOfModels * 460 * sizeof(double)); 
  
  localTree->EV_DNA          = (double *)malloc(localTree->NumberOfModels * 16 * sizeof(double));
  localTree->EV_AA           = (double *)malloc(localTree->NumberOfModels * 400 * sizeof(double));
  
  localTree->EI_DNA          = (double *)malloc(localTree->NumberOfModels * 12 * sizeof(double));
  localTree->EI_AA           = (double *)malloc(localTree->NumberOfModels * 380 * sizeof(double));  
  
  localTree->EIGN_DNA        = (double *)malloc(localTree->NumberOfModels * 3 * sizeof(double));
  localTree->EIGN_AA         = (double *)malloc(localTree->NumberOfModels * 19  * sizeof(double));
  
  localTree->frequencies_DNA = (double *)malloc(localTree->NumberOfModels * 4 * sizeof(double));
  localTree->frequencies_AA  = (double *)malloc(localTree->NumberOfModels * 20  * sizeof(double));
  
  localTree->initialRates_DNA = (double *)malloc(localTree->NumberOfModels * 5 * sizeof(double));
  localTree->initialRates_AA  = (double *)malloc(localTree->NumberOfModels * 190 * sizeof(double));
  
  localTree->xVector      = (double **)malloc(localTree->mxtips * sizeof(double *)); 
  
  localTree->gammaRates    = (double *)malloc(localTree->NumberOfModels * 4 * sizeof(double)); 
  localTree->invariants    = (double *)malloc(localTree->NumberOfModels * sizeof(double));
  localTree->model         = (int *)   malloc(localTree->originalCrunchedLength * sizeof(int)); 
 
  localTree->partitionData = (pInfo*)malloc(sizeof(pInfo) * localTree->NumberOfModels);  
  
  localTree->td[0].count = 0;
  localTree->td[0].ti    = (traversalInfo *)malloc(sizeof(traversalInfo) * localTree->mxtips);
 
  localTree->NumberOfCategories = tr->NumberOfCategories; 

  allocStrides(localTree);
 
 
#endif  

  printf("\nThis is RAxML Worker Pthread Number: %d\n", tid);
   
  while(1)
    {           
      /* 
	 old 
	 while (myCycle == jobCycle);
	 myCycle = !myCycle;
      */

      /* new */
      while (myCycle == threadJob);
      myCycle = threadJob;

#ifdef _LOCAL_DATA     
      strided_Bounds(tid, localTree->cdta->endsite,   n, &startIndex, &endIndex);
      strided_Bounds(tid, localTree->parsimonyLength, n, &parsimonyStartIndex, &parsimonyEndIndex);       
#else     
      calcBounds(tid, n, 0, tr->cdta->endsite, &startIndex, &endIndex);              
      calcBounds(tid, n, 0, tr->parsimonyLength, &parsimonyStartIndex, &parsimonyEndIndex);
#endif

      execFunction(tr, localTree, startIndex, endIndex, parsimonyStartIndex, parsimonyEndIndex, tid, n);
     
      barrierBuffer[tid] = 1;
    }

  return (void*)NULL;
}

static void startPthreads(tree *tr)
{  
  pthread_t *threads;
  pthread_attr_t attr;
  int rc, t;
  threadData *tData;


  /* TODO pthread_attr_getstackaddr and pthread_attr_setstackaddr */
  
  jobCycle        = 0;
  /* old */
  /* threadJob       = -1; */
  /* new */
  threadJob       = 0;
  
  printf("\nThis is the RAxML Master Pthread\n");
  
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_DETACHED);
  
  threads    = (pthread_t *)malloc(NumberOfThreads * sizeof(pthread_t));
  tData      = (threadData *)malloc(NumberOfThreads * sizeof(threadData));
  reductionBuffer          = (double *)malloc(sizeof(double) *  NumberOfThreads);
  reductionBufferTwo       = (double *)malloc(sizeof(double) *  NumberOfThreads);
  reductionBufferParsimony = (int *)malloc(sizeof(int) *  NumberOfThreads);
  barrierBuffer            = (int *)malloc(sizeof(int) *  NumberOfThreads);
  
#ifdef _LOCAL_DATA
  allocStrides(tr);
#endif

  for(t = 0; t < NumberOfThreads; t++)
    barrierBuffer[t] = 0;      

#ifndef _MAC
  pinThread2Cpu(0);
#endif

  for(t = 1; t < NumberOfThreads; t++)
    {
      tData[t].tr  = tr;
      tData[t].threadNumber = t;
      rc = pthread_create(&threads[t], &attr, likelihoodThread, (void *)(&tData[t]));
      if(rc)
	{
	  printf("ERROR; return code from pthread_create() is %d\n", rc);
	  exit(-1);
	}
    }
}



#endif


/*************************************************************************************************************************************************************/

typedef struct {
  double lh;
  int tree;
  double weight;
} elw;

static int elwCompare(const void *p1, const void *p2)
{
  elw *rc1 = (elw *)p1;
  elw *rc2 = (elw *)p2;
  
  double i = rc1->weight;
  double j = rc2->weight;
  
  if (i > j)
    return (-1);
  if (i < j)
    return (1);
  return (0);
}






static void computeELW(tree *tr, analdef *adef, char *bootStrapFileName)
{
  int 
    numberOfTrees = 0,   
    i, k;
  char ch; 

  /* 
     double 
     bestLH = unlikely,
     elwSum = 0.0;     
  */

  FILE *infoFile;
  int *originalRateCategories = (int*)malloc(tr->cdta->endsite * sizeof(int));      
  int *originalInvariant      = (int*)malloc(tr->cdta->endsite * sizeof(int));
  long startSeed;           
  double **lhs;
  double **lhweights;
  elw *bootweights;

  infoFile = fopen(infoFileName, "a");   

  initModel(tr, tr->rdta, tr->cdta, adef); 
  allocNodex(tr, adef); 

  INFILE = fopen(bootStrapFileName, "r");       
  while((ch = getc(INFILE)) != EOF)
    {
      if(ch == ';')
	numberOfTrees++;
    }	 
  rewind(INFILE);

  if(numberOfTrees < 2)
    {
      printf("Error, there is only one tree in file %s which you want to use to conduct an ELW test\n", bootStrapFileName);

      exit(-1);
    }
  
  printf("\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);
  fprintf(infoFile, "\n\nFound %d trees in File %s\n\n", numberOfTrees, bootStrapFileName);

  bootweights = (elw *)malloc(sizeof(elw) * numberOfTrees);

  lhs = (double **)malloc(sizeof(double *) * numberOfTrees);

  for(k = 0; k < numberOfTrees; k++)
    lhs[k] = (double *)malloc(sizeof(double) * adef->multipleRuns);

  lhweights = (double **)malloc(sizeof(double *) * numberOfTrees);

  for(k = 0; k < numberOfTrees; k++)
    lhweights[k] = (double *)malloc(sizeof(double) * adef->multipleRuns);
 

  treeReadLen(INFILE, tr, adef);      
  modOpt(tr, adef);
  rewind(INFILE);

  /*
    This is for testing only !
    for(i = 0; i < numberOfTrees; i++)
    {
      treeReadLen(INFILE, tr, adef);
      treeEvaluate(tr, 2.0);
      bootweights[i].lh = tr->likelihood;
    }
    rewind(INFILE);
  */

  printf("Model optimization, first Tree: %f\n", tr->likelihood);
  fprintf(infoFile, "Model optimization, first Tree: %f\n", tr->likelihood);
  

  memcpy(originalRateCategories, tr->cdta->rateCategory, sizeof(int) * tr->cdta->endsite);
  memcpy(originalInvariant,      tr->invariant,          sizeof(int) * tr->cdta->endsite);

  assert(adef->boot > 0);
  /* this is ugly, should be passed as param to computenextreplicate() */
  startSeed = adef->boot;
  

  for(i = 0; i < numberOfTrees; i++)
    {              
      treeReadLen(INFILE, tr, adef);      
      resetBranches(tr);
      adef->rapidBoot = startSeed;     
     
      for(k = 0; k < adef->multipleRuns; k++)
	{
	  computeNextReplicate(tr, adef, originalRateCategories, originalInvariant);

	  if(k == 0)
	    treeEvaluate(tr, 2.0);
	  else
	    treeEvaluate(tr, 0.5);
	  /*printf("%d %d %f\n", i, k, tr->likelihood);*/
	  lhs[i][k] = tr->likelihood;	 
	}          

      reductionCleanup(tr, adef, originalRateCategories, originalInvariant);
    }        

  

  for(k = 0; k < adef->multipleRuns; k++)
    {
      double best = unlikely;
      double sum = 0.0;

      for(i = 0; i < numberOfTrees; i++)
	if(lhs[i][k] > best)
	  best = lhs[i][k];

      for(i = 0; i < numberOfTrees; i++)
	lhweights[i][k] = exp(lhs[i][k] - best);

      for(i = 0; i < numberOfTrees; i++)
	sum += lhweights[i][k];

      for(i = 0; i < numberOfTrees; i++)
	lhweights[i][k] = lhweights[i][k] / sum;

    }
  
  
  
  for(i = 0; i < numberOfTrees; i++)
    {
      double sum = 0.0;
      
      for(k = 0; k < adef->multipleRuns; k++)
	sum += lhweights[i][k];

      bootweights[i].weight = sum / ((double)adef->multipleRuns);     
      bootweights[i].tree   = i;
    }

  qsort(bootweights, numberOfTrees, sizeof(elw), elwCompare);


  {
    double sum = 0.0;
    
    /*printf("Tree\t Posterior Probability \t Cumulative posterior probability \t Original Likelihood\n");*/
    printf("Tree\t Posterior Probability \t Cumulative posterior probability\n");
    fprintf(infoFile, "Tree\t Posterior Probability \t Cumulative posterior probability\n");
    for(i = 0; i < numberOfTrees; i++)
      {
	 sum += bootweights[i].weight;
	 /*printf("%d\t\t %f \t\t %f \t\t\t %f\n", bootweights[i].tree, bootweights[i].weight, sum,  bootweights[i].lh);*/
	 printf("%d\t\t %f \t\t %f\n", bootweights[i].tree, bootweights[i].weight, sum); 
	 fprintf(infoFile, "%d\t\t %f \t\t %f\n", bootweights[i].tree, bootweights[i].weight, sum); 
      }
  }

  free(originalRateCategories);
  free(originalInvariant);

  fclose(INFILE); 
  fclose(infoFile); 
  exit(0);
}



static void computeDistances(tree *tr, analdef *adef)
{
  int i, j, modelCounter;
  double z0[NUM_BRANCHES];
  double result[NUM_BRANCHES];
  double t;  
  char distanceFileName[1024];

  FILE 
    *out, 
    *infoFile = fopen(infoFileName, "a");

  strcpy(distanceFileName,         workdir);   
  strcat(distanceFileName,         "RAxML_distances.");
  strcat(distanceFileName,         run_id);

  out = fopen(distanceFileName, "w");

  modOpt(tr, adef);
    
  

  printf("\nLog Likelihood Score after parameter optimization: %f\n\n", tr->likelihood); 
  printf("\nComputing pairwise ML-distances ...\n");

  fprintf(infoFile, "\nLog Likelihood Score after parameter optimization: %f\n\n", tr->likelihood); 
  fprintf(infoFile, "\nComputing pairwise ML-distances ...\n");
  
  for(modelCounter = 0; modelCounter < tr->NumberOfModels; modelCounter++)
    z0[modelCounter] = defaultz;

  t = gettime();

  for(i = 1; i <= tr->mxtips; i++)
    for(j = i + 1; j <= tr->mxtips; j++)
      {        
	double z, x;
	
	makenewzGenericDistance(tr, 10, z0, result, i, j);

	if(tr->multiBranch) 
	  {
	    int k;
	   	    
	    for(k = 0, x = 0.0; k < tr->numBranches; k++)
	      {
		assert(tr->partitionContributions[k] != -1.0);
		assert(tr->fracchanges[k] != -1.0);
		z = result[k];
		if (z < zmin) 
		  z = zmin;      	 		
		x += (-log(z) * tr->fracchanges[k]) * tr->partitionContributions[k];
	      }	
	  }
	else
	  {
	    z = result[0];
	    if (z < zmin) 
	      z = zmin;      	 
	    x = -log(z) * tr->fracchange;
	  }
	      
	/*printf("%s-%s \t %f\n", tr->nameList[i], tr->nameList[j], x);*/
	fprintf(out, "%s-%s \t %f\n", tr->nameList[i], tr->nameList[j], x);
      }

  fclose(out);

  t = gettime() - t;

  printf("\nTime for pair-wise ML distance computation of %d distances: %f seconds\n", 
	 (tr->mxtips * tr->mxtips - tr->mxtips) / 2, t);
  printf("\nDistances written to file: %s\n", distanceFileName);

  fprintf(infoFile, "\nTime for pair-wise ML distance computation of %d distances: %f seconds\n", 
	  (tr->mxtips * tr->mxtips - tr->mxtips) / 2, t);
  fprintf(infoFile, "\nDistances written to file: %s\n", distanceFileName);
  
  fclose(infoFile);

  exit(0);
}

int main (int argc, char *argv[]) 
{   
  rawdata      *rdta;
  cruncheddata *cdta;
  tree         *tr;         
  analdef      *adef;  
   
#ifdef _USE_OMP
  {
    int tid;
#pragma omp parallel private(tid)
    {
      tid = omp_get_thread_num();      
  
      if (tid == 0) 
	{
	  NumberOfThreads = omp_get_num_threads();
	  printf("\nRAxML with OpenMP: Number of threads = %d\n\n",  NumberOfThreads);
	}
    }  
  }
#endif


#ifdef PARALLEL
  MPI_Init(&argc, &argv); 
  MPI_Comm_rank(MPI_COMM_WORLD, &processID);  
  MPI_Comm_size(MPI_COMM_WORLD, &numOfWorkers);       
  if(processID == 0)
    printf("\nThis is the RAxML MPI Master process\n");
  else
    printf("\nThis is the RAxML MPI Worker Process Number: %d\n", processID);
#else
  processID = 0;
#endif
 
  masterTime = gettime();            

  adef = (analdef *)malloc(sizeof(analdef));
  rdta = (rawdata *)malloc(sizeof(rawdata));
  cdta = (cruncheddata *)malloc(sizeof(cruncheddata));     
  tr   = (tree *)malloc(sizeof(tree));

  initAdef(adef); 
  get_args(argc,argv, adef, tr);      

  if(adef->model == M_PROTCAT || adef->model == M_GTRCAT)
    tr->rateHetModel = CAT;
  else
    {
      if(adef->useInvariant)
	tr->rateHetModel = GAMMA_I;
      else
	tr->rateHetModel = GAMMA;
    }
  
  /* 
     This is a very ugly numerical bug fix, that intends to avoid the unaesthetic phenomena
     that can occur during model param optimization due to the dependency between parameters 
     alpha and invar which are NOT independent from each other. 
     When using P-Invar set likelihood epsilon to a lower value!  

     TODO-MIX this is very ugly !

  */
          
  if(adef->useInvariant && adef->likelihoodEpsilon > 0.001)
    adef->likelihoodEpsilon = 0.001;         
  
  readData(adef, rdta, cdta, tr);   
  
  checkOutgroups(tr, adef);
  makeFileNames();    

  if(adef->useExcludeFile)
    {
      handleExcludeFile(tr, adef, rdta);
      exit(0);
    }
  
  if(adef->mode != SEQUENCE_SIMILARITY_FILTER)
    {      
 /*     checkSequences(tr, rdta, adef); */
    }
  else
    {     
      reduceBySequenceSimilarity(tr, rdta, adef);
      exit(0);
    }

  if(adef->mode == SPLIT_MULTI_GENE)
    {     
      splitMultiGene(tr, rdta);
      exit(0);
    }
  
  if(adef->mode == CHECK_ALIGNMENT)
    {
      printf("Alignment format can be read by RAxML \n");
      exit(0);
    }  
   			       
  makeweights(adef, rdta, cdta, tr);  
  makevalues(rdta, cdta, tr, adef);     
 
  if(adef->generateBS)
    {   
      generateBS(tr, adef);
      exit(0);
    }

#ifdef _USE_PTHREADS 
  startPthreads(tr);
#endif 

  if(adef->computeELW)    
    computeELW(tr, adef, bootStrapFile);    

  if(adef->boot)         
    makeboot(adef, tr);    
 
  initModel(tr, rdta, cdta, adef);                                                         
  
/*  printModelAndProgramInfo(tr, adef, argc, argv);  */
 
  if(adef->bootStopOnly > 0)
    {
      computeBootStopOnly(tr, adef, bootStrapFile);
      exit(0);
    } 

  switch(adef->mode)
    {   
    case DISTANCE_MODE:
      printf("RAxML computation of pair-wise distances\n");
      getStartingTree(tr, adef);
      computeDistances(tr, adef);
      break;
    case ARNDT_MODE:
      printf("OPT_ARNDT\n");
      getStartingTree(tr, adef);
      optimizeArndt(tr, adef);
      break;
    case MEHRING_ALGO:
      getStartingTree(tr, adef);
      determineSequencePosition(tr, adef);    
      break;
    case  PARSIMONY_ADDITION:
      getStartingTree(tr, adef);
      printStartingTree(tr, adef, TRUE); 
      break;
    case OPTIMIZE_RATES:    
      if(adef->computePerSiteLLs)
	computePerSiteLLs(tr, adef, bootStrapFile);
      else
	optimizeRatesOnly(tr, adef);     
      break;
    case TREE_EVALUATION:    
      getStartingTree(tr, adef);      
      if(adef->likelihoodTest)	
	computeLHTest(tr, adef, bootStrapFile);
      else
	{	
	  modOpt(tr, adef);	
	  printLog(tr, adef, TRUE);         
	  printResult(tr, adef, TRUE);     
	  break;     
	}
    case CALC_BIPARTITIONS:          
      calcBipartitions(tr, adef, tree_file, bootStrapFile);   
      break;
    case BIG_RAPID_MODE:
      if(adef->boot)	
	doBootstrap(tr, adef, rdta, cdta);
      else
	{         
	  if(adef->rapidBoot)
	    {
#ifdef _VINCENT
	      doAllInOneVincent(tr, adef);
#else
	      doAllInOne(tr, adef);
#endif
	    }
	  else     	    	     
	    doInference(tr, adef, rdta, cdta);          			     	   
	}
      break;
    default:
      assert(0);
    }

/*  finalizeInfoFile(tr, adef); */

#ifdef PARALLEL
  MPI_Finalize();
#endif

  return 0;
}
