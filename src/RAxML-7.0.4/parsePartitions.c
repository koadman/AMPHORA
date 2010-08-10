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

/*****************************FUNCTIONS FOR READING MULTIPLE MODEL SPECIFICATIONS************************************************/


extern char modelFileName[1024];
extern char excludeFileName[1024];
extern char proteinModelFileName[1024];


extern char inverseMeaningDNA[16];
extern char inverseMeaningPROT[23];
extern char seq_file[1024];

static boolean lineContainsOnlyWhiteChars(char *line)
{
  int i, n = strlen(line);

  if(n == 0)
    return TRUE;

  for(i = 0; i < n; i++)
    {
      if(!whitechar(line[i]))
	return FALSE;
    }
  return TRUE;
}


static int isNum(char c)
{
  
  return (c == '0' || c == '1' || c == '2' || c == '3' || c == '4' ||
	  c == '5' || c == '6' || c == '7' || c == '8' || c == '9');
}


static void skipWhites(char **ch)
{
  while(**ch == ' ' || **ch == '\t')
    *ch = *ch + 1;
}

static void analyzeIdentifier(char **ch, int modelNumber, tree *tr)
{
  char ident[2048] = "";
  char model[128] = "";
  char *protModels[10] = {"DAYHOFF", "DCMUT", "JTT", "MTREV", "WAG", "RTREV", "CPREV", "VT", "BLOSUM62", "MTMAM"};
  char *dnaModels[1]   = {"DNA"};
  char thisModel[1024];
  int i = 0, n, j;
  int containsComma = 0;

  while(**ch != '=')
    {
      if(**ch != ' ' && **ch != '\t')
	{
	  ident[i] = **ch;      
	  i++;
	}
      *ch = *ch + 1;
    }
  
  n = i;
  i = 0;
  
  for(i = 0; i < n; i++)
    if(ident[i] == ',') 
      containsComma = 1;

  if(!containsComma)
    {
      printf("Error, model file must have format: DNA or AA model, then a comma, and then the partition name\n");
      exit(-1);
    }
  else
    {
      boolean found = FALSE;
      i = 0;
      while(ident[i] != ',')
	{
	  model[i] = ident[i];
	  i++;
	}      
      
      /* AA */

      for(i = 0; i < 10 && !found; i++)
	{	
	  strcpy(thisModel, protModels[i]);
	  
	  if(strcmp(model, thisModel) == 0)
	    {
	      /*
		adef->protModels[modelNumber] = i;
		adef->protFreqs[modelNumber]  = 0;
		tr->dataType[modelNumber]     = AA_DATA;
	      */
	      
	      tr->partitionData[modelNumber].protModels = i;		  
	      tr->partitionData[modelNumber].protFreqs  = 0;
	      tr->partitionData[modelNumber].dataType   = AA_DATA;
	      found = TRUE;
	    }
	  	  
	  strcpy(thisModel, protModels[i]);
	  strcat(thisModel, "F");
	  
	  if(strcmp(model, thisModel) == 0)
	    {
	      /*
		adef->protModels[modelNumber] = i;
		adef->protFreqs[modelNumber]  = 1;
		tr->dataType[modelNumber]     = AA_DATA;
	      */

	      tr->partitionData[modelNumber].protModels = i;		  
	      tr->partitionData[modelNumber].protFreqs  = 1;
	      tr->partitionData[modelNumber].dataType   = AA_DATA;
	      found = TRUE;
	    }	  	 	  
	}
      
      if(!found)
	{	
	  strcpy(thisModel, dnaModels[0]);
	  
	  if(strcmp(model, thisModel) == 0)
	    {
	      /* 
		 adef->protModels[modelNumber] = -1;
		 adef->protFreqs[modelNumber]  = -1;
		 tr->dataType[modelNumber]     = DNA_DATA;
	      */
	      
	      tr->partitionData[modelNumber].protModels = -1;		  
	      tr->partitionData[modelNumber].protFreqs  = -1;
	      tr->partitionData[modelNumber].dataType   = DNA_DATA;
	      
	      found = TRUE;
	    }	  	  	  	 	  
	}

      if(!found)
	{
	  printf("ERROR: you specified the unknown protein or DNA model %s for partition %d\n", model, modelNumber);
	  exit(-1);
	}
           

      i = 0;
      while(ident[i++] != ',');      

      tr->partitionData[modelNumber].partitionName = (char*)malloc((n - i + 1) * sizeof(char));          

      j = 0;
      while(i < n)	
	tr->partitionData[modelNumber].partitionName[j++] =  ident[i++];

      tr->partitionData[modelNumber].partitionName[j] = '\0';                      
    }
}



static void setModel(int model, int position, int *a)
{
  if(a[position] == -1)
    a[position] = model;
  else
    {
      printf("ERROR trying to assign model %d to position %d \n", model, position);
      printf("while already model %d has been assigned to this position\n", a[position]);
      exit(-1);
    }      
}


static int myGetline(char **lineptr, int *n, FILE *stream)
{
  char *line, *p;
  int size, copy, len;
  int chunkSize = 256 * sizeof(char);

   if (*lineptr == NULL || *n < 2) 
    {
      line = (char *)realloc(*lineptr, chunkSize);
      if (line == NULL)
	return -1;
      *lineptr = line;
      *n = chunkSize;
    }

   line = *lineptr;
   size = *n;
  
   copy = size;
   p = line;
   
   while(1)
     {
       while (--copy > 0)
	 {
	   register int c = getc(stream);
	   if (c == EOF)
	     goto lose;
	   else
	     {
	       *p++ = c;
	       if(c == '\n' || c == '\r')	
		 goto win;
	     }
	 }

       /* Need to enlarge the line buffer.  */
       len = p - line;
       size *= 2;
       line = realloc (line, size);
       if (line == NULL)
	 goto lose;
       *lineptr = line;
       *n = size;
       p = line + len;
       copy = size - len;
     }
   
 lose:
  if (p == *lineptr)
    return -1;
  /* Return a partial line since we got an error in the middle.  */
 win:
  *p = '\0';
  return p - *lineptr;
}



void parsePartitions(analdef *adef, rawdata *rdta, tree *tr)
{
  FILE *f; 
  int numberOfModels = 0; 
  int nbytes = 0;
  char *ch;
  char *cc = (char *)NULL;
  char **p_names;
  int n, i, l;
  int lower, upper, modulo;
  char buf[256];
  int **partitions;
  int pairsCount;
  int as, j;
  int k; 

  f = fopen(modelFileName, "r");
  
  if (!f)
    {
      printf( "Could not open multiple model file: %s\n", modelFileName);
      exit(-1);
    }

 
  while(myGetline(&cc, &nbytes, f) > -1)
    {     
      if(!lineContainsOnlyWhiteChars(cc))
	{
	  numberOfModels++;
	}
      if(cc)
	free(cc);
      cc = (char *)NULL;
    }     
  
  rewind(f);
      
  p_names = (char **)malloc(sizeof(char *) * numberOfModels);
  partitions = (int **)malloc(sizeof(int *) * numberOfModels);
      
  /*if(adef->protModels == (int *)NULL)
    adef->protModels = (int *)malloc(sizeof(int) * numberOfModels);
    if(adef->protFreqs == (int *)NULL)
    adef->protFreqs = (int *)malloc(sizeof(int) * numberOfModels);

    tr->dataType = (int *)malloc(sizeof(int) * numberOfModels);*/
  
  tr->partitionData = (pInfo*)malloc(sizeof(pInfo) * numberOfModels);

      
  for(i = 0; i < numberOfModels; i++) 
    {
      /*
	adef->protModels[i] = adef->proteinMatrix;
	adef->protFreqs[i]  = adef->protEmpiricalFreqs;
	tr->dataType[i]     = -1;
      */

      tr->partitionData[i].protModels = adef->proteinMatrix;
      tr->partitionData[i].protFreqs  = adef->protEmpiricalFreqs;
      tr->partitionData[i].dataType   = -1;
    }

  for(i = 0; i < numberOfModels; i++)    
    partitions[i] = (int *)NULL;
    
  i = 0;
  while(myGetline(&cc, &nbytes, f) > -1)
    {          
      if(!lineContainsOnlyWhiteChars(cc))
	{
	  n = strlen(cc);	 
	  p_names[i] = (char *)malloc(sizeof(char) * (n + 1));
	  strcpy(&(p_names[i][0]), cc);
	  i++;
	}
      if(cc)
	free(cc);
      cc = (char *)NULL;
    }         

  for(i = 0; i < numberOfModels; i++)
    {           
      ch = p_names[i];     
      pairsCount = 0;
      skipWhites(&ch);
      
      if(*ch == '=')
	{
	  printf("Identifier missing prior to '=' in %s\n", p_names[i]);
	  exit(-1);
	}
      
      analyzeIdentifier(&ch, i, tr);
      ch++;
            
    numberPairs:
      pairsCount++;
      partitions[i] = (int *)realloc((void *)partitions[i], (1 + 3 * pairsCount) * sizeof(int));
      partitions[i][0] = pairsCount;
      partitions[i][3 + 3 * (pairsCount - 1)] = -1; 	
      
      skipWhites(&ch);
      
      if(!isNum(*ch))
	{
	  printf("%c Number expected in %s\n", *ch, p_names[i]);
	  exit(-1);
	}   
      
      l = 0;
      while(isNum(*ch))		 
	{
	  /*printf("%c", *ch);*/
	  buf[l] = *ch;
	  ch++;	
	  l++;
	}
      buf[l] = '\0';
      lower = atoi(buf);
      partitions[i][1 + 3 * (pairsCount - 1)] = lower;   
      
      skipWhites(&ch);
      
      /* NEW */
      
      if((*ch != '-') && (*ch != ','))
	{
	  if(*ch == '\0' || *ch == '\n' || *ch == '\r')
	    {
	      upper = lower;
	      goto SINGLE_NUMBER;
	    }
	  else
	    {
	      printf("'-' or ',' expected in %s\n", p_names[i]);
	      exit(-1);
	    }
	}	 
      
      if(*ch == ',')
	{	     
	  upper = lower;
	  goto SINGLE_NUMBER;
	}
      
      /* END NEW */
      
      ch++;   
      
      skipWhites(&ch);
      
      if(!isNum(*ch))
	{
	  printf("%c Number expected in %s\n", *ch, p_names[i]);
	  exit(-1);
	}    
      
      l = 0;
      while(isNum(*ch))
	{    
	  buf[l] = *ch;
	  ch++;	
	  l++;
	}
      buf[l] = '\0';
      upper = atoi(buf);     
    SINGLE_NUMBER:
      partitions[i][2 + 3 * (pairsCount - 1)] = upper;        	  
      
      if(upper < lower)
	{
	  printf("Upper bound %d smaller than lower bound %d for this partition: %s\n", upper, lower,  p_names[i]);
	  exit(-1);
	}
      
      skipWhites(&ch);
      
      if(*ch == '\0' || *ch == '\n' || *ch == '\r') /* PC-LINEBREAK*/
	{    
	  goto parsed;
	}
      
      if(*ch == ',')
	{	 
	  ch++;
	  goto numberPairs;
	}
      
      if(*ch == '\\')
	{
	  ch++;
	  skipWhites(&ch);
	  
	  if(!isNum(*ch))
	    {
	      printf("%c Number expected in %s\n", *ch, p_names[i]);
	      exit(-1);
	    }     
	  
	  l = 0;
	  while(isNum(*ch))
	    {
	      buf[l] = *ch;
	      ch++;	
	      l++;
	    }
	  buf[l] = '\0';
	  modulo = atoi(buf);      
	  partitions[i][3 + 3 * (pairsCount - 1)] = modulo; 	
	  
	  skipWhites(&ch);
	  if(*ch == '\0' || *ch == '\n' || *ch == '\r')
	    {	     
	      goto parsed;
	    }
	  if(*ch == ',')
	    {	       
	      ch++;
	      goto numberPairs;
	    }
	}  
      
      assert(0);
       
    parsed:
      i = i;
    }
  
  fclose(f);
 
  /*********************************************************************************************************************/ 

  for(i = 0; i <= rdta->sites; i++)
    tr->model[i] = -1;
  
  for(i = 0; i < numberOfModels; i++)
    {   
      as = partitions[i][0];     
      
      for(j = 0; j < as; j++)
	{
	  lower = partitions[i][1 + j * 3];
	  upper = partitions[i][2 + j * 3]; 
	  modulo = partitions[i][3 + j * 3];	
	 
	  if(modulo == -1)
	    {
	      for(k = lower; k <= upper; k++)
		setModel(i, k, tr->model);
	    }
	  else
	    {
	      for(k = lower; k <= upper; k += modulo)
		{
		  if(k <= rdta->sites)
		    setModel(i, k, tr->model);	      
		}
	    }
	}        
    }


  for(i = 1; i < rdta->sites + 1; i++)
    {
      
      if(tr->model[i] == -1)
	{
	  printf("ERROR: Alignment Position %d has not been assigned any model\n", i);
	  exit(-1);
	}      
    }  

  for(i = 0; i < numberOfModels; i++)
    {
      free(partitions[i]);
      free(p_names[i]);
    }
  
  free(partitions);
  free(p_names);    
    
  tr->NumberOfModels = numberOfModels;     
  
  if(adef->perGeneBranchLengths)
    {
      if(tr->NumberOfModels > NUM_BRANCHES)
	{
	  printf("You are trying to use %d partitioned models for an individual per-gene branch length estimate.\n", tr->NumberOfModels);
	  printf("Currently only %d are allowed to improve efficiency.\n", NUM_BRANCHES);
	  printf("\n");
	  printf("In order to change this please replace the line \"#define NUM_BRANCHES   %d\" in file \"axml.h\" \n", NUM_BRANCHES);
	  printf("by \"#define NUM_BRANCHES   %d\" and then re-compile RAxML.\n", tr->NumberOfModels);
	  exit(-1);
	}
      else
	{
	  tr->multiBranch = 1;
	  tr->numBranches = tr->NumberOfModels;
	}
    }
}

/*******************************************************************************************************************************/

void handleExcludeFile(tree *tr, analdef *adef, rawdata *rdta)
{
  FILE *f;  
  char ch, buf[256];
  int        
    j, value, i,
    state = 0,
    numberOfModels = 0,
    l = -1,
    excludeRegion   = 0,
    excludedColumns = 0,
    modelCounter    = 1;
  int
    *excludeArray, *countArray, *modelList;
  int
    **partitions;

  printf("\n\n");

  f = fopen(excludeFileName, "r");
  
  if (!f)
    {
      printf( "Could not open multiple model file: %s\n", excludeFileName);
      exit(-1);
    }

  while((ch = getc(f)) != EOF)
    {
      if(ch == '-')
	numberOfModels++;
    } 

  excludeArray = (int*)malloc(sizeof(int) * (rdta->sites + 1));
  countArray   = (int*)malloc(sizeof(int) * (rdta->sites + 1));
  modelList    = (int *)malloc((rdta->sites + 1)* sizeof(int));

  partitions = (int **)malloc(sizeof(int *) * numberOfModels);  
  for(i = 0; i < numberOfModels; i++)
    partitions[i] = (int *)malloc(sizeof(int) * 2);

  rewind(f);
  
  while((ch = getc(f)) != EOF)
    {     
      switch(state)
	{
	case 0: /* get first number */
	  if(!whitechar(ch))
	    {
	      if(!isNum(ch))
		{
		  printf("exclude file must have format: number-number [number-number]*\n");
		  exit(-1);
		}
	      l = 0;
	      buf[l++] = ch;
	      state = 1;
	    }
	  break;
	case 1: /*get the number or detect - */
	  if(!isNum(ch) && ch != '-')
	    {
	      printf("exclude file must have format: number-number [number-number]*\n");
	      exit(-1);
	    }
	  if(isNum(ch))
	    {
	      buf[l++] = ch;
	    }
	  else
	    {
	      buf[l++] = '\0';	     
	      value = atoi(buf);
	      partitions[excludeRegion][0] = value;
	      state = 2;
	    }
	  break;
	case 2: /*get second number */
	  if(!isNum(ch))
	    {
	      printf("exclude file must have format: number-number [number-number]*\n");
	      exit(-1);
	    }
	  l = 0;
	  buf[l++] = ch;
	  state = 3;
	  break;
	case 3: /* continue second number or find end */	 
	  if(!isNum(ch) && !whitechar(ch))
	    {
	      printf("exclude file must have format: number-number [number-number]*\n");
	      exit(-1);
	    }
	  if(isNum(ch))
	    {
	      buf[l++] = ch;
	    }
	  else
	    {	      
	      buf[l++] = '\0';	     
	      value = atoi(buf);
	      partitions[excludeRegion][1] = value;
	      excludeRegion++;
	      state = 0;
	    }
	  break;
	default:
	  assert(0);
	}
    }
     
  if(state == 3)
    {
      buf[l++] = '\0';     
      value = atoi(buf);
      partitions[excludeRegion][1] = value;
      excludeRegion++;
    }
  
  assert(excludeRegion == numberOfModels);

  for(i = 0; i <= rdta->sites; i++)
    {
      excludeArray[i] = -1;
      countArray[i] = 0;      
      modelList[i] = -1;
    }  

  for(i = 0; i < numberOfModels; i++)
    {
      int lower = partitions[i][0];
      int upper = partitions[i][1];

      if(lower > upper)
	{
	  printf("Misspecified exclude region %d\n", i);
	  printf("lower bound %d is greater than upper bound %d\n", lower, upper);
	  exit(-1);
	}

      if(lower == 0)
	{
	  printf("Misspecified exclude region %d\n", i);
	  printf("lower bound must be greater than 0\n");
	  exit(-1);
	}

      if(upper > rdta->sites)
	{
	  printf("Misspecified exclude region %d\n", i);
	  printf("upper bound %d must be smaller than %d\n", upper, (rdta->sites + 1));
	  exit(-1);
	}	
      for(j = lower; j <= upper; j++)
	{
	  if(excludeArray[j] != -1)
	    {
	      printf("WARNING: Exclude regions %d and %d overlap at position %d (already excluded %d times)\n", 
		     excludeArray[j], i, j, countArray[j]);
	    }
	  excludeArray[j] = i;
	  countArray[j]   =  countArray[j] + 1;	 
	}
    }

  for(i = 1; i <= rdta->sites; i++)
    {
      if(excludeArray[i] != -1)
	excludedColumns++;
      else
	{
	  modelList[modelCounter] = tr->model[i];
	  modelCounter++;
	}
    }

  printf("You have excluded %d out of %d columns\n", excludedColumns, rdta->sites);

  if(excludedColumns == rdta->sites)
    {
      printf("Error: You have excluded all sites\n");
      exit(-1);
    }

  if(adef->useMultipleModel && (excludedColumns > 0))
    {      
      char mfn[2048];
      FILE *newFile;

      strcpy(mfn, modelFileName);
      strcat(mfn, ".");
      strcat(mfn, excludeFileName);

      newFile = fopen(mfn, "w");

      printf("\nA partition file with analogous model assignments for non-excluded columns is printed to file %s\n", mfn);     
	      
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

  
  {
    FILE *newFile;
    char mfn[2048];
   

    strcpy(mfn, seq_file);
    strcat(mfn, ".");
    strcat(mfn, excludeFileName);
    
    newFile = fopen(mfn, "w");
    
    printf("\nAn alignment file with excluded columns is printed to file %s\n\n\n", mfn);
    
    fprintf(newFile, "%d %d\n", tr->mxtips, rdta->sites - excludedColumns);
    
    for(i = 1; i <= tr->mxtips; i++)
      {   
	char *tipI =  &(rdta->y[i][1]);
	fprintf(newFile, "%s ", tr->nameList[i]);
	
	for(j = 0; j < rdta->sites; j++)
	  {
	    if(excludeArray[j + 1] == -1)
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

  
  fclose(f);
  for(i = 0; i < numberOfModels; i++)
    free(partitions[i]);
  free(partitions);  
  free(excludeArray);
  free(countArray);
  free(modelList);
}


void parseProteinModel(analdef *adef)
{
  FILE *f; 
  int doublesRead = 0;
  int result = 0;
  int i, j;
  double acc = 0.0;

  assert(adef->userProteinModel);
  printf("User-defined prot mod %s\n", proteinModelFileName);

  adef->externalAAMatrix = (double*)malloc(420 * sizeof(double));

  f = fopen(proteinModelFileName, "r");
  
  if (!f)
    {
      printf( "Could not open external AA subsitution model file: %s\n", proteinModelFileName);
      exit(-1);
    }

  while(doublesRead < 420)
    {     
      result = fscanf(f, "%lf", &(adef->externalAAMatrix[doublesRead++]));           

      if(result == EOF)
	{
	  printf("Error protein model file must consist of exactly 420 entries \n");
	  printf("The first 400 entries are for the rates of the AA matrix, while the\n");
	  printf("last 20 should contain the empirical base frequencies\n");
	  printf("Reached End of File after %d entries\n", (doublesRead - 1));
	  exit(-1);
	}    
    }
       
  fclose(f);

  /* CHECKS */
  for(i = 0; i < 20; i++)
    for(j = 0; j < 20; j++)
      {
	if(i != j)
	  {
	    if(adef->externalAAMatrix[i * 20 + j] != adef->externalAAMatrix[j * 20 + i])
	      {
		printf("Error user-defined Protein model matrix must be symmetric\n");
		printf("Entry P[%d][%d]=%f at position %d is not equal to P[%d][%d]=%f at position %d\n", 
		       i, j,  adef->externalAAMatrix[i * 20 + j], (i * 20 + j),
		       j, i,  adef->externalAAMatrix[j * 20 + i], (j * 20 + i));
		exit(-1);
	      }
	  }
      }

  acc = 0.0;

  for(i = 400; i < 420; i++)    
    acc += adef->externalAAMatrix[i];         

  if((acc > 1.0 + 1.0E-6) || (acc <  1.0 - 1.0E-6))
    {
      printf("Base frequencies in user-defined AA substitution matrix do not sum to 1.0\n");
      printf("the sum is %1.80f\n", acc);
      exit(-1);
    }

}
