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


#ifdef PARALLEL
int numOfWorkers;
#endif

int processID;
infoList iList;
FILE   *INFILE;

int Thorough = 0;

char run_id[128] = "", 
  workdir[1024] = "", 
  seq_file[1024] = "", 
  tree_file[1024]="", 
  weightFileName[1024] = "", 
  modelFileName[1024] = "", 
  excludeFileName[1024] = "",
  bootStrapFile[1024] = "", 
  permFileName[1024] = "", 
  resultFileName[1024] = "", 
  logFileName[1024] = "", 
  checkpointFileName[1024] = "", 
  infoFileName[1024] = "", 
  randomFileName[1024] = "",   
  bootstrapFileName[1024] = "", 
  bipartitionsFileName[1024] = "",
  ratesFileName[1024] = "", 
  perSiteLLsFileName[1024] = "", 
  lengthFileName[1024] = "", 
  lengthFileNameModel[1024] = "",
  proteinModelFileName[1024] = "";

char *likelihood_key   = "likelihood",
  *ntaxa_key        = "ntaxa",
  *smoothed_key     = "smoothed";

char inverseMeaningDNA[16] = {'_', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', '-'};
char inverseMeaningPROT[23] = {'A','R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 
			       'T', 'W', 'Y', 'V', 'B', 'Z', '-'};
double masterTime;

const int protTipParsimonyValue[23] = {1 /* A*/, 2 /*R*/, 4/*N*/, 8/*D*/, 16/*C*/, 32/*Q*/, 64/*E*/, 128/*G*/, 
				       256/*H*/, 512 /*I*/, 1024/*L*/, 2048/*K*/, 4096/*M*/, 
				       8192 /*F*/, 16384/*P*/, 32768 /*S*/, 65535 /* T*/, 131072 /*W*/, 262144/*Y*/, 
				       524288 /*V*/, 12 /* N & D */, 96 /*Q & E*/, 1048575};

int partCount = 0;

int optimizeRatesInvocations = 1;  
int optimizeRateCategoryInvocations = 1;
int optimizeAlphaInvocations = 1;
int optimizeInvarInvocations = 1;

#ifdef _USE_OMP
volatile int             NumberOfThreads;
#endif


#ifdef _USE_PTHREADS
volatile int             jobCycle;
volatile int             threadJob;
volatile int             NumberOfThreads;
volatile double          *reductionBuffer;
volatile double          *reductionBufferTwo;
volatile int             *reductionBufferParsimony;
volatile int             *barrierBuffer;
#endif
