/* version 3.6. (c) Copyright 1993-2005 by the University of Washington.
   Written by Joseph Felsenstein, Akiko Fuseki, Sean Lamont, Andrew Keeffe,
   and Doug Buxton.
   Permission is granted to copy and use this program provided no fee is
   charged for it and provided that this copyright notice is not removed. */

#include "phylip.h"
#include "seq.h"

typedef enum {
  seqs, morphology, restsites, genefreqs
} datatype;

typedef enum {
  dna, rna, protein
} seqtype;


#ifndef OLDC
/* function prototypes */
void   getoptions(void);
void   seqboot_inputnumbers(void);
void   inputoptions(void);
char **matrix_char_new(long rows, long cols);
void   matrix_char_delete(char **mat, long rows);
void   seqboot_inputdata(void);
void   allocrest(void);
void   freerest(void);
void   allocnew(void);
void   freenew(void);
void   allocnewer(long newergroups, long newersites);
void   doinput(int argc, Char *argv[]);
void   bootweights(void);
void   writedata(long);
void   bootwrite(void);
void   freenewer(void);
/* function prototypes */
#endif

/*** Config vars ***/
/* Mutually exclusive booleans for boostrap type */
boolean bootstrap;

/* Bootstrap/jackknife sample frequency */
boolean regular = true;  /* Use 50% sampling with bootstrap/jackknife */
double fracsample = 0.5; /* ...or user-defined sample freq, [0..inf) */

/* Output format: mutually exclusive, none indicates PHYLIP */

boolean weights = false;/* Read weights file */
boolean all;             /* All alleles present in infile? */

boolean firstrep; /* TODO Must this be global? */
longer seed;

/* Filehandles and paths */
/* Usual suspects declared in phylip.c/h */
FILE *outcatfile, *outweightfile, *outmixfile, *outancfile, *outfactfile;
Char infilename[FNMLNGTH], outfilename[FNMLNGTH], catfilename[FNMLNGTH], outcatfilename[FNMLNGTH],
  weightfilename[FNMLNGTH], outweightfilename[FNMLNGTH], mixfilename[FNMLNGTH], outmixfilename[FNMLNGTH], ancfilename[FNMLNGTH], outancfilename[FNMLNGTH],
  factfilename[FNMLNGTH], outfactfilename[FNMLNGTH];
long sites, loci, maxalleles, groups, 
  nenzymes, reps, ws, blocksize, categs, maxnewsites;

datatype data;
seqtype seq;
steptr oldweight, where, how_many, mixdata, ancdata;

/* Original dataset */
/* [0..spp-1][0..sites-1] */
Char **nodep   = NULL;           /* molecular or morph data */

Char *factor = NULL;  /* factor[sites] - direct read-in of factors file */
long *factorr = NULL; /* [0..sites-1] => nondecreasing [1..groups] */

long *alleles = NULL;

/* Mapping with read-in weights eliminated
 * Allocated once in allocnew() */
long newsites;
long newgroups;
long *newwhere   = NULL;    /* Map [0..newgroups-1] => [1..newsites] */
long *newhowmany = NULL;    /* Number of chars for each [0..newgroups-1] */

/* Mapping with bootstrapped weights applied */
/* (re)allocated by allocnewer() */
long newersites, newergroups;
long *newerfactor  = NULL;  /* Map [0..newersites-1] => [1..newergroups] */
long *newerwhere   = NULL;  /* Map [0..newergroups-1] => [1..newersites] */
long *newerhowmany = NULL;  /* Number of chars for each [0..newergroups-1] */
long **charorder   = NULL;  /* Permutation [0..spp-1][0..newergroups-1] */
long **sppord      = NULL;  /* Permutation [0..newergroups-1][0..spp-1] */

void getoptions()
{
  /* interactively set options */
  long reps0;
  long inseed, inseed0, loopcount, loopcount2;
  Char ch;
  boolean done1;

  data = seqs;
  seq = dna;
  blocksize = 1;
  fracsample = 1.0;

  interleaved = true;
  loopcount = 0;
  initseed(&inseed, &inseed0, seed);
  /* These warnings only appear when user has changed an option that
   * subsequently became inapplicable */
}  /* getoptions */


void seqboot_inputnumbers()
{
  /* read numbers of species and of sites */
  long i;

  fscanf(infile, "%ld%ld", &spp, &sites);
  loci = sites;
  maxalleles = 1;
}  /* seqboot_inputnumbers */


void inputoptions()
{
  /* input the information on the options */
  long weightsum, maxfactsize, i, j, k, l, m;

  for (i = 1; i <= (sites); i++)
    factorr[i - 1] = i;
  for (i = 0; i < (sites); i++)
    oldweight[i] = 1;
  for (i = 0; i < (loci); i++)
    how_many[i] = 0;
  for (i = 0; i < (loci); i++)
    where[i] = 0; 
  for (i = 1; i <= (sites); i++) {
    how_many[factorr[i - 1] - 1]++;
    if (where[factorr[i - 1] - 1] == 0)
      where[factorr[i - 1] - 1] = i;
  }
  groups = factorr[sites - 1];
  newgroups = 0;
  newsites = 0;
  maxfactsize = 0;

  for(i = 0 ; i < loci ; i++){
    if(how_many[i] > maxfactsize){
      maxfactsize = how_many[i];
    }
  }
  maxnewsites = groups * maxfactsize;
  allocnew();
  for (i = 0; i < groups; i++) {
    if (oldweight[where[i] - 1] > 0) {
      newgroups++;
      newsites += how_many[i];
      newwhere[newgroups - 1] = where[i];
      newhowmany[newgroups - 1] = how_many[i];
    }
  }
}  /* inputoptions */


char **matrix_char_new(long rows, long cols)
{
  char **mat;
  long i;

  assert(rows > 0); assert(cols > 0);

  mat = (char **)Malloc(rows*sizeof(char *));
  for (i = 0; i < rows; i++)
    mat[i] = (char *)Malloc(cols*sizeof(char));

  return mat;
}


void matrix_char_delete(char **mat, long rows)
{
  long i;
  
  assert(mat != NULL);
  for (i = 0; i < rows; i++)
    free(mat[i]);
  free(mat);
}


double **matrix_double_new(long rows, long cols)
{
  double **mat;
  long i;

  assert(rows > 0); assert(cols > 0);
  
  mat = (double **)Malloc(rows*sizeof(double *));
  for (i = 0; i < rows; i++)
    mat[i] = (double *)Malloc(cols*sizeof(double));

  return mat;
}


void matrix_double_delete(double **mat, long rows)
{
  long i;
  
  assert(mat != NULL);
  for (i = 0; i < rows; i++)
    free(mat[i]);
  free(mat);
}


void seqboot_inputdata()
{
  /* input the names and sequences for each species */
  long i, j, k, l, m, n, basesread, basesnew=0;
  double x;
  Char charstate;
  boolean allread, done;

  nodep = matrix_char_new(spp, sites);
  j = nmlngth + (sites + (sites - 1) / 10) / 2 - 5;
  if (j < nmlngth - 1)
    j = nmlngth - 1;
  if (j > 37)
    j = 37;

  interleaved = (interleaved && ((data == seqs) || (data == restsites)));
  basesread = 0;
  allread = false;
  while (!allread) {
    /* eat white space -- if the separator line has spaces on it*/
    do {
      charstate = gettc(infile);
    } while (charstate == ' ' || charstate == '\t');
    ungetc(charstate, infile);
    if (eoln(infile)) 
      scan_eoln(infile);
    i = 1;
    while (i <= spp) {
      if ((interleaved && basesread == 0) || !interleaved)
        initname(i-1);
      j = interleaved ? basesread : 0;
      done = false;
      while (!done && !eoff(infile)) {
        if (interleaved)
          done = true;
        while (j < sites && !(eoln(infile) ||eoff(infile))) {
          charstate = gettc(infile);
          if (charstate == '\n' || charstate == '\t')
            charstate = ' ';
          if (charstate == ' ' ||
              (data == seqs && charstate >= '0' && charstate <= '9'))
            continue;
          uppercase(&charstate);
          j++;
          if (charstate == '.')
            charstate = nodep[0][j-1];
          nodep[i-1][j-1] = charstate;
        }
        if (interleaved)
          continue;
        if (j < sites) 
          scan_eoln(infile);
        else if (j == sites)
          done = true;
      }
      if (interleaved && i == 1)
        basesnew = j;
      scan_eoln(infile);
      if ((interleaved && j != basesnew) || ((!interleaved) && j != sites)){
        printf("\n\nERROR: sequences out of alignment at site %ld", j+1);
        printf(" of species %ld\n\n", i);
        exxit(-1);}
      i++;
    }
    if (interleaved) {
      basesread = basesnew;
      allread = (basesread == sites);
    } else
      allread = (i > spp);
  }
  if (!printdata)
    return;

}  /* seqboot_inputdata */


void allocrest()
{ /* allocate memory for bookkeeping arrays */

  oldweight = (steptr)Malloc(sites*sizeof(long));
  weight = (steptr)Malloc(sites*sizeof(long));
  where = (steptr)Malloc(loci*sizeof(long));
  how_many = (steptr)Malloc(loci*sizeof(long));
  factor = (Char *)Malloc(sites*sizeof(Char));
  factorr = (steptr)Malloc(sites*sizeof(long));
  nayme = (naym *)Malloc(spp*sizeof(naym));
}  /* allocrest */

void freerest()
{
  /* Free bookkeeping arrays */
  free(oldweight);
  free(weight);
  free(where);
  free(how_many);
  free(factor);
  free(factorr);
  free(nayme);
}


void allocnew(void)
{ /* allocate memory for arrays that depend on the lenght of the 
     output sequence*/
  /* Only call this function once */
  assert(newwhere == NULL && newhowmany == NULL); 
  newwhere = (steptr)Malloc(loci*sizeof(long));
  newhowmany = (steptr)Malloc(loci*sizeof(long));
} /* allocnew */


void freenew(void)
{ /* free arrays allocated by allocnew() */
  /* Only call this function once */
  assert(newwhere != NULL);
  assert(newhowmany != NULL);

  free(newwhere);
  free(newhowmany);
}


void allocnewer(long newergroups, long newersites)
{ /* allocate memory for arrays that depend on the length of the bootstrapped
     output sequence */
  /* Assumes that spp remains constant */
  static long curnewergroups = 0;
  static long curnewersites  = 0;

  long i;

  if (newerwhere != NULL) {
    if (newergroups > curnewergroups) {
      free(newerwhere);
      free(newerhowmany);
      for (i = 0; i < spp; i++)
        free(charorder[i]);
      newerwhere = NULL;
    }
    if (newersites > curnewersites) {
      free(newerfactor);
      newerfactor = NULL;
    }
  }

  if (charorder == NULL)
    charorder = (steptr *)Malloc(spp*sizeof(steptr));

  /* Malloc() will fail if either is 0, so add a dummy element */
  if (newergroups == 0)
    newergroups++;
  if (newersites == 0)
    newersites++;
  
  if (newerwhere == NULL) {
    newerwhere = (steptr)Malloc(newergroups*sizeof(long));
    newerhowmany = (steptr)Malloc(newergroups*sizeof(long));
    for (i = 0; i < spp; i++)
      charorder[i] = (steptr)Malloc(newergroups*sizeof(long));
    curnewergroups = newergroups;
  }
  if (newerfactor == NULL) {
    newerfactor = (steptr)Malloc(newersites*sizeof(long));
    curnewersites = newersites;
  }
}


void freenewer()
{
  /* Free memory allocated by allocnewer() */
  /* spp must be the same as when allocnewer was called */
  long i;

  if (newerwhere) {
    free(newerwhere);
    free(newerhowmany);
    free(newerfactor);
    for (i = 0; i < spp; i++)
      free(charorder[i]);
    free(charorder);
  }
}


void doinput(int argc, Char *argv[])
{ /* reads the input data */
  getoptions();
  seqboot_inputnumbers();
  allocrest();
  inputoptions();
  seqboot_inputdata();
}  /* doinput */


void bootweights()
{ /* sets up weights by resampling data */
  long i, j, k, blocks;
  double p, q, r;
  long grp = 0, site = 0;

  ws = newgroups;
  for (i = 0; i < (ws); i++)
    weight[i] = 0;
    blocks = fracsample * newgroups / blocksize;
    for (i = 1; i <= (blocks); i++) {
      j = (long)(newgroups * randum(seed)) + 1;
      for (k = 0; k < blocksize; k++) {
        weight[j - 1]++;
        j++;
        if (j > newgroups)
          j = 1;
      }
    }

  /* Count number of replicated groups */
  newergroups = 0;
  newersites  = 0;
  for (i = 0; i < newgroups; i++) {
    newergroups += weight[i];
    newersites  += newhowmany[i] * weight[i];
  }

  if (newergroups < 1) {
    fprintf(stdout, "ERROR: sampling frequency or number of sites is too small\n");
    exxit(-1);
  }
  
  /* reallocate "newer" arrays, sized by output groups:
   * newerfactor, newerwhere, newerhowmany, and charorder */
  allocnewer(newergroups, newersites);
  
  /* Replicate each group i weight[i] times */
  grp = 0;
  site = 0;
  for (i = 0; i < newgroups; i++) {
    for (j = 0; j < weight[i]; j++) {
      for (k = 0; k < newhowmany[i]; k++) {
        newerfactor[site] = grp + 1;
        site++;
      }
      newerwhere[grp] = newwhere[i];
      newerhowmany[grp] = newhowmany[i];
      grp++;
    }
  }
}  /* bootweights */

void writedata(long bs)
{
  /* write out one set of bootstrapped sequences */
  long i, j, k, l, m, n, n2;
  double x;
  FILE *outputfile;
  Char charstate;
  char filename[20];
  char bs_id[10];
  sprintf(bs_id,"%ld",bs);

  strcpy(filename, "bs.");  
  strcat(filename, bs_id);
  outputfile = fopen(filename,"w");
  sppord = (long **)Malloc(newergroups*sizeof(long *));
  for (i = 0; i < (newergroups); i++)
    sppord[i] = (long *)Malloc(spp*sizeof(long));
  for (j = 1; j <= spp; j++)
    sppord[0][j - 1] = j;
  for (i = 1; i < newergroups; i++) {
    for (j = 1; j <= (spp); j++)
      sppord[i][j - 1] = sppord[i - 1][j - 1];
  }
  fprintf(outputfile, "%5ld %5ld\n", spp, newersites);
  l = 1;
  /* When rewriting to PHYLIP, only convert interleaved <-> sequential
   * for molecular and restriction sites. */
  m = interleaved ? 60 : newergroups;
  do {
    if (m > newergroups)
      m = newergroups;
    for (j = 0; j < spp; j++) {
      n = 0;
      if ((l == 1)) {
        n2 = nmlngth;
        for (k = 0; k < n2; k++)
            putc(nayme[j][k], outputfile);
      }
      for (k = 0; k < nmlngth-n2; k++)
        fprintf(outputfile, " "); 
      fprintf(outputfile, " "); 
      for (k = l - 1; k < m; k++) {
        for (n2 = -1; n2 <= (newerhowmany[charorder[j][k]] - 2); n2++) {
          n++;
            charstate = nodep[sppord[charorder[j][k]][j] - 1]
                             [newerwhere[charorder[j][k]] + n2];
            putc(charstate, outputfile);
            if (n % 10 == 0 && n % 60 != 0)
              putc(' ', outputfile);
        }
      }
      putc('\n', outputfile);
    }
    if (interleaved) {
      if ((m <= newersites) && (newersites > 60))
        putc('\n', outputfile);
      l += 60;
      m += 60;
    }
  } while (interleaved && l <= newersites);
  for (i = 0; i < (newergroups); i++)
    free(sppord[i]);
  free(sppord);
  FClose(outputfile);
}  /* writedata */


void bootwrite()
{ /* does bootstrapping and writes out data sets */
  long i, j, rr, repdiv10;

  repdiv10 = reps / 10;
  if (repdiv10 < 1)
    repdiv10 = 1;
  firstrep = true;
  for (rr = 1; rr <= (reps); rr++) {
    bootweights();
    for (i = 0; i < spp; i++)
      for (j = 0; j < newergroups; j++)
        charorder[i][j] = j;
   
    writedata(rr);
  }
}  /* bootwrite */


int main(int argc, Char *argv[])
{  /* Read in sequences or frequencies and bootstrap or jackknife them */
  init(argc,argv);
  openfile(&infile, argv[1], "input file", "r", argv[0], infilename);
  reps = atol(argv[2]);
  doinput(argc, argv);
  bootwrite();
  
  freenewer();
  freenew();
  freerest();

  if (nodep)
    matrix_char_delete(nodep, spp);

  FClose(infile);
  
  return 0;
}
