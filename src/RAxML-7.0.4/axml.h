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
 *  Alexandros Stamatakis:"RAxML-VI-HPC: maximum likelihood-based phylogenetic analyses 
 *  with thousands of taxa and mixed models". 
 *  Bioinformatics 2006; doi: 10.1093/bioinformatics/btl446
 */

#include <assert.h>

#ifdef PARALLEL

#define COMPUTE_TREE 0
#define TREE         1
#define BS_TREE      2
#define ML_TREE      3
#define FINALIZE     4
#define JOB_REQUEST  5
#define PRINT_TREE   6
#define I_PRINTED_IT 7
#endif



#define smoothings     32         /* maximum smoothing passes through tree */
#define iterations     10         /* maximum iterations of iterations per insert */
#define newzpercycle   1          /* iterations of makenewz per tree traversal */
#define nmlngth        256        /* number of characters in species name */
#define deltaz         0.00001    /* test of net branch length change in update */
#define defaultz       0.9        /* value of z assigned as starting point */
#define unlikely       -1.0E300   /* low likelihood for initialization */
#define largeDouble    1.0E300    /* same as positive number */

#define SUMMARIZE_LENGTH -3
#define SUMMARIZE_LH     -2
#define NO_BRANCHES      -1

#define zmin       1.0E-15  /* max branch prop. to -log(zmin) (= 34) */
#define zmax (1.0 - 1.0E-6) /* min branch prop. to 1.0-zmax (= 1.0E-6) */
#define twotothe256  \
  115792089237316195423570985008687907853269984665640564039457584007913129639936.0
                                                     /*  2**256 (exactly)  */

#define minlikelihood  (1.0/twotothe256)
#define minusminlikelihood -minlikelihood

#define badRear         -1

#define NUM_BRANCHES   128

#define TRUE             1
#define FALSE            0

#define treeNone         0
#define treeNewick       1
#define treeProlog       2
#define treePHYLIP       3
#define treeMaxType      3
#define treeDefType  treePHYLIP

#define LIKELIHOOD_EPSILON 0.0000001

#define AA_SCALE 10.0
#define AA_SCALE_PLUS_EPSILON 10.001

/* ALPHA_MIN is critical -> numerical instability, eg for 4 discrete rate cats                    */
/* and alpha = 0.01 the lowest rate r_0 is                                                        */
/* 0.00000000000000000000000000000000000000000000000000000000000034878079110511010487             */
/* which leads to numerical problems Table for alpha settings below:                              */
/*                                                                                                */
/* 0.010000 0.00000000000000000000000000000000000000000000000000000000000034878079110511010487    */
/* 0.010000 yielded nasty numerical bugs in at least one case !                                   */
/* 0.020000 0.00000000000000000000000000000044136090435925743185910935350715027016962154188875    */
/* 0.030000 0.00000000000000000000476844846859006690412039180149775802624789852441798419292220    */
/* 0.040000 0.00000000000000049522423236954066431210260930029681736928018820007024736185030633    */
/* 0.050000 0.00000000000050625351310359203371872643495343928538368616365517027588794007897377    */
/* 0.060000 0.00000000005134625283884191118711474021861409372524676086868566926568746566772461    */
/* 0.070000 0.00000000139080650074206434685544624965062437960128249869740102440118789672851562    */
/* 0.080000 0.00000001650681201563587066858709818343436959153791576682124286890029907226562500    */
/* 0.090000 0.00000011301977332931251259273962858978301859735893231118097901344299316406250000    */
/* 0.100000 0.00000052651925834844387815526344648331402709118265192955732345581054687500000000    */


#define ALPHA_MIN    0.02 
#define ALPHA_MAX    1000.0

#define RATE_MIN     0.0000001
#define RATE_MAX     1000000.0

#define INVAR_MIN    0.0001
#define INVAR_MAX    0.9999

#define TT_MIN       0.0000001
#define TT_MAX       1000000.0

#define FREQ_MIN     0.000001 /* TO AVOID NUMERICAL PROBLEMS WHEN FREQ == 0 IN PARTITIONED MODELS, ESPECIALLY WITH AA */

#define MODEL_EPSILON 0.0001
#define ITMAX 100



#define SHFT(a,b,c,d)                (a)=(b);(b)=(c);(c)=(d);
#define SIGN(a,b)                    ((b) > 0.0 ? fabs(a) : -fabs(a))

#define ABS(x)    (((x)<0)   ?  (-(x)) : (x))
#define MIN(x,y)  (((x)<(y)) ?    (x)  : (y))
#define MAX(x,y)  (((x)>(y)) ?    (x)  : (y))
#define NINT(x)   ((int) ((x)>0 ? ((x)+0.5) : ((x)-0.5)))


#define PointGamma(prob,alpha,beta)  PointChi2(prob,2.0*(alpha))/(2.0*(beta))

#define programName        "RAxML"
#define programVersion     "7.0.4"
#define programDate        "April 2008"

#define  TREE_EVALUATION            0
#define  BIG_RAPID_MODE             1
#define  PARALLEL_MODE              2
#define  CALC_BIPARTITIONS          3
#define  SPLIT_MULTI_GENE           4
#define  CHECK_ALIGNMENT            5
#define  OPTIMIZE_RATES             6
#define  PARSIMONY_ADDITION         7
#define  SEQUENCE_SIMILARITY_FILTER 8
#define  MEHRING_ALGO               9
#define  ARNDT_MODE                 10
#define  DISTANCE_MODE              11

#define M_GTRCAT      1
#define M_GTRGAMMA    2
#define M_PROTCAT     5
#define M_PROTGAMMA   6

#define DAYHOFF    0
#define DCMUT      1
#define JTT        2
#define MTREV      3
#define WAG        4
#define RTREV      5
#define CPREV      6
#define VT         7
#define BLOSUM62   8
#define MTMAM      9
#define GTR        10

#define BITS_BYTE 8

#define GTRCAT           0
#define GTRCATMULT       1
#define PROTCAT          2
#define PROTCATMULT      3
#define GTRGAMMA         4
#define GTRGAMMAI        5
#define GTRGAMMAMULT     6
#define GTRGAMMAMULTI   7
#define PROTGAMMA        8
#define PROTGAMMAI       9
#define PROTGAMMAMULT    10
#define PROTGAMMAMULTI   11
#define PARSIMONY_PROT   12
#define PARSIMONY_DNA    13

/* bootstopping stuff */

#define BOOTSTOP_PERMUTATIONS 100
#define START_BSTOP_TEST      10

#define FC_THRESHOLD          99
#define FC_SPACING            50
#define FC_LOWER              0.99
#define FC_INIT               20

#define BC_THRESHOLD          0.05
#define BC_SPACING            100
#define BC_INIT               50

/* bootstopping stuff end */

#define AA_CAT                20
#define AA_GAMMA              80
#define DNA_CAT               4
#define DNA_GAMMA             16


#define TIP_TIP     0
#define TIP_INNER   1
#define INNER_INNER 2

#define SMALL_DATA  1
#define LARGE_DATA  2

#define BINARY_DATA 0
#define DNA_DATA    1
#define AA_DATA     2

#define DNA_RATES   5
#define AA_RATES    190

#define CAT         0
#define GAMMA       1
#define GAMMA_I     2


typedef  int boolean;


typedef struct {
  float val;
  int number;
} qtData;


typedef struct
{
  int  parsimonyScore;
  int  parsimonyState; 
} 
  parsimonyVector;


typedef struct ratec 
{
  double accumulatedSiteLikelihood;
  double rate;
} 
  rateCategorize;


typedef struct 
{  
  int tipCase;
  int pNumber;
  int qNumber;
  int rNumber;
  double qz[NUM_BRANCHES];
  double rz[NUM_BRANCHES];
} traversalInfo;
  
typedef struct 
{
  traversalInfo *ti;
  int count;
} traversalData;


/********** BOOTSTOPPING */




typedef struct {  
  int length;
  int freq1;
  int freq2;
  unsigned char *isSet;
  int  *entries; 
} bipList;

typedef struct {
  bipList *b;  
  int n;
  int count; 
  int treeVectorLength;   
  } BL;



/****************************/

typedef  struct noderec 
{
#ifdef _MULTI_GENE  
  struct noderec  *backs[NUM_BRANCHES];
  char            xs[NUM_BRANCHES];
#endif
  double           z[NUM_BRANCHES];
  struct noderec  *next;
  struct noderec  *back;
  int              number; 
  char             x;    
} 
  node, *nodeptr;
 
typedef struct 
  {
    double lh;
    int number;
  }
  info;

typedef struct bInf {
  double likelihood;  
  nodeptr node;
} bestInfo;

typedef struct iL {
  bestInfo *list;
  int n;
  int valid;
} infoList;

typedef struct {
  nodeptr p; /* not required for other lists*/
  int pNum;
  int qNum;
  int support;
  int length;/*not required for other lists*/
  int *entries;
} bList;


typedef  struct 
{
  int              numsp;      
  int              sites;      
  char             **y;  
  char             *y0;
  char             *yBUF;
  int              *wgt;        
  int              *wgt2;        
} rawdata;

typedef  struct {
  int             *alias;       /* site representing a pattern */  
  int             *aliaswgt;    /* weight by pattern */
  int             *rateCategory; 
  int              endsite;     /* # of sequence patterns */ 
  double          *patrat;      /* rates per pattern */
  double          *patratStored;
  double          *wr;          /* weighted rate per pattern */
  double          *wr2;         /* weight*rate**2 per pattern */
} cruncheddata;


typedef struct {
  int count;
  int *entries;
} qtList;


typedef struct {
  int   lower;
  int   upper;
  int   dataType;
  int   protModels;
  int   protFreqs;
  int   modelOffset;
  char *partitionName;
} pInfo;


typedef struct 
{
  nodeptr p;
  int freq;  
} insertionBranch;

typedef struct
{
  insertionBranch *ib; 
  int max;
  int reference;
} insertionPoints;

  
typedef  struct  {       
  pInfo            *partitionData;
  
  double           *likelihoodArray;
  int              *expArray;

  parsimonyVector  *parsimonyData; 

#ifdef _MULTI_GENE
  double        perPartitionLH[NUM_BRANCHES];
  int           doMulti;
  traversalData td[NUM_BRANCHES];
#else
  traversalData td[NUM_BRANCHES];
#endif    
 
  int              *dataVector;
  
  insertionBranch  *ib;
  insertionPoints  *ip;
  int numberOfTipsForInsertion;
  /* model-dependent stuff */ 

  
  
  double           *sumBuffer;
  double           *siteLL_Vector;
  double           coreLZ;
  int              modelNumber;
  int              multiBranch;
  int              numBranches;

#ifdef _LOCAL_DATA
      
  pInfo            *strided_partitionData;
  char             **strided_yVector;
  char             *strided_y0;

  
  int strideLength;
  int *strided_aliaswgt;  
  int *strided_invariant;
  int *strided_model;
  int *strided_rateCategory;
  int *strided_dataVector;

  double *strided_wr;
  double *strided_wr2;
  double *strided_patrat;
  double *strided_patratStored;
  double *strided_lhs;
  double *strided_siteLL_Vector;
#endif


#ifdef _MULTI_GENE
  nodeptr  *startVector;
  char    **tipMissing;
#endif
#ifdef _USE_PTHREADS 
  int              currentModel;
  int              mySpan;
  double           lower_spacing;
  double           upper_spacing;
  double           *lhs;   
#endif

  double           *tipVectorDNA;   
  double           *tipVectorAA;

  double           *EV_DNA;
  double           *EV_AA;

  double           *EI_DNA;
  double           *EI_AA;
  
  double           *EIGN_DNA;
  double           *EIGN_AA;

  double           *frequencies_DNA;
  double           *frequencies_AA;
  
  double           *initialRates_DNA;
  double           *initialRates_AA;
  
 
  /* the stuff below is shared among DNA and AA, span does 
     not change depending on datatype */

  double           *gammaRates;
  double           *alphas;
  double           *invariants;
  double           *fracchanges;  

  /* model stuff end */
  
  double           **xVector;
  char             **yVector;
  int              numberOfProteinPositions;
  int              numberOfNucleotidePositions;
  int              originalCrunchedLength;
  int              fullSites;
  int              *originalModel;
  int              *originalDataVector;
  int              *originalWeights;  


  double            partitionContributions[NUM_BRANCHES]; 
  double            fracchange;
  double            lhCutoff;
  double            lhAVG;
  unsigned long     lhDEC;              
  unsigned long     itCount;
  int               numberOfInvariableColumns;
  int               weightOfInvariableColumns;   
  int               likelihoodFunction;
  int               rateHetModel;



  double           startLH;
  double           endLH;
  double           likelihood;   
  double          *likelihoods;
  bList           *ML_Tree;
  int             *invariant;
  int              countML_Tree;
  int              numberOfTrees;
  node           **nodep;
  node            *start; 
  int              mxtips;
  int              *model; 

  int              *constraintVector;

  

  int              ntips;
  int              nextnode;
  int              NumberOfCategories;
  int              NumberOfModels;
  int              parsimonyLength;   
  int              checkPointCounter;
  int              treeID; 
  int              numberOfOutgroups;
  int             *outgroupNums;
  char           **outgroups;
  boolean          bigCutoff;  
  boolean          smoothed;
  boolean          rooted;
  boolean          grouped;
  boolean          constrained;
  boolean          doCutoff;
  boolean          mixedData;
  rawdata         *rdta;        
  cruncheddata    *cdta; 
  
  char **nameList;
  char *tree_string;
  int treeStringLength;
  int bestParsimony;
  double bestOfNode;
  nodeptr removeNode;
  nodeptr insertNode;

  double zqr[NUM_BRANCHES];
  double currentZQR[NUM_BRANCHES];

  double currentLZR[NUM_BRANCHES];
  double currentLZQ[NUM_BRANCHES];
  double currentLZS[NUM_BRANCHES];
  double currentLZI[NUM_BRANCHES];
  double lzs[NUM_BRANCHES];
  double lzq[NUM_BRANCHES];
  double lzr[NUM_BRANCHES];
  double lzi[NUM_BRANCHES];

} tree;


/***************************************************************/

typedef struct
{
  double z[NUM_BRANCHES];
  nodeptr p, q;
} 
  connectRELL, *connptrRELL;

typedef  struct 
{
  connectRELL     *connect;       
  int             *constraintVector;
  int             start;
  double          likelihood;
} 
  topolRELL;


typedef  struct 
{
  int max;
  int members;
  topolRELL **t; 
} 
  topolRELL_LIST;


/**************************************************************/



typedef struct conntyp {
    double           z[NUM_BRANCHES];           /* branch length */
    node            *p, *q;       /* parent and child sectors */
    void            *valptr;      /* pointer to value of subtree */
    int              descend;     /* pointer to first connect of child */
    int              sibling;     /* next connect from same parent */
    } connect, *connptr;

typedef  struct {
    double           likelihood;
  int              initialTreeNumber;
    connect         *links;       /* pointer to first connect (start) */
    node            *start;
    int              nextlink;    /* index of next available connect */
                                  /* tr->start = tpl->links->p */
    int              ntips;
    int              nextnode;
    int              scrNum;      /* position in sorted list of scores */
    int              tplNum;      /* position in sorted list of trees */
       
    boolean          smoothed;    /* branch optimization converged? */
    } topol;

typedef struct {
    double           best;        /* highest score saved */
    double           worst;       /* lowest score saved */
    topol           *start;       /* starting tree for optimization */
    topol          **byScore;
    topol          **byTopol;
    int              nkeep;       /* maximum topologies to save */
    int              nvalid;      /* number of topologies saved */
    int              ninit;       /* number of topologies initialized */
    int              numtrees;    /* number of alternatives tested */
    boolean          improved;
    } bestlist;

typedef  struct {  
  int              categories;
  int              model;
  int              bestTrav;
  int              max_rearrange;
  int              stepwidth;
  int              initial;
  boolean          initialSet;
  int              mode;
  long             boot;          
  long             rapidBoot;
  boolean          bootstrapBranchLengths;
  boolean          restart;
  boolean          useWeightFile;  
  boolean          useMultipleModel;
  boolean          constraint;
  boolean          grouping;
  boolean          randomStartingTree;
  boolean          categorizeGamma;
  boolean          useInvariant; 
#ifdef _VINCENT
  boolean          optimizeBSmodel;
#endif
  int            protEmpiricalFreqs;
  int            proteinMatrix;
  int            checkpoints;
  int            startingTreeOnly;
  int            useMixedModel;
  int            multipleRuns;
  long           parsimonySeed;
  int            multiBoot;
  boolean        reallyThoroughBoot;
  boolean        perGeneBranchLengths;
  boolean        likelihoodTest;
  boolean        outgroup;
  boolean        permuteTreeoptimize;
  boolean        allInOne;
  boolean        treeLength;
  boolean        computePerSiteLLs;
  boolean        generateBS;
  boolean        bootStopping;
  boolean        useExcludeFile;
  boolean        userProteinModel;
  boolean        rapidML_Addition;
  boolean        computeELW;
  boolean        computeDistance;
  int            bootStopOnly;
  double         likelihoodEpsilon;  
  double         sequenceSimilarity;    
  double         gapyness;
  int            similarityFilterMode;
  double         bootstopCutoff;
  double        *externalAAMatrix;
} analdef;





/****************************** FUNCTIONS ****************************************************/






extern double gettime ( void );
extern int gettimeSrand ( void );
extern double randum ( long *seed );
extern void getxnode ( nodeptr p );
extern void hookup ( nodeptr p, nodeptr q, double *z, int numBranches);
extern void hookupDefault ( nodeptr p, nodeptr q, int numBranches);
extern boolean whitechar ( int ch );
extern void makeboot ( analdef *adef, tree *tr );
extern void allocNodex ( tree *tr, analdef *adef );
extern void freeNodex ( tree *tr );
extern void errorExit ( int e );
extern void printResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printBootstrapResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printBipartitionResult ( tree *tr, analdef *adef, boolean finalPrint );
extern void printLog ( tree *tr, analdef *adef, boolean finalPrint );
extern void printStartingTree ( tree *tr, analdef *adef, boolean finalPrint );
extern void writeInfoFile ( analdef *adef, tree *tr, double t );
extern int main ( int argc, char *argv[] );
extern int countTips( nodeptr p, int numsp);
extern void calcBipartitions ( tree *tr, analdef *adef, char *bestTreeFileName, char *bootStrapFileName );
extern void makeVal ( char code, double *val );
extern void baseFrequenciesGTR ( rawdata *rdta, cruncheddata *cdta, tree *tr, analdef *adef );
extern void initReversibleGTR (tree *tr, analdef *adef, int model);
extern double LnGamma ( double alpha );
extern double IncompleteGamma ( double x, double alpha, double ln_gamma_alpha );
extern double PointNormal ( double prob );
extern double PointChi2 ( double prob, double v );
extern void makeGammaCats (int model, double *alphas, double *gammaRates);
extern void assignLikelihoodFunctions ( tree *tr, analdef *adef );
extern void initModel ( tree *tr, rawdata *rdta, cruncheddata *cdta, analdef *adef );
extern void doAllInOne ( tree *tr, analdef *adef );
#ifdef _VINCENT
extern void doAllInOneVincent( tree *tr, analdef *adef );
#endif

extern void doBootstrap ( tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta );
extern void doInference ( tree *tr, analdef *adef, rawdata *rdta, cruncheddata *cdta );
extern void resetBranches ( tree *tr );
extern void modOpt ( tree *tr, analdef *adef );
extern void modOptModel ( tree *tr, analdef *adef, int model );
extern void optimizeRateCategories ( tree *tr, int _maxCategories );
extern void specialRateCategories ( tree *tr, int _maxCategories, double modelEpsilon );
extern void quickAndDirtyOptimization ( tree *tr, analdef *adef );
extern int optimizeModel ( tree *tr, analdef *adef);
extern void optimizeAllRateCategories ( tree *tr );
extern void optimizeRatesOnly ( tree *tr, analdef *adef );
extern void parsePartitions ( analdef *adef, rawdata *rdta, tree *tr);
extern void computeBOOTRAPID (tree *tr, analdef *adef, long *radiusSeed);
extern void optimizeRAPID ( tree *tr, analdef *adef );
extern void thoroughOptimization ( tree *tr, analdef *adef, topolRELL_LIST *rl, int index );
extern int treeOptimizeThorough ( tree *tr, int mintrav, int maxtrav);

extern int checker ( tree *tr, nodeptr p );
extern int randomInt ( int n );
extern void makePermutation ( int *perm, int n, analdef *adef );
extern boolean tipHomogeneityChecker ( tree *tr, nodeptr p, int grouping );
extern void makeRandomTree ( tree *tr, analdef *adef );
extern void nodeRectifier ( tree *tr );
extern void makeParsimonyTree ( tree *tr, analdef *adef );
extern void makeParsimonyTreeIncomplete ( tree *tr, analdef *adef );

extern void tred2 ( double *a, const int n, const int np, double *d, double *e );
extern double pythag ( double a, double b );
extern void tqli ( double *d, double *e, int n, int np, double *z );
extern boolean initrav ( tree *tr, nodeptr p );
extern boolean initravDIST ( tree *tr, nodeptr p, int distance );
extern void initravPartition ( tree *tr, nodeptr p, int model );
extern boolean update ( tree *tr, nodeptr p );
extern boolean smooth ( tree *tr, nodeptr p );
extern boolean smoothTree ( tree *tr, int maxtimes );
extern boolean localSmooth ( tree *tr, nodeptr p, int maxtimes );
extern void initInfoList ( int n );
extern void freeInfoList ( void );
extern void insertInfoList ( nodeptr node, double likelihood );
extern boolean smoothRegion ( tree *tr, nodeptr p, int region );
extern boolean regionalSmooth ( tree *tr, nodeptr p, int maxtimes, int region );
extern nodeptr removeNodeBIG ( tree *tr, nodeptr p, int numBranches);
extern nodeptr removeNodeRestoreBIG ( tree *tr, nodeptr p );
extern boolean insertBIG ( tree *tr, nodeptr p, nodeptr q, int numBranches);
extern boolean insertRestoreBIG ( tree *tr, nodeptr p, nodeptr q );
extern void restoreTopologyOnly ( tree *tr, bestlist *bt );
extern boolean testInsertBIG ( tree *tr, nodeptr p, nodeptr q );
extern void addTraverseBIG ( tree *tr, nodeptr p, nodeptr q, int mintrav, int maxtrav );
extern int rearrangeBIG ( tree *tr, nodeptr p, int mintrav, int maxtrav );
extern void traversalOrder ( nodeptr p, int *count, nodeptr *nodeArray );
extern double treeOptimizeRapid ( tree *tr, int mintrav, int maxtrav, analdef *adef, bestlist *bt );
extern boolean testInsertRestoreBIG ( tree *tr, nodeptr p, nodeptr q );
extern void restoreTreeFast ( tree *tr );
extern int determineRearrangementSetting ( tree *tr, analdef *adef, bestlist *bestT, bestlist *bt );
extern void computeBIGRAPID ( tree *tr, analdef *adef );
extern void computeBIGRAPIDMULTIBOOT ( tree *tr, analdef *adef );
extern boolean treeEvaluate ( tree *tr, double smoothFactor );
extern boolean treeEvaluatePartition ( tree *tr, double smoothFactor, int model );
extern void saveTopolRELL ( tree *tr, topolRELL *tpl );
extern void restoreTopolRELL ( tree *tr, topolRELL *tpl );
extern void initTL ( topolRELL_LIST *rl, tree *tr, int n );
extern void freeTL ( topolRELL_LIST *rl , tree *tr);
extern void restoreTL ( topolRELL_LIST *rl, tree *tr, int n );
extern void resetTL ( topolRELL_LIST *rl );
extern void saveTL ( topolRELL_LIST *rl, tree *tr, int index );
extern void *tipValPtr ( nodeptr p );
extern int cmpTipVal ( void *v1, void *v2 );
extern topol *setupTopol ( int maxtips );
extern void freeTopol ( topol *tpl );
extern void saveTree ( tree *tr, topol *tpl );
extern boolean restoreTreeRecursive ( topol *tpl, tree *tr );
extern boolean restoreTree ( topol *tpl, tree *tr );
extern boolean restoreTopology ( topol *tpl, tree *tr );
extern int initBestTree ( bestlist *bt, int newkeep, int numsp );
extern void resetBestTree ( bestlist *bt );
extern boolean freeBestTree ( bestlist *bt );
extern int cmpSubtopol ( connptr p10, connptr p1, connptr p20, connptr p2 );
extern int cmpTopol ( void *tpl1, void *tpl2 );
extern int cmpTplScore ( void *tpl1, void *tpl2 );
extern int findInList ( void *item, void *list[], int n, int (* cmpFunc)(void *, void *) );
extern int findTreeInList ( bestlist *bt, tree *tr );
extern int saveBestTree ( bestlist *bt, tree *tr );
extern int recallBestTreeRecursive ( bestlist *bt, int rank, tree *tr );
extern int recallBestTree ( bestlist *bt, int rank, tree *tr );
extern int recallBestTopology ( bestlist *bt, int rank, tree *tr );
extern boolean readKeyValue ( char *string, char *key, char *format, void *value );
extern int lookupBipartition ( tree *tr, nodeptr p );
extern char *Tree2String ( char *treestr, tree *tr, nodeptr p, boolean printBranchLengths, boolean printNames, boolean printLikelihood, boolean rellTree, boolean finalPrint, analdef *adef, int perGene );
extern void printTreePerGene(tree *tr, analdef *adef, char *fileName, char *permission);
extern int treeFinishCom ( FILE *fp, char **strp );
extern int treeGetCh ( FILE *fp );
extern boolean treeLabelEnd ( int ch );
extern boolean treeGetLabel ( FILE *fp, char *lblPtr, int maxlen );
extern boolean treeFlushLabel ( FILE *fp );
extern int treeFindTipByLabel ( char *str, tree *tr );
extern int treeFindTipName ( FILE *fp, tree *tr );
extern void treeEchoContext ( FILE *fp1, FILE *fp2, int n );
extern boolean treeProcessLength ( FILE *fp, double *dptr );
extern int treeFlushLen ( FILE *fp );
extern boolean treeNeedCh ( FILE *fp, int c1, char *where );
extern boolean addElementLen ( FILE *fp, tree *tr, nodeptr p, boolean readBranchLengths );
extern int saveTreeCom ( char **comstrp );
extern boolean processTreeCom ( FILE *fp, tree *tr );
extern boolean treeReadLen ( FILE *fp, tree *tr, analdef *adef );
extern void treeReadTopologyOnly ( FILE *fp, tree *tr, analdef *adef, boolean readBranches );
extern boolean addElementLenMULT ( FILE *fp, tree *tr, nodeptr p, int partitionCounter );
extern boolean treeReadLenMULT ( FILE *fp, tree *tr, analdef *adef );
extern void getStartingTree ( tree *tr, analdef *adef );
extern double treeLength(tree *tr, int model);
extern void restoreTL_Light(topolRELL_LIST *rl, tree *tr, int n);
extern void computeBootStopOnly(tree *tr, analdef *adef, char *bootStrapFileName);
extern boolean bootStop(tree *tr, BL *b, int numberOfTrees, double *pearsonAverage, double bootstopCutoff);
extern double evaluatePartialGeneric (tree *, int i, double ki);
extern double evaluateGeneric (tree *tr, nodeptr p);
extern void newviewGeneric (tree *tr, nodeptr p);
extern void makenewzGeneric(tree *tr, nodeptr p, nodeptr q, double *z0, int maxiter, double *result);
extern void makenewzGenericDistance(tree *tr, int maxiter, double *z0, double *result, int taxon1, int taxon2);
extern double evaluatePartitionGeneric (tree *tr, nodeptr p, int model);
extern void newviewPartitionGeneric (tree *tr, nodeptr p, int model);
extern void evaluateGenericVector (tree *tr, nodeptr p, double *v);
extern void categorizeGeneric (tree *tr, nodeptr p);
extern double makenewzPartitionGeneric(tree *tr, nodeptr p, nodeptr q, double z0, int maxiter, int model);
extern boolean isTip(int number, int maxTips);
extern double *getLikelihoodArray(int number, int mxtips, double **xVector);
extern int *getScalingArray(int number, int endsite, int mxtips, int *scalingArray);
extern void computeTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches);


extern void   newviewIterative(tree *, int startIndex, int endIndex);

extern double evaluateIterative(tree *, int startIndex, int endIndex);
extern double evaluateIterativePartition(tree *, int lower, int upper, int model);

extern void newviewIterativePartition(tree *, int lower, int upper, int model);

extern void makenewzIterative(tree *, int startIndex, int endIndex);
extern void execCore(tree *, double *dlnLdlz, double *d2lnLdlz2, int lower, int upper, int model);
extern void makenewzIterativePartition(tree *, int startIndex, int endIndex, int model);
extern void execCorePartition(tree *, double *dlnLdlz, double *d2lnLdlz2, int lower, int upper, int model);
extern void determineFullTraversal(nodeptr p, tree *tr);
extern void optRateCat(tree *, int i, double lower_spacing, double upper_spacing, double *lhs);

extern int evaluateParsimonyIterative(tree *, int lower, int upper);
extern void newviewParsimonyIterative(tree *, int startIndex, int endIndex);

extern double evaluateGenericInitrav (tree *tr, nodeptr p);
extern double evaluateGenericInitravPartition(tree *tr, nodeptr p, int model);
extern void evaluateGenericVectorIterative(tree *, int startIndex, int endIndex);
extern void categorizeIterative(tree *, int startIndex, int endIndex);
extern void onlyInitrav(tree *tr, nodeptr p);
extern void onlyInitravPartition(tree *tr, nodeptr p, int model);

extern void fixModelIndices(tree *tr, analdef *adef, int endsite);
extern void calculateModelOffsets(tree *tr);
extern void gammaToCat(tree *tr, analdef *adef);
extern void catToGamma(tree *tr, analdef *adef);
extern void handleExcludeFile(tree *tr, analdef *adef, rawdata *rdta);

extern nodeptr findAnyTip(nodeptr p, int numsp);
extern void determineSequencePosition(tree *tr, analdef *adef);

extern void parseProteinModel(analdef *adef);

extern void computeFullTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int numBranches);

extern void optimizeArndt(tree *tr, analdef *adef);

extern void rapidML_Addition(tree *tr, analdef *adef);

extern void quickOpt(tree *tr, analdef *adef);

extern void computeNextReplicate(tree *tr, analdef *adef, int *originalRateCategories, int *originalInvariant);

extern void reductionCleanup(tree *tr, analdef *adef, int *originalRateCategories, int *originalInvariant);

#ifdef _MULTI_GENE
extern  boolean treeEvaluateMulti(tree *tr, double smoothFactor);
extern  void determineFullMultiTraversal(tree *tr);
extern  void computeMultiTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int model);
extern  void getxsnode (nodeptr p, int model);
extern  void computeFullMultiTraversalInfo(nodeptr p, traversalInfo *ti, int *counter, int maxTips, int model);
#endif


#ifdef _USE_PTHREADS

#define THREAD_NEWVIEW                0
#define THREAD_EVALUATE               1
#define THREAD_MAKENEWZ               2
#define THREAD_SUM_MAKENEWZ           3
#define THREAD_NEWVIEW_PARTITION      4
#define THREAD_EVALUATE_PARTITION     5
#define THREAD_SUM_MAKENEWZ_PARTITION 6
#define THREAD_MAKENEWZ_PARTITION     7
#define THREAD_RATE_CATS              8
#define THREAD_NEWVIEW_PARSIMONY      9
#define THREAD_EVALUATE_PARSIMONY     10
#define THREAD_EVALUATE_VECTOR        11
#define THREAD_CATEGORIZE             12

#ifdef _LOCAL_DATA

#define THREAD_PREPARE_PARSIMONY        13
#define THREAD_FINISH_PARSIMONY         14
#define THREAD_ALLOC_LIKELIHOOD         15
#define THREAD_FREE_LIKELIHOOD          16
#define THREAD_COPY_REVERSIBLE          17
#define THREAD_COPY_RATE_CATS           18
#define THREAD_COPY_GAMMA_RATES         19
#define THREAD_COPY_INVARIANTS          20
#define THREAD_COPY_INIT_MODEL          21
#define THREAD_NEXT_REPLICATE           22




extern void optRateCat_LOCAL(tree *localTree, int lower, int upper, double lower_spacing, double upper_spacing, double *lhs);
#endif

extern void masterBarrier(int jobType, tree *tr);
#endif
