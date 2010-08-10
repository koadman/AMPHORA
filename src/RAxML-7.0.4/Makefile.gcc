# Makefile August 2006 by Alexandros Stamatakis
# Makefile cleanup October 2006, Courtesy of Peter Cordes <peter@cordes.ca>

CC = gcc 
CFLAGS =  -O3 -fomit-frame-pointer -funroll-loops
#-Wall -pedantic -Wunused-parameter -Wredundant-decls  -Wreturn-type  -Wswitch-default -Wunused-value -Wimplicit  -Wimplicit-function-declaration  -Wimplicit-int -Wimport -pedantic-errors -Wunused  -Wunused-function  -Wunused-label -Wno-int-to-pointer-cast -Wbad-function-cast  -Wmissing-declarations -Wmissing-prototypes  -Wnested-externs  -Wold-style-definition -Wstrict-prototypes  -Wdeclaration-after-statement -Wpointer-sign -Wextra -Wmissing-noreturn -Wredundant-decls -Wunused -Wunused-function -Wunused-parameter -Wunused-value  -Wunused-variable

LIBRARIES = -lm 

RM = rm -f

objs    = axml.o raxmlParsimony.o rev_functions.o optimizeModel.o multiple.o searchAlgo.o topologies.o parsePartitions.o treeIO.o models.o bipartitionList.o rapidBootstrap.o evaluatePartialGeneric.o evaluateGeneric.o newviewGeneric.o makenewzGeneric.o  evaluateGenericVector.o categorizeGeneric.o 

all : raxmlHPC 

GLOBAL_DEPS = axml.h globalVariables.h

raxmlHPC : $(objs)
	$(CC) -o raxmlHPC $(objs) $(LIBRARIES) 

# Optimization flag -O3 will yield an erroneous code for rev_functions !

rev_functions.o : rev_functions.c $(GLOBAL_DEPS)
	$(CC) -c -o $@ $<

evaluatePartialGeneric.o : evaluatePartialGeneric.c $(GLOBAL_DEPS)
bipartitionList.o : bipartitionList.c $(GLOBAL_DEPS)
optimizeModel.o : optimizeModel.c $(GLOBAL_DEPS)
multiple.o : multiple.c $(GLOBAL_DEPS)
axml.o : axml.c $(GLOBAL_DEPS)
raxmlParsimony.o : raxmlParsimony.c $(GLOBAL_DEPS)
searchAlgo.o : searchAlgo.c $(GLOBAL_DEPS)
topologies.o : topologies.c $(GLOBAL_DEPS)
parsePartitions.o : parsePartitions.c $(GLOBAL_DEPS)
treeIO.o : treeIO.c $(GLOBAL_DEPS)
models.o : models.c $(GLOBAL_DEPS)
rapidBootstrap.o : rapidBootstrap.c $(GLOBAL_DEPS) 
#evaluatePartialGeneric.o : evaluatePartialGeneric.c $(GLOBAL_DEPS)
evaluateGeneric.o : evaluateGeneric.c $(GLOBAL_DEPS)
newviewGeneric.o : newviewGeneric.c $(GLOBAL_DEPS)
makenewzGeneric.o : makenewzGeneric.c $(GLOBAL_DEPS)
evaluateGenericVector.o : evaluateGenericVector.c $(GLOBAL_DEPS)
categorizeGeneric.o : categorizeGeneric.c $(GLOBAL_DEPS)

clean : 
	$(RM) *.o raxmlHPC
