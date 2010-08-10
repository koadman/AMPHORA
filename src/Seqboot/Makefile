# Makefile
#
# Generic Linux/Unix Makefile for PHYLIP 3.6.

PACKAGE=phylip
VERSION=3.67

# You should not need to change anything, though if you want you could
# change the first (noncomment) statements to some directory locations
# that might be more useful.
# The compressed tar archive phylip.tar.gz when uncompressed and extracted
# puts the source code into a directory ./src, and also makes two other
# directories ./exe and ./doc for the final executables and the documentation
# files.  Only change the EXEDIR settings if you want something different
# from that structure.
# If it causes trouble in compiling, the CFLAGS statement below may also need
# to be changed.
#
#  To use the PHYLIP v3.6 Makefile, type
#      make install       to compile the whole package, install
#                           the executables in $(EXEDIR), and
#                           remove the object files to save space
#      make all           to compile the whole package but not install it
#                           or remove the object files
#      make put           to move the executables into $(EXEDIR)
#      make clean         to remove all object files and executables from the
#                           current directory
#      make dnaml         to compile and link one program, (in this example,
#                           DnaML) and leave the executable and object files
#                           in the current directory (where the source code is).
#                           You will have to move the executable into the
#                           executables directory (e.g. "mv dnaml ../exe")
#                           Note that the program name should be lower case.
# 
# ----------------------------------------------------------------------------
#  (Starting here is the section where you may want to change things)
# ----------------------------------------------------------------------------
#
# these are the statements we have been talking about:
# one of the reasons for changing them would be to put the executables
# and/or the documentation on a different file system.
# The default configuration is to have within the overall PHYLIP
# directory three subdirectories:  "src" for source code, "exe" for the
# executables, and "doc" for the documentation files.
#
# the following specifies the directory where the executables will be placed
EXEDIR  = ../exe
#
# ----------------------------------------------------------------------------
#
# In the following statements (the ones that set CFLAGS, DFLAGS, LIBS
# and DLIBS, CC and DC) you should make sure each is set properly.
# Usually this will simply involve making sure that the proper statement
# has no "#" as its first character and that all other possibilities
# have "#" for their first character.
#
# ----------------------------------------------------------------------------
#
# This is the CFLAGS statement:
#
# if these statements say "-g"  and that causes trouble, remove the "-g"
# if the "-lX11" does not work (X compiling not working), you may want to
# remove that switch or fix it.
#
# Here are some possible CFLAGS statements:
#
#
#A minimal one
#CFLAGS =
#
# A basic one for debugging
CFLAGS  = -g 
#
# An optimized one for gcc
#CFLAGS  = -O3 -fomit-frame-pointer
#
# For some serious debugging using Gnu gcc
#
#CFLAGS=-g -Wall -Wmain -Wmissing-prototypes -Wreturn-type -Wstrict-prototypes  -Wunused -Werror -Wredundant-decls -Waggregate-return -Wcast-align -Wcomment
#
# For Digital Alpha systems with Compaq Tru64 Unix
# (however, be aware that this may cause floating-point problems in programs
#  like Dnaml owing to not using IEEE floating point standards).
#CFLAGS = -fast
#
# ----------------------------------------------------------------------------
#
# and here are some possible DFLAGS statements:
#
# A minimal one
#DFLAGS  = -DX $(CFLAGS)
#
# A basic one for debugging
#DFLAGS  = -g -DX
#
#
# For Gnu C++ for runs
#DFLAGS  = -I/usr/X11R6/include -O3 -DX -fomit-frame-pointer
#
# For Digital Alpha systems with Compaq Tru64 Unix
#DFLAGS = -DX -fast 
#
# for Linux with X Windows development packages installed
# or for MacOS X with X Windows installed
DFLAGS = $(CFLAGS) -DX -I/usr/X11R6/include
#
# ----------------------------------------------------------------------------
#
# These are the libraries for the CC and DC compiles, respectively
#
LIBS    = -lm
#
# if the Xlib library for the X windowing system is somewhere
# unexpected, you may have to change the path /usr/X11R6/lib in this one
#
# For gcc for Linux with X windows development packages installed
# or for MacOS X with X windows installed
DLIBS= -L/usr/X11R6/lib/ -lX11 -lXaw -lXt
#
# ----------------------------------------------------------------------------
#
# The next two assignments are the invocations of the compiler for the
# ordinary compiles and the tree-drawing programs, CC and DC
#
# This one specifies the "cc" C compiler
#CC        = cc
#
#  To use GCC instead, if it is not the compiler that "cc" invokes
CC        = gcc
#
# This one specifies the "cc" C compiler for the Draw programs
DC        = cc $(DFLAGS)
#
# To use GCC instead, if it is not the compiler that "cc" invokes
#DC        = gcc $(DFLAGS)
#
# ----------------------------------------------------------------------------
#  (After this point there should not be any reason to change anything)
# ----------------------------------------------------------------------------
#
#
# the list of programs
#
PROGS   =          clique consense contml contrast dnacomp dnadist \
                   dnainvar dnaml dnamlk dnamove dnapars dnapenny \
                   dolmove dollop dolpenny factor fitch gendist kitsch \
                   mix move neighbor pars penny proml promlk protdist \
                   protpars restdist restml retree seqboot treedist \
                   drawgram drawtree

#
# general commands
#

#  The first target it executed if you just type "make".  It tells you how to
#  use the Makefile.
#
help:
	@echo ""
	@echo " To use the PHYLIP v3.6 Makefile, type"
	@echo "     make install       to compile the whole package and install"
	@echo "                          the executables in $(EXEDIR), and then"
	@echo "                          remove the object files to save space"
	@echo "     make all           to compile the whole package but not install it"
	@echo "                          or remove the object files"
	@echo "     make put           to move the executables into $(EXEDIR)"
	@echo "     make clean         to remove all object files and executables from the"
	@echo "                          current directory"
	@echo "     make dnaml         to compile and link one program, (in this example,"
	@echo "                          Dnaml) and leave the executable and object files"
	@echo "                          in the current directory (where the source code is)."
	@echo "                          You will have to move the executable into the"
	@echo "                          executables directory (e.g. \"mv dnaml $(EXEDIR)\")"
	@echo "                          Note that the program name should be lower case."
	@echo " "

introduce:
	@echo "Building PHYLIP version $(VERSION)"

all:        introduce $(PROGS)
	@echo "Finished compiling."
	@echo ""

install:        all put clean
	@echo "Done."
	@echo ""

put:
	@echo "Installing PHYLIP v3.6 binaries in $(EXEDIR)"
	@mkdir -p $(EXEDIR)
	@cp $(PROGS) $(EXEDIR)
	@echo "Installing font files in $(EXEDIR)"
	@cp font* $(EXEDIR)
	@echo "Finished installation."
	@echo ""

clean:
	@echo "Removing object files to save space"
	@rm -f *.o
	@echo "Finished removing object files.  Now will remove"
	@echo "executable files from the current directory, but not from the"
	@echo "executables directory.  (If some are not here, the makefile"
	@echo "will terminate with an error message but this is not a problem)"
	@echo ""
	@echo "Removing executables from this directory"
	@rm -f $(PROGS)
	@echo "Finished cleanup."
	@echo ""

#
# compile object files shared between programs
# (make's implicit rule for %.o will take care of these)
#

phylip.o:     phylip.h
seq.o:        phylip.h seq.h
disc.o:       phylip.h disc.h
discrete.o:   phylip.h discrete.h
dollo.o:      phylip.h dollo.h
wagner.o:     phylip.h wagner.h
dist.o:       phylip.h dist.h
cont.o:       phylip.h cont.h
mlclock.o:    phylip.h mlclock.h
moves.o:      phylip.h moves.h
printree.o:   phylip.h printree.h

#
# compile the individual programs
#

clique.o:       disc.h phylip.h

clique:         clique.o disc.o phylip.o 
	$(CC) $(CFLAGS) $^ $(LIBS) -o clique

cons.o:         cons.h phylip.h 

consense.o:     cons.h phylip.h

consense:       consense.o phylip.o cons.o
	$(CC) $(CFLAGS) $^ $(LIBS) -o consense

contml.o:       cont.h phylip.h

contml:       contml.o cont.o phylip.o 
	$(CC) $(CFLAGS) $^ $(LIBS) -o contml

contrast.o:       contrast.c cont.h phylip.h
	$(CC) $(CFLAGS) -c contrast.c

contrast:       contrast.o cont.o phylip.o 
	$(CC) $(CFLAGS) contrast.o cont.o phylip.o $(LIBS) -o contrast

dnacomp.o:      dnacomp.c seq.h phylip.h
	$(CC) $(CFLAGS) -c dnacomp.c

dnacomp:      dnacomp.o seq.o phylip.o 
	$(CC) $(CFLAGS) dnacomp.o seq.o phylip.o $(LIBS) -o dnacomp

dnadist.o:      dnadist.c seq.h phylip.h
	$(CC) $(CFLAGS) -c dnadist.c

dnadist:      dnadist.o seq.o phylip.o 
	$(CC) $(CFLAGS) dnadist.o seq.o phylip.o $(LIBS) -o dnadist

dnainvar.o:      dnainvar.c seq.h phylip.h
	$(CC) $(CFLAGS) -c dnainvar.c

dnainvar:      dnainvar.o seq.o phylip.o 
	$(CC) $(CFLAGS) dnainvar.o seq.o phylip.o $(LIBS) -o dnainvar

dnaml.o:      dnaml.c seq.h phylip.h
	$(CC) $(CFLAGS) -c dnaml.c

dnaml:      dnaml.o seq.o phylip.o 
	$(CC) $(CFLAGS) dnaml.o seq.o phylip.o $(LIBS) -o dnaml

dnamlk.o:      dnamlk.c seq.h phylip.h mlclock.h printree.h
	$(CC) $(CFLAGS) -c dnamlk.c

dnamlk:      dnamlk.o seq.o phylip.o mlclock.o printree.o
	$(CC) $(CFLAGS) $^ $(LIBS) -o dnamlk

dnamove.o:      dnamove.c seq.h moves.h phylip.h
	$(CC) $(CFLAGS) -c dnamove.c

dnamove:      dnamove.o seq.o moves.o phylip.o 
	$(CC) $(CFLAGS) dnamove.o seq.o moves.o phylip.o $(LIBS) -o dnamove

dnapenny.o:      dnapenny.c seq.h phylip.h
	$(CC) $(CFLAGS) -c dnapenny.c

dnapenny:      dnapenny.o seq.o phylip.o 
	$(CC) $(CFLAGS) dnapenny.o seq.o phylip.o $(LIBS) -o dnapenny

dnapars.o:      dnapars.c seq.h phylip.h
	$(CC) $(CFLAGS) -c dnapars.c

dnapars:      dnapars.o seq.o phylip.o 
	$(CC) $(CFLAGS) dnapars.o seq.o phylip.o $(LIBS) -o dnapars

dolmove.o:       dolmove.c disc.h moves.h dollo.h phylip.h
	$(CC) $(CFLAGS) -c dolmove.c

dolmove:       dolmove.o disc.o moves.o dollo.o phylip.o 
	$(CC) $(CFLAGS) dolmove.o disc.o moves.o dollo.o phylip.o $(LIBS) -o dolmove

dollop.o:       dollop.c disc.h dollo.h phylip.h
	$(CC) $(CFLAGS) -c dollop.c

dollop:       dollop.o disc.o dollo.o phylip.o 
	$(CC) $(CFLAGS) dollop.o disc.o dollo.o phylip.o $(LIBS) -o dollop

dolpenny.o:       dolpenny.c disc.h dollo.h phylip.h
	$(CC) $(CFLAGS) -c dolpenny.c

dolpenny:       dolpenny.o disc.o dollo.o phylip.o 
	$(CC) $(CFLAGS) dolpenny.o disc.o dollo.o phylip.o $(LIBS) -o dolpenny

factor.o:       factor.c phylip.h
	$(CC) $(CFLAGS) -c factor.c

factor:       factor.o phylip.o 
	$(CC) $(CFLAGS) factor.o phylip.o $(LIBS) -o factor

fitch.o:        fitch.c dist.h phylip.h
	$(CC) $(CFLAGS) -c fitch.c

fitch:        fitch.o dist.o phylip.o 
	$(CC) $(CFLAGS) fitch.o dist.o phylip.o $(LIBS) -o fitch

gendist.o:      gendist.c phylip.h
	$(CC) $(CFLAGS) -c gendist.c

gendist:      gendist.o phylip.o 
	$(CC) $(CFLAGS) gendist.o phylip.o $(LIBS) -o gendist

kitsch.o:        kitsch.c dist.h phylip.h
	$(CC) $(CFLAGS) -c kitsch.c

kitsch:        kitsch.o dist.o phylip.o 
	$(CC) $(CFLAGS) kitsch.o dist.o phylip.o $(LIBS) -o kitsch

mix.o:        mix.c disc.h wagner.h phylip.h
	$(CC) $(CFLAGS) -c mix.c

mix:        mix.o disc.o wagner.o phylip.o 
	$(CC) $(CFLAGS) mix.o disc.o wagner.o phylip.o $(LIBS) -o mix

move.o:        move.c disc.h moves.h wagner.h phylip.h
	$(CC) $(CFLAGS) -c move.c

move:        move.o disc.o moves.o wagner.o phylip.o 
	$(CC) $(CFLAGS) move.o disc.o moves.o wagner.o phylip.o $(LIBS) -o move

neighbor.o:        neighbor.c dist.h phylip.h
	$(CC) $(CFLAGS) -c neighbor.c

neighbor:        neighbor.o dist.o phylip.o 
	$(CC) $(CFLAGS) neighbor.o dist.o phylip.o $(LIBS) -o neighbor

pars.o:   pars.c discrete.h phylip.h
	$(CC) $(CFLAGS) -c pars.c

pars: pars.o pars.c discrete.o discrete.c phylip.o 
	$(CC) $(CFLAGS) pars.o discrete.o phylip.o $(LIBS) -o pars

penny.o:  penny.c disc.h wagner.h phylip.h
	$(CC) $(CFLAGS) -c penny.c

penny:  penny.o disc.o wagner.o phylip.o 
	$(CC) $(CFLAGS) penny.o disc.o wagner.o  phylip.o $(LIBS) -o penny

proml.o:      proml.c seq.h phylip.h
	$(CC) $(CFLAGS) -c proml.c

proml:      proml.o seq.o phylip.o
	$(CC) $(CFLAGS) proml.o seq.o phylip.o $(LIBS) -o proml

promlk.o:      promlk.c seq.h phylip.h mlclock.h printree.h
	$(CC) $(CFLAGS) -c promlk.c

promlk:      promlk.o seq.o phylip.o mlclock.o printree.o
	$(CC) $(CFLAGS) $^ $(LIBS) -o promlk

protdist.o:      protdist.c seq.h phylip.h
	$(CC) $(CFLAGS) -c protdist.c

protdist:      protdist.o seq.o phylip.o 
	$(CC) $(CFLAGS) protdist.o seq.o phylip.o $(LIBS) -o protdist

protpars.o: protpars.c seq.h phylip.h
	$(CC) $(CFLAGS) -c protpars.c

protpars: protpars.o seq.o phylip.o 
	$(CC) $(CFLAGS) protpars.o seq.o phylip.o $(LIBS) -o protpars

restdist.o: restdist.c seq.h phylip.h
	 $(CC) $(CFLAGS) -c restdist.c

restdist: restdist.o seq.o phylip.o 
	$(CC) $(CFLAGS) restdist.o seq.o phylip.o $(LIBS) -o restdist

restml.o: restml.c seq.h phylip.h
	 $(CC) $(CFLAGS) -c restml.c

restml: restml.o seq.o phylip.o 
	$(CC) $(CFLAGS) restml.o seq.o phylip.o $(LIBS) -o restml

retree.o:       retree.c moves.h phylip.h
	$(CC) $(CFLAGS) -c retree.c

retree:       retree.o moves.o phylip.o 
	$(CC) $(CFLAGS) retree.o moves.o phylip.o $(LIBS) -o retree

seqboot.o:      seqboot.c phylip.h
	$(CC) $(CFLAGS) -c seqboot.c

seqboot:      seqboot.o seq.o phylip.o 
	$(CC) $(CFLAGS) seqboot.o seq.o phylip.o $(LIBS) -o seqboot

treedist.o:     treedist.c cons.h phylip.h
	$(CC) $(CFLAGS) -c treedist.c

treedist:     treedist.o phylip.o cons.o
	$(CC) $(CFLAGS) treedist.o cons.o phylip.o $(LIBS) -o treedist


# ----------------------------------------------------------------------------
#  The following section is used to build a PHYLIP distribution. All sources
#  and documentation files must be placed in the current directory.
#
#  Usage:
#  	make distdir	- Build the distribution dir phylip-<vers>/
#  	make dist	- Make a tarred and gzipped phylip-<vers>.tar.gz
# ----------------------------------------------------------------------------

DIST_COMMON = phylip.html 

DOCS= doc/clique.html doc/consense.html doc/contchar.html doc/contml.html \
      doc/contrast.html doc/discrete.html doc/distance.html doc/dnacomp.html \
      doc/dnadist.html doc/dnainvar.html doc/dnaml.html doc/dnamlk.html \
      doc/dnamove.html doc/dnapars.html doc/dnapenny.html doc/dollop.html \
      doc/dolmove.html doc/dolpenny.html doc/drawgram.html doc/draw.html \
      doc/drawtree.html doc/factor.html doc/fitch.html doc/gendist.html \
      doc/kitsch.html doc/main.html doc/mix.html doc/move.html \
      doc/neighbor.html doc/pars.html doc/penny.html doc/proml.html \
      doc/promlk.html doc/protdist.html doc/protpars.html doc/restdist.html \
      doc/restml.html doc/retree.html doc/seqboot.html doc/sequence.html \
      doc/treedist.html doc/phylip.gif

SOURCES= Makefile Makefile.bcc Makefile.cyg Makefile.msvc \
	 Makefile.osx compalpha.com linkmac clique.c cons.c consense.c cons.h cont.c \
	 cont.h contml.c contrast.c disc.c disc.h discrete.c discrete.h dist.c \
	 dist.h dnacomp.c dnadist.c dnainvar.c dnaml.c dnamlk.c dnamove.c \
	 dnapars.c dnapenny.c dollo.c dollo.h dollop.c dolmove.c dolpenny.c \
	 interface.c interface.h io.h kitsch.c macface.c macface.h mix.c move.c \
	 moves.c moves.h neighbor.c newmove.c pars.c penny.c \
	 phylip.c phylip.h proml.c promlk.c protdist.c protpars.c restdist.c \
	 restml.c retree.c seqboot.c seq.c seq.h treedist.c wagner.c wagner.h \
	 mlclock.c mlclock.h printree.c printree.h

MAC= \
		Info.plist.in boot.icns clique.icns command.in consense.icns \
	contml.icns contrast.icns disc.icns dist.icns dna.icns dnacomp.icns \
	dnadist.icns dnainvar.icns dnaml.icns dnamlk.icns dnamove.icns \
	dnapars.icns dnapenny.icns dollo.icns dollop.icns dolmove.icns \
	dolpenny.icns drawgram.icns drawtree.icns factor.icns fitch.icns \
	gendist.icns kitsch.icns mac.sit mix.icns move.icns neighbor.icns \
	pars.icns penny.icns proml.icns promlk.icns protdist.icns protein.icns \
	protpars.icns restdist.icns restml.icns restrict.icns retree.icns \
	seqboot.icns treedist.icns

ICONS= 		boot.ico clique.ico clique.rc clique.rcb consense.ico \
		consense.rc consense.rcb contml.ico contml.rc contml.rcb \
		contrast.ico contrast.rc contrast.rcb disc.ico dist.ico dna.ico \
		dnacomp.rc dnacomp.rcb dnadist.rc dnadist.rcb dnainvar.rc \
		dnainvar.rcb dnaml.rc dnaml.rcb dnamlk.rc dnamlk.rcb dnamove.rc \
		dnamove.rcb dnapars.rc dnapars.rcb dnapenny.rc dnapenny.rcb \
		dollo.ico dollop.rc dollop.rcb dolmove.rc dolmove.rcb \
		dolpenny.rc dolpenny.rcb drawgram.ico drawgram.rc drawgram.rcb \
		drawtree.ico drawtree.rc drawtree.rcb factor.rc factor.rcb \
		fitch.rc fitch.rcb gendist.ico gendist.rc gendist.rcb kitsch.rc \
		kitsch.rcb mix.rc mix.rcb move.rc move.rcb neighbor.rc \
		neighbor.rcb pars.rc pars.rcb penny.rc penny.rcb proml.rc \
		proml.rcb promlk.rc promlk.rcb protdist.rc protdist.rcb \
		protein.ico protpars.rc protpars.rcb restdist.rc restdist.rcb \
		restml.rc restml.rcb restrict.ico retree.ico retree.rc \
		retree.rcb seqboot.rc seqboot.rcb treedist.ico treedist.rc \
		treedist.rcb

FONTS= font1 font2 font3 font4 font5 font6

DISTDIR=$(PACKAGE)-$(VERSION)
dist_SRCDIR=$(DISTDIR)/src
dist_DOCDIR=$(DISTDIR)/doc
dist_EXEDIR=$(DISTDIR)/exe
MACICONDIR=src/mac

SHELL=bash

# We use this target to create a tarred and gzipped distribution of PHYLIP
dist: distdir
	-chmod -R a+r $(DISTDIR)
	tar chozf $(DISTDIR).tar.gz $(DISTDIR)
	-rm -rf $(DISTDIR)

# This target creates the distribution directory
distdir: $(DIST_COMMON) $(DOCS) $(SOURCES)
	-rm -rf $(DISTDIR)
	mkdir $(DISTDIR) && \
	mkdir $(dist_EXEDIR) && \
	mkdir $(dist_DOCDIR) && \
	mkdir $(dist_SRCDIR)
	mkdir $(dist_SRCDIR)/mac
	mkdir $(dist_SRCDIR)/icons
	for i in $(DIST_COMMON) ; do \
	  cp -r $$i $(DISTDIR) ; \
	done
	for i in $(DOCS) ; do \
	  cp -r $$i $(dist_DOCDIR) ; \
	done
	for i in $(SOURCES) ; do \
	  cp -r $$i $(dist_SRCDIR) ; \
	done
	for i in $(MAC) ; do \
	  cp -r mac/$$i $(dist_SRCDIR)/mac ; \
	done
	for i in $(ICONS) ; do \
	  cp -r icons/$$i $(dist_SRCDIR)/icons ; \
	done
	for i in $(FONTS) ; do \
	  cp -r $$i $(dist_SRCDIR) ; \
	done

# This target untars the dist and checks that it can be compiled and remade
distcheck: dist
	-rm -rf $(DISTDIR)
	tar xzf $(DISTDIR).tar.gz
	cd $(DISTDIR)/$(SRCDIR) \
	  && make all
	-rm -rf $(DISTDIR)
	@echo "$(DISTDIR).tar.gz is ready for distribution"

# Makefile
