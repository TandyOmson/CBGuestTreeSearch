
#  nab configuration file.
#  Created on Wed Jul 17 16:35:09 +08 2024 via ./configure --openblas

CC=gcc -Wall
CXX=g++
FC=gfortran
FLIBS_ARCH=-lgfortran -w
LDFLAGS=

###############################################################################

# (1)  Location of the installation

BASEDIR=/home/spine/nabc
BINDIR=$(BASEDIR)/bin
LIBDIR=$(BASEDIR)/lib
INCDIR=$(BASEDIR)/include

###############################################################################

#  (2) Flags that depend on OS type

SHARED_SUFFIX=.so
MAKE_SHARED=-shared
LM=-lm

###############################################################################

#  (3)  CC compilers

CFLAGS=-DBINTRAJ -I/home/spine/nabc/include    
CNOOPTFLAGS=-g -O0
COPTFLAGS=-g -Ofast -mtune=native

###############################################################################

#  (4)  other flags:

FLIBS= -lsff -lnabc -larpack -llapack -lblas -lnetcdf  -lgfortran -w
CIFTARGET=skip
AR=    ar rv
RANLIB=ranlib
MV=mv
CP=cp
MAKE=make --no-print-directory
VB=@

