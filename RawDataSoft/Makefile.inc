#####################################################################
#
# Set global compile / link flags here
#
#####################################################################

VERSION = 1
SUBVER  = 0

TOP := $(dir $(lastword $(MAKEFILE_LIST)))
LIBDIR = $(TOP)lib
INCDIR = $(TOP)include

DOXYGEN    = doxygen

ARCH       = $(shell root-config --arch)
ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)
ROOTGLIBS  = $(shell root-config --glibs)
ROOTEVLIBS = $(shell root-config --evelibs)

CXX        = g++
CXXFLAGS   = -g -Wall -fPIC
LD         = g++

CXXFLAGS  += $(ROOTCFLAGS) -I $(INCDIR) # -D WITHLIBEVENT2
LIBS       = $(ROOTLIBS) $(SYSLIBS) -lMinuit 
GLIBS      = $(ROOTGLIBS) $(SYSLIBS)
LIBEV2     = #-levent -levent_pthreads # to enable with WITHLIBEVENT
EVLIBS     = $(ROOTEVLIBS) -lVMC -lEG -lGeomPainter -lFFTW -pthread $(LIBEV2)
LDFLAGS    = -Wl,--no-as-needed $(LIBS) $(GLIBS) $(EVLIBS)
SOFLAGS    = -shared

# libevent2 and threads


AR         = ar
ARFLAGS    = rcsv
