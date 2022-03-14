CXXFLAGS:=$(shell root-config --cflags) -I$(PYTHIA8)/include  \
    $(shell root-config --libs) -L$(PYTHIA8)/lib -lpythia8

# .SUFFIXES:      .o .cxx

# pythia8_simple: pythia8_simple.o

