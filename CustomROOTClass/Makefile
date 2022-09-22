# Makefile for ROOT and/or PYTHIA scripts

# DON'T TOUCH THE FOLLOWING LINE
CXXFLAGS:=$(shell root-config --cflags) -I$(PYTHIA8)/include $(shell root-config --libs) -L$(PYTHIA8)/lib -lpythia8

##### The following section allows for ROOT I/O of custom classes:
##### This is needed if you want to fill a Tree with a custom Event class for instance

ROOTPATH:=-I$(ROOTSYS)/include 		# Path to ROOT header files (TFile.h, ...)
PWD=$(shell pwd)					# Present Working Directory
HEADERS=classtest.h					# Headers containing the classes
SOURCES=classtest.cpp				# Source files accompanying the headers
LINKDEF=linkdef.h 					# linkdef file containing pragma statements concerning which classes to include in the dictionary
DICT=EventDict.cpp					# Name of the dictionary to be created, can pretty much be anything

# Only do this if test.cpp exists. Executing `make test` targets this specifically. Explanation of the recipe below.
# remove all the '@' at the start of the lines for debugging
test: test.cpp
	@echo -n Making dictionary with rootcling...\  

	@rootcling -f $(DICT) $(ROOTPATH) $(HEADERS) $(LINKDEF)

	@echo Success!
	@echo -n Making library libEvent.so...\ 

	@g++ -shared -o libEvent.so -fPIC -std=c++17 $(ROOTPATH) $(DICT) $(SOURCES) $(root-config --ldflags --libs)
	
	@echo Success!
	@echo -n Creating executable and linking libraries...\ 

	@g++ -o $@ $^ -L. -lEvent -Wl,-rpath=$(PWD) $(CXXFLAGS)

	@echo Success!
	@echo Looks like everything worked! Cleaning up...
	
	@rm $(DICT)

# The first line creates a dictionary EventDict.cpp from the headers supplied. 
#	It only creates dictionaries for the classes specified in linkdef.h 
# The second line makes a library libEvent.so from the source files of the headers and the dictionary created in the first step.
# 	I don't know exactly what the -fPIC and root-config flags do, but I'm pretty sure they are needed
# 	The -std=c++17 flag forces c++17, which is required by ROOT.
#	-I... flag specifies the location of other header files needed (I.E. ROOT headers such as TFile.h, ...)
# The third makes the executable and links the required libraries. 
#	-L. says that it should look in the current directory for libraries as well (hence the dot) when compiling
#	-lEvent corresponds to libEvent.so (-lXXX --> libXXX.so)
#	-Wl,rpath=$(PWD) is similar to -L., except it tells the executable to look in the current directory for libraries.
#		This is important, otherwise you need to manually set the path after every compile
# 	$(CXXFLAGS) are kinda obvious I think.