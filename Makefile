# Makefile for ROOT and/or PYTHIA scripts

# DON'T TOUCH THE FOLLOWING LINE
CXXFLAGS:=$(shell root-config --cflags) $(shell root-config --libs)

##### The following section allows for ROOT I/O of custom classes:
##### This is needed if you want to fill a Tree with a custom Event class for instance

CPATH:=/home/rik/cernbox/PhD/code# 	# path to code repository

ifeq (stbc,$(findstring stbc,$(shell hostname)))
CPATH:=/user/rspijker/project/code
$(info stbc node detected! Changing path accordingly.)
endif

LIB:=$(CPATH)/lib#
SRC:=$(CPATH)/src#
INCLUDE:=$(CPATH)/include#
ROOTPATH:=$(ROOTSYS)/include# 		# Path to ROOT header files (TFile.h, ...)
HEADERS := SmallTrack.h SmallEvent.h#				# Headers containing the classes
FILES := SmallTrack.cpp SmallEvent.cpp#				# Source files accompanying the headers
SOURCES := $(FILES:%.cpp=$(SRC)/%.cpp)			# prepend source dir to filenames.
LINKDEF := linkdef.h# 					# linkdef file containing pragma statements concerning which classes to include in the dictionary
DICT:=EventDict.cpp#					# Name of the dictionary to be created, can pretty much be anything

# Only do this if tree_handler_v2.cpp exists. Executing `make tree_handler_v2` targets this specifically. Explanation of the recipe below.
# remove all the '@' at the start of the lines for debugging
tree_handler_v2: tree_handler_v2.cpp
	@echo -n Making dictionary with rootcling...\  

	@rootcling -f $(LIB)/$(DICT) -I$(ROOTPATH) -I$(INCLUDE) $(HEADERS) $(LINKDEF)

	@echo Success!
	@echo -n Making library libEvent.so...\ 

	@g++ -shared -o $(LIB)/libEvent.so -fPIC -std=c++17 -I$(ROOTPATH) -I$(INCLUDE) -I$(SRC) $(LIB)/$(DICT) $(SOURCES) $(root-config --ldflags --libs)
	
	@echo Success!
	@echo -n Creating executable and linking libraries...\ 

	@g++ -o $@ $^ -L$(LIB) -lEvent -Wl,-rpath=$(LIB) $(CXXFLAGS) -I$(INCLUDE)

	@echo Success!
	@echo Looks like everything worked! Cleaning up...
	
	@rm $(LIB)/$(DICT)

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