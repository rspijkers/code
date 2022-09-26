#include "SmallTrack.h"

// Default constructor
SmallTrack::SmallTrack(){
    pdg = 0;
    pT = 0;
    eta = 0;
    phi = 0;
}

// clear all info
void SmallTrack::Clear(){
	pdg = 0;
    pT = 0;
    eta = 0;
    phi = 0;;
}

// Destructor
SmallTrack::~SmallTrack(){
	Clear();
}
