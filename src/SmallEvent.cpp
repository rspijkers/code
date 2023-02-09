#include "SmallEvent.h"

// Default constructor
SmallEvent::SmallEvent(){
    eventId = 0;
    Ntracks4p0 = 0;
	Ntracks2p4 = 0;
	Ntracks0p8 = 0;
	pTssbar = 0;
    Candidates = {};
}

// clear all info
void SmallEvent::Clear(){
	eventId = 0;
	Ntracks4p0 = 0;
	Ntracks2p4 = 0;
	Ntracks0p8 = 0;
	pTssbar = 0;
	Candidates.clear();
}

// Destructor
SmallEvent::~SmallEvent(){
	Clear();
}