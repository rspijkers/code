#include "SmallEvent.h"

// Default constructor
SmallEvent::SmallEvent(){
    eventId = 0;
    Ntracks = 0;
	pTssbar = 0;
    Candidates = {};
}

// clear all info
void SmallEvent::Clear(){
	eventId = 0;
	Ntracks = 0;
	pTssbar = 0;
	Candidates.clear();
}

// Destructor
SmallEvent::~SmallEvent(){
	Clear();
}