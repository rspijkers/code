#include "SmallEvent.h"

// Default constructor
SmallEvent::SmallEvent(){
    eventId = 0;
    Ntracks = 0;
    Tracks = {};
    Candidates = {};
}

// clear all info
void SmallEvent::Clear(){
	eventId = 0;
	Ntracks = 0;
	Tracks.clear();
	Candidates.clear();
}

// Destructor
SmallEvent::~SmallEvent(){
	Clear();
}