#ifndef SMALLEVENT_H
#define SMALLEVENT_H
#include "SmallTrack.h"
#include <vector>

class SmallEvent {
    private:
        Int_t eventId; 
        Int_t Ntracks;
        std::vector<SmallTrack> Tracks;
        std::vector<SmallTrack> Candidates;

    public:
        SmallEvent();
        ~SmallEvent();
        void Clear();

        // setters
        void setEventId(Int_t i) {eventId = i;}
        void setNtracks(Int_t n) {Ntracks = n;}
        void addTrack(SmallTrack track) {Tracks.push_back(track);}
        void addCandidate(SmallTrack cand) {Candidates.push_back(cand);}
        void setCandidates(std::vector<SmallTrack> cands) {Candidates = cands;}

        // getters
        Int_t getEventId() const {return eventId;}
        Int_t getNtracks() const {return Ntracks;}
        std::vector<SmallTrack> getCandidates() const {return Candidates;}

        ClassDef(SmallEvent,1); 
};
#endif