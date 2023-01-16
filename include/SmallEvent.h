#ifndef SMALLEVENT_H
#define SMALLEVENT_H
#include "SmallTrack.h"
#include <vector>

class SmallEvent {
    private:
        Int_t eventId; 
        Int_t Ntracks;
        Double_t pTssbar;
        std::vector<SmallTrack> Candidates;

    public:
        SmallEvent();
        ~SmallEvent();
        void Clear();

        // setters
        void setEventId(Int_t i) {eventId = i;}
        void setNtracks(Int_t n) {Ntracks = n;}
        void setpTssbar(Double_t pT) {pTssbar = pT;}
        void addCandidate(SmallTrack cand) {Candidates.push_back(cand);}
        void setCandidates(std::vector<SmallTrack> cands) {Candidates = cands;}

        // getters
        Int_t getEventId() const {return eventId;}
        Int_t getNtracks() const {return Ntracks;}
        Double_t getpTssbar() const {return pTssbar;}
        std::vector<SmallTrack> getCandidates() const {return Candidates;}

        ClassDef(SmallEvent,1); 
};
#endif