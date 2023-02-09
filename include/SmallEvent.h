#ifndef SMALLEVENT_H
#define SMALLEVENT_H
#include "SmallTrack.h"
#include <vector>

class SmallEvent {
    private:
        Int_t eventId; 
        Int_t Ntracks4p0;
        Int_t Ntracks2p4;
        Int_t Ntracks0p8;
        Double_t pTssbar;
        std::vector<SmallTrack> Candidates;

    public:
        SmallEvent();
        ~SmallEvent();
        void Clear();

        // setters
        void setEventId(Int_t i) {eventId = i;}
        void setNtracks4p0(Int_t n) {Ntracks4p0 = n;}
        void setNtracks2p4(Int_t n) {Ntracks2p4 = n;}
        void setNtracks0p8(Int_t n) {Ntracks0p8 = n;}
        void setpTssbar(Double_t pT) {pTssbar = pT;}
        void addCandidate(SmallTrack cand) {Candidates.push_back(cand);}
        void setCandidates(std::vector<SmallTrack> cands) {Candidates = cands;}

        // getters
        Int_t getEventId() const {return eventId;}
        Int_t getNtracks4p0() const {return Ntracks4p0;}
        Int_t getNtracks2p4() const {return Ntracks2p4;}
        Int_t getNtracks0p8() const {return Ntracks0p8;}
        Double_t getpTssbar() const {return pTssbar;}
        std::vector<SmallTrack> getCandidates() const {return Candidates;}

        ClassDef(SmallEvent,1); 
};
#endif