#include <iostream>
// #include "classtest.h"
#include "TTree.h"
#include "TFile.h"

class TestEvent {
    public:
        Int_t eventId; // eventId
        TestEvent() {}
        TestEvent(int i) {eventId = i;}
        ClassDef(TestEvent,1); // TestEvent class
};

int test(){
    TFile *file = new TFile("testoutput.root", "RECREATE");
    TTree *tree = new TTree("events", "event");
    TestEvent ev;
    tree->Branch("event", &ev);
    for (int i = 0; i < 10; i++){
        ev = TestEvent(i);
        tree->Fill();
    }
	file->Write();
    std::cout << ev.eventId << std::endl;
    return 0;
}

int main(){
    return test();
}