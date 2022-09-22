#include "TClass.h"

class TestEvent {
    public:
        Int_t eventId; // eventId
        TestEvent() {eventId = 9999;}
        TestEvent(int i);
        ~TestEvent() {;}
        ClassDef(TestEvent,1); // TestEvent class
};