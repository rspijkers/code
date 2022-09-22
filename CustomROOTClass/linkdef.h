// All Linkdef.h files start with this preamble
#ifdef __CINT__
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link C++ nestedclasses;

// List all your classes that require Dictionary generation:
#pragma link C++ class TestEvent+;
// #pragma link C++ class TTree::Branch<TestEvent>+;

#endif