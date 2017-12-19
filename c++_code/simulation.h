#ifndef _SIMULATION_H_
#define _SIMULATION_H_

// -----------------------------------------------------------------------------
//                           INTERFACE TO OpenMM
// -----------------------------------------------------------------------------
// These four functions and an opaque structure are used to interface our main
// program with OpenMM without the main program having any direct interaction
// with the OpenMM API. This is a clean approach for interfacing with any MD
// code, although the details of the interface routines will differ. This is
// still just "locally written" code and is not required by OpenMM.

struct AtomType;
struct MyAtomInfo;
struct MyOpenMMData;
static MyOpenMMData* 
myInitializeOpenMM( const MyAtomInfo    atoms[],
                    double              temperature,
                    double              frictionInPerPs,
                    double              stepSizeInFs,
                    double              boxEdgeLengthInNm, 
                    std::string&        platformName); 
static void          myStepWithOpenMM(MyOpenMMData*, int numSteps);
static void          myGetOpenMMState(MyOpenMMData*, bool wantEnergy,
                                      double& time, double& energy, 
                                      MyAtomInfo atoms[]);
static void          myTerminateOpenMM(MyOpenMMData*);
static void
myWritePDBFrame(int frameNum, double timeInPs, double energyInKcal, 
                const MyAtomInfo atoms[]); 


#endif
