/* -----------------------------------------------------------------------------
 *               OpenMM(tm) HelloEthane example in C++ (June 2009)
 * -----------------------------------------------------------------------------
 * This is a complete, self-contained "hello world" example demonstrating 
 * GPU-accelerated simulation of a system with both bonded and nonbonded forces, 
 * using ethane (H3-C-C-H3) in a vacuum as an example. A multi-frame PDB file is 
 * written to stdout which can be read by VMD or other visualization tool to 
 * produce an animation of the resulting trajectory.
 *
 * Pay particular attention to the handling of units in this example. Incorrect
 * handling of units is a very common error; this example shows how you can
 * continue to work with Amber-style units of Angstroms and kCals while correctly
 * communicating with OpenMM in nanometers and kJoules.
 * -------------------------------------------------------------------------- */

#include <iostream>
#include <cstdio>
#include <string>
#include <vector>
#include <math.h>
#include <typeinfo>
#include "Gcmc.h"
#include "OpenMM.h"

// -----------------------------------------------------------------------------
//                                 MOCK MD CODE
// -----------------------------------------------------------------------------
// The code starting here and through main() below is meant to represent in 
// simplified form some pre-existing molecular dynamics code, which defines its 
// own data structures for force fields, the atoms in this simulation, and the 
// simulation parameters, and takes care of recording the trajectory. All this 
// has nothing to do with OpenMM; the OpenMM-dependent code comes later and is 
// clearly marked below.
// -----------------------------------------------------------------------------

//                   MODELING AND SIMULATION PARAMETERS
const double StepSizeInFs        = 2;       // integration step size (fs)
const double ReportIntervalInFs  = 2000;      // how often to generate PDB frame (fs)
const double SimulationTimeInPs  = 40;     // total simulation time (ps)
static const bool   WantEnergy   = true;

const double Temperature         = 88.9407; // Kelvins
const double FrictionInPerPs     = 1.0;     // collisions per picosecond
const double CutoffDistanceInAng = 8.5125;  // Angstroms
const double BoxEdgeLengthInNm   = 4.0;     //nms
const double VolumeInNm3          =pow(BoxEdgeLengthInNm,3); //nm3

const double R     = 8.3144598;             //ideal gas constant in J/(K*Mol)
const double Beta  = 1.0/(Temperature*R);//reciprocal temperature 
const double Mu    = -7.7946;               //chemical potential in kilojoule/mole
const double WaveLength = 0.030802;         //thermal wavelength in nms
const double WaveLengthCube = pow(WaveLength,3.0);


//                            FORCE FIELD DATA
// For this example we're using a tiny subset of the Amber99 force field.
// We want to keep the data in the original unit system to avoid conversion
// bugs; this requires conversion on the way in and out of OpenMM.

// Amber reduces nonbonded forces between 1-4 bonded atoms.
const double Coulomb14Scale      = 0.5;
const double LennardJones14Scale = 0.5;

const double Charge              = 0.0;     //charge for argon 
const double Sigma               = 0.3405;  //sigma for argon in nanometers
const double Epsilon             = 0.997967;//epsilon for argon in kJ/mol

//                             GCMC PARAMETERS
const int    M     = 2;                     //number of intermediate states + 1
const int    UpperLimit = 30.0;             //maximum num of atoms
const int    LowerLimit = 1.0;              //minimum num of atoms
const int    MdSteps = 10;                  //ms steps

//                            other related information
const std::string   RealName = "H";
const std::string   GhostName = "Hg";
const std::string   IntermediateName = "Hi";


struct AtomType {
    double mass, charge, sigmaInNms, epsilonInkJPerMol;
} atomType[] = {/*0 H*/ 39.9481, 0.0, 0.340500, 0.997967};
const int H = 0;

//                                MOLECULE DATA
// Now describe an ethane molecule by listing its atoms, bonds, angles, and 
// torsions. We'll provide a default configuration which centers the molecule 
// at (0,0,0) with the C-C bond along the x axis.

// Use this as an "end of list" marker so that we do not have to count; let the
// computer do that!
const int EndOfList=-1;

struct MyAtomInfo
{   int type; const char* pdb; double initPosInAng[3]; double posInAng[3];} 
atoms[] = 
    {/*0*/H,       " H  ",       {14.143,  8.683, 13.348},     {0,0,0},
     /*1*/H,       " Hg ",       {17.510,  9.476, 15.579},     {0,0,0},
     /*2*/H,       " Hg ",       {11.315, 11.646,  9.369},     {0,0,0}, 
     /*3*/H,       " Hg ",       {11.729, 15.249, 10.719},     {0,0,0},
     /*4*/H,       " Hg ",       {19.967, 10.689, 12.443},     {0,0,0},
     /*5*/H,       " Hg ",       {17.660,  7.968, 11.859},     {0,0,0}, 
     /*6*/H,       " Hg ",       {12.722, 12.157, 12.867},     {0,0,0},
     /*7*/H,       " Hg ",       {14.949,  9.732,  9.687},     {0,0,0},
     /*8*/H,       " Hg ",       {18.510, 13.489, 15.431},     {0,0,0},
     /*9*/H,       " Hg ",       {15.045, 11.846,  6.650},     {0,0,0}, 
     /*10*/H,      " Hg ",       {24.143,  8.683, 13.348},     {0,0,0},
     /*11*/H,      " Hg ",       {27.510,  9.476, 15.579},     {0,0,0},
     /*12*/H,      " Hg ",       {21.315, 11.646,  9.369},     {0,0,0}, 
     /*13*/H,      " Hg ",       {21.729, 15.249, 10.719},     {0,0,0},
     /*14*/H,      " Hg ",       {29.967, 10.689, 12.443},     {0,0,0},
     /*15*/H,      " Hg ",       {27.660,  7.968, 11.859},     {0,0,0},
     /*16*/H,      " Hg ",       {22.722, 12.157, 12.867},     {0,0,0}, 
     /*17*/H,      " Hg ",       {24.949,  9.732,  9.687},     {0,0,0},
     /*18*/H,      " Hg ",       {28.510, 13.489, 15.431},     {0,0,0},
     /*19*/H,      " Hg ",       {25.045, 11.846,  6.650},     {0,0,0}, 
     /*20*/H,      " Hg ",       {34.143,  8.683, 13.348},     {0,0,0},
     /*21*/H,      " Hg ",       {37.510,  9.476, 15.579},     {0,0,0},
     /*22*/H,      " Hg ",       {31.315, 11.646,  9.369},     {0,0,0}, 
     /*23*/H,      " Hg ",       {31.729, 15.249, 10.719},     {0,0,0},
     /*24*/H,      " Hg ",       {39.967, 10.689, 12.443},     {0,0,0},
     /*25*/H,      " Hg ",       {37.660,  7.968, 11.859},     {0,0,0},
     /*26*/H,      " Hg ",       {32.722, 12.157, 12.867},     {0,0,0}, 
     /*27*/H,      " Hg ",       {34.949,  9.732,  9.687},     {0,0,0},
     /*28*/H,      " Hg ",       {38.510, 13.489, 15.431},     {0,0,0},
     /*29*/H,      " Hg ",       {35.045, 11.846,  6.650},     {0,0,0}, 
     EndOfList};

// -----------------------------------------------------------------------------
//                           INTERFACE TO OpenMM
// -----------------------------------------------------------------------------
// These four functions and an opaque structure are used to interface our main
// program with OpenMM without the main program having any direct interaction
// with the OpenMM API. This is a clean approach for interfacing with any MD
// code, although the details of the interface routines will differ. This is
// still just "locally written" code and is not required by OpenMM.

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
myWritePDBFrame(int frameNum, double timeInPs, double energyInKJ, 
                const MyAtomInfo atoms[]); 

//                               PDB FILE WRITER
// Given state data, output a single frame (pdb "model") of the trajectory.
static void
myWritePDBFrame(int frameNum, double timeInPs, double energyInKJ, 
                const MyAtomInfo atoms[]) 
{
    // Write out in PDB format -- printf is so much more compact than formatted cout.
    printf("MODEL     %d\n", frameNum);
    printf("REMARK 250 time=%.3f ps; energy=%.3f kJ/mole\n", 
           timeInPs, energyInKJ);
    //for (int n=0; atoms[n].type != EndOfList; ++n)
    //    printf("ATOM  %5d %4s ETH     1    %8.3f%8.3f%8.3f  1.00  0.00\n", 
    //        n+1, atoms[n].pdb, 
    //        atoms[n].posInAng[0], atoms[n].posInAng[1], atoms[n].posInAng[2]);
    //printf("ENDMDL\n");
}

// -----------------------------------------------------------------------------
//                           OpenMM-USING CODE
// -----------------------------------------------------------------------------
// The OpenMM API is visible only at this point and below. Normally this would
// be in a separate compilation module; we're including it here for simplicity.
// -----------------------------------------------------------------------------

// Suppress irrelevant warnings from Microsoft's compiler.
#ifdef _MSC_VER
    #pragma warning(disable:4996)   // sprintf is unsafe 
#endif

using OpenMM::Vec3; // so we can just say "Vec3" below

// This is our opaque "handle" class containing all the OpenMM objects that
// must persist from call to call during a simulation. The main program gets 
// a pointer to one of these but sees it as essentially a void* since it 
// doesn't know the definition of this class.
struct MyOpenMMData {
    MyOpenMMData() : system(0), context(0), integrator(0), nonbond(0) {}
    ~MyOpenMMData() {delete context; delete integrator; delete system;}
    OpenMM::System*         system;
    OpenMM::Integrator*     integrator;
    OpenMM::Context*  context;
    OpenMM::NonbondedForce* nonbond;
};


// -----------------------------------------------------------------------------
//                           Argon MAIN PROGRAM
// -----------------------------------------------------------------------------
int main() {
    // ALWAYS enclose all OpenMM calls with a try/catch block to make sure that
    // usage and runtime errors are caught and reported.
    try {
        std::string   platformName;

        // Set up OpenMM data structures; returns OpenMM Platform name.
        MyOpenMMData* omm = myInitializeOpenMM(atoms, Temperature, 
                            FrictionInPerPs, StepSizeInFs, 
                            BoxEdgeLengthInNm, platformName);

        // Set up GCMC initial condition
        int numReal = 1;

        std::vector<int> realList,ghostList;
        for (int i = 0; i < numReal; i++)
            realList.push_back(i);
        for (int i = numReal; i < UpperLimit; i++){
            ghostList.push_back(i);
            omm->nonbond->setParticleParameters(i,0.0,Sigma,0.0);
        }
        omm->nonbond->setReactionFieldDielectric(0.0);
        omm->nonbond->setUseDispersionCorrection(false);
        omm->nonbond->updateParametersInContext(*omm->context);

        double*  lamda           = new double(M+1);
        double*  histogram       = new double(UpperLimit*M + 1);
        double*  eta             = new double(UpperLimit*M + 1);
        double*  energy          = new double(UpperLimit*M + 1);
        double*  freeEnergy      = new double(UpperLimit*M + 1);
  
        Gcmc *simulation = new Gcmc(Temperature, Beta, Mu, VolumeInNm3, BoxEdgeLengthInNm, 
                                    WaveLengthCube, M, numReal, UpperLimit - numReal, realList,
                                    ghostList, omm->system, omm->context, omm->integrator,
                                    omm->nonbond, MdSteps, UpperLimit, LowerLimit, Sigma,
                                    Epsilon, Charge, lamda, histogram, energy, freeEnergy,eta);
        std::cout << simulation->lamda[1] << std::endl;    
        simulation->step();
        double time, Energy;
        const OpenMM::State state = omm->context->getState(OpenMM::State::Positions);
        const std::vector<Vec3>& positions = state.getPositions();
        bool xx = simulation->checkClusterCriteria(positions);
        std::cout << positions[0][0] << xx << std::endl;
        myGetOpenMMState(omm, WantEnergy, time, Energy, atoms);
        simulation->step();
        myGetOpenMMState(omm, WantEnergy, time, Energy, atoms);
        simulation->step();
        myGetOpenMMState(omm, WantEnergy, time, Energy, atoms);
  
 
        //Gcmc gcmc(Temperature, Beta, Mu, VolumeInNm3, BoxEdgeLengthInNm, 
        //          WaveLengthCube, M, numReal, UpperLimit - numReal, realList,
        //          ghostList, omm->system, omm->context, omm->integrator,
        //          omm->nonbond, MdSteps, UpperLimit, LowerLimit, Sigma,
        //          Epsilon, Charge);

        // Run the simulation:
        //  (1) Write the first line of the PDB file and the initial configuration.
        //  (2) Run silently entirely within OpenMM between reporting intervals.
        //  (3) Write a PDB frame when the time comes.
        printf("REMARK  Using OpenMM platform %s\n", platformName.c_str());
        
        

//        const int NumSilentSteps = (int)(ReportIntervalInFs / StepSizeInFs + 0.5);
//        for (int frame=1; ; ++frame) {
//            double time, energy;
//            myGetOpenMMState(omm, WantEnergy, time, energy, atoms);
//            //myWritePDBFrame(frame, time, energy, atoms);
//
//            if (time >= SimulationTimeInPs)
//                break;
//            
//
//            myStepWithOpenMM(omm, NumSilentSteps);
//        } 
 
        // Clean up OpenMM data structures.
        //delete lamda,histogram,energy,freeEnergy;
        myTerminateOpenMM(omm);

        return 0; // Normal return from main.
    }

    // Catch and report usage and runtime errors detected by OpenMM and fail.
    catch(const std::exception& e) {
        printf("EXCEPTION: %s\n", e.what());
        return 1;
    }
}

// -----------------------------------------------------------------------------
//                      INITIALIZE OpenMM DATA STRUCTURES
// -----------------------------------------------------------------------------
// We take these actions here:
// (1) Load any available OpenMM plugins, e.g. Cuda and Brook.
// (2) Allocate a MyOpenMMData structure to hang on to OpenMM data structures
//     in a manner which is opaque to the caller.
// (3) Fill the OpenMM::System with the force field parameters we want to
//     use and the particular set of atoms to be simulated.
// (4) Create an Integrator and a Context associating the Integrator with
//     the System.
// (5) Select the OpenMM platform to be used.
// (6) Return the MyOpenMMData struct and the name of the Platform in use.
//
// Note that this function must understand the calling MD code's molecule and
// force field data structures so will need to be customized for each MD code.
static MyOpenMMData* 
myInitializeOpenMM( const MyAtomInfo    atoms[],
                    double              temperature,
                    double              frictionInPerPs,
                    double              stepSizeInFs,
                    double              boxEdgeLengthInNm, 
                    std::string&        platformName) 
{
    // Load all available OpenMM plugins from their default location.
    OpenMM::Platform::loadPluginsFromDirectory
       ("/home/xli646/local/openmm/lib/plugins");
  
    // Allocate space to hold OpenMM objects while we're using them.
    MyOpenMMData* omm = new MyOpenMMData();

    // Create a System and Force objects within the System. Retain a reference
    // to each force object so we can fill in the forces. Note: the System owns
    // the force objects and will take care of deleting them; don't do it yourself!
    OpenMM::System&                 system      = *(omm->system = new OpenMM::System());
    OpenMM::NonbondedForce&         nonbond     = *(omm->nonbond = new OpenMM::NonbondedForce());
    system.addForce(&nonbond);
    
    // Specify the atoms and their properties:
    //  (1) System needs to know the masses.
    //  (2) NonbondedForce needs charges,van der Waals properties (in MD units!).
    //  (3) Collect default positions for initializing the simulation later.
    std::vector<Vec3> initialPosInNm;
    for (int n=0; atoms[n].type != EndOfList; ++n) {
        const AtomType& atype = atomType[atoms[n].type];
        system.addParticle(atype.mass);
        nonbond.addParticle(atype.charge,
                            atype.sigmaInNms, 
                            atype.epsilonInkJPerMol); 
        // Convert the initial position to nm and append to the array.
        const Vec3 posInNm(atoms[n].initPosInAng[0] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[1] * OpenMM::NmPerAngstrom,
                           atoms[n].initPosInAng[2] * OpenMM::NmPerAngstrom);
        initialPosInNm.push_back(posInNm);
    }
   
    //create periodic box
    nonbond.setNonbondedMethod(OpenMM::NonbondedForce::CutoffPeriodic);
    nonbond.setCutoffDistance(CutoffDistanceInAng * OpenMM::NmPerAngstrom);
    system.setDefaultPeriodicBoxVectors(Vec3(boxEdgeLengthInNm,0,0),
                                        Vec3(0,boxEdgeLengthInNm,0), 
                                        Vec3(0,0,boxEdgeLengthInNm));        

    // Choose an Integrator for advancing time, and a Context connecting the
    // System with the Integrator for simulation. Let the Context choose the
    // best available Platform. Initialize the configuration from the default
    // positions we collected above. Initial velocities will be zero.
    omm->integrator = new OpenMM::LangevinIntegrator(temperature, frictionInPerPs, 
                                                     stepSizeInFs * OpenMM::PsPerFs);
    omm->integrator->setConstraintTolerance(0.00001);
    omm->context    = new OpenMM::Context(*omm->system, *omm->integrator);
    omm->context->setPositions(initialPosInNm);

    platformName = omm->context->getPlatform().getName();
    return omm;
}


// -----------------------------------------------------------------------------
//                     COPY STATE BACK TO CPU FROM OPENMM
// -----------------------------------------------------------------------------
static void
myGetOpenMMState(MyOpenMMData* omm, bool wantEnergy, 
                 double& timeInPs, double& energyInKJ,
                 MyAtomInfo atoms[])
{
    int infoMask = 0;
    infoMask = OpenMM::State::Positions;
    if (wantEnergy) {
        infoMask += OpenMM::State::Velocities; // for kinetic energy (cheap)
        infoMask += OpenMM::State::Energy;     // for pot. energy (expensive)
    }
    // Forces are also available (and cheap).

    const OpenMM::State state = omm->context->getState(infoMask);
    timeInPs = state.getTime(); // OpenMM time is in ps already

    // Copy OpenMM positions into atoms array and change units from nm to Angstroms.
    const std::vector<Vec3>& positionsInNm = state.getPositions();
    for (int i=0; i < (int)positionsInNm.size(); ++i)
        for (int j=0; j < 3; ++j)
            atoms[i].posInAng[j] = positionsInNm[i][j] * OpenMM::AngstromsPerNm;

    // If energy has been requested, obtain it and convert from kJ to kcal.
    energyInKJ = 0;
    std::cout << state.getPotentialEnergy() << "  " << state.getKineticEnergy() << std::endl;
    if (wantEnergy) 
        energyInKJ = (state.getPotentialEnergy() + state.getKineticEnergy());
}


// -----------------------------------------------------------------------------
//                     TAKE MULTIPLE STEPS USING OpenMM 
// -----------------------------------------------------------------------------
static void 
myStepWithOpenMM(MyOpenMMData* omm, int numSteps) {
    omm->integrator->step(numSteps);
}

// -----------------------------------------------------------------------------
//                     DEALLOCATE OpenMM OBJECTS
// -----------------------------------------------------------------------------
static void 
myTerminateOpenMM(MyOpenMMData* omm) {
    delete omm;
}





