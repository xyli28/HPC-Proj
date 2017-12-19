#include "OpenMM.h" 

#ifndef _GCMC_H_
#define _GCMC_H_


class Gcmc
{
public:
    //members
    double           temp;
    double           beta;
    double           mu;
    double           volume;
    double           boxEdgeLength;
    double           waveLengthCube;
    int              M;
    double           numReal;
    double           numGhost;
    std::vector<int> realList;
    std::vector<int> ghostList;
    OpenMM::System*  system;
    OpenMM::Context* context;
    OpenMM::Integrator* integrator;
    OpenMM::NonbondedForce* nonbond;
    int              mdSteps;
    int              upperLimit;
    int              lowerLimit;
    double           sigma;
    double           epsilon;
    double           charge;
    
    double           randIndex = 0.0;
    double           refIndex = 0.0;
    double           numInsr = 0.0;
    double           numDel = 0.0;
    double           sucInsr = 0.0;
    double           sucDel = 0.0;
    double           m = 0.0;           
    double           cutoff = 1.5*sigma;//cluster criteria
    double           f = 1.0;           //wang-laudau initial modification factor
    double           nIn = 0.0;         //coordination of one atom
    double           vIn = 4.0/3.0*M_PI*pow(cutoff,3);//excluded volume for each atom
    std::map<int,std::vector<int>> bond;
  
    double           *eta = new double[upperLimit*M+1];   
    double           *lamda = new double[M+1];
    double           *histogram = new double[upperLimit*M + 1];
    double           *energy = new double[upperLimit*M + 1];
    double           *freeEnergy = new double[upperLimit*M + 1];

public:
    //constructors    
    Gcmc() = default;
    Gcmc(double temp,double beta,double mu, double volume, double boxEdgeLength,
         double waveLengthCube, int M, double numReal, double numGhost, std::vector<int> realList,
         std::vector<int> ghostList, OpenMM::System *system, OpenMM::Context* context, 
         OpenMM::Integrator* integrator, OpenMM::NonbondedForce* nonbond, int mdSteps,
         int upperLimit, int lowerLimit, double sigma, double epsilon, double charge,
         double* lamda, double* histogram, double* energy, double* freeEnergy, double* eta);
    ~Gcmc();

public: 
    //functions
    void updateForce(double scaling);
    double getDistance(const std::vector<OpenMM::Vec3>& positions, int atomA, int atomB);
    
    //std::map<int>& generateGraph(std::vector<OpenMM::Vec3>& positions)
     
    void step();

    void generateGraph(const std::vector<OpenMM::Vec3>& positions);

    int clusterSizebyDFS(int start);
     
    bool checkClusterCriteria(const std::vector<OpenMM::Vec3>& positions);

    void  setRandomPositions();

    bool getRandomAtoms();

    void setRandomVelocities();

    double getInserProbability();

    double getDelProbability();
	
    double getAverageNum();

    void rezeroState();

    void rezeroHistogram();

    double* getEnergy();

    void rezeroEnergy();

    double* getFreeEnergy();
};

#endif


