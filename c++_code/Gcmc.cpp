#include <algorithm>
#include "Gcmc.h"
#include "OpenMM.h"

Gcmc::Gcmc(double temp,double beta,double mu, double volume, double boxEdgeLength,
    double waveLengthCube, int M, double numReal, double numGhost, std::vector<int> realList,
    std::vector<int> ghostList, OpenMM::System *system, OpenMM::Context* context, 
    OpenMM::Integrator* integrator, OpenMM::NonbondedForce* nonbond, int mdSteps,
    int upperLimit, int lowerLimit, double sigma, double epsilon, double charge,
    double* lamda, double* histogram, double* energy, double* freeEnergy, double* eta):
    temp(temp),beta(beta),mu(mu),volume(volume),boxEdgeLength(boxEdgeLength),
    waveLengthCube(waveLengthCube),M(M),numReal(numReal),numGhost(numGhost),realList(realList),
    ghostList(ghostList), system(system), context(context), integrator(integrator),
    nonbond(nonbond), mdSteps(mdSteps), upperLimit(upperLimit),
    lowerLimit(lowerLimit), sigma(sigma), epsilon(epsilon),charge(charge), lamda(lamda),
    histogram(histogram), energy(energy), freeEnergy(freeEnergy), eta(eta){

    for (int i = 0; i < upperLimit*M + 1; i++)
        eta[upperLimit*M + 1] = 0.0;

    for (int i = 0; i < M; i++)
        lamda[i] = i*i*1.0/(M*M);
    lamda[M] = 1.0;
 
    for (int i = 0; i < upperLimit*M + 1; i++)
        histogram[upperLimit*M + 1] = 0.0;

    for (int i = 0; i < upperLimit*M + 1; i++)
        energy[upperLimit*M + 1] = 0.0;

    for (int i = 0; i < upperLimit*M + 1; i++)
        freeEnergy[upperLimit*M + 1] = 0.0;
}

Gcmc::~Gcmc(){
    delete lamda;
    delete histogram;
    delete energy;
    delete freeEnergy;
}

void Gcmc::updateForce(double scaling){
    nonbond->setParticleParameters(randIndex,charge*scaling,
        sigma, epsilon*pow(scaling,2.0));
    nonbond->updateParametersInContext(*context);
}

double Gcmc::getDistance(const std::vector<OpenMM::Vec3>& positions, int atomA, int atomB){
    double dSq = 0.0;
    double d = 0.0;
    OpenMM::Vec3 dVec = positions[atomA] - positions[atomB];
    dSq = pow((dVec[0]-boxEdgeLength*round(dVec[0]/boxEdgeLength)),2);
    dSq = pow((dVec[1]-boxEdgeLength*round(dVec[1]/boxEdgeLength)),2);
    dSq = pow((dVec[2]-boxEdgeLength*round(dVec[2]/boxEdgeLength)),2);
    d = sqrt(dSq);
    return d;
}

void Gcmc::generateGraph(const std::vector<OpenMM::Vec3>& positions){
    bond.clear();    
    int size = realList.size();    

    //for (int i = 0; i < size; i++)
    //    bond[i] = *new std::vector<int>();

    for (int i = 0; i < size; i++)
        for (int j = i; j < size; j++ )
            if (getDistance(positions,i,j) < cutoff){
                bond[i].push_back(j);
                bond[j].push_back(i);
            }
}

int Gcmc::clusterSizebyDFS(int start){
    std::vector<int> visited;
    std::vector<int> stack = {start};
  
    while (!stack.empty()){
        int vertex = stack.back();
        stack.pop_back();
        if (std::find(visited.begin(),visited.end(),vertex) != visited.end()){
            visited.push_back(vertex);
            for (int i = 0; i < bond[vertex].size(); i++)
                if (std::find(visited.begin(),visited.end(),bond[vertex][i]) 
                    != visited.end())
                    stack.push_back(bond[vertex][i]);
        }
    }
    return visited.size();
}
    
bool Gcmc::checkClusterCriteria(const std::vector<OpenMM::Vec3>& positions){
    generateGraph(positions);
    int size = clusterSizebyDFS(realList[0]);
    if (size = realList.size())
        return true;
    else
        return false;
}

void  Gcmc::setRandomPositions(){
    const std::vector<OpenMM::Vec3>& positions = 
        context->getState(OpenMM::State::Positions).getPositions();
    refIndex = realList[rand()%realList.size()];
 
    nIn = 0.0;
    for (int i = 0; i< realList.size(); i++)
        if (realList[i] != refIndex)
            if (getDistance(positions,i,refIndex) < cutoff)
                nIn += 1;
    
    OpenMM::Vec3 randVec = *new OpenMM::Vec3;
    do {
        randVec[0] = (rand()/double(RAND_MAX)*2-1)*cutoff;
        randVec[1] = (rand()/double(RAND_MAX)*2-1)*cutoff;
        randVec[2] = (rand()/double(RAND_MAX)*2-1)*cutoff;
    } while (sqrt(pow(randVec[0],2)+pow(randVec[1],2)
        +pow(randVec[2],2)) <= cutoff);

    std::vector<OpenMM::Vec3> newPositions(positions);
    newPositions[randIndex] = newPositions[randIndex] + randVec;
    context->setPositions(newPositions);
}
   
bool Gcmc::getRandomAtoms(){
    const std::vector<OpenMM::Vec3>& positions =
    context->getState(OpenMM::State::Positions).getPositions(); 

    bool goodAtoms = false;
    nIn = 0.0;
    std::vector<int> randList;
    refIndex = realList[rand()%realList.size()];
    for (int i = 0; i < realList.size(); i++)
        if (realList[i] != refIndex)
            if (getDistance(positions,i,refIndex) < cutoff){
                nIn += 1.0;
                randList.push_back(i);
            }
    randIndex = randList[rand()%randList.size()];
    realList.erase(std::remove(realList.begin(),realList.end(),
        randIndex),realList.end());      
 
    if (!checkClusterCriteria(positions)){
        realList.push_back(randIndex); 
        return false;
    } else
        return true;
}  

void Gcmc::setRandomVelocities(){
    const std::vector<OpenMM::Vec3>& velocities =
        context->getState(OpenMM::State::Velocities).getVelocities();         
    context->setVelocitiesToTemperature(temp,rand());
    const std::vector<OpenMM::Vec3>& newVelocities =
        context->getState(OpenMM::State::Velocities).getVelocities();         
    std::vector<OpenMM::Vec3> newNewVelocities(velocities);
    newNewVelocities[randIndex] = newVelocities[randIndex];
    context->setVelocities(newVelocities);
} 

void Gcmc::step(){
    integrator->step(mdSteps);
}

double Gcmc::getInserProbability(){
    return sucInsr/numInsr;
}

double Gcmc::getDelProbability(){
    return sucDel/numDel;
}

void Gcmc::rezeroState(){
    numInsr = 0.0;
    numDel  = 0.0;
    sucInsr = 0.0;
    sucInsr = 0.0;
    rezeroHistogram();
}

void Gcmc::rezeroHistogram(){
    for (int i = 0; i < upperLimit*M + 1; i++){
        histogram[i] = 0.0;
        energy[i] = 0.0;
    }
}

double* Gcmc::getEnergy(){
    for (int i = lowerLimit*M; i < upperLimit*M + 1; i++)
        if (histogram[i] > 0.0)
            energy[i] = energy[i]/histogram[i];
    return energy;
}

void Gcmc::rezeroEnergy(){
    for (int i = lowerLimit*M; i < upperLimit*M + 1; i++)
        energy[i] = 0.0;
}

double* Gcmc::getFreeEnergy(){
    for (int i = lowerLimit*M; i < upperLimit*M + 1; i++)
        freeEnergy[i] = -(log(histogram[i]/histogram[lowerLimit*M])-eta[i]+
            eta[lowerLimit*M])/beta;
    return freeEnergy;
}
