#include <iostream>
#include <algorithm>
#include <omp.h>
#include "Gcmc.h"
#include "OpenMM.h"

Gcmc::Gcmc(double temp,double beta,double mu, double volume, double boxEdgeLength,
    double waveLengthCube, int M, int numReal, int numGhost, std::vector<int> realList,
    std::vector<int> ghostList, OpenMM::System *system, OpenMM::Context* context, 
    OpenMM::Integrator* integrator, OpenMM::NonbondedForce* nonbond, int mdSteps,
    int upperLimit, int lowerLimit, double sigma, double epsilon, double charge) :
    temp(temp),beta(beta),mu(mu),volume(volume),boxEdgeLength(boxEdgeLength),
    waveLengthCube(waveLengthCube),M(M),numReal(numReal),numGhost(numGhost),realList(realList),
    ghostList(ghostList), system(system), context(context), integrator(integrator),
    nonbond(nonbond), mdSteps(mdSteps), upperLimit(upperLimit),
    lowerLimit(lowerLimit), sigma(sigma), epsilon(epsilon),charge(charge){

    eta = new double[upperLimit*M+1];   
    lamda = new double[M+1];
    histogram = new double[upperLimit*M + 1];
    energy = new double[upperLimit*M + 1];
    freeEnergy = new double[upperLimit*M + 1];


    for (int i = 0; i < M; i++)
        lamda[i] = i*i*1.0/(M*M);
    lamda[M] = 1.0;
 
    for (int i = 0; i < upperLimit*M + 1; i++){
        histogram[i] = 0.0;
        eta[i] = 0.0;
        freeEnergy[i] = 0.0;
        energy[i] = 0.0;
    }
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

//#pragma omp parallel
{
//#pragma omp for  
    for (int i = 0; i < size; i++)
        for (int j = i + 1; j < size; j++ )
            if (getDistance(positions,realList[i],realList[j]) < cutoff){
//#pragma omp critical
                bond[realList[i]].push_back(realList[j]);
//#pragma omp critical
                bond[realList[j]].push_back(realList[i]);
            }
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
#pragma omp parallel
{
#pragma omp for 
    for (int i = 0; i< realList.size(); i++)
        if (realList[i] != refIndex)
            if (getDistance(positions,realList[i],refIndex) < cutoff)
#pragma omp atomic
                nIn += 1;
}   
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
#pragma omp parallel
{
#pragma omp for
    for (int i = 0; i < realList.size(); i++)
        if (realList[i] != refIndex)
            if (getDistance(positions,realList[i],refIndex) < cutoff){
#pragma omp atomic
                nIn += 1.0;
#pragma omp critical
                randList.push_back(realList[i]);
            }
}
    int d = realList.size();    
    realList.erase(std::remove(realList.begin(),realList.end(),
        randIndex),realList.end());    
    if (d!=realList.size()+1)
 
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

void Gcmc::insertion(){
    if (numReal == upperLimit){
    } else {
        numInsr += 1.0;
        double oldPotential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
        if (m == 0){
            int randNumber = rand()%ghostList.size();
            randIndex = ghostList[randNumber];
            ghostList.erase(ghostList.begin()+randNumber);
            setRandomPositions();
        }
        updateForce(lamda[m+1]);
        double newPotential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
        double dU = newPotential - oldPotential;
        double acc = exp(-beta*(dU-(lamda[m+1]-lamda[m])*mu)+eta[numReal*M+m+1]-eta[numReal*M+m])
                     * pow(vIn*numReal/waveLengthCube,lamda[m+1]-lamda[m])
                     * pow((numReal+lamda[m])*(nIn+lamda[m]),lamda[m])/
                     pow((numReal+lamda[m+1])*(nIn+lamda[m+1]),lamda[m+1]);
        acc = std::min(acc,1.0);
        if (rand()/double(RAND_MAX) < acc){
            sucInsr += 1.0;
            if (m == 0){
                m = 1;
                setRandomVelocities();
            } else if (m == M-1) {
                m = 0;
                realList.push_back(randIndex);
                numReal += 1;
                numGhost -= 1;
            } else {
                m += 1;
            }
        } else {
            updateForce(lamda[m]);
            if (m == 0)
                ghostList.push_back(randIndex);
        }       
    }
}

void Gcmc::deletion(){
    if (numReal == lowerLimit && m == 0){
    } else {
        numDel += 1;
        double oldPotential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
       
        if (m == 0){
            if (getRandomAtoms())  
                updateForce(lamda[M-1]);
            else
                return; 
        } else 
            updateForce(lamda[m-1]);
        double newPotential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
        double dU = newPotential - oldPotential;
        double acc = 1.0;
        if (m == 0)
            acc = exp(-beta*(dU-(lamda[M-1]-1.0)*mu)+eta[numReal*M+m-1]-eta[numReal*M+m])
                  *pow((numReal-1)*vIn/waveLengthCube,lamda[M-1]-1.0)
                  *pow(numReal*nIn,1.0)/
                  pow((numReal-1.0+lamda[M-1])*(nIn-1.0+lamda[M-1]),lamda[M-1]);
        else
            acc = exp(-beta*(dU-(lamda[m-1]-lamda[m])*mu)+eta[numReal*M+m-1]-eta[numReal*M+m])
                  *pow(numReal*vIn/waveLengthCube,lamda[m-1]-lamda[m])
                  *pow((numReal+lamda[m])*(nIn+lamda[m]),lamda[m])/
                  pow((numReal+lamda[m-1])*(nIn+lamda[m-1]),lamda[m-1]);
        acc = std::min(1.0,acc);
        if (rand()/double(RAND_MAX) < acc){
            sucDel += 1;
            if (m == 0){
                m = M-1;
                numReal -= 1;
                nIn -= 1;
                numGhost += 1;
            } else if (m == 1){
                m = 0;
                ghostList.push_back(randIndex);
            } else 
                m -= 1;
        } else 
            if (m == 0){
                updateForce(1.0);
                realList.push_back(randIndex);
            } else
                updateForce(lamda[m]); 
    }
}
                 
                
 
void Gcmc::step(){
    
    if (rand()/double(RAND_MAX) < 0.5)
        insertion();
    else
        deletion();

    histogram[numReal*M+m] += 1.0;
    energy[numReal*M+m] += context->getState(OpenMM::State::Energy).getPotentialEnergy();
    eta[numReal*M+m] -= f;
    //integrator->step(mdSteps);
}

double Gcmc::getInsrProbability(){
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
    rezeroEnergy();
}

void Gcmc::rezeroHistogram(){
    for (int i = 0; i < upperLimit*M + 1; i++){
        histogram[i] = 0.0;
        energy[i] = 0.0;
    }
}

void Gcmc::updateWangLandauFactor(){
    double sumHistogram = 0.0;
    for (int i = lowerLimit*M; i < upperLimit*M; i++)
        sumHistogram += histogram[i];
    double aveHistogram = sumHistogram/((upperLimit-lowerLimit)*M+1);
    int count = 0;
    for (int i = lowerLimit*M; i < upperLimit*M; i++)
        if (histogram[i] > 0.8*aveHistogram)
            count++;  
    if (count == (upperLimit-lowerLimit)*M+1 ){
        f = 0.5*f;
        rezeroHistogram();
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
