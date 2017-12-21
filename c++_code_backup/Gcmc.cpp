#include <iostream>
#include <algorithm>
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
    dSq += pow((dVec[0]-boxEdgeLength*round(dVec[0]/boxEdgeLength)),2);
    dSq += pow((dVec[1]-boxEdgeLength*round(dVec[1]/boxEdgeLength)),2);
    dSq += pow((dVec[2]-boxEdgeLength*round(dVec[2]/boxEdgeLength)),2);
    d = sqrt(dSq);
    return d;
}

void Gcmc::generateGraph(const std::vector<OpenMM::Vec3>& positions){
    bond.clear();    
    int size = realList.size();    
    std::cout << "generate graph" << std::endl;
    //for (int i = 0; i < size; i++)
    //    bond[i] = *new std::vector<int>();

    for (int i = 0; i < size; i++)
        for (int j = i + 1; j < size; j++ )
            if (getDistance(positions,realList[i],realList[j]) <= cutoff){
                bond[realList[i]].push_back(realList[j]);
                bond[realList[j]].push_back(realList[i]);
            }
    std::cout << "print bond" << std::endl;
    for (int i = 0; i < size; i++)
        std::cout << bond[realList[i]].size() << std::endl;
    std::cout << "bond end" << std::endl;
}

int Gcmc::clusterSizebyDFS(int start){
    std::vector<int> visited = {10000};
    std::vector<int> stack = {start};
  
    while (!stack.empty()){
        int vertex = stack.back();
        stack.pop_back();
        if (std::find(visited.begin(),visited.end(),vertex) == visited.end()){
            visited.push_back(vertex);
            for (int i = 0; i < bond[vertex].size(); i++)
                if (std::find(visited.begin(),visited.end(),bond[vertex][i]) 
                    == visited.end())
                    stack.push_back(bond[vertex][i]);
        }
    }
    std::cout << realList.size() << "vs" << visited.size()-1 << std::endl;
    return visited.size()-1;
}
    
bool Gcmc::checkClusterCriteria(const std::vector<OpenMM::Vec3>& positions){
    generateGraph(positions);
    int size = clusterSizebyDFS(realList[0]);
    if (size == realList.size())
        return true;
    else
        return false;
}

void Gcmc::setRandomPositions(){
    const std::vector<OpenMM::Vec3>& positions = 
        context->getState(OpenMM::State::Positions).getPositions();
    refIndex = realList[rand()%realList.size()];
    std::cout << "get distance after insert" << getDistance(positions,refIndex,randIndex) << std::endl; 
    std::cout << "check positions" << positions[randIndex][0] << " "<<positions[refIndex][0]<< std::endl; 
    std::cout << "check positions" << positions[randIndex][1] << " "<<positions[refIndex][1]<< std::endl; 
    std::cout << "check positions" << positions[randIndex][2] << " "<<positions[refIndex][2]<< std::endl; 
 
    nIn = 0.0;
    for (int i = 0; i< realList.size(); i++)
        if (realList[i] != refIndex)
            if (getDistance(positions,realList[i],refIndex) <= cutoff)
                nIn += 1;
    double x,y,z;
    do {
        x = (2*rand()/double(RAND_MAX)-1);
        y = (2*rand()/double(RAND_MAX)-1);
        z = (2*rand()/double(RAND_MAX)-1);
    } while (sqrt(pow(x,2)+pow(y,2)
        +pow(z,2)) > 1);
    
    std::cout << "insertDistance" <<sqrt(pow(x,2)+pow(y,2)  +pow(z,2)) << std::endl;
    std::cout << "ref/random " << refIndex << "  " << randIndex << std::endl; 
    x *= cutoff;    
    y *= cutoff;    
    z *= cutoff;    

    std::cout << x <<" " <<y << " " << z << std::endl; 
    
     
    std::vector<OpenMM::Vec3> newPositions;
    for (int i = 0; i < positions.size(); i++)
        newPositions.push_back(positions[i]);
    double a = newPositions[refIndex][0];
    double b = newPositions[refIndex][1];
    double c = newPositions[refIndex][2];
    newPositions[randIndex][0] = x + a;
    newPositions[randIndex][1] = y + b;
    newPositions[randIndex][2] = z + c;
    
    context->setPositions(newPositions);
    const std::vector<OpenMM::Vec3>& newnewpositions = 
        context->getState(OpenMM::State::Positions).getPositions();
    std::cout << "get distance after insert" << getDistance(positions,refIndex,randIndex) << std::endl; 
    std::cout << "get distance after insert" << getDistance(newPositions,refIndex,randIndex) << std::endl; 
    std::cout << "get distance after insert" << getDistance(newnewpositions,refIndex,randIndex) << std::endl; 
    std::cout << "check positions" << positions[randIndex][0] << " "<<positions[refIndex][0]<< std::endl; 
    std::cout << "check positions" << newPositions[randIndex][0] << " " << newPositions[refIndex][0]<<std::endl; 
    std::cout << "check positions" << newnewpositions[randIndex][0] << " " << newnewpositions[refIndex][0]<<std::endl; 
    std::cout << "check positions" << positions[randIndex][1] << " "<<positions[refIndex][1]<< std::endl; 
    std::cout << "check positions" << newPositions[randIndex][1] << " " << newPositions[refIndex][1]<<std::endl; 
    std::cout << "check positions" << newnewpositions[randIndex][1] << " " << newnewpositions[refIndex][1]<<std::endl; 
    std::cout << "check positions" << positions[randIndex][2] << " "<<positions[refIndex][2]<< std::endl; 
    std::cout << "check positions" << newPositions[randIndex][2] << " " << newPositions[refIndex][2]<<std::endl; 
    std::cout << "check positions" << newnewpositions[randIndex][2] << " " << newnewpositions[refIndex][2]<<std::endl; 
    for (int i = 0; i < realList.size(); i++)
        std::cout <<" "<< realList[i];
    std::cout << std::endl;
    newPositions.clear(); 
}
   
bool Gcmc::getRandomAtoms(){
    const std::vector<OpenMM::Vec3>& positions =
        context->getState(OpenMM::State::Positions).getPositions(); 

    bool goodAtoms = false;
    nIn = 0.0;
    std::vector<int> randList;
    
    std::cout << "reallistsize "<<realList.size() << std::endl; 
    refIndex = realList[rand()%realList.size()];
    for (int i = 0; i < realList.size(); i++)
        if (realList[i] != refIndex)
            if (getDistance(positions,realList[i],refIndex) <= cutoff){
                std::cout <<"distance " << getDistance(positions,realList[i],refIndex)<< std::endl;
                nIn += 1.0;
                randList.push_back(realList[i]);
            } else 
                std::cout <<"distance " << getDistance(positions,realList[i],refIndex)<< std::endl;
   
    std::cout << "randlistsize "<<randList.size() << std::endl; 
    randIndex = randList[rand()%randList.size()];
    realList.erase(std::remove(realList.begin(),realList.end(),
        randIndex),realList.end()); 
    std::cout << "reallistsize "<<realList.size() << std::endl;   
 
    if (!checkClusterCriteria(positions)){
        realList.push_back(randIndex); 
        std::cout << "reallistsize "<<realList.size() << std::endl;
        std::cout << "cluster fall apart"<< std::endl;
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
    context->setVelocities(newNewVelocities);
}

void Gcmc::insertion(){
    
    std::cout <<"insertion" << randIndex <<std::endl;
    std::cout << "m and lamdam" << m << " " << lamda[m] << std::endl;
    if (numReal == upperLimit){
    } else {
        numInsr += 1.0;
        double oldPotential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
        std::cout <<"old " <<oldPotential << std::endl;
        if (m == 0){
            int randNumber = rand()%ghostList.size();
            randIndex = ghostList[randNumber];
            ghostList.erase(ghostList.begin()+randNumber);
            setRandomPositions();
        }
        for (int i = 0; i < realList.size(); i++)
            std::cout <<" "<< realList[i];
        std::cout << std::endl; 
        updateForce(lamda[m+1]);
        double newPotential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
        std::cout <<"new " <<newPotential << std::endl;
        double dU = newPotential - oldPotential;

        std::cout << "dU " << dU << std::endl;
        double acc = exp(-beta*1000*(dU-(lamda[m+1]-lamda[m])*mu)+eta[numReal*M+m+1]-eta[numReal*M+m])
                     * pow(vIn*numReal/waveLengthCube,lamda[m+1]-lamda[m])
                     * pow((numReal+lamda[m])*(nIn+lamda[m]),lamda[m])/
                     pow((numReal+lamda[m+1])*(nIn+lamda[m+1]),lamda[m+1]);
        std::cout << "acc " << acc << std::endl;
        acc = std::min(acc,1.0);
        if (rand()/double(RAND_MAX) < acc){
            std::cout <<"suc" << std::endl;
            sucInsr += 1.0;
            if (m == 0){
                m = 1;
                setRandomVelocities();
                const std::vector<OpenMM::Vec3>& newnewpositions = 
                    context->getState(OpenMM::State::Positions).getPositions();
                std::cout << "get distance after suc of insertion" << getDistance(newnewpositions,refIndex,randIndex) << std::endl; 
            } else if (m == M-1) {
                m = 0;
                realList.push_back(randIndex);
                numReal += 1;
                numGhost -= 1;
            } else {
                m += 1;
            }
        } else {
            std::cout <<"fail" << std::endl;
            std::cout <<"do it"<< std::endl;
            double potential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
            std::cout << "before dorecover force after insertion " << potential << std::endl;
            updateForce(lamda[m]);
            std::cout <<"done it"<< std::endl;
            potential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
            std::cout << "after done dorecover force after insertion " << potential << std::endl;
            std::cout << "m and lamdam" << m << " " << lamda[m] << std::endl;
            if (m == 0)
                ghostList.push_back(randIndex);
            potential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
            std::cout << "1recover force after insertion " << potential << std::endl;
        }       
    }
}

void Gcmc::deletion(){
    std::cout <<"deletion" << randIndex << std::endl;
    if (numReal == lowerLimit && m == 0){
    } else {
        numDel += 1;
        double oldPotential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
        std::cout <<"old " <<oldPotential << std::endl;
       
        if (m == 0){
            if (getRandomAtoms())  
                updateForce(lamda[M-1]);
            else{
                std::cout << "can we get here" << std::endl;
                return; 
            }
        } else 
            updateForce(lamda[m-1]);
        double newPotential = context->getState(OpenMM::State::Energy).getPotentialEnergy();

        std::cout <<"new " <<newPotential << std::endl;
        double dU = newPotential - oldPotential;
        std::cout << "dU " << dU << std::endl;
        double acc = 1.0;
        if (m == 0)
            acc = exp(-beta*1000*(dU-(lamda[M-1]-1.0)*mu)+eta[numReal*M+m-1]-eta[numReal*M+m])
                  *pow((numReal-1)*vIn/waveLengthCube,lamda[M-1]-1.0)
                  *pow(numReal*nIn,1.0)/
                  pow((numReal-1.0+lamda[M-1])*(nIn-1.0+lamda[M-1]),lamda[M-1]);
        else
            acc = exp(-beta*1000*(dU-(lamda[m-1]-lamda[m])*mu)+eta[numReal*M+m-1]-eta[numReal*M+m])
                  *pow(numReal*vIn/waveLengthCube,lamda[m-1]-lamda[m])
                  *pow((numReal+lamda[m])*(nIn+lamda[m]),lamda[m])/
                  pow((numReal+lamda[m-1])*(nIn+lamda[m-1]),lamda[m-1]);
        std::cout << "acc " << acc << std::endl;
        acc = std::min(1.0,acc);
        if (rand()/double(RAND_MAX) < acc){
            std::cout <<"suc" << std::endl;
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
        } else{ 
            std::cout <<"fail" << std::endl;
            if (m == 0){
                updateForce(1.0);
                double potential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
                std::cout << "1recover force after deletion " << potential << std::endl;
                realList.push_back(randIndex);
            } else{
                updateForce(lamda[m]); 
                double potential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
                std::cout << "2recover force after deletion " << potential << std::endl;
            }
        }   
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
    double potential = context->getState(OpenMM::State::Energy).getPotentialEnergy();
    std::cout <<"after step potential " << potential << std::endl;
    std::cout <<"after step, check cluster again" << std::endl;
    for (int i = 0; i < realList.size(); i++)
        std::cout <<" "<< realList[i];
    std::cout << std::endl; 
    const std::vector<OpenMM::Vec3>& positions = 
        context->getState(OpenMM::State::Positions).getPositions();
    if (checkClusterCriteria(positions))
        std::cout << "fine" << std::endl;
    else
        std::cout << "not fine" << std::endl;
    std::cout << "randIndex after step" << randIndex << std::endl;
    std::cout << "final distance after step" << " "<< refIndex <<" " << randIndex <<" "<< getDistance(positions,randIndex,refIndex) << std::endl;
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
