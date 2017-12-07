#use sterling definition of cluster, 1.5 sigma as cutoff, no innerShell
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from math import exp
from math import log
from math import sqrt
from math import pi
from random import random
from random import seed
from random import randrange
from random import randint
import numpy as np
from rotation import random_rotation_matrix
from randvec import random_vector_inside_sphere
from randvec import random_vector_inside_shell

class Gcmc(object):
    """Gcmc provides a object to run hybrid GCMC/MD simulation using openMM. 
    
    """

    def __init__(self,parameters):
        """Intialize an gcmc object using the parameters passed by parameters list
         
        Parameters
        ---------
        parameters:
            parameter list
        """
        #input variables
        self.temp = parameters['temp']                     #temperature
        self.beta = parameters['beta']                     #reciprocal temperature beta 
        self.mu = parameters['mu']                         #chemical potential
        self.vol = parameters['vol']                       #volume
        self.box = parameters['box']                       #simulation box, 3D tuple
        self.atoms = parameters['atoms']                   #atom list
        self.wavelengthCube = parameters['wavelengthCube'] #thermal de Broglie wavelength cube
        self.M = parameters['M']                           #amount of intermediate states
        self.numReal = parameters['numReal']               #number of real particles
        self.numGhost = parameters['numGhost']             #number of ghost particles
        self.realList = parameters['realList']             #list of real particles
        self.ghostList = parameters['ghostList']           #list of ghost particles
        self.realName = parameters['realName']             #real name
        self.intermediateName = parameters['intermediateName']#intermediate name
        self.ghostName = parameters['ghostName']           #ghost name     
        self.simulation = parameters['simulation']         #simulation object        
        self.nonbondedforce = parameters['nonbondedforce'] #nonbonded force object
        self.customnonbondedforce = parameters['customnonbondedforce']#custom nonbonded force object
        self.mdStep = parameters['mdStep']                 #MD steps between MC steps
        self.upperLimit = parameters['upperLimit']         #upper limit for cluster size
        self.lowerLimit = parameters['lowerLimit']         #lower linit for cluster size
               
        #nonbonded interaction parameters
        self.sigma = parameters['sigma']                      
        self.epsilon = parameters['epsilon']     
        self.charge = parameters['charge']         
        
        #initiate other variables
        self.randIndex = 0                   #index of random particle to insert/delete
        self.refIndex = 0                    #index of reference atom
        self.numInsr = 0                     #num of insertion step
        self.numDel = 0                      #num of deletion step
        self.sucInsr = 0                     #num of successful insertion
        self.sucDel = 0                      #num of successful deletion
        self.m = 0                           #label of intermediate state, m=0 
                                             #means interaction is 0, m = M means the 
                                             #interaction is fully on
        self.lamda = []                      #scaling parameters for intermediate states
        #bias/coupling parameter for intermediate states
        self.eta = [0.0 for x in range(self.upperLimit*self.M+1)]
        #histogram for intermediate states
        self.histogram = [0.0 for x in range(self.upperLimit*self.M+1)]                                
        self.sumNumReal = 0                  #sum of of the num of particles at m == 0 state
        self.numSum = 0                      #num of m == 0 state
        
        #freeEnergy 
        self.freeEnergy = [0.0 for x in range(self.upperLimit*self.M+1)]                 
        self.energy = [0.0 * kilojoule/mole for x in range(self.upperLimit*self.M+1)]            
        self.cutoff = 1.5*self.sigma
        self.f = 1.0                         #wang-landau initial modification factor

        self.nIn = 0
        self.vIn = 4.0/3.0*pi*pow(self.cutoff,3)*nanometers**3
        for i in range(self.M):
            self.lamda.append(i*i*1.0/(self.M*self.M))
        self.lamda.append(1.0)
    
    def updateForce(self,scaling):
        """Funtion to scale the interaction of the selected atom
        
        Parameters
        ---------
        scaling:
            scaling parameter

        """
        self.nonbondedforce.setParticleParameters(self.randIndex,
            self.charge*scaling,self.sigma,self.epsilon*scaling**2)
        
        if scaling == 0.0:
            self.customnonbondedforce.setParticleParameters(self.refIndex,[0.0])
            self.customnonbondedforce.setParticleParameters(self.randIndex,[0.0])
        else:
            self.customnonbondedforce.setParticleParameters(self.refIndex,[1-scaling])
            self.customnonbondedforce.setParticleParameters(self.randIndex,[1-scaling])

        self.nonbondedforce.updateParametersInContext(self.simulation.context)
        self.customnonbondedforce.updateParametersInContext(self.simulation.context)

    def DFS(self,graph,start):
        """

        """
        visited = set()
        stack = [start]
        while stack:
            vertex = stack.pop()
            if vertex not in visited:
                visited.add(vertex)
                stack.extend(graph[vertex] - visited)
        return visited 

    def generateGraph(self,positions):
        """

        """
        bond = {}

        for atom in self.realList:
            bond[atom] = set()

        for atomA in self.realList:
            for atomB in self.realList:
                if atomA < atomB:
                    distance = self.getDistance(positions[atomA],positions[atomB])
                    if distance < self.cutoff:
                        bond[atomA].add(atomB)
                        bond[atomB].add(atomA)
        return bond 
    
    def checkClusterCriteria(self,positions):
        """

        """
        #set up the graph
        bond = self.generateGraph(positions)

        #using DFS-like method to find one cluster    
        size = len(self.DFS(bond,self.realList[0]))
        if size == len(self.realList):
            return True
        else:
            return False

    def findClusters(self,positions):
        """

        """
        #set up the graph 
        clusters = []
        num = 0
        bond = self.generateGraph(positions)

        realList = list(self.realList)
        while (num != self.numReal):
            cluster = self.DFS(bond,realList[0])
            num += len(cluster)  
            clusters.append(cluster)
            for atom in cluster:
                realList.remove(atom)

        if len(clusters) > 2:
            print ("that's messy, now we have to deal with 3 clusters")

        return clusters

    def setRandomPosition(self):
        """Function to set random position for new inserted atoms

        """
        positions = self.simulation.context.getState(getPositions
            =True).getPositions()

        self.refIndex = self.realList[randrange(self.numReal)]

        #calculate coordination number of the reference atom
        nIn = 0
        for atom in self.realList:
            if atom != self.refIndex:
                distance = self.getDistance(positions[self.refIndex],positions[atom])
                if distance < self.cutoff:
                    nIn += 1       

        vec = random_vector_inside_sphere()
        displace = Vec3(vec[0],vec[1],vec[2])*self.cutoff*nanometers
        positions[self.randIndex] = positions[self.refIndex]+displace
        
        self.simulation.context.setPositions(positions)
        return nIn 

    def getRandomAtoms(self):
        """

        """
        positions = self.simulation.context.getState(getPositions
            =True).getPositions()

        goodAtoms = False
        nIn = 0
        while (not goodAtoms):
            nIn = 0
            goodAtoms = True
            randList = []
            self.refIndex = self.realList[randrange(self.numReal)]
            for atom in self.realList:
                if atom != self.refIndex:
                    distance = self.getDistance(positions[self.refIndex],positions[atom])
                    if distance < self.cutoff:
                        nIn += 1
                        randList.append(atom)
            self.randIndex = randList.pop(randrange(nIn))
            self.realList.remove(self.randIndex)

            #make sure after deletion, it's still a cluster
            if (not self.checkClusterCriteria(positions)):
                goodAtoms = False

            if (not goodAtoms):
                self.realList.append(self.randIndex)
                return -1 
        #while (not goodAtoms):
        #    goodAtoms = True
        #    refList = []
        #    self.randIndex = self.realList.pop(randrange(self.numReal))
        #    for atom in self.realList:
        #        distance = self.getDistance(positions[self.randIndex],positions[atom])
        #        if distance < self.cutoff:
        #            refList.append(atom)
        #    self.refIndex = refList.pop(randrange(len(refList)))

        #    #make sure after deletion, it's still a cluster
        #    if (not self.checkClusterCriteria(positions)):
        #        goodAtoms = False

        #    if (not goodAtoms):
        #        self.realList.append(self.randIndex)

        return nIn 

    def setRandomVelocity(self):
        """Function to ser random velocity for new inserted atoms

        """
        velocities = self.simulation.context.getState(getVelocities
            =True).getVelocities()
        #self.simulation.context.setVelocitiesToTemperature(
        #    self.temp/kelvin,randint(1, 10000))
        self.simulation.context.setVelocitiesToTemperature(
            self.temp/kelvin,randint(1, 10000))
        new_velocities = self.simulation.context.getState(getVelocities
            =True).getVelocities()
        velocities[self.randIndex] = new_velocities[self.randIndex]
        self.simulation.context.setVelocities(velocities)


    def changeNameReal(self):
        """change the name of the atoms to real name

        """
        self.atoms[self.randIndex].name = self.realName
 
    def changeNameGhost(self):
        """change the name of the atoms to ghost name

        """
        self.atoms[self.randIndex].name = self.ghostName

    def changeNameIntermediate(self):
        """chaneg the name of the atoms to intermediate name

        """
        self.atoms[self.randIndex].name = self.intermediateName
        
    def insertion(self):
        """Function to insert

        """
        if self.numReal == self.upperLimit:
            pass
        else:
            self.numInsr += 1
            oldPotential = self.simulation.context.getState(getEnergy=True).getPotentialEnergy()
            #need to choose new atom to insert if m==0, otherwise just turn on more interaction
            if self.m == 0:
                #set random positions if insertion starting from noexisted ghost particle
                self.randIndex = self.ghostList.pop(randrange(self.numGhost))
                self.nIn = self.setRandomPosition()

            #update force
            self.updateForce(self.lamda[self.m+1])
        
            newPotential = self.simulation.context.getState(getEnergy=True).getPotentialEnergy()
            #calculate the acceptance ratio
            dU = newPotential - oldPotential
            acc = min(1, exp(-self.beta*(dU - (self.lamda[self.m+1]-self.lamda[self.m])*self.mu) +
                  self.eta[(self.numReal*self.M+self.m)+1] - self.eta[self.numReal*self.M+self.m]) 
                  * pow(self.vIn*self.numReal/self.wavelengthCube,self.lamda[self.m+1] - self.lamda[self.m])
                  * pow((self.numReal+self.lamda[self.m])*(self.nIn+self.lamda[self.m]),self.lamda[self.m])/
                  pow((self.numReal+self.lamda[self.m+1])*(self.nIn+self.lamda[self.m+1]),self.lamda[self.m+1]))             
            #acc = min(1, exp(-self.beta*(dU - (self.lamda[self.m+1]-self.lamda[self.m])*self.mu) + 
            #      self.eta[(self.numReal*self.M+self.m)+1] - self.eta[self.numReal*self.M+self.m]) * pow(self.vol/self.wavelengthCube,
            #      self.lamda[self.m+1] - self.lamda[self.m]) * pow(self.numReal+self.lamda[self.m],
            #      self.lamda[self.m])/pow(self.numReal + self.lamda[self.m+1],self.lamda[self.m+1]))        
            #decide whether to accept or reject
            if random() < acc:
                self.sucInsr += 1
                if self.m == 0:
                    self.m = 1
                    #set random velocity for the new particle
                    self.setRandomVelocity()
                    self.changeNameIntermediate()
                elif self.m == self.M-1:
                    self.m = 0
                    self.realList.append(self.randIndex)
                    self.changeNameReal()
                    self.numReal += 1
                    self.numGhost -= 1
                else:
                    self.m += 1
            else:
                self.updateForce(self.lamda[self.m])
                if self.m == 0:
                    self.ghostList.append(self.randIndex)

    def deletion(self):
        """Function to delete

        """
        #if there is no atoms, just pass
        if self.numReal == self.lowerLimit and self.m == 0:
            pass
        else:
            self.numDel += 1
            oldPotential = self.simulation.context.getState(getEnergy=True).getPotentialEnergy()
            #need to seletct atom to delete if m==0
            if self.m == 0:
                self.nIn = self.getRandomAtoms()
                if (self.nIn == -1):
                    return 0 
                #update force
                self.updateForce(self.lamda[self.M-1])
            else:
                self.updateForce(self.lamda[self.m-1]) 
            newPotential = self.simulation.context.getState(getEnergy=True).getPotentialEnergy() 
            #calculate the acceptance ratio
            dU = newPotential - oldPotential
            if self.m == 0:
                acc = min(1, exp(-self.beta*(dU - (self.lamda[self.M-1]-1.0)*self.mu)+ self.eta[(self.numReal*self.M+self.m)-1] -
                      self.eta[self.numReal*self.M+self.m]) 
                      * pow((self.numReal-1)*self.vIn/self.wavelengthCube,self.lamda[self.M-1]-1.0)
                      * pow(self.numReal*self.nIn,1.0)/
                      pow((self.numReal-1.0+self.lamda[self.M-1])*(self.nIn-1.0+self.lamda[self.M-1]),self.lamda[self.M-1]))
                #acc = min(1, exp(-self.beta*(dU - (self.lamda[self.M-1]-1.0)*self.mu)+ self.eta[(self.numReal*self.M+self.m)-1] - 
                #      self.eta[self.numReal*self.M+self.m]) * pow(self.vol/self.wavelengthCube,self.lamda[self.M-1]-1.0) * 
                #      pow(self.numReal,1.0)/pow(self.numReal-1.0+self.lamda[self.M-1],self.lamda[self.M-1]))
            else:
                acc = min(1, exp(-self.beta*(dU - (self.lamda[self.m-1]-self.lamda[self.m])*self.mu) + 
                      self.eta[(self.numReal*self.M+self.m)-1] - self.eta[self.numReal*self.M+self.m]) 
                      * pow(self.numReal*self.vIn/self.wavelengthCube,self.lamda[self.m-1]-self.lamda[self.m]) 
                      * pow((self.numReal+self.lamda[self.m])*(self.nIn+self.lamda[self.m]),self.lamda[self.m])/
                      pow((self.numReal+self.lamda[self.m-1])*(self.nIn+self.lamda[self.m-1]),self.lamda[self.m-1]))  
                #acc = min(1, exp(-self.beta*(dU - (self.lamda[self.m-1]-self.lamda[self.m])*self.mu) + 
                #      self.eta[(self.numReal*self.M+self.m)-1] - self.eta[self.numReal*self.M+self.m]) * pow(self.vol/self.wavelengthCube,
                #      self.lamda[self.m-1] - self.lamda[self.m]) * pow(self.numReal+self.lamda[self.m],
                #      self.lamda[self.m])/pow(self.numReal + self.lamda[self.m-1],self.lamda[self.m-1]))
            #decide whether to accept or reject
            if random() < acc:
                self.sucDel += 1
                if self.m == 0:
                    self.m = self.M-1
                    self.numReal -= 1
                    self.nIn -= 1
                    self.numGhost += 1
                    self.changeNameIntermediate()
                elif self.m == 1:
                    self.m = 0 
                    self.ghostList.append(self.randIndex)
                    self.changeNameGhost()
                else: 
                    self.m -= 1
            else: 
                if self.m == 0: 
                    self.updateForce(1.0)
                    self.realList.append(self.randIndex)
                else: 
                    self.updateForce(self.lamda[self.m])
        
    def step(self):
        """run a gcmc step, insert/delete

        """
        if random() >= 0.5:
            self.insertion()
        else: 
            self.deletion()
        #calculate statiscals
        if self.m == 0:
            self.sumNumReal += self.numReal
            self.numSum += 1

        self.histogram[self.numReal*self.M+self.m] += 1.0
        self.energy[self.numReal*self.M+self.m] += self.simulation.context.getState(getEnergy=True).getPotentialEnergy() 
        self.eta[self.numReal*self.M+self.m] -= self.f 
        #self.simulation.step(self.mdStep)      
        
        #After md, if it's a real state, use DFS to check cluster critria,
        #if the cluster falls apart, add a constraint to connect closest atoms
        #in two clusters and running 100/1000 md steps, keep checking cluster 
        #critria and add constraint until the it satisfies the cluster critria
        #if self.m == 0:  
        #    positions = self.simulation.context.getState(getPositions
        #        =True).getPositions()
        #    while (not self.checkClusterCriteria(positions)):
        #        print("fall apart, I can fix it")
        #        clusters = self.findClusters(positions)
        #        distance = 20.0
        #        indexA = 0
        #        indexB = 0
        #        for atomA in clusters[0]:
        #            for atomB in clusters[1]:
        #                newDistance = self.getDistance(positions[atomA],positions[atomB])
        #                if newDistance < distance:
        #                    distance = newDistance
        #                    indexA = atomA
        #                    indexB = atomB
        #        print (indexA)
        #        print (indexB)
        #        print (distance)
        #        self.customnonbondedforce.setParticleParameters(indexA,[0.1])
        #        self.customnonbondedforce.setParticleParameters(indexB,[0.1])
        #        self.customnonbondedforce.updateParametersInContext(self.simulation.context)
        #        self.simulation.step(100)
        #        self.customnonbondedforce.setParticleParameters(indexA,[0.0])
        #        self.customnonbondedforce.setParticleParameters(indexB,[0.0])
        #        self.customnonbondedforce.updateParametersInContext(self.simulation.context)
        #        positions = self.simulation.context.getState(getPositions
        #            =True).getPositions()
   
    def getInsrPrabability(self):
        """get insertion probability

        """
        return self.sucInsr*1.0/self.numInsr
        
    def getDelPrabability(self):
        """get deletion probability

        """
        return self.sucDel*1.0/self.numDel
   
    def getAverageNum(self):
        """get average num of particles, only count m=0 state,
        state of half interaction is not included into calculation

        """
        return self.sumNumReal*1.0/self.numSum
         
     
    def rezeroStat(self):
        """after equilibrium, to rezero the some statics

        """
        self.numInsr = 0
        self.numDel = 0
        self.sucInsr = 0
        self.sucDel = 0
        self.sumNumReal = 0
        self.numSum = 0
        self.rezeroHistogram() 
        
             
    def rezeroHistogram(self):
        """rezero histogram

        """
        for i in range(self.upperLimit*self.M+1):
            self.histogram[i] = 0.0
            self.energy[i] = 0.0 * kilojoule/mole

    def updateWangLandauFactor(self):
        """wang-landu sampling to optimize the bias, you have to call 
        this funtion to optimize in main program

        """
        sumHistogram = 0.0
        for i in range(self.lowerLimit*self.M,self.upperLimit*self.M+1):
            sumHistogram += self.histogram[i]
        aveHistogram = sumHistogram*1.0/((self.upperLimit-self.lowerLimit)*self.M+1)
        count = 0
        for i in range(self.lowerLimit*self.M,self.upperLimit*self.M+1):
            if (self.histogram[i] > 0.8*aveHistogram):
                count += 1
        if (count == ((self.upperLimit-self.lowerLimit)*self.M+1)):
            self.f = 0.5*self.f
            self.rezeroHistogram()

    def getEnergy(self):   
        """get energy for each state

        """
        energy = [0.0 * kilojoule/mole for x in range(self.upperLimit*self.M+1)]                 
        for i in range(self.lowerLimit*self.M,self.upperLimit*self.M+1):
            if self.histogram[i] > 0.0:
                energy[i] = self.energy[i]/self.histogram[i] 
        return energy 

    def getFreeEnergy(self):
        """get free energy for each state

        """ 
        for i in range(self.lowerLimit*self.M,self.upperLimit*self.M+1):
            self.freeEnergy[i] = -(log(self.histogram[i]*1.0/self.histogram[self.lowerLimit*self.M])-self.eta[i]+self.eta[self.lowerLimit*self.M])/self.beta
        return self.freeEnergy 

    def getDistance(self,positionA,positionB):
        """

        """        
        distanceSquare = 0.0
        distance = 0.0
        distanceVector =  (positionA-positionB)/nanometers
        for i in range(3):
            distanceSquare += (distanceVector[i]-self.box[i]*round(distanceVector[i]/self.box[i]))**2
        distance = sqrt(distanceSquare)
        return distance
