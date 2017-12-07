from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
from math import exp
from random import seed
import numpy as np
from gcmc_restraint_check import Gcmc

print ("GCMC for LiF ion pair")
#variable definition
R = 8.3144598*joule/(kelvin*mole)
temp = 80.4186*kelvin
beta = 1.0/(R*temp)
mu = -7.7964*kilojoule/mole          #chemical potential 
wavelength = 0.030802*nanometer         #thermal de broglie wavelength
wavelengthCube = wavelength**3
numReal = 1                               #num of real ion pairs in initial state
M = 2                                    #number of intermediate states + 1
charge = 0.0
sigma = 0.3405
epsilon = 0.997967
realName = 'H'
ghostName = 'Hg'
intermediateName = 'Hi'
upperLimit = 30                            #maximum num of ion pairs
lowerLimit = 1                            #minimum num of ion pairs 

#set down the initial state and simulation parameters
platform = Platform.getPlatformByName('CUDA')
properties = {'CudaPrecision':'mixed'}
pdb = PDBFile('argon_gcmc.pdb')
forcefield = ForceField('argon.xml','spce.xml')
atoms = list(pdb.topology.atoms())
system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffPeriodic,
         nonbondedCutoff=1.53225*nanometer, constraints = AllBonds, rigidWater = True)

#create custom nonbonded force
restraint = CustomNonbondedForce("A*b1*b2*max(0,r-rmax)^2") 
restraint.addGlobalParameter("A", 1000000.0*kilojoules_per_mole/nanometers**2)
restraint.addGlobalParameter("rmax", sigma*1.5*nanometers)
restraint.addPerParticleParameter("b")
setH = set()
for atom in atoms:
    para = [0.0]
    restraint.addParticle(para)
    setH.add(atom.index)
restraint.addInteractionGroup(setH,setH)
restraint.setCutoffDistance(2.0)
restraint.setNonbondedMethod(2)
system.addForce(restraint)

#get nonbondedforce object from system to switch some real particles to ghost particles
forces = system.getForces()
for force in forces:
    if type(force) is NonbondedForce:
        nonbondedforce = force
    if type(force) is CustomNonbondedForce:
        custombondedforce = force
nonbondedforce.setReactionFieldDielectric(0.0)
nonbondedforce.setUseDispersionCorrection(False)

#create two lists for real molecules and ghost molecules and switch some molecules 
#to ghost particles
realList = []
ghostList = []
for i in range(numReal):
    realList.append(i)
numGhost = 30-numReal
for i in range(numReal,30):
    ghostList.append(i)
    nonbondedforce.setParticleParameters(i,0.0,0.3405,0.0)

#generate simulation object and equilibrate
integrator = LangevinIntegrator(temp, 1/picosecond, 0.002*picoseconds)
simulation = Simulation(pdb.topology, system, integrator,platform)
simulation.context.setPositions(pdb.positions)
simulation.minimizeEnergy()
#simulation.reporters.append(PDBReporter('gcmc_out_%d-%d.pdb' %(lowerLimit,upperLimit), 100000000000000))
#simulation.reporters.append(StateDataReporter(stdout, 1, step=True,
#                            potentialEnergy=True, temperature=True))
simulation.step(10000)

box = simulation.topology.getUnitCellDimensions().value_in_unit(nanometers)
vol = box[0]*box[1]*box[2]*nanometer**3


#number of gcmc steps
numGcmc = 2000000
count = 0
#set random generator
seed()

#create gcmc object
parameters = {}
parameters['temp'] = temp
parameters['beta'] = beta
parameters['mu'] = mu
parameters['vol'] = vol
parameters['box'] = box
parameters['atoms'] = atoms
parameters['wavelengthCube'] = wavelengthCube
parameters['M'] = M
parameters['numReal'] = numReal
parameters['numGhost'] = numGhost
parameters['realList'] = realList
parameters['ghostList'] = ghostList
parameters['realName'] = realName
parameters['ghostName'] = ghostName
parameters['intermediateName'] = intermediateName
parameters['simulation'] = simulation
parameters['nonbondedforce'] = nonbondedforce
parameters['customnonbondedforce'] = restraint
parameters['mdStep'] = 10
parameters['sigma'] = sigma
parameters['charge'] = charge
parameters['epsilon'] = epsilon
parameters['upperLimit'] = upperLimit
parameters['lowerLimit'] = lowerLimit
gcmc = Gcmc(parameters)


#simulation.reporters.append(StateDataReporter(stdout, 1, step=True,
#                            potentialEnergy=True, temperature=True))
#simulation.reporters.append(DCDReporter('output_%d-%d.dcd' %(lowerLimit,upperLimit), 100))
#equilibration
#j = 0;
while (exp(gcmc.f) > 1.001):
    for i in range(3000):
        #print (j*1000+i)
        #state = simulation.context.getState(getPositions=True)
        #simulation.reporters[0].report(simulation, state)
        gcmc.step()
        print ("numReal %f numGhost %f" %(gcmc.numReal+gcmc.lamda[gcmc.m],gcmc.numGhost-gcmc.lamda[gcmc.m])) 
    print("#histogram")
    print(gcmc.histogram)
    print("eta")
    print(gcmc.eta)
    print("energy")
    print(gcmc.getEnergy())
    #every 1000 steps, update wanglaudafactor
    gcmc.updateWangLandauFactor()
    print("wang-landau factor")
    print(gcmc.getInsrPrabability())
    print(gcmc.getDelPrabability())
    print(exp(gcmc.f))
    #j += 1;

gcmc.rezeroStat()

#production
gcmc.f = 0.0
for i in range(numGcmc):
    gcmc.step()
    print ("numReal %f numGhost %f" %(gcmc.numReal+gcmc.lamda[gcmc.m],gcmc.numGhost-gcmc.lamda[gcmc.m])) 
print(gcmc.getAverageNum())
print(gcmc.getInsrPrabability())
print(gcmc.getDelPrabability())
print(gcmc.eta)
print(gcmc.histogram)
x = gcmc.getFreeEnergy()
print(x)
for i in range(61):
    print (x[i])
print (" ")
x = gcmc.getEnergy()
print(x)
for i in range(61):
    print (x[i])

