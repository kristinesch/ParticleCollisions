import numpy as np
from simulation import *

def getVelocities(Parr):
    velocities=[]
    for p in Parr:
        velocities.append(np.sqrt(p[VX]*p[VX]+p[VY]*p[VY]))
    return velocities


def histogram(velocities):
    plt.hist(velocities,bins=101,density=True)


#NB 2 degrees of freedom!!!!!!!!!!!!!!!!!!
def HistogramMB(velocities,m,filename): #filename is a string containing the name to save the file as
    v=np.linspace(0,4)
    averageVelocity=np.sum(velocities)/len(velocities)
    #kT=(averageVelocity*averageVelocity*m/2) #Boltzmann constant times Temperature
    kT=V0*V0*m/2
    mb=(m*v/kT)*np.exp(-(m*(v*v)/(2*kT))) #Maxwell-Boltzmann distribution

    #plotting
    fig,ax=plt.subplots(1,1)
    fig.suptitle("Speed distribution")
    ax.set_xlabel("Speed")
    ax.set_ylabel("Number of particles")
    ax.hist(velocities,bins=51,density=True)
    ax.plot(v,mb)
    fig.savefig(filename)
    plt.show()
    return averageVelocity

#returns average kinetic energy and velocity for the particles in the input array
def averageKineticEnergyAndVelocity(parr): #assumes same mass for whole array
    velocities=getVelocities(parr)
    averageVelocity=np.sum(velocities)/len(velocities)
    kineticEnergy=averageVelocity*averageVelocity*0.5*parr[1,M]
    return kineticEnergy, averageVelocity

def averageKineticEnergy(parr):
    velocities=getVelocities(parr)
    averageVelocity=np.sum(velocities)/len(velocities)
    kineticEnergy=averageVelocity*averageVelocity*0.5*parr[1,M]
    return kineticEnergy

def avrgEnergyP(p): #if p is only one particle
    velocitySquared=p[VX]*p[VX]+p[VY]*p[VY]
    kineticEnergy=velocitySquared*0.5*p[M]
    return kineticEnergy