import numpy as np
from simulation import *


#get the absoute values of the velocity for each particle in Parr
def getVelocities(Parr):
    velocities=[]
    for p in Parr:
        velocities.append(np.sqrt(p[VX]*p[VX]+p[VY]*p[VY]))
    return velocities

#plotting histogram of velocities
def histogram(velocities,filename):#filename is a string containing the name to save the file as
    fig,ax=plt.subplots(1,1)
    fig.suptitle("Speed distribution")
    ax.set_xlabel("Speed")
    ax.set_ylabel("Number of particles")
    ax.hist(velocities,bins=41, range=(0,2))
    fig.savefig(filename)
    plt.show()


#NB 2 degrees 
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
    ax.hist(velocities,bins=91,density=True)
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

"""
#Delete?
def averageKineticEnergy(parr):
    velocities=getVelocities(parr)
    averageVelocity=np.sum(velocities)/len(velocities)
    kineticEnergy=averageVelocity*averageVelocity*0.5*parr[1,M]
    return kineticEnergy
"""


def totKinEnergy(parr):
    kineticEnergy=(parr[:,VX]*parr[:,VX]+parr[:,VY]*parr[:,VY])*parr[:,M]*0.5
    return np.sum(kineticEnergy)
    

def craterSize(initialParr,Parr):
    stillParticles=0 #for counting
    N=len(Parr) #number of particles
    for i in range(N): #counting how many particles havent moved
        if (initialParr[i,X]==Parr[i,X] and initialParr[i,Y]==Parr[i,Y]):
            stillParticles+=1
    return (N-stillParticles)/N #percentage of particles that are part of the crater


def plotParameterVsCraterSize(filename,parameter,figurename): #plots data from file "filename.txt". parameter=name of parameter(string) #figurename=name of file figure is saved to
    #loading data
    data=np.load(filename) #filename is a saved 2d numpy array

    #plotting
    fig, ax=plt.subplots(1,1)
    ax.plot(data[0],data[1])
    ax.suptitle("Crater size vs "+parameter)
    ax.set_xlabel("Crater size")
    ax.set_ylabel(parameter)
    fig.savefig(figurename)
    plt.show()


