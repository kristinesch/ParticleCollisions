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


def histogram2(velocitiesm0,velocities4m0,filename):
    fig,ax=plt.subplots(1,2)
    fig.suptitle("Speed distribution")
    ax[0].set_xlabel("Speed")
    ax[0].set_ylabel("Number of particles")
    ax[1].set_xlabel("Speed")
    ax[1].set_ylabel("Number of particles")
    ax[0].hist(velocitiesm0,bins=41,range=(0,4))
    ax[0].set_title("mass m0")
    ax[1].hist(velocities4m0,bins=41,range=(0,4))
    ax[1].set_title("mass 4m0")
    plt.tight_layout()
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
    ax.set_ylabel("Probability density")
    ax.hist(velocities,bins=41,density=True)
    ax.plot(v,mb)
    fig.savefig(filename)
    plt.show()
    return averageVelocity


#returns average speed and kinetic energy of the particles in parr
def avrgKinEnergyAndVelocity(parr):
    velocity=parr[:,VX]*parr[:,VX]+parr[:,VY]*parr[:,VY]
    avrgVel=np.sum(velocity)/len(parr)
    kineticEnergy=(parr[:,VX]*parr[:,VX]+parr[:,VY]*parr[:,VY])*parr[:,M]*0.5
    avrgEk=np.sum(kineticEnergy)/len(parr)
    return avrgEk,avrgVel

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


