from simulation import *  
from functions import *

#average collision number=number of collisions/number of particles

RC=1


"""PROBLEM 1"""
#N=number of particles, n=number of iterations,m0=mass, interval=number of collisions between data sampling
def problem1(N1,m0,n,interval,start):
    velocities=[]
    #initialize mass array
    masses1=np.full(N1,m0)

    #initialize particle array
    parr1=initParticleArray(N1,masses1)

    #run simulation
    parrData=simulationData(parr1,n,interval,start,RC)

    #Get velocity data, and plot histogram with Maxwell-Boltzmann distribution
    for parr in parrData: #add velocity data to list
        velocities.extend(getVelocities(parr))

    HistogramMB(velocities,m0,"problem1")
        


"""PROBLEM 2"""

def problem2(N2,m0,n):
    N2=1000
    m0=1
    n=10000

    #initialize mass array
    masses2=np.full(N2,m0)
    x=int(N2/2)
    for i in range(x):
        masses2[x+i]=4*m0

    #Initialize particle array
    particleArray2=initParticleArray(N2,masses2)

    #initial plot
    initialVelocities=getVelocities(particleArray2)
    histogram(initialVelocities)

    #run simulation
    finalParr=simulation(particleArray2,n,RC)

    #Get arrays corresponding to the light and heavy particles
    Parrm0=finalParr[:x]
    Parr4m0=finalParr[x:]

    #Calculate average kinetic energies and velocities
    kineticEnergym0,averageVelocitym0=averageKineticEnergyAndVelocity(Parrm0)
    kineticEnergy4m0,averageVelocity4m0=averageKineticEnergyAndVelocity(Parr4m0)

    #printing results
    print("Average velocity for particles with mass m0: ",averageVelocitym0)
    print("Average kinetic energy for particles with mass m0",kineticEnergym0)
    print("Average velocity for particles with mass 4m0: ",averageVelocity4m0)
    print("Average kinetic energy for particles with mass 4m0",kineticEnergy4m0)



    

def problem3(N3,m0,n,interval,start):

    #initialize mass array
    masses3=np.full(N3,m0)
    x=int(N3/2)
    for i in range(x):
        masses3[x+i]=4*m0

    #lists for storing data 
    kinEm0=[]
    kinE4m0=[]
    velm0=[]
    vel4m0=[]

    #initialize array
    particleArray3=initParticleArray(N3,masses3)
    #run simulation and collect data
    data=simulationData(particleArray3,n,interval,start,RC)

    #calculate average velocity and kinetic energy and save to lists
    for parr in data:
        Parrm0=parr[:x]
        Parr4m0=parr[x:]
        Ekm0,vm0=averageKineticEnergyAndVelocity(Parrm0)
        Ek4m0,v4m0=averageKineticEnergyAndVelocity(Parr4m0)
        kinE4m0.append(Ek4m0)
        kinEm0.append(Ekm0)
        velm0.append(vm0)
        vel4m0.append(v4m0)

    #plotting kinetic energy:
    timesteps=np.arange(0,n+interval,interval)
    fig1,ax1=plt.subplots(1,1)
    fig1.suptitle("Average kinetic energy")
    ax1.set_xlabel("Number of collisions")
    ax1.set_ylabel("Kinetic energy")
    ax1.plot(timesteps,kinE4m0,label="Kinetic Energy 4m0")
    ax1.plot(timesteps,kinEm0,label="Kinetic Energy m0")
    plt.legend()
    fig1.savefig("KineticEneryPlot")
    plt.show()

    #plotting velocities:
    fig2,ax2=plt.subplots(1,1)
    fig2.suptitle("Average velocity")
    ax2.set_xlabel("Number of collisions")
    ax2.set_ylabel("Velocity")
    ax2.plot(timesteps,velm0,label="Velocity m0")
    ax2.plot(timesteps,vel4m0,label="Velocity 4m0")
    ax2.legend()
    fig2.savefig("VelocityPlot")
    plt.show()




#problem1(2000,1,20000,5000,10000)
problem2(500,1,5000)

problem3(500,1,5000,100,0)