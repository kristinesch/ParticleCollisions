from simulation import *  
from functions import *
from scipy import stats

#average collision number=number of collisions/number of particles




"""PROBLEM 1"""
#N=number of particles, n=number of iterations,m0=mass, interval=number of collisions between data sampling
def problem1(N1,m0,n,interval,start,rc,r0,v0):
    rc=1
    velocities=[]
    velCorr=[]
    #initialize mass array
    masses1=np.full(N1,m0)

    #initialize particle array
    parr1=initParticleArray(N1,masses1,r0,m0,v0,rc)
    initVel=getVelocities(parr1)
    histogram(initVel,"initHistp1")

    #run simulation
    parrData=simulationData(parr1,n,interval,start,rc)

    #Get velocity data, and plot histogram with Maxwell-Boltzmann distribution
    for parr in parrData: #add velocity data to list
        velCorr.append(getVelocities(parr))
        velocities.extend(getVelocities(parr))

    #plotting histogram and MB 
    HistogramMB(velocities,m0,"problem1")
    
    #checking correlation
    print("Pearson correlation coefficient: ",stats.pearsonr(velCorr[0],velCorr[1]))
    fig,ax=plt.subplots(1,1)
    ax.scatter(velCorr[0],velCorr[1],"bo")
    ax.set_xlabel("velocity data sample 1")
    ax.set_ylabel("velocity data sample 2")
    plt.suptitle("Correlation for "+str(N1)+" particles and "+str(interval)+" collisions")
    fig.savefig("correlationPlot")
    plt.show()

    
        


"""PROBLEM 2"""

def problem2(N2,m0,n,r0,v0,rc):
    N2=1000
    m0=1
    n=10000

    #initialize mass array
    masses2=np.full(N2,m0)
    x=int(N2/2)
    for i in range(x):
        masses2[x+i]=4*m0

    #Initialize particle array
    particleArray2=initParticleArray(N2,masses2,r0,m0,v0,rc)

    #initial plot
    initialVelocities=getVelocities(particleArray2)
    histogram(initialVelocities,"initHistp2")

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



    

def problem3(N3,m0,n,interval,start,r0,v0,rc):

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
    particleArray3=initParticleArray(N3,masses3,r0,m0,v0,rc)
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
    ax1.plot(timesteps,kinE4m0,label="Kinetic Energy of particles with mass 4m0")
    ax1.plot(timesteps,kinEm0,label="Kinetic Energy of particles with mass m0")
    ax1.legend()
    fig1.savefig("KineticEneryPlot")
    plt.show()

    #plotting velocities:
    fig2,ax2=plt.subplots(1,1)
    fig2.suptitle("Average velocity")
    ax2.set_xlabel("Number of collisions")
    ax2.set_ylabel("Velocity")
    ax2.plot(timesteps,velm0,label="Velocity of particles with mass m0")
    ax2.plot(timesteps,vel4m0,label="Velocity of particles with mass 4m0")
    ax2.legend()
    fig2.savefig("VelocityPlot")
    plt.show()



N1=3000
m0=1
n1=30000
interval1=10000
start1=10000
rc=1
r0=0.001 #must be large enough so that particle particle collisions happen often enough???
v0=1


problem1(N1,m0,n1,interval1,start1,rc,r0,v0)
#problem2(N2,m0,n,r0,v0,rc)
#problem3(N3,m0,n,interval,start,r0,v0,rc)


#problem1(2000,1,20000,5000,10000)
#problem2(500,1,5000)

#problem3(500,1,5000,100,0)