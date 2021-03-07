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
    #correlation plot
    fig,ax=plt.subplots(1,1)
    ax.plot(velCorr[0],velCorr[1],"bo")
    ax.set_xlabel("velocity data sample 1")
    ax.set_ylabel("velocity data sample 2")
    plt.suptitle("Correlation for "+str(N1)+" particles and "+str(interval)+" collisions")
    fig.savefig("correlationPlot")
    plt.show()

    
        


"""PROBLEM 2"""

def problem2(N2,m0,n,r0,v0,rc):

    #initialize mass array
    masses2=np.full(N2,m0)
    x=int(N2/2)
    for i in range(x):
        masses2[x+i]=4*m0

    #Initialize particle array
    particleArray2=initParticleArray(N2,masses2,r0,m0,v0,rc)

    #initial plot
    initParrm0=particleArray2[:x]
    initParr4m0=particleArray2[x:]
    histogram2(getVelocities(initParrm0),getVelocities(initParr4m0),"p2start")
    

    #run simulation
    finalParr=simulation(particleArray2,n,rc)

    #Get arrays corresponding to the light and heavy particles
    Parrm0=finalParr[:x]
    Parr4m0=finalParr[x:]
    #final plot
    histogram2(getVelocities(initParrm0),getVelocities(initParr4m0),"p2end")
    #HistogramMB(getVelocities(initParrm0),m0,"p2m0MB")
    #HistogramMB(getVelocities(initParr4m0),4*m0,"p24m0MB")


    #Calculate average kinetic energies and velocities
    kineticEnergym0,averageVelocitym0=avrgKinEnergyAndVelocity(Parrm0)
    kineticEnergy4m0,averageVelocity4m0=avrgKinEnergyAndVelocity(Parr4m0)

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
    kinE=[]
    velAll=[]
    #initialize array
    particleArray3=initParticleArray(N3,masses3,r0,m0,v0,rc)

    #run simulation and collect data
    data=simulationData(particleArray3,n,interval,start,rc)


    #calculate average velocity and kinetic energy and save to lists
    for parr in data:
        Parrm0=parr[:x]
        Parr4m0=parr[x:]
        Ekm0,vm0=avrgKinEnergyAndVelocity(Parrm0)
        Ek4m0,v4m0=avrgKinEnergyAndVelocity(Parr4m0)
        Ek,vel=avrgKinEnergyAndVelocity(parr)
        kinE4m0.append(Ek4m0)
        kinEm0.append(Ekm0)
        velm0.append(vm0)
        vel4m0.append(v4m0)
        kinE.append(Ek)
        velAll.append(vel)

    print(kinE)
    
    timesteps=np.arange(0,n+interval,interval)
    fig,ax=plt.subplots(2,1,sharex=True)
    fig.suptitle("Average kinetic energy and velocity, rc="+str(rc))
    #plotting kinetic energy:
    ax[0].set_ylabel("Average kinetic energy")
    ax[0].plot(timesteps,kinE4m0,label="Particles with mass 4m0")
    ax[0].plot(timesteps,kinEm0,label="Particles with mass m0")
    ax[0].plot(timesteps,kinE,label="All particles")
    ax[0].legend()

    #plotting velocities:
    ax[1].set_xlabel("Number of collisions")
    ax[1].set_ylabel("Average velocity")
    ax[1].plot(timesteps,velm0,label="Particles with mass m0")
    ax[1].plot(timesteps,vel4m0,label="Particles with mass 4m0")
    ax[1].plot(timesteps,velAll,label="Particles with mass 4m0")
    ax[1].legend()
    fig.savefig("VelocityEkPlot"+str(rc)+".png")
    plt.show()



N1=3000
m0=1
n1=30000
interval1=10000
start1=10000
rc=1
r0=0.001 #must be large enough so that particle particle collisions happen often enough???
v0=1


#problem1(N1,m0,n1,interval1,start1,rc,r0,v0)

N2=3000
n2=10000
problem2(N2,m0,n2,r0,v0,rc)


N3=3000
n3=30000
interval3=300
start3=0
# rcList=[0.8,0.9,1]
# for rc3 in rcList:
#     problem3(N3,m0,n3,interval3,start3,r0,v0,rc3)




#problem1(2000,1,20000,5000,10000)
#problem2(500,1,5000)

#problem3(500,1,5000,100,0)