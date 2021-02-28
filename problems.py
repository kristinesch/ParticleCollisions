from simulation import *  

#average collision number=number of collisions/number of particles

#FUNKER AV EN ELLER ANNEN GRUNN BARE I TERMINAL, WHAT THE ACTUAL FUCK
def getVelocities(Parr):
    velocities=[]
    for p in Parr:
        velocities.append(np.sqrt(p[VX]*p[VX]+p[VY]*p[VY]))
    return velocities

def histogram(velocities):
    plt.hist(velocities,bins=101,density=True)

    #kin=3/2kT
    #T=(kin/k*2/3)

#NB 2 degrees of freedom!!!!!!!!!!!!!!!!!!
def HistogramMB(velocities,m):
    v=np.linspace(0,2)
    averageVelocity=np.sum(velocities)/len(velocities)
    kT=averageVelocity*averageVelocity*m/2 #Boltzmann constant times Temperature
    mb=(m*v/kT)*np.exp(-(m*(v*v)/(2*kT))) #Maxwell-Boltzmann distributuin
    plt.hist(velocities,bins=41,density=True)
    plt.plot(v,mb)
    plt.show()
    return averageVelocity

#returns average kinetic energy and velocity for the particles in the input array
def averageKineticEnergyAndVelocity(parr): #assumes same mass for whole array
    velocities=getVelocities(parr)
    averageVelocity=np.sum(velocities)/len(velocities)
    kineticEnergy=averageVelocity*0.5*parr[1,M]
    return kineticEnergy, averageVelocity



# TestArray=simulation(100,1000)
# TestV=getVelocities(TestArray)
# HistogramMB(TestV,1)

def problem1():
    N1=1000
    m0=1

    masses1=np.zeros(N1)
    for i in range(N1):
        masses1[i]=m0

    parr1=initParticleArray(1000,masses1)
    parrData=simulationData(parr1,10000,2000)
    for parr in parrData:
        vel=getVelocities(parr)
        #print("VELOCITIES\n",vel)
        HistogramMB(vel,m0)
        


#TestArray=simulation(100,1000)
#TestV=getVelocities(TestArray)



"""PROBLEM 2"""

def problem2():
    N2=1000
    m0=1

    #initialize mass array
    masses2=np.zeros(N2)
    x=int(N2/2)
    for i in range(x):
        masses2[i]=m0
    for i in range(x):
        masses2[x+i]=4*m0

    print(masses2)
    particleArray2=initParticleArray(N2,masses2)

    #initial plot
    initialVelocities=getVelocities(particleArray2)
    histogram(initialVelocities)

    finalParr=simulation(particleArray2,N2,10000)
    Parrm0=finalParr[:x]
    velm0=getVelocities(Parrm0)

    Parr4m0=finalParr[x:]
    vel4m0=getVelocities(Parr4m0)

    averageVelocitym0=HistogramMB(velm0,m0)
    averageVelocity4m0=HistogramMB(vel4m0,4*m0)
    kineticEnergym0=0.5*m0*averageVelocitym0
    kineticEnergy4m0=0.5*4*m0*averageVelocity4m0
    print("Average velocity for particles with mass m0: ",averageVelocitym0)
    print("Average kinetic energy for particles with mass m0",kineticEnergym0)
    print("Average velocity for particles with mass 4m0: ",averageVelocity4m0)
    print("Average kinetic energy for particles with mass 4m0",kineticEnergy4m0)



    

def problem3():
    N3=1000
    m0=1
    n=15000 #number of total collisions
    interval=100 #intervals between data sampling

    #initialize mass array
    masses3=np.zeros(N3)
    x=int(N3/2)
    for i in range(x):
        masses3[i]=m0
    for i in range(x):
        masses3[x+i]=4*m0

    #lists for storing data 
    kinEm0=[]
    kinE4m0=[]
    velm0=[]
    vel4m0=[]

    particleArray3=initParticleArray(N3,masses3)

    data=simulationData(particleArray3,n,interval)

    for parr in data:
        Parrm0=parr[:x]
        Parr4m0=parr[x:]
        Ekm0,vm0=averageKineticEnergyAndVelocity(Parrm0)
        Ek4m0,v4m0=averageKineticEnergyAndVelocity(Parr4m0)
        kinE4m0.append(Ek4m0)
        kinEm0.append(Ekm0)
        velm0.append(vm0)
        vel4m0.append(v4m0)
    
    timesteps=np.arange(0,n+interval,interval)
    fig,ax=plt.subplots(1,1)
    ax.plot(timesteps,kinE4m0,label="Kinetic Energy 4m0")
    ax.plot(timesteps,kinEm0,label="Kinetic Energy m0")
    ax.plot(timesteps,velm0,label="Velocity m0")
    ax.plot(timesteps,vel4m0,label="Velocity 4m0")
    ax.legend()
    fig.savefig("KineticEnergyAndVelocity")
    plt.show()






problem3()