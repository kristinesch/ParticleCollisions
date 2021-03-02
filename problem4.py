from simulation import *
from functions import *

def initParticlesForProjectileImpact(N,m0,r0,m,r,v0): #m0 is the mass of the small particles, N=number of small particles
    #v0=downward speed of projectile particle
    print("Initializing particle array...")
    parr=np.zeros((N+1,7))

    #initializing ground particles:

    occupiedPositions=[] #list of occupied positions
    for i in range(N): #assign random values for positions
        Overlapping=True
        while(Overlapping==True): #checking if position is already occupied
                pos=np.array([random.uniform(r0,1-r0),random.uniform(r0,0.5-r0)]) #positions in the lower half of the box
                Overlapping=False
                for posi in occupiedPositions:
                    if (4*r0*r0>np.dot((posi-pos),(posi-pos))):
                        Overlapping=True
                        break

        parr[i,X]=pos[0]
        parr[i,Y]=pos[1]
        
        #add values for velocity, mass and radius
        parr[i,VX]=0
        parr[i,VY]=0
        parr[i,M]=m0
        parr[i,R]=r0

        occupiedPositions.append(pos)
        print((i*100)/N,"%")


    #initializing projectile particle:
    parr[N,X]=0.5
    parr[N,Y]=0.75
    parr[N,VX]=0
    parr[N,VY]=-v0
    parr[N,M]=m
    parr[N,R]=r

    print("Initialization done")
    return parr #returns initialized particle array

def plotCraterFormation(parr,r0,r):
    N=len(parr)
    print(N)
    groundParticles=parr[:N-1]
    #print(groundParticles)
    projectile=parr[N-1]
    plt.axis([Xmin,Xmax,Ymin,Ymax])
    plt.scatter(groundParticles[:,X],groundParticles[:,Y],c="b")
    plt.scatter(projectile[X],projectile[Y],c="r",s=r*r*500*500)
    plt.show()


def problem4Visualization(N,n,interval,start,m0,r0,m,r,v0,RC):
    parr=initParticlesForProjectileImpact(N,m0,r0,m,r,v0)
    data=simulationData(parr,n,interval,start,RC)

    for d in data:
        print("ok")
        plotCraterFormation(d,r0,r)
    print("done")

def craterSimulation(N,interval,m0,r0,m,r,v0,RC):

    #initialization
    parr=initParticlesForProjectileImpact(N,m0,r0,m,r,v0)
    initialEnergy=averageKineticEnergy(parr[:N-1])+avrgEnergyP(parr[N]) #ground particles start with zero kinetic energy!!
    print("initE",avrgEnergyP(parr[N]))
    events=firstCollisions(parr)
    Time=0.0 #set start time to zero

    #run simulation
    print("start")
    while((averageKineticEnergy(parr[:N-1])+avrgEnergyP(parr[N]))>=0.1*initialEnergy): #continue until 10% or less than initial energy
        print((averageKineticEnergy(parr[:N-1])+avrgEnergyP(parr[N]))/initialEnergy)
        Time = runSimulation(interval,events,parr,Time,RC) #run simulation and update time
    print("done")
    return parr
 

    


"""values for the different parameters"""
N=1000 #number of ground particles
n=10000 #number of collisions in simulation
interval=100 #interval between data sampling
start=0 #start of data sampling
m0=1 #mass of ground particles
r0=np.sqrt(1/(N*4*math.pi)) #radius of mass particles
m=m0*25 #mass of projectile partile (as stated in problem text)
r=r0*5 #radius of projectile particle
v0=5 #initial velocity of projectile particle
RC=0.5



#problem4Visualization(N,n,interval,start,m0,r0,m,r,v0,RC)

craterSimulation(N,interval,m0,r0*0.9,m,r,v0,RC)