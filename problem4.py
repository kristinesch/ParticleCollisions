from simulation import *
from functions import *


#NOTE TO SELF: might be more effetive to write a function that just gives the positions, since these take time to calculate

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

    fig,ax=plt.subplots(1,1)
    ax.axis([Xmin,Xmax,Ymin,Ymax])
    ax.scatter(groundParticles[:,X],groundParticles[:,Y],c="b")
    ax.scatter(projectile[X],projectile[Y],c="r",s=r*r*500*500)
    fig.savefig("p4Visualization")
    plt.show()


def problem4Visualization(N,n,interval,start,m0,r0,m,r,v0,RC):
    parr=initParticlesForProjectileImpact(N,m0,r0,m,r,v0)
    data=simulationData(parr,n,interval,start,RC)

    for d in data:
        print("ok")
        plotCraterFormation(d,r0,r)
    print("done")

def craterVisualization(parr,n,interval,start,RC):
    # Xmin=0
    # Xmax=1
    # Ymin=0
    # Ymax=1
    #parr=initParticlesForProjectileImpact(N,m0,r0,m,r,v0)
    data=simulationData(parr,n,interval,start,RC)
    fig,ax=plt.subplots(2,2)

    parr1=data[0]
    groundParticles1=parr1[:N-1]
    projectile1=parr1[N-1]
    ax[0,0].axis([Xmin,Xmax,Ymin,Ymax])
    ax[0,0].scatter(groundParticles1[:,X],groundParticles1[:,Y],c="b",s=0.01)
    ax[0,0].scatter(projectile1[X],projectile1[Y],c="r",s=r*r*50*50)

    parr2=data[1]
    groundParticles2=parr2[:N-1]
    projectile2=parr2[N-1]
    ax[0,1].axis([Xmin,Xmax,Ymin,Ymax])
    ax[0,1].scatter(groundParticles2[:,X],groundParticles2[:,Y],c="b",s=0.01)
    ax[0,1].scatter(projectile2[X],projectile2[Y],c="r",s=r*r*50*50)

    parr3=data[2]
    groundParticles3=parr3[:N-1]
    projectile3=parr3[N-1]
    ax[1,0].axis([Xmin,Xmax,Ymin,Ymax])
    ax[1,0].scatter(groundParticles3[:,X],groundParticles3[:,Y],c="b",s=0.01)
    ax[1,0].scatter(projectile3[X],projectile3[Y],c="r",s=r*r*50*50)

    parr4=data[3]
    groundParticles4=parr4[:N-1]
    projectile4=parr4[N-1]
    ax[1,1].axis([Xmin,Xmax,Ymin,Ymax])
    ax[1,1].scatter(groundParticles4[:,X],groundParticles4[:,Y],c="b",s=0.01)
    ax[1,1].scatter(projectile4[X],projectile4[Y],c="r",s=r*r*50*50)

    fig.savefig("p4Visualization")
    plt.show()



def craterSimulation(parr,N,interval,RC):

    #initialization
    initialEnergy=totKinEnergy(parr) #ground particles start with zero kinetic energy!!
    print("initE",initialEnergy)
    events=firstCollisions(parr)
    Time=0.0 #set start time to zero
    #run simulation
    print("start")
    while((totKinEnergy(parr))>=0.1*initialEnergy): #continue until 10% or less than initial energy
        print(totKinEnergy(parr)/(initialEnergy))
        Time = runSimulation(interval,events,parr,Time,RC) #run simulation and update time
    print("done")
    return parr

#notes: projectil energy decreases the whole time
#ground particle energy increases, because they get energy from the projectile
 





def massCraterSize(initialParr,masslist,N,interval,RC): #masslist is a list of masses to plot for #filename=name of file data is saved to
    #initialize array
    #initialParr=initParticlesForProjectileImpact(N,m0,r0,1,r,v0) #set m to 1, can change later
    print(masslist)
    craterData=[]
    for mass in masslist:

        #initialization
        parr=initialParr.copy() #copy intitial particle array
        parr[N,M]=mass
        print("current mass",parr[N,M])
        initialEnergy=totKinEnergy(parr)

        events=firstCollisions(parr)
        Time=0.0 #set start time to zero

        newInterval=int(interval*(mass/m0)) #to make more appropriate intervals
        print("interval",newInterval)

        #run simulation
        print("start")
        while((totKinEnergy(parr))>=0.1*initialEnergy): #continue until 10% or less than initial energy
            print(totKinEnergy(parr)/(initialEnergy))
            Time = runSimulation(newInterval,events,parr,Time,RC) #run simulation and update time
        crater=craterSize(initialParr,parr) #calculate crater size
        craterData.append(crater) 
    npCrater=np.array(craterData)
    np.save("massCraterData.npy",masslist,npCrater) #save data to binary file


    


"""values for the different parameters"""
N=500 #number of ground particles
n=8000 #number of collisions in simulation
interval=10 #interval between data sampling
start=0 #start of data sampling
m0=1#mass of ground particles
r0=np.sqrt(1/(N*4*math.pi)) #radius of mass particles
m=m0*10 #mass of projectile partile (as stated in problem text)
r=r0*5 #radius of projectile particle
v0=5 #initial velocity of projectile particle
RC=0.5



#problem4Visualization(N,n,interval,start,m0,r0,m,r,v0,RC)

#0.9 for faster init!!!!!


massList=np.linspace(m0,30*m0,5)

"""Initializing particle array"""
#parr=initParticlesForProjectileImpact(N,m0,r0,m,r,v0)
#np.save("particleArray.npy",parr) #save array as binary file

initParticleArray=np.load("particleArray.npy")
print("ok")
newArray=initParticleArray.copy()
newArray[N,M]=m



#craterSimulation(newArray,N,100,RC)

craterVisualization(initParticleArray,n,2000,0,RC)

#massCraterSize(initParticleArray,massList,N,interval,RC)

#plotParameterVsCraterSize("MassCrater.npy","Projectile Mass","MassVsCraterPlot")


#code for plotting for different parameters: take parr as input
