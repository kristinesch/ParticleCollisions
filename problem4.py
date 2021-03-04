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



def craterSimulation(parr,N,interval,m0,r0,m,r,v0,RC):

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
 





def massCraterSize(initialParr,masslist,filename,N,interval,m0,r0,r,v0,RC): #masslist is a list of masses to plot for #filename=name of file data is saved to
    #initialize array
    #initialParr=initParticlesForProjectileImpact(N,m0,r0,1,r,v0) #set m to 1, can change later

    file=open(filename,"w") #file to write to

    print(masslist)
    for mass in masslist:

        #initialization
        parr=initialParr.copy() #copy intitial particle array
        parr[N,M]=mass
        print("current mass",parr[N,M])
        initialEnergy=totKinEnergy(parr)

        events=firstCollisions(parr)
        Time=0.0 #set start time to zero

        newInterval=int(interval*(mass/m0))
        print("interval",newInterval)
        #run simulation
        print("start")
        while((totKinEnergy(parr))>=0.1*initialEnergy): #continue until 10% or less than initial energy
            print(totKinEnergy(parr)/(initialEnergy))
            Time = runSimulation(newInterval,events,parr,Time,RC) #run simulation and update time
            

        crater=craterSize(initialParr,parr)
        file.write(str(crater)) #writing data to file
        print("data saved")
    file.close()


    


"""values for the different parameters"""
N=1000 #number of ground particles
n=10000 #number of collisions in simulation
interval=100 #interval between data sampling
start=0 #start of data sampling
m0=1 #mass of ground particles
r0=np.sqrt(1/(N*4*math.pi)) #radius of mass particles
m=m0*15 #mass of projectile partile (as stated in problem text)
r=r0*5 #radius of projectile particle
v0=5 #initial velocity of projectile particle
RC=0.5



#problem4Visualization(N,n,interval,start,m0,r0,m,r,v0,RC)

#0.9 for faster init!!!!!


#massList=np.linspace(m0,40*m0,10)

"""Initializing particle array"""
#parr=initParticlesForProjectileImpact(N,m0,r0,m,r,v0)
#np.save("particleArray.npy",parr) #save array as binary file

initParticleArray=np.load("particleArray.npy")
print("ok")

craterSimulation(initParticleArray,N,interval,m0,r0,m,r,v0,RC)


#massCraterSize(initialParr,massList,"massVSCraterSize.txt",N,interval,m0,r0*0.9,r,v0,RC)


#code for plotting for different parameters: take parr as input
