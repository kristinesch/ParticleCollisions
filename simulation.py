from matplotlib import numpy as np 
from matplotlib import pyplot as plt
from random import seed
import random
import math
import heapq
from numba import jit

#%load_ext line_profiler
random.seed(1)

"""Colission count:
The first element for each particle list is the particles current collision count
Update the particles collision count for every collision
Collision count used to check if valid collision
Each collision event is saved as a list where the two last elements are the collision counts for the involved particle(s)
If two particles: both their collision count at the time of collision is saved
If wall collisoin: first entry is the collsiion count of the particle, last entry doesnt matter 

"""

"""event queue structure"""
#A queue of events. Each event is a list containg time untill next collision(deltaT)
#index i of particle, index j of wall or particle, collision count of particle at time of collision
#if j=-2: vertical wall. if j=-1: horizontal wall. Else: particle collision


#index consts:
Dt=0 #Time of collision
I=1 #index of particle 1 in the particle array
J=2 #index of particle 2
CCI=3 #collision count of particle 1
CCJ=4 #collision count of particle 2
#PL=5 #particle list: list over all other particles that have registered a collision with this particle REMOVE?


"""Particles"""
#Saved in array, one column corresponds to one particle

#index constants:
C=0 #collision count
X=1 #position (X,Y)
Y=2
VX=3 #velocity (VX,VY)
VY=4
R=5 #radius
M=6 #mass


#REMOVE?????????????
"""constants""" #NOTE: might change into parameters in the functions...
RC=1  #restitution coefficient (degree of inelasticity)
Xmin=0
Xmax=1  #max. x-value of position
Ymin=0
Ymax=1
r=0.0001 #Should be much less than 1!!! 
m=1  #just picked some values for r and m. Might change them later
V0=1 #Absolute value of velocity


"""Iitialization of particle array"""
#N is the number of particles, the other arguments are the max. and min values for the position coords
#!!!!!!Assumes equal radii for all particles!!!!!!!!!!

#initializes the particle array, with random posistions
def initParticleArray(N,masses,r,m,v0,rc,xmin=0,xmax=1,ymin=0,ymax=1): #masses is an array of masses
    print("Initializing particle array...")
    parr=np.zeros((N,7))
    occupiedPositions=[] #list of occupied positions

    for i in range(N): #assign random values for positions
        Overlapping=True
        while(Overlapping==True): #checking if position is already occupied
                pos=np.array([random.uniform(r,1-r),random.uniform(r,1-r)])
                Overlapping=False
                for posi in occupiedPositions:
                    if (4*r*r>np.dot((posi-pos),(posi-pos))):
                        Overlapping=True
                        break

        parr[i,X]=pos[0]
        parr[i,Y]=pos[1]
        
        #velocities
        theta=random.uniform(0,2*math.pi)
        parr[i,VX]=V0*math.cos(theta)
        parr[i,VY]=V0*math.sin(theta)
        parr[i,M]=masses[i] #masses
        parr[i,R]=r

        occupiedPositions.append(pos)
        print((i*100)/N,"%")
    print("Initialization done")
    return parr #returns initialized particle array




"""Plotting/Visualization"""

#plotting the particle positions (used for testing)
def plotParticlePositions(parr):
    plt.axis([Xmin,Xmax,Ymin,Ymax])
    plt.plot(parr[:,X],parr[:,Y],"ro")
    plt.show()

#used for testing
def plotPositionsAndVelocities(parr):
    plt.plot(parr[:,X],parr[:,Y],"ro",markersize=r*200)
    plt.plot([parr[:,X],parr[:,X]+parr[:,VX]],[parr[:,Y],parr[:,Y]+parr[:,VY]]) #testing
    plt.axis([Xmin,Xmax,Ymin,Ymax])
    plt.plot([0,1],[1,1]) #helping line
    plt.show()



"""FUNCTIONS"""

#moves system forward in time: update positions: (x+vx*dt) 
def moveForwardInTime(parr,dt):
    parr[:,X]=parr[:,X]+parr[:,VX]*dt
    parr[:,Y]=parr[:,Y]+parr[:,VY]*dt



#updating the velocity values after collision
#the input parr is the particle array

def verticalCollision(parr,i,dt,RC): #i=particle index, dt = remaining time until the collision
    moveForwardInTime(parr,dt)
    #updating velocities
    parr[i,VX]=-RC*parr[i,VX]
    parr[i,VY]=RC*parr[i,VY]
    #update collision count
    parr[i,C]+=1

def horizontalCollision(parr,i,dt,RC): 
    moveForwardInTime(parr,dt)
    parr[i,VX]=RC*parr[i,VX]
    parr[i,VY]=-RC*parr[i,VY]
    #update collision count
    parr[i,C]+=1




#if particle collides with vertical wall, add event to heapq
#@jit(nopython = True)
def nextVerticalCollision(parr,i,events,time): #i is the particle index, events is heapq of events, time=current time
    p=parr[i]
    if p[VX]>0:
        deltaT=(1-p[R]-p[X])/p[VX]
        e=[deltaT+time,i,-2,parr[i,C]]
        heapq.heappush(events,e)
    elif p[VX]<0:
        deltaT=(p[R]-p[X])/p[VX]
        e=[deltaT+time,i,-2,parr[i,C]] 
        heapq.heappush(events,e)


#@jit(nopython = True)
def nextHorizontalCollision(parr,i,events,time):
    p=parr[i]
    if p[VY]>0:
        deltaT=(1-p[R]-p[Y])/p[VY]
        e=[deltaT+time,i,-1,parr[i,C]] 
        heapq.heappush(events,e)
    elif p[VY]<0:
        deltaT=(p[R]-p[Y])/p[VY]
        e=[deltaT+time,i,-1,parr[i,C]]
        heapq.heappush(events,e)
 



#optimalize calculation speed
@jit(nopython = True)
def computeStuff(i,j,parr):
    Rij=parr[i,R]+parr[j,R] #radiuses added together
    deltaX=np.array([(parr[j,X]-parr[i,X]),(parr[j,Y]-parr[i,Y])])
    deltaV=np.array([(parr[j,VX]-parr[i,VX]),(parr[j,VY]-parr[i,VY])]) #calculation as given in problem text
    VdotX=(np.dot(deltaV,deltaX))
    d=(VdotX)*(VdotX)-(np.dot(deltaV,deltaV))*((np.dot(deltaX,deltaX))-Rij*Rij)
    return deltaX,deltaV,d,VdotX 


#@jit(nopython = True)
def particlesColTime(i,j,parr,events,time): #p1 and p2 is the index of the particles rows in parr 
    deltaX, deltaV, d,VdotX=computeStuff(i,j,parr)

    if (d>0 and VdotX<0):
        #calculate new time
        #deltaT=calculateDt(deltaV,deltaX,d)
        #print(np.dot(deltaV,deltaX))
        deltaT=-(np.dot(deltaV,deltaX)+np.sqrt(d))/np.dot(deltaV,deltaV)
        e=[deltaT+time,i,j,parr[i,C],parr[j,C]]  #last 2 entries are the current collisioncount for particles i and
        heapq.heappush(events,e)  #add event to heapq




def particleCollision(parr,i,j,dt,RC):  #i and j are particle indeces
    #BEFORE collision:
    vi=[parr[i,VX],parr[i,VY]] #saving original velocity!
    vj=[parr[j,VX],parr[j,VY]]
    #move system forward in time: 
    moveForwardInTime(parr,dt)
    #update velocities for colliding particles:
    mi=parr[i,M]
    mj=parr[j,M]    
    Rij=parr[i,R]+parr[j,R]
    newDeltaV=[(parr[j,VX]-parr[i,VX]),(parr[j,VY]-parr[i,VY])]
    newDeltaX=[(parr[j,X]-parr[i,X]),(parr[j,Y]-parr[i,Y])]
    #new velocities calculated from given formula:
    newVi=vi+np.dot(((1+RC)*(mj/(mi+mj))*(np.dot(newDeltaV,newDeltaX))/pow(Rij,2)),newDeltaX)
    newVj=vj-np.dot(((1+RC)*(mi/(mi+mj))*(np.dot(newDeltaV,newDeltaX))/pow(Rij,2)),newDeltaX)
    #updating the particle array:
    parr[i,VX]=newVi[0]
    parr[i,VY]=newVi[1]
    parr[j,VX]=newVj[0]
    parr[j,VY]=newVj[1]
    #updating collsion count for both particles:
    parr[i,C]+=1
    parr[j,C]+=1




"""Initialization"""

def firstCollisions(parr): #takes in the particle array, returns the heapq events
    events=[]
    N=len(parr)
    for i in range(N): #for row i in parr (particle i)
        if(N>1):
            for j in range(N):
                if i!=j: #if different particles:
                    particlesColTime(i,j,parr,events,0) 
        nextHorizontalCollision(parr,i,events,0)
        nextVerticalCollision(parr,i,events,0)
    return events


#@jit(nopython = True)
def nextCollisions(parr,events,i,time): #calculating the next collisions for particle i 
    N=len(parr)
    if(N>1):
        for j in range(N):
            if (i!=j): #if different particles:     
                particlesColTime(i,j,parr,events,time) 
    nextHorizontalCollision(parr,i,events,time)
    nextVerticalCollision(parr,i,events,time)



"""LOOP"""

def runSimulation(n,events,parr, t0,RC): #n is the number of valid iterations/collisions, t0=start time
    k=1 #for counting
    Time=t0 #current time
    while k<n:
        #first let the first valid collision happen:
        #pop removes the event from the event list
        cc=heapq.heappop(events) # cc = currentCollision
        CollisionTime=cc[Dt] #DT=time+dt
        dt=CollisionTime-Time #time until next collision

        if (parr[cc[I],C]==cc[CCI]): #if valid for 1st particle

            if (cc[J]==-2): #if collision with vertical wall
                verticalCollision(parr,cc[I],dt,RC) 
                Time=CollisionTime #updating time
                nextCollisions(parr,events,cc[I],CollisionTime) #calculate new collisions for the colliding particle
                k+=1

            elif (cc[J]==-1): #if collision with horizontal wall
                horizontalCollision(parr,cc[I],dt,RC)
                Time=CollisionTime #updating time
                nextCollisions(parr,events,cc[I],CollisionTime) #calculate new collisions for the colliding particle:
                k+=1

            elif (parr[cc[J],C]==cc[CCJ]): #if valid for second particle
                    particleCollision(parr,cc[I],cc[J],dt,RC)
                    nextCollisions(parr,events,cc[J],CollisionTime) #calculate new next collisions
                    Time=CollisionTime #updating time
                    nextCollisions(parr,events,cc[I],CollisionTime) #calculate new collisions for the colliding particle
                    k+=1

            if (k%1000==0): #just to keep track of how much time is left
                print("collision number",k)
            
    return Time     #return current time at the end of the simulation    

                     

#just runs the simulation for a given particle array
def simulation(particleArray,n,RC): #particles, number of iterations in simulation
    eventQueue=firstCollisions(particleArray)
    runSimulation(n,eventQueue,particleArray,0,RC)
    return particleArray

#runs the simulation for a given particle array, but saves data at given interval
def simulationData(particleArray,n,interval,start,RC): #start=number of collisions before starting sampling, interval=number of collisions between sampling
    eventQueue=firstCollisions(particleArray)
    data=[]
    print("start")
    if (start==0): #copy initial array only if start of sampling set to zero
        data.append(particleArray.copy())
    Time=0.0 #set start time to zero
    for i in range(int(n/interval)+1):
        Time = runSimulation(interval,eventQueue,particleArray, Time,RC) #run simulation and update time
        if ((i+1)*interval>=start): #save data
            print("data saved")
            data.append(particleArray.copy()) 
    return data


"""TESTING STUFF"""
# test=simulation(1000,5000)
# print(test)
# plotPositionsAndVelocities(test)

#p1=simulation(100,100) 

#%lprun -f test x = test(10, 10)

#parrTest=initParr(100)
#eventTest=firstCollisions(parrTest)

#test(3,3)
#%lprun -f firstCollisions x = firstCollisions(parrTest)

#%lprun -f test x = test(1000,1000)


#treigDrit(1,1,parrTest)

#%lprun -f runSimulation x = runSimulation(100,eventTest,parrTest)

#%lprun -f particlesColTime x = particlesColTime(1,1,parrTest,eventTest)
#%lprun -f treigDrit x = treigDrit(1,1,parrTest)


