
#problem: koden funker ikke med jit????? får index error?

from matplotlib import numpy as np 
from matplotlib import pyplot as plt
#from random import random
from random import seed
import random
import math
import heapq
from numba import jit

#%load_ext line_profiler
random.seed(1)


#hva med kollisjonene som forsvinner pga invalidity? det Sigrid forklarte

"""Colission count:
the first element for each particle list is the particles current collision count
Updtae the particles collisoion count for every collision
Collision count used to check if valid collision
Each collision event is saved as a list where the two last elements are the collision counts for the involved particle(s)
If two particles: both their collision count at the time of collision is saved
If wall collisoin: first entry is the collsiion count of the particle, last entry doesnt matter (no last entry needed?)

"""

"""event queue structure"""
#A queue of events. Each event is a list containg time untill next collision(deltaT)
#index i of particle, index j of wall or particle, collision count of particle at time of collision
#if j=-2: vertical wall. if j=-1: horizontal wall. Else: particle collision


#index consts:
Dt=0 #Time of collision
I=1 #index of particle 1 in the particle array
J=2 #index of article 2
CCI=3 #collision count of particle 1
CCJ=4 #collision count of particle 2
PL=5 #particle list: list over all other particles that have registered a collision with this particle

"""constants"""
 #number of particles
RC=1 #restitutuion coefficient (degree of inelasticity)
Xmin=0
Xmax=1  #max. x-value of position
Ymin=0
Ymax=1
r=0.001 #Should be much less than 1!!! m
m=1  #foreløpige verdier for m og r. Hvis de skal variere, kan det fikses senere
V0=1


"""Particles"""
#Saved in array, one column corresponds to one particle

#index constants:
C=0 #collision count
X=1
Y=2
VX=3
VY=4
R=5
M=6


#Iitialization of particle array
#N is the number of particles, the other arguments are the max. and min values for the position coords

#NB skrevet for like radiuser!

#initializes the particle array, with random posistions
def initParticleArray(N,xmin=0,xmax=1,ymin=0,ymax=1):
    print("Initializing particle array...")
    parr=np.zeros((N,7))
    occupiedPositions=[] #list of occupied positions

    for i in range(N): #gir randome verdier til alle hastigheter og posisjoner (ikke r eller m)
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
        parr[i,M]=m 
        parr[i,R]=r

        occupiedPositions.append(pos)
        print((i*100)/N,"%")
    print("Initialization done")
    return parr




"""Plotting/Visualization"""

def plotPositionsAndVelocities(parr):
    plt.plot(parr[:,X],parr[:,Y],"ro",markersize=r*200)
    plt.plot([parr[:,X],parr[:,X]+parr[:,VX]],[parr[:,Y],parr[:,Y]+parr[:,VY]]) #testing
    plt.axis([Xmin,Xmax,Ymin,Ymax])
    plt.plot([0,1],[1,1]) #hjelpelinje
    plt.show()

"""Functions"""

def moveForwardInTime(parr,dt):
    parr[:,X]=parr[:,X]+parr[:,VX]*dt
    parr[:,Y]=parr[:,Y]+parr[:,VY]*dt


#updating the velocity values after collision
#the input parr is the particle array

def verticalCollision(parr,i,dt): #i is the particle index
    moveForwardInTime(parr,dt)
    #updating velocities
    parr[i,VX]=-RC*parr[i,VX]
    parr[i,VY]=RC*parr[i,VY]
    #update collision count
    parr[i,C]+=1

def horizontalCollision(parr,i,dt): 
    moveForwardInTime(parr,dt)
    parr[i,VX]=RC*parr[i,VX]
    parr[i,VY]=-RC*parr[i,VY]
    parr[i,C]+=1



#NB: DENNE KAN NOK OPTIMALISERES!!!!

#if particle collides with vertical wall, add event to heapq
#@jit(nopython = True)
def nextVerticalCollision(parr,i,time): #i is the particle index, events is heapq of events
    p=parr[i]
    if p[VX]>0:
        deltaT=(1-p[R]-p[X])/p[VX]
        e=[deltaT+time,i,-2,parr[i,C]]
        return e
    elif p[VX]<0:
        deltaT=(p[R]-p[X])/p[VX]
        e=[deltaT+time,i,-2,parr[i,C]] #SISTE ELEMENT HER ER CURRENT PARTICLE COUNT
        #nb går utifra at lengden på hvert element ikke har noe å si...
        return e
        #print("deltaT",deltaT,"time",time,"\n")

#@jit(nopython = True)
def nextHorizontalCollision(parr,i,time):
    p=parr[i]
    if p[VY]>0:
        deltaT=(1-p[R]-p[Y])/p[VY]
        e=[deltaT+time,i,-1,parr[i,C]] #første element er tiden kollisjonen skjer ved, altså ikke tid til neste kollisjon!
        return e
    elif p[VY]<0:
        deltaT=(p[R]-p[Y])/p[VY]
        e=[deltaT+time,i,-1,parr[i,C]]
        return e
 


#GANSKE JALLA Å BARE TA UT DELER AV FUNKSJONEN FOR Å JITTE DE MEN WHATEVER

@jit(nopython = True)
def treigDrit(i,j,parr):
    Rij=parr[i,R]+parr[j,R]  #radiusene lagt sammen LITT USIKKER PÅ DETTE
    deltaX=np.array([(parr[j,X]-parr[i,X]),(parr[j,Y]-parr[i,Y])])
    deltaV=np.array([(parr[j,VX]-parr[i,VX]),(parr[j,VY]-parr[i,VY])])
    VdotX=(np.dot(deltaV,deltaX))
    d=(VdotX)*(VdotX)-(np.dot(deltaV,deltaV))*((np.dot(deltaX,deltaX))-Rij*Rij)
    return deltaX,deltaV,d,VdotX

"""
@jit(nopython = True)
def calculateDt(deltaV,deltaX,d):
    return -(np.dot(deltaV,deltaX)+np.sqrt(d))/np.dot(deltaV,deltaV)
    """

#@jit(nopython = True)
#FEIL svar når partiklene er inni hverandre
#@jit(nopython = True)
def particlesColTime(i,j,parr,time): #p1 and p2 is the index of the particles rows in parr 
    #deltaX, deltaV, d,VdotX=treigDrit(i,j,parr)
    Rij=parr[i,R]+parr[j,R]  #radiusene lagt sammen LITT USIKKER PÅ DETTE
    deltaX=np.array([(parr[j,X]-parr[i,X]),(parr[j,Y]-parr[i,Y])])
    deltaV=np.array([(parr[j,VX]-parr[i,VX]),(parr[j,VY]-parr[i,VY])])
    VdotX=(np.dot(deltaV,deltaX))
    d=(VdotX)*(VdotX)-(np.dot(deltaV,deltaV))*((np.dot(deltaX,deltaX))-Rij*Rij)

    if (d>0 and VdotX<0):
        #calculate new time
        deltaT=-(np.dot(deltaV,deltaX)+np.sqrt(d))/np.dot(deltaV,deltaV)
        e=[deltaT+time,i,j,parr[i,C],parr[j,C]]  #last 2 entries are the current collisioncount for particles i and
        return e  



#NB flytter systemet frem i tid OG oppdaterer hastigheter
def particleCollision(parr,i,j,dt):  #i and j are particle indeces
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
                    events.append(particlesColTime(i,j,parr,0))
        events.append(nextHorizontalCollision(parr,i,0))
        events.append(nextVerticalCollision(parr,i,0))
        #print(events)
    eventQueue=[]
    for e in events:
        if e!=None:
            #print(e)
            eventQueue.append(e)
    heapq.heapify(eventQueue)
    return eventQueue #return a heapqueue of the events

#@jit(nopython = True)
def nextCollisions(parr,i,time): #TAR INN PARTIKKEL-INDEKS!! og tid!
    events=[]
    N=len(parr)
    if(N>1):
        for j in range(N):
            if (i!=j): #if different particles:     #BURDE FINNE NOE MER Å SJEKKE HER
                #FOR Å IKKE KALLE PÅ TREIGDIRT UNØDVENDIG
                events.append(particlesColTime(i,j,parr,0))
        events.append(nextHorizontalCollision(parr,i,0))
        events.append(nextVerticalCollision(parr,i,0))
        
    return events




"""LOOP"""

def runSimulation(n,eventQueue,parr): #n is the number of valid iterations/collisions
    print("DRITT")
    k=1 #for counting
    Time=0 #current time
    print("start loop")
    while k<n:
        #first let the first valid collision happen:
        #pop removes the event from the event list
        cc=heapq.heappop(eventQueue) # cc = currentCollision
        #print("cc",cc)
        CollisionTime=cc[Dt] #DT=time+dt
        dt=CollisionTime-Time #time until next collision
        newEvents=[]
        #print("DT",DT,"time",Time,"dt",dt)
        print("P",[cc[I]],cc)
        if (parr[cc[I],C]==cc[CCI]): #if valid for 1st particle
            print("COLLISION")
            if (cc[J]==-2): #if collision with vertical wall
                verticalCollision(parr,cc[I],dt) #DT-time er tid til neste kollisjon
                Time=CollisionTime #updating time
                newEvents.append(nextCollisions(parr,cc[I],CollisionTime)) #calculate new collisions for the colliding particle:
                k+=1
                #print("Collision number: ",k," Time: ",Time,"\n")
            elif (cc[J]==-1): #if collision with horizontal wall
                horizontalCollision(parr,cc[I],dt)
                Time=CollisionTime #updating time
                newEvents.append(nextCollisions(parr,cc[I],CollisionTime)) #calculate new collisions for the colliding particle:
                k+=1
                #print("Collision number: ",k," Time: ",Time,"\n")
            elif (parr[cc[J],C]==cc[CCJ]): #if valid for second particle
                    particleCollision(parr,cc[I],cc[J],dt)
                    #newEvents.append(nextCollisions(parr,cc[J],CollisionTime)) #with updated time
                    Time=CollisionTime #updating time
                    newEvents.append(nextCollisions(parr,cc[I],CollisionTime)) #calculate new collisions for the colliding particle:
                    k+=1
                    #print("Collision number: ",k," Time: ",Time,"\n")
            #for i in range(len(parr)):
                #newEvents.append(nextCollisions(parr,i,CollisionTime))
            if (k%1000==0):
                print("collision number",k)
            for elist in newEvents:
                for e in elist:
                    if e!=None:
                        heapq.heappush(eventQueue,e)
            
            
            

            #plotting for each iteration        
            #plotPositionsAndVelocities(parr)
            

        


def simulation(N,n): #particles, number of iterations in simulation
    particleArray=initParticleArray(N)
    eventQueue=firstCollisions(particleArray)
    #plotPositionsAndVelocities(parrTest)
    runSimulation(n,eventQueue,particleArray)
    #plotPositionsAndVelocities(parrTest)
    print("Done")
    return particleArray


test=simulation(1000,5000)
print(test)
plotPositionsAndVelocities(test)



# parrTest=initParticleArray(50)
# eventTest=firstCollisions(parrTest)
# %lprun -f runSimulation x = runSimulation(100,eventTest,parrTest)



