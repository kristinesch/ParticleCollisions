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
def initParticleArray(N,masses,xmin=0,xmax=1,ymin=0,ymax=1): #masses is an array of masses
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
        parr[i,M]=masses[i]
        parr[i,R]=r

        occupiedPositions.append(pos)
        print((i*100)/N,"%")
    print("Initialization done")
    return parr




"""Plotting/Visualization"""
def plotParticlePositions(parr):
    plt.axis([Xmin,Xmax,Ymin,Ymax])
    plt.plot(parr[:,X],parr[:,Y],"ro",s=5)


def plotPositionsAndVelocities(parr):
    #print(parr)
    #fig, ax = plt.subplots()
    plt.plot(parr[:,X],parr[:,Y],"ro",markersize=r*200)
    #print("VY",parr[:,X].tolist(),parr[:,Y].tolist(),parr[:,VX].tolist(),parr[:,VY].tolist())
    plt.plot([parr[:,X],parr[:,X]+parr[:,VX]],[parr[:,Y],parr[:,Y]+parr[:,VY]]) #testing
    #ax.quiver(parr[:,X].tolist(),parr[:,Y].tolist(),parr[:,VX].tolist(),parr[:,VY].tolist())
    #
    plt.axis([Xmin,Xmax,Ymin,Ymax])
    plt.plot([0,1],[1,1]) #hjelpelinje
    plt.show()

"""Functions"""

def moveForwardInTime(parr,dt):
    #move system forward in time: update positions: (x+vx*dt) and time????
    #print("Old Position:  X:",parr[:,X],"Y: ",parr[:,Y],"dt",dt)
    parr[:,X]=parr[:,X]+parr[:,VX]*dt
    parr[:,Y]=parr[:,Y]+parr[:,VY]*dt
    #print("New Position",parr[:,X],parr[:,Y])


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
    #update collision count
    parr[i,C]+=1



#NB: DENNE KAN NOK OPTIMALISERES!!!!

#if particle collides with vertical wall, add event to heapq
#@jit(nopython = True)
def nextVerticalCollision(parr,i,events,time): #i is the particle index, events is heapq of events
    p=parr[i]
    if p[VX]>0:
        deltaT=(1-p[R]-p[X])/p[VX]
        e=[deltaT+time,i,-2,parr[i,C]]
        heapq.heappush(events,e)
    elif p[VX]<0:
        deltaT=(p[R]-p[X])/p[VX]
        e=[deltaT+time,i,-2,parr[i,C]] #SISTE ELEMENT HER ER CURRENT PARTICLE COUNT
        #nb går utifra at lengden på hvert element ikke har noe å si...
        heapq.heappush(events,e)
        #print("deltaT",deltaT,"time",time,"\n")

#@jit(nopython = True)
def nextHorizontalCollision(parr,i,events,time):
    p=parr[i]
    if p[VY]>0:
        deltaT=(1-p[R]-p[Y])/p[VY]
        e=[deltaT+time,i,-1,parr[i,C]] #første element er tiden kollisjonen skjer ved, altså ikke tid til neste kollisjon!
        heapq.heappush(events,e)
        #print("deltaT",deltaT,"time",time,"\n")
        #print(time)
    elif p[VY]<0:
        deltaT=(p[R]-p[Y])/p[VY]
        #print("deltaT",deltaT,"time",time,"\n")
        #print(time)
        e=[deltaT+time,i,-1,parr[i,C]]
        heapq.heappush(events,e)
 


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
def particlesColTime(i,j,parr,events,time): #p1 and p2 is the index of the particles rows in parr 
    deltaX, deltaV, d,VdotX=treigDrit(i,j,parr)

    if (d>0 and VdotX<0):
        #calculate new time
        #deltaT=calculateDt(deltaV,deltaX,d)
        #print(np.dot(deltaV,deltaX))
        deltaT=-(np.dot(deltaV,deltaX)+np.sqrt(d))/np.dot(deltaV,deltaV)
        e=[deltaT+time,i,j,parr[i,C],parr[j,C]]  #last 2 entries are the current collisioncount for particles i and
        heapq.heappush(events,e)  #add event to heapq
        #print("particlecoltime",deltaT)

        #(parr[i,PL]).append(j)
        #(parr[j, PL]).append(i)



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
                    particlesColTime(i,j,parr,events,0) 
        nextHorizontalCollision(parr,i,events,0)
        nextVerticalCollision(parr,i,events,0)
    return events

#@jit(nopython = True)
def nextCollisions(parr,events,i,time): #TAR INN PARTIKKEL-INDEKS!! og tid!
    N=len(parr)
    if(N>1):
        for j in range(N):
            if (i!=j): #if different particles:     #BURDE FINNE NOE MER Å SJEKKE HER
                #FOR Å IKKE KALLE PÅ TREIGDIRT UNØDVENDIG
                particlesColTime(i,j,parr,events,time) 
    nextHorizontalCollision(parr,i,events,time)
    nextVerticalCollision(parr,i,events,time)



"""LOOP"""

def runSimulation(n,events,parr, t0): #n is the number of valid iterations/collisions
    k=1 #for counting
    Time=t0 #current time
    #print("start loop")
    while k<n:
        #first let the first valid collision happen:
        #pop removes the event from the event list
        cc=heapq.heappop(events) # cc = currentCollision
        #print("cc",cc)
        CollisionTime=cc[Dt] #DT=time+dt
        dt=CollisionTime-Time #time until next collision
        #print("DT",DT,"time",Time,"dt",dt)
        #print("P",[cc[I]],cc)
        if (parr[cc[I],C]==cc[CCI]): #if valid for 1st particle
            if (cc[J]==-2): #if collision with vertical wall
                verticalCollision(parr,cc[I],dt) #DT-time er tid til neste kollisjon
                Time=CollisionTime #updating time
                nextCollisions(parr,events,cc[I],CollisionTime) #calculate new collisions for the colliding particle:
                k+=1
                #print("Collision number: ",k," Time: ",Time,"\n")
            elif (cc[J]==-1): #if collision with horizontal wall
                horizontalCollision(parr,cc[I],dt)
                Time=CollisionTime #updating time
                nextCollisions(parr,events,cc[I],CollisionTime) #calculate new collisions for the colliding particle:
                k+=1
                #print("Collision number: ",k," Time: ",Time,"\n")
            elif (parr[cc[J],C]==cc[CCJ]): #if valid for second particle
                    particleCollision(parr,cc[I],cc[J],dt)
                    nextCollisions(parr,events,cc[J],CollisionTime) #with updated time

                    Time=CollisionTime #updating time
                    nextCollisions(parr,events,cc[I],CollisionTime) #calculate new collisions for the colliding particle:
                    k+=1
                    #print("Collision number: ",k," Time: ",Time,"\n")
            if (k%1000==0):
                print("collision number",k)
            
    return Time        
            #plotting for each iteration        
            #plotPositionsAndVelocities(parr)
            

        


def simulation(particleArray,N,n): #particles, number of iterations in simulation
    eventQueue=firstCollisions(particleArray)
    #plotPositionsAndVelocities(parrTest)
    runSimulation(n,eventQueue,particleArray,0)
    #plotPositionsAndVelocities(parrTest)
    print("Done")
    return particleArray

def simulationData(particleArray,n,interval):
    eventQueue=firstCollisions(particleArray)
    data=[]
    data.append(particleArray.copy())
    Time=0.0
    for i in range(int(n/interval)):
        Time = runSimulation(interval,eventQueue,particleArray, Time)
        data.append(particleArray.copy())
        print(i*interval/n,"%")
    return data

"""TESTING"""
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


