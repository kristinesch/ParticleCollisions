from matplotlib import pyplot as plt
from problem1 import *


def HistogramMB(v,m):
    mb=MaxwellBoltzman(v,m)
    plt.hist(v,density=1,bins=11)
    plt.plot(v,mb)
    plt.show()

TestArray=simulation(500,1000)
TestV=getVelocities(TestArray)
HistogramMB(TestV,1)
print("M",m)




