#!/usr/bin/env python

import numpy as np 
import matplotlib.pyplot as plt
import random as random

from numpy.random import rand


#initialized values
N = 25
latticeDimension = [N,N]
sweeps_no = 1000
T = np.arange(0.1,3.8,0.1)
beta = 1.0/T
kb = 1.0 # 1.38*10**-23
t_steps = 100 #number of points for temperature

z1 = sweeps_no*N*N #dividing factor

#Initialized energies, magnetization, specific heat and the magnetic susceptibility is an array of zeros 
Energy       = np.zeros(t_steps);   Magnetization  = np.zeros(t_steps)
SpecificHeat = np.zeros(t_steps);   Susceptibility = np.zeros(t_steps)


#creating the initial lattice with random spin orientations
def lattice(N):
    arr = np.random.rand(N,N)
    arr = np.where(arr <= 0.5, 1,-1)
    return arr
    
#function allowing us to determine the H value based on the spins of nearest neighbours
#def nearestneighbours(lattice, N, x, y):
#   return lattice(N)[(x - 1) % N, y] + lattice(N)[(x + 1) % N, y] + \
#       lattice(N)[x, (y - 1) % N] + lattice(N)[x, (y + 1) % N]

#converting the lattice function into a variable we can extract data from       
lat = lattice(N)
    
        
#main ising model     
J = 1    
t = 0
switch = 0
energies = 0
energycount = []
res = []
#fig = plt.figure()

def monte(state, beta):
    '''Monte Carlo move using Metropolis algorithm '''
    for k in range(N):
        for l in range(N):
                #defining random points
                i = np.random.randint(0, N)
                j = np.random.randint(0, N)
                site =  state[i, j] #spin state of point
                neighbours = state[(i - 1) % N, j] + state[(i + 1) % N, j] + \
       state[i, (j - 1) % N] + state[i, (j + 1) % N] #spin state of neighbours
                energychange = 2*site*neighbours
                #accepting the flip if the energy change is less than zero
                if energychange < 0: 
                    site *= -1
                #may also accept flip if the flipping probability is greater than a randomly generated number
                elif np.random.random() < np.exp(-energychange*beta):
                    site *= -1
                state[i, j] = site
    return state
    #fig.canvas.draw()
    #plt.pause(0.01)    
    #image = plt.imshow(lat, cmap='gray')
    #plt.show() 
    


def EnergyCal(state):
    '''Calculating the energy of the spin state'''
    Energy = 0
    for i in range(len(state)):
        for j in range(len(state)):
            site = state[i,j]
            neighbours = state[(i - 1) % N, j] + state[(i + 1) % N, j] + \
       state[i, (j - 1) % N] + state[i, (j + 1) % N] #spin state of neighbours
            Energy += -neighbours*site 
    return Energy/4. #prevents counting multiple times



def Magnetization(state):
    '''computes the magnetization of all sites in the lattice'''
    mag = np.sum(state)
    return mag


    
for m in range(len(T)):
    #initialized values for the energies and the magnetizations at zero       
    E1 = M1 = E2 = M2 = 0
    state = lattice(N) 
    newbeta = 1.0/T[m]
    
    
    for i in range(sweeps_no):         # equilibrate
        monte(state, newbeta)           # Monte Carlo moves

    for i in range(sweeps_no):
        monte(state, newbeta)           
        Ene = EnergyCal(state)     # calculate the energy
        Mag = Magnetization(state)        # calculate the magnetisation

        E1 = E1 + Ene
        M1 = M1 + Mag
        #M2 = M2 + Mag*Mag 
        E2 = E2 + Ene*Ene

        Energy[m]         = z1*E1
        #Magnetization[m]  = z1*M1
f = plt.figure()
plt.plot(T, Energy)
plt.show()
    
image = plt.imshow(lat, cmap='gray')
plt.show() 
