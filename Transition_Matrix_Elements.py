#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import matplotlib.pyplot as plt
import numpy as np
from math import pi
from qiskit import *
from qiskit import QuantumCircuit,execute,Aer,IBMQ
from qiskit.compiler import transpile,assemble
from matplotlib import style
from qiskit.circuit.library import QFT
from qiskit.tools.monitor import job_monitor
from itertools import repeat
from qiskit import IBMQ, Aer, transpile, assemble,execute
from qiskit import ClassicalRegister, QuantumRegister
from qiskit.visualization import plot_histogram, plot_bloch_multivector
from qiskit.compiler import transpile,assemble
import operator
import time 
import random
from numpy.random import seed
from numpy.random import randint 
import numpy.linalg as linalg
import copy


# ##=================== KE function ===================##
# 
# For 3 qubits [q0, q1, q2], we first use q0 & q1 simulantously to undergo the kinetic circuit. 
# 
# Then, we use q1 & q2. The pattern goes on for bigger systems.
# 
# Typically,
# 
# 1) q0 is a control for X and q1 is a traget
# 
# 2) Rz(- pi /2) on q0
# 
# 3) q1 is a control for X and q0 is a traget
# 
# 4) Ry(- hopping angle) on q0
# 
# 5) q1 is a control for X and q0 is a traget
# 
# 6) Ry(+ hopping angle) on q0
# 
# 7) Rz(+ pi/2) on q0
# 
# 8) q0 is a control for X and q1 is a traget
# 
# 
# This makes sense because the kinetic term is just hopping from one location to the next.
# So, in general for n qubits, we will connect the i-th qubit with the i-th + 1 qubit..
# 
# It would be function of qc, hopping angle (hopping parameter h*dt), and the number of qubits.

# In[ ]:


#qc stand for quatum computer but it's actually just a neat way for psi or quantum state
#it's the key variable in any system' it's the block upon which this system is applied
def Kinetic_Energy(qc, hopping_angle, num_qubits):
    for i in range(num_qubits - 1):
        qc.cx(i, i + 1)
        qc.rz(-np.pi/2.0, i)
        qc.cx(i + 1, i)
        qc.ry(-hopping_angle, i)
        qc.cx(i + 1, i)
        qc.ry(hopping_angle, i)
        qc.rz(np.pi/2.0, i)
        qc.cx(i, i + 1)

    
# psi = QuantumCircuit(3)
# Kinetic_Energy(psi, 0.002, 3)
# psi.draw(output="mpl")


# ##=================== PE function ===================##
# 
# Acoording to how far the ith site is from origin, the potential acts ~ x^2 
# This fact translates to ((L - 1 - 2 i)/2)^2 being the elements on the diagonal
# 'L' is the total number of sites while 'i' is the site number.
# It would be function of qc, potential_angle = -l*dt (where l = m w^2 /2) and the number of qubits.
# 

# In[ ]:


def Potential_Energy(qc, potential_angle, num_qubits): 
    for i in range(num_qubits):
        if i == int ((num_qubits-1)/2):
            continue
        r = (i-((num_qubits-1)/2))
        qc.p(potential_angle*r**2,i)

    
# # sanity check
# psi=QuantumCircuit(3)
# Potential_Energy(psi, -0.001, 3)
# psi.draw(output="mpl")


# Now, we aready to write a full time-evolution Hamiltonain quantum circuit
# This Hamiltonian is of one particle in a simple hamrnoic oscillator potential
# The time evoultion of any Hamiltonian is a unitary opertion 
# We are using frist order Lee_Torotter-Suzuki formula 

# In[ ]:


##=================== Complete SHO N times ===================##

# 'time_steps' defines the number of time steps 'dt' to be taken 
# so that the system evolves with total time T = time_steps * dt 

def Uintary_Time_Evolution_SHO(qc, hopping_angle, potential_angle, time_steps, num_qubits):
    for i in range(time_steps):
        Kinetic_Energy(qc, hopping_angle*0.5, num_qubits)
        Potential_Energy(qc, potential_angle, num_qubits)
        Kinetic_Energy(qc, hopping_angle*0.5, num_qubits)


# # Sanity check. This is a deep circuit
# psi = QuantumCircuit(3)
# Uintary_Time_Evolution_SHO(psi, 0.002, -0.001, 10, 3)
# psi.draw(output="mpl")


# Now, let's add electromagnetic (E1) interaction to this system.
# Here, we are building the Unitary Time Evolution circuit for the E1 transition.
# The idea is to encode the entire electromagnetic field into one qubit, 
# which we refer to as an 'auxiliary qubit.' 
# Why only one? Because it has a narrow spectrum for interaction,
# and the interaction occurs within a cavity, 
# so we are only interested in photon modes of 1 (indicating the presence of a photon) 
# and photon modes of 0 (indicating no photon). It only takes one qubit!
# 
# We are going to sandwich the Unitary Time Evolution (U) circuit of the E1 transition between two X-gates.
# Why? Because using two X-gates transforms our circuit from its current form into the required form. Please check Nillsen&Chung (2nd Ed. Ch 4 section 3) for details. 
# Additionally, we need to apply U once and then U† (the conjugate transpose of U) afterward since everything must remain unitary

# In[ ]:


#### ====================================== the E1 circuit  ====================================== ###

# beta, the E1 parameter or photon term (beta = alpha*R*dt = gamma/2)

def E1_circuit(qc, beta, num_qubits):
    for i in range(num_qubits-1):
        R = (i-((num_qubits-2)/2))  #R is the distance [L-(L-1)/2]
        qc.cx(num_qubits-1,i)     
        qc.cx(i, num_qubits-1)     
        #this is the start of the Uc circuit, an X gate.. and the Z-angle could also be pi/2
        qc.rz(-np.pi, i)   
        qc.cx(num_qubits-1, i)
        qc.ry(beta*R, i)
        qc.cx(num_qubits-1, i)
        qc.ry(-beta*R, i)
        qc.rz(np.pi, i)  #could be -pi/2
        qc.cx(i, num_qubits-1)
        #here is the end of the Uc circuit.
        qc.cx(num_qubits-1, i)



# # Do U Want To See This? 
# psi=QuantumCircuit(3+1) #you gotta add 1 becuase the whole thing is controled by the photon which is the nth+1 qubit
# E1_circuit(psi, 0.002, 3)
# psi.draw(output="mpl")


# Now, this following function applies the Unitary Time Evolution full Hamiltonian (Kinetic + Potential + E1)
# Let's call it Unitary_SHO_E1 

# In[ ]:


##==================## THE FULL SHO + E1 FUNCTION ##==================##

def Unitary_SHO_E1(qc, hopping_angle, potential_angle, time_steps, beta, num_qubits):
    for i in range(time_steps):
        Kinetic_Energy(qc, hopping_angle*0.5, num_qubits-1)
        Potential_Energy(qc, potential_angle*0.5, num_qubits-1)
        #photon energy term (E is the photon kinetic energy)
        qc.rz(-E*0.5, num_qubits-1) 
        #photon E1 term 
        E1_circuit(qc, beta, num_qubits)
        qc.rz(-E*0.5, num_qubits-1)
        Potential_Energy(qc, potential_angle*0.5, num_qubits-1)
        Kinetic_Energy(qc, hopping_angle*0.5, num_qubits-1)



# In[ ]:


#Sanity Check; let's chose some physical parameters to see this system in action
a = 1.0/197.33   #lattice spacing = 1 fermi
w = 8*a          #omega, the harmonic potential frequancy = 8 MeV = 8*a in lattice units
dt = 0.01/w      #time step
E = 1*w*dt       #epsilon, the photon kinetic energy

# psi = QuantumCircuit(3)
# Unitary_SHO_E1(psi, 0.002, -0.001, 10, 0.002, 2)
# psi.draw(output="mpl") 


# Now, let's design a function that takes care of plugging in the initial state 
# as well as elvolving this state under this Unitary_SHO_E1 circuit
# the initial state is a actually a tensor product of the one of the SHO states & one of the photon states
# for exmaple we could start from the ground state of the SHO particle with no photon
# so it'd be |ground ⊗ 0 > ... where |0> represents the state on zero photon
# therefore this circuit takes literally two states
# meaning we have to initialize two states 
# one state for the SHO particle itself (it'd be a 2^n component state)
# and another state for the photon mode and it's either (0, 1) or (1, 0) repressenting no phton or a photon respectively
# this is why this circuit has the two labels psi (for particle state) & phi (photon state)
# let's call it two level circuit

# In[ ]:


## ==================  SHO + E1 Circuit ============== ##
#This circuit does the E1 transition on a particle in SHO potential

def two_level_circuit(psi, phi, hopping_angle, potential_angle, time_steps, beta, num_qubits):
    
    #circuit of n quantum bits and n classical bits
    qc = QuantumCircuit(num_qubits, num_qubits)
    
    #append psi on the first n-1 bits
    qc.append(psi, list(range(0, num_qubits-1))) 
    
    #append a photon on the nth + 1 qubit
    qc.append(phi, list(range(num_qubits-1, num_qubits)))
    
    #apply U(H = KE + PE + E1)
    Unitary_SHO_E1(qc, hopping_angle, potential_angle, time_steps, beta, num_qubits) 
    
    return qc


# In[ ]:


## ========== Let's see how this circuit looks like ========== ##

# ================= Defining the Variables=======================
a = 1.6/197.33   #lattice spacing

w = 8*a    #omega, the harmonic potential frequancy

dt = 0.01/w   # time step

m_lat = 940*a   # particle mass in lattice units

L = 7   #Number of lattice sites (particle sites)

num_qubits = L + 1 

h = 0.5/m_lat     # hopping parameter in lattice units

l = 0.5*m_lat*w**2    #the SHO potential energy factor (1/2 m w^2)


##Angles
hopping_angle = h*dt    #kientic energy angle

potential_angle = -l*dt    #potntial energy angle

E = 1*w*dt     #epsilon, the photon kinetic energy

q = 1 #charge

sigma = 11*a

e = np.sqrt(4*np.pi/137)

beta = -dt*q*e*sigma**2  #the photon rotation angle

w1 = 2*q*e*sigma**2/(np.sqrt(2*m_lat*w))

T = 0.65*np.pi/w1   #total time

Total_time = 5*np.pi/w

time_steps = int(Total_time/dt)


# ##=================== Testing the SHO + E1 Circuit ===================##
# Now, let's get some rabi oscillations 
# for now I will provide the exact eigen state by hand 
# However, there is a code that computes the specturm of the lattice Hamiltonian evolution 

# In[ ]:


# This is the ground state on the SHO particle
# here there should be a function the provides the exact ground state of the system!
# This the photon state of zero (no photon)
photon_state = [0, 1]

#let's create circuits for these states
psi = QuantumCircuit(num_qubits-1)
phi = QuantumCircuit(1)

## the initialized states: 
psi.initialize(ground_state)
phi.initialize(photon_state)


# In[ ]:


#make a list of zeros and a one anywhere else .. ex: |000100> or |000001>  
key = '0' * (num_qubits-2)   
keys = [0]*(num_qubits-1) 

for i in range(num_qubits-1): 
   kkeys = "0" + key[:i] + "1" + key[i:] #inserting 0 at the 1st bit and 1 anywhere else, once at a time
   keys[i] = kkeys  #the list of keys


# In[ ]:


#the number of trails
trials = 5000

#let's run the whole ting many times
iterations = 5000    #the number of times we run the whole thing
rabi = np.zeros((iterations+1, 2))   #2D array where time = i *steps  & it is (f+1)x2 array


# In[ ]:


def rabi_oscillations(num_qubits, keys, iterations):
    for i in range(1, iterations+1):
        Total = np.pi/w1*i
        time_steps = int(Total/dt)
        rabi[i, 0] = dt*time_steps*w
        result =  two_level_circuit(psi, phi, hopping_angle, potential_angle, time_steps, beta, num_qubits)
        
        temp = copy.deepcopy(result)
        for k in range(num_qubits):
            temp.measure(k, k)
        simulator = Aer.get_backend('aer_simulator')
        
        job = execute(temp, simulator, shots = trials, seed = 3555)
        result = job.result()
        counts = result.get_counts(temp)
        
        for j in range(num_qubits-1):
            if counts.get(keys[j]): #this should be the defined list of keys
                rabi[i, 1] += counts.get(keys[j])/trials
        return rabi
    
        #the probability of absorbation is stored run-by-run in the 2D array 


# In[ ]:


# Run it!
start_time = time.time()

rabi = rabi_oscillations(num_qubits, keys, iterations)

end_time = time.time()
execution_time = end_time - start_time
print(f"Execution Time: {execution_time:.6f} seconds")
    
##### Probability plot #####
# lis = list(range(0, f))
# x = np.array(lis)

plt.rcParams["figure.figsize"] = (8,6)


# plot the functions
plt.scatter(rabi[:, 0], rabi[:, 1], marker='o', c = 'red', label='Probability Curve')


plt.xlabel('wt')
plt.ylabel('Probability')
plt.title('Probability Oscillations as Sin^2 funstion of time')
plt.grid(True)
plt.legend()
plt.show()


# In[ ]:




