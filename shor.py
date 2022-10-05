import random
import math
from qiskit import *
from qiskit.circuit import *
from qiskit.extensions import *
from qiskit.circuit.library import *
from qiskit.extensions.simulator.snapshot import snapshot
from qiskit.quantum_info.operators import Operator
from qiskit.extensions.simulator.snapshot import snapshot
from qiskit.compiler import transpile
from qiskit.transpiler import PassManager
from qiskit.transpiler.passes import Unroller
from scipy import optimize
import matplotlib.pyplot as plt
import numpy as np
from qiskit import QuantumCircuit, Aer, transpile, assemble
from qiskit.visualization import plot_histogram
from math import gcd
from numpy.random import randint
import pandas as pd
from fractions import Fraction
from matplotlib.pyplot import plot,show

########################################
###           PARTIE OUTILS          ###
########################################

def GCD(x, y):#calcule le pgcd de deux nombres entiers positifs
   while(y):
       x, y = y, x % y
   return x

def gateMult(a,p,N,n):#construit les portes U
    nn = 2 ** n
    M = [[0 for x in range(nn)] for i in range(nn)]
    for x in range(2**n):
        res = x
        if x<N:
            res = a**p*x%N
        M[res][x] = 1
    U = Operator(M)
    return(UnitaryGate(U))

########################################
###         PARTIE QUANTIQUE         ###
########################################

def searchPeriod(a, N):
    size_eig = len(bin(N))-2 # nombre de qubits
    size_phi = 5 # qubits pour alilenter les U
    eig = QuantumRegister(size_eig, name="eig")
    phi = QuantumRegister(size_phi, name="phi")
    ceig = ClassicalRegister(size_eig, name="ceig")
    qc = QuantumCircuit(eig,phi,ceig)

    # Initialisation
    for i in range(0,size_eig):
        qc.initialize(0,i)
    qc.initialize(1,size_eig)
    for i in range(size_eig+1,size_eig+size_phi):
        qc.initialize(0,i)

    # H
    for i in range(0,size_eig):
        qc.h(eig[i])
        
    # U^2^n
    for i in range(0,size_eig):
        qc.append(gateMult(a,1,N,size_phi).power(2**(size_eig-i-1)).control(),[eig[size_eig-1-i],phi[0],phi[1],phi[2],phi[3],phi[4]])

    # QFT
    qc.append(QFT(size_eig).inverse(),eig)

    # mesures
    for i in range(0, size_eig):
        qc.measure(eig[i],ceig[i])

    backend = BasicAer.get_backend('qasm_simulator')
    job = execute(qc, backend, shots=1000)
    d = job.result().get_counts() #On récupère les résulats

    measured_phases = [] # on y conserve les valeurs des phases
    for output in d:
        decimal = int(output, 2)  # Conversion en entier base 10
        phase = decimal/(2**size_eig)  # Trouve la pahse correspondante
        measured_phases.append(phase)

    rows = []
    for phase in measured_phases:
        frac = Fraction(phase).limit_denominator(N) # phase = s / r, donc on approxime r
        rows.append(frac.denominator)
    
    r = max(rows, key = rows.count) #on récupère le r déviné le plus souvent
    return r

N = 15
r = 1

########################################
###         PARTIE CLASSIQUE         ###
########################################
while True:
    a = random.randint(2,N-1)
    b = 0
    pgcd = GCD(a,N)
    if pgcd != 1:# le pgcd est une solution triviale
        a = pgcd
        b = int(N/a)
        break
    else:
        r = searchPeriod(a,N)#partie quantique
        if not r % 2 == 1:
            if not math.pow(a,r/2)%N == -1 %N:
                a = GCD(math.pow(a,r/2)+1,N)
                b = GCD(math.pow(a,r/2)-1,N)
                if a * b == N and a!= N and b!=N:# sinon la periode est erronée
                   break
    
    
print("solutions : " + str(a) + " et " + str(b))
