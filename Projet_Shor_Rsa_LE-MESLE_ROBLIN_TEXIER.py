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

#algorithme d'Euclide étendu
#pour 2 entiers a et b, calcule :
# - leur PGCD (r)
# - un couple de coefficients de Bézouts (u et v)
# Dans notre cas, a et b seront premiers entre eux, donc leur PGCD sera toujours 1, et u et v deviennent les inverses de a modulo b et de b modulo a, respectivement
def euclideEtendu(a, b) :
    r, u, v, rp, up, vp = a, 1, 0, b, 0, 1

    while rp != 0 :
        q = r // rp
        r, u, v, rp, up, vp = rp, up, vp, r-(q*rp), u-(q*up), v-(q*vp)
    return (r, u, v)

#Génère la liste des nombres premiers jusqu'à n
def eratosthene(n):
    L = [ i for i in range(2,n+1) ]
    P = [ ]
    while len(L) != 0:
        P.append(L[0])
        i = L[0]
        for k in L:
            if k % i == 0:
                del(L[L.index(k)])
    return P

########################################
###            PARTIE RSA            ###
########################################

#La fonction de génération de clé prend en entrée 2 nombres premiers p et q, ainsi qu'un exposant e, lui aussi nombre premier
# on renvoit d, la clé privée, ainsi que n, la clé publique
def keysgen(p, q, e) :
    n = p*q
    phi = (p-1)*(q-1)
    c, d, dd = euclideEtendu(e, phi)
    return(d % phi, n)

#Pour chiffrer un message, on prend le message m, et la clé publique n et e
#On renvoit le message chiffré
def chiffrement(m, n, e) :
    return pow(m, e, n)

#Pour déchiffrer un message, on prend en entrée le message chiffré x, la clé privée d et n
#On renvoit le message déchiffré
def dechiffrement(x, n, d) :
    return pow(x, d, n)

#Génère 3 nombres premiers utilisables pour effectuer un chiffrement RSA
def genPrime():
    p = eratosthene(100000)
    e = p[-(int)(len(p)*0.2):]
    a = random.randint(0,len(e))
    b = a
    while b == a :
        b = random.randint(0,len(e))
    c = a
    while ( c == a or c == b ) or ((e[a]-1)*(e[b]-1) < e[c]) or (GCD((e[a]-1)*(e[b]-1) , e[c]) != 1) :
        c = random.randint(0,len(e))
    return e[a],e[b],e[c]

#Génère 3 nombres premiers utilisables par RSA mais assez petits
#pour être utilisables par Shor en un temps humain sur un ordinateur classique
#simulant du calcul quantique
def genLittlePrime():
    while True:
        try:
            e = eratosthene(15)
            a = random.randint(0,len(e))
            b = a
            while b == a :
                b = random.randint(0,len(e))
            c = a
            while ( c == a or c == b ) or ((e[a]-1)*(e[b]-1) < e[c]) or (GCD((e[a]-1)*(e[b]-1) , e[c]) != 1):
                c = random.randint(0,len(e))
            break
        except:
            e=[]
    return e[a],e[b],e[c]

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

########################################
###         PARTIE CLASSIQUE         ###
########################################
def solve(N):
    while True:
        a = random.randint(2,N-1)
        b = 0
        pgcd = GCD(a,N)
        if pgcd != 1:# le pgcd est une solution triviale
            a = pgcd
            b = int(N/a)
            print("Solutions triviales trouvées")
            break
        else:
            r = searchPeriod(a,N)#partie quantique
            if not r % 2 == 1:
                if not math.pow(a,r/2)%N == -1 %N:
                    a = GCD(math.pow(a,r/2)+1,N)
                    b = GCD(math.pow(a,r/2)-1,N)
                    if a * b == N and a!= N and b!=N:# sinon la periode est erronée
                       print("Solutions trouvées par calcul de phase")
                       break
    return a,b

########################################
###           PARTIE MAIN            ###
########################################

#création des clés
N=100
while N > 25:
    p,q,e = genLittlePrime()
    priv, pub = keysgen(p,q,e)
    N = p * q    
print("Premiers générés : ",p,q,e)
print("Clé publique : ",pub)
print("Clé privée : ", priv)

print()
#test du codage RSA
message = random.randint(2,5)
messageC = chiffrement(message, pub, e)
messageD = dechiffrement(messageC, pub, priv)
print("Message original : ",message)
print("Message codé : ", messageC)
print("Message décodé : ",messageD)

print()
#cassage avec Shor
a,b = solve(pub)
print("Solutions : " + str(a) + " et " + str(b))
privCasse, pub = keysgen(a,b,e)
print("Clé privée cassée : ",privCasse)
print("Message décodé avec clé cassée : ", dechiffrement(messageC, pub, privCasse))

