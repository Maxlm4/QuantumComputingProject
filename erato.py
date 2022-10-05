import random

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

p = eratosthene(1000)
e = p[-50:]

a = random.randint(0,len(e))
b = a
while b == a :
    b = random.randint(0,len(e))
c = a
while ( c == a or c == b ) and ((e[a]-1)*(e[b]-1) > e[c]):
    c = random.randint(0,len(e))


print(e[a])
print(e[b])
print((e[a]-1)*(e[b]-1))
print(e[c])
