#algorithme d'Euclide étendu
#pour 2 entiers a et b, calcule :
# - leur PGCD (r)
# - un couple de coéficients de Bézouts (u et v)
# Dans notre cas, a et b seront premiers entre eux, donc leur PGCD sera toujours 1, et u et v deviennent les inverses de a modulo b et de b modulo a, respectivement
def euclideEtendu(a, b) :
    r, u, v, rp, up, vp = a, 1, 0, b, 0, 1

    while rp != 0 :
        q = r // rp
        r, u, v, rp, up, vp = rp, up, vp, r-(q*rp), u-(q*up), v-(q*vp)
    return (r, u, v)

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

