### Imports ###

import numpy as np



### Exercice 2 ###

# 1)

GL = np.array([[0, 1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0],
               [0, 0, 0, 1, 0, 1, 0, 0, 0, 1, 0, 0],
               [0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0],
               [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1],
               [0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
               [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]], dtype=float)


# 2)

B = 263     # B = 263 W.K^(-1).m^(-2)
l = 140     # λ = 140 W.K^(-1).m^(-1)

def conductance(S, c=True, L=0, l=l, B=B) :
    return S * (B if c else (l / L if L > 0 else 0))

S10 = ((15 - 5) * 10**(-3)) * ((30 - 20) * 10**(-3))
GL[1, 9] = conductance(S10)     # GL_2_10 = 2,630 * 10^(-2) W.K^(-1)


# 3)

def distance(x1, y1, x2, y2) :
    return np.sqrt((x2 - x1)**2 + (y2 - y1)**2)

x3 = ((15 - 5) / 2) * 10**(-3)
y3 = ((60 - 30) / 2) * 10**(-3)
x4 = ((35 - 25) / 2) * 10**(-3)
y4 = ((70 - 20) / 2) * 10**(-3)
e = 6 * 10**(-3)
GL[2, 3] = conductance(((80-18)*(10**(-3)))*(e/2), False, distance(x3,y3,x4,y4))     # GL_3_4 = 2,604 W.K^(-1)


# 4)

GL[0, 3] = 0    # GL_1_4 non-défini


# 5)

S1 = ((18 - 0) * 10**(-3))**2
GL[0, 4] = conductance(S1, False, e)     # GL_1_5 = 7,560 W.K^(-1)


### Exercice 3 ###

# 1)

def cartToSpher(x, y, z) :
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    phi = np.arctan2(y / x)
    return r, theta, phi

def spherToCart(r, theta, phi) :
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z
