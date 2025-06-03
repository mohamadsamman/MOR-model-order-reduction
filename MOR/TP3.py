### Imports ###

import numpy as np
import pandas as pd
import pykrige as pk



### Exercice 1 ###

# 1)

def produitVectoriel(u, v) :
    u = np.array(u)
    v = np.array(v)
    return np.cross(u, v)

def aireTriangle(A, B, C) :
    A = np.array(A)
    B = np.array(B)
    C = np.array(C)
    return np.linalg.norm(produitVectoriel(B-A, C-A)) / 2

def isInTriangle(A, B, C, X) :
    return np.isclose(aireTriangle(A,B,C), aireTriangle(X,B,C)+aireTriangle(A,X,C)+aireTriangle(A,B,X))

def estVraie(b) :
    print("Faux" if not(b) else "Vrai")

A  = (   0,    0, 0)
B  = (   0,    1, 0)
C  = (   1,    0, 0)
X1 = ( 0.5,  0.5, 0)
print("\nA =\n", A, "\nB =\n", B, "\nC =\n", C)
print("\nX1 =\n", X1)
estVraie(isInTriangle(A, B, C, X1))

X2 = (0.51, 0.51, 0)
print("\nX2 =\n", X2)
estVraie(isInTriangle(A, B, C, X2))


# 2)

with open("test_temperature_triangle.csv", encoding="utf-8-sig") as f :
    data1 = np.genfromtxt(f, delimiter=";", dtype=float, autostrip=True)    # 4 × n1 (températures stockées dans la ligne n°4, ce n'est pas très malin de séparer en deux...)
    n1 = np.shape(data1)[1]     # n1 = 9
    f.close()


# 3)

def distance(A, B) :
    A = np.array(A)
    B = np.array(B)
    return float(sum((B - A)**2)**0.5)

def getDataBis(X=X1, sort=True, data=data1) :
    dataBis = np.block([[                                             data],
                        [distance(p, X) for p in np.transpose(data[:3, :])],    # 1 × n1 (distances stockées dans la ligne n°5)
                        [                     np.arange(np.shape(data)[1])]])   # 1 × n1 (indices stockées dans la ligne n°6)
    return dataBis[:, np.argsort(dataBis[4,:])] if sort else dataBis            # tri par ordre croissant selon la distance entre le point et X1 (lorsque cela est souhaité)

dataBis1 = getDataBis()


# 4) et 5)

def getPointsTriangle(dataBis=dataBis1, X=X1) :
    pointsTriangle = np.transpose(dataBis[:, :3])                                                               # 3 × 6 (températures stockées dans la colonne n°4)
    pointsTriangle = np.transpose(np.block([[np.transpose(pointsTriangle)],
                                            [aireTriangle(X, pointsTriangle[1,:3], pointsTriangle[2,:3]),
                                             aireTriangle(pointsTriangle[0,:3], X, pointsTriangle[2,:3]),
                                             aireTriangle(pointsTriangle[0,:3], pointsTriangle[1,:3], X)]]))    # 3 × 7 (sous-aires stockées dans la colonne n°7)
    return pointsTriangle

pointsTriangle1 = getPointsTriangle()


# 6)

def tempPond(pointsTriangle=pointsTriangle1) :
    return sum(pointsTriangle[:, 3] * pointsTriangle[:, 6]) / sum(pointsTriangle[:, 6])

X1Temp = tempPond()
print("\nX1Temp =", X1Temp)


# 7)

def getAllDataBis(X=X2, sort=True, data=data1) :
    dataBis = getDataBis(X=X, sort=sort, data=data)
    pointsTriangle = getPointsTriangle(dataBis=dataBis, X=X)
    XTemp = tempPond(pointsTriangle=pointsTriangle)
    return dataBis, pointsTriangle, XTemp

dataBis2, pointsTriangle2, X2Temp = getAllDataBis()
print("\nX2Temp =", X2Temp)

X3 = (0.49, 0.49, 0)
print("\nX3 =\n", X3)
estVraie(isInTriangle(A, B, C, X3))
dataBis3, pointsTriangle3, X3Temp = getAllDataBis(X3)
print("\nX3Temp =", X3Temp)

print("\ndata1 =\n", data1, "\n "+str(np.shape(data1)), sep="")     # aucune critique pertinente...



### Exercice 2 ###

# 1)

data2 = pd.read_csv("donnees_temperature.csv", dtype=float, encoding="utf-8-sig")   # n × 3 (températures stockées dans la colonne n°3, c'est toujours aussi ballot de séparer en deux...)


# 2)

def distanceMatrix(p=np.array(data2.iloc[:,:2])) :
    n = np.shape(p)[0]
    return np.array([[distance(p[i,:], p[j,:]) if i <= j else 0 for j in range(n)] for i in range(n)]), n

M, n2 = distanceMatrix()    # n2 = 100


# 3)


