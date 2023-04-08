#author Gerbaud FLorent
#Algorithme LU
#03/03/2023

import numpy as np
import math
import matplotlib.pyplot as plt
from copy import deepcopy

def Dichotomie(f,a,b,eps,Nmax):
    
    residu=[]
    if(f(a)*f(b)>=0):
        print("La méthode n'est pas applicable sur cette fonction ou sur cet intervalle")
    k=0
    c=(a+b)/2
    fc=f(c)
    while(abs(fc)>eps and k<Nmax):
        residu.append(abs(fc))
        k=k+1
        if(f(a)*f(c)<0):
            b=c
        elif(f(c)*f(b)<0) :
            a=c
        else:
            print("la méthode n'est plus applicable")
        c=(a+b)/2
        fc=f(c)
    return (c,residu,k)

def PointFixe(F,x0,eps,Nmax):
    residu=[]
    x=x0
    k=0
    fx=F(x)-x
    while(abs(fx)>eps and k<Nmax):
        residu.append(abs(fx))
        x=F(x) #x_{n+1}=F(x_n)
        fx=F(x)-x #approche sol
        k=k+1
    return x,k,residu

def Newton(fDf,x0,eps,Nmax):
    return PointFixe(fDf,x0,eps,Nmax)
        
def f(x):
    return np.power(x,2)-2
def Df(x):
    return 2.0*x
def fDf(x):
    return x-(f(x)/Df(x))

def F(x):
    return (2.0*x+2.0/x)/3.0


eps=1e-6
x0=1/2
print("x0: ",x0)


#_____________________________________________ Test Dichotomie _____________________________________________#
(c,r,k)=Dichotomie(f,-7,5,eps,300)
print(c,r[-1],k)
#_____________________________________________ Test PointFixe _____________________________________________#
(x,k,r0)=PointFixe(F,x0,1e-6,300)
print(x,k)
#_____________________________________________ Test Newton _____________________________________________#
(x1,k1,r1)=Newton(fDf,x0,1e-6,30000)
print(x1,k1)

# Partie Satelite


#_____________________________________________ Newton in N dim _____________________________________________#

def NewtonDimN(fonc, jac, x0, eps, Nmax,a,b,c,d,h):
    x = x0
    k = 0
    fx = fonc(x,x0,a,b,c,d,h)
    dfx = jac(x,a,b,c,d,h)
    while(np.linalg.norm(fx) > eps and k < Nmax):
        k=k+1 #incrément de k
        x=x-np.dot(np.linalg.inv(dfx),fx) # calcul du nouveau x
        fx=fonc(x,x0,a,b,c,d,h)
        dfx=jac(x,a,b,c,d,h) # mise ajour de la jacobienne
    return x

