import numpy as np
import sympy
from mpmath import *
from matplotlib import pyplot as mp

def f(x,y):
    return y+x-x**2+1
def euler(xi,yi,h,m):
    u , v = [],[]
    for i in range(1,m):
        k1= h*f(xi,yi)
        k2 = h*f(xi+h,yi+k1)
        yi=yi+0.5*(k1+k2)
        xi= xi+h
        v=v+[yi]
        u=u+[xi]
    return [u,v]
def real(u):
    vreal = []
    fun=odefun(lambda x,y:y,0,1)
    for a in u:
        vreal = vreal+[fun(a)]
    return [vreal]

h= 0.1
m= 20
xi, yi = 0,1
[u,v]= euler(xi,yi,h,m)
[real]= real(u)
lx=np.linspace(0,3,20)
ly=np.linspace(0,20,20)
A,B =np.meshgrid(lx,ly)
U=1
mp.quiver(A,B,U,lx)
mp.xlabel("x")
mp.ylabel("f(x)")
mp.plot(u,v,'or', label = "euler")
mp.plot(u,real,'bs', label = "real")
mp.legend()
mp.show()
