from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

def f(x,y,z):
    return x+y+z
def g(x,y,z):
    return -x+y-z
def sistemas(xi,yi,zi,h,m):
    u,v,w = [],[],[]
    for i in range(1,m):
        k1y=h*f(xi,yi,zi)
        k1z=h*g(xi,yi,zi)
        k2y=h*f(xi+h,yi+k1y,zi+k1z)
        k2z=h*g(xi+h,yi+k1y,zi+k1z)
        yi=yi+0.5*(k1y+k2y)
        v=v+[yi]
        zi=zi+0.5*(k1z+k2z)
        w=w+[zi]
        xi= xi+h
        u=u+[xi]
    return [u,v,w]

[x,y,z]=sistemas(0,1,2,0.001,10000)
[x1,y1,z1]=sistemas(0,1,2,0.1,100)
figura_3d = plt.figure()
ax = figura_3d.gca(projection = '3d')
ax.view_init(45, -35)
plt.plot(x,y,z,'or', label = "h= 0,001 m= 10000")
plt.plot(x1,y1,z1,'bs', label = "h= 0,1 m= 100")
plt.legend()
plt.show()
