import numpy as np
from sympy import integrate
from sympy.abc import a, x, y
from sympy import Symbol
from sympy import diff
import math
import matplotlib.pyplot as plt

def Real(x):
    return x/2+x**2-(11/20)*math.exp(2*x)+7/4
def Euler(f1,n,h):
  x = Symbol('x')
  y = Symbol('y')

  derivada1x = diff(2*y-(2*x**2)+x-3,x)
  derivada2x = diff(derivada1x,x)
  derivadaxy = diff(derivada1x,y)

  derivada1y = diff(2*y-(2*x**2)+x-3,y)
  derivada2y = diff(derivada1x,y)

  realL=[]
  eulerL=[]
  heunL=[]
  kuttaL=[]
  taylor=[]
  puntos=[]
  x=0
  y1,y2,y3,y4=1.2,1.2,1.2,1.2

  i= 0
  aux=y1
  yreal=Real(x)
  realL.append(yreal)
  eulerL.append(y1)
  heunL.append(y2)
  kuttaL.append(y3)
  taylor.append(y4)
  puntos.append(i);
  y=y1
  error= float(((yreal-y1)/yreal) *100)
  error2=float(((yreal-y2)/yreal) *100)
  error3 = float(((yreal-y3)/yreal) *100)
  error4 = float(((yreal-y4)/yreal) *100)
  print(("| {:^15} _ {:^15} _ {:^15} _ {:^15} _ {:^15} _ {:^15} _ {:^15} _ {:^15} _ {:^15} _ {:^15} |").format("---------------","---------------","---------------","---------------","---------------","---------------","---------------","---------------","---------------","---------------"))
  print(("| {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^15} | {:^15} |").format("i","y-real","Euler","Heun","Runge-K","Taylor","errorEu(%)","errorHeun(%)","errorRunge(%)","errorTay(%)"))
  print(("| {:^15} + {:^15} + {:^15} + {:^15} + {:^15} + {:^15} + {:^15} + {:^15} + {:^15} + {:^15} |").format("---------------","---------------","---------------","---------------","---------------","---------------","---------------","---------------","---------------","---------------"))
  print(("| {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} |").format(i,yreal,y1,y2,y3,y4,error,error2,error3,error4))
  i=h
  while i<n:
      puntos.append(i);
      y1=aux+ ( h * f1((i-h),aux) )  #EULER
      y2 = y2 + ((f1((i-h),y2) + f1(i+h/1.1,y1))/2) * h #HEUN
      k1 = f1(i-h , y3)
      k2 = f1((i-h)+h/2 , y3 + (k1/2)*h)
      k3 = f1((i-h)+h/2 , y3 + (k2/2)*h)
      k4 = f1(i,y3 + k3*h)
      y3 = y3 + (1/6)*(k1+2*k2+2*k3+k4)*h #RUNGE KUTTA
      x=i-h
      y = y4
      #y4 = f1(i-h,y4) + eval(str(derivada1x))*h + eval(str(derivada1y))*0+(1/2)*((eval(str(derivada2y))*h*h)+(2*eval(str(derivadaxy))*h)+(eval(str(derivada2y))))
      y4 = y4 + ((f1(x,y4) + eval(str(derivada1x)) * (h/2))*h)

      yreal=Real(x)
      aux=y1
      error=((yreal-y1)/yreal) *100
      error2=((yreal-y2)/yreal) *100
      error3=(((yreal-y3)/yreal) *100)
      error4 = (((yreal-y4)/yreal) *100)

      realL.append(yreal)
      eulerL.append(y1)
      heunL.append(y2)
      kuttaL.append(y3)
      taylor.append(y4)
      print(("| {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} | {:^15.2f} |").format(i,yreal,y1,y2,y3,y4,error,error2,error3,error4))
      i+=h
  print(("| {:^15} + {:^15} + {:^15} + {:^15} + {:^15} + {:^15} + {:^15} + {:^15} + {:^15} + {:^15} |").format("---------------","---------------","---------------","---------------","---------------","---------------","---------------","---------------","---------------","---------------"))
  plt.plot(puntos,realL,'r--',label='Real')
  plt.plot(puntos,eulerL,label='Euler')
  plt.plot(puntos,heunL,label='Heun')
  plt.plot(puntos,kuttaL,label='Runge-Kutta')
  plt.plot(puntos,taylor,label='Taylor')
  plt.legend()
  plt.grid();
  plt.title("Comparacion para 2*y-(2*x**2)+x-3 , con %d puntos y un paso de %.1f"%(n,h))
  plt.show()
  return y1



x = Symbol('x')
y = Symbol('y')
f2= integrate(2*y-(2*x**2)+x-3,x)
f2=str(f2)+'+1.2'
x = 0
y = 1.2
#print ("--> ",eval(f2))
Euler(lambda x,y: 2*y-(2*x**2)+x-3 ,10,0.1)
print()
#print(Euler(lambda x:  12*x**2 -2*x**3 - 20*x+ 8.5,10,0.5,f2))
