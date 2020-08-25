# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 18:11:08 2020

@author: ignacio
"""

import scipy as sp
from scipy.integrate import odeint

#Parametros
cm = 0.01
inch = 2.54*cm
g = 9.81


#Coeficientes de arrastre
p = 1.225
cd = 0.47
D = 8.5*inch
r = D/2
A = sp.pi*r**2
CD = 0.5*p*cd*A

#Masa
m = 15.

#Viento
Viento = [0, 10.0, 20.0]    


for V in Viento:
    
    def bala(z,t):
        zp = sp.zeros(4)
    
        zp[0] = z[2]
        zp[1] = z[3]
    
        v = z[2:4]    
        v[0] = v[0] - V
        v2 = sp.dot(v,v)
        vnorm = sp.sqrt(v2)
        FD = -CD * v2 * (v/ vnorm)
        zp[2] = FD[0] / m
        zp[3] = FD[1] / m - g
    
        return zp
    
    # vector de tiempo
    t = sp.linspace(0, 30, 1001)

    #Parte en el origen y tiene vx =

    vi = 100*1000./3600.
    z0 = sp.array([0, 0, vi, vi])

    sol = odeint(bala, z0, t)
    
    import matplotlib.pylab as plt

    x = sol[:,0]
    y = sol[:,1]
    name=(f"V = {V} m/s")    
    plt.figure(1)
    plt.plot(x,y,label=name)
    plt.ylabel("Y(m)")
    plt.xlabel("X(m)")
    plt.grid(True)
    plt.legend()


plt.xlim(0, 150)
plt.ylim(0, 50)
plt.savefig("Plott_P1_Entrega1")  
#plt.show()     
    



   
    

    

    
