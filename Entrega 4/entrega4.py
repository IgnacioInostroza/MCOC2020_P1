# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 08:04:43 2020

@author: ignacio
"""


from scipy.integrate import odeint
import numpy as np
import math
import matplotlib.pylab as plt


m = 1 #kg
f = 1.0 #Hz
e = 0.2
w = 2*(math.pi)*f
k = m*(w**2)
c = 2*e*w*m

def eulerint(zp, z0, t, Nsubdivisiones):
    Nt= len(t)
    Ndim = len(z0)
    
    z = np.zeros((Nt, Ndim)) 
    z[0, :] = z0
    
    for i in range (1, Nt):
        t_anterior = t[i-1]
        
        dt = (t[i] - t[i-1])/Nsubdivisiones
        
        z_temp = z[i-1, :]*1.0
        
        for k in range(Nsubdivisiones):
            z_temp += dt * zp(z_temp, t_anterior + k*dt)
        z[i, :] = z_temp
        
    return z    

# mx" + cx' + kx = 0        
# y' = -(k/m)*x - (c/m)*y
# x' = y

def zp(z, t):
    zp = np.zeros(2)
    
    zp[0] = z[1]
    zp[1] = -(k/m)*z[0] - (c/m)*z[1]
        
    return zp

 

z0 = np.array([0,1]) #vector inicial
t = np.linspace(0, 4., 100)

y = c/(2*m)
I_analitica = np.sqrt(z0[0]**2 + ((z0[1] + y*z0[0])/w)**2) 
z_real = I_analitica * np.exp(-y*t) * np.sin(w*t)

sol = odeint(zp, z0, t)
z_odeint = sol[:, 0]

sol = eulerint(zp, z0, t, 1)
z_euler_NSD1 = sol[:,0]

sol = eulerint(zp, z0, t, 10)
z_euler_NSD10 = sol[:,0]

sol = eulerint(zp, z0, t, 100)
z_euler_NSD100 = sol[:,0]


plt.figure()
plt.plot(t,z_real, color = "black", linewidth = 2, label = "real")
plt.plot(t,z_odeint, "b", label="odeint")
plt.plot(t,z_euler_NSD1, ":", color = "green", label="eulerint_ NSD = 1")
plt.plot(t,z_euler_NSD10, ":", color = "red", label="eulerint_NSD = 10")
plt.plot(t,z_euler_NSD10,":", color = "orange", label="eulerint_NSD = 100")

plt.xlabel("t(s)")
plt.ylabel("x(t) (m)")


plt.grid(True)
plt.legend()
plt.savefig("Plott_Entrega4")
plt.show



                     