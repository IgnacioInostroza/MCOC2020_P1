# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 19:53:37 2020

@author: ignacio
"""


from scipy.integrate import odeint
import math
import numpy as np

#Datos
hr=3600
km=1e3
G = 6.67408e-11 #m/s**2
m = 5.972e24 #kg
radio =6371*km #km
omega=-7.2921150e-5
HB=700*km

Fgmax = G*m/radio**2


zp = np.zeros(6)    
def satelite(z,t):
    c = math.cos(omega*t)
    s = math.sin(omega*t)
    R = np.array([[c, s, 0],[-s, c, 0],[0, 0, 1]]) 
    R1 = omega * np.array([[-s, c, 0],[-c, -s, 0],[0, 0, 0]]) 
    R2 = (omega**2)* np.array([[-c, -s, 0],[s, -c, 0],[0, 0, 0]])
    
    z1 = z[0:3] 
    z2 = z[3:6]
    
    r2 = np.dot(z1,z1)
    r = np.sqrt(r2)
    
    Fg=(-G*m/r**2)*(R@(z1/r))
    
    zp[0:3] = z2
    zp[3:6] = R.T@(Fg - (2*(R1@z2) + (R2@z1)))    
        
    return zp

from datetime import datetime
ti = "2020-07-29T22:59:42.000000"
ti = ti.split("T")
ti = "{} {}".format(ti[0], ti[1])
ti = datetime.strptime(ti, '%Y-%m-%d %H:%M:%S.%f')

tf = "2020-07-31T00:59:42.000000"
tf = tf.split("T")
tf = "{} {}".format(tf[0], tf[1])
tf = datetime.strptime(tf, '%Y-%m-%d %H:%M:%S.%f')

deltaT = (tf - ti).seconds


x_i = 1665352.492338
y_i = -6756092.122778
z_i = 1281662.400464

vx_i = -1871.276354
vy_i = 934.835139
vz_i = 7304.321490


x_f = 18066.885210
y_f = -3451953.355853
z_f = -6183737.173452

vx_f = -2407.237376
vy_f = -6269.867602
vz_f = 3494.919702

t = np.linspace(0, deltaT, 9361)

z0 = np.array([x_i, y_i, z_i, vx_i, vy_i, vz_i])

sol = odeint(satelite, z0, t)

x = sol[:, :]

pos_final = np.array([x_f, y_f, z_f, vx_f, vy_f, vz_f]) - sol[-1]
   
deriva = np.sqrt(((x_f) - sol[-1,0])**2 + ((y_f) - sol[-1, 1])**2 + ((z_f) - sol[-1, 2])**2)

print (deriva)



  


