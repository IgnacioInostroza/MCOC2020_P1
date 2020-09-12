# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 19:53:37 2020

@author: ignacio
"""


from scipy.integrate import odeint
import math
import numpy as np
import datetime as dt
#from leer_eof import leer_eof
from time import perf_counter as cron
import xml
import xml.etree.ElementTree as ET
from numpy import zeros
from sys import argv



def utc2time(utc, ut1, EOF_datetime_format = "%Y-%m-%dT%H:%M:%S.%f"):
	t1 = dt.datetime.strptime(ut1,EOF_datetime_format)
	t2 = dt.datetime.strptime(utc,EOF_datetime_format)
	return (t2 - t1).total_seconds()


def leer_eof(fname):
	tree = ET.parse(fname)
	root = tree.getroot()

	Data_Block = root.find("Data_Block")		
	List_of_OSVs = Data_Block.find("List_of_OSVs")

	count = int(List_of_OSVs.attrib["count"])

	t = zeros(count)
	x = zeros(count)
	y = zeros(count)
	z = zeros(count)
	vx = zeros(count)
	vy = zeros(count)
	vz = zeros(count)

	set_ut1 = False
	for i, osv in enumerate(List_of_OSVs):
		UTC = osv.find("UTC").text[4:]
		
		x[i] = osv.find("X").text   #conversion de string a double es implicita
		y[i] = osv.find("Y").text
		z[i] = osv.find("Z").text
		vx[i] = osv.find("VX").text
		vy[i] = osv.find("VY").text
		vz[i] = osv.find("VZ").text

		if not set_ut1:
			ut1 = UTC
			set_ut1 = True

		t[i] = utc2time(UTC, ut1)

	return t, x, y, z, vx, vy, vz


tti = cron()

archivo = argv[1]
#archivo = "S1A_OPER_AUX_POEORB_OPOD_20200819T121202_V20200729T225942_20200731T005942.EOF"
at, ax, ay, az, avx, avy, avz = leer_eof(archivo)

Pred = archivo.replace(".EOF",".PRED")

#Datos
hr=3600
km=1000.
km3=(1000.)**3
km5=(1000.)**5
km6=(1000.)**6
v = 398600441500000.000000 # m**3/s**2 : G*(masa de la tierra)
G = 6.67408e-11 #m/s**2
m = 5.972e24 #kg
radio = 6378136.3  
omega=-7.2921150e-5
HB=700.*km



Fgmax = v/radio**2

    
def satelite(z,t):
    c = math.cos(omega*t)
    s = math.sin(omega*t)
    Rp = np.array([[c, s, 0],[-s, c, 0],[0, 0, 1]])
    Rp_T = np.array([[c, -s, 0],[s, c, 0],[0, 0, 1]]) 
    R1 = omega * np.array([[-s, c, 0],[-c, -s, 0],[0, 0, 0]]) 
    R2 = (omega**2)* np.array([[-c, -s, 0],[s, -c, 0],[0, 0, 0]])
    
    x = z[0:3] 
    xp = z[3:6]
    
    r2 = np.dot(x,x)
    r = np.sqrt(r2)
    
    xs = Rp@x
    rnorm = xs / r
    Fg = -v/r**2 * rnorm
    
    #z2 = xs[2]**2
    #rflat = xs[0]**2 + xs[1]**2
    from sympy import sqrt
    
   
    
    Fx = -5.26658496e+25*xs[0]*xs[2]**2*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.5) - 3.0*xs[0]*(2.63329248e+25*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 8.7776416e+24)*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-2.5) - 4.0*xs[0]*(-6.54783215e+29*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) + 3.92869929e+29*xs[2]/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2))*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.0) - 9.0*xs[0]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-5.5)*(-1.11145500984375e+64*xs[2]**8/(xs[0]**2 + xs[1]**2 + xs[2]**2)**4 + 2.074716018375e+64*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 1.1969515490625e+64*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 2.17627554375e+63*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 6.04520984375e+61) - 8.0*xs[0]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-5.0)*(-1.6124802676125e+57*xs[2]**7/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(7/2) + 2.6047758169125e+57*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) - 1.1839890076875e+57*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) + 1.315543341875e+56*xs[2]/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2)) - 7.0*xs[0]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.5)*(2.0907567883125e+50*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 2.8510319840625e+50*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 9.503439946875e+49*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 4.52544759375e+48) - 6.0*xs[0]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.0)*(-7.544857768875e+42*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) + 8.38317529875e+42*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) - 1.796394706875e+42*xs[2]/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2)) - 5.0*xs[0]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.5)*(-4.67333269375e+36*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 4.0057137375e+36*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 4.0057137375e+35) + (1.964349645e+30*xs[0]*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) - 3.92869929e+29*xs[0]*xs[2]/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2))*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-2.0) + (1.8693330775e+37*xs[0]*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 8.011427475e+36*xs[0]*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2)*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-2.5) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.5)*(8.89164007875e+64*xs[0]*xs[2]**8/(xs[0]**2 + xs[1]**2 + xs[2]**2)**5 - 1.244829611025e+65*xs[0]*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**4 + 4.78780619625e+64*xs[0]*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 4.3525510875e+63*xs[0]*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.0)*(1.12873618732875e+58*xs[0]*xs[2]**7/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(9/2) - 1.30238790845625e+58*xs[0]*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(7/2) + 3.5519670230625e+57*xs[0]*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) - 1.315543341875e+56*xs[0]*xs[2]/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2)) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.5)*(-1.2544540729875e+51*xs[0]*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**4 + 1.140412793625e+51*xs[0]*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 1.900687989375e+50*xs[0]*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.0)*(3.7724288844375e+43*xs[0]*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(7/2) - 2.514952589625e+43*xs[0]*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) + 1.796394706875e+42*xs[0]*xs[2]/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2))
    Fy = -5.26658496e+25*xs[1]*xs[2]**2*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.5) - 3.0*xs[1]*(2.63329248e+25*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 8.7776416e+24)*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-2.5) - 4.0*xs[1]*(-6.54783215e+29*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) + 3.92869929e+29*xs[2]/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2))*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.0) - 9.0*xs[1]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-5.5)*(-1.11145500984375e+64*xs[2]**8/(xs[0]**2 + xs[1]**2 + xs[2]**2)**4 + 2.074716018375e+64*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 1.1969515490625e+64*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 2.17627554375e+63*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 6.04520984375e+61) - 8.0*xs[1]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-5.0)*(-1.6124802676125e+57*xs[2]**7/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(7/2) + 2.6047758169125e+57*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) - 1.1839890076875e+57*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) + 1.315543341875e+56*xs[2]/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2)) - 7.0*xs[1]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.5)*(2.0907567883125e+50*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 2.8510319840625e+50*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 9.503439946875e+49*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 4.52544759375e+48) - 6.0*xs[1]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.0)*(-7.544857768875e+42*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) + 8.38317529875e+42*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) - 1.796394706875e+42*xs[2]/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2)) - 5.0*xs[1]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.5)*(-4.67333269375e+36*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 4.0057137375e+36*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 4.0057137375e+35) + (1.964349645e+30*xs[1]*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) - 3.92869929e+29*xs[1]*xs[2]/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2))*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-2.0) + (1.8693330775e+37*xs[1]*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 8.011427475e+36*xs[1]*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2)*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-2.5) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.5)*(8.89164007875e+64*xs[1]*xs[2]**8/(xs[0]**2 + xs[1]**2 + xs[2]**2)**5 - 1.244829611025e+65*xs[1]*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**4 + 4.78780619625e+64*xs[1]*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 4.3525510875e+63*xs[1]*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.0)*(1.12873618732875e+58*xs[1]*xs[2]**7/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(9/2) - 1.30238790845625e+58*xs[1]*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(7/2) + 3.5519670230625e+57*xs[1]*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) - 1.315543341875e+56*xs[1]*xs[2]/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2)) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.5)*(-1.2544540729875e+51*xs[1]*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**4 + 1.140412793625e+51*xs[1]*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 1.900687989375e+50*xs[1]*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.0)*(3.7724288844375e+43*xs[1]*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(7/2) - 2.514952589625e+43*xs[1]*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) + 1.796394706875e+42*xs[1]*xs[2]/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2))
    Fz = -3.0*xs[2]*(2.63329248e+25*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 8.7776416e+24)*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-2.5) - 4.0*xs[2]*(-6.54783215e+29*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) + 3.92869929e+29*xs[2]/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2))*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.0) - 9.0*xs[2]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-5.5)*(-1.11145500984375e+64*xs[2]**8/(xs[0]**2 + xs[1]**2 + xs[2]**2)**4 + 2.074716018375e+64*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 1.1969515490625e+64*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 2.17627554375e+63*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 6.04520984375e+61) - 8.0*xs[2]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-5.0)*(-1.6124802676125e+57*xs[2]**7/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(7/2) + 2.6047758169125e+57*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) - 1.1839890076875e+57*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) + 1.315543341875e+56*xs[2]/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2)) - 7.0*xs[2]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.5)*(2.0907567883125e+50*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 2.8510319840625e+50*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 9.503439946875e+49*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 4.52544759375e+48) - 6.0*xs[2]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.0)*(-7.544857768875e+42*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) + 8.38317529875e+42*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) - 1.796394706875e+42*xs[2]/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2)) - 5.0*xs[2]*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.5)*(-4.67333269375e+36*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 4.0057137375e+36*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2) - 4.0057137375e+35) + (-5.26658496e+25*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 5.26658496e+25*xs[2]/(xs[0]**2 + xs[1]**2 + xs[2]**2))*(xs[0]**2 + xs[1]**2 + xs[2]**2)**(-1.5) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.5)*(8.89164007875e+64*xs[2]**9/(xs[0]**2 + xs[1]**2 + xs[2]**2)**5 - 2.1339936189e+65*xs[2]**7/(xs[0]**2 + xs[1]**2 + xs[2]**2)**4 + 1.72361023065e+65*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 5.223061305e+64*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 4.3525510875e+63*xs[2]/(xs[0]**2 + xs[1]**2 + xs[2]**2)) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-4.0)*(1.12873618732875e+58*xs[2]**8/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(9/2) - 2.431124095785e+58*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(7/2) + 1.6575846107625e+58*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) - 3.68352135725e+57*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) + 1.315543341875e+56/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2)) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.5)*(-1.2544540729875e+51*xs[2]**7/(xs[0]**2 + xs[1]**2 + xs[2]**2)**4 + 2.3948668666125e+51*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 1.3304815925625e+51*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 1.900687989375e+50*xs[2]/(xs[0]**2 + xs[1]**2 + xs[2]**2)) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-3.0)*(3.7724288844375e+43*xs[2]**6/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(7/2) - 6.2873814740625e+43*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) + 2.6945920603125e+43*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) - 1.796394706875e+42/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2)) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-2.5)*(1.8693330775e+37*xs[2]**5/(xs[0]**2 + xs[1]**2 + xs[2]**2)**3 - 2.670475825e+37*xs[2]**3/(xs[0]**2 + xs[1]**2 + xs[2]**2)**2 + 8.011427475e+36*xs[2]/(xs[0]**2 + xs[1]**2 + xs[2]**2)) + (xs[0]**2 + xs[1]**2 + xs[2]**2)**(-2.0)*(1.964349645e+30*xs[2]**4/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(5/2) - 2.357219574e+30*xs[2]**2/(xs[0]**2 + xs[1]**2 + xs[2]**2)**(3/2) + 3.92869929e+29/sqrt(xs[0]**2 + xs[1]**2 + xs[2]**2))
    
    sum_JNs = (-1)*np.array([Fx, Fy, Fz]) # se necesita (-grad)
    
    zp = np.zeros(6) 
    
    zp[0:3] = xp
    zp[3:6] = Rp_T@(Fg + sum_JNs -(2*(R1@xp) + (R2@x)))    
        
    return zp
    

#from datetime import datetime
#ti = "2018-08-14T22:59:42.000000"
#ti = ti.split("T")
#ti = "{} {}".format(ti[0], ti[1])
#ti = datetime.strptime(ti, '%Y-%m-%d %H:%M:%S.%f')

#tf = "2018-08-16T00:59:42.000000"
#tf = tf.split("T")
#tf = "{} {}".format(tf[0], tf[1])
#tf = datetime.strptime(tf, '%Y-%m-%d %H:%M:%S.%f')

#deltaT = (tf - ti).seconds

#x_i = 1665352.492338
#y_i = -6756092.122778
#z_i = 1281662.400464

#vx_i = -1871.276354
#vy_i = 934.835139
#vz_i = 7304.321490


#x_f = 18066.885210
#y_f = -3451953.355853
#z_f = -6183737.173452

#vx_f = -2407.237376
#vy_f = -6269.867602
#vz_f = 3494.919702

#t = np.linspace(0, deltaT, 9361)

#z0 = np.array([x_i, y_i, z_i, vx_i, vy_i, vz_i])

#sol = odeint(satelite, z0, t)

t = at
x0 = radio + HB
vt = 6820.*3.6 #m/s

z0 = np.array([ax[0], ay[0], az[0], avx[0], avy[0], avz[0]])


sol = odeint(satelite, z0, t )

x = sol[:,0]
y = sol[:,1]
z = sol[:,2]


delta = at[-1]
t = np.linspace(0, delta, len(at))
t1=cron()

sol = odeint(satelite, z0, t)
t2=cron()


deltaT= t2-t1
#print(deltaT)

x = sol[:,0]
y = sol[:,1]
z = sol[:,2]
vx = sol[:,3]
vy = sol[:,4]
vz = sol[:,5]


grad_x = np.gradient(vx,at)
grad_y = np.gradient(vy,at)
grad_z = np.gradient(vz,at)

grad_ax = np.gradient(avx,at)
grad_ay = np.gradient(avy,at)
grad_az = np.gradient(avz,at)

deriva = np.sqrt((x[-1] - ax[-1])**2 + (y[-1] - ay[-1])**2 + (z[-1] - az[-1])**2)

#print (f"{deriva}")

ttf = cron()
ttotal = ttf-tti


with open(Pred,"w") as fout:
    Nt = len(t)
    fout.write("<?xml version=\"1.0\"?>\n")
    fout.write("<Earth_Explorer_File>\n")
    fout.write("  <Earth_Explorer_Header>\n")
    fout.write("    <Fixed_Header>\n")
    fout.write(f"      <File_Name>{archivo}</File_Name>\n")
    fout.write("      <File_Description>Precise Orbit Ephemerides (POE) Orbit File</File_Description>\n")
    fout.write("      <Notes></Notes>\n")
    fout.write("      <Mission>Sentinel-1B</Mission>\n")
    fout.write("      <File_Class>OPER</File_Class>\n")
    fout.write("      <File_Type>AUX_POEORB</File_Type>\n")
    fout.write("      <Validity_Period>\n")
    fout.write("        <Validity_Start>UTC=2020-08-03T22:59:42</Validity_Start>\n")
    fout.write("        <Validity_Stop>UTC=2020-08-05T00:59:42</Validity_Stop>\n")
    fout.write("      </Validity_Period>\n")
    fout.write("      <File_Version>0001</File_Version>\n")
    fout.write("      <Source>\n")
    fout.write("        <System>OPOD</System>\n")
    fout.write("        <Creator>OPOD</Creator>\n")
    fout.write("        <Creator_Version>0.0</Creator_Version>\n")
    fout.write("        <Creation_Date>UTC=2020-08-24T11:11:52</Creation_Date>\n")
    fout.write("      </Source>\n")
    fout.write("    </Fixed_Header>\n")
    fout.write("    <Variable_Header>\n")
    fout.write("      <Ref_Frame>EARTH_FIXED</Ref_Frame>\n")
    fout.write("      <Time_Reference>UTC</Time_Reference>\n")
    fout.write("    </Variable_Header>\n")
    fout.write("  </Earth_Explorer_Header>\n")
    fout.write("<Data_Block type=\"xml\">\n")
    fout.write(f"  <List_of_OSVs count=\"{len(at)}\">\n")
    for i in range(len(at)):
        obj = dt.datetime(2020,8,3,22,59,42,000000)
        fecha = (obj + dt.timedelta(seconds=t[i])).strftime("%Y-%m-%dT%H:%M:%S.%f")
        fout.write("    <OSV>\n")
        fout.write(f"      <UTC>UTC={fecha}</UTC>\n")
        fout.write("      <Absolute_Orbit>+22765</Absolute_Orbit>\n")
        fout.write(f"      <X unit=\"m\">{x[i]}</X>\n")
        fout.write(f"      <Y unit=\"m\">{y[i]}</Y>\n")
        fout.write(f"      <Z unit=\"m\">{z[i]}</Z>\n")
        fout.write(f"      <VX unit=\"m/s\">{vx[i]}</VX>\n")
        fout.write(f"      <VY unit=\"m/s\">{vy[i]}</VY>\n")
        fout.write(f"      <VZ unit=\"m/s\">{vz[i]}</VZ>\n")
        fout.write("      <Quality>NOMINAL</Quality>\n")
        fout.write("    </OSV>\n")
    fout.write("  </List_of_OSVs>\n")
    fout.write("</Data_Block>\n</")
    fout.write("Earth_Explorer_File>")       

