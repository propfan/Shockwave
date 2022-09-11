# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 20:19:13 2022

@author: pedro
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

#Variables iniciales
#p0=rho0*R*T0
M1=5
p1=101325
T1=288
beta_punto2=np.deg2rad(21.5)

#Polynomial from NASA 200-20000 K
#Values for N2
a=np.array([[2.21037150E+04,-3.81846182E+02,6.08273836E+00,-8.53091441E-03,1.38464619E-05,-9.62579362E-09,2.51970581E-12],
            [5.87712406E+05,-2.23924907E+03,6.06694922E+00,-6.13968550E-04,1.49180668E-07,-1.92310549E-11,1.06195439E-15],
            [8.31013916E+08,-6.42073354E+05,2.02026464E+02,-3.06509205E-02,2.48690333E-06,-9.70595411E-11,1.43753888E-15]])
Mm=28.013/1000

def R_g(Mm):
    return 8.31451/Mm

def Gamma(T,Mm):
    Rg=R_g(Mm)
    global cp
    if (T >= 200) & (T<=1000):
        cp=Rg*(a[0,0]*T**-2+a[0,1]*T**-1+a[0,2]+a[0,3]*T+a[0,4]*T**2+a[0,5]*T**3+a[0,6]*T**4)
    elif (T > 1000) & (T<=6000):
        cp=Rg*(a[1,0]*T**-2+a[1,1]*T**-1+a[1,2]+a[1,3]*T+a[1,4]*T**2+a[1,5]*T**3+a[1,6]*T**4)
    elif (T > 6000) & (T<=20000):
        cp=Rg*(a[2,0]*T**-2+a[2,1]*T**-1+a[2,2]+a[2,3]*T+a[2,4]*T**2+a[2,5]*T**3+a[2,6]*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def cp_TtT(T,Tt,Mm):
    Rg=R_g(Mm)
    global AA
    if (Tt >= 200) & (Tt<=1000):
        AA=Rg*(-a[0,0]*(1/Tt-1/T)+a[0,1]*np.log(Tt/T)+a[0,2]*(Tt-T)+a[0,3]*(1/2)*(Tt**2-T**2)+a[0,4]*(1/3)*(Tt**3-T**3)+a[0,5]*(1/4)*(Tt**4-T**4)+a[0,6]*(1/5)*(Tt**5-T**5))    
    elif (Tt > 1000) & (Tt<=6000):
        AA=Rg*(-a[1,0]*(1/Tt-1/T)+a[1,1]*np.log(Tt/T)+a[1,2]*(Tt-T)+a[1,3]*(1/2)*(Tt**2-T**2)+a[1,4]*(1/3)*(Tt**3-T**3)+a[1,5]*(1/4)*(Tt**4-T**4)+a[1,6]*(1/5)*(Tt**5-T**5))    
    elif (Tt > 6000) & (Tt<=20000):
        AA=Rg*(-a[2,0]*(1/Tt-1/T)+a[2,1]*np.log(Tt/T)+a[2,2]*(Tt-T)+a[2,3]*(1/2)*(Tt**2-T**2)+a[2,4]*(1/3)*(Tt**3-T**3)+a[2,5]*(1/4)*(Tt**4-T**4)+a[2,6]*(1/5)*(Tt**5-T**5))    
    return AA

def funShockI(M1,p1,T1, beta):
    #Velocidad Tangencial
    Rg=R_g(Mm)
    g1=Gamma(T1,Mm)
    a1=(g1*Rg*T1)**0.5
    u1=M1*a1
    rho1=p1/(Rg*T1)
    w1=u1*np.cos(beta)
    
    #Velocidad Normal
    u1n=u1*np.sin(beta)
    M1n=u1n/a1
    
    #Temperatura remanso
    func_ht = lambda Tt : -u1**2+2*cp_TtT(T1, Tt,Mm)
    Tt_initial_guess = (1+M1n**2*(g1-1)*0.5)*T1
    Tt_s=fsolve(func_ht, Tt_initial_guess)
    
    #Resolución T2
    func_T2= lambda T2 : -Rg*T1/u1n-u1n+Rg*T2/(2*cp_TtT(T2,Tt_s,Mm)-w1**2)**0.5+(2*cp_TtT(T2,Tt_s,Mm)-w1**2)**0.5
    T2_initial_guess = (1+M1n**2*(g1-1)*0.5)*T1*0.85
    T2_s=fsolve(func_T2, T2_initial_guess)
    
    u2n=(2*cp_TtT(T2_s, Tt_s,Mm)-w1**2)**0.5
    u2=(u2n**2+w1**2)**0.5
    
    Theta_2rad=beta-np.arcsin(u2n/u2)
    Theta_2deg=np.rad2deg(Theta_2rad)
    
    g2=Gamma(T2_s,Mm)
    a2=(g2*Rg*T2_s)**0.5
    M2=u2/a2
    M2n=u2n/a2
    
    #Saltos
    T21=T2_s/T1
    u21n=u2n/u1n
    rho21=1/u21n
    p2=T2_s*Rg*rho1*rho21
    p21=p2/p1
    rho2=rho21*rho1
    return p21, Theta_2deg, M2, T2_s, p2,rho21,u2,rho2


BETA_in = []
T2_s= []
R21_s= []
p21_s= []
Theta_deg_s=[]
for beta1 in np.arange(np.deg2rad(5),np.deg2rad(19.4), 0.001):   
  
    p21,Theta_deg,M2,T2,p2,R21,u2_s,ro2=funShockI (M1,p1,T1, beta1)
    
    BETA_in = np.append(BETA_in,beta1)
    T2_s= np.append(T2_s,T2)
    R21_s= np.append(R21_s,R21)
    p21_s= np.append(p21_s,p21)
    Theta_deg_s= np.append(Theta_deg_s,Theta_deg)
    print(Theta_deg)

#Plot
plt.plot(BETA_in,T2_s, color='green', label="Salto T")
plt.plot(BETA_in,R21_s,color='blue', label="Salto \u03C1")
plt.plot(BETA_in,p21_s, color='black', label="Salto p")
plt.plot(BETA_in,Theta_deg_s, color='red', label="Theta [º]")

plt.grid(color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('M1')
plt.ylabel('Salto de variables')
plt.legend(ncol=3)
plt.show()
