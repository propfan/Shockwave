# -*- coding: utf-8 -*-
"""
Created on Sun Sep 11 16:55:20 2022

@author: pedro
"""

import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

#Variables iniciales
#p0=rho0*R*T0
T1=300

Mm=28.013/1000
Tv=3390

def R_g(Mm):
    return 8.31451/Mm

#----------**Oscillator APROX**----------
def Thetav(T,Tv):
    Thetav=Tv/(2*T)  
    return Thetav

def T_INV(Thetav,Tv):
    T=Tv/(2*Thetav)  
    return T

def Thetavt_delta0(T, Tv,M):
    Thetavt=Thetav(T,Tv)/(1+(1/5)*M**2) 
    return Thetavt

def z(Thetav):
    z_thetav=(Thetav/np.sinh(Thetav))**2
    return z_thetav

def gamma(Thetav,delta):
    gamma=(7+2*delta*z(Thetav))/(5+2*delta*z(Thetav))
    return gamma

def rho_jump (gM1,gM2):
    Rho_jump=(gM1/gM2)*((1+gM2)/(1+gM1))
    return Rho_jump

def p_jump (gM1,gM2):
    p_jump=((1+gM1)/(1+gM2))
    return p_jump

#----------**Oscillator APROX**----------


def Oscillator_NormalShock (M1,T1):
    Rg=R_g(Mm)
    
    #----------**Oscillator APROX Delta=0**----------
    #Condiciones zona1
    Thetav1=Thetav(T1,Tv)
    g1_d0=gamma(Thetav1,0)
    g1M1_d0=g1_d0*M1**2
    
    #Solucion Thetavt
    Thetavt_d0=Thetavt_delta0(T1, Tv,M1)

    #Solucion Thetav2
    func_Thetav2_d0 = lambda Thetav2d0: -Thetav2d0/Thetav1+(g1M1_d0/(7*(-1+Thetav2d0/Thetavt_d0)))*((1+7*(-1+Thetav2d0/Thetavt_d0))/(1+g1M1_d0))**2
    Thetav2_d0_initial_guess = Thetav(T1*((-1+7*M1**2)*(5+M1**2))/(36*M1**2),Tv)
    Thetav2_d0=fsolve(func_Thetav2_d0, Thetav2_d0_initial_guess)

    g2_d0=gamma(Thetav2_d0,0)
    M2_d0=(7*(-1+Thetav2_d0/Thetavt_d0)/g2_d0)**0.5
    g2M2_d0=g2_d0*M2_d0**2

    #Solución Saltos
    T2_d0=T_INV(Thetav2_d0,Tv)
    R21_d0=rho_jump(g1M1_d0,g2M2_d0)
    p21_d0=p_jump(g1M1_d0,g2M2_d0)
    T21_d0=T2_d0/T1
    #----------**Oscillator APROX Delta=0**----------
    #----------**Oscillator APROX Delta=1**----------
    #Condiciones zona 1
    g1=gamma(Thetav1,1)
    g1M1=g1*M1**2

    #Solucion Thetavt
    func_Thetavt = lambda Thetavt: -g1M1+7*(-1+Thetav1/Thetavt)+2*Thetav1*(1/np.tanh(Thetavt)-1/np.tanh(Thetav1))
    Thetavt_initial_guess = Thetavt_d0
    Thetavt_s=fsolve(func_Thetavt, Thetavt_initial_guess)
    
    #Solucion Thetav2
    Thetav2_initial_guess = Thetav2_d0
    M2_initial_guess = M2_d0
    
    def f(variables):
        #Resolución Thetav2,M2 para delta=1
        (Thetav2,M2)= variables
        f1=-Thetav2/Thetav1+(g1M1/(gamma(Thetav2,1)*M2**2))*((1+gamma(Thetav2,1)*M2**2)/(1+g1M1))**2
        f2=-gamma(Thetav2,1)*M2**2+7*(-1+Thetav2/Thetavt_s)+2*Thetav2*(1/np.tanh(Thetavt_s)-1/np.tanh(Thetav2))
        return [f1,f2]
    
    (Thetav2_s,M2_s)=fsolve(f, (Thetav2_initial_guess,M2_initial_guess))

    #Calculo del resto de variables
    g2=gamma(Thetav2_s,1)
    M2=((7*(-1+Thetav2_s/Thetavt_s)+2*Thetav2_s*(1/np.tanh(Thetavt_s)-1/np.tanh(Thetav2_s)))/g2)**0.5
    g2M2=g2*M2**2

    #Solución Saltos
    T2=T_INV(Thetav2_s,Tv)
    R21=rho_jump(g1M1,g2M2)
    p21=p_jump(g1M1,g2M2)
    T21=T2/T1
    #----------**Oscillator APROX Delta=1**----------
    return T21_d0,R21_d0,p21_d0,T21,R21,p21

M1_in = []
T21_s_d0= []
R21_s_d0= []
p21_s_d0= []
T21_s= []
R21_s= []
p21_s= []
for M1 in np.arange(2,10, 0.1):   
    T21_d0,R21_d0,p21_d0,T21,R21,p21=Oscillator_NormalShock (M1,T1)
    
    M1_in = np.append(M1_in,M1)
    T21_s_d0= np.append(T21_s_d0,T21_d0)
    R21_s_d0= np.append(R21_s_d0,R21_d0)
    p21_s_d0= np.append(p21_s_d0,p21_d0)
    T21_s= np.append(T21_s,T21)
    R21_s= np.append(R21_s,R21)
    p21_s= np.append(p21_s,p21)

#Plot
plt.plot(M1_in,T21_s, color='green', label="Salto T")
plt.plot(M1_in,T21_s_d0, color='green', linestyle='--', label="Salto T para d0")
plt.plot(M1_in,R21_s,color='blue', linewidth=0.5, label="Salto \u03C1")
plt.plot(M1_in,R21_s_d0,color='blue', linewidth=0.5, linestyle='--', label="Salto \u03C1 para d0")
plt.plot(M1_in,p21_s, color='black', label="Salto p")
plt.plot(M1_in,p21_s_d0, color='black', linestyle='--', label="Salto p para d0")

plt.grid(color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('M1')
plt.ylabel('Salto de variables')
plt.legend(ncol=2)
plt.show()
