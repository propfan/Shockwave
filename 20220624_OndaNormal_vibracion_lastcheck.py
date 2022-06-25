# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 16:46:46 2022

@author: pedro
"""

import numpy as np
from scipy.optimize import fsolve


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

def Rho_21 (gM1,gM2):
    Rho_jump=(gM1/gM2)*((1+gM2)/(1+gM1))
    return Rho_jump

def p_21 (gM1,gM2):
    p_jump=((1+gM1)/(1+gM2))
    return p_jump

#Condiciones iniciales
M1=6
T1=300
R=8.31451

#Propiedades del gas
#Def para N2
Mm=28.013/1000
Tv=3390

Rg=R/Mm

#Solución para delta=0
#Condiciones zona1
Thetav1=Thetav(T1,Tv)
g1_d0=gamma(Thetav1,0)
g1M1_d0=g1_d0*M1**2

Thetavt_d0=Thetavt_delta0(T1, Tv,M1)

func_Thetav2_d0 = lambda Thetav2d0: -Thetav2d0/Thetav1+(g1M1_d0/(7*(-1+Thetav2d0/Thetavt_d0)))*((1+7*(-1+Thetav2d0/Thetavt_d0))/(1+g1M1_d0))**2
Thetav2_d0_initial_guess = Thetav(T1*((-1+7*M1**2)*(5+M1**2))/(36*M1**2),Tv)
Thetav2_d0=fsolve(func_Thetav2_d0, Thetav2_d0_initial_guess)

g2_d0=gamma(Thetav2_d0,0)


M2_d0=(7*(-1+Thetav2_d0/Thetavt_d0)/g2_d0)**0.5
g2M2_d0=g2_d0*M2_d0**2

T2_d0=T_INV(Thetav2_d0,Tv)
rho_21_d0=Rho_21(g1M1_d0,g2M2_d0)
p_21_d0=p_21(g1M1_d0,g2M2_d0)
T_21_d0=T2_d0/T1

print('Aproximación despreciando la vibración')
print('M2 para delta=0',M2_d0)
print('T2 para delta=0',T2_d0)
print('Salto de presiones p21 para delta=0',p_21_d0)
print('Salto de densidades rho21 para delta=0',rho_21_d0)
print('Salto de temperaturas T21 para delta=0',T_21_d0)

#HASTA AQUI CODIGO COMPROBADO PARA DELTA=0-----------------------------------------------

#Cálculo para Delta=1
#Condiciones zona 1
g1=gamma(Thetav1,1)
g1M1=g1*M1**2

func_Thetavt = lambda Thetavt: -g1M1+7*(-1+Thetav1/Thetavt)+2*Thetav1*(1/np.tanh(Thetavt)-1/np.tanh(Thetav1))
Thetavt_initial_guess = Thetavt_d0
Thetavt_s=fsolve(func_Thetavt, Thetavt_initial_guess)

Thetav2_initial_guess = Thetav2_d0
M2_initial_guess = M2_d0

def f(variables):
    (Thetav2,M2)= variables
    f1=-Thetav2/Thetav1+(g1M1/(gamma(Thetav2,1)*M2**2))*((1+gamma(Thetav2,1)*M2**2)/(1+g1M1))**2
    f2=-gamma(Thetav2,1)*M2**2+7*(-1+Thetav2/Thetavt_s)+2*Thetav2*(1/np.tanh(Thetavt_s)-1/np.tanh(Thetav2))
    return [f1,f2]

(Thetav2_s,M2_s)=fsolve(f, (Thetav2_initial_guess,M2_initial_guess))

"""
#Calculo a cholon
func_Thetav2 = lambda Thetav2:-Thetav2/Thetav1+(g1M1/(7*(-1+Thetav2/Thetavt_s)+2*Thetav2*(1/np.tanh(Thetavt_s)-1/np.tanh(Thetav2))))*((1+7*(-1+Thetav2/Thetavt_s)+2*Thetav2*(1/np.tanh(Thetavt_s)-1/np.tanh(Thetav2)))/(1+g1M1))**2
Thetav2_initial_guess = -Thetav2_d0*1.0
Thetav2_s=fsolve(func_Thetav2, Thetav2_initial_guess)
"""


#Calculo del resto de variables
g2=gamma(Thetav2_s,1)

#Por aqui es posible que se encuentre el error
M2=((7*(-1+Thetav2_s/Thetavt_s)+2*Thetav2_s*(1/np.tanh(Thetavt_s)-1/np.tanh(Thetav2_s)))/g2)**0.5
g2M2=g2*M2**2

T2=T_INV(Thetav2_s,Tv)
rho_21=Rho_21(g1M1,g2M2)
p_21=p_21(g1M1,g2M2)
T_21=T2/T1
print('Aproximación incluyendo la vibración')
print('M2 para delta=1',M2)
print('T2 para delta=1',T2)
print('Salto de presiones p21 para delta=1',p_21)
print('Salto de densidades rho21 para delta=1',rho_21)
print('Salto de temperaturas T21 para delta=1',T_21)
