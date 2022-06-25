# -*- coding: utf-8 -*-
"""
Created on Sat Jun 25 18:51:36 2022

@author: pedro
"""
import numpy as np
from scipy.optimize import fsolve

def Gamma_realH2(T):
    R=8.31451
    #Def para H2
    Mm=28.013/1000
    Rg=R/Mm
    if (T >= 200) & (T <= 1000):
        cp=Rg*(3.53100E+00-1.23660E-04*T-5.02999E-7*T**2+2.43530E-09*T**3-1.40881E-12*T**4)
    else:
        cp=Rg*(2.95257e0+1.39690E-3*T-4.92631E-07*T**2+7.89010E-11*T**3-4.60755E-15*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def cp_TtT(T,Tt):
    R=8.31451
    #Def para H2
    Mm=28.013/1000
    Rg=R/Mm
    if (T >= 200) & (T <= 1000):
        AA=+Rg*(3.53100E+00*(Tt-T)-1.23660E-04*0.5*(Tt**2-T**2)-5.02999E-7*(1/3)*(Tt**3-T**3)+2.43530E-09*(1/4)*(Tt**4-T**4)-1.40881E-12*(1/5)*(Tt**5-T**5))  
    else:
        AA=+Rg*(2.95257e0*(Tt-T)+1.39690E-3*0.5*(Tt**2-T**2)-4.92631E-07*(1/3)*(Tt**3-T**3)+7.89010E-11*(1/4)*(Tt**4-T**4)-4.60755E-15*(1/5)*(Tt**5-T**5))  
    return AA

#Variables iniciales
#p0=rho0*R*T0
M1=10
p1=1e5
T1=1300
beta1=np.radians(10)

R=8.31451
#Propiedades del gas
#Def para H2
Mm=28.013/1000
Rg=R/Mm

#Velocidad Tangencial
gamma1=Gamma_realH2(T1)
a1=(gamma1*Rg*T1)**0.5
u1=M1*a1
rho1=p1/(Rg*T1)
w1=u1*np.cos(beta1)

#Velocidad Normal
u1n=u1*np.sin(beta1)
M1n=u1n/a1

#Temperatura remanso
func_ht = lambda Tt : -0.5*u1**2+ cp_TtT(T1, Tt)
Tt_initial_guess = (1+M1n**2*(gamma1-1)*0.5)*T1
Tt_s=fsolve(func_ht, Tt_initial_guess)

#Resolución T2
func_T2= lambda T2 : -Rg*T1/u1-u1+Rg*T2/(2*cp_TtT(T2, Tt_s))**0.5+(2*cp_TtT(T2, Tt_s))**0.5
T2_initial_guess = Tt_s*0.99
T2_s=fsolve(func_T2, T2_initial_guess)

u2n=(2*cp_TtT(T2_s, Tt_s))**0.5
u2=(u2n**2+w1**2)**0.5

#Saltos
T21=T2_s/T1
u21n=u2n/u1n
rho21=1/u21n
p2=(Gamma_realH2(T2_s)*Rg*T2_s*rho21*rho1)/Gamma_realH2(T2_s)
p21=p2/p1

gamma2=Gamma_realH2(T2_s)
a2=(gamma2*Rg*T2_s)**0.5
M2n=u2n/a2
M2=u2/a2

print('Solution for Tt',Tt_s)
print('Solution for T2',T2_s)

print('M1',M1)
print('M1n',M1n)
print('M2',M2)
print('M2n',M2n)

print('Jump results')
print('Jump T21',T21)
print('Jump rho_21',rho21)
print('Jump p21',p21)
