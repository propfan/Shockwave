# -*- coding: utf-8 -*-
"""
Created on Sun Jun 26 19:22:35 2022

@author: pedro
"""
import numpy as np
from scipy.optimize import fsolve

def f(variables):
    (P, R, T, A)= variables
    f1 = -P+(1 + A)*R*T 
    f2 = -P + 1+g1*(1-1/R)*M1**2 
    f3 = -T+(4-1/R-2*A*Theta_ion)/((4-R)*(1+A))
    f4 = -A**2/(1-A)+Bi*T**(3/2)*np.exp(-Theta_ion/T)/R
    return [f1,f2,f3,f4]


#Inputs
T1=300 
p1=666
M1=30

#Constant
R_cte=8.31451           #Ideal gas constant[J/molK]
kb= 1.380648E-23        #Boltzmann constant[J/K]
kv_ev=8.617333262E-5    #Boltzmann constant[eV/K]
h = 1.05457E-34         #Planck constant[J/s]
me=9.1093837E-31        #Electron mass [kg]
#Aproximacion - Stability of ionizing shock waves in monatomic gases-
Gi=1E-10                


#Propiedades del gas monoatómico
#Def para O
Mm=15.9994/1000
Ii=13.6181 #[eV]

"""
#Ar
Mm=0.03995
Ii=15.7596
"""
Rg=R_cte/Mm

g1=5/3
a1=(Rg*g1*T1)**0.5

Ti=Ii/kv_ev
#Non-dimension T
Theta_ion=Ti/T1

Bi=Gi*(me/(2*np.pi*h**2))**(3/2)*(kb*T1)**(5/2)/p1

print(Theta_ion)
print(Bi)

(Pjump,Rjump, Tjump, alpha)=fsolve(f, (3000, 12, 50,0.8))

#Gas monoatomico ionizado
g2=5/3
M2=M1/(Rjump*g2)

print('M1=',M1) 
print('M2=',M2)
print('P21=',Pjump)
print('T21=',Tjump)
print('rho21=',Rjump)
print('Alpha Dissotiation=',alpha)