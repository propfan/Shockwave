# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 01:06:27 2022

@author: pedro
"""
import numpy as np
from scipy.optimize import fsolve

def f(variables):
  (P, R, T, A)= variables
  f1 = -P+(1 + A)*R*T 
  f2 = -P + 1+(7/5)*(1-1/R)*M1**2 
  f3 = -T+(6-1/R-2*A*bd-2*(1-A)*bv/(np.exp(bv/T)-1))/(2*(A+3)-(1+ A)*R) 
  f4 = -A**2/(1-A)+B*np.exp(-bd/T)*(T**0.5)*(1-np.exp(-bv/T))/R
  #Aproximación
  #f4 = -A**2/(1-A)+R*B*np.exp(-bd/T)
  return [f1,f2,f3,f4]

#Inputs
T1=300 
p1=1E5
u1=6600

#-----------***-----------
#Propiedades del gas
#Def para O2
Mm=31.9988/1000
Tv=2270
Td=59500
#-----------***-----------

#Constant
R_cte=8.31451
G=25/3 
m= 2.6567*10**-26 
tr=2.08 
kb= 1.380648*10**-23 
h = 1.05457*10**-34

Rg=R_cte/Mm
rho1=p1/(Rg*T1)

a1=(Rg*7/5*T1)**0.5
bv = Tv/T1 
bd = Td/T1 
M1 = u1/a1
B=(G*m*tr*T1**0.5*(np.pi*m*kb/h**2)**(3/2))/rho1

print('B=',B) 
print('bv=',bv) 
print('bd=',bd)

(Pjump,Rjump, Tjump, alpha)=fsolve(f, (100, 11, 25,0.8))

alpha_r=alpha*(1-alpha)/(2-alpha)
alpha_T=-alpha_r*(0.5+bd*(1-(1+bv/bd)*sympy.exp(bv/Tjump))/(Tjump*(1-sympy.exp(-bv/Tjump))))

evib=(bv/Tjump)/(-1+sympy.exp(bv/Tjump))
a21=((5*Tjump/7)*(1+alpha+alpha_r+(1+alpha+alpha_T)*(2*(1+alpha)-alpha_r*(1-evib+2*bd/Tjump))/(5+alpha+2*(1-alpha)*evib**2*sympy.exp(-bv/Tjump)+alpha_T*(1-2*evib+2*bd/Tjump))))**0.5

M2=M1/(Rjump*a21)

print('M1=',M1) 
print('M2=',M2)
print('P21=',Pjump)
print('T21=',Tjump)
print('rho21=',Rjump)
print('Alpha Dissotiation=',alpha)
