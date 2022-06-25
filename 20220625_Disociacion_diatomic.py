# -*- coding: utf-8 -*-
"""
Created on Sun Jun 12 01:06:27 2022

@author: pedro
"""
import numpy as np
import sympy 
from sympy import Symbol, nsolve, exp, sqrt

P = Symbol('P') 
T = Symbol('T') 
R = Symbol('R')
A = Symbol('A')

#Inputs
T1=300 
p1=1E5

#Constant
R_cte=8.31451


#Propiedades del gas
#Def para O2
Mm=31.9988/1000
Rg=R_cte/Mm

rho1=p1/(Rg*T1)
a1=(Rg*7/5*T1)**0.5

#Def para 02
#-----------***-----------
bv = 2270/T1 
bd = 59500/T1 
M1 = 6600/a1

G=25/3 
m= 2.6567*10**-26 
tr=2.08 
kb= 1.380648*10**-23 
h = 1.05457*10**-34
B=(G*m*tr*T1**0.5*(np.pi*m*kb/h**2)**(3/2))/rho1
#-----------***-----------

print('B=',B) 
print('bv=',bv) 
print('bd=',bd)

f1 = -P+(1 + A)*R*T 
f2 = -P + 1+(7/5)*(1-1/R)*M1**2 
f3 = -T+(6-1/R-2*A*bd-2*(1-A)*bv/(sympy.exp(bv/T)-1))/(2*(A+3)-(1+ A)*R) 
f4 = -A**2/(1-A)+B*sympy.exp(-bd/T)*(T**0.5)*(1-sympy.exp(-bv/T))/R
#Aproximación
#f4 = -A**2/(1-A)+R*B*sympy.exp(-bd/T)

Pjump,Rjump, Tjump, alpha=nsolve((f1, f2, f3, f4), (P, R, T, A), (100, 11, 25,0.8))

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




