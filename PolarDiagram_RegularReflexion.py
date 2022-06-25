# -*- coding: utf-8 -*-
"""
Created on Mon May 23 20:34:18 2022

@author: pedro
"""
import numpy as np
from matplotlib import pyplot as plt

#Condiciones de entrada 
Mi=2.8
gamma=1.4
beta_i=35*np.pi/180

#Construccion de la polar incidente
def funShockI(g, M, b):
    #Se definen parametros onda incidente
    Mn=M*np.sin(b)
    P=1+(2*g)*(Mn**2-1)/(g+1)

    Theta_rad=np.arctan((2*(Mn**2-1)/np.tan(b))/(M**2*(g+np.cos(2*b))+2))
    Theta_i_d=Theta_rad*180/np.pi 

    return P, Theta_i_d

#Construccion de la polar reflejada
def funShockR(g, M1, b1):
    M1n=M1*np.sin(b1)    
    M2=((((g-1)*M1n**2+2)**2+(g+1)**2*M1n**2*M1**2*(np.cos(b1))**2)/((2*g*M1n**2-(g-1))*((g-1)*M1n**2+2)))**0.5
    Theta1=np.arctan((2*(M1n**2-1)/np.tan(b1))/(M1**2*(g+np.cos(2*b1))+2))
    
    M2n=((M1n**2+2/(g-1))/(M1n**2*(2*g/(g-1))-1))**0.5
    
    #Se introducen valores de la zona reflejada
    Pr_0=1/(1+(2*g)*(M2n**2-1)/(g+1))
    return M2, Pr_0, Theta1


def funSchock2r(g,M2,b2):
    M2n=M2*np.sin(b2)
    Pr=(1+(2*g)*(M2n**2-1)/(g+1))
    #hasta aqui benne    
    Theta2=np.arctan((2*(M2n**2-1)/np.tan(b2))/(M2**2*(g+np.cos(2*b2))+2))
    Theta2_d=Theta2*180/np.pi
    
    return Pr, Theta2_d


#Contrucción de la rama positiva y negatuva onda incidente
b_p =np.linspace(np.arcsin(1/Mi),np.pi/2, 100)
b_n =np.linspace(-np.pi/2,-np.arcsin(1/Mi), 100)

#Llamada función de onda incidente
P_p, Theta_d_p= funShockI(gamma,Mi,b_p)
P_n, Theta_d_n= funShockI(gamma,Mi,b_n)

#Llamada función de onda reflejada
M2, Pr_0, Theta1=funShockR(gamma, Mi, beta_i)

#Construcción de la rama negativa onda reflejada
b_r =np.linspace(-np.pi/2, -np.arcsin(1/M2), 100)
Pr, Theta2_d=funSchock2r(gamma, M2, b_r)

print(M2, Pr_0, Theta1*180/np.pi)
#Plot onda reflejada
plt.plot(Theta_d_p, P_p, color='red')
plt.plot(Theta_d_n, P_n, color='blue')
plt.plot(Theta1*180/np.pi+Theta2_d, Pr*Pr_0, color='black')

plt.grid(color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('\u03B8')
plt.ylabel('p/p\u2080')
plt.show()
