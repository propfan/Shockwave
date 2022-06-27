import numpy as np
from matplotlib import pyplot as plt

#JANNAF function
def Gamma_realN2_1000(T):
    #Def para N2
    Mm=28.013/1000
    #For <-1000K
    a1=2.210371497E+04
    a2=-3.818461820E+02
    a3=6.082738360E+00
    a4=-8.530914410E-03
    a5=1.384646189E-05
    a6=-9.625793620E-09
    a7=2.519705809E-12
    
    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)

    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_realN2p1000(T):
    #Def para N2
    Mm=28.013/1000
    #FOR >1000K<6000
    a1=5.877124060E+05
    a2=-2.239249073E+03
    a3=6.066949220E+00
    a4=-6.139685500E-04
    a5=1.491806679E-07
    a6=-1.923105485E-11
    a7=1.061954386E-15
    
    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_realN2p6000(T):
    #Def para N2
    Mm=28.013/1000
    #FOR >1000K
    a1=8.310139160E+08
    a2=-6.420733540E+05
    a3=2.020264635E+02
    a4=-3.065092046E-02
    a5=2.486903333E-06
    a6=-9.705954110E-11
    a7=1.437538881E-15
    
    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

#Diatomic harmonic oscillator functions 
def Thetav(T,Tv):
    Thetav=Tv/(2*T)  
    return Thetav

def z(Thetav):
    z_thetav=(Thetav/np.sinh(Thetav))**2
    return z_thetav

def gammav(Thetav,delta):
    gamma=(7+2*delta*z(Thetav))/(5+2*delta*z(Thetav))
    return gamma

#Propiedades del gas
#Def para N2
Tv=3390

#Construcción gráficas gamma
T_1000 = np.linspace(200, 999, 500)
Tp1000 = np.linspace(1000, 5999, 1000)
Tp6000 = np.linspace(6000, 20000, 1000)
"""
#Solución para JANNAF
g_janaf= Gamma_realN2(T)


#Solución para delta=0
g1_d0=gamma
#Solución para delta=1
g1_d1=gamma
"""
#Plot resultados
plt.plot(T_1000,Gamma_realN2_1000(T_1000), color='red')
plt.plot(T_1000,gammav(Thetav(T_1000,Tv),0), color='blue')
plt.plot(T_1000,gammav(Thetav(T_1000,Tv),1),color='black')

plt.plot(Tp1000,Gamma_realN2p1000(Tp1000), color='red')
plt.plot(Tp1000,gammav(Thetav(Tp1000,Tv),0), color='blue')
plt.plot(Tp1000,gammav(Thetav(Tp1000,Tv),1),color='black')

plt.plot(Tp6000,Gamma_realN2p6000(Tp6000), color='red')
plt.plot(Tp6000,gammav(Thetav(Tp6000,Tv),0), color='blue')
plt.plot(Tp6000,gammav(Thetav(Tp6000,Tv),1),color='black')


plt.grid(color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('T [K]')
plt.ylabel('\u03B3')
plt.show()
