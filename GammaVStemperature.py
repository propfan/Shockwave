import numpy as np
from matplotlib import pyplot as plt

#JANNAF function
def Gamma_realN2_1000(T):
    #Def para N2
    Mm=28.013/1000
    #For <-1000K
    a1=	2.21037150E+04
    a2=	-3.81846182E+02
    a3=	6.08273836E+00
    a4=	-8.53091441E-03
    a5=	1.38464619E-05
    a6=	-9.62579362E-09
    a7=	2.51970581E-12
    
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
    
    a1=	5.87712406E+05
    a2=	-2.23924907E+03
    a3=	6.06694922E+00
    a4=	-6.13968550E-04
    a5=	1.49180668E-07
    a6=	-1.92310549E-11
    a7=	1.06195439E-15
   
    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_realN2p6000(T):
    #Def para N2
    Mm=28.013/1000
    #FOR >6000K
    a1=	8.31013916E+08
    a2=	-6.42073354E+05
    a3=	2.02026464E+02
    a4=	-3.06509205E-02
    a5=	2.48690333E-06
    a6=	-9.70595411E-11
    a7=	1.43753888E-15
    
    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_realO2_1000(T):
    #Def para O2
    Mm=31.9988/1000
    #For <-1000K
    a1=	-3.42556342E+04
    a2=	4.84700097E+02
    a3=	1.11901096E+00
    a4=	4.29388924E-03
    a5=	-6.83630052E-07
    a6=	-2.02337270E-09
    a7=	1.03904002E-12
   
    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)

    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_realO2p1000(T):
    #Def para O2
    Mm=31.9988/1000
    #FOR >1000K<6000
    a1=	-1.03793902E+06
    a2=	2.34483028E+03
    a3=	1.81973204E+00
    a4=	1.26784758E-03
    a5=	-2.18806799E-07
    a6=	2.05371957E-11
    a7=	-8.19346705E-16
   
    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_realO2p6000(T):
    #Def para O2
    Mm=31.9988/1000
    #FOR >6000K
    a1=	4.97529430E+08
    a2=	-2.86610687E+05
    a3=	6.69035225E+01
    a4=	-6.16995902E-03
    a5=	3.01639603E-07
    a6=	-7.42141660E-12
    a7=	7.27817577E-17

    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_realAir_1000(T):
    #Def para Air
    Mm=28.9651159/1000
    #For <-1000K
    a1=	1.00995016E+04
    a2=	-1.96827561E+02
    a3=	5.00915511E+00
    a4=	-5.76101373E-03
    a5=	1.06685993E-05
    a6=	-7.94029797E-09
    a7=	2.18523191E-12

    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)

    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_realAirp1000(T):
    #Def para Air
    Mm=28.9651159/1000
    #FOR >1000K<6000
    a1=	2.41521443E+05
    a2=	-1.25787460E+03
    a3=	5.14455867E+00
    a4=	-2.13854179E-04
    a5=	7.06522784E-08
    a6=	-1.07148349E-11
    a7=	6.57780015E-16
    
    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_realH2_1000(T):
    #Def para H2
    Mm=2.01588/1000
    #For <-1000K
    a1=	4.07832321E+04
    a2=	-8.00918604E+02
    a3=	8.21470201E+00
    a4=	-1.26971446E-02
    a5=	1.75360508E-05
    a6=	-1.20286027E-08
    a7=	3.36809349E-12

    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)

    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_realH2p1000(T):
    #Def para H2
    Mm=2.01588/1000
    #FOR >1000K<6000
    a1=	5.60812801E+05
    a2=	-8.37150474E+02
    a3=	2.97536453E+00
    a4=	1.25224912E-03
    a5=	-3.74071619E-07
    a6=	5.93662520E-11
    a7=	-3.60699410E-15
    
    #Constant
    R=8.31451
    Rg=R/Mm
    cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_realH2p6000(T):
    #Def para H2
    Mm=2.01588/1000
    #FOR >6000K
    a1=	4.96688412E+08
    a2=	-3.14754715E+05
    a3=	7.98412188E+01
    a4=	-8.41478921E-03
    a5=	4.75324835E-07
    a6=	-1.37187349E-11
    a7=	1.60546176E-16
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
TvN2=3390
TvO2=2270
#Construcción gráficas gamma
T_1000 = np.linspace(200, 999, 500)
Tp1000 = np.linspace(1000, 5999, 1000)
Tp6000 = np.linspace(6000, 20000, 1000)


#Plot resultados

plt.plot(T_1000,Gamma_realAir_1000(T_1000), color='green',linestyle='-', label="Air")
plt.plot(Tp1000,Gamma_realAirp1000(Tp1000), color='green',linestyle='-')

plt.plot(T_1000,Gamma_realH2_1000(T_1000), color='black',linestyle='dashdot', label="H\u2082")
plt.plot(Tp1000,Gamma_realH2p1000(Tp1000), color='black',linestyle='dashdot')
plt.plot(Tp6000,Gamma_realH2p6000(Tp6000), color='black',linestyle='dashdot')

plt.plot(T_1000,Gamma_realO2_1000(T_1000), color='red',linestyle=':', label="O\u2082")
plt.plot(Tp1000,Gamma_realO2p1000(Tp1000), color='red',linestyle=':')
plt.plot(Tp6000,Gamma_realO2p6000(Tp6000), color='red',linestyle=':')

plt.plot(T_1000,Gamma_realN2_1000(T_1000), color='blue',linestyle='--', label="N\u2082")
plt.plot(Tp1000,Gamma_realN2p1000(Tp1000), color='blue',linestyle='--')
plt.plot(Tp6000,Gamma_realN2p6000(Tp6000), color='blue',linestyle='--')

"""
plt.plot(T_1000,gammav(Thetav(T_1000,TvN2),0), color='blue',linestyle='--', label="\u03B4=0")
plt.plot(Tp1000,gammav(Thetav(Tp1000,TvN2),0), color='blue',linestyle='--')
plt.plot(Tp6000,gammav(Thetav(Tp6000,TvN2),0), color='blue',linestyle='--')

plt.plot(T_1000,gammav(Thetav(T_1000,TvN2),1),color='black',linestyle='--', label="\u03B4=1 for N\u2082")
plt.plot(Tp1000,gammav(Thetav(Tp1000,TvN2),1),color='black',linestyle='--')
plt.plot(Tp6000,gammav(Thetav(Tp6000,TvN2),1),color='black',linestyle='--')

plt.plot(T_1000,gammav(Thetav(T_1000,TvO2),1),color='black', label="\u03B4=1 for O\u2082")
plt.plot(Tp1000,gammav(Thetav(Tp1000,TvO2),1),color='black')
plt.plot(Tp6000,gammav(Thetav(Tp6000,TvO2),1),color='black')
"""

plt.grid(color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('T [K]')
plt.ylabel('\u03B3 [-]')
plt.legend(ncol=4)

plt.xticks(ticks=[600,2000,4000,7000,9000,16000,20000], labels=[600,2000,4000,7000,9000,16000,20000])
plt.show()
