import numpy as np
from matplotlib import pyplot as plt

#JANNAF function
def Gamma_N2(T):
     #Def para N2
    Mm=28.013/1E3
    a=np.array([[2.21037150E+04,-3.81846182E+02,6.08273836E+00,-8.53091441E-03,1.38464619E-05,-9.62579362E-09,2.51970581E-12],
                [5.87712406E+05,-2.23924907E+03,6.06694922E+00,-6.13968550E-04,1.49180668E-07,-1.92310549E-11,1.06195439E-15],
                [8.31013916E+08,-6.42073354E+05,2.02026464E+02,-3.06509205E-02,2.48690333E-06,-9.70595411E-11,1.43753888E-15]])
    
    R=8.31451
    Rg=R/Mm
    if ((T >= 200) & (T<=1000)):
        cp=Rg*(a[0,0]*T**-2+a[0,1]*T**-1+a[0,2]+a[0,3]*T+a[0,4]*T**2+a[0,5]*T**3+a[0,6]*T**4)
    elif ((T > 1000) & (T<=6000)):
        cp=Rg*(a[1,0]*T**-2+a[1,1]*T**-1+a[1,2]+a[1,3]*T+a[1,4]*T**2+a[1,5]*T**3+a[1,6]*T**4)
    elif ((T > 6000) & (T<=20000)):
        cp=Rg*(a[2,0]*T**-2+a[2,1]*T**-1+a[2,2]+a[2,3]*T+a[2,4]*T**2+a[2,5]*T**3+a[2,6]*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_O2(T):
     #Def para O2
    Mm=31.9988/1E3
    a=np.array ([[-3.42556342E+04,4.84700097E+02,1.11901096E+00,4.29388924E-03,-6.83630052E-07,-2.02337270E-09,1.03904002E-12],
                 [-1.03793902E+06,2.34483028E+03,1.81973204E+00,1.26784758E-03,-2.18806799E-07,2.05371957E-11,-8.19346705E-16],
                 [4.97529430E+08,-2.86610687E+05,6.69035225E+01,-6.16995902E-03,3.01639603E-07,-7.42141660E-12,7.27817577E-17]])
    R=8.31451
    Rg=R/Mm
    if ((T >= 200) & (T<=1000)):
        cp=Rg*(a[0,0]*T**-2+a[0,1]*T**-1+a[0,2]+a[0,3]*T+a[0,4]*T**2+a[0,5]*T**3+a[0,6]*T**4)
    elif ((T > 1000) & (T<=6000)):
        cp=Rg*(a[1,0]*T**-2+a[1,1]*T**-1+a[1,2]+a[1,3]*T+a[1,4]*T**2+a[1,5]*T**3+a[1,6]*T**4)
    elif ((T > 6000) & (T<=20000)):
        cp=Rg*(a[2,0]*T**-2+a[2,1]*T**-1+a[2,2]+a[2,3]*T+a[2,4]*T**2+a[2,5]*T**3+a[2,6]*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_H2(T):
     #Def para H2
    Mm=2.01588/1E3
    a=np.array([[4.07832321E+04,-8.00918604E+02,8.21470201E+00,-1.26971446E-02,1.75360508E-05,-1.20286027E-08,3.36809349E-12],
                [5.60812801E+05,-8.37150474E+02,2.97536453E+00,1.25224912E-03,-3.74071619E-07,5.93662520E-11,-3.60699410E-15],
                [4.96688412E+08,-3.14754715E+05,7.98412188E+01,-8.41478921E-03,4.75324835E-07,-1.37187349E-11,1.60546176E-16]])    
    
    R=8.31451
    Rg=R/Mm
    if ((T >= 200) & (T<=1000)):
        cp=Rg*(a[0,0]*T**-2+a[0,1]*T**-1+a[0,2]+a[0,3]*T+a[0,4]*T**2+a[0,5]*T**3+a[0,6]*T**4)
    elif ((T > 1000) & (T<=6000)):
        cp=Rg*(a[1,0]*T**-2+a[1,1]*T**-1+a[1,2]+a[1,3]*T+a[1,4]*T**2+a[1,5]*T**3+a[1,6]*T**4)
    elif ((T > 6000) & (T<=20000)):
        cp=Rg*(a[2,0]*T**-2+a[2,1]*T**-1+a[2,2]+a[2,3]*T+a[2,4]*T**2+a[2,5]*T**3+a[2,6]*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def Gamma_Air(T):
     #Def para Air
    Mm=28.9651159/1e3
    a=np.array([[1.00995016E+04,-1.96827561E+02,5.00915511E+00,-5.76101373E-03,1.06685993E-05,-7.94029797E-09,2.18523191E-12],
                [2.41521443E+05,-1.25787460E+03,5.14455867E+00,-2.13854179E-04,7.06522784E-08,-1.07148349E-11,6.57780015E-16]])
    
    R=8.31451
    Rg=R/Mm
    if ((T >= 200) & (T<=1000)):
        cp=Rg*(a[0,0]*T**-2+a[0,1]*T**-1+a[0,2]+a[0,3]*T+a[0,4]*T**2+a[0,5]*T**3+a[0,6]*T**4)
    else:
        cp=Rg*(a[1,0]*T**-2+a[1,1]*T**-1+a[1,2]+a[1,3]*T+a[1,4]*T**2+a[1,5]*T**3+a[1,6]*T**4)
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

Tin=[]
Tin2=[]
GammaN2_J=[]
GammaAir_J=[]
GammaH2_J=[]
GammaO2_J=[]
GammaN2_d0=[]
GammaN2_d1=[]
GammaO2_d1=[]

for T in range (200,20000, 100):
   Tin= np.append(Tin,T)
   GammaN2_J= np.append(GammaN2_J,Gamma_N2(T))
   GammaO2_J= np.append(GammaO2_J,Gamma_O2(T))
   GammaH2_J= np.append(GammaH2_J,Gamma_H2(T))
   GammaN2_d0= np.append(GammaN2_d0,gammav(Thetav(T,TvN2),0))
   GammaN2_d1= np.append(GammaN2_d1,gammav(Thetav(T,TvN2),1))
   GammaO2_d1= np.append(GammaO2_d1,gammav(Thetav(T,TvO2),1))

for T2 in range (200,6000, 100):
   Tin2= np.append(Tin2,T2)
   GammaAir_J= np.append(GammaAir_J,Gamma_Air(T2))


plt.plot(Tin,GammaN2_d0,color='blue',linestyle='--', label="\u03B4=0")
plt.plot(Tin,GammaN2_d1,color='black',linestyle='--', label="\u03B4=1 for N\u2082")
plt.plot(Tin,GammaO2_d1,color='black', label="\u03B4=1 for O\u2082")

plt.plot(Tin2,GammaAir_J, color='green',linestyle='-', label="Air")
plt.plot(Tin,GammaN2_J, color='blue',linestyle='--', label="N\u2082")
plt.plot(Tin,GammaO2_J, color='red',linestyle=':', label="O\u2082")
plt.plot(Tin,GammaH2_J, color='black',linestyle='dashdot', label="H\u2082")
plt.grid(color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('T [K]')
plt.ylabel('\u03B3 [-]')
plt.legend(ncol=4)
plt.show()
