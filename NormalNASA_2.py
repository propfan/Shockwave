import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

#Polynomial from NASA 200-20000 K
#Values for N2
a=np.array([[2.21037150E+04,-3.81846182E+02,6.08273836E+00,-8.53091441E-03,1.38464619E-05,-9.62579362E-09,2.51970581E-12],
   [5.87712406E+05,-2.23924907E+03,6.06694922E+00,-6.13968550E-04,1.49180668E-07,-1.92310549E-11,1.06195439E-15],
   [8.31013916E+08,-6.42073354E+05,2.02026464E+02,-3.06509205E-02,2.48690333E-06,-9.70595411E-11,1.43753888E-15]])



def Gamma(T):
     #Def para N2
    Mm=28.013/1000
   
    R=8.31451
    Rg=R/Mm
    if np.where((T >= 200) & (T<=1000)):
        cp=Rg*(a[0,0]*T**-2+a[0,1]*T**-1+a[0,2]+a[0,3]*T+a[0,4]*T**2+a[0,5]*T**3+a[0,6]*T**4)
    elif np.where((T > 1000) & (T<=6000)):
        cp=Rg*(a[1,0]*T**-2+a[1,1]*T**-1+a[1,2]+a[1,3]*T+a[1,4]*T**2+a[1,5]*T**3+a[1,6]*T**4)
    elif np.where((T > 6000) & (T<=20000)):
        cp=Rg*(a[2,0]*T**-2+a[2,1]*T**-1+a[2,2]+a[2,3]*T+a[2,4]*T**2+a[2,5]*T**3+a[2,6]*T**4)

    gammaT=(1-Rg/cp)**-1
    return gammaT

def cp(T):
     #Def para N2
    Mm=28.013/1000
   
    R=8.31451
    Rg=R/Mm
    if np.where((T >= 200) & (T<=1000)):
        cp=Rg*(a[0,0]*T**-2+a[0,1]*T**-1+a[0,2]+a[0,3]*T+a[0,4]*T**2+a[0,5]*T**3+a[0,6]*T**4)
    elif np.where((T > 1000) & (T<=6000)):
        cp=Rg*(a[1,0]*T**-2+a[1,1]*T**-1+a[1,2]+a[1,3]*T+a[1,4]*T**2+a[1,5]*T**3+a[1,6]*T**4)
    elif np.where((T > 6000) & (T<=20000)):
        cp=Rg*(a[2,0]*T**-2+a[2,1]*T**-1+a[2,2]+a[2,3]*T+a[2,4]*T**2+a[2,5]*T**3+a[2,6]*T**4)
    return cp


def cp_TtT(T,Tt):
     #Def para N2
    Mm=28.013/1000
    
    R=8.31451
    Rg=R/Mm
    if np.where((Tt >= 200) & (Tt<=1000)):
        AA=2*Rg*(-a[0,0]*(1/Tt-1/T)+a[0,1]*np.log(Tt/T)+a[0,2]*(Tt-T)+a[0,3]*(1/2)*(Tt**2-T**2)+a[0,4]*(1/3)*(Tt**3-T**3)+a[0,5]*(1/4)*(Tt**4-T**4)+a[0,6]*(1/5)*(Tt**5-T**5))    
    elif np.where((Tt > 1000) & (Tt<=6000)):
        AA=2*Rg*(-a[1,0]*(1/Tt-1/T)+a[1,1]*np.log(Tt/T)+a[1,2]*(Tt-T)+a[1,3]*(1/2)*(Tt**2-T**2)+a[1,4]*(1/3)*(Tt**3-T**3)+a[1,5]*(1/4)*(Tt**4-T**4)+a[1,6]*(1/5)*(Tt**5-T**5))    
    elif np.where((Tt > 6000) & (Tt<=20000)):
        AA=2*Rg*(-a[2,0]*(1/Tt-1/T)+a[2,1]*np.log(Tt/T)+a[2,2]*(Tt-T)+a[2,3]*(1/2)*(Tt**2-T**2)+a[2,4]*(1/3)*(Tt**3-T**3)+a[2,5]*(1/4)*(Tt**4-T**4)+a[2,6]*(1/5)*(Tt**5-T**5))    
    return AA

def funShock_perfect(M):
    #Se definen parametros onda incidente
    g=1.4
    P=(7*M**2-1)/6
    T=(7*M**2-1)*(5+M**2)/(36*M**2)
    R=(6*M**2)/(M**2+5)
    M2=((M1**2+5)/(7*M**2-1))**0.5
    return P, T, R, M2

#INPUTS
T1=300
p1=1e5
#Definición para N2
Mm=28.013/1000

M2nasa=[]
gnasa=[]
Pjump_d0=[]
Tjump_d0=[]
Rjump_d0=[]
Pjump=[]
Tjump=[]
Rjump=[]
M1_in=[]
for M1 in np.arange(2,10,0.1):   
    gamma1=Gamma(T1)
    
    R=8.31451
    Rg=R/Mm
    
    a1=(gamma1*Rg*T1)**0.5
    u1=M1*a1
    rho1=p1/(Rg*T1)
    
    #Temperatura remanso
    func_ht = lambda Tt : -u1**2+ cp_TtT(T1, Tt)
    Tt_initial_guess = (1+M1**2*(gamma1-1)*0.5)*T1
    Tt_s=fsolve(func_ht, Tt_initial_guess)
    
    #Resolución T2 & u2n
    func_T2 = lambda T2 : -Rg*T1/u1-u1+Rg*T2/(cp_TtT(T2,Tt_s))**0.5+(cp_TtT(T2,Tt_s))**0.5
    T2_initial_guess = Tt_s*0.99999
    T2_s=fsolve(func_T2, T2_initial_guess)
    
    u2=(cp_TtT(T2_s,Tt_s))**0.5
    a2=(Gamma(T2_s)*Rg*T2_s)**0.5
    M2=u2/a2
       
    p21=(p1+rho1*u1**2*(1-u2/u1))/p1
    rho21=u1/u2
    T21=T2_s/T1
    
    (p_21_d0,T_21_d0,rho_21_d0,M2_d0)=funShock_perfect(M1)
    
    M1_in=np.append(M1_in,M1)
                                  
    Pjump_d0=np.append(Pjump_d0,p_21_d0)
    Tjump_d0=np.append(Tjump_d0,T_21_d0)
    Rjump_d0=np.append(Rjump_d0,rho_21_d0)
    Pjump=np.append(Pjump,p21)
    Tjump=np.append(Tjump,T21)
    Rjump=np.append(Rjump,rho21)
    
    gnasa=np.append(gnasa,gamma1)
    M2nasa=np.append(M2nasa,M2)
    
"""
#plt.plot(M1_in,Pjump_d0, color='blue', label="N\u2082 perfecto") 
#plt.plot(M1_in,Pjump, color='orange',label="N\u2082 NASA")

#plt.plot(M1_in,Rjump_d0, color='blue', label="N\u2082 perfecto") 
#plt.plot(M1_in,Rjump, color='orange',label="N\u2082 NASA")

plt.plot(M1_in,Tjump_d0, color='blue', label="N\u2082 perfecto") 
plt.plot(M1_in,Tjump, color='orange',label="N\u2082 NASA")
"""
plt.plot(M1_in,Pjump/Pjump_d0, color='black', label="p jump") 
plt.plot(M1_in,Rjump/Rjump_d0, color='blue', label="\u03C1 jump") 
plt.plot(M1_in,Tjump/Tjump_d0, color='green', label="T jump") 
plt.plot(M1_in,gnasa/1.4, color='cyan',linestyle='--', label="\u03B3") 
plt.plot(M1_in,M2nasa/M2_d0, color='red',linestyle=':', label="M\u2082") 

plt.grid(color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('M1 [-]')
plt.ylabel('NASA imperfecto/gas perfecto') #para alpha \u03B1  rho \u03C1
plt.yticks(ticks=[0,1,3,5,7,9,11], labels=[0,1,3,5,7,9,11])
plt.legend(ncol=3)
plt.show()
