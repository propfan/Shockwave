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

def cp_TtT(T,Tt):
     #Def para N2
    Mm=28.013/1000
    
    R=8.31451
    Rg=R/Mm
    if(Tt>= 200) & (Tt<=1000):
        AA=2*Rg*(-a[0,0]*(1/Tt-1/T)+a[0,1]*np.log(Tt/T)+a[0,2]*(Tt-T)+a[0,3]*(1/2)*(Tt**2-T**2)+a[0,4]*(1/3)*(Tt**3-T**3)+a[0,5]*(1/4)*(Tt**4-T**4)+a[0,6]*(1/5)*(Tt**5-T**5))    
    elif (Tt> 1000) & (Tt<=6000):
        AA=2*Rg*(-a[1,0]*(1/Tt-1/T)+a[1,1]*np.log(Tt/T)+a[1,2]*(Tt-T)+a[1,3]*(1/2)*(Tt**2-T**2)+a[1,4]*(1/3)*(Tt**3-T**3)+a[1,5]*(1/4)*(Tt**4-T**4)+a[1,6]*(1/5)*(Tt**5-T**5))    
    elif (Tt> 6000) & (Tt<=20000):
        AA=2*Rg*(-a[2,0]*(1/Tt-1/T)+a[2,1]*np.log(Tt/T)+a[2,2]*(Tt-T)+a[2,3]*(1/2)*(Tt**2-T**2)+a[2,4]*(1/3)*(Tt**3-T**3)+a[2,5]*(1/4)*(Tt**4-T**4)+a[2,6]*(1/5)*(Tt**5-T**5))    
    return AA

def funShock_imperfect(M, beta1, T1 ,p1):
    #Definición para N2
    Mm=28.013/1000
    
    #Se definen parametros onda incidente
    M1n=M1*np.sin(beta1)
    gamma1=Gamma(T1)
    
    R=8.31451
    Rg=R/Mm
    
    a1=(gamma1*Rg*T1)**0.5
    u1=M1*a1
    rho1=p1/(Rg*T1)
    u1n=u1*np.sin(beta1)
    w1=u1*np.cos(beta1)
    
    #Temperatura remanso
    func_ht = lambda Tt : -u1**2+ cp_TtT(T1, Tt)
    Tt_initial_guess = (1+M1**2*(gamma1-1)*0.5)*T1
    Tt_s=fsolve(func_ht, Tt_initial_guess)
    
    #Resolución T2 & u2n
    func_T2 = lambda T2 : -Rg*T1/u1n-u1n+Rg*T2/(cp_TtT(T2,Tt_s)-w1**2)**0.5+(cp_TtT(T2,Tt_s)-w1**2)**0.5
    T2_initial_guess= ((7*M1n**2-1)*(M1n**2+5)/(36*M1n**2))*T1
    (T2_s)=fsolve(func_T2, T2_initial_guess)
    
    
    u2=(cp_TtT(T2_s,Tt_s))**0.5
    a2=(Gamma(T2_s)*Rg*T2_s)**0.5
    M2=u2/a2
    
    u2n=(u2**2-w1**2)**0.5
    M2N=u2n/a2
    
    p21=(p1+rho1*u1n**2*(1-u2n/u1n))/p1
    
    Theta_rad=beta1-np.arctan(np.tan(beta1)*(u2n/u1n))
    Theta_deg=np.rad2deg(Theta_rad)
    return p21, Theta_deg

def funShock_perfect(M, b):
    #Se definen parametros onda incidente
    g=1.4
    Mn=M*np.sin(b)
    P=1+(2*g)*(Mn**2-1)/(g+1)
    Theta_rad=np.arctan((2*(Mn**2-1)/np.tan(b))/(M**2*(g+np.cos(2*b))+2))
    Theta_deg=np.degrees(Theta_rad)
    return P, Theta_deg

#INPUTS
M1=3.5
T1=300
p1=1e5


#Contrucción de la rama positiva y negatuva onda incidente
b_p =np.linspace(np.arcsin(1/M1),np.pi/2, 100)
b_p_imp =np.linspace(np.arcsin(1/M1),np.pi/2, 100)

#Gases calorificamente perfectos
P_p_perfect, Theta_d_p_perfect= funShock_perfect(M1,b_p)

#Gases calorificamente imperfectos
P_p_imperfect, Theta_d_p_imperfect= funShock_imperfect(M1, b_p_imp, T1 ,p1)

#Plot
plt.plot(Theta_d_p_perfect, P_p_perfect, color='red')

plt.plot(Theta_d_p_imperfect, P_p_imperfect, color='blue')

plt.grid(color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('\u03B8')
plt.ylabel('p/p\u2080')
plt.show()
