import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

def Gamma(T):
     #Def para N2
    Mm=28.013/1000
    
    R=8.31451
    Rg=R/Mm
    if np.where((T > 1000) & (T<6000)):
        #FOR [1000 K,6000 K]
        a1=	5.87712406E+05
        a2=	-2.23924907E+03
        a3=	6.06694922E+00
        a4=	-6.13968550E-04
        a5=	1.49180668E-07
        a6=	-1.92310549E-11
        a7=	1.06195439E-15
        cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)
    elif np.where(T > 6000):
        #FOR >6000K
        a1=	8.31013916E+08
        a2=	-6.42073354E+05
        a3=	2.02026464E+02
        a4=	-3.06509205E-02
        a5=	2.48690333E-06
        a6=	-9.70595411E-11
        a7=	1.43753888E-15
        cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)

    elif np.where((T > 200) & (T<=1000)):
        #For <-1000K
        a1=	2.21037150E+04
        a2=	-3.81846182E+02
        a3=	6.08273836E+00
        a4=	-8.53091441E-03
        a5=	1.38464619E-05
        a6=	-9.62579362E-09
        a7=	2.51970581E-12
        cp=Rg*(a1*T**-2+a2*T**-1+a3+a4*T+a5*T**2+a6*T**3+a7*T**4)
        
    gammaT=(1-Rg/cp)**-1
    return gammaT

def cp_TtT(T,Tt):
    R=8.31451
    #Def para N2
    Mm=28.013/1000
    Rg=R/Mm
    if np.where((T > 1000) & (T<6000)):
        #FOR [1000 K,6000 K]
        a1=	5.87712406E+05
        a2=	-2.23924907E+03
        a3=	6.06694922E+00
        a4=	-6.13968550E-04
        a5=	1.49180668E-07
        a6=	-1.92310549E-11
        a7=	1.06195439E-15
        
        return Rg*(-a1*(1/Tt-1/T)+a2*np.log(Tt/T)+a3*(Tt-T)+a4*(1/2)*(Tt**2-T**2)+a5*(1/3)*(Tt**3-T**3)+a6*(1/4)*(Tt**4-T**4)+a7*(1/5)*(Tt**5-T**5))
    elif np.where(T > 6000):
        #FOR >6000K
        a1=	8.31013916E+08
        a2=	-6.42073354E+05
        a3=	2.02026464E+02
        a4=	-3.06509205E-02
        a5=	2.48690333E-06
        a6=	-9.70595411E-11
        a7=	1.43753888E-15
        
        return Rg*(-a1*(1/Tt-1/T)+a2*np.log(Tt/T)+a3*(Tt-T)+a4*(1/2)*(Tt**2-T**2)+a5*(1/3)*(Tt**3-T**3)+a6*(1/4)*(Tt**4-T**4)+a7*(1/5)*(Tt**5-T**5))  
    elif np.where((T > 200) & (T<=1000)):
        #For <-1000K
        a1=	2.21037150E+04
        a2=	-3.81846182E+02
        a3=	6.08273836E+00
        a4=	-8.53091441E-03
        a5=	1.38464619E-05
        a6=	-9.62579362E-09
        a7=	2.51970581E-12
        
        return Rg*(-a1*(1/Tt-1/T)+a2*np.log(Tt/T)+a3*(Tt-T)+a4*(1/2)*(Tt**2-T**2)+a5*(1/3)*(Tt**3-T**3)+a6*(1/4)*(Tt**4-T**4)+a7*(1/5)*(Tt**5-T**5))      


def funShock_imperfect(M, beta1, T1 ,p1):
    #Definición para N2
    Mm=28.013/1000
    
    #Se definen parametros onda incidente
    Mn=M*np.sin(beta1)
    gamma1=Gamma(T1)
    
    R=8.31451
    Rg=R/Mm
    
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

    u21n=u2n/u1n
    rho21=1/u21n
    p21=(p1+u1n**2-u2n*u1n*rho1)/p1
    
    Theta_rad=beta1-np.arcsin(u2n/u2)
    Theta_deg=np.degrees(Theta_rad)

    return p21, Theta_deg


def funShock_perfect(M, b):
    #Se definen parametros onda incidente
    g=1.4
    Mn=M*np.sin(b)
    P=1+(2*g)*(Mn**2-1)/(g+1)
    Theta_rad=np.arctan((2*(Mn**2-1)/np.tan(b))/(M**2*(g+np.cos(2*b))+2))
    Theta_deg=np.degrees(Theta_rad)
    return P, Theta_deg

#Condiciones de entrada 
M1=2
T1=300.0
p1=1e5


#Contrucción de la rama positiva y negatuva onda incidente
b_p =np.linspace(np.arcsin(1/M1),np.pi/2, 100)
b_n =np.linspace(-np.pi/2,-np.arcsin(1/M1), 100)

b_p_imp =np.linspace(0.2,np.pi/2, 300)
#b_n_imp =np.linspace(-1,0, 300)

#Llamada función de onda incidente
#Gases calorificamente perfectos
P_p_perfect, Theta_d_p_perfect= funShock_perfect(M1,b_p)
P_n_perfect, Theta_d_n_perfect= funShock_perfect(M1,b_n)

#Gases calorificamente imperfectos
P_p_imperfect, Theta_d_p_imperfect= funShock_imperfect(M1, b_p_imp, T1 ,p1)
#P_n_imperfect, Theta_d_n_imperfect= funShock_imperfect(M1, b_n_imp, T1 ,p1)

#Plot
plt.plot(Theta_d_p_perfect, P_p_perfect, color='red')
plt.plot(Theta_d_n_perfect, P_n_perfect, color='red')

plt.plot(Theta_d_p_imperfect, P_p_imperfect, color='blue')
#plt.plot(Theta_d_n_imperfect, P_n_imperfect, color='black')

plt.grid(color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('\u03B8')
plt.ylabel('p/p\u2080')
plt.show()
