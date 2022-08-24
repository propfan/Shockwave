import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import fsolve

#Variables iniciales
#p0=rho0*R*T0
M1=5
p1=1e5
T1=300
beta_punto2=np.deg2rad(35)
#Gamma gases perfectos
gamma=1.4

#Polynomial from NASA 200-20000 K
#Values for N2
a=np.array([[2.21037150E+04,-3.81846182E+02,6.08273836E+00,-8.53091441E-03,1.38464619E-05,-9.62579362E-09,2.51970581E-12],
            [5.87712406E+05,-2.23924907E+03,6.06694922E+00,-6.13968550E-04,1.49180668E-07,-1.92310549E-11,1.06195439E-15],
            [8.31013916E+08,-6.42073354E+05,2.02026464E+02,-3.06509205E-02,2.48690333E-06,-9.70595411E-11,1.43753888E-15]])
Mm=28.013/1000

def R_g(Mm):
    return 8.31451/Mm

def Gamma(T,Mm):
    Rg=R_g(Mm)
    global cp
    if (T >= 200) & (T<=1000):
        cp=Rg*(a[0,0]*T**-2+a[0,1]*T**-1+a[0,2]+a[0,3]*T+a[0,4]*T**2+a[0,5]*T**3+a[0,6]*T**4)
    elif (T > 1000) & (T<=6000):
        cp=Rg*(a[1,0]*T**-2+a[1,1]*T**-1+a[1,2]+a[1,3]*T+a[1,4]*T**2+a[1,5]*T**3+a[1,6]*T**4)
    elif (T > 6000) & (T<=20000):
        cp=Rg*(a[2,0]*T**-2+a[2,1]*T**-1+a[2,2]+a[2,3]*T+a[2,4]*T**2+a[2,5]*T**3+a[2,6]*T**4)
    gammaT=(1-Rg/cp)**-1
    return gammaT

def cp_TtT(T,Tt,Mm):
    Rg=R_g(Mm)
    global AA
    if (Tt >= 200) & (Tt<=1000):
        AA=Rg*(-a[0,0]*(1/Tt-1/T)+a[0,1]*np.log(Tt/T)+a[0,2]*(Tt-T)+a[0,3]*(1/2)*(Tt**2-T**2)+a[0,4]*(1/3)*(Tt**3-T**3)+a[0,5]*(1/4)*(Tt**4-T**4)+a[0,6]*(1/5)*(Tt**5-T**5))    
    elif (Tt > 1000) & (Tt<=6000):
        AA=Rg*(-a[1,0]*(1/Tt-1/T)+a[1,1]*np.log(Tt/T)+a[1,2]*(Tt-T)+a[1,3]*(1/2)*(Tt**2-T**2)+a[1,4]*(1/3)*(Tt**3-T**3)+a[1,5]*(1/4)*(Tt**4-T**4)+a[1,6]*(1/5)*(Tt**5-T**5))    
    elif (Tt > 6000) & (Tt<=20000):
        AA=Rg*(-a[2,0]*(1/Tt-1/T)+a[2,1]*np.log(Tt/T)+a[2,2]*(Tt-T)+a[2,3]*(1/2)*(Tt**2-T**2)+a[2,4]*(1/3)*(Tt**3-T**3)+a[2,5]*(1/4)*(Tt**4-T**4)+a[2,6]*(1/5)*(Tt**5-T**5))    
    return AA

#------******------
#Código gases perfectos
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
#------******------
#------******------
#Contrucción de la rama positiva y negatuva onda incidente
b_p =np.linspace(np.arcsin(1/M1),beta_punto2, 200)

#Llamada función de onda incidente
P_p, Theta_d_p= funShockI(gamma,M1,b_p)

#Llamada función de onda reflejada
M2_perf, Pr_0, Theta1=funShockR(gamma, M1, beta_punto2)

#Construcción de la rama negativa onda reflejada
b_r =np.linspace(-0.825, -np.arcsin(1/M2_perf), 200)
Pr, Theta2_d=funSchock2r(gamma, M2_perf, b_r)

print(M2_perf, Pr_0, Theta1*180/np.pi)
#------******------


#Construccion de la polar incidente
def funShockI(T1,M1, beta):
    #Velocidad Tangencial
    Rg=R_g(Mm)
    g1=Gamma(T1,Mm)
    a1=(g1*Rg*T1)**0.5
    u1=M1*a1
    rho1=p1/(Rg*T1)
    w1=u1*np.cos(beta)
    
    #Velocidad Normal
    u1n=u1*np.sin(beta)
    M1n=u1n/a1
    
    #Temperatura remanso
    func_ht = lambda Tt : -u1**2+2*cp_TtT(T1, Tt,Mm)
    Tt_initial_guess = (1+M1n**2*(g1-1)*0.5)*T1
    Tt_s=fsolve(func_ht, Tt_initial_guess)
    
    #Resolución T2
    func_T2= lambda T2 : -Rg*T1/u1n-u1n+Rg*T2/(2*cp_TtT(T2,Tt_s,Mm)-w1**2)**0.5+(2*cp_TtT(T2,Tt_s,Mm)-w1**2)**0.5
    T2_initial_guess = (1+M1n**2*(g1-1)*0.5)*T1*0.9
    T2_s=fsolve(func_T2, T2_initial_guess)
    
    u2n=(2*cp_TtT(T2_s, Tt_s,Mm)-w1**2)**0.5
    u2=(u2n**2+w1**2)**0.5
    
    Theta_2rad=beta-np.arcsin(u2n/u2)
    Theta_2deg=np.rad2deg(Theta_2rad)
    
    g2=Gamma(T2_s,Mm)
    a2=(g2*Rg*T2_s)**0.5
    M2=u2/a2
    M2n=u2n/a2
    
    #Saltos
    T21=T2_s/T1
    u21n=u2n/u1n
    rho21=1/u21n
    p2=T2_s*Rg*rho1*rho21
    p21=p2/p1

    return p21, Theta_2deg, M2, T2_s, p2

#Construccion Onda reflejada
def funShockR(T1,M1,p1,p01,theta1,beta):
    #Velocidad Tangencial
    Rg=R_g(Mm)
    g1=Gamma(T1,Mm)
    a1=(g1*Rg*T1)**0.5
    u1=M1*a1
    rho1=p1/(Rg*T1)
    w1=u1*np.cos(beta)
    
    #Velocidad Normal
    u1n=u1*np.sin(beta)
    M1n=u1n/a1
    
    #Temperatura remanso
    func_ht = lambda Tt : -u1**2+2*cp_TtT(T1, Tt,Mm)
    Tt_initial_guess = (1+M1n**2*(g1-1)*0.5)*T1
    Tt_s=fsolve(func_ht, Tt_initial_guess)
    
    #Resolución T2
    func_T2= lambda T2 : -Rg*T1/u1n-u1n+Rg*T2/(2*cp_TtT(T2,Tt_s,Mm)-w1**2)**0.5+(2*cp_TtT(T2,Tt_s,Mm)-w1**2)**0.5
    T2_initial_guess = (1+M1n**2*(g1-1)*0.5)*T1*0.9
    T2_s=fsolve(func_T2, T2_initial_guess)
    
    u2n=(2*cp_TtT(T2_s, Tt_s,Mm)-w1**2)**0.5
    u2=(u2n**2+w1**2)**0.5
    
    Theta_2rad=-(beta-np.arcsin(u2n/u2))
    Theta_2deg=(theta1+np.rad2deg(Theta_2rad))
    
    g2=Gamma(T2_s,Mm)
    a2=(g2*Rg*T2_s)**0.5
    M2=u2/a2
    M2n=u2n/a2
    
    #Saltos
    T21=T2_s/T1
    u21n=u2n/u1n
    rho21=1/u21n
    p2=T2_s*Rg*rho1*rho21
    p21=p2/p1*p01

    return p21, Theta_2deg, M2, T2_s, p2

#Onda incidente
beta_in = []
p21_s= []
Theta_1deg_s=[]

for beta1 in np.arange(np.arcsin(1/M1),beta_punto2, 0.001): 
    p21, Theta_1deg, M2_no, T2_no, p2_no =funShockI(T1,M1, beta1)
    beta_in = np.append(beta_in,np.rad2deg(beta1))
    p21_s= np.append(p21_s,p21)
    Theta_1deg_s= np.append(Theta_1deg_s,Theta_1deg)
p21_beta1, Theta_1deg_point2, M2, T2, p2_beta1 = funShockI(T1,M1, beta_punto2)
print(M2,p21_beta1,Theta_1deg_point2,T2)

#Onda reflejada
beta_re = []
p32_s= []
Theta_2deg_s=[]

for beta2 in np.arange(np.arcsin(1/M2),np.deg2rad(46), 0.001): 
    p32, Theta_2deg, M2_no, T2_no, p2_no =funShockR(T2,M2,p21_beta1*p1,p21_beta1,Theta_1deg_point2, beta2)
    beta_re = np.append(beta_re,np.rad2deg(beta2))
    p32_s= np.append(p32_s,p32)
    Theta_2deg_s= np.append(Theta_2deg_s,Theta_2deg)

#Plot
#Gases calorificamente imperfectos
plt.plot(Theta_1deg_s,p21_s, color='red', linestyle='dashdot',linewidth=0.5)
plt.plot(Theta_2deg_s,p32_s, color='black', linestyle='dashdot',linewidth=0.5)

#Gases perfectos
plt.plot(Theta_d_p, P_p, color='blue', linestyle='dotted',linewidth=0.5)
plt.plot(Theta1*180/np.pi+Theta2_d, Pr*Pr_0, color='blue', linestyle='dotted',linewidth=0.5)

plt.grid(color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('\u03B8')
plt.ylabel('p/p\u2080')
plt.show()
