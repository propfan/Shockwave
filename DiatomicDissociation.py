import numpy as np
from scipy.optimize import fsolve
from matplotlib import pyplot as plt

def f(variables):
  (P, R, T, A)= variables
  f1 = -P+(1 + A)*R*T 
  f2 = -P + 1+(7/5)*(1-1/R)*M1**2 
  f3 = -T+(6-1/R-2*A*bd-2*(1-A)*bv/(np.exp(bv/T)-1))/(2*(A+3)-(1+ A)*R)  #-T+(6-1/R-2*A*bd)/(8-R*(1+A))
  f4 = -A**2/(1-A)+B*np.exp(-bd/T)*(T**0.5)*(1-np.exp(-bv/T))/R
  #Aproximación
  #f4 = -A**2/(1-A)+R*B*np.exp(-bd/T)
  return [f1,f2,f3,f4]

#Inputs
T1=300 
p1=1E5
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

Pjump_dis=[]
Tjump_dis=[]
Rjump_dis=[]
alpha_dis=[]

Pjump_perf=[]
Tjump_perf=[]
Rjump_perf=[]
alpha_perf=[]

M1_in=[]
for M1 in np.arange(3,20, 1):

    Rg=R_cte/Mm
    rho1=p1/(Rg*T1)
    
    a1=(Rg*7/5*T1)**0.5
    bv = Tv/T1 
    bd = Td/T1 
    B=(G*m*tr*T1**0.5*(np.pi*m*kb/h**2)**(3/2))/rho1
    
    print('B=',B) 
    print('bv=',bv) 
    print('bd=',bd)
    
    (Pjump,Rjump, Tjump, alpha)=fsolve(f, ((7*M1**2-1)/(6), (6*M1**2)/(M1**2+5), (7*M1**2-1)*(M1**2+5)/(36*M1**2),0.5))
    
    alpha_r=alpha*(1-alpha)/(2-alpha)
    alpha_T=-alpha_r*(0.5+bd*(1-(1+bv/bd)*np.exp(bv/Tjump))/(Tjump*(1-np.exp(-bv/Tjump))))
    
    evib=(bv/Tjump)/(-1+np.exp(bv/Tjump))
    a21=((5*Tjump/7)*(1+alpha+alpha_r+(1+alpha+alpha_T)*(2*(1+alpha)-alpha_r*(1-evib+2*bd/Tjump))/(5+alpha+2*(1-alpha)*evib**2*np.exp(-bv/Tjump)+alpha_T*(1-2*evib+2*bd/Tjump))))**0.5
    
    M2=M1/(Rjump*a21)
   
    M1_in=np.append(M1_in,M1)
    Tjump_dis=np.append(Tjump_dis,Tjump)
    Rjump_dis=np.append(Rjump_dis,Rjump)
    alpha_dis=np.append(alpha_dis,alpha)
    Pjump_dis=np.append(Pjump_dis,Pjump)

    Tjump_perf=np.append(Tjump_perf,(7*M1**2-1)*(M1**2+5)/(36*M1**2))
    Rjump_perf=np.append(Rjump_perf,(6*M1**2)/(M1**2+5))
    Pjump_perf=np.append(Pjump_perf,(7*M1**2-1)/(6))
    
    print('M1=',M1) 
    print('M2=',M2)
    print('P21=',Pjump)
    print('T21=',Tjump)
    print('rho21=',Rjump)
    print('Alpha Dissotiation=',alpha)

#plt.plot(M1_in,Tjump_dis, color='blue', label="O\u2082 dis") 
#plt.plot(M1_in,Tjump_perf, color='orange',label="O\u2082 perfecto")

#plt.plot(M1_in,Rjump_dis, color='blue', label="O\u2082 dis") 
#plt.plot(M1_in,Rjump_perf
#plt.plot(M1_in,alpha_dis)
plt.plot(M1_in,Pjump_dis, color='blue', label="O\u2082 dis") 
plt.plot(M1_in,Pjump_perf, color='orange',label="O\u2082 perfecto")

plt.grid(color='grey', linestyle='--', linewidth=0.5)
plt.xlabel('M1 [-]')
plt.ylabel('Salto p [-]') #para alpha \u03B1  rho \u03C1
plt.legend(ncol=3)
plt.show()
