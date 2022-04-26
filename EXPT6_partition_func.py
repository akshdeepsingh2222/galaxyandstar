import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import math as mt
K=1.3806e-23 # J/K  Boltzmann constant
h=6.626e-34  # J S planck's constant 
m=3.332e-27  # mass of H2 molecule in Kg
Na=6.022e23  #Avagadro No.
N=1*Na
def user_simpson(func,a,b,n):
    h = float((b-a)/n)
    result = (1)*(func(a)+func(b))
    for i in range(1,n,2):
        result+= 4*(func(a+i*h))
    for j in range(2,n-1,2):
        result+=2*(func(a+j*h))
    result*=h/3
    return result
def deriv(x,y):
    der=[]
    for i in range(len(x)-1):
        p=(y[i+1]-y[i])/(x[i+1]-x[i])
        der.append(p)
    return der
def Z(V,Temp):  
    mat=[]
    for i in range(len(V)): 
        mat_temp=[] # for diff value of T and for same V 
        V1=V[i]
        for j in range(len(Temp)):
            T1=Temp[j]
            def Z_integrand(x):  # x = nj
                power= ((h**2)*(x**2))/(8*m*(V1**(2/3))*K*T1)
                x=(np.pi/2)*(x**2)*(np.exp(-power))
                return x
            val=user_simpson(Z_integrand, 1e-4,1e11,100)
            mat_temp.append(mt.log(val))
        mat.append((mat_temp))
    return np.array(mat) 
def Z_analytical(V,Temp):  
    mat=[]
    for i in range(len(V)): 
        mat_temp=[] # for diff value of T and for same V 
        V1=V[i]
        for j in range(len(Temp)):
            T1=Temp[j]
            val=V1*(T1**1.5)*((2*np.pi*m*K)/(h**2))**(1.5)
            mat_temp.append(mt.log(val))
        mat.append((mat_temp))
    return np.array(mat) 

####################### func VS Temperature for constant Volume ######
T_pts=np.linspace(150,450,100)
V=np.linspace(20e-3,50e-3,4)
mat=Z(V,T_pts) #   4* 100

Temp=np.linspace(150,450,4)
V_pts=np.linspace(20e-3,50e-3,100)
mat_anal=Z_analytical(V, Temp) #########
mat_con_t=Z(V_pts,Temp).T
##################################################
for i in range(len(V)):   # ln(Z) vs T for constant Volume
    plt.plot(T_pts,mat[i],label="V="+str(V[i]))
    plt.plot(Temp,mat_anal[i],'.')
plt.grid()
plt.title("ln(Z) vs T for constant Volume")
plt.xlabel("Temp")
plt.ylabel("ln(Z(V,T))")
plt.legend()
plt.show()
mat_anal_trans=mat_anal.T
for j in range(len(Temp)):
    plt.plot(V_pts,mat_con_t[j],label="T="+str(Temp[j]))
    plt.plot(V,mat_anal_trans[j],'.')
plt.grid()
plt.title("ln(Z) vs V for constant Temperature")
plt.xlabel("Volume")
plt.ylabel("ln(Z(V,T))")
plt.legend()
plt.show()
############## PRESSURE ################
mat_con_t=Z(V,T_pts).T
P=[]
for i in range(len(T_pts)-1):  #P vs T for constant T
    P.append((N*K*T_pts[i])*np.array(deriv(V,mat_con_t[i])))
P=np.array(P)   # P[i]= P(T_i)
# print(P.size)
# np.reshape(len(T_pts)-1,len(V))
# print(P.size)
for i in range(len(Temp)):
    plt.plot(V[:len(V)-1],P[15*i],label="T="+str(T_pts[15*i]))
plt.grid()
plt.title("P vs V for constant T")
plt.xlabel("V")
plt.ylabel("Pressure")
plt.legend()
plt.show()

for j in range(len(V)-1): #P vs T for constant Volume
    plt.plot(T_pts[:len(T_pts)-1],P.T[j],label="V="+str(V[j]))
plt.grid()
plt.title("P vs T for constant Volume")
plt.xlabel("Temperature")
plt.ylabel("Pressure")
plt.legend()
plt.show()
########## INTERNAL ENERGY ###########
U_temporary=[]
# V=[20e-3]
V=np.array(V)
mat=Z(V,T_pts)
for j in range(len(V)):  
    U_temporary.append(deriv(T_pts,mat[j]))
U_temporary=np.array(U_temporary)
# print(U_temporary)
U_temporary=U_temporary.T
# U_temporary.reshape(len(T_pts)-1,len(V))
# print(U_temporary)
U=[]
T_const=N*K*np.square(T_pts)
# print(T_const)
for l in range(len(T_pts)-1):
    g=(T_const[l])*U_temporary[l]
    U.append(g)
U=np.array(U) 
U=U.T 

# U.reshape(len(V),len(T_pts)-1)
# # print(U)
for i in range(1):
    plt.plot(T_pts[:len(T_pts)-1],U[i],label="V="+str(V[i]))
plt.grid()
plt.title("U vs T for constant Volume")
plt.xlabel("Temperature")
plt.ylabel("Internal Energy")
plt.legend()
plt.show()

Slope,intercept=np.polyfit(T_pts[:len(T_pts)-1],U[0],1)
print("slope of U vs t at cons V:",Slope/N)  # Specific heat per particle


###############################ENTROPY ###########
U_trans=U.T
term1=[]
for i in range(len(T_pts)-1):
    term1.append(U_trans[i]/T_pts[i])

term1=np.array(term1)
term1=term1.T
mat=Z(V,T_pts)
term2=[]
con_val=np.log(N)
for j in range(len(V)):
    val=mat[j]
    val[:]=[number + (1-con_val) for number in val]
    term2.append(val)
term2=np.array(term2)
term2=term2*(N*K)
term2=term2[...,:-1]
S=term1+term2
for i in range(len(V)):
    plt.plot(T_pts[:len(T_pts)-1],S[i],label=str(V[i]))
plt.grid()
plt.title("S vs T for constant Volume")
plt.xlabel("Temperature")
plt.ylabel("ENTROPY")
plt.legend()
plt.show()

for i in range(len(V)):
    plt.plot(V,S.T[15*i],label=str(T_pts[15*i]))
plt.grid()
plt.title("S vs V for constant Temp")
plt.xlabel("Volume")
plt.ylabel("Entropy")
plt.legend()
plt.show()





