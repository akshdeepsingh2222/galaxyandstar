import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import math as mt
K=8.617e-5 # eV/K  Boltzmann constant
h=4.135e-15  # eV S planck's constant 
m=3.332e-27  # mass of H2 molecule in Kg
Na=6.022e23  #Avagadro No.
N=10
n_level_system=2
Ej=np.arange(0,n_level_system,1)  #in eV

def deriv(x,y):
    size=len(x)
    der=[]
    for i in range(len(x)-1):  # 98 99
        p=(y[i+1]-y[i])/(x[i+1]-x[i])
        der.append(p)
    val=(y[size-1]-y[size-2])/(x[size-1]-x[size-2])  #99  100
    der.append(val)
    return der

def Z(Temp,Ej):  
    mat=[] # for diff value of T and for same V     
    for j in range(len(Temp)):
        T=Temp[j]
        def Zustands_summe(x):  # x = Ej
            Gj=1
            sum1=0
            for i in range(len(x)):
                power=(Ej[i]/(K*T))
                sum1+=Gj*(np.exp(-power))
            return sum1
        val=Zustands_summe(Ej)
        mat.append((val))
    return np.array(mat) 

def NjbyN(Temp,Ej):
    mat=Z(Temp,Ej)
    Nj=[]
    for i in range(len(Ej)):
        Nj_const_Ej=[]
        for j in range(len(Temp)):
            val=(1/mat[j])*np.exp(-Ej[i]/(K*Temp[j]))
            Nj_const_Ej.append(val)
        Nj.append(Nj_const_Ej)
    return Nj

def internal_energy(T,Ej):
    T1=T
    U_temporary=[]
    mat=Z(T1,Ej)
    U_temporary=np.array(deriv(T1,np.log(mat)))
    U=[]
    T_const=N*K*np.square(T1)
    for l in range(len(T1)):
        g=(T_const[l])*U_temporary[l]
        U.append(g)
    U=np.array(U) 
    U=U.T 
    return U

def Free_Energy(T1,Ej):
    Part=np.log(Z(T1,Ej))
    Part=Part.T
    F=[]
    for i in range(len(Part)):
        F.append(-N*K*T1[i]*Part[i])
    F=np.array(F)
    F=F.T
    return F

def Entropy(T,Ej):
    T1=T
    F=Free_Energy(T1, Ej)
    y=deriv(T1, F)
    y=np.array(y)
    y=-y
    return y

##########################################
# TEMPERATURE RANGES 
m1=4
# m2=0.25
no_of_pts=50
T1=np.linspace(1,m1*(n_level_system/(K)),no_of_pts)
# T2=np.linspace(1,m2*(n_level_system/(K)),no_of_pts)
# T2=np.linspace(1,(n_level_system/(K)),no_of_pts)
# T_star=T1/T2
T_star=np.linspace(0,m,no_of_pts)

######################### Z ##############
Partition_fn=Z(T1,Ej)
# Partition_fn2=Z(T2,Ej)
plt.plot(T1,Partition_fn)
plt.grid(b=True,which='both',axis="both")
plt.minorticks_on()
plt.ylabel("Z")
plt.xlabel("Temperature")
plt.title("Z v/s T")
plt.legend()
plt.show()
################### Nj/N    ##########
Temp=T1
Nj=NjbyN(Temp, Ej)
for i in range(len(Ej)):
    plt.scatter(Temp,Nj[i],label="Ej="+str(Ej[i]))
plt.grid(b=True,which='both',axis="both")
plt.minorticks_on()
plt.ylabel("Nj/N")
plt.xlabel("Temperature")
plt.legend()
plt.show()
# ####################INTERNAL ENERGY ################

U=internal_energy(T1, Ej)
plt.plot(T1,U)
plt.grid(b=True,which='both',axis="both")
plt.minorticks_on()
plt.ylabel("Internal Energy")
plt.xlabel("Temperature")
plt.title("U v/s T")
plt.show()
# ########################## FREE ENERGY #############

F=Free_Energy(T1, Ej)
plt.plot(T1,F)
plt.grid(b=True,which='both',axis="both")
plt.minorticks_on()
plt.ylabel("Free energy")
plt.xlabel("Temperature")
plt.title("F v/s T")
plt.show()
# # ###################### ENTROPY ##########

y= Entropy(T1, Ej) 
print(y)
plt.plot(T1,y)
plt.grid(b=True,which='both',axis="both")
plt.minorticks_on()
plt.ylabel("Entropy")
plt.xlabel("Temperature")
plt.title("S v/s T")
plt.show()
################################
# fig, axs = plt.subplots(2)
# fig.suptitle('Plots')
# axs[0].plot(T2, Partition_fn2)
# axs[0].set_ylabel("Entropy")
# axs[0].minorticks_on()
# axs[0].grid(visible=True,axis="both",which="both")
# axs[1].plot(T1, Partition_fn)
# axs[1].set_ylabel("Entropy")
# axs[1].minorticks_on()
# axs[1].grid(visible=True,axis="both",which="both")
# axs[1].set_xlabel("Temperature")
# plt.show()
##########################Entropy##########
# def entropy(T,Ej):
#     T1=T
#     Partition=Z(T1,Ej)
#     Term1=[]
#     U=internal_energy(T1, Ej)
#     U=U.T
#     Temp=T1
#     for i in range(len(Temp)):
#         Term1.append(U[i]/Temp[i])
#     Term1=np.array(Term1)
#     Term1=Term1.T
#     Term2=N*K*(np.log(Partition)-(1-np.log(N)))
#     Entropy=(Term1+Term2)
#     return Entropy
# Entropy = entropy(T1, Ej)
# plt.plot(T1,Entropy)
# plt.grid(b=True,which='both',axis="both")
# plt.minorticks_on()
# plt.ylabel("Entropy")
# plt.xlabel("Temperature")
# plt.title("S v/s T")
# plt.show()