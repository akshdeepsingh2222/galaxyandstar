import numpy as np
import matplotlib.pyplot as plt

K=1.3806e-23 # J/K  Boltzmann constant
h=6.626e-34  # J S planck's constant 
c=3e8 #m/s
def G(x):   # energy density 
    g=[]
    for i in range(len(x)):
        g.append(np.pi*x[i]**2)
    return g

def F_rj(x):# energy state for 
    func=[]
    for i in range(len(x)):
        func.append(x[i]**2)
    return func

def F_p(x):   # energy state for planck
    func=[]
    for i in range(len(x)):
        y=(x[i]**3)/(np.exp(x[i])-1)
        func.append(y)
    return func

def energy_state_rj(x,t):   # energy density
    E_star=K*t
    l_star=h*c/(E_star)
    constant=(8*np.pi*E_star)/(l_star**3)
    y=np.array(F_rj(x))
    L=constant*y
    return L

def energy_state_p(x,t):
    E_star=K*t
    l_star=h*c/(E_star)
    constant=(8*np.pi*E_star)/(l_star**3)
    y=np.array(F_p(x))
    L=constant*y
    return L

###########################################
power=np.linspace(10,30,100)
v_all=10**(power)
lo=1e-10  # 1 angstrom
v_o=c/(2*lo)
x_freq=v_all/v_o
visible_range=np.linspace(4e14,8e14,100)
x_freq_visible=visible_range/v_o

plt.plot(x_freq_visible,G(x_freq_visible))
# plt.xlabel("G(x)")
# plt.ylabel("x (visible range)")
plt.title( "Reduced density of states v/s dimless freq (visible)")
plt.grid()
plt.show()

plt.plot(x_freq,G(x_freq),'o-')
plt.grid()
plt.title(" Reduced density of states v/s dimless freq(complete)")
# plt.xlabel("G(x)")  
# plt.ylabel("x (complete range)")
plt.show()

plt.plot(np.log10(x_freq),np.log10(G(x_freq)))
plt.grid()
plt.title("log G(x) vs log(x) for whole range")
plt.show()
####################################
x=np.linspace(1e-2,12,100)
plt.plot(x,F_rj(x))
plt.title("average energy RJ v/s x")
plt.grid()
plt.show()
plt.plot(x,F_p(x))
plt.title("average energy Planck v/s x")
plt.grid()
plt.show()
#####################$
T=[12,15,18]
T=np.array(T)
T=T*1e2
for i in range (len(T)):
    plt.plot(x,energy_state_rj(x, T[i]),label=str(T[i]))
plt.grid()
plt.title("energy density RJ v/s x")
plt.legend()
plt.show()    

for i in range (len(T)):
    plt.plot(x,energy_state_p(x, T[i]),label=str(T[i]))
plt.title("energy density Planck v/s x")
plt.grid()
plt.legend()
plt.show()  


for i in range (len(T)):
    nu=(x*(K*T[i]))/h
    plt.plot(nu,energy_state_rj(nu, T[i]),label=str(T[i]))
plt.grid()
plt.title("spectral energy density RJ v/s nu")
plt.legend()
plt.show()    

for i in range (len(T)):
    nu=(x*(K*T[i]))/h
    plt.plot(nu,energy_state_p(x, T[i]),label=str(T[i]))
plt.grid()
plt.title("spectral energy density Planck v/s nu")
plt.legend()
plt.show()    



