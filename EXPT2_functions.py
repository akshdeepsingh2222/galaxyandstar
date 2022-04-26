import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# def Max_boltzmann(Vel,mass,T):
#     K=1.381e-23       # Boltzmann Constant
#     a_kg=1.66054e-27  # Conversion factor for amu to kg 
#     probab_density=[]
#     for i in range(len(Vel)):
#         z=(-0.5*mass*a_kg*(Vel[i])**(2))/(K*T)
#         y=(((mass*a_kg)/(2*np.pi*K*T))**(3/2))*(4*np.pi*(Vel[i])**(2)*np.exp(z))
#         probab_density.append(y)
#     return probab_density
def Max_boltzmann(E,T,E_f):
    K=1.381e-23       # Boltzmann Constant
    E_f = 1.38e-23   # J/K : constant:fermi_energy_level
    probab_density=[]
    for i in range(len(E)):
        z=(E[i]-E_f)/(K*T)
        y= 1/(np.exp(z))
        probab_density.append(y)
    return probab_density

def Fermi_dirac(E,T,E_f):
    K=1.381e-23       # Boltzmann Constant
    E_f = 1.38e-23   # J/K : constant:fermi_energy_level
    probab_density=[]
    for i in range(len(E)):
        z=(E[i]-E_f)/(K*T)
        y= 1/(np.exp(z)+1)
        probab_density.append(y)
    return probab_density

def Bose_einstein(E,T,E_f):
    K=1.381e-23       # Boltzmann Constant
    E_f = 1.38e-23   # J/K : constant:fermi_energy_level
    probab_density=[]
    for i in range(len(E)):
        z=(E[i]-E_f)/(K*T)
        y= 1/(np.exp(z)-1)
        probab_density.append(y)
    return probab_density


K=1.381e-23        # J/K : constant:fermi_energy_level
E_f = 1.38e-23
T=0.1*(E_f/K)
# T=1
E=np.linspace(0,5,200)*E_f
prob=Bose_einstein(E, T, E_f)
prob1=Fermi_dirac(E, T, E_f)
prob2=Max_boltzmann(E,T,E_f)
plt.plot(E,prob,label='Bose-Einstein')
# plt.plot(E,prob2,label='Maxwell Boltzmann')
plt.plot(E,prob1,label='fermidirac')
plt.legend()
plt.grid()
plt.show()
#################################################################
# mass=1.008     # in AMU
# T=10             # Average temp in Kelvin
# V=np.arange(0,1500,1)   #Vel in M/s

# prob=Max_boltzmann(V, mass, T)
# # prob1=Max_boltzmann(V, 6.941, T)
# plt.plot(V,prob)
# # plt.plot(V,prob1)
# plt.grid()
# plt.show()
# #############################



# K=1.381e-23        # J/K : constant:fermi_energy_level
# E_f = 1.38e-23
# T=0.1*(E_f/K)
# E=np.linspace(0,2*E_f,200)
# prob=Fermi_dirac(E, T, E_f)
# plt.plot(E,prob)
# plt.grid()
# plt.show()
# ##############################


# K=1.381e-23        # J/K : constant:fermi_energy_level
# E_f = 1.38e-23
# T=0.01*(E_f/K)
# E=np.linspace(0,2*E_f,200)
# prob=Bose_einstein(E, T, E_f)
# plt.plot(E,prob)
# plt.grid()
# plt.show()


