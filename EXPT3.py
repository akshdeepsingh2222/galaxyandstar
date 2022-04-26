import numpy as np
import matplotlib.pyplot as plt
from sympy import *
deg_freedom=3
Na=6.02214076e23   #avagadro No.
h=6.626e-34        # J Sec Planck's constant
K=1.38064852e-23   #m2 kg s-2 K-1
R=8.3145           #J/mol K
sigma=1e-5

def f(x,sigma,a):
    # x=np.around(x,5)
    h=(x-a)**2/((2*sigma))
    integrand=(np.exp(-h))/np.sqrt(2*np.pi*sigma)
    integrand =np.absolute(integrand)
    return integrand 
def user_simpson(func,a,b,n):
    h = float((b-a)/n)
    result = (1)*(func(a)+func(b))
    for i in range(1,n,2):
        result+= 4*(func(a+i*h))
    for j in range(2,n-1,2):
        result+=2*(func(a+j*h))
    result*=h/3
    return result

def G_einstein(nu):
    G=[]
    for i in range(len(nu)):
        x=deg_freedom*Na*f(nu[i],sigma,1)
        G.append(x)
    return G
def G_debye(nu):
    G=[]
    V=[]
    for i in range(len(nu)):
        x=(9*Na*(nu[i]**2))/(nu_d)
        V.append(x)
    
    for j in range(len(V)):
        if nu[j] >= 1:
          G.append(0)
        else:
            G.append(V[j])
    return G

nu_e=1e12      # s^-1 for liquid (Einstein frequency)
nu_d=1e12        # s^-1 arbitrary
x=np.linspace(0,2,100)

y=(G_einstein(x))
y=np.array(y)
y=y/(3*Na)
y1=G_debye(x)
plt.plot(x*(nu_e/nu_d),y,label="Einstein")
# plt.plot(x,y1,'.',label="Debye")
plt.xlabel("V/Vx")
plt.ylabel("G(v/vx)")
plt.title("Density of States")
plt.grid()
plt.legend()
plt.show()
#####################################
theta_e=(h*nu_e)/K
theta_d=(h*nu_d)/K
def Cv_einstein(x):
    x=1/x
    y=[]
    for i in range(len(x)):
        z=((np.exp(x[i]))/(np.square(np.exp(x[i])-1)))*x[i]**2
        y.append(z)
    return y

def Cv_integrand(x):
    y=(x**3)/(np.exp(x)-1)
    return y
def Cv_debye(y):
    H=[]
    for k in range(len(y)):
        x=1/y[k]
        term1=(-3*x)/(np.exp(x)-1)
        term2=(12/(x**3))
        term3=user_simpson(Cv_integrand, 1e-10, x,1000)
        Cv=(term1+term2*term3)
        H.append(Cv)
    return H
def dulong_petit(x):
    y=[]
    for i in range(len(x)):
        y.append(1)
    return y
T=Symbol('T')
def U_einstein():
    b=3*Na*h*nu_e
    a=(h*nu_e)/K
    # a=1
    # b=1
    func=b/(exp(a/T)-1)
    z=diff(func,T)*(1/(3*R))
    return z
V=Symbol('v')
def U_debye():
    a=9*Na*h
    b=nu_d
    c=h/K
    j=(a/b**3)*(V**3/(exp(c*(V/T))-1))
    integ=integrate(j,(V,0,nu_d))
    diff_debye=diff(integ,T)
    return diff
L=U_einstein()
P=U_debye()
print("Total internal energy in case of Einstein distribution:")
print(L)
print("Total internal energy in case of Debye distribution:")
print(P)
U_e=[]
U_d=[]
for i in range(len(x)):
    U_e.append(L.subs(T,x[i]*theta_e))
    # U_d.append(P.subs(T,x[i]))
U_e=np.array(U_e)
U_d=np.array(U_d)
x=np.linspace(0,2,100)
plt.plot(x,Cv_einstein(x),label="Einstein")
plt.plot(x,Cv_debye(x),label="Debye")
plt.plot(x,U_e,'r.',label="diff(U_e)")
plt.plot(x,dulong_petit(x),label="Dulong Petit")
plt.xlabel("T/theta_x")
plt.ylabel("Cv/3R")
plt.legend()
plt.grid()
plt.show()

# plt.plot(x,np.real(U_d))
# plt.show()




