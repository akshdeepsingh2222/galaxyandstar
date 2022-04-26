import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

hc=1.98644568e-25  #J.m
k= 1.380649e-23 #joule per kelvin # Boltzmann constant
n=100
# b=0
# b_new=0

def user_simpson(func,a,b,n):
    h = float((b-a)/n)
    result = (1)*(func(a)+func(b))
    for i in range(1,n,2):
        result+= 4*(func(a+i*h))
    for j in range(2,n-1,2):
        result+=2*(func(a+j*h))
    result*=h/3
    return result
    
def Fp(x):
    y=(x**3)/(np.exp(x)-1)
    return y


x=np.linspace(1e-2,12,101)
z=[]
for i in range(len(x)):
    z.append(Fp(x[i]))
z=np.array(z)
plt.plot(x,z)
plt.grid()
plt.title("Fp v/s x")
plt.show()
Zo = np.amax(z)  # max of Fp
index1=np.where(z == Zo)  # index of Zo
Xp=x[index1]
print("max Fp:",Zo,"|| index of Fp max",index1,"|| corresponding X",Xp)

integration=user_simpson(Fp,1e-4,12,n)
expected_integration=integration/2
y=[]
for i in range(len(x)):
    y.append(user_simpson(Fp,1e-4,x[i],n))
count=0 
y_acc_values=[]
while True:    
    tol=0.5     # in percentage
    tol_value=(tol/100)*expected_integration
    extreme_val1= expected_integration+tol_value
    extreme_val2= expected_integration -tol_value
    if y[count]<=extreme_val1 and y[count]>= extreme_val2:
       y_acc_values.append(y[count])
       count+=1
       
    elif y[count]>extreme_val1:
        break
    else :
        # print("Value of new index:",count)
        count+=1 

y_acc_values=np.array(y_acc_values)
if y_acc_values.size == 0:
    print("no values matched !")
    print("either increase the nodal points of x or increase the tolerance ")

else:
    y_accurate_diff=y_acc_values - expected_integration
    y_min_diff=np.amin(y_accurate_diff)
    ind=np.where(y_min_diff==y_accurate_diff)
    y_most_acc=y_acc_values[ind]
    n_index=np.where(y == y_most_acc)
    new_index=int(n_index[0])
    new_perc_error=(abs(y[new_index]-expected_integration)/expected_integration)*100
            
    Xp_new=x[new_index]
    print("expected integration value(I/2):",expected_integration)  
    print("calculated integration value at",new_index,"index is " ,y[new_index])  
    print("percentage error in integration :",new_perc_error,"%")  
    print("new Xp:",Xp_new)
    
    
    
    b=float((hc)/(k*Xp))    # kelvin metre
    b_new=(hc)/(k*Xp_new)    # kelvin metre
    print("value of b:",b) 
    print("value of new b:",b_new) 
perc_err_b=(abs(b-0.003)/0.003)*100
perc_err_b_new=(abs(b_new-0.003)/0.003)*100
print("percentage error in b :",perc_err_b,"%")
print("percentage error in b new :",perc_err_b_new ,"%")
print("#################################")

exp_val=((np.pi)**4)/15
Ip=user_simpson(Fp,1e-4,np.max(x),100)
print("numerical result of Ip:",Ip)
print("expected value of Ip",exp_val) # (pi^4)/15

##########################################
def C(T):  # factor for temperature
    p=[]
    for i in range(len(T)):    
        eo=k*T[i]
        lo=hc/eo
        h=(8*np.pi*eo)/(lo**3)
        p.append(h)
    return p

def U(T):   # total energy density
    x1=C(T)
    x1=np.array(x1)
    x1=x1*(Ip)
    x1=x1.tolist()
    return x1
def F(T):   # radiant flux
    c=3e8  # m/s
    s=U(T)
    s=np.array(s)
    s=s*(c/4)
    s=s.tolist()
    return s
t=(np.arange(100,10000,500))*1000
T=t.tolist()
plt.plot(T,F(T),'ro-')
plt.grid()
plt.title("Radiant flux v/s T")
plt.ylabel("F(T)")
plt.xlabel("T")
plt.show()

log_temp=np.log(T)
log_F=np.log(F(T))

plt.plot(log_temp,log_F,"go-")
plt.title("log radiant flux v/s log temp")
plt.grid()
plt.show()

Slope,intercept=np.polyfit(log_temp,log_F,1)
print("slope of log data:",Slope)
print("intercept of log data:",intercept)

sigma=np.exp(intercept)
print("Stefan-Boltzmann Constant :",sigma)  # W m^(-2) K^(-4)

l=5.65628953367244e-08#0,13
l1=5.6622892536659206e-08#0,12