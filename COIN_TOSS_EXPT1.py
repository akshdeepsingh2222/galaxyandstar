import numpy as np
import pandas as pd
import random
import matplotlib.pyplot as plt
from scipy.stats import binom
def coin_toss(Nc,Nt):
    coin=['HEAD','TAIL']
    all_coin=[]
    for k in range(Nc):
        toss_one_coin=[]
        for i in range(Nt):
            x=random.choice(coin)
            toss_one_coin.append(x)
        all_coin.append(toss_one_coin)
    
    
    all_coin = np.array(all_coin)
    all_coin=all_coin.T
    # print(all_coin)
    all_coin = all_coin.tolist()
    occurrences=[]
    for i in range(Nt):
        occurrences.append(all_coin[i].count("HEAD"))
    
    # print(occurrences)
    possible_outcome=np.arange(Nc+1)
    # print(possible_outcome)
    frequency=[]
    for j in range(len(possible_outcome)):
        y=possible_outcome[j]
        count=0
        for l in range(len(occurrences)):
            if y == occurrences[l]:
                count+=1
        frequency.append(count)
    
    # print(frequency)
    frequency=np.array(frequency)
    occurrences=np.array(occurrences)
    probability=frequency/Nt
    return possible_outcome,occurrences,probability

def binom_dist(n, p):
    prob_dist = []
    for i in range(n+1):
        prob_dist.append(binom.pmf(i, n, p))

    return prob_dist

########## SAME NO. OF COINS AND DIFF NO. OF TRIALS    ###############
Nc=13
trials=[10*Nc,100*Nc,500*Nc,1000*Nc,10000*Nc]
xvals=np.arange(Nc+1)
# trials=np.arange(10,10000,1000)
# print(trials)
binomial_dist = binom_dist(Nc, 0.5)
for i in range(len(trials)):
    x,occ,y=coin_toss(Nc, trials[i])
    plt.plot(x,y,'.-',label=trials[i])
plt.xlabel("no. of possible outcome")
plt.ylabel("Probability")
plt.title("for different no. of trials and same no. of coins")

plt.plot(xvals, binomial_dist, label = "Binomial Distribution", color = "black", zorder = 10)
plt.grid()
plt.legend()
plt.show()

########## SAME NO. OF TRIALS AND DIFF NO. OF COINS    ###############

Nc=np.arange(1,10,2)
trials=1000

for i in range(len(Nc)):
    x,occ,y=coin_toss(Nc[i], trials)
    plt.plot(x,y,'.-',label=Nc[i])
plt.xlabel("no. of possible outcome")
plt.ylabel("Probability")
plt.title("for different no. of coins and same no. of trials")
plt.grid()
plt.legend()
plt.show()
#################### CUMULATIVE PROBABILITY ##################

Nt=10000
trial=np.arange(1,Nt+1,1)
Nc=4
poss,occurrence,prob = coin_toss(Nc, Nt)
cumm_prob=[]
cumm_outcome=[]
temp=0
for j in range(len(occurrence)):
    temp1=occurrence[j]+temp
    temp=temp1
    cumm_outcome.append(temp1)
cumm_outcome=np.array(cumm_outcome)
cumm_total=np.arange(1,Nt+1,1)

cumm_total=cumm_total*Nc
# print(cumm_total)
prob=cumm_outcome/cumm_total
prob_tails=1-prob
# print(prob_tails)
plt.plot(trial,prob,label='HEADS')
plt.plot(trial,prob_tails,label='TAILS')
plt.title("cumulative prob V/S trials upto 10000")
plt.legend()
plt.grid()
plt.show()
############################################################
plt.plot(trial[:500],prob[:500],label='HEADS')
plt.plot(trial[:500],prob_tails[:500],label='TAILS')
plt.title("cumulative prob V/S trials upto 500")
plt.legend()
plt.grid()
plt.show()
############################
# slope=[]
# for i in range(len(prob)-1):
#     x2=prob[i+1]
#     x1=prob[i]
#     y2=trial[i+1]
#     y1=trial[i]
#     slope.append((y2-y1)/(x2-x1))    

# plt.plot(trial[:Nt-1],slope)
# plt.show()
# #############################################################
prob=np.array(prob)
mean=np.mean(prob)
fluc= abs((prob) - (mean))
# y=[]
# for i in range(Nt):
#     y.append(1/(trial[i]**(0.75)))
    
plt.plot(trial[:Nt-1500],fluc[:Nt-1500])
# plt.plot(trial,y)
plt.grid()
plt.show()

























