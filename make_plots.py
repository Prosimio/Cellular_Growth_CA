# -*- coding: utf-8 -*-
"""
Created on Sun Jul  1 00:55:12 2018

@author: Prosimios
"""
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt

def save_obj(obj, name, folder ):
    with open(folder+'/'+ name + '.pkl', 'wb') as f:
        pkl.dump(obj, f, pkl.HIGHEST_PROTOCOL)

def load_obj(name, folder ):
    with open(folder+'/' + name + '.pkl', 'rb') as f:
        return pkl.load(f)

all_ratios = load_obj('all_ratios', 'data')

GRR = ['1:1','20:19','10:9','8:7','4:3','3:2','2:1']

mean_ratios = np.zeros(len(GRR))
std_ratios = np.zeros(len(GRR))

for j in range(len(GRR)):
    last_ratios = []
    plt.figure()
    for i in range(len(all_ratios[j])):
        
        plt.plot(all_ratios[j][i])
        last_ratios.append(all_ratios[j][i][-1])
        #mean_ratio += all_ratios[j][i][-1]
    
    plt.title('Growth rate ratio '+GRR[j])
    plt.ylabel('cell type ratio')
    plt.xlabel('check point step number')
    #plt.show()
    plt.savefig('plots/'+ 'ratio ' +str(j)+'.png', transparent = True)
    
    last_ratios = np.asarray(last_ratios)
    mean_ratios[j]= np.mean(last_ratios)
    std_ratios[j] = np.std(last_ratios)
    print(mean_ratios[j])

expected = [1/2, 20/39, 10/19, 8/15, 4/7, 3/5, 2/3]
GRR_num = [1/1,20/19,10/9,8/7,4/3,3/2,2/1]
plt.figure()
#X =np.arange(8)[1:]
X = GRR_num
plt.plot(X, mean_ratios, 'bo', label= 'observed ratio')
plt.plot(X, expected, 'ro', label = 'growth ratio')
plt.errorbar(X, mean_ratios, yerr=std_ratios)
plt.ylabel('cell type ratio')
plt.legend(loc = 'lower right')
plt.xticks([])
#plt.xticks(['1:1','10:9','4:3','3:2'])
plt.xticks(X, GRR, rotation = 80)

plt.savefig('plots/'+ 'ratios_obs_exp2.png', transparent=True)
