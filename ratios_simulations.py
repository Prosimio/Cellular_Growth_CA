# -*- coding: utf-8 -*-
"""
Created on Sat Jun 30 21:55:18 2018

@author: Prosimios
"""
import pickle as pkl

def save_obj(obj, name, folder ):
    with open(folder+'/'+ name + '.pkl', 'wb') as f:
        pkl.dump(obj, f, pkl.HIGHEST_PROTOCOL)

def load_obj(name, folder ):
    with open(folder+'/' + name + '.pkl', 'rb') as f:
        return pkl.load(f)

GRR = ['1:1','20:19','10:9','8:7','4:3','3:2','2:1']
all_sims = []

filename_path = 'D:/Dropbox/Dia dia/Sistemas complejos IIQ3763/Isaac/Semestral project/Codes/Cellular_Growth_CA/celular_growth_ca.py'
work_dir = 'D:/Dropbox/Dia dia/Sistemas complejos IIQ3763/Isaac/Semestral project/Codes/Cellular_Growth_CA'

for i in range(len(GRR)):
    
    argums = '150 100 15000 ' + GRR[i]
    runfile(filename_path, args = argums , wdir= work_dir)
    ratio_values = load_obj('ratios', 'data' )
    all_sims.append(ratio_values)

save_obj(all_sims, 'all_ratios', 'data')
