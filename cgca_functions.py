# -*- coding: utf-8 -*-
"""
Created on Sat Jun 30 10:21:50 2018

@author: Prosimios
"""
import numpy as np
import matplotlib.pyplot as plt
import time
import matplotlib

#modify some matplotlib parameters to manage the images for illustrator
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

def show_grid(grid_array):
    plt.figure()
    plt.imshow(grid_array, cmap=plt.cm.gray)
    plt.show()
    
def select_cell(grid):
    # 2 = growing cell
    # 1 = stationary cell
    # 0 = empty space
    
    #this rule makes each cell divide only one time per step
    g_index =  np.nonzero(grid == 2) # growth index = where cell value == 2
    
    #choose a random cell in the dividing state
    index_pos = int(np.random.rand(1)[0]*g_index[0].shape[0])
    m = g_index[0][index_pos]
    n = g_index[1][index_pos]

    
    #save the cell grid index positions
    cell_index = [m,n]
        
    return(cell_index)

def check_nbhd(grid, cell_index):
    #chek free spaces in the neighbourhood
    # fs: array
    #    index of free spaces in the neighborhood
    m = cell_index[0]
    n = cell_index[1]
    
    #define the neighborhood
    nb = grid[m-1:m+2,n-1:n+2]  #nb = neighborhood
    
    #define the free spaces in nb
    fs = np.where(nb == 0)
    
    return(fs)
    
def nb_prob(grid, cell_index, prob_dist = 'contact_linear'):
    # assign division probabilities based on empty space cell contacts
    
    # prob_dist: uniform - contact_linear - contact_exp
    #    contact linear is the default or if another thing is written
    
    # fs: array
    #    index of free spaces in the neighborhood
    # return
    # prob: list
    #    list with the [0,1] probability partition limit of each free space
    #    e.g. prob = [0.23, 0.81, 1] --> second cell has bigger probability
    
    m = cell_index[0]
    n = cell_index[1]
    
    #define neighborhood
    nb = grid[m-1:m+2,n-1:n+2]  #nb = neighborhood
    
    #define the free spaces in nb
    fs = np.where(nb == 0)
    #define cell spaces in bn
    cs = np.where(nb != 0)
    
    fs_num = len(fs[0])
    prob = np.zeros(fs_num)
    contacts = np.zeros(fs_num)
    
    if prob_dist != 'uniform':
        
        # if prob_dist is something different from the options, contact_linear is the default
        for i in range(fs_num):
            mg = m + fs[0][i] - 1  #-1 to convert [0 1 2] to [-1 0 1]
            ng = n + fs[1][i] - 1

            i_nb = grid[mg-1:mg+2,ng-1:ng+2] # i position neighborhood

            occup = np.where(i_nb != 0)
            contacts[i] = len(occup[0])  #save the number of contacts of this position

        if prob_dist == 'contact_exp':
            contacts = np.exp(contacts)
            
    else:
        contacts = np.ones(fs_num)  #assign uniform values
        
    total = sum(contacts)

    prob[0] = (contacts[0]/total)

    for i in range(1,fs_num):
        prob[i] = prob[i-1]+contacts[i]/total
        
    return(prob)


def cell_divition_uniform(grid, cell_index, fs):
    # uniform neighborhood divition probability
    #fs: free neighborhood spaces
    
    m = cell_index[0]
    n = cell_index[1]
    
    if len(fs[0]) == 1:
        grid[m,n] = 1     # then, that position will not divide again

    #grown over an empty position
    #new_pos = int(np.random.rand(1)[0]*fs[0].shape[0])
    new_pos = int(np.random.rand(1)[0]*fs[0].shape[0]) #new pos in the neighbour matrix
    m_new = m + fs[0][new_pos] - 1  #-1 to convert [0 1 2] to [-1 0 1]
    n_new = n + fs[1][new_pos] - 1
    grid[m_new, n_new] = 2     # crates the new cell
        
    ncell_index = [m_new ,n_new]
    
    return(grid, ncell_index)

def cell_divition(grid, cell_index, fs, fs_proba):
    #fs: free neighborhood spaces
    #fs_proba: free spaces growth probabilities
    
    m = cell_index[0]
    n = cell_index[1]
    
    if len(fs[0]) == 1:
        grid[m,n] = 1     # then, that position will not divide again

    #grown over an empty position

    rand_val = np.random.rand(1)[0]
    # find the first position which is bigger than rand_val
    new_pos = np.where( (fs_proba > rand_val) == True )[0][0]  #new pos in the neighbour matrix
   
    m_new = m + fs[0][new_pos] - 1  #-1 to convert [0 1 2] to [-1 0 1]
    n_new = n + fs[1][new_pos] - 1
    grid[m_new, n_new] = 2     # crates the new cell
        
    ncell_index = [m_new ,n_new]
    
    return(grid, ncell_index)

def initial_plasmids(grid, pattern_num = 0, num_plas = 2, max_copy = 4):
    # grid: initial grid
    c_index =  np.nonzero(grid) # c_index, 
    #cell_number = c_index[0].shape[0]
    
    gs = grid.shape
    pattern = np.zeros((gs[0],gs[1],max_copy)) #initialize the pattern array
    
    # add different patterns
    if pattern_num == 0:  #random plasmid pattern
        for i in range(c_index[0].shape[0]): #assign a random plasmid pattern to each cell position
            pattern[c_index[0][i],c_index[1][i],:] = ((num_plas +1 )*np.random.rand(max_copy)).astype(int) 
                                                            #num_plas +1 to add "no-plasmid" state
            
    elif pattern_num == 1:
        pattern = np.ones((grid.shape))
        
    return(pattern)
    
def role_divideFlag(plasmids):
    #plasmids: cell plasmids vector
    
    max_plasmids = plasmids.shape[0]
    num_plasmids = np.nonzero(plasmids)[0].shape[0]
    divisor = max_plasmids*1.1  #arbitrary defined to make division(max_plasmids number) < 1
       
    # make a cuadratic function of probabilities
    probability = (num_plasmids/divisor)**2 
    #if a cell has no plasmids --> will not divide
    
    if np.random.rand(1) < probability:
        return(1) # divide
    else:
        return(0) # not divide
    
    #Probability tables
    #plasmid_nums = np.arange(max_plasmids +1)
    #probability = (plasmid_nums/divisor)**2 

def create_image(grid, plasgrid):
    im_s = plasgrid.shape
    aux_imR = np.zeros((im_s[0],im_s[1],im_s[2]))
    aux_imG = np.zeros((im_s[0],im_s[1],im_s[2]))
    
    for i in range(im_s[2]):
        aux_imR[:,:,i] = 1*(plasgrid[:,:,i]==1)
        aux_imG[:,:,i] = 1*(plasgrid[:,:,i]==2) 
    
    aux_imR = np.sum(aux_imR,axis=2)
    aux_imG = np.sum(aux_imG,axis=2)
    aux_transparency = 0.5*(grid[:,:]==1) + 1*(grid[:,:]==2)
    
    # create the image
    im_grid = np.zeros((im_s[0],im_s[1],im_s[2]))
    im_grid[:,:,0] = np.multiply(np.ones((im_s[0],im_s[1])),aux_imR)
    im_grid[:,:,1] = np.multiply(np.ones((im_s[0],im_s[1])),aux_imG)
    #im_grid[:,:,2] = np.ones((100,100))*250
    im_grid[:,:,3] = np.multiply(np.ones((im_s[0],im_s[1])),aux_transparency)
     # stationary cell -> transparency = 0.5)
    
    return(im_grid)

def create_image2(grid, plasgrid):
    im_s = plasgrid.shape
    aux_imR = np.zeros((im_s[0],im_s[1],im_s[2]))
    aux_imG = np.zeros((im_s[0],im_s[1],im_s[2]))
    
    for i in range(im_s[2]):
        aux_imR[:,:,i] = 1*(plasgrid[:,:,i]==1)
        aux_imG[:,:,i] = 1*(plasgrid[:,:,i]==2) 
    
    aux_imR = np.multiply(1*(np.sum(aux_imR,axis=2)>0),50*(grid[:,:]==1)) + 1*(np.sum(aux_imR,axis=2)>0)
    aux_imG = np.multiply(1*(np.sum(aux_imG,axis=2)>0),50*(grid[:,:]==1)) + 1*(np.sum(aux_imG,axis=2)>0)
    
    # create the image
    im_grid = np.zeros((im_s[0],im_s[1],3))
    im_grid[:,:,0] = np.multiply(np.ones((im_s[0],im_s[1])),aux_imR)
    im_grid[:,:,1] = np.multiply(np.ones((im_s[0],im_s[1])),aux_imG)
    #im_grid[:,:,2] = np.ones((100,100))*250
    
    return(im_grid)
    
def plasmid_gProb(g_ratio=[1,1], p_types = [1,2]):
    #define a growtn probability (=velocity) based on the plasmids
    #g_ratio: ratio of growth rate between genotypes (i.e. plasmids)
    #p_types: plasmids types or labels
    
    #built the probability class vector
    cat_len = len(g_ratio)
    probs = np.zeros(cat_len)
    
    denominator = sum(g_ratio)

    probs[0] = g_ratio[0]/denominator

    for i in range(1,cat_len):
        probs[i] = probs[i-1]+g_ratio[i]/denominator
    
    return(probs)
    
    
def plasm_g_test(plasmids,probs):
        #perform the probability test
    rand_val = np.random.rand(1)    #random value
    growth = False
    
    ptypes = np.unique(plasmids[np.where(plasmids>0)]) #plasmid types in the vector
    
    if ptypes.size == 1:  #if has one type of plasmid
        
        #determine the first position > random value
        pos = np.where( (probs > rand_val) == True )[0][0]
        
        growth_type = pos + 1  #to transform position to the corresponding plasmid type
        #found = np.where(plasmids == ptype)[0]
        
        
        #ptotal = np.where(plasmids > 0)[0]  #total num plasmids
        #pdif =  ptotal.size - found.size  
        
        if growth_type == ptypes:
            growth = True
        #if found.size>0:
        
    else:       #if has more than one type of plasmid
        
        mean_prob = probs[-1]/probs.size  #mean plasmid probability
        
        if rand_val < mean_prob:
            growth = True
        
        
    return(growth)

def cell_ratio(plasmgrid, ptype = [1,2]):
    
    c_num_plasm = np.sum(plasmgrid>0, axis=2) #number of plasmids in each grid

    plasm_sum = np.sum(plasmgrid, axis = 2)
    divition = np.divide(plasm_sum,c_num_plasm)
    
    #total = np.sum(np.isnan(divition) == False, axis = (0,1)) #it include cells with mix plasmids
    found = np.zeros(len(ptype))
    total = 0
    for i in range(len(ptype)):
        found[i] = len(np.where(divition == ptype[i])[0])
        total += found[i]
    
    ratio = found[0]/total
        
    return(ratio)

def plasmid_update(plasmids, pos_index):
    #plasmids: vector with plasmids. e.g [0,1,1,0,2]
    state = 2 # cell state = growing state
    
    plasmids_pos = np.nonzero(plasmids)
    empty_pos = np.where(plasmids == 0)
    
    num_plas = plasmids_pos[0].shape[0]
    
    if num_plas == 0:
        #it means no plasmid in the cell
        state = 1   #to not evaluate this cell in the loop again
        
    elif num_plas == plasmids.shape[0]:
        #it means all plasmids positions are full
        return(plasmids, state)
    
    else:
        copied_pos = np.random.randint(num_plas)
        plasmids[empty_pos[0][0]] = plasmids[plasmids_pos[0][copied_pos]]
        #copy the plasmid in the first free space
    
    return(plasmids, state)   

def divide_plasmids(plasmids):
    #plasmids: cell plasmids
    p_size = plasmids.size
    
    mother_p = np.zeros(p_size)
    child_p = np.zeros(p_size)
    
    np.random.shuffle(plasmids) #shuffle the plasmids
    
    if (p_size & 1) == 1: #odd case
        
        #sum a random value to choose which cell keep more plasmids
        rand_val = np.random.rand(1)
        half_p = int(p_size/2 + rand_val) 
        
    else:  #even case
        half_p = int(p_size/2)
    
    mother_p[:half_p] = plasmids[:half_p]
    child_p[half_p:]= plasmids[half_p:] 
    
    return(mother_p, child_p)
    
def initial_pattern(grid, pattern_num):
    
    pattern = {}  #initiate initial pattern dictionary
    
    # add different patterns
    pattern[0] = np.array([[2]])
    pattern[1] = np.array([[0, 0, 2, 0, 0],[0,2,2,2,0],[2,2,1,2,2],[0,2,2,2,0],[0,0,2,0,0]])
    pattern[2] = np.ones((2,35))*2
    
    #make elements which are not in the border to be = 1
    fixed_pat = pattern[pattern_num]
    
    
    #put the pattern in the grid
    gs = grid.shape
    
    m0 = int(gs[0]/2)
    n0 = int(gs[1]/2)
    
    ps = fixed_pat.shape
    mpm = int(ps[0]/2)
    npm = int(ps[1]/2)
            
    for i in range(ps[0]):
        for j in range(ps[1]):
            m = m0 + (i - mpm)
            n = n0 + (j - npm)
            
            grid[m,n] = fixed_pat[i,j]


    return(grid)



