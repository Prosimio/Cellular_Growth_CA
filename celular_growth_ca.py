# -*- coding: utf-8 -*-
"""
Created on Sat Jun 30 10:52:30 2018

@author: Prosimios
"""

import numpy as np
import matplotlib.pyplot as plt
import time
import sys
import cgca_functions as cgca

def main():
    #sys.argv[1] = grid size value   e.g. 200
    #sys.argv[2] = number of simulation  e.g. 10 
    #sys.argv[3] = steps per simulation  e.g. 50000 
    
    #initialize input variables 
    #(in the future will be better to use a file with them as input)
    
    grid_size = 150
    sim_num = 1
    steps = 15000
    
    if len(sys.argv) > 1:
        grid_size = sys.argv[1]
    else:
        grid_size = int(input("Please enter value of grid size (it will be size [n,n]): n="))
    
    if len(sys.argv) > 2:
        sim_num = sys.argv[2]
    else:
        sim_num = int(input("Please enter the number of simulations: "))
        
    if len(sys.argv) > 3:
        steps  = sys.argv[3]
    else:
        steps  = int(input("Please enter the number of steps per simulation: "))
    
    ######################################
    ###  OTHER THINGS SET BY THE USER  ###

    prob_dist = 'contact_exp' # grid selection probability, used on nb_prob()
    
    #to save the last figure
#    filename = 'Segregation\\ratios\\image_%03d.jpg'
    save_im_final = True  #put True to save them
    filename_j = 'finalSatate_%04d.jpg'
    
    
    #to save images of a sequence
    every = 100   #every how many save images and/or computations
    save_im_step = False  #put True to save them
    filename_i = 'im_%04d.jpg'
    
    #####################################        
    
    all_ratios = []
    
    
    if save_im_final == True or save_im_step == True:
        plt.figure()
        
    ##########################
    ## Run the simulations  ##  
    
    for j in range(sim_num):
        
        time1 = time.clock()
        
        ######################################
        ###  OTHER THINGS SET BY THE USER  ###
        
        initial_pattern_choise =  1
        g_ratio= [1,1] 
        # g_ratio = [2,1] --> plasmid 1 divide twice fast than plasmid 2

        
        #####################################
        
        # create a initial empty grid
        grid = np.zeros((grid_size,grid_size))
    
        # define the initial grid and plasmid pattern
        grid = cgca.initial_pattern(grid, initial_pattern_choise)
        plasm_grid = cgca.initial_plasmids(grid)
    
        # Uncoment Show the initial state
    #    plt.figure(2)
    #    im_grid = create_image(grid, plasm_grid)
    #    plt.title('step 0')       
    #    plt.imshow(im_grid)

            
        #counter of step saved images
        count = 0
    
        #determine the growth probability of each genotype
        plas_probs = cgca.plasmid_gProb(g_ratio)
        
        ratios = []  #to store the cell type ratios
        
        
        for i  in range(steps):
    
            #select a random growing cell
            cell_pos = cgca.select_cell(grid)
            free_nb = cgca.check_nbhd(grid, cell_pos)
    
            if free_nb[0].shape[0] > 0: # go if there is a place in the neighborhood
    
                plasmids = plasm_grid[cell_pos[0], cell_pos[1],:] #get its plasmids
                c_growth = cgca.plasm_g_test(plasmids, plas_probs)
    
                if c_growth == True: 
    
                    #update its plasmids and cell state, n:new
                    n_plasmids, n_state = cgca.plasmid_update(plasmids, cell_pos)
    
                    plasm_grid[cell_pos[0], cell_pos[1],:] = n_plasmids
                    grid[cell_pos[0], cell_pos[1]] = n_state
                     #state will not be evaluated before role_divide
                     #role_divide function shouldn´t allow divition of that cell
    
                    divide_flag = cgca.role_divideFlag(n_plasmids)
    
                    #perform the divition if flag changed
                    if divide_flag != 0:
                        #assign a cell to a new position
                        free_proba = cgca.nb_prob(grid, cell_pos, prob_dist)
                        grid, nCell_pos = cgca.cell_divition(grid, cell_pos, free_nb, free_proba) 
    
                        #split the mother plasmids
                        m_plasmids, c_plasmids = cgca.divide_plasmids(n_plasmids)
    
                        #assign mother and child plasmids
                        plasm_grid[cell_pos[0], cell_pos[1],:] = m_plasmids
                        plasm_grid[nCell_pos[0], nCell_pos[1],:] = c_plasmids
            else:
                grid[cell_pos[0],cell_pos[1]] = 1
    
            #store values and images
            if i%every == 0 or i == steps-1:
                ratios.append(cgca.cell_ratio(plasm_grid))
                
                if save_im_step == True:
                    plt.title('step '+ str(i+1))
                    im_grid = cgca.create_image2(grid, plasm_grid)
                    plt.imshow(im_grid)
                    plt.savefig(filename_i%(count), transparent=True)
                    count += 1
                    
            #Plot the result
            if i == steps-1 and save_im_final == True:
                
                im_grid = cgca.create_image2(grid, plasm_grid)
                plt.imshow(im_grid)
                plt.title('step '+ str(i+1))
                plt.savefig(filename_j%(j), transparent=True)

                
        all_ratios.append(np.asarray(ratios))
        
        elapsed = time.clock() - time1
        print('elapsed time: ' + str(elapsed)+ ' sec')
    return(all_ratios)
  
if __name__== "__main__":
  main()