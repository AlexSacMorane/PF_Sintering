#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import random, skfmm, math

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def check_overlap(L_radius_grains, L_pos_grains):
    '''
    Determine if grains are overlapping.

    Used in Insert_Grains() function.
    '''
    overlap = False
    # real - real
    L_overlap = []
    for i_g in range(len(L_radius_grains)-1):
        # radius of grains
        radius_i = L_radius_grains[i_g]
        # position of grains
        pos_i = L_pos_grains[i_g]
        for j_g in range(i_g+1, len(L_radius_grains)):
            # radius of grains
            radius_j = L_radius_grains[j_g]
            # position of grains
            pos_j = L_pos_grains[j_g]
            # check distance
            if np.linalg.norm(pos_i-pos_j)<radius_i+radius_j:
                overlap = True
                L_overlap.append((i_g, j_g))
    return overlap, L_overlap

#-------------------------------------------------------------------------------

def generate_microstructure(dict_user):
    '''
    Insert n_grains grains in the domain. The grains are circle defined by a radius (uniform distribution).
    The position of the grains is randomly set, avoiding overlap between particules.
    A maximum number of tries is done per grain insertion.

    Map of etai and c are generated.
    '''
    # Initialize the mesh lists
    L_x = np.arange(dict_user['x_min'], dict_user['x_max'] +0.1*dict_user['d_mesh'], dict_user['d_mesh'])
    L_y = np.arange(dict_user['y_min'], dict_user['y_max'] +0.1*dict_user['d_mesh'], dict_user['d_mesh'])

    # initialize the list of grain maps
    L_M_eta = []

    # Insert grains
    L_pos_grains = []
    L_radius_grains = []
    i_grain = 0
    # check conditions
    while i_grain < dict_user['n_grains'] :
        # Random radius of the grain
        R_try = max(dict_user['mean_R']*(1+dict_user['var_R']*(random.random()-0.5)*2), dict_user['d_mesh']*5)
        # Random position of the grain center
        x_try = random.uniform(dict_user['x_min']+dict_user['n_margins']*dict_user['d_mesh']+R_try, dict_user['x_max']-dict_user['n_margins']*dict_user['d_mesh']-R_try)
        y_try = random.uniform(dict_user['y_min']+dict_user['n_margins']*dict_user['d_mesh']+R_try, dict_user['y_max']-dict_user['n_margins']*dict_user['d_mesh']-R_try)
        # Save grain
        L_pos_grains.append(np.array([x_try, y_try]))
        L_radius_grains.append(R_try)
        # initialize the map of the grain
        L_M_eta.append(np.zeros((len(L_y),len(L_x))))
        # prepare next grains
        i_grain = i_grain + 1

    # compute the configuration
    for i_steps in range(1, dict_user['n_steps']+1):
        print('Increase radius step',i_steps,'/',dict_user['n_steps'])

        # compute tempo radius at this step
        L_radius_grains_step = []
        for radius in L_radius_grains:
            L_radius_grains_step.append(radius*i_steps/dict_user['n_steps'])

        # check if there is no overlap
        overlap, L_overlap = check_overlap(L_radius_grains_step, L_pos_grains)
        while overlap:
            # save old positions
            L_pos_grains_old = L_pos_grains.copy()
            # iterate on overlap list to move problematic grains
            overlap = L_overlap[0]
            # get indices
            i_g = overlap[0]
            j_g = overlap[1]
            # get radius
            r_i = L_radius_grains_step[i_g]
            r_j = L_radius_grains_step[j_g]
            # get displacement vector (i->j)
            u_ij = np.array(L_pos_grains[j_g] - L_pos_grains[i_g])
            u_ij = u_ij/np.linalg.norm(u_ij)
            # move grain
            L_pos_grains[i_g] = L_pos_grains_old[j_g] - u_ij*(r_i+r_j)
            L_pos_grains[j_g] = L_pos_grains_old[i_g] + u_ij*(r_i+r_j)
            
            # check position of grain after displacement
            for i_grain in range(len(L_pos_grains)):
                # - x limit
                if L_pos_grains[i_grain][0] < dict_user['x_min'] + L_radius_grains_step[i_grain] + dict_user['n_margins']*dict_user['d_mesh']:
                    L_pos_grains[i_grain][0] = dict_user['x_min'] + (1+random.random()*2)*L_radius_grains_step[i_grain] + dict_user['n_margins']*dict_user['d_mesh']
                # + x limit
                if dict_user['x_max'] - L_radius_grains_step[i_grain] - dict_user['n_margins']*dict_user['d_mesh'] < L_pos_grains[i_grain][0]:
                    L_pos_grains[i_grain][0] = dict_user['x_max'] - (1+random.random()*2)*L_radius_grains_step[i_grain] - dict_user['n_margins']*dict_user['d_mesh']
                # - y limit
                if L_pos_grains[i_grain][1] < dict_user['y_min'] + L_radius_grains_step[i_grain] + dict_user['n_margins']*dict_user['d_mesh']:
                    L_pos_grains[i_grain][1] = dict_user['y_min'] + (1+random.random()*2)*L_radius_grains_step[i_grain] + dict_user['n_margins']*dict_user['d_mesh']
                # + y limit
                if dict_user['y_max'] - L_radius_grains_step[i_grain] - dict_user['n_margins']*dict_user['d_mesh'] < L_pos_grains[i_grain][1]:
                    L_pos_grains[i_grain][1] = dict_user['y_max'] - (1+random.random()*2)*L_radius_grains_step[i_grain] - dict_user['n_margins']*dict_user['d_mesh']

            # look if overlap exists
            overlap, L_overlap = check_overlap(L_radius_grains_step, L_pos_grains)

    # increase the radius to ensure at least one contact per grain
    for i_grain in range(len(L_pos_grains)-1):
        # compute the distance to the closer grain
        min_distance = None
        for j_grain in range(len(L_pos_grains)):
            if i_grain != j_grain:
                if min_distance == None:
                    distance = np.linalg.norm(np.array(L_pos_grains[i_grain]) - L_pos_grains[j_grain])
                    min_distance = distance
                    min_radius = L_radius_grains[j_grain]
                else :
                    distance = np.linalg.norm(np.array(L_pos_grains[i_grain]) - L_pos_grains[j_grain])    
                    if distance < min_distance:
                        min_distance = distance
                        min_radius = L_radius_grains[j_grain]
        # adapt the radius to ensure contact
        L_radius_grains[i_grain] = min_distance - min_radius

    print('\ncompute maps')
    # Initialize the arrays
    M_c = np.zeros((len(L_y),len(L_x)))
    M_etas_plot = np.zeros((len(L_y),len(L_x)))
    # iterate on grains
    for i_grain in range(len(L_pos_grains)):
        x_grain = L_pos_grains[i_grain][0]
        y_grain = L_pos_grains[i_grain][1]
        Center_grain = np.array([x_grain, y_grain])
        r_grain = L_radius_grains[i_grain]
        # find the nearest node of the center
        L_search = list(abs(np.array(L_x-x_grain)))
        i_x_center = L_search.index(min(L_search))
        L_search = list(abs(np.array(L_y-y_grain)))
        i_y_center = L_search.index(min(L_search))
        # compute the number of node (depending on the radius)
        n_nodes = int(r_grain/(L_x[1]-L_x[0]))+4
        for i_x in range(max(0,i_x_center-n_nodes),min(i_x_center+n_nodes+1,len(L_x))):
            for i_y in range(max(0,i_y_center-n_nodes),min(i_y_center+n_nodes+1,len(L_y))):
                x = L_x[i_x]
                y = L_y[i_y]
                Point = np.array([x, y])
                distance = np.linalg.norm(Point-Center_grain)
                # Update map etas
                if distance <= r_grain:
                    L_M_eta[i_grain][-1-i_y, i_x] = 1
                    M_c[-1-i_y, i_x] = M_c[-1-i_y, i_x] + 1
                    if M_etas_plot[-1-i_y, i_x] == 0: # do not erase data
                        M_etas_plot[-1-i_y, i_x] = i_grain+1
                elif distance >= r_grain:
                    L_M_eta[i_grain][-1-i_y, i_x] = 0

    # print maps
    fig, (ax1, ax2) = plt.subplots(2,1,figsize=(16,9))
    im = ax1.imshow(M_etas_plot, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]))
    im = ax2.imshow(M_c, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]))
    fig.tight_layout()
    fig.savefig('output/IC_maps_bin.png')
    plt.close(fig)

    # prepare map for plot
    M_etas_plot = np.zeros((len(L_y),len(L_x)))
    # iterate on the grains
    for i_grain in range(len(L_M_eta)):
        # adapt eta maps
        L_M_eta[i_grain] = L_M_eta[i_grain] - 0.5

        # compute the signed distance functions
        sd_eta = skfmm.distance(L_M_eta[i_grain], dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))
        
        # compute the phase field variables
        for i_x in range(len(L_x)):
            for i_y in range(len(L_y)):
                if sd_eta[i_y, i_x] > dict_user['w_int']/2: # inside the grain
                    L_M_eta[i_grain][i_y, i_x] = 1
                    if M_etas_plot[-1-i_y, i_x] == 0: # do not erase data
                        M_etas_plot[-1-i_y, i_x] = i_grain + 1
                elif sd_eta[i_y, i_x] < -dict_user['w_int']/2: # outside the grain
                    L_M_eta[i_grain][i_y, i_x] = 0
                else : # in the interface
                    L_M_eta[i_grain][i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_eta[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))
                    if M_etas_plot[-1-i_y, i_x] == 0: # do not erase data
                        M_etas_plot[-1-i_y, i_x] = i_grain + L_M_eta[i_grain][i_y, i_x]

    # adapt the concentration map
    M_c = M_c - 0.5

    # compute the signed distance functions
    sd_c = skfmm.distance(M_c, dx = np.array([L_x[1]-L_x[0],L_y[1]-L_y[0]]))
        
    # compute the phase field variables
    for i_x in range(len(L_x)):
        for i_y in range(len(L_y)):
            if sd_c[i_y, i_x] > dict_user['w_int']/2: # inside the grain
                M_c[i_y, i_x] = 1
            elif sd_c[i_y, i_x] < -dict_user['w_int']/2: # outside the grain
                M_c[i_y, i_x] = 0
            else : # in the interface
                M_c[i_y, i_x] = 0.5*(1+math.cos(math.pi*(-sd_c[i_y, i_x]+dict_user['w_int']/2)/(dict_user['w_int'])))

    # Plot maps
    fig, ((ax1),(ax2)) = plt.subplots(2,1,figsize=(9,25))
    # parameters
    title_fontsize = 30
    # psi
    im = ax1.imshow(M_etas_plot, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]))
    fig.colorbar(im, ax=ax1)
    ax1.set_title(r'Map of $\eta$s',fontsize = title_fontsize)
    # phi
    im = ax2.imshow(M_c, interpolation = 'nearest', extent=(L_x[0],L_x[-1],L_y[0],L_y[-1]))
    fig.colorbar(im, ax=ax2)
    ax2.set_title(r'Map of c',fontsize = title_fontsize)
    fig.savefig('output/IC_maps_pf.png')
    plt.close(fig)

    # save in dicts
    dict_user['L_x'] = L_x
    dict_user['L_y'] = L_y
    dict_user['L_M_eta'] = L_M_eta
    dict_user['M_c'] = M_c

    print('write data')
    # iterate on grains
    for i_grain in range(len(L_M_eta)):
        # Write phase variables
        file_to_write_etai = open('data/eta'+str(i_grain)+'.txt','w')
        # x
        file_to_write_etai.write('AXIS X\n')
        line = ''
        for x in dict_user['L_x']:
            line = line + str(x)+ ' '
        line = line + '\n'
        file_to_write_etai.write(line)
        # y
        file_to_write_etai.write('AXIS Y\n')
        line = ''
        for y in dict_user['L_y']:
            line = line + str(y)+ ' '
        line = line + '\n'
        file_to_write_etai.write(line)
        # data
        file_to_write_etai.write('DATA\n')
        for l in range(len(dict_user['L_y'])):
            for c in range(len(dict_user['L_x'])):
                file_to_write_etai.write(str(L_M_eta[i_grain][-1-l][c])+'\n')
        # close
        file_to_write_etai.close()

    # write the concentration
    file_to_write_c = open('data/c.txt','w')
    # x
    file_to_write_c.write('AXIS X\n')
    line = ''
    for x in dict_user['L_x']:
        line = line + str(x)+ ' '
    line = line + '\n'
    file_to_write_c.write(line)
    # y
    file_to_write_c.write('AXIS Y\n')
    line = ''
    for y in dict_user['L_y']:
        line = line + str(y)+ ' '
    line = line + '\n'
    file_to_write_c.write(line)
    # data
    file_to_write_c.write('DATA\n')
    for l in range(len(dict_user['L_y'])):
        for c in range(len(dict_user['L_x'])):
            file_to_write_c.write(str(M_c[-1-l][c])+'\n')
    # close
    file_to_write_c.close()