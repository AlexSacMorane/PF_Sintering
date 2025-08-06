#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, label
from pathlib import Path
import os, shutil, time, random, skfmm, math

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def create_folder(name):
    '''
    Create a new folder.
    '''
    if Path(name).exists():
        shutil.rmtree(name)
    os.mkdir(name)

#-------------------------------------------------------------------------------

def index_to_str(j):
  '''
  An integer is converted to a float with 3 components
  '''
  if j < 10:
      j_str = '00'+str(j)
  elif 10 <= j and j < 100:
      j_str = '0'+str(j)
  else :
      j_str = str(j)
  return j_str

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

#-------------------------------------------------------------------------------

def write_input(dict_user):
    '''
    Write an input file from a template.
    '''
    print('\nwrite input')
    file_to_write = open('PF_Sintering.i','w')
    file_to_read = open('PF_Sintering_template.i','r')
    lines = file_to_read.readlines()
    file_to_read.close()

    j = 0
    for line in lines:
        j = j + 1
        if j == 4:
            line = line[:-1] + ' ' + str(len(dict_user['L_x'])) + '\n'
        if j == 5:
            line = line[:-1] + ' ' + str(len(dict_user['L_y'])) + '\n'
        if j == 7:
            line = line[:-1] + ' ' + str(dict_user['L_x'][0]) + '\n'
        if j == 8:
            line = line[:-1] + ' ' + str(dict_user['L_x'][-1]) + '\n'
        if j == 9:
            line = line[:-1] + ' ' + str(dict_user['L_y'][0]) + '\n'
        if j == 10:
            line = line[:-1] + ' ' + str(dict_user['L_y'][-1]) + '\n'
        if j == 15:
            line = ''
            for i_grain in range(len(dict_user['L_M_eta'])):
                line = line + '\t[./eta'+str(i_grain)+']\n'+\
                              '\t\torder = FIRST\n'+\
                              '\t\tfamily = LAGRANGE\n'+\
                              '\t\toutputs = exodus\n'+\
                              '\t\t[./InitialCondition]\n'+\
                              '\t\t\ttype = FunctionIC\n'+\
                              '\t\t\tfunction = eta'+str(i_grain)+'_txt\n'+\
                              '\t\t[../]\n'+\
                              '\t[../]\n'
        if j == 30:
            line = ''
            for i_grain in range(len(dict_user['L_M_eta'])):
                # prepare coupled_variables str
                coupled_variables = "'"
                for j_grain in range(len(dict_user['L_M_eta'])):
                    if j_grain != i_grain:
                        coupled_variables = coupled_variables + 'eta'+str(j_grain) + ' '
                coupled_variables = coupled_variables + "c'"
                # write line
                line = line + '\t# Order parameter eta'+str(i_grain)+'\n'+\
                              '\t[./deta'+str(i_grain)+'dt]\n'+\
                              '\t\ttype = TimeDerivative\n'+\
                              '\t\tvariable = eta'+str(i_grain)+'\n'+\
                              '\t[../]\n'+\
                              '\t[./ACBulk_eta'+str(i_grain)+']\n'+\
                              '\t\ttype = AllenCahn\n'+\
                              '\t\tvariable = eta'+str(i_grain)+'\n'+\
                              '\t\tcoupled_variables = '+coupled_variables+'\n'+\
                              '\t\tmob_name = L\n'+\
                              '\t\tf_name = f_tot\n'+\
                              '\t[../]\n'+\
                              '\t[./ACInterface_eta'+str(i_grain)+']\n'+\
                              '\t\ttype = ACInterface\n'+\
                              '\t\tvariable = eta'+str(i_grain)+'\n'+\
                              '\t\tmob_name = L\n'+\
                              '\t\tkappa_name = kappa_eta\n'+\
                              '\t[../]\n'
        if j == 37:
            # prepare coupled_variables str
            coupled_variables = "'"
            for i_grain in range(len(dict_user['L_M_eta'])):
                coupled_variables = coupled_variables + 'eta'+str(i_grain) + ' '
            coupled_variables = coupled_variables + "'"
            # write line
            line = line[:-1] + ' ' + coupled_variables + '\n'
        if j == 61:
            line = line[:-1] + "'" + str(dict_user['L']) + ' ' + str(dict_user['kappa_eta']) + ' ' +\
                                     str(dict_user['M']) + ' ' + str(dict_user['kappa_c']) + "'\n"
        if j == 67:
            # prepare coupled_variables str
            coupled_variables = "'"
            for i_grain in range(len(dict_user['L_M_eta'])):
                coupled_variables = coupled_variables + 'eta'+str(i_grain) + ' '
            coupled_variables = coupled_variables + "c'"
            # write line
            line = line[:-1] + ' ' + coupled_variables + '\n'
        if j == 69:
            line = line[:-1] + "'" + str(dict_user['A']) + ' ' + str(dict_user['B']) + "'\n"
        if j == 70:
            # prepare expression str
            expression = "'A*(c^2)*((1-c)^2) + B*(c^2 + 6*(1-c)*("
            for i_grain in range(len(dict_user['L_M_eta'])):
                expression = expression + 'eta'+str(i_grain) + '^2+'
            expression = expression[:-1] + ') - 4*(2-c)*('   
            for i_grain in range(len(dict_user['L_M_eta'])):
                expression = expression + 'eta'+str(i_grain) + '^3+' 
            expression = expression[:-1] + ') + 3*(' 
            for i_grain in range(len(dict_user['L_M_eta'])):
                expression = expression + 'eta'+str(i_grain) + '^2+' 
            expression = expression[:-1] + ")^2)'"     
            # write line
            line = line[:-1] + ' ' + expression + '\n'
        if j == 78:
            line = ''
            for i_grain in range(len(dict_user['L_M_eta'])):
                line = line + '\t[eta'+str(i_grain)+'_txt]\n'+\
                              '\t\ttype = PiecewiseMultilinear\n'+\
                              '\t\tdata_file = data/eta'+str(i_grain)+'.txt\n'+\
                              '\t[../]\n'
        if j == 104 or j == 105 or j == 108 or j == 109:
            line = line[:-1] + ' ' + str(dict_user['crit_res']) + '\n'
        if j == 112:
            line = line[:-1] + ' ' + str(dict_user['n_ite_max']) + '\n'
        if j == 116:
            line = line[:-1] + ' ' + str(dict_user['dt_PF']) + '\n'      
        if j == 121:
            line = ''
            for i_grain in range(len(dict_user['L_M_eta'])):
                line = line + '\t[eta'+str(i_grain)+'_pp]\n'+\
                              '\t\ttype = ElementAverageValue\n'+\
                              '\t\tvariable = eta'+str(i_grain)+'\n'+\
                              '\t[../]\n'
        if j == 131:
            # prepare variable str
            variable = "'"
            for i_grain in range(len(dict_user['L_M_eta'])):
                variable = variable + 'eta'+str(i_grain) + ' '
            variable = variable[:-1]+"'"  
            # write line
            line = line[:-1] + ' ' + variable + '\n'
        if j == 150:
            # prepare show str
            show = "'c_pp "
            for i_grain in range(len(dict_user['L_M_eta'])):
                show = show + 'eta'+str(i_grain) + '_pp '
            show = show[:-1]+"'"  
            # write line
            line = line[:-1] + ' ' + show + '\n'

        file_to_write.write(line)
    file_to_write.close()

#-------------------------------------------------------------------------------
# Parameter
#-------------------------------------------------------------------------------

# domain 
x_min = -30
x_max =  30
y_min = -30
y_max =  30
n_margins = 6 # number of pixels as margins on the size

# mesh
n_mesh = 300 # number of element in x direction of the mesh
             # the number of nodes is n_mesh+1
d_mesh = min(x_max-x_min, y_max-y_min)/n_mesh # size of the mesh element

# description of the PSD
n_grains = 9
mean_R = 15/2
var_R = 0.2 # -

# parameters for IC
n_steps = 30 # steeply increase the radius

# Description of the phase field variables
A_free_energy = 16 # the energy barrier value used for free energies description
B_free_energy = 1 # the energy barrier value used for free energies description
n_int = 6 # number of mesh in the interface
w_int = d_mesh*n_int # the interface thickness
kappa_eta = 0.5 # gradient coefficient for free energies (eta)
kappa_c = 1 # gradient coefficient for free energies (c)
L = 1 # Mobility value used for free energies (eta)
M = L # Mobility value used for free energies (c)

# computing information
n_proc = 4 # number of processor used
crit_res = 1e-4 # convergence criteria on residual

# PF time parameters
dt_PF = 0.01  # time step
n_ite_max = 100 # maximum number of iteration

# compute performances
tic = time.perf_counter()

#------------------------------------------------------------------------------
# create dict
#------------------------------------------------------------------------------

dict_user = {
    'x_min': x_min,
    'x_max': x_max,
    'y_min': y_min,
    'y_max': y_max,
    'n_margins': n_margins,
    'd_mesh': d_mesh,
    'n_grains': n_grains,
    'mean_R': mean_R,
    'var_R': var_R,
    'n_steps': n_steps,
    'w_int': w_int,
    'L': L,
    'kappa_eta': kappa_eta,
    'M': M,
    'kappa_c': kappa_c,
    'A': A_free_energy,
    'B': B_free_energy,
    'n_proc': n_proc,
    'crit_res': crit_res,
    'n_ite_max': n_ite_max,
    'dt_PF': dt_PF
}

#-------------------------------------------------------------------------------
# Plan simulation
#-------------------------------------------------------------------------------

create_folder('output')
create_folder('data')

#-------------------------------------------------------------------------------
# Generation and save of the microstructure
#-------------------------------------------------------------------------------

generate_microstructure(dict_user)

#-------------------------------------------------------------------------------
# write the input file
#-------------------------------------------------------------------------------

write_input(dict_user)

#-------------------------------------------------------------------------------
# run the simulation
#-------------------------------------------------------------------------------

print('\n run simulation')
os.system('mpiexec -n '+str(dict_user['n_proc'])+' ~/projects/moose/modules/phase_field/phase_field-opt -i PF_Sintering.i')


raise ValueError('stop')

#-------------------------------------------------------------------------------
# sort files
#-------------------------------------------------------------------------------

os.rename('PF_Sintering_csv.csv','output/PF_Sintering_csv.csv')
os.rename('PF_Sintering_out.e','output/PF_Sintering_out.e')
os.rename('PF_Sintering.i','output/PF_Sintering.i')
