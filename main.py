#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import os, shutil, time

# Own 
from CreateIC import *

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
n_ite_max = 1 # maximum number of iteration

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
