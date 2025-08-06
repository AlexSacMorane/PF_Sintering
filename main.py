#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    import numpy as np
    import matplotlib.pyplot as plt
    from pathlib import Path
    import os, shutil, time, pickle

    # Own 
    from CreateIC import *
    from WriteI import *
    from SortFiles import *
    from PostProcessing import *

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
# MAIN code
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # compute performances
    tic = time.perf_counter()

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
    n_grains = 10
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

    # PF time parameters
    dt_PF = 0.01  # time step
    n_ite_max = 20 # maximum number of iteration

    # computing information
    n_proc = 4 # number of processor used
    crit_res = 1e-4 # convergence criteria on residual

    # sorting files
    reduce_vtk = True # reduce or not the number of vtk
    n_vtk_max = 10 # if reduced, maximal number of vtk files

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
        'dt_PF': dt_PF,
        'reduce_vtk': reduce_vtk,
        'n_vtk_max': n_vtk_max,
    }

#-------------------------------------------------------------------------------
# Plan simulation
#-------------------------------------------------------------------------------

    create_folder('output')
    create_folder('output/vtk')
    create_folder('data')

#-------------------------------------------------------------------------------
# Generation and save of the microstructure
#-------------------------------------------------------------------------------

    generate_microstructure(dict_user)

#-------------------------------------------------------------------------------
# write the input file
#-------------------------------------------------------------------------------

    write_input_pf(dict_user)

#-------------------------------------------------------------------------------
# run the simulation
#-------------------------------------------------------------------------------

    print('\nrun simulation')
    os.system('mpiexec -n '+str(dict_user['n_proc'])+' ~/projects/moose/modules/phase_field/phase_field-opt -i PF_Sintering.i')

#-------------------------------------------------------------------------------
# sort files
#-------------------------------------------------------------------------------

    # move
    os.rename('PF_Sintering_csv.csv','output/PF_Sintering_csv.csv')
    os.rename('PF_Sintering_out.e','output/PF_Sintering_out.e')
    os.rename('PF_Sintering.i','output/PF_Sintering.i')
    Sort_vtk(dict_user)
    # remove
    shutil.rmtree('data')

#-------------------------------------------------------------------------------
# save dict
#-------------------------------------------------------------------------------

    with open('output/dict_user', 'wb') as handle:
        pickle.dump(dict_user, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # compute performances
    tac = time.perf_counter()
    hours = (tac-tic)//(60*60)
    minutes = (tac-tic - hours*60*60)//(60)
    seconds = int(tac-tic - hours*60*60 - minutes*60)
    print("\nSimulation time : "+str(hours)+" hours "+str(minutes)+" minutes "+str(seconds)+" seconds")
    print('Simulation ends')

#-------------------------------------------------------------------------------
# pp
#-------------------------------------------------------------------------------
      
    # see PostProcessing.py for the called functions
    pp(dict_user)
    