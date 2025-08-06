#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

# Own
from main import index_to_str

#-------------------------------------------------------------------------------
# Functions
#-------------------------------------------------------------------------------

def pp(dict_user):
    '''
    Main function for the post-processing.    
    '''
    # read and plot the .csv data.
    ReadCSV(dict_user)

    # initialize the dict
    dict_tempo = {
        'L_M_Neck': []
    }
    # iterate on the .pvtu files
    for i in range(dict_user['last_j']+1):

        # read the .pvtu
        ReadPVTU(dict_user, dict_tempo, 'output/vtk/PF_Sintering_other_'+index_to_str(i)+'.pvtu')

        # compute the neck area
        ComputeNeck(dict_tempo)

    # plot the neck area
    PlotNeck(dict_tempo)

#-------------------------------------------------------------------------------

def ReadCSV(dict_user):
    '''
    Read the csv file generated with Moose (made by postprocessors).
    '''
    # read file
    f = open('output/PF_Sintering_csv.csv', "r")
    lines = f.readlines()
    f.close()
    # init data
    time_pp = []
    c_pp = []
    L_eta_pp = []
    for i_grain in range(dict_user['n_eta']):
        L_eta_pp.append([])

    # iterate on lines
    for line in lines[1:]:
        line = line.replace("\n", "")
        data = line.split(',')
        # read data
        time_pp.append(float(data[0]))
        c_pp.append(float(data[1]))
        for i_grain in range(dict_user['n_eta']):
            L_eta_pp[i_grain].append(data[2+i_grain])

    # plot time-mass
    fig, ax1 = plt.subplots(1,1,figsize=(16,9))
    ax1.plot(time_pp, c_pp, linewidth=6)
    ax1.set_xlabel('time (-)', fontsize=25)
    ax1.set_ylabel('mean mass concentration (-)', fontsize=25)
    ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
    fig.tight_layout()
    fig.savefig('output/evol_time_mass.png')
    plt.close(fig)
    # output
    print('conservation of the mass:')
    print('mass min:', min(c_pp), '(delta ', int(100*(np.mean(c_pp)-min(c_pp))/np.mean(c_pp)),'% of mean value)')
    print('mass max:', max(c_pp), '(delta ', int(100*(max(c_pp)-np.mean(c_pp))/np.mean(c_pp)),'% of mean value)')

    # plot time-etas
    fig, ax1 = plt.subplots(1,1,figsize=(16,9))
    for i_grain in range(dict_user['n_eta']):
        ax1.plot(time_pp, L_eta_pp[i_grain], linewidth=6)
    ax1.set_xlabel('time (-)', fontsize=25)
    ax1.set_ylabel('mean grain concentration (-)', fontsize=25)
    ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
    fig.tight_layout()
    fig.savefig('output/evol_time_etas.png')
    plt.close(fig)

#-------------------------------------------------------------------------------

def ReadPVTU(dict_user, dict_tempo, filename):
    '''
    Read the PVTU file generated with Moose.
    '''
    # read the .pvtu file 
    dataset = pv.read(filename)

    # extract the points
    L_points = dataset.points

    # extract the concentration 
    L_c = dataset['c']
    
    # extract the etas
    L_L_eta = []
    for i_grain in range(dict_user['n_eta']):
        L_L_eta.append(dataset['eta'+str(i_grain)])

    # save
    dict_tempo['L_points'] = L_points
    dict_tempo['L_c'] = L_c
    dict_tempo['L_L_eta'] = L_L_eta

#-------------------------------------------------------------------------------

def ComputeNeck(dict_tempo):
    '''
    Compute the neck area.

    Defined as the integral of the function eta_i*eta_j.
    '''
    # initialization
    M_Neck = np.zeros((len(dict_tempo['L_L_eta']), len(dict_tempo['L_L_eta'])))
    # iterate on the pairs
    for i_grain in range(len(dict_tempo['L_L_eta'])-1):
        for j_grain in range(i_grain+1,len(dict_tempo['L_L_eta'])):
            L_eta_i_j = np.array(dict_tempo['L_L_eta'][i_grain])*np.array(dict_tempo['L_L_eta'][j_grain])
            M_Neck[i_grain, j_grain] = np.mean(L_eta_i_j)
    # save
    dict_tempo['L_M_Neck'].append(M_Neck)

#-------------------------------------------------------------------------------

def PlotNeck(dict_tempo):
    '''
    Plot the neck area.

    Defined as the integral of the function eta_i*eta_j.
    '''
    # open figure
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    # iterate on the pairs
    for i_grain in range(len(dict_tempo['L_L_eta'])-1):
        for j_grain in range(i_grain+1,len(dict_tempo['L_L_eta'])):
            L_Neck_i_j = []
            # extract the pair with time
            for ite in range(len(dict_tempo['L_M_Neck'])):
                L_Neck_i_j.append(dict_tempo['L_M_Neck'][ite][i_grain, j_grain])
            # plot
            ax1.plot(L_Neck_i_j, linewidth=6)
    # close figure
    ax1.set_xlabel('iteration (-)', fontsize=25)
    ax1.set_ylabel('neck area (-)', fontsize=25)
    ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
    fig.tight_layout()
    fig.savefig('output/evol_ite_neck.png')
    plt.close(fig)

