#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv
import pickle, porespy, skimage

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

    print('\nread pvtu')
    # initialize the dict
    dict_tempo = {
        'L_M_Neck': [],
        'L_freeSurface': [],
        'L_M0': [],
        'L_M1': [],
        'L_M2': [],
        'L_M3': [],
        'map': None
    }
    # iterate on the .pvtu files
    for i in range(dict_user['last_j']+1):

        # read the .pvtu
        ReadPVTU(dict_user, dict_tempo, 'output/vtk/PF_Sintering_other_'+index_to_str(i)+'.pvtu')

        # Biswas, 2018
        # compute the neck area
        ComputeNeck(dict_tempo)
        # compute the free surface area
        ComputeFreeSurface(dict_tempo)
    
        # Guevel, 2022
        RebuildMap(dict_user, dict_tempo)
        ComputeMorphometers(dict_tempo)

    # plot final map
    PlotMap(dict_user, dict_tempo)

    # Biswas, 2018
    # plot the neck area
    PlotNeck(dict_tempo)
    # plot the free surface area
    PlotFreeSurface(dict_tempo)

    # Guevel, 2022
    PlotMorphometers(dict_tempo)

#-------------------------------------------------------------------------------

def ReadCSV(dict_user):
    '''
    Read the csv file generated with Moose (made by postprocessors).
    '''
    print('\nread csv')

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
            L_eta_pp[i_grain].append(float(data[2+i_grain]))

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

#-------------------------------------------------------------------------------

def ComputeFreeSurface(dict_tempo):
    '''
    Compute the free surface area.

    Defined as the integral of the function c_int (=1 for  c in [0.45, 0.55]).
    '''
    # pp data
    L_c_int = []
    for c in dict_tempo['L_c']:
        if 0.45 <= c and c <= 0.55:
            L_c_int.append(1)
        else :
            L_c_int.append(0)
    # save
    dict_tempo['L_freeSurface'].append(np.mean(L_c_int))
    
#-------------------------------------------------------------------------------

def PlotFreeSurface(dict_tempo):
    '''
    Plot the free surface area.

    Defined as the integral of the function c_int (=1 for  c in [0.45, 0.55]).
    '''
    # open figure
    fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
    ax1.plot(dict_tempo['L_freeSurface'], linewidth=6)
    # close figure
    ax1.set_xlabel('iteration (-)', fontsize=25)
    ax1.set_ylabel('neck area (-)', fontsize=25)
    ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
    fig.tight_layout()
    fig.savefig('output/evol_ite_freesurface.png')
    plt.close(fig)

#-------------------------------------------------------------------------------

def RebuildMap(dict_user, dict_tempo):
    '''
    Rebuild the map from the output of the pyvista.
    '''
    # initialization
    M_void = np.zeros((len(dict_user['L_y']), len(dict_user['L_x'])))
    M_solid = np.zeros((len(dict_user['L_y']), len(dict_user['L_x'])))

    # the map is not know
    if dict_tempo['map'] == None:
        # init map
        map = []
        # iterate on the points
        for i_point in range(len(dict_tempo['L_points'])):
            # search node in the mesh
            L_search = list(abs(np.array(dict_user['L_x']-list(dict_tempo['L_points'][i_point])[0])))
            i_x = L_search.index(min(L_search))
            L_search = list(abs(np.array(dict_user['L_y']-list(dict_tempo['L_points'][i_point])[1])))
            i_y = L_search.index(min(L_search))
            # save map
            map.append([i_x, i_y])
            # rebuild maps
            if dict_tempo['L_c'][i_point] < 0.5:
                M_void[-1-i_y, i_x] = 1
            else:
                M_solid[-1-i_y, i_x] = 1
        # save map
        dict_tempo['map'] = map

    # the map is know
    else:
        # iterate on the points
        for i_point in range(len(dict_tempo['L_points'])):
            # read the map
            i_x = dict_tempo['map'][i_point][0]
            i_y = dict_tempo['map'][i_point][1]
            # rebuild maps
            if dict_tempo['L_c'][i_point] < 0.5:
                M_void[-1-i_y, i_x] = 1
            else:
                M_solid[-1-i_y, i_x] = 1
    
    # extract maps (exclude void on the sides)
    # find i_x_min
    i_x_min = 0
    while np.max(M_solid[:, i_x_min]) == 0:
        i_x_min = i_x_min + 1
    # find i_x_max
    i_x_max = M_solid.shape[1]-1
    while np.max(M_solid[:, i_x_max]) == 0:
        i_x_max = i_x_max - 1
    # find i_y_min
    i_y_min = 0
    while np.max(M_solid[i_y_min, :]) == 0:
        i_y_min = i_y_min + 1
    # find i_y_max
    i_y_max = M_solid.shape[0]-1
    while np.max(M_solid[i_y_max, :]) == 0:
        i_y_max = i_y_max - 1
    
    # save
    dict_tempo['M_void'] = M_void.copy()
    dict_tempo['M_void_extracted'] = M_void.copy()[i_y_min: i_y_max+1, i_x_min: i_x_max+1]
    dict_tempo['M_solid'] = M_solid.copy()
    dict_tempo['M_solid_extracted'] = M_solid.copy()[i_y_min: i_y_max+1, i_x_min: i_x_max+1]

#-------------------------------------------------------------------------------

def PlotMap(dict_user, dict_tempo):
    '''
    plot the current configuration.
    '''
    # Plot maps
    fig, (ax1) = plt.subplots(1,2,figsize=(16,9))
    # parameters
    ax1.imshow(dict_tempo['M_solid'], interpolation = 'nearest', extent=(dict_user['L_x'][0],dict_user['L_x'][-1],dict_user['L_y'][0],dict_user['L_x'][-1]))
    fig.savefig('output/map_c.png')
    plt.close(fig)


#-------------------------------------------------------------------------------

def ComputeMorphometers(dict_tempo):
    '''
    Compute the morphometers.

    see Guével, 2022
    '''
    # M0 porosity
    M0 = porespy.metrics.porosity(dict_tempo['M_void_extracted'])
    # M1 perimeter
    M1 = skimage.measure.perimeter(dict_tempo['M_solid'])
    # M2 grain size
    M2 = 0
    # M3 Euler
    M3 = skimage.measure.euler_number(dict_tempo['M_solid'])

    # save
    dict_tempo['L_M0'].append(M0)
    dict_tempo['L_M1'].append(M1)
    dict_tempo['L_M2'].append(M2)
    dict_tempo['L_M3'].append(M3)

#-------------------------------------------------------------------------------

def PlotMorphometers(dict_tempo):
    '''
    Plot the morphometers.

    see Guevel, 2022.
    '''
    # open figure
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(2*16,2*9))
    # M0
    ax1.plot(dict_tempo['L_M0'], linewidth=6)
    ax1.set_xlabel('iteration (-)', fontsize=25)
    ax1.set_ylabel('M0 porosity (-)', fontsize=25)
    ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
    # M1
    ax2.plot(dict_tempo['L_M1'], linewidth=6)
    ax2.set_xlabel('iteration (-)', fontsize=25)
    ax2.set_ylabel('M1 perimeter (-)', fontsize=25)
    ax2.tick_params(axis='both', labelsize=20, width=3, length=3)
    # M2
    ax3.plot(dict_tempo['L_M2'], linewidth=6)
    ax3.set_xlabel('iteration (-)', fontsize=25)
    ax3.set_ylabel('M2 grain size (-)', fontsize=25)
    ax3.tick_params(axis='both', labelsize=20, width=3, length=3)
    # M3
    ax4.plot(dict_tempo['L_M3'], linewidth=6)
    ax4.set_xlabel('iteration (-)', fontsize=25)
    ax4.set_ylabel('M3 euler (-)', fontsize=25)
    ax4.tick_params(axis='both', labelsize=20, width=3, length=3)
    # close
    fig.tight_layout()
    fig.savefig('output/evol_ite_morphometers.png')
    plt.close(fig)

#-------------------------------------------------------------------------------
# MAIN code
#-------------------------------------------------------------------------------

if __name__ == '__main__':

    # load dict_user
    with open('output/dict_user', 'rb') as handle:
        dict_user = pickle.load(handle)

    # call pp function
    pp(dict_user)