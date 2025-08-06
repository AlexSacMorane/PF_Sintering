#-------------------------------------------------------------------------------
# Librairies
#-------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import pyvista

#-------------------------------------------------------------------------------
# Fonction
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
# User
#-------------------------------------------------------------------------------

# output
template_namefile_pvtu = "input2_other_"
n_file = 39

# tracker
L_x_0 = []
L_y_0 = []
L_x_1 = []
L_y_1 = []
L_dx = []
L_dy = []
L_volume = []

#-------------------------------------------------------------------------------
# Read files
#-------------------------------------------------------------------------------

# iterate on files
for i_file in range(n_file):
    print(i_file+1, '/', n_file)

    # Read
    dataset = pyvista.read(template_namefile_pvtu+index_to_str(i_file)+'.pvtu')

    # Extract points
    points = dataset.points # shape: (N_points, 3)

    # Extract eta0 and eta1
    eta0s = dataset['eta0'] # shape: (N_points, )
    eta1s = dataset['eta1'] # shape: (N_points, )

    # Extract c
    cs = dataset['c'] # shape: (N_points, )

    # compute centers of eta0 and eta1
    center0 = np.array([0, 0, 0])
    center1 = np.array([0, 0, 0])
    # find the limits
    x_min = None
    # compute porosity
    # iterate on the points
    for i in range(len(points)):
        # compute centers
        center0 = center0 + float(eta0s[i])*np.array(points[i])
        center1 = center1 + float(eta1s[i])*np.array(points[i])
        # find limits
        if cs[i]>0.5:
            if x_min == None:
                x_min = points[i][0]
                x_max = points[i][0]
                y_min = points[i][1]
                y_max = points[i][1]
            else:
                if points[i][0] < x_min:
                    x_min = points[i][0]
                if x_max < points[i][0]:
                    x_max = points[i][0]
                if points[i][1] < y_min:
                    y_min = points[i][1]
                if y_max < points[i][1]:                
                    y_max = points[i][1]
    # compute centers
    center0 = center0/float(np.sum(eta0s))
    center1 = center1/float(np.sum(eta1s))

    # save
    L_x_0.append(center0[0])
    L_y_0.append(center0[1])
    L_x_1.append(center1[0])
    L_y_1.append(center1[1])
    L_dx.append(abs(center0[0]-center1[0]))
    L_dy.append(abs(center0[1]-center1[1]))
    L_volume.append((x_max-x_min)*(y_max-y_min))

#-------------------------------------------------------------------------------
# Plot
#-------------------------------------------------------------------------------

fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))
# absolute x coordinates
ax1.plot(L_x_0, linewidth=6)
ax1.plot(L_x_1, linewidth=6)
# absolute y coordinates
ax2.plot(L_y_0, linewidth=6)
ax2.plot(L_y_1, linewidth=6)
# close 
ax1.set_xlabel(r'iteration (-)', fontsize=25)
ax1.set_ylabel(r"x coordinate of the center (-)", fontsize=25)
ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
ax2.set_xlabel(r'iteration (-)', fontsize=25)
ax2.set_ylabel(r"y coordinate of the center (-)", fontsize=25)
ax2.tick_params(axis='both', labelsize=20, width=3, length=3)
fig.tight_layout()
fig.savefig('coordinate.png')
plt.close(fig)

fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,9))
# delta x coordinates
ax1.plot(L_dx, linewidth=6)
# delta y coordinates
ax2.plot(L_dy, linewidth=6)
# close 
ax1.set_xlabel(r'iteration (-)', fontsize=25)
ax1.set_ylabel(r"delta x coordinate of the centers (-)", fontsize=25)
ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
ax2.set_xlabel(r'iteration (-)', fontsize=25)
ax2.set_ylabel(r"delta y coordinate of the centers (-)", fontsize=25)
ax2.tick_params(axis='both', labelsize=20, width=3, length=3)
fig.tight_layout()
fig.savefig('delta_coordinate.png')
plt.close(fig)

fig, (ax1) = plt.subplots(1,1,figsize=(16,9))
ax1.plot(L_volume, linewidth=6)
ax1.set_xlabel(r'iteration (-)', fontsize=25)
ax1.set_ylabel(r"volume (-)", fontsize=25)
ax1.tick_params(axis='both', labelsize=20, width=3, length=3)
fig.tight_layout()
fig.savefig('volume.png')
plt.close(fig)
