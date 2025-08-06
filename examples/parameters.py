# librairies
import numpy as np
import matplotlib.pyplot as plt

# function 
def f_M(Dvol, Dvap, Dsuf, Dgb, c, etas):
    '''
    Compute the mobility M from different parameters
    '''
    M_vol = Dvol*c**3*(10-15*c+6*c**2)
    M_vap = Dvap*(1-(c**3*(10-15*c+6*c**2)))
    M_suf = Dsuf*c*(1-c)
    M_gb = Dgb*etas
    return M_vol + M_vap + M_suf + M_gb

# input 
Dvol = 0.01 
Dvap = 0.001
Dsuf = 4.0
Dgb = 0.4

L_c = np.linspace(0, 1, 20)
L_etas = np.linspace(0, 1, 20)

M_M = np.zeros((len(L_c), len(L_etas)))
min_M = None

# output
for i_c in range(len(L_c)):
    for i_etas in range(len(L_etas)):
        M_M[i_c, i_etas] = f_M(Dvol, Dvap, Dsuf, Dgb, L_c[i_c], L_etas[i_etas])

        if min_M == None:
            min_M = M_M[i_c, i_etas]
            c_etas_min = (L_c[i_c], L_etas[i_etas])
            max_M = M_M[i_c, i_etas]
            c_etas_max = (L_c[i_c], L_etas[i_etas])
        else:  
            if M_M[i_c, i_etas] < min_M:
                min_M = M_M[i_c, i_etas]
                c_etas_min = (L_c[i_c], L_etas[i_etas])
            if max_M < M_M[i_c, i_etas]:
                max_M = M_M[i_c, i_etas]
                c_etas_max = (L_c[i_c], L_etas[i_etas])

print(min_M, c_etas_min)
print(max_M, c_etas_max)

im = plt.imshow(M_M)
plt.colorbar(im)
plt.show()


