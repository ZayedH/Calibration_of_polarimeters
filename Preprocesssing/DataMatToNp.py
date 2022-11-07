import MinimizeRapportEigenValues as m
import numpy as np
import sys
from scipy.io import loadmat
import matplotlib.pyplot as plt
import time


sys.path.insert(1, 'AandW_pixel_calibration')
import Simulation as sim


# path = os.getcwd()
# print(path)
a = loadmat("Preprocesssing/D_matrix/B0.mat")
b = loadmat("Preprocesssing/D_matrix/E1_PSG.mat")
c = loadmat("Preprocesssing/D_matrix/E2_PSG.mat")
d = loadmat("Preprocesssing/D_matrix/E3_PSG.mat")
B_0 = a['B0']
B_1 = b['E1_PSG']
B_2 = c['E2_PSG']
B_3 = d['E3_PSG']

size_matrix = np.shape(B_0)[0]

### Extract submatrix for NaiveCalibration and Recursive###
n = 256

b0 = B_0[0:n, 0:n]
b1 = B_1[0:n, 0:n]
b2 = B_2[0:n, 0:n]
b3 = B_3[0:n, 0:n]

A = np.zeros((n, n, 4, 4))
W = np.zeros((n, n, 4, 4))

# Matrices M de Muller
M_Air = sim.M_Air

M_Pol0 = sim.f_Polar

M_Pol90 = sim.M_Pol90

M_Ret30 = sim.M_Ret30

####___Extract submatrices to do the mean over a central window___###

sizewindow =40  ## a window aroud the center
center = size_matrix//2

# normalisation is just doing nothing
b0_moy = np.einsum(
    'ijkl->kl', B_0[center-sizewindow:center+sizewindow, center-sizewindow:center+sizewindow]) 
b1_moy = np.einsum(
    'ijkl->kl', B_1[center-sizewindow:center+sizewindow, center-sizewindow:center+sizewindow])
b2_moy = np.einsum(
    'ijkl->kl', B_2[center-sizewindow:center+sizewindow, center-sizewindow:center+sizewindow])
b3_moy = np.einsum(
    'ijkl->kl', B_3[center-sizewindow:center+sizewindow, center-sizewindow:center+sizewindow])


t_Icp_Ic_Is=m.t_Icp_Ic_Is(m.computeMuller(b0_moy ,b2_moy))

t=t_Icp_Ic_Is[0]
Icp=t_Icp_Ic_Is[1]
Ic=t_Icp_Ic_Is[2]
Is=t_Icp_Ic_Is[3]

M_Ret0_real=(t/2)*np.array([[1, Icp, 0, 0],
                 [Icp, 1, 0, 0],
                 [0, 0, Ic, Is],
                 [0, 0, -Is, Ic]])

M_Pol0_real=(t/2)*np.array([[1, 1, 0, 0],
                 [1, 1, 0, 0],
                 [0, 0, 0, 0],
                 [0, 0, 0, 0]])

lamda_16_lamda_15=m.Find_real(90,30,M_Air,M_Pol0_real,M_Pol0_real,M_Ret0_real,b0_moy,b1_moy,b2_moy,b3_moy)
# plt.figure(figsize=(8,6),
#            facecolor='w')
# img= plt.imshow(lamda_16_lamda_15[2])
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
X, Y = np.meshgrid(lamda_16_lamda_15[0], lamda_16_lamda_15[1])
surf = ax.plot_surface(X, Y,lamda_16_lamda_15[2])

plt.show()

print(lamda_16_lamda_15[3])
        


        # t1 = time.perf_counter()
        # calibrationAandW()
        # t2 = time.perf_counter()
        # rec(0,n-1,0,n-1,n)
        # t3 = time.perf_counter()
        # print("naive: ",t2 - t1)
        # print("Rec: ",t3 - t2)
