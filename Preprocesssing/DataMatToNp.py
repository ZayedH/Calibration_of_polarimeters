import MinimizeRapportEigenValues as m
import numpy as np
import sys
from scipy.io import loadmat


sys.path.insert(1, 'AandW_pixel_calibration')
import Simulation as sim


a = loadmat("Preprocesssing/D_matrix/B0.mat")
b = loadmat("Preprocesssing/D_matrix/E1_PSG.mat")
c = loadmat("Preprocesssing/D_matrix/E2_PSG.mat")
d = loadmat("Preprocesssing/D_matrix/E3_PSG.mat")
B_0 = a['B0']
B_1 = b['E1_PSG']
B_2 = c['E2_PSG']
B_3 = d['E3_PSG']

size_matrix_x = np.shape(B_0)[0]
size_matrix_y = np.shape(B_0)[1]

shape=np.shape(B_0)

A = np.zeros(shape)
W = np.zeros(shape)

##Parameters of the Polarisor oriented at 0
t_Pol0 = np.zeros((size_matrix_x,size_matrix_y))
Icp_Pol0 = np.zeros((size_matrix_x,size_matrix_y))
Ic_Pol0 = np.zeros((size_matrix_x,size_matrix_y))
Is_Pol0 = np.zeros((size_matrix_x,size_matrix_y))

##Parameters of the Polarisor oriented at 90
t_Pol90 = np.zeros((size_matrix_x,size_matrix_y))
Icp_Pol90 = np.zeros((size_matrix_x,size_matrix_y))
Ic_Pol90 = np.zeros((size_matrix_x,size_matrix_y))
Is_Pol90 = np.zeros((size_matrix_x,size_matrix_y))

##Parameters of the Retardator oriented at 30
t_Ret30 = np.zeros((size_matrix_x,size_matrix_y))
Icp_Ret30 = np.zeros((size_matrix_x,size_matrix_y))
Ic_Ret30 = np.zeros((size_matrix_x,size_matrix_y))
Is_Ret30 = np.zeros((size_matrix_x,size_matrix_y))

M_Air = sim.M_Air


sizewindow = 4  
center_x = size_matrix_x//2
center_y = size_matrix_y//2


b0_moy = np.einsum(
    'ijkl->kl', B_0[center_x-sizewindow:center_x+sizewindow, center_y -sizewindow:center_y +sizewindow]) 
b1_moy = np.einsum(
    'ijkl->kl', B_1[center_x-sizewindow:center_x+sizewindow, center_y -sizewindow:center_y +sizewindow])
b2_moy = np.einsum(
    'ijkl->kl', B_2[center_x-sizewindow:center_x+sizewindow, center_y -sizewindow:center_y +sizewindow])
b3_moy = np.einsum(
    'ijkl->kl', B_3[center_x-sizewindow:center_x+sizewindow, center_y -sizewindow:center_y +sizewindow])


M_Pol0_moy = m.ComputeMullerWithoutRotation_moy(b0_moy , b1_moy)
M_Pol90exp_moy = m.ComputeMullerWithoutRotation_moy(b0_moy , b2_moy)

M_Ret30exp_moy = m.ComputeMullerWithoutRotation_moy(b0_moy , b3_moy)


lamda_16_lamda_15=m.Find_real(90,45,M_Air,M_Pol0_moy,M_Pol90exp_moy,M_Ret30exp_moy,b0_moy,b1_moy,b2_moy,b3_moy)

np.save("A_and_W_storage/lamda_16_lamda_15.npy",lamda_16_lamda_15[2])
np.save("A_and_W_storage/thetaP_x.npy",lamda_16_lamda_15[0])
np.save("A_and_W_storage/thetaR_y.npy",lamda_16_lamda_15[1])
np.save("A_and_W_storage/couple_min.npy",lamda_16_lamda_15[3])


thetaP=lamda_16_lamda_15[3][0] 
thetaR=lamda_16_lamda_15[3][1] 



