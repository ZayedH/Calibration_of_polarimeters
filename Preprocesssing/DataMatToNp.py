import numpy as np
import sys
from scipy.io import loadmat

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

n = 256

b0 =B_0[0:n, 0:n]
b1 =B_1[0:n, 0:n]
b2 =B_2[0:n, 0:n]
b3 =B_3[0:n, 0:n]

A = np.zeros((n, n, 4, 4))
W = np.zeros((n, n, 4, 4))

# Matrices M de Muller
M_Air = sim.M_Air

M_Pol0 = sim.f_Polar

M_Pol90 = sim.M_Pol90

M_Ret30 = sim.M_Ret30




# t1 = time.perf_counter()
# calibrationAandW()
# t2 = time.perf_counter()
# rec(0,n-1,0,n-1,n)
# t3 = time.perf_counter()
# print("naive: ",t2 - t1)
# print("Rec: ",t3 - t2)