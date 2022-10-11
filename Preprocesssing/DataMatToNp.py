import time
import numpy as np
import sys
from scipy.io import loadmat

sys.path.insert(1, 'AandW_pixel_calibration')
from Calibration_A import *

sys.path.insert(1, 'AandW_complete_calibration')
from calibrationAandW import *


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

t1 = time.perf_counter()
A  , W = calibrationAandW(M_Air, M_Pol0, M_Pol90, M_Ret30, B_0, B_1, B_2, B_3)
t2 = time.perf_counter()
print(np.shape(A))
print(t2 - t1)