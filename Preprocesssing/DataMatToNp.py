
import numpy as np
from scipy.io import loadmat


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
