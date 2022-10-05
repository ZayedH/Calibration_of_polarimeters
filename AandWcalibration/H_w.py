import numpy as np


def H_w(M1, M2, B1, B2):
    invM1 = np.linalg(M1)
    invB2 = np.linalg(B2)
