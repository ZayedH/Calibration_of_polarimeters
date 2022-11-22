import numpy as np
from K_A import K_A


def Calibration_A(M_0, M_1, M_2, M_3, B_0, B_1, B_2, B_3):
    KA = K_A(M_0, M_1, M_2, M_3, B_0, B_1, B_2, B_3)
    eigenvalues, eigenvectors = np.linalg.eigh(KA)
    v = eigenvectors[:, 0]
    # v0=np.max(v)
    # v = (1/v[0])*v

    return v.reshape((4, 4))





