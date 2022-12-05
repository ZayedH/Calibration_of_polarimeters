import numpy as np
from K_W import K_W


def Calibration_W(M_0, M_1, M_2, M_3, B_0, B_1, B_2, B_3):
    KW = K_W(M_0, M_1, M_2, M_3, B_0, B_1, B_2, B_3)
    # print("Condition Value :",np.linalg.cond( KW))
    eigenvalues, eigenvectors = np.linalg.eigh(KW)
    v = eigenvectors[:, 0]
    # v0=np.max(v)
    # v = (1/v0)*v
    # print("eigenvalues :" ,eigenvalues)
    lamda_16_lamda_15=eigenvalues[0]/eigenvalues[1]
    # print("eigenvalues normalized :" ,lamda_16_lamda_15)

    return [v.reshape((4, 4)),lamda_16_lamda_15]


