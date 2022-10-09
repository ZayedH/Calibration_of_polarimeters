import numpy as np
from K_A import K_A
import Simulation as sim


def Calibration_A(M_0, M_1, M_2, M_3, B_0, B_1, B_2, B_3):
    KA = K_A(M_0, M_1, M_2, M_3, B_0, B_1, B_2, B_3)
    eigenvalues, eigenvectors = np.linalg.eigh(KA)
    v = eigenvectors[:, 0]
    v = (1/v[0])*v

    return v.reshape((4, 4))


# Application
# Matrices M de Muller
M_Air = sim.M_Air

M_Pol0 = sim.f_Polar

M_Pol90 = sim.M_Pol90

M_Ret30 = sim.M_Ret30


# Matrices B des mesures
B_Air = sim.MBN

B_Pol0 = sim.MBN_CalPol0

B_Pol90 = sim.MBN_CalPol90

B_Ret30 = sim.MBN_CalRet30


print("##### K_A eigenvalues ##########")
print(Calibration_A(M_Air, M_Pol0, M_Pol90,
      M_Ret30, B_Air, B_Pol0, B_Pol90, B_Ret30))
