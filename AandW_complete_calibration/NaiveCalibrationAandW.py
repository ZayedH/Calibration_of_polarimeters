import time
import sys
from AandW_pixel_calibration.Calibration_A import M_Pol0
sys.path.insert(1, 'AandW_pixel_calibration')
from Calibration_A import Calibration_A
from Calibration_W import Calibration_W

sys.path.insert(1, 'Preprocesssing')
from DataMatToNp import *
from MinimizeRapportEigenValues import *



def calibrationAandW():
    for i in range(n):
        for j in range(n):
            b0_pixel = b0[i, j, :]
            b1_pixel = b1[i, j, :]
            b2_pixel = b2[i, j, :]
            b3_pixel = b3[i, j, :]

            ##Calcul des matrices de Muller
            M_Pol0 = ComputeMullerWithoutRotation(b0_pixel , b1_pixel)
            M_Pol0 = f_Rotation(thetaP[i]*np.pi/180)@ComputeMullerWithoutRotation(b0_pixel , b2_pixel)@f_Rotation(-thetaP[i]*np.pi/180)
            M_Pol0 = f_Rotation(thetaR[i]*np.pi/180)@ComputeMullerWithoutRotation(b0_pixel , b3_pixel)@f_Rotation(-thetaR[i]*np.pi/180)

            # Make pixel calibration
            A_pixel = Calibration_A(
                M_Air, M_Pol0, M_Pol90, M_Ret30, B0_pixel, B1_pixel, B2_pixel, B3_pixel)
            W_pixel = Calibration_W(
                M_Air, M_Pol0, M_Pol90, M_Ret30, B0_pixel, B1_pixel, B2_pixel, B3_pixel)
            A[i, j] = A_pixel
            W[i, j] = W_pixel
            # print('i ', i)
            # print('j ', j)

    return 


t1 = time.perf_counter()
calibrationAandW()
t2 = time.perf_counter()
print("naive: ",t2 - t1)



