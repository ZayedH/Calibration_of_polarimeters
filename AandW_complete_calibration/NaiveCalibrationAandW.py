import sys
import numpy as np
sys.path.insert(1, 'AandW_pixel_calibration')
from Calibration_A import Calibration_A
from Calibration_W import Calibration_W

sys.path.insert(1, 'Preprocesssing')
from DataMatToNp import *
from MinimizeRapportEigenValues import *



def calibrationAandW():
    for i in range(size_matrix_x):
        for j in range(size_matrix_y):
            B0_pixel = B_0[i, j, :]
            B1_pixel = B_1[i, j, :]
            B2_pixel = B_2[i, j, :]
            B3_pixel = B_3[i, j, :]

            if(np.linalg.det( B0_pixel)!=0):
                M_Pol0 = ComputeMullerWithoutRotation(B0_pixel , B1_pixel)
                M_Pol90 = f_Rotation(thetaP*np.pi/180)@ComputeMullerWithoutRotation(B0_pixel , B2_pixel)@f_Rotation(-thetaP*np.pi/180)
                M_Ret30 = f_Rotation(thetaR*np.pi/180)@ComputeMullerWithoutRotation(B0_pixel , B3_pixel)@f_Rotation(-thetaR*np.pi/180)

                # Make pixel calibration
                A_pixel = Calibration_A(
                    M_Air, M_Pol0, M_Pol90, M_Ret30, B0_pixel, B1_pixel, B2_pixel, B3_pixel)
                W_pixel = Calibration_W(
                    M_Air, M_Pol0, M_Pol90, M_Ret30, B0_pixel, B1_pixel, B2_pixel, B3_pixel)
                A[i, j] = A_pixel
                W[i, j] = W_pixel[0]

    return 

calibrationAandW()


np.save("A_and_W_storage/A.npy" , A)
np.save("A_and_W_storage/W.npy" , W)




