import sys
sys.path.insert(1, 'AandW_pixel_calibration')
from Calibration_W import Calibration_W
from Calibration_A import Calibration_A
import numpy as np




def calibrationAandW(M_0, M_1, M_2, M_3, B_0, B_1, B_2, B_3):
    n = np.shape(B_0)[0]
    A = np.zeros((n, n, 4, 4))
    W = np.zeros((n, n, 4, 4))
    for i in range(n):
        for j in range(n):
            B0_pixel = B_0[i, j, :]
            B1_pixel = B_1[i, j, :]
            B2_pixel = B_2[i, j, :]
            B3_pixel = B_3[i, j, :]
            # Make pixel calibration
            A_pixel = Calibration_A(
                M_0, M_1, M_2, M_3, B0_pixel, B1_pixel, B2_pixel, B3_pixel)
            W_pixel = Calibration_W(
                M_0, M_1, M_2, M_3, B0_pixel, B1_pixel, B2_pixel, B3_pixel)
            A[i, j, :] = A_pixel
            W[i, j, :] = W_pixel
            print('i ', i)
            print('j ', j)

    return A, W
