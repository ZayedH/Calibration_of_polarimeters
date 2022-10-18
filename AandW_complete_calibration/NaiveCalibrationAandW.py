import time
import sys
sys.path.insert(1, 'AandW_pixel_calibration')
from Calibration_A import Calibration_A
from Calibration_W import Calibration_W

sys.path.insert(1, 'Preprocesssing')
from DataMatToNp import *



def calibrationAandW():
    for i in range(n):
        for j in range(n):
            B0_pixel = b0[i, j, :]
            B1_pixel = b1[i, j, :]
            B2_pixel = b2[i, j, :]
            B3_pixel = b3[i, j, :]
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



