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

            
            M_Pol0 = M_Pol0_moy
            M_Pol90_ = M_Pol90exp_moy
            M_Ret30_ = M_Ret30exp_moy

            if(np.linalg.det( B0_pixel)!=0):
                # print(np.linalg.cond( B0_pixel))
                C1=ComputeCmatrix(B0_pixel,B1_pixel)
                C2=ComputeCmatrix(B0_pixel,B2_pixel)
                C3=ComputeCmatrix(B0_pixel,B3_pixel)
                if(np.linalg.cond(C1)<30):
                    M_Pol0 = ComputeMullerWithoutRotation(C1)
                t_sur2 = M_Pol0[0,0]
                t_Pol0[i,j] = 2*t_sur2
                Icp_Pol0[i,j] = M_Pol0[0,1]/(t_sur2)
                Ic_Pol0[i,j] = M_Pol0[2,2]/(t_sur2)
                Is_Pol0[i,j] = M_Pol0[2,3]/(t_sur2)
                if(np.linalg.cond(C2)<30):
                    M_Pol90_ = ComputeMullerWithoutRotation(C2)
                t_sur2 = M_Pol90_[0,0]
                t_Pol90[i,j] = 2*t_sur2
                Icp_Pol90[i,j] = M_Pol90_[0,1]/(t_sur2)
                Ic_Pol90[i,j] = M_Pol90_[2,2]/(t_sur2)
                Is_Pol90[i,j] = M_Pol90_[2,3]/(t_sur2)
                M_Pol90 = f_Rotation(thetaP*np.pi/180)@M_Pol90_@f_Rotation(-thetaP*np.pi/180)
                if(np.linalg.cond(C3)<30):
                    M_Ret30_ = ComputeMullerWithoutRotation(C3)
                t_sur2 = M_Ret30_[0,0]
                t_Ret30[i,j] = 2*t_sur2
                Icp_Ret30[i,j] = M_Ret30_[0,1]/(t_sur2)
                Ic_Ret30[i,j] = M_Ret30_[2,2]/(t_sur2)
                Is_Ret30[i,j] = M_Ret30_[2,3]/(t_sur2)
                M_Ret30 = f_Rotation(thetaR*np.pi/180)@M_Ret30_@f_Rotation(-thetaR*np.pi/180)

                # Make pixel calibration
                A_pixel = Calibration_A(
                    M_Air, M_Pol0, M_Pol90, M_Ret30, B0_pixel, B1_pixel, B2_pixel, B3_pixel)
                W_pixel = Calibration_W(
                    M_Air, M_Pol0, M_Pol90, M_Ret30, B0_pixel, B1_pixel, B2_pixel, B3_pixel)
                A[i, j] = A_pixel
                W[i, j] = W_pixel[0]

    return 

calibrationAandW()

##Save A and W
np.save("A_and_W_storage/A.npy" , A)
np.save("A_and_W_storage/W.npy" , W)

##Save parameters of the Polarisor oriented at 0
np.save("parameters_storage/t_Pol0.npy" , t_Pol0)
np.save("parameters_storage/Icp_Pol0.npy" , Icp_Pol0)
np.save("parameters_storage/Ic_Pol0.npy" , Ic_Pol0)
np.save("parameters_storage/Is_Pol0.npy" , Is_Pol0)

##Save parameters of the Polarisor oriented at 90
np.save("parameters_storage/t_Pol90.npy" , t_Pol90)
np.save("parameters_storage/Icp_Pol90.npy" , Icp_Pol90)
np.save("parameters_storage/Ic_Pol90.npy" , Ic_Pol90)
np.save("parameters_storage/Is_Pol90.npy" , Is_Pol90)

##Save parameters of the Retardator oriented at 30
np.save("parameters_storage/t_Ret30.npy" , t_Ret30)
np.save("parameters_storage/Icp_Ret30.npy" , Icp_Ret30)
np.save("parameters_storage/Ic_Ret30.npy" , Ic_Ret30)
np.save("parameters_storage/Is_Ret30.npy" , Is_Ret30)





