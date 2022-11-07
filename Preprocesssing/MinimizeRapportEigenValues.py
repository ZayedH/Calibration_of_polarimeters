import numpy as np
import sys
sys.path.insert(1, 'AandW_pixel_calibration')
import Calibration_W as w
# We search ThetaP for polarizer and ThetaR for retarder


def computeMuller(b0, b):
    b0inv = np.linalg.inv(b0)
    M_similar = b@b0inv
    eigenvalues, eigenvectors = np.linalg.eig(M_similar)
    return eigenvalues


def t_Icp_Ic_Is(eigenvalues):
    "We assume that it respects the theoretical form"
    real = np.real(eigenvalues)
    imag = np.imag(eigenvalues)

    t = real[0]+real[1]
    Icp = (real[1]-real[0])/t
    Ic = (real[2]+real[3])/t  # ???????????????????????
    Is = (imag[2]-imag[3])/t

    return t, Icp, Ic, Is  # ok even >=1 because there is a different in the articl


def f_Rotation(a): return np.array([[1, 0, 0, 0],
                                    [0, np.cos(2*a), np.sin(2*a), 0],
                                    [0, -np.sin(2*a), np.cos(2*a), 0],
                                    [0, 0, 0, 1]])

def Find_real(thetaP,thetaR,M_0, M_1, M_2, M_3, B_0, B_1, B_2, B_3):
    m=50
    thetaP_x=np.linspace(thetaP-10,thetaP+10,m)
    thetaR_y=np.linspace(thetaR-5,thetaR+5,m)
    lamda_16_lamda_15=np.zeros((m,m))
    min =4
    couple=0
    for i in range(m):
        for j in range(m):
            M_2_thetaP=f_Rotation(thetaP_x[i]*np.pi/180)@M_2@f_Rotation(-thetaP_x[i]*np.pi/180)
            M_3_thetaR=f_Rotation(thetaR_y[j]*np.pi/180)@M_3@f_Rotation(-thetaR_y[j]*np.pi/180)
            lamda_16_lamda_15[i][j]= w.Calibration_W(M_0, M_1, M_2_thetaP,M_3_thetaR, B_0, B_1, B_2, B_3)[1]
            if(min>lamda_16_lamda_15[i][j]):
                min=lamda_16_lamda_15[i][j]
                couple=[thetaP_x[i],thetaR_y[j]]

    return thetaP_x,thetaR_y,lamda_16_lamda_15,couple

