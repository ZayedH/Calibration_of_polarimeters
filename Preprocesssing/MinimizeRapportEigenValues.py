import numpy as np
import sys
sys.path.insert(1, 'AandW_pixel_calibration')
import Calibration_W as w



def ComputeEigenvaluesCmatrix(b0, b):
    b0inv = np.linalg.inv(b0)
    M_similar = b0inv@b  
    eigenvalues, eigenvectors = np.linalg.eig(M_similar)

    return eigenvalues

def OrganizeEigenSamp(Eig, Lim):
    """ Input arguments:
    # - Eig: eigenvalues to sort
    # - Lim: tolerance of sorting. Above the limit any inequality is ignored
    #        and considered as a noise effect rather than a physical effect in
    #        the sample
    # + orEig: sorted eigenvalues

    # From the ECM theory we have, under certain sort, we have that:          # 
    # mu1 = tau*(1+cos(2*psi))           mu2 = tau*(1-cos(2*psi))           #
    # mu3 = 2*tau*sin(2*Psi)*exp(1i*del) mu4 = conj(mu(3))                  #
    # This program tries to mantain a fixed sort of the eigenvalues in order  #
    # to ensure the validity of the above mentioned equations                 #                
    """
    # Utility definitions
    N      = np.size(Eig) # number of eigenvalues to analize
    orEig  = np.zeros((N,1) , dtype=complex)  # predefinition of output argument
    conEig = np.conj(Eig)   # complex conjugate of the eigenvalues
    imEig  = abs(np.imag(Eig))   # imaginary part of the eigenvalues
    Ni     = sum(imEig > 0)   # Number of complex eigenvalues. Tipically 
                            # there are two (non nul retardance case) or 
                            # any (null retardance case). Anyhow due to 
                            # noise different cases may appear                 
                            
    equals =    False       # preassumption that there's no equalities 
    i      = 0          # analised eigenvalues counter 

    while(i <= N-2 and (not equals)): # loop to seek the related eigenvalues
        aux = Eig[i] # Reference
        j   = i+1
        while(j <= N-1 and (not equals)): # Loop to compare the reference with the rest 
                                # of eigenvalues. It stops once the a seeking
                                # condition is satisfied or all the 
                                # eigenvalues are done. The seeking condition
                                # deps on the amount of complex  
                                # eigenvalues.
                                
            # Stablishment of the seeking condition
            if(Ni == 0) : # All eigenvalues are real
                cond = abs(aux - Eig[j]) < Lim # Two equal eigevalues under
                                                # certain tolerance

            elif(Ni == 2) : # Two complex eigenvalues
                cond = abs(aux - conEig[j]) == 0  # Two complex conjugate 
                                                    # eigenvalues
            
            elif(Ni == 4) :
                cond1 = abs(aux - conEig[j]) == 0      # Two complex conjugate
                cond2 = abs(np.imag(aux)/np.real(aux)) > Lim # Apreciable retardance
                cond = cond1 and cond2

            if(cond): # Condition satisfied
                if(np.imag(Eig[i]) >= 0) : # Following convention
                                    # of positive retardance)
                    orEig[3] = Eig[i]              
                    orEig[2] = Eig[j]               

                elif(np.imag(Eig[i]) < 0) :
                    orEig[2] = Eig[i]
                    orEig[3] = Eig[j]
                equals = True # Eigenvalues 3 and 4 found
                # Delete 3 and 4. Two remaining eigenvalues to sort
                Eig = np.delete(Eig , (i,j))
                
            j = j+1
        i = i+1

    if(equals) : # sort the two remaining eigenvalues
        if(Eig[1] >= Eig[0]) : # convention: psi between 0 y pi/4
            orEig[0] = Eig[0]
            orEig[1] = Eig[1]
        else:
            orEig[0] = Eig[1]
            orEig[1] = Eig[0]

    elif(Ni == 0) : # All eigenvalues are real but no relation was found            
        orEig = np.sort(Eig)[::-1] # Sort from mayor to minor
        
    return orEig.reshape((N,1))
    
             

def Compute_t_Icp_Ic_Is(eigenvalues):
    "We assume that it respects the theoretical form"
    # max=np.max(np.abs(eigenvalues)) 
    t = np.real(eigenvalues[0][0] + eigenvalues[1][0])
    Icp = np.real((eigenvalues[1][0]  - eigenvalues[0][0]))/t # We have to fix a sign for Icp
    Ic = np.real((eigenvalues[2][0] + eigenvalues[3][0]))/t 
    Is = np.imag((eigenvalues[3][0] - eigenvalues[2][0]))/t      # We have to fix a sign for Is
    max=np.max(np.abs([t,Icp,Ic,Is]))
    return t/max, Icp/max, Ic/max, Is/max  # ok even >=1 because there is a different in the article
    
def Compute_t_Icp_Ic_Is_moy(eigenvalues):
    "We assume that it respects the theoretical form"
    # max=np.max(np.abs(eigenvalues)) 
    t = np.real(eigenvalues[0][0] + eigenvalues[1][0])
    Icp = np.real((eigenvalues[1][0]  - eigenvalues[0][0]))/t # We have to fix a sign for Icp
    Ic = np.real((eigenvalues[2][0] + eigenvalues[3][0]))/t 
    Is = np.imag((eigenvalues[3][0] - eigenvalues[2][0]))/t      # We have to fix a sign for Is
    max=np.max(np.abs([t,Icp,Ic,Is]))
    return t, Icp, Ic, Is  # ok even >=1 because there is a different in the article

def ComputeMullerWithoutRotation(b0, b):
    eigenvalues = ComputeEigenvaluesCmatrix(b0 , b)
    t_Icp_Ic_Is = Compute_t_Icp_Ic_Is(OrganizeEigenSamp(eigenvalues , 0.008))
    t = t_Icp_Ic_Is[0]
    Icp = t_Icp_Ic_Is[1]
    Ic = t_Icp_Ic_Is[2]
    Is = t_Icp_Ic_Is[3]

    Muller_WithoutRotation = (t/2)*np.array([[1, Icp, 0, 0],
                                             [Icp, 1, 0, 0],
                                             [0, 0, Ic, Is],
                                            [0, 0, -Is, Ic]])

    return Muller_WithoutRotation


def ComputeMullerWithoutRotation_moy(b0, b):
    eigenvalues = ComputeEigenvaluesCmatrix(b0 , b)
    t_Icp_Ic_Is = Compute_t_Icp_Ic_Is_moy(OrganizeEigenSamp(eigenvalues , 0.008))
    t = t_Icp_Ic_Is[0]
    Icp = t_Icp_Ic_Is[1]
    Ic = t_Icp_Ic_Is[2]
    Is = t_Icp_Ic_Is[3]

    Muller_WithoutRotation = (t/2)*np.array([[1, Icp, 0, 0],
                                             [Icp, 1, 0, 0],
                                             [0, 0, Ic, Is],
                                            [0, 0, -Is, Ic]])

    return Muller_WithoutRotation


def f_Rotation(a): return np.array([[1, 0, 0, 0],
                                    [0, np.cos(2*a), np.sin(2*a), 0],
                                    [0, -np.sin(2*a), np.cos(2*a), 0],
                                    [0, 0, 0, 1]])

def Find_real(thetaP,thetaR,M_0, M_1, M_2, M_3, B_0, B_1, B_2, B_3):
    m=50
    thetaP_x=np.linspace(thetaP-30,thetaP+30,m)
    thetaR_y=np.linspace(thetaR-30,thetaR+30,m)
    lamda_16_lamda_15=np.zeros((m,m))
    min =40
    couple=0
    for i in range(m):
        for j in range(m):
            M_2_thetaP=f_Rotation(thetaP_x[i]*np.pi/180)@M_2@f_Rotation(-thetaP_x[i]*np.pi/180)
            M_3_thetaR=f_Rotation(thetaR_y[j]*np.pi/180)@M_3@f_Rotation(-thetaR_y[j]*np.pi/180)
            lamda_16_lamda_15[i][j]= np.log(np.sqrt(np.abs(w.Calibration_W(M_0, M_1, M_2_thetaP,M_3_thetaR, B_0, B_1, B_2, B_3)[1])))
            if(min>lamda_16_lamda_15[i][j]):
                min=lamda_16_lamda_15[i][j]
                couple=[thetaP_x[i],thetaR_y[j]]
    couple.append(min)

    return [thetaP_x,thetaR_y,lamda_16_lamda_15,couple]

