import numpy as np

# une méthode naïve pour calculer H_w


def H_w(M1, M2, B1, B2):
    invM1 = np.linalg.inv(M1)
    invB2 = np.linalg.inv(B2)
    B1_B2 = invB2@B2  # (B1)^-1*B2
    M1_M2 = invM1@M2  # (M1)^-1*M2
    zero = np.zeros((4, 4))
    for i in range(4):
        for j in range(4):
            zero[i, j] = 1
            if(i == 0 and j == 0):
                W = (M1_M2@zero-zero@M1_M2).reshape((16, 1))
            else:
                W = np.concatenate(
                    (W, (M1_M2@zero-zero@B1_B2).reshape((16, 1))), axis=1)
    return W


# I = np.identity(4)
# I[0][3] = 4
# I[3][0] = 4
# wI = np.ones((4, 4))
# print(H_w(I, wI, wI, I))
