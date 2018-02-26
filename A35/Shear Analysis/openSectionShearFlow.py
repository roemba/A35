import numpy as np
from math import *


class openSectionShearFlow:

    """
    This is only per cell, there are two cells in total.

    Inputs: I_zz, I_yy, I_zy, [[z_boom, y_boom, boom_area, k, pnltype, l, pnltype, m, pnltype]]
    Outputs: q_b
    """
    @staticmethod
    def calculation(I_zz, I_yy, I_zy, V_z, V_y, boomArray):

        arrayAz = []
        arrayAy = []

        for i in range(0, len(boomArray)):
            arrayAz.append(boomArray[i][2]*boomArray[i][0])
            arrayAy.append(boomArray[i][2]*boomArray[i][1])

        print arrayAz
        print arrayAy

        sumAz = np.sum(arrayAz)
        sumAy = np.sum(arrayAy)

        p1 = (((I_zz*V_z) - (I_zy*V_y))/(I_zz*I_yy - (I_zy**2)))*sumAz
        p2 = (((I_yy*V_y) - (I_zy*V_z))/(I_zz*I_yy - (I_zy**2)))*sumAy

        q_b = -1.*p1 - p2

        return q_b
