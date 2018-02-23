import numpy as np
from math import *


class openSectionShearFlow:

    """
    This is only per cell, there are two cells in total.

    Inputs: I_zz, I_yy, I_zy, [[z_boom, y_boom], boom_area, [k, l, m]]
    Outputs: q_b
    """
    @staticmethod
    def calculation(I_zz, I_yy, I_zy, V_z, V_y, boomArray):

        p1 = (((I_zz*V_z) - (I_zy*V_y))/(I_zz*I_yy - (I_zy**2)))*np.sum((boomArray[1]*boomArray[0][0]))
        p2 = (((I_yy*V_y) - (I_zy*V_z))/(I_zz*I_yy - (I_zy**2)))*np.sum((boomArray[1]*boomArray[0][1]))

        q_b = -1.*p1 - p2

        return q_b
