import numpy as np
from math import *

class closedSectionShearFlow:
    """
    Inputs: I_zz, I_yy, I_zy, V_z, V_y, L_c, L_t, h, G, t_sk, t_sp, boomArray, openShearArray
    boomArray has index: [z, y, area, k, k_type, l, l_type, m, m_type, isStringerBool]
    openShearArray has index: [cellNumber, startBoom, endBoom, q_b]

    Outputs: q_b, z_sc
    """
    @staticmethod
    def calculation(I_zz, I_yy, I_zy, V_z, V_y, h, t_sk, t_sp, boomArray):



        return q_b, z_sc
