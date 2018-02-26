import numpy as np
from math import *

class closedSectionShearFlow:

    """
    Inputs: I_zz, I_yy, I_zy, V_z, V_y, L_c, L_t, h, G, t_sk, t_sp,  k, k_type, l, l_type, m, m_type

    Outputs: q_b
    """
    @staticmethod
    def calculation(I_zz, I_yy, I_zy, V_z, V_y, L_c, L_t, h, G, t_sk, t_sp, boomArray):

        #CHECK IF USING t_sp
        rateOfTwist = (1./(np.pi*((h/2.)**2)))*(((q_s01*L_c)/(t_sk*G))+(((q_s01 - q_s02)*L_t)/(t_sp*G)))
        rateOfTwist =

        return 0.
    """
    
    """
    @staticmethod
    def geometry():

        return 0
