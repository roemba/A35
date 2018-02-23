import numpy as np

"""

Area Moment Of Inertia (output in m^4)
Essentially Steiner Term

"""

class numbericalMOI:




    """


    Inputs: z_centroid, y_centroid, [z_boom, y_boom, boom_area, k, l, m]
    Outputs: I_zz, I_yy, I_zy
    """
    @staticmethod
    def getMOI(z_c, y_c, boomArray):

        arrayI_zz = boomArray[2]*((boomArray[1]-y_c)**2.)
        arrayI_yy = boomArray[2]*((boomArray[0]-z_c)**2.)
        arrayI_zy = boomArray[2]*((boomArray[1]-y_c)*(boomArray[0]-z_c))
        #print arrayI_zz
        #print arrayI_yy

        I_zz = np.sum(arrayI_zz)
        I_yy = np.sum(arrayI_yy)

        #I_zy, be careful with this one
        I_zy = np.sum(arrayI_zy)
        
        # Due to symmetry
        #I_zy = 0.

        return I_zz, I_yy, I_zy
