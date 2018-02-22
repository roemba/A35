import numpy as np

"""

Area Moment Of Inertia (output in m^4)
call

"""

class areaMomentOfInertia:

    """
    Inputs z_centroid, y_centroid, z_boom, y_boom, boom_area
    """
    @staticmethod
    def getMOI(z_c, y_c, z, y, b_a):

        #distance from centroid
        d =

        I_zz = b_a*((y-y_c)**2.)
        I_yy = b_a*((z-z_c)**2.)

        return I_zz, I_yy