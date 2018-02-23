import numpy as np

"""

Area Moment Of Inertia (output in m^4)
Essentially Steiner Term

"""

class numericalMOI:

    """


    Inputs: z_centroid, y_centroid, [z_boom, y_boom, boom_area, [k, pnltype], [l, pnltype], [m, pnltype]]
    Outputs: I_zz, I_yy, I_zy
    """
    @staticmethod
    def getMOI(z_c, y_c, boomArray):

        arrayI_zz = []
        arrayI_yy = []
        arrayI_zy = []

        for i in range(0, len(boomArray)):
            arrayI_zz.append(boomArray[i][2]*((boomArray[i][1]-y_c)**2.))
            arrayI_yy.append(boomArray[i][2]*((boomArray[i][0]-z_c)**2.))
            arrayI_zy.append(boomArray[i][2]*((boomArray[i][1]-y_c)*(boomArray[i][0]-z_c)))

        I_zz = np.sum(arrayI_zz)
        I_yy = np.sum(arrayI_yy)

        #I_zy, be careful with this one
        I_zy = np.sum(arrayI_zy)

        # Due to symmetry
        #I_zy = 0.

        return I_zz, I_yy, I_zy

