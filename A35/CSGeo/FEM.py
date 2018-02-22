
import numpy as np

class FEM:
    """
    Borrowed stringerPos function from CSGeometry to have this standalone for call and output

    Output array output has z and y coordinates. pitch and boom area are
    assumed constants; in loop the boom area may be modified to match
    z and y are positive according to reference frame (z in dir of flight, y upwards)
    Reference point for axis system is at Leading Edge.
    """

    """
    Produce output in form: [z_pos, y_pos, [k, l (, m)]]
    """

    @staticmethod
    def crossSection(h, C_a, n_st):
        radius = h /2.
        long = C_a - radius
        semi = radius * np.pi
        hypo = np.sqrt(radius**2 + long**2)
        angle = np.arctan(radius / long)
        length = semi + 2. * hypo
        pitch = length / n_st
        return radius, long, semi, hypo, angle, length, pitch

    @staticmethod
    def discretization(n_sector_1, n_sector_2, n_sector_4, C_a, h, n_st):
        """
        Sector 1: TE to negative y connection spar
        Sector 2: Semi-circular LE
        Sector 3: positive y connection spar to TE: equal spacing to Sector 1
        Sector 4: Spar
        Fixed point 1: TE
        Fixed point 2: negative y spar-skin connection
        Fixed point 3: positive y spar-skin connection
        """

        # get stringer pitch values
        stringerPitch = FEM.crossSection(h, C_a, n_st)[-1]

        # Now check each sector value



        return outputArray

