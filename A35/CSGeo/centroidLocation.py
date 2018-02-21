"""

Expect stringer position array [n, z, y, angle] x (n_st)
compute

"""

import numpy as np
from CSGeometry import CSGeometry as CS

class Centroid:
    """
    lol

    From CSGeometry use all functions (called with import as CS)

    nomenclature: _xyz => to/wrt axis xyz

    Computes centroid with LE->TE positive z

    """

    def computeCentroid(self, C_a, h, t_sk, t_sp, t_st, w_st, h_st, n_st):
        # stringerCoords is an array in the form [n, z, y, angle] x (n_st)
        # centroid_z_pos to be found, assume symmetry about z: centroid_y_pos = 0
        Area = 0.
        Vol = 0.

        # Spar
        tempArea, tempOff = CS.spar(h, t_sp, t_sk)
        Area += tempArea
        Vol += tempArea * tempOff

        # Circular part
        tempArea, tempOff = CS.circleSkin(h, t_sk)
        Area += tempArea
        Vol += tempArea * tempOff

        # Straight edges
        tempArea, tempOff = CS.straightSkin(C_a, h, t_sk)
        Area += 2. * tempArea
        Vol += 2. * tempArea * tempOff

        # Stringers (yey)
        tempArea, tempOff = CS.stringer(w_st, h_st, t_st, t_sk)
        stringerCoords = CS.stringerPos(h, C_a, n_st, t_sk)
        for stringer in stringerCoords:
            #stringerNumber = stringer[0]           # not used
            stringer_z = stringer[1]
            #stringer_y = stringer[2]               # not used
            stringer_angle = stringer[3]

            # centroid distance offset
            addOff = np.sin(stringer_angle) * tempOff
            dist_z = stringer_z - addOff
            # Area, 'volume'
            Area += tempArea
            Vol += tempArea * dist_z

        centroid_z_pos = Vol / Area
        return centroid_z_pos

