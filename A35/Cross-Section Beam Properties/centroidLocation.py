"""

Expect stringer position array [n, z, y, angle] x (n_st)
compute

"""

import numpy as np

class Centroid:
    """
    lol

    From CSGeometry.py get

    nomenclature: _xyz => to/wrt axis xyz

    Computes centroid with LE->TE positive z

    """

    def stringer(self, w_st, h_st, t_st):
        stArea = (w_st + h_st) * t_st
        offsetFromSkin = (h_st * h_st + w_st * t_st) * t_st / (2 * stArea)
        return stArea, offsetFromSkin

    def straightSkin(self, C_a, h, t_sk):
        radius = h/2
        long = C_a - radius
        hypo = np.sqrt(long^2 + radius^2)
        straightArea = hypo * t_sk
        #dist_z = radius/2           # distance to z axis from centroid ----- due to symmetry irrelevant
        dist_y = long/2 + radius    # distance to y axis from centroid
        return straightArea, dist_y

    def circleSkin(self, h, t_sk):
        radius = h/2
        circArea = np.pi * radius * t_sk
        dist_y = radius - (2 * radius / np.pi)   # check this?
        return circArea, dist_y

    def spar(self, h, t_sp):
        spArea = h * t_sp
        dist_y = h/2
        return spArea, dist_y


    # combine all..

    def computeCentroid(self, C_a, h, t_sk, t_sp, t_st, w_st, h_st, stringerCoords):
        # stringerCoords is an array in the form [n, z, y, angle] x (n_st)
        # centroid_z_pos to be found, assume symmetry about z: centroid_y_pos = 0
        Area = 0.
        Vol = 0.

        # Spar
        tempArea, tempOff = Centroid.spar(h, t_sp)
        Area += tempArea
        Vol += tempArea * tempOff

        # Circular part
        tempArea, tempOff = Centroid.circleSkin(h, t_sk)
        Area += tempArea
        Vol += tempArea * tempOff

        # Straight edges
        tempArea, tempOff = Centroid.straightSkin(C_a, h, t_sk)
        Area += 2 * tempArea
        Vol += 2 * tempArea * tempOff

        # Stringers (yey)
        tempArea, tempOff = Centroid.stringer(w_st, h_st, t_st)
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

