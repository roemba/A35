
"""
Stringer locations in [z,y, angle] locations
(angle is rotation wrt z-axis for centroid, MOI)
With [0, 0] at LE of aileron

inputs:
h
C_a
n_st
"""

import numpy as np

class CSGeometry:
    """

    Take h, C_a, n_st to compute stringer coordinates

    csLength computes total length of skin cross-section

    """

    def stringerPos(self, h, C_a, n_st):
        radius = h/2.
        triangleLongSection = C_a - radius
        semiCircle = radius * np.pi
        hypo = np.sqrt(radius**2 + triangleLongSection**2)        # Compute length straight part of aileron
        angle = np.arctan(radius / triangleLongSection)         # Compute half-angle
        length = semiCircle + 2.*hypo                           # Total length cross-section

        centerSpacing = length / n_st                           # stringer pitch
        finalArray = np.zeros((n_st, 4))
        for i in n_st:
            distFromTE = centerSpacing * (i + 0.5)

            if distFromTE > hypo:
                distFromTE -= hypo

                if distFromTE > semiCircle:
                    distFromTE -= semiCircle
                    y_loc = radius - distFromTE * np.sin(angle)
                    z_loc = radius + distFromTE * np.cos(angle)
                    outAngle = angle

                else:   # distFromTE < semiCircle
                    tempAngle = distFromTE / radius
                    y_loc = 0. - np.cos(tempAngle)
                    z_loc = radius - np.sin(tempAngle)
                    outAngle = -1 * tempAngle           # negative for simplified centroid & AMOI calc

            else:       # distFromTE < hypo
                y_loc = 0. - distFromTE * np.sin(angle)
                z_loc = C_a - distFromTE * np.cos(angle)
                outAngle = angle

            finalArray[i] = [i, z_loc, y_loc, outAngle]
        return finalArray

    def stringer(self, w_st, h_st, t_st):
        stArea = (w_st + h_st) * t_st
        offsetFromSkin = (h_st * h_st + w_st * t_st) * t_st / (2. * stArea)
        return stArea, offsetFromSkin

    def straightSkin(self, C_a, h, t_sk):
        radius = h/2.
        long = C_a - radius
        hypo = np.sqrt(long**2 + radius**2)
        straightArea = hypo * t_sk
        #dist_z = radius/2           # distance to z axis from centroid ----- due to symmetry irrelevant
        dist_y = long/2. + radius    # distance to y axis from centroid
        return straightArea, dist_y

    def circleSkin(self, h, t_sk):
        radius = h/2.
        circArea = np.pi * radius * t_sk
        dist_y = radius - (2. * radius / np.pi)   # check this?
        return circArea, dist_y

    def spar(self, h, t_sp):
        spArea = h * t_sp
        dist_y = h/2.
        return spArea, dist_y

