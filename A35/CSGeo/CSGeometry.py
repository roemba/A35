
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
import matplotlib.pyplot as plt

class CSGeometry:
    """

    Take h, C_a, n_st to compute stringer coordinates

    csLength computes total length of skin cross-section

    """
    @staticmethod
    def stringerPos(h, C_a, n_st, t_sk):
        radius = h/2. - t_sk
        triangleLongSection = C_a - radius
        semiCircle = (radius) * np.pi
        hypo = np.sqrt(radius**2 + triangleLongSection**2)      # Compute length straight part of aileron
        angle = np.arctan(radius / triangleLongSection)         # Compute half-angle
        length = semiCircle + 2.*hypo                           # Total length cross-section

        centerSpacing = length / n_st                           # stringer pitch
        print 'stringer spacing: ', centerSpacing
        finalArray = np.zeros((int(n_st), 4))
        for i in range(int(n_st)):
            distFromTE = centerSpacing * (float(i) + 0.5)

            if distFromTE > hypo:
                distFromTE -= hypo

                if distFromTE > semiCircle:
                    distFromTE -= semiCircle
                    y_loc = radius - distFromTE * np.sin(angle)
                    z_loc = radius + distFromTE * np.cos(angle)
                    outAngle = angle

                else:   # distFromTE < semiCircle
                    tempAngle = distFromTE / radius
                    z_loc = radius - radius * np.sin(tempAngle)
                    y_loc = - radius * np.cos(tempAngle)


                    outAngle = -1. * tempAngle           # negative for simplified centroid & AMOI calc

            else:       # distFromTE < hypo
                y_loc = 0. - distFromTE * np.sin(angle)
                z_loc = C_a - distFromTE * np.cos(angle)
                outAngle = angle

            finalArray[i] = [i, z_loc, y_loc, outAngle]

        #plt.plot(finalArray[:,1], finalArray[:,2])
        #plt.show()
        return finalArray

    @staticmethod
    def stringer(w_st, h_st, t_st, t_sk):
        stArea = (w_st + h_st) * t_st
        offsetFromSkin = ((h_st * h_st + w_st * t_st) * t_st / (2. * stArea)) + t_sk
        return stArea, offsetFromSkin

    @staticmethod
    def straightSkin(C_a, h, t_sk):
        radius = h/2.
        long = C_a - radius
        hypo = np.sqrt(long**2 + radius**2)
        straightArea = hypo * t_sk
        #dist_z = radius/2           # distance to z axis from centroid ----- due to symmetry irrelevant
        dist_y = long/2. + radius    # distance to y axis from centroid
        return straightArea, dist_y

    @staticmethod
    def circleSkin(h, t_sk):
        radius = h/2.
        circArea = np.pi * radius * t_sk
        dist_y = radius - (2. * radius / np.pi)   # check this?
        return circArea, dist_y

    @staticmethod
    def spar(h, t_sp, t_sk):
        spArea = (h - (2. * t_sk)) * t_sp
        dist_y = h/2.
        return spArea, dist_y

