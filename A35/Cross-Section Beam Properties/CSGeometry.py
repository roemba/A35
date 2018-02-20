
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

class stringerPositions:
    """

    Take h, C_a, n_st to compute stringer coordinates

    csLength computes total length of skin cross-section

    """

    def Pos(self, h, C_a, n_st):
        radius = h/2
        triangleLongSection = C_a - radius
        semiCircle = radius * np.pi
        hypo = np.sqrt(radius^2 + triangleLongSection^2)        # Compute length straight part of aileron
        angle = np.atan(radius / triangleLongSection)           # Compute half-angle
        length = semiCircle + 2*hypo                            # Total length cross-section

        centerSpacing = length / n_st                           # stringer pitch
        finalArray = np.zeros((n_st, 4))
        for i in n_st:
            distFromTE = centerSpacing * i

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
                    outAngle = -1 * tempAngle           # negative for simplified centroid calc

            else:       # distFromTE < hypo
                y_loc = 0. - distFromTE * np.sin(angle)
                z_loc = C_a - distFromTE * np.cos(angle)
                outAngle = angle

            finalArray[i] = [i, z_loc, y_loc, outAngle]



