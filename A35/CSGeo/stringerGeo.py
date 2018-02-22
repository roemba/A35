
import numpy as np

class stringerGeo:
    """
    Borrowed stringerPos function from CSGeometry to have this standalone for call and output

    Tested and true. Array output has [pitch, z, y, boom area] for z and y coordinates. pitch and boom area are
    assumed constants; in loop the boom area may be modified to match
    z and y are positive according to reference frame (z in dir of flight, y upwards)
    Reference point for axis system is at Leading Edge.
    """

    @staticmethod
    def boomArea(w_st, h_st, t_st): return (w_st + h_st) * t_st

    @staticmethod
    def stringerPos(h, C_a, n_st, w_st, h_st, t_st):
        radius = h / 2.
        triangleLongSection = C_a - radius
        semiCircle = (radius) * np.pi
        hypo = np.sqrt(radius ** 2 + triangleLongSection ** 2)  # Compute length straight part of aileron
        angle = np.arctan(radius / triangleLongSection)         # Compute half-angle
        length = semiCircle + 2. * hypo                         # Total length cross-section

        centerSpacing = length / n_st                           # stringer pitch

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

                else:                                               # distFromTE < semiCircle
                    tempAngle = distFromTE / radius
                    z_loc = radius - radius * np.sin(tempAngle)
                    y_loc = - radius * np.cos(tempAngle)

                    outAngle = -1. * tempAngle                      # negative for simplified centroid & AMOI calc

            else:                                                   # distFromTE < hypo
                y_loc = 0. - distFromTE * np.sin(angle)
                z_loc = C_a - distFromTE * np.cos(angle)
                outAngle = angle

            finalArray[i] = [centerSpacing, - z_loc, y_loc, stringerGeo.boomArea(w_st, h_st, t_st)]    #, outAngle]

        return finalArray
