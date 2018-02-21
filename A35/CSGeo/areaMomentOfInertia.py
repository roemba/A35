import numpy as np
from CSGeometry import CSGeometry as g

"""

Area Moment Of Inertia (output in m^4)
call

"""

class areaMomentOfInertia:

    """
    Use subscript b for base, subscript s for steiner

    Order:
    Stringer:
    genericStringer: compute AMOI in terms of geometry of stringer only
    rotationEffect: compute rotation effect on stringer (apply to genericStringer)
    steinerTerm: compute steinerTerms
    totalStringer: combine rotationEffects and steinerTerm

    Spar:
    sparMOI: compute AMOI and steinerTerms

    Skin:
    circleMOI: compute AMOI and steinerTerms
    straightMOI: compute AMOI and steinerTerms

    Total:
    getAMOI: Add all AMOI
    """

    @staticmethod
    def genericStringer(w_st, h_st, t_st, t_sk):

        cy = g.stringer(w_st, h_st, t_st, t_sk)[1]

        I_zz_b = (((t_st)**3)*w_st)/12. + (((h_st-t_st)**3)*t_st)/12.
        I_yy_b = (((t_st**3.)*(h_st-t_st))/12.) + (((w_st**3)*t_st)/12.)

        #Steiner Term
        I_zz_s = (((cy-(t_st/2.))**2)*t_st*w_st) + ((((h_st-t_st)/2.)-cy)**2)*(h_st-t_st)*t_st
        I_yy_s = 0.

        I_zz = I_zz_b + I_zz_s
        I_yy = I_yy_b + I_yy_s

        return I_zz, I_yy

    @staticmethod
    def rotationEffect(I_zz, I_yy, angle):
        return I_zz, I_yy, angle

    @staticmethod
    def steinerTerm(area, zpos, ypos, cy, angle):
        return I_zz, I_yy

    @staticmethod
    def sparMOI(t_sp, h, t_sk, z_centroid):
        area, offset = g.spar(h, t_sp, t_sk)
        height = h - 2. * t_sk
        # Spar in own centroid at [t_sp, 0]
        I_zz_b = area * height**2
        I_yy_b = area * t_sp**2

        # Steiner term towards global centroid
        I_zz_s = 0.             # symmetric airfoil -> both centroids coincide with chord (z-axis)
        I_yy_s = area * (offset + z_centroid)**2        # z_centroid should be negative. check?


        # Totals
        I_zz = I_zz_b + I_zz_s
        I_yy = I_yy_b + I_yy_s

        return I_zz, I_yy

    @staticmethod
    def circleMOI(h, t_sk, z_centroid):
        # Approach: Take MOI of thin-walled circle and divide by 2,
        # for I_yy: compute steiner term from center of circle to centroid. (math behind it verified)
        area = g.circleSkin(h, t_sk)[0]
        I_zz = np.pi * t_sk * h**3 / 8. / 2.    # is final, centered about z-axis

        I_yy_b = I_zz
        I_yy_s = area * (h/2. + z_centroid)**2

        I_yy = I_yy_b + I_yy_s

        return I_zz, I_yy

    @staticmethod
    def getAMOI():
        return I_zz, I_yy