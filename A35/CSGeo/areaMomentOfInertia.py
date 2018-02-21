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
    steinerTerm: compute stringer steinerTerms
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

        area, cy = g.stringer(w_st, h_st, t_st, t_sk)

        I_zz_b = (((t_st)**3)*w_st)/12. + (((h_st-t_st)**3)*t_st)/12.
        I_yy_b = (((t_st**3.)*(h_st-t_st))/12.) + (((w_st**3)*t_st)/12.)

        #Steiner Term
        I_zz_s = (((cy-(t_st/2.))**2)*t_st*w_st) + ((((h_st-t_st)/2.)-cy)**2)*(h_st-t_st)*t_st
        I_yy_s = 0.

        I_zz = I_zz_b + I_zz_s
        I_yy = I_yy_b + I_yy_s

        return I_zz, I_yy, area, cy

    @staticmethod
    def rotationEffect(I_zz_old, I_yy_old, angle):
        # I_zy = 0, hence ignored.
        # Sourced calcresource.com/moment-of-inertia-rotation.html
        # Verified by hand to match
        I_zz = (I_zz_old + I_yy_old) / 2. + (I_zz_old - I_yy_old) * np.cos(2.*angle)/2.
        I_yy = (I_zz_old + I_yy_old) / 2. - (I_zz_old - I_yy_old) * np.cos(2.*angle)/2.
        return I_zz, I_yy

    @staticmethod
    def steinerTerm(area, zcoord, ycoord, cy, angle, z_centroid):
        # z_centroid should be negative!
        # zcoord is positive for all instances
        y = ycoord - cy * np.cos(angle)
        z = zcoord - cy * np.sin(angle) + z_centroid        # sign may be incorrect. check

        I_zz = area * y**2
        I_yy = area * z**2

        return I_zz, I_yy

    @staticmethod
    def stringerMOI(w_st, h_st, t_st, t_sk, zcoord, ycoord, angle, z_centroid):
        # zcoord should be positive
        # z_centroid should be negative
        I_zz_temp, I_yy_temp, area, cy = areaMomentOfInertia.genericStringer(w_st, h_st, t_st, t_sk)
        I_zz_b, I_yy_b = areaMomentOfInertia.rotationEffect(I_zz_temp, I_yy_temp, angle)
        I_zz_s, I_yy_s = areaMomentOfInertia.steinerTerm(area, zcoord, ycoord, cy, angle, z_centroid)

        I_zz = I_zz_b + I_zz_s
        I_yy = I_yy_b + I_yy_s

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
    def straightMOI(C_a, h, t_sk, z_centroid):
        # Approach: for I_zz take straight line length 2xhypo, once
        # for I_zz take straight line length 2xhypo, once. centered at length/2
        # z_centroid is negative
        area, zDistToLE = g.straightSkin(C_a, h, t_sk)
        radius = h/2.
        long = C_a - radius
        doubleHypo = 2. * np.sqrt(radius**2 + long**2)
        angle = np.arctan(radius / long)

        I_zz_temp = doubleHypo * t_sk**3 / 12.
        I_yy_temp = doubleHypo**3 * t_sk / 12.

        I_zz, I_yy_b = areaMomentOfInertia.rotationEffect(I_zz_temp, I_yy_temp, angle)
        I_yy_s = area * (zDistToLE + z_centroid)**2

        I_yy = I_yy_b + I_yy_s

        return I_zz, I_yy

    @staticmethod
    def getAMOI(C_a, h, w_st, h_st, n_st, t_st, t_sk, t_sp, z_centroid):
        I_zz = 0.
        I_yy = 0.

        # Spar
        I_zz += areaMomentOfInertia.sparMOI(t_sp, h, t_sk, z_centroid)[0]
        I_yy += areaMomentOfInertia.sparMOI(t_sp, h, t_sk, z_centroid)[1]

        # Skin
        I_zz += areaMomentOfInertia.circleMOI(h, t_sk, z_centroid)[0] + areaMomentOfInertia.straightMOI(C_a, h, t_sk, z_centroid)[0]
        I_yy += areaMomentOfInertia.circleMOI(h, t_sk, z_centroid)[1] + areaMomentOfInertia.straightMOI(C_a, h, t_sk, z_centroid)[1]

        # Stringers
        stringerCoords = g.stringerPos(h, C_a, n_st, t_sk)
        for stringer in stringerCoords:
            # idx [0] not used
            z_loc = stringer[1]
            y_loc = stringer[2]
            angle = stringer[3]

            I_zz += areaMomentOfInertia.stringerMOI(w_st, h_st, t_st, t_sk, z_loc, y_loc, angle, z_centroid)[0]
            I_yy += areaMomentOfInertia.stringerMOI(w_st, h_st, t_st, t_sk, z_loc, y_loc, angle, z_centroid)[1]
        
        return I_zz, I_yy