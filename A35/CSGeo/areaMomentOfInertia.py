import numpy as np
import CSGeometry.CSGeometry as g

"""

Area Moment Of Inertia (output in m^4)
call

"""

class areaMomentOfInertia:

    """
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
    straightPart: compute AMOI and steinerTerms

    """


    def genericStringer(self, w_st, h_st, t_st):

        a, cy = g.stringer(w_st, h_st, t_st)

        I_zz_b = (((t_st)**3)*w_st)/12. + (((h_st-t_st)**3)*t_st)/12.
        I_yy_b = (((t_st**3.)*(h_st-t_st))/12.) + (((w_st**3)*t_st)/12.)

        #Steiner Term
        I_zz_s = (((cy-(t_st/2.))**2)*t_st*w_st) + ((((h_st-t_st)/2.)-cy)**2)*(h_st-t_st)*t_st
        I_yy_s = 0.

        I_zz = I_zz_b + I_zz_s
        I_yy = I_yy_b + I_yy_s

        return I_zz, I_yy, cy

    def rotationEffect(self, I_zz, I_yy, angle):


    def getAMOI(self):




        return I_zz, I_yy