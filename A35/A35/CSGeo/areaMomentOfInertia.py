import numpy as np
import CSGeometry.CSGeometry as g

"""

Area Moment Of Inertia (output in m^4)
call

"""

class areaMomentOfInertia:

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

    def rotationEffects(self, I_zz, I_yy, angle):

        I_zz_r = ((I_zz + I_yy)/2) + ((I_zz - I_yy)/2)*np.cos(2*angle)
        I_yy_r = ((I_zz + I_yy)/2) - ((I_zz - I_yy)/2)*np.sin(2*angle)

        return I_zz_r, I_yy_r

    def getAMOI(self):




        return I_zz, I_yy