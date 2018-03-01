import numpy as np
from beamTheory.main import beamTheory as bm
from math import *

class closedSectionShearFlow:
    """
    Inputs: V_y, C_a, h, t_sk, t_sp, z_sc, boomArray, openShearArray
    boomArray has index: [z, y, area, k, k_type, l, l_type, m, m_type, isStringerBool]
    openShearArray has index: [cellNumber, startBoom, endBoom, q_b]

    Outputs: q_s01, q_s02, d_theta_d_x
    """
    @staticmethod
    def perpDistShearCenter(z1, y1, z2, y2, z_sc):
        den = np.sqrt((y1 - y2) ** 2 + (z1 - z2) ** 2)
        if z1 == z2:
            p = abs(z1 - z_sc)
        elif y1 == y2:
            p = abs(y1)
        else:
            gradient = (y1 - y2) / (z1 - z2)
            c1 = y1 - gradient*z1
            c2 = z_sc/gradient
            z_perp = (c2-c1)/(gradient + 1/gradient)
            y_perp = -(1/gradient)*z_perp + c2
            p = np.sqrt((z_perp - z_sc)**2. + y_perp**2.)

        l = den
        return p, l


    @staticmethod
    def calculation(T, C_a, h, t_sk, t_sp, z_sc, boomArray, openShearArray, G):

        # Areas of cells
        A1 = np.pi * h**2 / 8.
        A2 = (C_a - h / 2.) * h / 2.

        # loopArray structure:
        # [ cellNumber, sum(l_i/t_i), sum_spar(l_i/t_i), sum(q_bi*l_i/t_i), sum(p*l_i*q_bi), sum_spar(p*l_i*q_bi) ]
        loopArray = np.zeros((2, 6))

        # Is global, not specific to one cell, still returned in array (only use cell 2)
        p_l_q = 0.

        # Loop cells -> loop panels -> check cell match -> get p, find t
        for cell in range(1, 3):    # cells 1, 2
            # (re)set summation values
            l_t = 0.
            l_t_spar = 0.
            q_bi = 0.           # is actually q_bi * l_i / t_i
            p_l_q_spar = 0.     # We only need the one for cell 1; hence recompute
            for panel in openShearArray:
                cellNumber, startBoom, endBoom, q_b, panelIndex = panel

                if int(cellNumber) == cell:
                    # Find boom values for associated booms
                    boom1 = boomArray[int(startBoom)]
                    boom2 = boomArray[int(endBoom)]

                    # p: perp. dist. from line of application
                    # l: length of panel
                    p, l = closedSectionShearFlow.perpDistShearCenter(boom1[0], boom2[0], boom1[1], boom2[1], z_sc)

                    # Find index of connection and get t based upon string
                    # First prevent mistakes if boom[:3] could
                    boomTemp = boom1[3:3+6:2]
                    if boom1[np.where(boomTemp == int(endBoom))[0] + 1] == 'spar':
                        if cell == 1 and p > h/2.:
                            p *= -1.
                        t = t_sp
                        l_t_spar += l / t
                        p_l_q_spar += p * l * q_b
                    else:
                        t = t_sk

                    l_t += l / t
                    q_bi += q_b * l / t
                    p_l_q += p * l * q_b

            loopArray[cell - 1] = [cell, l_t, l_t_spar, q_bi, p_l_q, p_l_q_spar]

        # 3x3 matrix for finding qs01, qs02, d(theta)/dx
        matrix = np.array([[loopArray[0, 1],    -loopArray[0, 2],   -A1*2.*G],
                           [-loopArray[1, 2],   loopArray[1, 1],    -A2*2.*G ],
                           [2. * A1,             2. * A2,             0.   ]])
        ans = np.array([[-loopArray[0, 3]],
                        [-loopArray[1, 3]],
                        [-loopArray[1, 4] + loopArray[0, 5] + T]])

        q_s01, q_s02, d_theta_d_x = np.linalg.solve(matrix, ans)
        return q_s01[0], q_s02[0], d_theta_d_x[0], loopArray


    # Idea behind torsionUpdate: run loop to find convergent value for Torsion and shear center
    # Might need to make this a for loop for a set amount of iterations to see what happens..
    @staticmethod
    def torsionUpdate(T0, V_y, C_a, h, t_sk, t_sp, G, z_sc, boomArray, openShearArray,
                      q, P, A_y, A_z, B_y, B_z, C_y, R_z, theta, x_3, x_2, x_1, x_a, x):
        if T0 == 0.:
            T0 = 0.000000001
        run = True
        while run:
            q_s01, q_s02, d_theta_d_x, loopArray = closedSectionShearFlow.calculation\
                                                            (V_y, C_a, h, t_sk, t_sp, z_sc, boomArray, openShearArray, G)

            A1 = np.pi * h / 2.
            A2 = (C_a - h / 2.) * h / 2.

            l_t1 = loopArray[0, 1]
            l_t2 = loopArray[1, 1]

            cellFactor1 = A1 ** 2 * l_t1
            cellFactor2 = A2 ** 2 * l_t2

            z_sc = - V_y * 4 * G * d_theta_d_x * (cellFactor1 + cellFactor2)

            T = bm.calculateTorqueForX(q, P, R_z, A_z, A_y, B_z, B_y, C_y, C_a, h, theta, x_3, x_2, x_1, x_a, z_sc, x)
            if abs((T/T0)*100. - 100.) > 1.:
                run = False
            else:
                T0 = T

        return T, z_sc, q_s01, q_s02, d_theta_d_x

