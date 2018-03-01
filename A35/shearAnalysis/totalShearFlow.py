
import numpy as np

class shearFlowAndDeflection:
    """


    """
    # max value locations
    crossSectionNumber = None
    max_q = 0.
    panelNumber = 0.

    @staticmethod
    def crossSectionMaxQ(q_s01, q_s02, boomArray, openShearArray):
        # This function should be called for each rib location
        # From beamTheory get: T, V_y, V_z
        # From FEM get: boomArray
        # From openSectionShearFlow get: openShearArray
        # from closedSectionShearFlow get: q_s01, q_s02, d_theta_d_x

        # output: startBoom, endBoom, q_s

        i = 0

        totalShearArray = np.zeros((np.shape(openShearArray)[0], 3))
        for cell in range(1, 3):
            for panel in openShearArray:
                cellNumber, startBoom, endBoom, q_bi, panelIndex = panel
                if cell == cellNumber:
                    boom1 = boomArray[int(startBoom)]
                    boomTemp = boom1[3:3 + 6:2]
                    # Part exclusively for spar; only part with multiple q_s0
                    tempBoomIndex = np.where(boom1 == int(endBoom))
                    if boom1[tempBoomIndex[0] + 1] == 'spar' and cell == 1:
                       # same but different direction for each cell -> take cell 1 only
                        q_s = q_bi + q_s01 - q_s02
                        totalShearArray[i] = [startBoom, endBoom, q_s]
                    else:   # shear flow in skin
                        if cell == 1: q_s0 = q_s01
                        elif cell == 2: q_s0 = q_s02
                        q_s = q_bi + q_s0
                        totalShearArray[i] = [startBoom, endBoom, q_s]
                    i += 1

        max_q = np.max(np.abs(totalShearArray[:, 2]))       # max shear flow in cross section
        max_q_index = np.where(max_q == np.abs(totalShearArray[:, 2]))
        abc = totalShearArray[max_q_index[0]]
        idxMax = abc[0, 0]
        idx2Max = abc[0, 1]

        return max_q, idxMax, idx2Max
        #   output format: [max_q, min_q, startBoom_max_q, endBoom_max_q, startBoom_min_q, endBoom_min_q]


    @staticmethod
    def bendingDeflectionInY(E, I_zz, theta, aero_q, x_1, x_2, x_3, x_a, A_y, A_z, B_y, B_z, C_y, k1, k2, R_z, P, x):
        cos = np.cos(theta)
        sin = np.sin(theta)
        div = - E * I_zz
        qContribution = - aero_q * x**4 * cos
        rotatedA = A_y * cos + A_z * sin
        rotatedB = B_y * cos + B_z * sin
        rotatedC = C_y * cos
        rotatedP = P * sin
        rotatedR = R_z * sin
        stepx1 = (x - x_1)**3 * rotatedA / 6.
        stepx2 = (x - x_2)**3 * rotatedB / 6.
        stepx3 = (x - x_3)**3 * rotatedC / 6.
        stepa1 = (x - x_2 + x_a / 2.) * rotatedR / 6.
        stepa2 = -(x - x_2 - x_a / 2.) * rotatedP / 6.

        deflection = k1 * x + k2 + qContribution
        if x > x_1:
            deflection += stepx1
        if x > x_2:
            deflection += stepx2
        if x > x_3:
            deflection += stepx3
        if x > x_2 - x_a / 2.:
            deflection += stepa1
        if x > x_2 + x_a / 2.:
            deflection += stepa2

        deflection /= div
        return deflection

    @staticmethod
    def bendingDeflectionInZ(E, I_yy, theta, aero_q, x_1, x_2, x_3, x_a, A_y, A_z, B_y, B_z, C_y, R_z, P, x):


