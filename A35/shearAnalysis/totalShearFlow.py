
import numpy as np

class shearFlowAndDeflection:
    """


    """
    # max value locations
    crossSectionNumber = None
    max_q = 0.
    panelNumber = 0.

    @staticmethod
    def crossSectionMaxQ(q_s01, q_s02, boomArray, openShearArray, t_sk, t_sp):
        # This function should be called for each rib location
        # From beamTheory get: T, V_y, V_z
        # From FEM get: boomArray
        # From openSectionShearFlow get: openShearArray
        # from closedSectionShearFlow get: q_s01, q_s02, d_theta_d_x

        # output: startBoom, endBoom, q_s

        i = 0

        totalShearArray = np.zeros((np.shape(openShearArray)[0], 4))
        for cell in range(1, 3):
            for panel in openShearArray:
                cellNumber, startBoom, endBoom, q_bi, panelIndex = panel
                if cell == cellNumber:
                    boom1 = boomArray[int(startBoom)]
                    # Part exclusively for spar; only part with multiple q_s0
                    tempBoomIndex = np.where(boom1 == int(endBoom))
                    if boom1[tempBoomIndex[0] + 1] == 'spar' and cell == 1:
                       # same but different direction for each cell -> take cell 1 only
                        q_s = q_bi + q_s01 - q_s02
                        totalShearArray[i] = [startBoom, endBoom, q_s, q_s / t_sp]
                    else:   # shear flow in skin
                        if cell == 1: q_s0 = q_s01
                        elif cell == 2: q_s0 = q_s02
                        q_s = q_bi + q_s0
                        totalShearArray[i] = [startBoom, endBoom, q_s, q_s / t_sk]
                    i += 1

        max_q = np.max(np.abs(totalShearArray[:, 2]))       # max shear flow in cross section
        max_q_index = np.where(max_q == np.abs(totalShearArray[:, 2]))
        abc = totalShearArray[max_q_index[0]]
        idxQMax = abc[0, 0]
        idx2QMax = abc[0, 1]

        max_shear_stress = np.max(np.abs(totalShearArray[:, 3]))       # max shear flow in cross section
        max_shear_stress_index = np.where(max_shear_stress == np.abs(totalShearArray[:, 3]))
        abc = totalShearArray[max_shear_stress_index[0]]
        idxStressMax = abc[0, 0]
        idx2StressMax = abc[0, 1]

        return max_q, idxQMax, idx2QMax, max_shear_stress, idxStressMax, idx2StressMax
        #   output format: [max_q, min_q, startBoom_max_q, endBoom_max_q, startBoom_min_q, endBoom_min_q]


    @staticmethod
    def bendingDeflectionInY(E, I_zz, theta, aero_q, x_1, x_2, x_3, x_a, A_y, A_z, B_y, B_z, C_y, R_z, k1, k2, P, x):
        cos = np.cos(theta)
        sin = np.sin(theta)
        div = - E * I_zz
        qContribution = - aero_q * x**4 * cos / 24.
        rotatedA = A_y * cos + A_z * sin
        rotatedB = B_y * cos + B_z * sin
        rotatedC = C_y * cos
        rotatedP = P * sin
        rotatedR = R_z * sin
        stepx1 = (x - x_1)**3 * rotatedA / 6.
        stepx2 = (x - x_2)**3 * rotatedB / 6.
        stepx3 = (x - x_3)**3 * rotatedC / 6.
        stepa1 = (x - x_2 + x_a / 2.)**3 * rotatedR / 6.
        stepa2 = -(x - x_2 - x_a / 2.)**3 * rotatedP / 6.

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
    def integrationConstantsY(theta, aero_q, x_1, x_2, x_a, A_y, A_z, R_z):
        cos = np.cos(theta)
        sin = np.sin(theta)

        deflectionAt2 = - ( (aero_q / 24.) * (x_2 ** 4) * sin ) \
                        - ( (x_2 - x_1) ** 3 * (- A_y * sin + A_z * cos) / 6. ) \
                        - ( (x_a / 2.) ** 3 * (R_z * cos) / 6. )

        deflectionAt1 = - aero_q * x_1 ** 4 * sin / 24.

        matrix = np.array([[x_1, 1], [x_2, 1]])
        ans = np.array([[deflectionAt1], [deflectionAt2]])

        result = np.linalg.solve(matrix, ans)
        k1 = result[0]
        k2 = result[1]

        return k1, k2

    @staticmethod
    def bendingDeflectionInZ(E, I_yy, theta, aero_q, x_1, x_2, x_3, x_a, A_y, A_z, B_y, B_z, C_y, R_z, k1, k2, P, x):
        cos = np.cos(theta)
        sin = np.sin(theta)

        deflection = ( aero_q * x**4 * sin / 24. ) + k1 * x + k2
        if x > x_1:
            deflection += (x - x_1)**3 * (- A_y * sin + A_z * cos) / 6.
        if x > x_2:
            deflection += (x - x_2)**3 * (- B_y * sin + B_z * cos) / 6.
        if x > x_3:
            deflection += (x - x_3)**3 * (- C_y * sin) / 6.
        if x > x_2 - x_a / 2.:
            deflection += (x - x_2 + x_a / 2.)**3 * R_z * cos / 6.
        if x > x_2 + x_a / 2.:
            deflection -= (x - x_2 - x_a / 2.)**3 * P * cos / 6.

        deflection /= (- E * I_yy)
        return deflection

    @staticmethod
    def maxNormalStress(boomArray, M_y, M_z, I_yy, I_zz):
        maxStress = 0.
        minStress = 0.
        for boom in range(np.shape(boomArray)[0]):
            z = boomArray[boom, 0]
            y = boomArray[boom, 1]
            stress = M_z * y / I_zz + M_y * z / I_yy
            if abs(stress) > maxStress:
                maxStress = abs(stress)
                maxStressBoom = boom
            if stress < minStress:
                minStress = stress
                minStressBoom = boom

        if maxStress == 0.:
            maxStressBoom = 0
        if minStress == 0.:
            minStressBoom = 0

        return maxStress, minStress, maxStressBoom, minStressBoom


    @staticmethod
    def rotatePoint(theta, x0, y0):
        x = x0 * np.cos(theta) - y0 * np.sin(theta)
        y = y0 * np.cos(theta) + x0 * np.sin(theta)

        return x, y
