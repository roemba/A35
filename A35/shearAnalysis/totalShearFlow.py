
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
