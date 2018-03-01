
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

        totalShearArray = np.zeros((np.shape(openShearArray)[0], 2))
        for cell in range(1, 3):
            for panel in range(np.shape(openShearArray)[0]):
                cellNumber, startBoom, endBoom, q_bi, panelIndex = panel
                if cell == cellNumber:
                    boom1 = boomArray[int(startBoom)]
                    boomTemp = boom1[3:3 + 6:2]
                    # Part exclusively for spar; only part with multiple q_s0
                    if boom1[np.where(boomTemp == int(endBoom))[0] + 1] == 'spar':
                        if cell == 1:   # same but different direction for each cell -> take cell 1 only
                            q_s = q_bi + q_s01 - q_s02
                            totalShearArray[panelIndex] = [startBoom, endBoom, q_s]
                    else:   # shear flow in skin
                        if cell == 1: q_s0 = q_s01
                        elif cell == 2: q_s0 = q_s02
                        q_s = q_bi + q_s0
                        totalShearArray[panelIndex] = [startBoom, endBoom, q_s]

        max_q = np.max(totalShearArray[:, 2])       # max shear flow in cross section
        idxMax = totalShearArray[np.where(max_q == totalShearArray[:, 2])[0]][0]
        min_q = np.min(totalShearArray[:, 2])       # min shear flow in cross section
        idxMin = totalShearArray[np.where(min_q == totalShearArray[:, 2])[0]][1]

        return max_q, min_q, idxMax, idxMax + 1, idxMin, idxMin + 1
        #   output format: [max_q, min_q, startBoom_max_q, endBoom_max_q, startBoom_min_q, endBoom_min_q]
