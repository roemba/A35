import numpy as np

class openSectionShearFlow:
    """
    This is only per cell, there are two cells in total.

    Inputs: I_zz, I_yy, I_zy, V_z, V_y, boomArray
    boomArray has index: [z, y, area, k, k_type, l, l_type, m, m_type, isStringerBool]
    Outputs: q_b
    """

    @staticmethod
    def perBoomCalc(I_zz, I_yy, I_zy, V_z, V_y, z, y, area):
        '''
        This part function computes the z and y component shear flows for a specific boom
        inputs z, y, area are thus boom-specific
        '''
        zArea = z * area
        yArea = y * area

        div = I_zz * I_yy - I_zy**2

        zIzz = I_zz * V_z
        zIzy = - I_zy * V_y
        yIyy = I_yy * V_y
        yIzy = - I_zy * V_z

        q1 = zArea * (zIzz + zIzy) / div
        q2 = yArea * (yIyy + yIzy) / div

        return - q1 - q2

    @staticmethod
    def totOpenCalc(I_zz, I_yy, I_zy, V_z, V_y, boomArray, n_sector_1, n_sector_2, n_sector_4):
        '''
         This part function uses perBoomCalc to find the shear flows in each section
         taking the first link in each cell (by convention) to have zero shear.
         outArray has form: [cellNumber, startBoom, endBoom, q_b]
         Output is dependent on the currently evaluated cell.
        '''
        # Cell 2 (Trailing edge and spar)
        for panelIndex in range(2 * (n_sector_1 + 1) + n_sector_4 + 1):
            cellNumber = 2
            # Sector 1 (skin)
            if panelIndex == 0:
                q_b = 0.
                startBoom = 0
                endBoom = 1
                outArray = np.array([cellNumber, startBoom, endBoom, q_b])
                continue
            elif panelIndex <= n_sector_1:
                startBoom = panelIndex
                endBoom = panelIndex + 1

            # Sector 4 (spar)
            elif panelIndex <= n_sector_1 + n_sector_4 + 1:
                if panelIndex == n_sector_1 + 1:
                    startBoom = panelIndex
                    endBoom = n_sector_1 * 2 + n_sector_2 + 2
                elif panelIndex == n_sector_1 + 1 + n_sector_4:
                    startBoom = n_sector_1 + n_sector_2 + 1 + panelIndex
                    endBoom = n_sector_1 + n_sector_1 + 2
                else:
                    startBoom = n_sector_1 + n_sector_2 + 1 + panelIndex
                    endBoom = startBoom + 1

            # Sector 3 (skin)
            else:
                startBoom = panelIndex + n_sector_2 - 2
                endBoom = startBoom + 1
                if panelIndex == 2 * (n_sector_1 + 1) + n_sector_4:
                    endBoom = 0    # connected to boom 0 (TE)

            z, y, area = boomArray[startBoom, :3]
            # negative because clockwise instead of counter-clockwise
            q_b = - openSectionShearFlow.perBoomCalc(I_zz, I_yy, I_zy, V_z, V_y, z, y, area)
            if panelIndex == 1:
                q_b = outArray[-1]
            else:
                q_b += outArray[-1, -1] # Because it's supposed to be a cumulative thing

            np.vstack([outArray, [cellNumber, startBoom, endBoom, q_b]])

        # Cell 1 (Leading edge and spar)
        for panelIndex in range(n_sector_2 + n_sector_4 + 2):
            cellNumber = 1
            panelOffset = n_sector_1 + 1
            # Sector 2 (skin)
            if panelIndex == 0:
                q_b = 0.
                startBoom = panelOffset
                endBoom = 1 + panelOffset
                # Setting to 0 again since different cell, building upon this one now
                np.vstack([outArray, [cellNumber, startBoom, endBoom, q_b]])
                continue
            elif panelIndex <= n_sector_2:
                startBoom = panelIndex + panelOffset
                endBoom = startBoom + 1

            # Sector 4 (spar)
            else:
                if panelIndex == n_sector_2 + 1:
                    startBoom = panelIndex + panelOffset
                    endBoom = 2 * panelOffset + n_sector_2 + n_sector_4

                else:
                    startBoom = 2 * panelOffset + 2 * n_sector_2 + n_sector_4 + 1 - panelIndex
                    endBoom = startBoom - 1
                    if panelIndex == n_sector_2 + n_sector_4 + 1:
                        endBoom = panelOffset

            z, y, area = boomArray[startBoom, :3]
            q_b = openSectionShearFlow.perBoomCalc(I_zz, I_yy, I_zy, V_z, V_y, z, y, area)
            q_b += outArray[-1, -1]

            np.vstack([outArray, [cellNumber, startBoom, endBoom, q_b]])

        return outArray
        # outArray has the form [cellNumber, startBoom, endBoom, q_b], q_b is adjusted for counterclockwise positive.
