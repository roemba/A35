import numpy as np

class beamTheory:
    def displacementStepFunctions(self, E, I, d_3, d_1, x_3, x_2, x_1):
        b = x_3 - x_1
        a = x_2 - x_1

        coefficientMatrix = np.matrix([[-1., ((1./6.)*(b-a)**3. - (b*(a**2.))/6. + (1./6.)*a**3.)/((1./6.)*b**3. - (b*(a**2.))/6.), 0.],
                                    [1., -1., 1.],
                                    [0., -a, b]])

        resultMatrix = np.matrix([[(E * I * d_3 - E * I * d_1 * ((b / a) - 1.)) / ((1. / 6.) * b ** 3. - (b * (a ** 2.)) / 6.)],
                                  [0.],
                                  [0.]])
        solutionMatrix = np.linalg.solve(coefficientMatrix, resultMatrix)
        solutionArray = np.array(solutionMatrix)

        print "A_y = {a_y}, C_y = {c_y}, B_y = {B_y}".format(
            a_y=round(solutionArray[0, :][0], 2), c_y=round(solutionArray[2, :][0], 2),
            B_y=round(solutionArray[1, :][0], 2))

        return solutionMatrix

    def staticEquilibrium(self, c_a, l_a, x_1, x_2, x_3, x_a, P, q, a_y, c_y, h):
        coefficientMatrix = np.matrix([[1., 0., 0., 0.],
                                       [0., 1., 1., 1.],
                                       [0., 0., x_2-x_1, x_a/2.],
                                       [0.25*c_a, 0., 0., -h/2.]
                                    ])
        #[0] = B_y, [1] = B_z, [2] = A_z, [3] = R_z

        resultMatrix = np.matrix([[q*l_a - a_y - c_y],
                                  [P],
                                  [-P*(x_a/2)],
                                  [-0.25*c_a*(a_y+c_y)-P*(h/2.)]
                                  ])

        solutionMatrix = np.linalg.solve(coefficientMatrix, resultMatrix)

        solutionArray = np.array(solutionMatrix)
        print "B_y = {b_y}, B_z = {b_z}, A_z = {a_z}, R_z = {r_z}".format(
            b_y=round(solutionArray[0, :][0], 2), b_z=round(solutionArray[1, :][0], 2),
            a_z=round(solutionArray[2, :][0], 2), r_z=round(solutionArray[3, :][0], 2))

        return solutionMatrix
