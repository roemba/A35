import numpy as np

class beamTheory:
    def staticEquilibrium(self, E, I, c_a, l_a, x_1, x_2, x_3, x_a, P, q, h, d_3, d_1):
        coefficientMatrix = np.matrix([[1., 0., 1., 0., 1., 0., 0., 0.],
                                       [0., 1., 0., 1., 0., 1., 0., 0.],
                                       [0., x_2-x_1, 0., 0., 0., x_a/2., 0., 0.],
                                       [-(x_2-x_1), 0., 0., 0., x_3-x_2, 0., 0., 0.],
                                       [0., 0., 0., 0., 0., h/2., 0., 0.],
                                       [0., 0., 0., 0., 0., 0., x_1, 1.],
                                       [(1./6.)*(x_2-x_1)**3., 0., 0., 0., 0., 0., x_2, 1.],
                                       [(1./6.)*(x_3-x_1)**3., 0., (1./6.)*(x_3-x_2)**3., 0., 0., 0., x_3, 1.]

                                    ])

        resultMatrix = np.matrix([[q*l_a],
                                  [P],
                                  [-P*(x_a/2.)],
                                  [0.],
                                  [P*(h/2.)-q*l_a*0.25*c_a],
                                  [(q/24.)*x_1**4. - E*I*d_1],
                                  [(q/24.)*x_2**4.],
                                  [(q/24.)*x_3**4. - E*I*d_3]
                                  ])

        solutionMatrix = np.linalg.solve(coefficientMatrix, resultMatrix)

        solutionArray = np.array(solutionMatrix)
        print "A_y = {a_y}, A_z = {a_z}, B_y = {b_y}, B_z = {b_z}, C_y = {c_y}, R_z = {r_z}, k_1 = {k_1}, k_2 = {k_2}".format(
            a_y=round(solutionArray[0, :][0], 2), a_z=round(solutionArray[1, :][0], 2),
            b_y=round(solutionArray[2, :][0], 2), b_z=round(solutionArray[3, :][0], 2),
            c_y=round(solutionArray[4, :][0], 2), r_z=round(solutionArray[5, :][0], 2),
            k_1=round(solutionArray[6, :][0], 2), k_2=round(solutionArray[7, :][0], 2))

        return solutionMatrix

    def bendingMomentAndShearForcesForX(self, A_y, A_z, B_y, B_z, C_y, R_z, x_3, x_2, x_1, x_a, q, P, x):

        m_zz = q*(1./2.)*x**2.
        v_y = q*x
        m_yy = 0.
        v_z = 0.

        if x >= x_1:
            m_zz -= A_y*(x-x_1)
            v_y -= A_y
            m_yy -= A_z*(x-x_1)
            v_z -= A_z

        if x >= (x_2 - (x_a/2.)):
            m_yy -= R_z*(x-x_2+(x_a/2.))
            v_z -= R_z

        if x >= x_2:
            m_zz -= B_y*(x-x_2)
            v_y -= B_y
            m_yy -= B_z*(x-x_2)
            v_z -= B_z

        if x >= (x_2 + (x_a/2.)):
            m_yy += P*(x-x_2-(x_a/2.))
            v_z += P

        if x >= x_3:
            m_zz -= C_y*(x-x_3)
            v_y -= C_y

        return [m_yy, m_zz, v_y, v_z]
