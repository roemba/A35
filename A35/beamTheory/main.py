import numpy as np

class beamTheory:
    @staticmethod
    def staticEquilibrium(E, I_zz, I_yy, c_a, l_a, x_1, x_2, x_3, x_a, P, q, h, d_3, d_1, delta):
        coefficientMatrix = np.matrix([[1., 0., 1., 0., 1., 0., 0., 0., 0., 0., 0.],
                                       [0., 1., 0., 1., 0., 1., 0., 0., 0., 0., 0.],
                                       [0., x_2-x_1, 0., 0., 0., x_a/2., 0., 0., 0., 0., 0.],
                                       [-(x_2-x_1), 0., 0., 0., x_3-x_2, 0., 0., 0., 0., 0., 0.],
                                       [0., 0., 0., 0., 0., (h/2.)*np.cos(delta), 0., 0., 0., 0., 0.],
                                       [0., 0., 0., 0., 0., 0., x_1, 1., 0., 0., 0.],
                                       [(1./6.)*(x_2-x_1)**3., 0., 0., 0., 0., 0., x_2, 1., 0., 0., 0.],
                                       [(1./6.)*(x_3-x_1)**3., 0., (1./6.)*(x_3-x_2)**3., 0., 0., 0., x_3, 1., 0., 0., 0.],
                                       [0., 0., 0., 0., 0., 0., 0., 0., x_1, 1., 0.],
                                       [0., (1./6.)*(x_2-x_1)**3., 0., 0., 0., (1./6.)*(x_a/2.)**3., 0., 0., x_2, 1., 0.],
                                       [0., (1./6.)*(x_3-x_1)**3., 0., (1./6.)*(x_3-x_2)**3, 0., (1./6.)*(x_3 - x_2 + (x_a/2.))**3., 0., 0., x_3, 1., E*I_yy],
                                    ])

        resultMatrix = np.matrix([[q*l_a],
                                  [P],
                                  [-P*(x_a/2.)],
                                  [0.],
                                  [P*(h/2.)*np.cos(delta)+q*l_a*(0.25*c_a-(h/2.))*np.cos(delta)],
                                  [(q/24.)*x_1**4. - E*I_zz*d_1],
                                  [(q/24.)*x_2**4.],
                                  [(q/24.)*x_3**4. - E*I_zz*d_3],
                                  [0],
                                  [0],
                                  [P*(1./6.)*(x_3 - x_2 - (x_a/2.))**3]
                                  ])

        solutionMatrix = np.linalg.solve(coefficientMatrix, resultMatrix)

        solutionArray = np.array(solutionMatrix)
        print "A_y = {a_y}, A_z = {a_z}, B_y = {b_y}, B_z = {b_z}, C_y = {c_y}, R_z = {r_z}, k_1_y = {k_1_y}, k_2_y = {k_2_y}, k_1_z = {k_1_z}, k_2_z = {k_2_z}, delta_c_z = {delta_c_z}".format(
            a_y=round(solutionArray[0, :][0], 2), a_z=round(solutionArray[1, :][0], 2),
            b_y=round(solutionArray[2, :][0], 2), b_z=round(solutionArray[3, :][0], 2),
            c_y=round(solutionArray[4, :][0], 2), r_z=round(solutionArray[5, :][0], 2),
            k_1_y=round(solutionArray[6, :][0], 2), k_2_y=round(solutionArray[7, :][0], 2),
            k_1_z=round(solutionArray[8, :][0], 2), k_2_z=round(solutionArray[9, :][0], 2),
            delta_c_z=solutionArray[10, :][0])

        return solutionMatrix

    @staticmethod
    def bendingMomentAndShearForcesForX(A_y, A_z, B_y, B_z, C_y, R_z, x_3, x_2, x_1, x_a, q, P, c_a, h, theta, x):

        m_zz = q*(1./2.)*x**2.
        v_y = q*x
        m_yy = 0.
        v_z = 0.
        T = q*x*(0.25*c_a-(h/2.))*np.cos(theta)

        if x >= x_1:
            m_zz -= A_y*(x-x_1)
            v_y -= A_y
            m_yy -= A_z*(x-x_1)
            v_z -= A_z

        if x >= (x_2 - (x_a/2.)):
            m_yy -= R_z*(x-x_2+(x_a/2.))
            v_z -= R_z
            T -= R_z*(h/2.)*np.cos(theta)

        if x >= x_2:
            m_zz -= B_y*(x-x_2)
            v_y -= B_y
            m_yy -= B_z*(x-x_2)
            v_z -= B_z

        if x >= (x_2 + (x_a/2.)):
            m_yy += P*(x-x_2-(x_a/2.))
            v_z += P
            T += P*(h/2.)

        if x >= x_3:
            m_zz -= C_y*(x-x_3)
            v_y -= C_y

        return [m_yy, m_zz, v_y, v_z]

    @staticmethod
    def calculateTorqueForX(q, P, R_z, A_z, A_y, B_z, B_y, C_y, c_a, h, theta, x_3, x_2, x_1, x_a, x_sc, x):
        T = q*x*(0.25*c_a-x_sc)*np.sin(theta)

        if x >= x_1:
            T -= (A_z*np.cos(theta) + A_y*np.sin(theta))*((h/2.) - x_sc)

        if x >= (x_2 - (x_a/2.)):
            T -= R_z*np.cos(theta)*((h/2.) - x_sc) + R_z*np.sin(theta)*(h/2.)

        if x >= x_2:
            T -= (B_z*np.cos(theta) + B_y*np.sin(theta))*((h/2.) - x_sc)

        if x >= (x_2 + (x_a/2.)):
            T += P*np.cos(theta)*((h/2.) - x_sc) + P*np.sin(theta)*(h/2.)

        if x >= x_3:
            T -= C_y*np.sin(theta)*((h/2.) - x_sc)

        return T
