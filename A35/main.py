from parameters import parameters as pm
from beamTheory.main import beamTheory
from Utilities.coordinateSwap import coordinateSwap
from CSGeo.centroidLocation import Centroid
from CSGeo.areaMomentOfInertia import areaMomentOfInertia as AMOI
from CSGeo.boomAreas import boomAreas as ba
from CSGeo.numericalMOI import numericalMOI
from CSGeo.FEM import FEM
import matplotlib as mplt
import matplotlib.pyplot as plt
import numpy as np

#Replace these values with initial condition
# i_zz_cs = 0.000014925648
# i_yy_cs = 0.000014925648
# i_zy_cs = 0.
initial2DBoomArray = FEM.boomPositions(5, 5, 5, pm.chord, pm.height, pm.stiffenernumber,
                                       pm.stiffenerheight, pm.stiffenerwidth,
                                       pm.stiffenerthickness)
stringer_area = (pm.stiffenerheight-pm.stiffenerthickness)*pm.stiffenerthickness + pm.stiffenerthickness*pm.stiffenerwidth
print pm.verticaldisplacementhinge1, pm.verticaldisplacementhinge3

"""
boomArray = [z, y, A, [k, connection_type], [l, type], [(m), type]]

"""

# Setting reference coordinate system at LE
# (z coincides with chord, positive in direction of flight)
# (y positive upward flight direction)
centroid_z = - Centroid.computeCentroid(pm.chord, pm.height, pm.skinthickness,
                                        pm.sparthickness, pm.stiffenerthickness,
                                        pm.stiffenerwidth, pm.stiffenerheight,
                                        pm.stiffenernumber)
centroid_y = 0.

i_zz_cs, i_yy_cs, i_zy_cs = numericalMOI.getMOI(centroid_z, centroid_y, initial2DBoomArray)

n_of_crossections = 500
xtab = np.linspace(0., pm.span, num=n_of_crossections)

boomArray3D = []
for x in xtab:
    boomArray3D.append([initial2DBoomArray, [i_zz_cs, i_yy_cs, i_zy_cs]])

#Start loop
ytab = []
i_zz_cs_old = 0.00000000000000001
i = 0
while abs((i_zz_cs/i_zz_cs_old)*100. - 100.) > 1.:
#while i <= 100:
    i_zz_cs_old = i_zz_cs
    i += 1
    print "Loop: " + str(i)
    reactionForces = beamTheory.staticEquilibrium(pm.youngsmodulus, i_zz_cs, pm.chord,
                                                  pm.span,
                                                  pm.xlocation1, pm.xlocation2,
                                                  pm.xlocation3, pm.d12, pm.actuatorload,
                                                  pm.aerodynamicload,
                                                  pm.height, pm.verticaldisplacementhinge3,
                                                  pm.verticaldisplacementhinge1, pm.maxupwarddeflection).A1

    i_zz_sum = 0.
    i_yy_sum = 0.
    i_zy_sum = 0.
    #Starting the iterations over all the cross sections
    for index in xrange(xtab.shape[0]):
        boomArray = boomArray3D[index][0]
        x = xtab[index]

        externalForcesArray = beamTheory.bendingMomentAndShearForcesForX(reactionForces[0], reactionForces[1],
                                                                         reactionForces[2],
                                                                         reactionForces[3], reactionForces[4],
                                                                         reactionForces[5],
                                                                         pm.xlocation3, pm.xlocation2,
                                                                         pm.xlocation1, pm.d12,
                                                                         pm.aerodynamicload,
                                                                         pm.actuatorload, pm.chord,
                                                                         pm.height, pm.maxupwarddeflection,
                                                                         x)

        torque = beamTheory.calculateTorqueForX(pm.aerodynamicload, pm.actuatorload, reactionForces[5],
                                                reactionForces[1], reactionForces[0], reactionForces[3], reactionForces[2],
                                                reactionForces[4], pm.chord, pm.height,
                                                pm.maxupwarddeflection,
                                                pm.xlocation3, pm.xlocation2, pm.xlocation1,
                                                pm.d12,
                                                pm.height / 2., x)

        m_zz_cs, m_yy_cs = coordinateSwap.APtoCS(externalForcesArray[1], externalForcesArray[0],
                                                 pm.maxupwarddeflection)
        v_z_cs, v_y_cs = coordinateSwap.APtoCS(externalForcesArray[3], externalForcesArray[2],
                                               pm.maxupwarddeflection)
        internalForcesArray = [m_yy_cs, m_zz_cs, v_y_cs, v_z_cs]

        boomAreaClass = ba(boomArray, internalForcesArray[0], internalForcesArray[1], i_zz_cs, i_yy_cs)
        boomArray3D[index][0] = boomAreaClass.calculateBoomAreas(pm.sparthickness, pm.skinthickness, stringer_area)

        #Calculate MOI for each cross section
        boomArray3D[index][1][0], boomArray3D[index][1][1], boomArray3D[index][1][2] = numericalMOI.getMOI(centroid_z, centroid_y, boomArray3D[index][0])
        i_zz_sum += boomArray3D[index][1][0]
        i_yy_sum += boomArray3D[index][1][1]
        i_zy_sum += boomArray3D[index][1][2]

    i_zz_cs = i_zz_sum/len(xtab)
    i_yy_cs = i_yy_sum/len(xtab)
    i_zy_cs = i_zy_sum/len(xtab)
    ytab.append(i_zz_cs)

plotxtable = range(i)
fig1 = plt.figure()
ax1 = fig1.add_subplot(111)
ax1.plot(plotxtable, ytab)
ax1.set_xlabel("n of iterations")
fig1.show()
print "Result I_zz_cs:", i_zz_cs, "Result I_yy_cs:", i_yy_cs, "Result I_zy_cs:", i_zy_cs
print

#Calculate all the shear forces and bending moments using analytical MoI
i_zz_cs = 0.00001153484934
i_yy_cs = 0.00006369185039

reactionForces = beamTheory.staticEquilibrium(pm.youngsmodulus, i_zz_cs, pm.chord,
                                              pm.span,
                                              pm.xlocation1, pm.xlocation2,
                                              pm.xlocation3, pm.d12, pm.actuatorload,
                                              pm.aerodynamicload,
                                              pm.height, pm.verticaldisplacementhinge3,
                                              pm.verticaldisplacementhinge1, pm.maxupwarddeflection).A1

ytab = []
for index in xrange(xtab.shape[0]):
    boomArray = boomArray3D[index][0]
    x = xtab[index]

    externalForcesArray = beamTheory.bendingMomentAndShearForcesForX(reactionForces[0], reactionForces[1],
                                                                     reactionForces[2],
                                                                     reactionForces[3], reactionForces[4],
                                                                     reactionForces[5],
                                                                     pm.xlocation3, pm.xlocation2,
                                                                     pm.xlocation1, pm.d12,
                                                                     pm.aerodynamicload,
                                                                     pm.actuatorload, pm.chord,
                                                                     pm.height, pm.maxupwarddeflection,
                                                                     x)

    torque = beamTheory.calculateTorqueForX(pm.aerodynamicload, pm.actuatorload, reactionForces[5],
                                            reactionForces[1], reactionForces[0], reactionForces[3], reactionForces[2],
                                            reactionForces[4], pm.chord, pm.height,
                                            pm.maxupwarddeflection,
                                            pm.xlocation3, pm.xlocation2, pm.xlocation1,
                                            pm.d12,
                                            pm.height / 2., x)

    m_zz_cs, m_yy_cs = coordinateSwap.APtoCS(externalForcesArray[1], externalForcesArray[0],
                                             pm.maxupwarddeflection)
    v_z_cs, v_y_cs = coordinateSwap.APtoCS(externalForcesArray[3], externalForcesArray[2],
                                           pm.maxupwarddeflection)
    internalForcesArray = [m_yy_cs, m_zz_cs, v_y_cs, v_z_cs]
    ytab.append(internalForcesArray)

    boomAreaClass = ba(boomArray, internalForcesArray[0], internalForcesArray[1], i_zz_cs, i_yy_cs)
    boomArray3D[index][0] = boomAreaClass.calculateBoomAreas(pm.sparthickness, pm.skinthickness, stringer_area)

npYArray = np.array(ytab)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_title("Bending moments for " + str(n_of_crossections) + " cross sections")
ax2.plot(xtab, npYArray[:, 0], label="M_yy_cs")
ax2.plot(xtab, npYArray[:, 1], label="M_zz_cs")
ax2.grid(b=True, which='both', color='0.65', linestyle='-')
ax2.legend()
ax2.set_xlabel("x")
fig2.show()

fig3 = plt.figure()
ax3 = fig3.add_subplot(111)
ax3.set_title("Shear forces for " + str(n_of_crossections) + " cross sections")
ax3.plot(xtab, npYArray[:, 2], label="V_y_cs")
ax3.plot(xtab, npYArray[:, 3], label="V_z_cs")
ax3.legend()
ax3.set_xlabel("x")
fig3.show()
