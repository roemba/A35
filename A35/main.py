from parameters import parameters as pm
from beamTheory.main import beamTheory
from Utilities.coordinateSwap import coordinateSwap
from CSGeo.centroidLocation import Centroid
from CSGeo.boomAreas import boomAreas as ba
from CSGeo.numericalMOI import numericalMOI
from CSGeo.FEM import FEM
from shearAnalysis.openSectionShearFlow import openSectionShearFlow
from shearAnalysis.closedSectionShearFlow import closedSectionShearFlow
from shearAnalysis.totalShearFlow import shearFlowAndDeflection
import matplotlib.pyplot as plt
import numpy as np

saveFigs = False

# ------------------------------------------------------------------------------------
# Create a boom mesh
# ------------------------------------------------------------------------------------
initial2DBoomArray = FEM.boomPositions(5, 5, 5, pm.chord, pm.height, pm.stiffenernumber,
                                       pm.stiffenerheight, pm.stiffenerwidth,
                                       pm.stiffenerthickness)
stringer_area = (pm.stiffenerheight-pm.stiffenerthickness)*pm.stiffenerthickness + pm.stiffenerthickness*pm.stiffenerwidth

# ------------------------------------------------------------------------------------
# Draw figure of boom locations
# ------------------------------------------------------------------------------------
stringerBoom = np.array([0., 0.])
for boom in initial2DBoomArray:
    if boom[2] != 0:
        stringerBoom = np.vstack([stringerBoom, [boom[0], boom[1]]])

fig = plt.figure()
ax = fig.add_subplot(111)
ax.hlines(0, -0.6, 0.05)
ax.vlines(-0.225/2, -0.225/2 - 0.05, 0.225/2 + 0.05)
ax.vlines(0, -0.225/2 - 0.05, 0.225/2 + 0.05)
ax.set_xlabel("$z$ ($m$)")
ax.set_ylabel("$y$ ($m$)")
ax.grid()
ax.scatter(initial2DBoomArray[:, 0], initial2DBoomArray[:, 1], label="Boom")
ax.scatter(stringerBoom[:, 0], stringerBoom[:, 1], marker='*', c='red', label="Stringer")
ax.legend()
fig.show()

"""
boomArray = [z, y, A, [k, connection_type], [l, type], [(m), type]]

"""

# ------------------------------------------------------------------------------------
# Calculating the centroid analytically
#
# Setting reference coordinate system at LE
# (z coincides with chord, positive in direction of flight)
# (y positive upward flight direction)
# ------------------------------------------------------------------------------------
centroid_z = - Centroid.computeCentroid(pm.chord, pm.height, pm.skinthickness,
                                        pm.sparthickness, pm.stiffenerthickness,
                                        pm.stiffenerwidth, pm.stiffenerheight,
                                        pm.stiffenernumber)
centroid_y = 0.

# ------------------------------------------------------------------------------------
# Calculate the initial moment of inertia for the iteration
# ------------------------------------------------------------------------------------
i_zz_cs, i_yy_cs, i_zy_cs = numericalMOI.getMOI(centroid_z, centroid_y, initial2DBoomArray)

n_of_crossections = 500
riblocations = [pm.xlocation1, pm.xlocation2 - (pm.d12 / 2.), pm.xlocation2 + (pm.d12 / 2.), pm.xlocation3]

xtab = np.linspace(0., pm.span, num=n_of_crossections)
xtab = np.concatenate((xtab, riblocations))
xtab = np.sort(xtab)

boomArray3D = []
for x in xtab:
    boomArray3D.append([initial2DBoomArray, [i_zz_cs, i_yy_cs, i_zy_cs]])

# ------------------------------------------------------------------------------------
# Start iterating the moment of inertia until convergence has been achieved
# ------------------------------------------------------------------------------------
ytab = []
i_zz_cs_old = 0.00000000000000001
i = 0
while abs((i_zz_cs/i_zz_cs_old)*100. - 100.) > 1.:
#while i <= 10:
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

        boomAreaClass = ba(boomArray, internalForcesArray[0], internalForcesArray[1], i_zz_cs, i_yy_cs, centroid_z, centroid_y)
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

# ------------------------------------------------------------------------------------
# Start final calculation of reaction forces using analytical MoI and shear center
# ------------------------------------------------------------------------------------
i_zz_cs = 0.00001153484934
i_yy_cs = 0.00006369185039
z_sc = 0.0634 + pm.height/2.

reactionForces = beamTheory.staticEquilibrium(pm.youngsmodulus, i_zz_cs, pm.chord,
                                              pm.span,
                                              pm.xlocation1, pm.xlocation2,
                                              pm.xlocation3, pm.d12, pm.actuatorload,
                                              pm.aerodynamicload,
                                              pm.height, pm.verticaldisplacementhinge3,
                                              pm.verticaldisplacementhinge1, pm.maxupwarddeflection).A1

crosssections = []
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
                                            z_sc, x)

    m_zz_cs, m_yy_cs = coordinateSwap.APtoCS(externalForcesArray[1], externalForcesArray[0],
                                             pm.maxupwarddeflection)
    v_z_cs, v_y_cs = coordinateSwap.APtoCS(externalForcesArray[3], externalForcesArray[2],
                                           pm.maxupwarddeflection)
    internalForcesArray = [m_yy_cs, m_zz_cs, v_y_cs, v_z_cs, torque]
    crosssections.append(internalForcesArray)

    boomAreaClass = ba(boomArray, internalForcesArray[0], internalForcesArray[1], i_zz_cs, i_yy_cs, centroid_z, centroid_y)
    boomArray3D[index][0] = boomAreaClass.calculateBoomAreas(pm.sparthickness, pm.skinthickness, stringer_area)

# ------------------------------------------------------------------------------------
# Draw figures of bending, shear and torque
# ------------------------------------------------------------------------------------

npYArray = np.array(crosssections)
fig2 = plt.figure()
ax2 = fig2.add_subplot(111)
ax2.set_title("Bending moments for $n = " + str(n_of_crossections) + "$")
ax2.set_xlim(0, pm.span)
ax2.plot(xtab, npYArray[:, 0], label="$M_{y_{cs}}$")
ax2.plot(xtab, npYArray[:, 1], label="$M_{z_{cs}}$")
ax2.grid(b=True, which='both', color='0.65', linestyle='-')
ax2.legend()
ax2.set_xlabel("Span ($m$)")
ax2.set_ylabel("Moment ($Nm$)")
fig2.show()

fig3 = plt.figure(figsize=(7,4.8))
ax3 = fig3.add_subplot(111)
ax3.set_title("Shear forces for $n = " + str(n_of_crossections) + "$")
ax3.set_xlim(0, pm.span)
ax3.plot(xtab, npYArray[:, 2], label="$V_{y_{cs}}$")
ax3.plot(xtab, npYArray[:, 3], label="$V_{z_{cs}}$")
ax3.grid(b=True, which='both', color='0.65', linestyle='-')
ax3.legend()
ax3.set_xlabel("Span ($m$)")
ax3.set_ylabel("Force ($N$)")
fig3.show()

fig4 = plt.figure()
ax4 = fig4.add_subplot(111)
ax4.set_title("Torque for $n = " + str(n_of_crossections) + "$")
ax4.set_xlim(0, pm.span)
ax4.plot(xtab, npYArray[:, 4], label="$T_{CS}$")
ax4.grid(b=True, which='both', color='0.65', linestyle='-')
ax4.legend()
ax4.set_xlabel("Span ($m$)")
ax4.set_ylabel("Torque ($Nm$)")
fig4.show()

# ------------------------------------------------------------------------------------
# Start calculating shear flows
# ------------------------------------------------------------------------------------

numberOfBoomsArray = FEM.discretization(5, 5, 5, pm.chord, pm.height, pm.stiffenernumber)

openShearFlowCrossSections = []
closedShearFlowCrossSections = []
displacements = []
normalStress = []
fig5 = plt.figure(figsize=(8, 4.8))
for index in xrange(xtab.shape[0]):
    boomArray = boomArray3D[index][0]
    x = xtab[index]

    openShearFlow = openSectionShearFlow.totOpenCalc(i_zz_cs, i_yy_cs, 0., crosssections[index][3], crosssections[index][2],
                                                     boomArray, int(numberOfBoomsArray[0][0]), int(numberOfBoomsArray[1][0]),
                                                     int(numberOfBoomsArray[3][0]))

    openShearFlowCrossSections.append(openShearFlow)


    q_s01, q_s02, d_theta_d_x, loopArray = closedSectionShearFlow.calculation(crosssections[index][4], pm.chord, pm.height, pm.skinthickness, pm.sparthickness,
                                                         z_sc, boomArray, openShearFlow, pm.shearmodulus)

    q_max, idqmax, id2qmax, max_shear_stress, idxstressMax, idx2stressMax = shearFlowAndDeflection.crossSectionMaxQ(q_s01, q_s02, boomArray, openShearFlow, pm.skinthickness, pm.sparthickness)

    closedShearFlow = [q_s01, q_s02, d_theta_d_x, q_max, idqmax, id2qmax, max_shear_stress, idxstressMax, idx2stressMax]
    closedShearFlowCrossSections.append(closedShearFlow)

    d_y = shearFlowAndDeflection.bendingDeflectionInY(pm.youngsmodulus, i_zz_cs, pm.maxupwarddeflection,
                                                         pm.aerodynamicload, pm.xlocation1, pm.xlocation2,
                                                         pm.xlocation3, pm.d12, reactionForces[0], reactionForces[1],
                                                         reactionForces[2], reactionForces[3], reactionForces[4],
                                                         reactionForces[5], reactionForces[6], reactionForces[7],
                                                         pm.actuatorload, x)

    k1_z, k2_z = shearFlowAndDeflection.integrationConstantsY(pm.maxupwarddeflection,
                                                              pm.aerodynamicload, pm.xlocation1, pm.xlocation2,
                                                              pm.d12, reactionForces[0], reactionForces[1], reactionForces[5])

    d_z = shearFlowAndDeflection.bendingDeflectionInZ(pm.youngsmodulus, i_yy_cs, pm.maxupwarddeflection,
                                                      pm.aerodynamicload, pm.xlocation1, pm.xlocation2,
                                                      pm.xlocation3, pm.d12, reactionForces[0], reactionForces[1],
                                                      reactionForces[2], reactionForces[3], reactionForces[4],
                                                      reactionForces[5], k1_z, k2_z,
                                                      pm.actuatorload, x)
    # Normal stresses
    maxStress, minStress, maxStressBoomIdx, minStressBoomIdx = \
                                        shearFlowAndDeflection.maxNormalStress(boomArray, crosssections[index][0],
                                                                            crosssections[index][1], i_yy_cs, i_zz_cs)
    # Taken min separately, since compression has additional interesting effects..
    # max takes abs(min) as well
    normalStress.append([maxStress, minStress])
    displacements.append([d_y, d_z])

    if index == 5:
        ax5 = fig5.add_subplot(111)
        ax5.hlines(0, -0.6, 0.05)
        ax5.vlines(-0.225/2, -0.225/2 - 0.05, 0.225/2 + 0.05)
        ax5.vlines(0, -0.225/2 - 0.05, 0.225/2 + 0.05)
        ax5.set_xlabel("$z$ ($m$)")
        ax5.set_ylabel("$y$ ($m$)")
        ax5.grid()
        ax5.scatter(initial2DBoomArray[:, 0], initial2DBoomArray[:, 1], label="Boom")
        ax5.scatter(stringerBoom[:, 0], stringerBoom[:, 1], marker='*', c='red', label="Stringer")
        ax5.legend()
        for i in xrange(initial2DBoomArray.shape[0]):
            ax5.annotate(i, (initial2DBoomArray[i, 0], initial2DBoomArray[i, 1]))
        for shearArray in openShearFlow:
            startBoomIndex = int(shearArray[1])
            startBoom = boomArray[startBoomIndex]
            endBoomIndex = int(shearArray[2])
            endBoom = boomArray[endBoomIndex]
            ax5.plot([startBoom[0], endBoom[0]], [startBoom[1], endBoom[1]])

        fig5.show()

npClosedShearFlowCrossSections = np.array(closedShearFlowCrossSections)
fig6 = plt.figure()
ax6 = fig6.add_subplot(111)
ax6.set_title("Rate of twist for $n = " + str(n_of_crossections) + "$")
ax6.set_xlim(0, pm.span)
ax6.plot(xtab, npClosedShearFlowCrossSections[:, 2], label=r"$\frac{d\theta}{d x}_{CS}}$")
ax6.grid(b=True, which='both', color='0.65', linestyle='-')
ax6.legend()
ax6.set_xlabel("Span ($m$)")
ax6.set_ylabel(r"$\frac{d\theta}{d x}$ ($rad m^{-1}$)")
fig6.show()

fig7 = plt.figure(figsize=(8, 4.8))
ax7 = fig7.add_subplot(111)
ax7.set_title("Maximum shear flow")
ax7.set_xlim(0, pm.span)
for index in xrange(len(riblocations)):
    lb = ""
    if index == 0:
        lb = "Rib 1-4"
    riblocation = riblocations[index]
    ax7.vlines(riblocation, 0, np.max(npClosedShearFlowCrossSections[:, 3]), label=lb, linestyles="dashed")
    xtab_index = np.where(xtab == riblocation)[0]
    ax7.annotate('%.2e' % round(float(npClosedShearFlowCrossSections[xtab_index, 3]), 2), (riblocation, npClosedShearFlowCrossSections[xtab_index, 3] + 3000.))
ax7.plot(xtab, npClosedShearFlowCrossSections[:, 3], label="$q_{max_{CS}}$")
ax7.grid(b=True, which='both', color='0.65', linestyle='-')
ax7.legend()
ax7.set_xlabel("Span ($m$)")
ax7.set_ylabel("Shear flow ($Nm^{-1}$)")
fig7.show()

fig12 = plt.figure(figsize=(8, 4.8))
ax12 = fig12.add_subplot(111)
ax12.set_title("Maximum shear stress")
ax12.set_xlim(0, pm.span)
for index in xrange(len(riblocations)):
    lb = ""
    if index == 0:
        lb = "Rib 1-4"
    riblocation = riblocations[index]
    ax12.vlines(riblocation, 0, np.max(npClosedShearFlowCrossSections[:, 6]), label=lb, linestyles="dashed")
    xtab_index = np.where(xtab == riblocation)[0]
    ax12.annotate('%.2e' % round(float(npClosedShearFlowCrossSections[xtab_index, 6]), 2), (riblocation, npClosedShearFlowCrossSections[xtab_index, 6] + 3000000.))
ax12.plot(xtab, npClosedShearFlowCrossSections[:, 6], label=r"$\tau_{max_{CS}}$")
ax12.grid(b=True, which='both', color='0.65', linestyle='-')
ax12.legend()
ax12.set_xlabel("Span ($m$)")
ax12.set_ylabel("Shear stress ($Pa$)")
fig12.show()

npDisplacements = np.array(displacements)
fig8 = plt.figure()
ax8 = fig8.add_subplot(111)
ax8.set_title("Displacement in $y$")
ax8.set_xlim(0, pm.span)
ax8.plot(xtab, npDisplacements[:, 0], label=r"$\delta_{CS}}$")
ax8.grid(b=True, which='both', color='0.65', linestyle='-')
ax8.legend()
ax8.set_xlabel("Span ($m$)")
ax8.set_ylabel("Displacement ($m$)")
fig8.show()

npNormStress = np.array(normalStress)
fig10 = plt.figure()
ax10 = fig10.add_subplot(111)
ax10.set_title("Maximum Normal Stress")
ax10.set_xlim(0, pm.span)
ax10.annotate("Max = " + '%.2e' % np.max(npNormStress[:, 0]), (xtab[np.where(npNormStress[:,0] == np.max(npNormStress[:, 0]))[0]], np.max(npNormStress[:, 0])))
ax10.plot(xtab, npNormStress[:, 0], label=r"$\sigma_{n_{CS}}$")
ax10.grid(b=True, which='both', color='0.65', linestyle='-')
ax10.legend()
ax10.set_xlabel("Span ($m$)")
ax10.set_ylabel("Normal stress ($Pa$)")
fig10.show()

actuator1Index = np.where(xtab == (pm.xlocation2 - pm.d12/2.))[0][0]
# LE_x0, LE_y0 = shearFlowAndDeflection.rotatePoint(-pm.maxupwarddeflection, pm.height/2., 0.)
# TE_x0, TE_y0 = shearFlowAndDeflection.rotatePoint(-pm.maxupwarddeflection, pm.height/2. - pm.chord, 0.)
# coordinates_backwards = []
# for index in xrange(actuator1Index, 0, -1):
#     x0 = xtab[index]
#     d_y = displacements[index][0]
#     x1 = xtab[index - 1]
#     d_theta_d_x0 = closedShearFlowCrossSections[index][2]
#     theta_1 = abs(x1-x0)*d_theta_d_x0
#
#     LE_x0, LE_y0 = shearFlowAndDeflection.rotatePoint(-(theta_1), LE_x0, LE_y0)
#     TE_x0, TE_y0 = shearFlowAndDeflection.rotatePoint(-(theta_1), TE_x0, TE_y0)
#     coordinates_backwards.append([x0, LE_x0, LE_y0+d_y, TE_x0, TE_y0+d_y])

LE_x0, LE_y0 = shearFlowAndDeflection.rotatePoint(-pm.maxupwarddeflection, pm.height/2., 0.)
TE_x0, TE_y0 = shearFlowAndDeflection.rotatePoint(-pm.maxupwarddeflection, pm.height/2. - pm.chord, 0.)
coordinates_forwards = []
for index in xrange(0, xtab.shape[0] - 1, 1):
    x0 = xtab[index]
    d_y = displacements[index][0]
    x1 = xtab[index + 1]
    d_theta_d_x0 = closedShearFlowCrossSections[index][2]
    theta_1 = abs(x1-x0)*d_theta_d_x0

    LE_x0, LE_y0 = shearFlowAndDeflection.rotatePoint(-(theta_1), LE_x0, LE_y0)
    TE_x0, TE_y0 = shearFlowAndDeflection.rotatePoint(-(theta_1), TE_x0, TE_y0)
    coordinates_forwards.append([x0, LE_x0, LE_y0+d_y, TE_x0, TE_y0+d_y])

npCoordinates_forwards = np.array(coordinates_forwards)
# npCoordinates_backwards = np.array(coordinates_backwards)
fig9 = plt.figure()
ax9 = fig9.add_subplot(111)
ax9.set_title("LE and TE Displacement in $y$ with initial deflection, bending and twist")
ax9.set_xlim(0, pm.span)
ax9.vlines(pm.xlocation2 - pm.d12/2., -0.10, 0.30, linestyles="dashed", label="Actuator 1")
ax9.plot(npCoordinates_forwards[:, 0], npCoordinates_forwards[:, 2], label=r"LE Deflection $\delta_{LE_{AP}}$", color="b")
# ax9.plot(npCoordinates_backwards[:, 0], npCoordinates_backwards[:, 2], color="b")
ax9.plot(npCoordinates_forwards[:, 0], npCoordinates_forwards[:, 4], label=r"TE Deflection $\delta_{TE_{AP}}$", color="r")
# ax9.plot(npCoordinates_backwards[:, 0], npCoordinates_backwards[:, 4], color="r")
ax9.grid(b=True, which='both', color='0.65', linestyle='-')
ax9.legend()
ax9.set_xlabel("Span ($m$)")
ax9.set_ylabel("Displacement ($m$)")
fig9.show()

fig11 = plt.figure()
ax11 = fig11.add_subplot(111)
ax11.set_title("LE and TE Displacement in $z$ with initial deflection, bending and twist")
ax11.set_xlim(0, pm.span)
ax11.vlines(pm.xlocation2 - pm.d12/2., -0.5, 0.30, linestyles="dashed", label="Actuator 1")
ax11.plot(npCoordinates_forwards[:, 0], npCoordinates_forwards[:, 1], label=r"LE Deflection $\delta_{LE_{AP}}$", color="b")
# ax11.plot(npCoordinates_backwards[:, 0], npCoordinates_backwards[:, 1], color="b")
ax11.plot(npCoordinates_forwards[:, 0], npCoordinates_forwards[:, 3], label=r"TE Deflection $\delta_{TE_{AP}}$", color="r")
# ax11.plot(npCoordinates_backwards[:, 0], npCoordinates_backwards[:, 3], color="r")
ax11.grid(b=True, which='both', color='0.65', linestyle='-')
ax11.legend()
ax11.set_xlabel("Span ($m$)")
ax11.set_ylabel("Displacement ($m$)")
fig11.show()

if saveFigs:
    fig.savefig("boom_locations.png")
    fig2.savefig("bending_moments.png")
    fig3.savefig("shear_forces.png")
    fig4.savefig("torque.png")
    fig5.savefig("panel_connections.png")
    fig6.savefig("rate_of_twist.png")
    fig7.savefig("q_max.png")
    fig8.savefig("y_displacement_cs.png")
    fig9.savefig("LE_TE_displacement_y.png")
    fig10.savefig("Max_normal_stress.png")
    fig11.savefig("LE_TE_displacement_z.png")
    fig12.savefig("Max_shear_stress.png")
