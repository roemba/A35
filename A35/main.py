from parameters import parameters
from beamTheory.main import beamTheory
from Utilities.coordinateSwap import coordinateSwap
from CSGeo.centroidLocation import Centroid
from CSGeo.areaMomentOfInertia import areaMomentOfInertia as AMOI
import matplotlib as mplt
import matplotlib.pyplot as plt
from CSGeo import *
import numpy as np

reactionForces = beamTheory.staticEquilibrium(parameters.youngsmodulus, 0.000014925648, 0.000014925648, parameters.chord,
                                              parameters.span,
                                              parameters.xlocation1, parameters.xlocation2,
                                              parameters.xlocation3, parameters.d12, parameters.actuatorload,
                                              parameters.aerodynamicload,
                                              parameters.height, parameters.verticaldisplacementhinge3,
                                              parameters.verticaldisplacementhinge1, parameters.maxupwarddeflection).A1

ytab = []
xtab = np.linspace(0., parameters.span, num=500)
# xtab = [0.0765]
for x in xtab:
    externalForcesArray = beamTheory.bendingMomentAndShearForcesForX(reactionForces[0], reactionForces[1],
                                                                     reactionForces[2],
                                                                     reactionForces[3], reactionForces[4],
                                                                     reactionForces[5],
                                                                     parameters.xlocation3, parameters.xlocation2,
                                                                     parameters.xlocation1, parameters.d12,
                                                                     parameters.aerodynamicload,
                                                                     parameters.actuatorload, parameters.chord,
                                                                     parameters.height, parameters.maxupwarddeflection,
                                                                     x)

    torque = beamTheory.calculateTorqueForX(parameters.aerodynamicload, parameters.actuatorload, reactionForces[5],
                                            reactionForces[1], reactionForces[0], reactionForces[3], reactionForces[2],
                                            reactionForces[4], parameters.chord, parameters.height,
                                            parameters.maxupwarddeflection,
                                            parameters.xlocation3, parameters.xlocation2, parameters.xlocation1,
                                            parameters.d12,
                                            parameters.height / 2., x)

    m_zz_cs, m_yy_cs = coordinateSwap.APtoCS(externalForcesArray[1], externalForcesArray[0],
                                             parameters.maxupwarddeflection)
    v_z_cs, v_y_cs = coordinateSwap.APtoCS(externalForcesArray[3], externalForcesArray[2],
                                           parameters.maxupwarddeflection)
    internalForcesArray = [m_yy_cs, m_zz_cs, v_y_cs, v_z_cs]
    ytab.append(torque)

# plt.plot(xtab, ytab)
# plt.show()

# Setting centroid position for simple calling etc.
# Setting reference coordinate system at LE
# (z coincides with chord, positive in direction of flight)
# (y positive upward flight direction)
centroid_z = - Centroid.computeCentroid(parameters.chord, parameters.height, parameters.skinthickness, parameters.sparthickness, parameters.stiffenerthickness, parameters.stiffenerwidth, parameters.stiffenerheight, parameters.stiffenernumber)
centroid_y = 0.
print centroid_z

# Area moment of inertias
I_zz, I_yy = AMOI.getAMOI(parameters.chord, parameters.height, parameters.stiffenerwidth, parameters.stiffenerheight, parameters.stiffenernumber, parameters.stiffenerthickness, parameters.skinthickness, parameters.sparthickness, centroid_z)
print "I_zz", I_zz * 1000.**4
print "I_yy", I_yy * 1000.**4
>>>>>>> 7f06039c9ea557ed80c007e956d494254998ee94
