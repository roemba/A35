from parameters import parameters
from beamTheory.main import beamTheory
import matplotlib as mplt
import matplotlib.pyplot as plt
from CSGeo import *
import numpy as np

beam = beamTheory()
reactionForces = beam.staticEquilibrium(parameters.youngsmodulus, 0.000014925648, parameters.chord, parameters.span,
                                        parameters.xlocation1, parameters.xlocation2,
                                        parameters.xlocation3, parameters.d12, parameters.actuatorload,
                                        parameters.aerodynamicload,
                                        parameters.height, parameters.verticaldisplacementhinge3,
                                        parameters.verticaldisplacementhinge1).A1

ytab = []
xtab = np.linspace(0., parameters.span, num=500)
for x in xtab:
    internalForcesArray = beam.bendingMomentAndShearForcesForX(reactionForces[0], reactionForces[1], reactionForces[2],
                                                               reactionForces[3], reactionForces[4], reactionForces[5],
                                                               parameters.xlocation3, parameters.xlocation2,
                                                               parameters.xlocation1, parameters.d12,
                                                               parameters.aerodynamicload, parameters.actuatorload, x)
    ytab.append(internalForcesArray[3])

# plt.plot(xtab, ytab)
# plt.show()
