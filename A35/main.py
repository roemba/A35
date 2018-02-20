
from parameters import parameters
from beamTheory.main import beamTheory
import numpy as np

beam = beamTheory()
verticalLoads = beam.displacementStepFunctions(parameters.youngsmodulus, 0.00001240432379,
                                               parameters.verticaldisplacementhinge3,
                                               parameters.verticaldisplacementhinge1, parameters.xlocation3,
                                               parameters.xlocation2, parameters.xlocation1)

beam.staticEquilibrium(parameters.chord, parameters.span, parameters.xlocation1, parameters.xlocation2,
                       parameters.xlocation3, parameters.d12, parameters.actuatorload, parameters.aerodynamicload,
                       np.array(verticalLoads)[0, :][0], np.array(verticalLoads)[2, :][0], parameters.height)
