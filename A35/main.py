
from parameters import parameters
from beamTheory.main import beamTheory
from CSGeo import *
import numpy as np

beam = beamTheory()
beam.staticEquilibrium(parameters.youngsmodulus, 0.000014925648, parameters.chord, parameters.span, parameters.xlocation1, parameters.xlocation2,
                       parameters.xlocation3, parameters.d12, parameters.actuatorload, parameters.aerodynamicload,
                       parameters.height, parameters.verticaldisplacementhinge3, parameters.verticaldisplacementhinge1)