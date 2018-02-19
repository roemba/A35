
from parameters import parameters
from beamTheory.main import beamTheory

beam = beamTheory()
beam.staticEquilibrium(parameters.chord, parameters.span, parameters.xlocation1, parameters.xlocation2, parameters.xlocation3, parameters.d12, parameters.actuatorload, parameters.aerodynamicload)