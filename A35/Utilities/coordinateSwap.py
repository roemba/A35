import numpy as np
from math import *

class coordinateSwap:
    
    def APtoCS(self, z_AP, y_AP, theta):
        
        z_CS = z_AP*np.cos(theta) - y_AP*np.sin(theta)
        y_CS = z_AP*np.sin(theta) + y_AP*np.cos(theta)
        
        return z_CS, y_CS
    
    def CStoAP(self, z_CS, y_CS, theta):
        
        z_AP = y_CS*np.sin(theta) + z_CS*np.cos(theta)
        y_AP = y_CS*np.cos(theta) - z_CS*np.sin(theta)
        
        return z_AP, y_AP