import math

class boomAreas:
    boomArray = []
    m_y = 0.
    m_z = 0.
    I_yy = 0.
    I_zz = 0.

    """
    boomArray = [z, y, A, k, connection_type, l, connection_type, (m), connection_type, isStringer]
    where z = z-coordinate,
    y = y-coordinate,
    A = boom area,
    k = index of first connected boom
    l = index of second connected boom
    m = (optional) index of third connected boom
    connection_type = "skin" or "spar"
    
    """
    def __init__(self, boomArray, m_y, m_z, I_yy, I_zz):
        self.boomArray = boomArray
        self.m_y = m_y
        self.m_z = m_z
        self.I_yy = I_yy
        self.I_zz = I_zz

    def calculateBoomAreas(self, spar_thickness, skin_thickness, stringer_area):
        for index in xrange(len(self.boomArray)):
            boom = self.boomArray[index]
            if boom[9]:
                boomArea = stringer_area
            else:
                boomArea = 0.

            if self.m_z != 0.:
                ownDistanceToNeutralAxis = self.getDistanceToNeutralAxis(index)

            for i in [3, 5, 7]:
                if boom[i] is not None:
                    distanceToConnectedBoom = self.getSkinDistance(index, boom[i])

                    ratio = 0.
                    if self.m_z != 0.:
                        connectedBoomDistanceToNeutralAxis = self.getDistanceToNeutralAxis(boom[i])
                        if abs(ownDistanceToNeutralAxis) >= 0.000001 and abs(connectedBoomDistanceToNeutralAxis) >= 0.000001:
                            ratio = connectedBoomDistanceToNeutralAxis/ownDistanceToNeutralAxis
                        else:
                            ratio = 0.

                    if boom[i+1] == "spar":
                        boomArea += self.calculateBoomArea(spar_thickness, distanceToConnectedBoom, ratio)
                    else:
                        boomArea += self.calculateBoomArea(skin_thickness, distanceToConnectedBoom, ratio)

            self.boomArray[index][2] = boomArea

        return self.boomArray


    def getSkinDistance(self, index_boom_1, index_boom_2):
        boom1 = self.boomArray[index_boom_1]
        boom2 = self.boomArray[index_boom_2]
        return math.sqrt((boom2[0]-boom1[0])**2. + (boom2[1]-boom1[1])**2.)

    def getDistanceToNeutralAxis(self, index_boom):
        boom = self.boomArray[index_boom]
        z = boom[0]
        y = boom[1]
        return abs((self.I_yy/self.I_zz)*y - (self.m_y/self.m_z)*z)/math.sqrt((self.m_y/self.m_z)**2. + (self.I_yy/self.I_zz)**2.)

    def calculateBoomArea(self, thickness, distance, distance_ratio):
        return ((thickness*distance)/6.)*(2. + distance_ratio)
