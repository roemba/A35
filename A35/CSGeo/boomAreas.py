import math

class boomAreas:
    boomArray = []
    m_y = 0.
    m_z = 0.
    I_yy = 0.
    I_zz = 0.

    """
    boomArray = [z, y, A, [k, connection_type], [l, type], [(m), type]]
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

    def calculateBoomAreas(self, spar_thickness, skin_thickness):
        for index, boom in self.boomArray:
            boomArea = 0.
            ownDistanceToNeutralAxis = self.getDistanceToNeutralAxis(index)

            for i in range(3, 6):
                if boom[i][0] is not None:
                    distanceToConnectedBoom = self.getSkinDistance(index, boom[i][0])
                    connectedBoomDistanceToNeutralAxis = self.getDistanceToNeutralAxis(boom[i][0])

                    print connectedBoomDistanceToNeutralAxis/ownDistanceToNeutralAxis
                    if boom[i][1] == "spar":
                        boomArea += self.calculateBoomArea(spar_thickness, distanceToConnectedBoom,
                                                           connectedBoomDistanceToNeutralAxis/ownDistanceToNeutralAxis)
                    else:
                        boomArea += self.calculateBoomArea(skin_thickness, distanceToConnectedBoom,
                                                           connectedBoomDistanceToNeutralAxis/ownDistanceToNeutralAxis)

            print index, boomArea
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
        return ((self.I_yy/self.I_zz)*y - (self.m_y/self.m_z)*z)/math.sqrt((self.m_y/self.m_z)**2. + (self.I_yy/self.I_zz)**2.)

    def calculateBoomArea(self, thickness, distance, distance_ratio):
        return ((thickness*distance)/6.)*(2. + distance_ratio)
