
import numpy as np

class FEM:
    """
    Borrowed stringerPos function from CSGeometry to have this standalone for call and output

    Output array output has z and y coordinates. pitch and boom area are
    assumed constants; in loop the boom area may be modified to match
    z and y are positive according to reference frame (z in dir of flight, y upwards)
    Reference point for axis system is at Leading Edge.
    """

    """
    Produce output in form: [z_pos, y_pos, [k, l (, m)]]
    """

    @staticmethod
    def crossSection(h, C_a, n_st):
        radius = h / 2.
        longS = C_a - radius
        semi = radius * np.pi
        hypo = np.sqrt(radius**2 + longS**2)
        angle = np.arctan(radius / longS)
        length = semi + 2. * hypo
        pitch = length / n_st
        return radius, longS, semi, hypo, angle, length, pitch



    @staticmethod
    def generalSector(n_booms, stringerPitch, startDist, endDist, length, minBooms, sectorNumber):
        """
        General sector boom spacing function that places booms centered about stringers
        with two additional booms in the boundary panels

        :return:                n_boom, normPitch, startPitch, endPitch, stringerAmount
        """
        if sectorNumber != '4':
            stringerAmount = int(length / stringerPitch)

            # First solve issues with less than minimum booms
            if n_booms - 4 > minBooms:
                if n_booms < stringerAmount + 4:
                    n_booms = stringerAmount + 4  # Add four for two additional at each corner case
                    print 'Sector', sectorNumber, 'does not have the minimum amount of booms; it is now:', n_booms
            else:
                n_booms = minBooms + 4
                print 'Sector', sectorNumber, 'does not have the minimum amount of booms; it is now:', n_booms

            # Then solve such that booms are placed at the stringers
            key = False
            if (n_booms - 4) % stringerAmount != 0:
                n_booms -= (n_booms - 4) % stringerAmount
                key = True
            if ((n_booms - 4) / stringerAmount) % 2 == 0:
                n_booms -= stringerAmount
                key = True
            if key == True:
                print 'Simple program is not advanced enough for complicated mesh-making; ' \
                      'reducing Sector', sectorNumber, 'boom amount to:', n_booms

            # From first to last boom, spacing is equal at: stringerPitch * stringerAmount / (n_booms - 4)
            # this changes for the first set to     startDist / ((boomsPerStringer + 3) / 2)
            # and for the last set it changes to    endDist / ((boomsPerStringer + 3) / 2)

            boomsPerStringer = (n_booms - 4) / stringerAmount

            normPitch = stringerPitch * stringerAmount / (n_booms - 4)
            startPitch = startDist / ((boomsPerStringer + 3) / 2)
            endPitch = endDist / ((boomsPerStringer + 3) / 2)

        else:
            normPitch = length / n_booms
            startPitch = 0
            endPitch = 0
            stringerAmount = 0

        return n_booms, normPitch, startPitch, endPitch, stringerAmount



    @staticmethod
    def sector_1(n_booms, stringerPitch, hypo):
        """
        Sector 1, 3: TE to y connection spar
        """

        name = '1, 3'
        minBooms = 0
        stringerAmount = int(hypo / stringerPitch)
        endDist = hypo - stringerPitch * (stringerAmount - 0.5)
        startDist = stringerPitch / 2

        n_booms, normPitch, startPitch, endPitch, stringerAmount = FEM.generalSector(n_booms, stringerPitch, startDist,
                                                                                     endDist, hypo, minBooms, name)

        print 'sector 1, 3'
        print n_booms, normPitch, startPitch, endPitch, stringerAmount
        return n_booms, normPitch, startPitch, endPitch, stringerAmount



    @staticmethod
    def sector_2(n_booms, stringerPitch, hypo, semi, stringersIn1):
        """
        Sector 2: Semi-circular LE
        """
        name = '2'
        minBooms = 17   # small angle approximation requires 18 panels for 180 degrees

        # now check for spacing things
        # As with sector 1, 3: require skin-booms symmetrically about stringers
        # first: get stringers in section...
        startDist = stringerPitch * (stringersIn1 + 1) - hypo
        stringerAmount = int((semi + startDist) / stringerPitch)
        endDist = startDist     # Due to symmetry

        n_booms, normPitch, startPitch, endPitch, stringerAmount = FEM.generalSector(n_booms, stringerPitch, startDist,
                                                                                         endDist, semi + startDist, minBooms, name)

        print 'sector 2'
        print n_booms, normPitch, startPitch, endPitch, stringerAmount
        return n_booms, normPitch, startPitch, endPitch, stringerAmount



    @staticmethod
    def sector_4(n_booms, height):
        """
        Sector 4: Spar
        """
        name = '4'
        minBooms = 0
        stringerPitch = None

        startDist = 0
        endDist = 0

        n_booms, normPitch, startPitch, endPitch, stringerAmount = FEM.generalSector(n_booms, stringerPitch, startDist,
                                                                                     endDist, height, minBooms, name)

        print 'sector 4'
        print n_booms, normPitch, startPitch, endPitch, stringerAmount
        return n_booms, normPitch, startPitch, endPitch, stringerAmount



    @staticmethod
    def discretization(n_sector_1, n_sector_2, n_sector_4, C_a, h, n_st):
        """
        Sector 1: TE to negative y connection spar
        Sector 2: Semi-circular LE
        Sector 3: positive y connection spar to TE: equal spacing to Sector 1
        Sector 4: Spar
        Fixed boom 1: TE
        Fixed boom 2: negative y spar-skin connection
        Fixed boom 3: positive y spar-skin connection
        """

        # get stringer pitch values
        radius, longS, semi, hypo, angle, length, stringerPitch = FEM.crossSection(h, C_a, n_st)

        # Decision time:
        # Each stringer has an amount of skin-booms on each side, equal spacing within (pitch/2) of itself
        # Each boundary condition (except spar?) contains two additional booms

        outputArray = np.zeros((4, 5))

        # Check each sector value to match spacing of stringers as equally as possible
        # Sectors 1, 3:n_booms, normPitch, startPitch, endPitch, stringerAmount

        n_sector_1, normPitch_1, startPitch_1, endPitch_1, stringersIn1 = FEM.sector_1(n_sector_1, stringerPitch, hypo)
        n_sector_3, normPitch_3, endPitch_3, startPitch_3, stringersIn3 =\
            n_sector_1, normPitch_1, startPitch_1, endPitch_1, stringersIn1

        # Sector 2:
        n_sector_2, normPitch_2, startPitch_2, endPitch_2, stringersIn2 = \
            FEM.sector_2(n_sector_2, stringerPitch, hypo, semi, stringersIn1)

        # Sector 4:
        n_sector_4, normPitch_4, startPitch_4, endPitch_4, stringersIn4 = FEM.sector_4(n_sector_4, h)


        # Sector 1, 3:
        outputArray[0] = n_sector_1, normPitch_1, startPitch_1, endPitch_1, stringersIn1
        outputArray[2] = n_sector_3, normPitch_3, startPitch_3, endPitch_3, stringersIn3

        outputArray[1] = n_sector_2, normPitch_2, startPitch_2, endPitch_2, stringersIn2

        # Sector 4:
        outputArray[3] = n_sector_4, normPitch_4, startPitch_4, endPitch_4, stringersIn4

        print ''
        return outputArray

print FEM.discretization(18, 17, 5, 0.547, 0.225, 17)
#b = FEM.discretization(77, 52, 5, 0.547, 0.225, 17)


