
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
        radius = h /2.
        longS = C_a - radius
        semi = radius * np.pi
        hypo = np.sqrt(radius**2 + longS**2)
        angle = np.arctan(radius / longS)
        length = semi + 2. * hypo
        pitch = length / n_st
        return radius, longS, semi, hypo, angle, length, pitch

    @staticmethod
    def sector_1(n_booms, stringerPitch, hypo):
        """
        Sector 1: TE to negative y connection spar
        """
        stringerAmount = int(hypo / stringerPitch)

        if n_booms < stringerAmount:  # check that boom amount is at least the amount of stringers
            n_booms = stringerAmount
            print 'Sectors 1 and 3 do not have the minimum amount of booms; it is now:', n_booms

        # note that the total distance to be divided should be (stringerAmount - 1) * stringerPitch + stringerPitch/2
        # which leaves out the residual towards boom in the spar:
        lastDist = hypo - stringerPitch * (stringerAmount - 0.5)

        # check if n_sector_1 % stringerAmount
        if n_booms % stringerAmount != 0:
            n_booms -= n_booms % stringerAmount
        if (n_booms / stringerAmount) % 2 == 0:
            n_booms -= stringerAmount
            print 'Simple program not advanced enough for complicated mesh-making; reducing ' \
                  'sector 1 and 3 boom amount to:', n_booms

        # Where until last boom (near spar), spacing is equal at: stringerPitch * stringerAmount / n_sector_1
        # this changes for the last set to ((n_sector_1 / stringerAmount) - 1) lastDist / 2
        normPitch = stringerPitch * stringerAmount / n_booms
        lastPitch = lastDist / ((n_booms / stringerAmount - 1) / 2)

        print 'sector 1, 3'
        print n_booms, normPitch, lastPitch
        return n_booms, normPitch, lastPitch, stringerAmount

    @staticmethod
    def sector_2(n_booms, stringerPitch, hypo, semi, stringersIn1):
        """
        Sector 2: Semi-circular LE
        """
        startDist = stringerPitch * (stringersIn1 + 1) - hypo

        # small angle approximation is required for linearization of inter-boom skin
        # which leads to a minimum of 18 skin panels..
        # check number of booms to equal 17 or more first:
        if n_booms < 17:
            n_booms = 17
            print 'To allow linearization of panels between the booms, 18 skin panels are required;' \
                  'the minimum boom count is thus:', n_booms

        # now check for spacing things
        # As with sector 1, 3: require skin-booms symmetrically about stringers
        # first: get stringers in section...
        stringerAmount = int((semi + startDist) / stringerPitch)

        if n_booms % stringerAmount != 0:
            n_booms -= n_booms % stringerAmount
        if (n_booms / stringerAmount) % 2 == 0:
            n_booms -= stringerAmount
            print 'Simple program not advanced enough for complicated mesh-making; reducing ' \
                  'sector 2 boom amount to:', n_booms

        normPitch = stringerPitch * stringerAmount / n_booms
        startEndPitch = startDist / ((n_booms / stringerAmount - 1) / 2)

        print 'sector 2'
        print n_booms, normPitch, startEndPitch

        return n_booms, normPitch, startEndPitch


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

        # Check each sector value to match spacing of stringers as equally as possible
        # Sectors 1, 3:
        n_sector_1, normPitch_1, lastPitch_1, stringersIn1 = FEM.sector_1(n_sector_1, stringerPitch, hypo)

        # Sector 2:
        n_sector_2, normPitch_2, lastPitch_2 = FEM.sector_2(n_sector_2, stringerPitch, hypo, semi, stringersIn1)



        # Hey guys, this shit isn't finished yet. far from it even




        outputArray = np.array([])
        print ''
        return outputArray

a = FEM.discretization(18, 17, 0, 0.547, 0.225, 17)
b = FEM.discretization(77, 52, 0, 0.547, 0.225, 17)


