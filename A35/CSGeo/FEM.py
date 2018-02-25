
import numpy as np
import matplotlib.pyplot as plt

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

            if key:
                print 'Simple program is not advanced enough for complicated mesh-making; ' \
                      'reducing Sector', sectorNumber, 'boom amount to:', n_booms

            # From first to last boom, spacing is equal at: stringerPitch * stringerAmount / (n_booms - 4)
            # this changes for the first set to     startDist / (boomsPerStringer + 5)
            # and for the last set it changes to    endDist / (boomsPerStringer + 5)

            boomsPerStringer = (n_booms - 4) / stringerAmount

            normPitch = stringerPitch * stringerAmount / (n_booms - 4)
            startPitch = 2 * startDist / (boomsPerStringer + 5)
            endPitch = 2 * endDist / (boomsPerStringer + 5)

        else:       # Because sector 4 (spar) doesn't have any stringers, and equal spacing is the current go-to
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

   #     print 'sector 1, 3'
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
        endDist = startDist     # Due to symmetry

        n_booms, normPitch, startPitch, endPitch, stringerAmount = FEM.generalSector(n_booms, stringerPitch, startDist,
                                                                                     endDist, semi + startDist,
                                                                                     minBooms, name)

        # Due to discretization accuracy issues, the return statements for pitch are changed to radians
        radius = semi / np.pi
        normAngle = normPitch / radius
        normPitch = normAngle

        return n_booms, normPitch, startDist, endDist, stringerAmount



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

     #   print 'sector 4'
      #  print n_booms, normPitch, startPitch, endPitch, stringerAmount
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

        # Sector 2:     # sector 2 'start' and 'end' pitch modified to dist instead
        n_sector_2, normPitch_2, startPitch_2, endPitch_2, stringersIn2 = \
            FEM.sector_2(n_sector_2, stringerPitch, hypo, semi, stringersIn1)

        # Sector 4:
        n_sector_4, normPitch_4, startPitch_4, endPitch_4, stringersIn4 = FEM.sector_4(n_sector_4, h)


        # Sector 1, 3:
        outputArray[0] = n_sector_1, normPitch_1, startPitch_1, endPitch_1, stringersIn1
        outputArray[2] = n_sector_3, normPitch_3, startPitch_3, endPitch_3, stringersIn3
        outputArray[1] = n_sector_2, normPitch_2, startPitch_2, endPitch_2, stringersIn2
        outputArray[3] = n_sector_4, normPitch_4, startPitch_4, endPitch_4, stringersIn4

#        print ''
 #       print outputArray
  #      print ''
        # Quick ref: output format: [   sectorNumber, normPitch, startPitch, endPitch, stringerAmount  ]
        return outputArray



    # Let's get this party started
    @staticmethod
    def boomPositions(n_sector_1, n_sector_2, n_sector_4, C_a, h, n_st, h_st, w_st, t_st):
        """
        Takes array of boom pitch values; spits out coordinates compliant with assumptions
        :return: [z_pos, y_pos, area, [k, l, m]]
        """

        # From choices:     boom idx 0                              TE          ->  [-C_a, h]
        #                   boom idx n_sector_1 + 1                 Spar neg. y ->  [-h/2, -h/2]
        #                   boom idx n_sector_2 + n_sector_1 + 2    Spar pos. y ->  [-h/2, h/2]
        #
        # values k, l are almost entirely direct index neighbours in outArray return

        # First step: get boom pitch values:
        discrArray = FEM.discretization(n_sector_1, n_sector_2, n_sector_4, C_a, h, n_st)

        # Geometry values:
        radius, longS, semi, hypo, angle, length, stringerPitch = FEM.crossSection(h, C_a, n_st)
        stringerArea = (h_st + w_st) * t_st

        # Array creation
        # outArray = np.vstack([A, [z_pos, y_pos, area, k, k_type, l, l_type, m, m_type]])
        # boom idx 0:
        m, m_type = None, None
        outArray = np.array(
            [-C_a, int(0), int(0), int(1), 'skin', int(2 * n_sector_1 + n_sector_2 + 2), 'skin', m, m_type])

        # sector looping:
        for sector in range(4):
            # Fixed is the value m, since the only non-'N/A' values are in fixed points 2 and 3
            booms, normPitch, startPitch, endPitch, stringers = discrArray[sector]
            if stringers == 0:
                normBoomNumber = 0  # There is one case where stringers = 0
            else:
                normBoomNumber = ((booms - 4) / stringers - 1) / 2  # Follows assumption

            interStringerBooms = 2 * normBoomNumber  # just like the name
            isStringer = False
            stringerTicker = 0  # ticks through each cycle, is used to determine when isStringer is flipped for area
            collPitch = 0.      # reset per sector: collective boom pitch since
            startBooms = True  # very obvious
            endBooms = False
            if startPitch == 0.: startBooms = False
            isFirstStringer = True

            # The following variables are for rework of sector 2 discretization/spacing without influencing all else
            # This rework makes the code look like crazy
            sec2Pitch = 0.
            if booms % 2 != 0:
                tempBooms = (booms - 1) / 2  # amount of booms in one quadrant minus boom in LE
                notStartBooms = ((stringers - 1) / 2) * (interStringerBooms + 1)
            else:
                tempBooms = booms / 2
                notStartBooms = (stringers / 2) * (interStringerBooms + 1) - normBoomNumber
            normAngleSpacing = ((semi / 2.) - startPitch) / radius / notStartBooms
            startEndSpacing = startPitch / radius / (normBoomNumber + 2)
            semi2Ticker = 0

            for boom in range(int(booms)):
                if startBooms or ((not startBooms) and isFirstStringer):
                    pitch = startPitch
                    if isStringer and isFirstStringer:
                        isFirstStringer = False
                elif endBooms:
                    pitch = endPitch
                else:
                    pitch = normPitch
                if sector + 1 == 4:     # prevent difficult code requirements
                    pitch = normPitch

                if sector + 1 == 4:   # getting weird stuff; prevention
                    area = int(0)
                elif isStringer:
                    area = stringerArea
                else:
                    area = int(0)

                # Step per iter
                collPitch += pitch

                if sector + 1 == 1:  # TE -> lower spar intersect
                    z = collPitch * np.cos(angle) - C_a
                    y = - collPitch * np.sin(angle)
                    print pitch, collPitch, z, y

                # Major rework happened here; messy is a result
                if sector + 1 == 2:
                    if boom <= tempBooms:       # to ensure that this implies the first half only
                        if tempBooms - notStartBooms == boom + 1: pitchUsed = startEndSpacing
                        else: pitchUsed = normAngleSpacing            # spacing between spar, and first stringer
                        sec2Pitch += pitchUsed
                        if booms % 2 != 0 and boom == tempBooms:    # hardcoding is navigates around issue
                            z = 0.
                            y = 0.
                        else:
                            z = radius * (np.sin(sec2Pitch) - 1.)
                            y = radius * np.cos(sec2Pitch)

                    else:                       # Second half literally copies first half; works due to symmetry y
                        semi2Ticker += 1  # to work towards index of referenced values
                        if booms % 2 != 0:
                            symmCopy = outArray[- 2 * semi2Ticker]
                        else:
                            symmCopy = outArray[- 2 * semi2Ticker + 1]

                        z = symmCopy[0]
                        y = -1 * symmCopy[1]

                if sector + 1 == 3:  # upper spar intersect -> TE
                    z = - radius - collPitch * np.cos(angle)
                    y = radius - collPitch * np.sin(angle)

                if sector + 1 == 4:  # lower spar intersect -> upper spar intersect
                    z = - radius
                    y = - radius + collPitch

                # end of second nested loop:
                stringerTicker += 1
                if isStringer:
                    isStringer = False
                    stringerTicker = 0

                if boom + 1 == normBoomNumber + 2:
                    startBooms = False
                    isStringer = True  # because the assumption is how it is designed
                if boom + 1 == booms - (normBoomNumber + 2):
                    if not (sector == 4):  # Let's prevent weird issues shall we?
                        endBooms = True

                if (not (startBooms or endBooms)) and stringerTicker == interStringerBooms:
                    isStringer = True

                # Let's get the neighbouring boom index
                k = int(np.shape(outArray)[0] - 1)  # takes amount of rows. THIS is the issue if there is any
                l = int(k + 2)

                if sector + 1 != 4:
                    k_type = 'skin'
                    l_type = 'skin'

                else:
                    k_type = 'spar'
                    l_type = 'spar'

                if sector + 1 == 3 and boom == range(int(booms)): l = 0
                outArray = np.vstack([outArray, [z, y, area, k, k_type, l, l_type, m, m_type]])

    ###################################################################################################################
    #################################################### SEPARATOR ####################################################
    ###################################################################################################################

            # How else do we get those spar caps in there?
            if sector + 1 == 1:
                outArray = np.vstack(
                    [outArray, [-h / 2, -h / 2, 0, n_sector_1, 'skin', n_sector_1 + 2, 'skin', - n_sector_4, 'spar']])
            if sector + 1 == 2:
                outArray = np.vstack(
                    [outArray, [-h / 2, h / 2, 0, n_sector_2 + n_sector_1 + 1, 'skin', n_sector_2 + n_sector_1 + 3, 'skin', -1, 'spar']])

        # Not sure why this is necessary, but it is.
        for boomIndex in range(len(outArray)):
            for itemIndex in [0, 1, 2, 3, 5, 7]:
                item = outArray[boomIndex, itemIndex]
                if item is not None:
                    if float(item) == 0.0:
                        outArray[boomIndex, itemIndex] = int(item)
                    else:
                        outArray[boomIndex, itemIndex] = float(item)

        return outArray


# OUTPUT -- [   z   ,   y   ,   area    ,   ... neighbours ...  ]

testArray = FEM.boomPositions(5, 20, 5, 0.547, 0.225, 17, 0.015, 0.02, 0.0012)
#print 'testarray', testArray
#print testArray[12+19]
#print '\n spar', testArray[-5:]
stringerBoom = np.array([0.05, 0.])

partArray = testArray[11:12+20, :3]#[12:12+19, :3]
for boom in testArray:
    if boom[2] != 0:
        stringerBoom = np.vstack([stringerBoom, [boom[0], boom[1]]])


#print testArray[3, 0:3]

#print '\n stringerBoom'
#print stringerBoom

z = testArray[12:12+19, 0]
zStringer = stringerBoom[:, 0]
y = testArray[12:12+19, 1]
yStringer = stringerBoom[:, 1]

z = testArray[:, 0]
y = testArray[:, 1]

plt.hlines(0, -0.6, 0.05)
plt.vlines(-0.225/2, -0.225/2 - 0.05, 0.225/2 + 0.05)
plt.vlines(0, -0.225/2 - 0.05, 0.225/2 + 0.05)
plt.scatter(z, y)
plt.scatter(stringerBoom[:, 0], stringerBoom[:, 1], marker='*', c='red')
plt.show()

# Current problems:
# Sectors 1, 3 do not display correctly as they should in scatter plot
# Sector 2 shows spacing incorrect -> potential fix: take (n_booms / 2) - 1 and div over spacing, set fix at LE?
