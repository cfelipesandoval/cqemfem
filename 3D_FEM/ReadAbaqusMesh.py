from operator import itemgetter

def readAbq(fName):
    # [nodes, elms, bndry] = readAbq(fName)
    #
    #  INPUT: fName               ; filepath
    # OUTPUT: [nodes, elms, bndry]; [nodes] stores xyz-coords of each node
    #                               [elms] stores indices of [nodes] in first 3 cols, relative permittivity value in last col (CUBIT BLOCK NAME MUST BE "eprX")
    #                               [bndry] stores indices of [nodes]
    #
    # ERROR CODES:  Abq1 # Unrecognized element type: only trimesh and boundaries accepted.
    #               Abq2 # Unrecognized material property: name block "eprX" for permittivity of X, else do not rename for boundaries.
    #               Abq3 # Node numbers are not compressed (fix??): enter "compress all" into Cubit command line before exporting.


    file = open(fName, encoding='utf-8')
    nodes = []
    elements = []
    bndry = []
    cmprsChk = 1                                    # verify node numbers are compressed

    while 1:                                        # move to NODES
        if (file.readline()).find("*NODE") != -1:
            break
            
    while 1:                                        # read out NODES
        currLine = file.readline()
        if currLine.find("*") != -1:                # terminate if end reached
            break
        else:                                       # construct [nodes]
            currLine = currLine.split(",")
            if cmprsChk != int(currLine[0]):
                raise Exception("code:Abq3")
            nodes.append([float(currLine[0]),
                          float(currLine[1]),
                          float(currLine[2]),
                          float(currLine[3])])
            cmprsChk += 1

    while 1:                                        # move to ELEMENTS
        if currLine.find("E L E M E N T S") != -1:
            currLine = file.readline()
            break
        currLine = file.readline()
    
    while currLine.find("**") == -1 and len(        # read out ELEMENTS
        currLine):
        if currLine.find("ELEMENT,") != -1:         # element property line catch
            isBndry = currLine.split(",")
            typeMat = isBndry[2].split("\n")
            typeMat = typeMat[0].strip()
            isBndry = isBndry[1].strip()
            match isBndry:                          # DEAL WITH CUBIT IMPORT MISBEHAVIOR
                case "TYPE=C3D4" | "TYPE=B21" | "TYPE = B31":                    
                    isBndry = 0                      # elmDim for how many nodes needed in [elms] and [bndry]
                case "TYPE=STRI3":
                    isBndry = 1
                case _:
                    raise Exception("code:Abq1")
            if typeMat.find("epr") != -1:           # extracts relative permittivity
                typeMat = typeMat.split("epr")
                typeMat = (typeMat.pop()).split("_")
                temp = float(typeMat.pop(0))
                try:                                # scuffed try-except to catch any decimals
                    typeMat = temp + float(
                        typeMat[0]) / 10
                except:
                    typeMat = temp
                # isBndry = 0
            elif typeMat.find("EB") == -1:
                raise Exception("code:Abq2")
            # else:
            #     isBndry = 1
        else:
            currLine = currLine.split(",")
            if isBndry:                             # constructs [bndry]: first elm is a flag if 1D boundary
                bndry.append([int(currLine[1]),
                              int(currLine[2]),
                              int(currLine[3])])
            else:                                   # constructs [elements]: last elm is a flag if 2D mesh
                elements.append([int(currLine[1]),
                                 int(currLine[2]),
                                 int(currLine[3]),
                                 int(currLine[4])]) 
                
            # output different dimension objects in different arrays??
        currLine = file.readline()
    
    file.close()
    # bndry.sort(key=itemgetter(2))   # sort from 2D to 3D bndry
    return [nodes, elements, bndry]