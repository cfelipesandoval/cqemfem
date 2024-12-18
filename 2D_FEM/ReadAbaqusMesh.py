def readAbq(fName):
    # [nodes, elms, bndry] = readAbq(fName)
    #  INPUT:   fName            ; filepath
    # OUTPUT:   [nodes, elms, bndry]; [nodes] stores xyz-coords of each node
    #                                 [elms] stores indices of [nodes] in first 3 cols, relative permittivity value in last col (CUBIT BLOCK NAME MUST BE "epr=X")
    #                                 [bndry] stores indices of [nodes]

    file = open(fName, encoding='utf-8')
    nodes = []
    elements = []
    bndry = []

    while 1:                                        # move to NODES
        if (file.readline()).find("*NODE") != -1:
            break
            
    while 1:                                        # read out NODES
        currLine = file.readline()
        if currLine.find("*") != -1:                # terminate if end reached
            break
        else:                                       # construct [nodes]
            currLine = currLine.split(",")
            nodes.append([float(currLine[0]),
                          float(currLine[1]),
                          float(currLine[2]),
                          float(currLine[3])])

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
            match isBndry:                          # cases may change
                case "TYPE=STRI3":
                    isBndry = 0
                case "TYPE=B21":
                    isBndry = 1
                case _:
                    print("ERROR! code:Abq1")       # Unrecognized element type: only tris and boundaries accepted.
                    return None
            if typeMat.find("epr") != -1:           # extracts relative permittivity supporting decimals
                typeMat = typeMat.split("epr")
                typeMat = (typeMat.pop()).split("_")
                temp = float(typeMat.pop(0))
                try:                                # scuffed try-except to catch any decimals
                    typeMat = temp + float(typeMat[0]) / 10
                except:
                    typeMat = temp
            elif typeMat.find("EB") == -1:
                print("ERROR! code:Abq2")           # Unrecognized material property: name block "epr=X" for permittivity
                return None
        else:
            currLine = currLine.split(",")
            if isBndry:                             # constructs [bndry]
                bndry.append([int(currLine[1]),
                              int(currLine[2])])
            else:                                   # constructs [elements]
                elements.append([typeMat,
                                 int(currLine[1]),
                                 int(currLine[2]),
                                 int(currLine[3])])
        currLine = file.readline()
    
    file.close()
    return [nodes, elements, bndry]

def readAbaqusMesh2D(fName):
    # [nodes, elms, bndry] = readAbaqusMesh2D(fName)
    #  INPUT:   fName               ; filepath
    # OUTPUT:   [nodes, elms, bndry]; [nodes] stores xyz-coords in 3 cols + arbitrary node number in last col (-1 flag for bndry)
    #                                 [elms] stores indices of [nodes] for each triangle vertex in 3 cols + block number in last col
    #                                 [bndry] stores boundary edge nodes 

    file = open(fName, encoding='utf-8')
    nodes = []
    elms = []
    bndry = []
    idx = 1                                         # arbitrary node indexing
    bNum = 0                                        # block number: IDEALLY EPR

    while 1:                                        # move to NODES
        if (file.readline()).find("*NODE") != -1:
            break
            
    while 1:                                        # read out NODES
        currLine = file.readline()
        if currLine.find("*") != -1:                # terminate if end reached
            break
        else:
            currLine = currLine.split(",")
            nodes.append([float(currLine[1]),       # construct [nodes]
                          float(currLine[2]),
                          float(currLine[3])])
    
    while 1:                                        # move to ELEMENTS
        currLine = file.readline()
        if currLine.find("E L E M E N T S") != -1:
            currLine = file.readline()
            break
              
    while currLine.find("**") == -1:                # read out ELEMENTS
        if currLine.find("ELEMENT,") != -1:         # element property line catch
            isBndry = currLine.split(",")
            typeMat = isBndry[2].split("\n")        # extract block number
            typeMat = typeMat[0].strip()
            isBndry = isBndry[1].strip()            # determines if bndry
            match isBndry:                          # cases may change
                case "TYPE=STRI3":
                    isBndry = 0
                case "TYPE=B21":
                    isBndry = 1
                case _:
                    print("ERROR! code:Abq1")       # Unrecognized element type: only tris and boundaries accepted.
                    return None
            if typeMat.find("epr") != -1:           # assigns ARBITRARY block number
                bNum += 1
        else:
            currLine = currLine.split(",")
            currLine.pop(0)
            if isBndry:                             # construct [bndry]
                bndry.append([int(currLine[0]),
                              int(currLine[1])])
            else:
                elms.append([int(currLine[0]),      # construct [elms]
                             int(currLine[1]),
                             int(currLine[2]),
                             bNum])
        currLine = file.readline()

    for lcv in range(len(nodes)):                   # append flag/arbitrary indexing to [nodes]
        if lcv + 1 in bndry:
            nodes[lcv].append(-1)
        else:
            nodes[lcv].append(idx)
            idx += 1

    file.close()
    return [nodes, elms, bndry]