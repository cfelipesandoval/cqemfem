from operator import itemgetter
import re
import sys
import time

def readAbq_ports(fName):
    # [nodes, elms, bndry, ports] = readAbq(fName)
    #
    #  INPUT: fName               ; filepath
    # OUTPUT: [nodes, elms, bndry]; [nodes] stores xyz-coords of each node, with first element as global node number
    #                               [elms] stores indices of [nodes] in first 3 cols, relative permittivity value in last col (CUBIT BLOCK NAME MUST BE "eprX")
    #                               [bndry] stores indices of [nodes]
    #                               [ports] stores indices of [nodes] along the lumped ports
    #
    # ERROR CODES:  Abq1 # Unrecognized element type: only curves B21/B31, tris STRI3, and tets C3D4 accepted.
    #               Abq2 # Unrecognized material property: name block "eprX" for permittivity of X, else do not rename for boundaries.
    #               Abq3 # Node numbers are not compressed: enter "compress all" into Cubit command line before exporting.
    #
    # NEW CODE DESIGNED FOR 3D

    file = open(fName, encoding='utf-8')
    nodes = []                                          # safe initializations
    elements = []
    bndry = []
    ports = []
    cmprsChk = 1                                        # verify node numbers are compressed

    while 1:                                            # move to NODES
        if (file.readline()).find("*NODE") != -1:
            break
            
    while 1:                                            # read out NODES
        currLine = file.readline()
        if currLine.find("*") != -1:                    # terminate if end reached
            break
        else:                                           # construct [nodes]
            currLine = currLine.split(",")
            if cmprsChk != int(currLine[0]):
                raise Exception("Node numbers are not compressed. code:Abq3")
            nodes.append([float(currLine[0]),           # global node index
                          float(currLine[1]),           # coordinates
                          float(currLine[2]),
                          float(currLine[3])])
            cmprsChk += 1

    while 1:                                            # move to ELEMENTS
        if currLine.find("E L E M E N T S") != -1:
            currLine = file.readline()
            break
        currLine = file.readline()
    
    while currLine.find("**") == -1 and len(currLine):  # read out ELEMENTS
        if currLine.find("ELEMENT,") != -1:             # element property line catch
            isBndry = currLine.split(",")               # extracting necessary information
            typeMat = isBndry[2].split("\n")
            typeMat = typeMat[0].strip()
            isBndry = isBndry[1].strip()
            portNum = 0
            match isBndry:                              # change if cubit decides to export a different element
                case "TYPE=C3D4" | "TYPE=B21" | "TYPE = B31":
                    elmDim = 1                          # elmDim for how many nodes needed in [elms] and [bndry]
                case "TYPE=STRI3":
                    elmDim = 0
                case _:
                    raise Exception("Unrecognized element type. code:Abq1")
            # extract material property
            if typeMat.find("epr") != -1 or typeMat.find("mpr") != -1:
                epr = re.compile("epr(\d+\.\d+|\d+)")   # difficult to search manually
                mpr = re.compile("mpr(\d+\.\d+|\d+)")   # hence regex
                typeMat = [float(epr.search(typeMat).group(1)),
                           float(mpr.search(typeMat).group(1))]
                isBndry = 0
            elif typeMat.find("Port") != -1:
                port = re.compile("Port(\d+\.\d+|\d+)")
                portNum = float(port.search(typeMat).group(1))
                isBndry = 0
            elif typeMat.find("EB") == -1:
                raise Exception("Unrecognized material property. code:Abq2")
            else:                                       # unnamed blocks are bndry blocks
                isBndry = 1
        else:
            currLine = currLine.split(",")
            if isBndry:                                 # constructs [bndry]: first elm is a flag if 1D boundary
                bndry.append([int(currLine[1]),
                              int(currLine[2]),
                              -1 if elmDim else int(currLine[3])])
            elif portNum != 0:
                ports.append([int(currLine[1]),
                              int(currLine[2]),
                              -1 if elmDim else int(currLine[3]),
                              int(portNum)])
            else:                                       # constructs [elements]: last elm is a flag if 2D mesh
                elements.append([int(currLine[0]),      # Cubit elm num (needed for AMR)
                                 int(currLine[1]),      # node indices
                                 int(currLine[2]),
                                 int(currLine[3]),
                                 int(currLine[4]) 
                                 if elmDim else -1,     # differentiate 2D/3D
                                 typeMat[0],            # electric permittivity
                                 typeMat[1]])           # magnetic permeability
        currLine = file.readline()
    
    bndry.sort(key=itemgetter(2))                       # sort from 2D to 3D bndry
    file.close()
    return [nodes, elements, bndry, ports]

def readAbq(fName):
    # [nodes, elms, bndry] = readAbq(fName)
    #
    #  INPUT: fName               ; filepath
    # OUTPUT: [nodes, elms, bndry]; [nodes] stores xyz-coords of each node, with first element as global node number
    #                               [elms] stores indices of [nodes] in first 3 cols, relative permittivity value in last col (CUBIT BLOCK NAME MUST BE "eprX")
    #                               [bndry] stores indices of [nodes]
    #
    # ERROR CODES:  Abq1 # Unrecognized element type: only curves B21/B31, tris STRI3, and tets C3D4 accepted.
    #               Abq2 # Unrecognized material property: name block "eprX" for permittivity of X, else do not rename for boundaries.
    #               Abq3 # Node numbers are not compressed: enter "compress all" into Cubit command line before exporting.
    #
    # NEW CODE DESIGNED FOR 3D

    file = open(fName, encoding='utf-8')
    nodes = []                                          # safe initializations
    elements = []
    bndry = []
    cmprsChk = 1                                        # verify node numbers are compressed

    while 1:                                            # move to NODES
        if (file.readline()).find("*NODE") != -1:
            break
            
    while 1:                                            # read out NODES
        currLine = file.readline()
        if currLine.find("*") != -1:                    # terminate if end reached
            break
        else:                                           # construct [nodes]
            currLine = currLine.split(",")
            if cmprsChk != int(currLine[0]):
                raise Exception("Node numbers are not compressed. code:Abq3")
            nodes.append([float(currLine[0]),           # global node index
                          float(currLine[1]),           # coordinates
                          float(currLine[2]),
                          float(currLine[3])])
            cmprsChk += 1

    while 1:                                            # move to ELEMENTS
        if currLine.find("E L E M E N T S") != -1:
            currLine = file.readline()
            break
        currLine = file.readline()
    
    while currLine.find("**") == -1 and len(currLine):  # read out ELEMENTS
        if currLine.find("ELEMENT,") != -1:             # element property line catch
            isBndry = currLine.split(",")               # extracting necessary information
            typeMat = isBndry[2].split("\n")
            typeMat = typeMat[0].strip()
            isBndry = isBndry[1].strip()
            match isBndry:                              # change if cubit decides to export a different element
                case "TYPE=C3D4" | "TYPE=B21" | "TYPE = B31":
                    elmDim = 1                          # elmDim for how many nodes needed in [elms] and [bndry]
                case "TYPE=STRI3":
                    elmDim = 0
                case _:
                    raise Exception("Unrecognized element type. code:Abq1")
            # extract material property
            if typeMat.find("epr") != -1 or typeMat.find("mpr") != -1:
                epr = re.compile("epr(\d+\.\d+|\d+)")   # difficult to search manually
                mpr = re.compile("mpr(\d+\.\d+|\d+)")   # hence regex
                typeMat = [float(epr.search(typeMat).group(1)),
                           float(mpr.search(typeMat).group(1))]
                isBndry = 0
            elif typeMat.find("EB") == -1:
                raise Exception("Unrecognized material property. code:Abq2")
            else:                                       # unnamed blocks are bndry blocks
                isBndry = 1
        else:
            currLine = currLine.split(",")
            if isBndry:                                 # constructs [bndry]: first elm is a flag if 1D boundary
                bndry.append([int(currLine[1]),
                              int(currLine[2]),
                              -1 if elmDim else int(currLine[3])])
            else:                                       # constructs [elements]: last elm is a flag if 2D mesh
                elements.append([int(currLine[0]),      # Cubit elm num (needed for AMR)
                                 int(currLine[1]),      # node indices
                                 int(currLine[2]),
                                 int(currLine[3]),
                                 int(currLine[4]) 
                                 if elmDim else -1,     # differentiate 2D/3D
                                 typeMat[0],            # electric permittivity
                                 typeMat[1]])           # magnetic permeability
        currLine = file.readline()
    
    bndry.sort(key=itemgetter(2))                       # sort from 2D to 3D bndry
    file.close()
    return [nodes, elements, bndry]

def readAbaqusMesh2D(fName):
    # [nodes, elms, bndry] = readAbaqusMesh2D(fName)
    #  INPUT:   fName               ; filepath
    # OUTPUT:   [nodes, elms, bndry]; [nodes] stores xyz-coords in 3 cols + arbitrary node number in last col (-1 flag for bndry)
    #                                 [elms] stores indices of [nodes] for each triangle vertex in 3 cols + block number in last col
    #                                 [bndry] stores boundary edge nodes 
    #
    # REPLICATES MATLAB OUTPUT

    file = open(fName, encoding='utf-8')
    nodes = []
    elms = []
    bndry = []
    idx = 1                                         # arbitrary node indexing
    bNum = 0                                        # block number
    temp = []                                       # temporary storage

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
                chk = 1 if int(currLine[0]) < int(currLine[1]) else 0
                bndry.append([int(currLine[0]) if chk else int(currLine[1]),
                              int(currLine[1]) if chk else int(currLine[0])])
            else:
                elms.append([int(currLine[0]),      # construct [elms]
                             int(currLine[1]),
                             int(currLine[2]),
                             bNum])
        currLine = file.readline()

    for lcv in range(len(bndry)):                   # prepare for next step
        temp.append(bndry[lcv][0])
        temp.append(bndry[lcv][1])

    for lcv in range(len(nodes)):                   # append flag/arbitrary indexing to [nodes]
        if lcv + 1 in temp:
            nodes[lcv].append(-1)
        else:
            nodes[lcv].append(idx)
            idx += 1

    file.close()
    bndry.sort(key=itemgetter(0))
    return [nodes, elms, bndry]

# DO NOT USE FF: ABSOLUTELY NONFUNCTIONAL
def meshRefine(fName, cubitPath, elmsList, E_eigs, B_eigs):
    # [n, e, b] = meshRefine(fName, cubitPath, elmsList, E_eigs, B_eigs)
    # 
    # INPUTS MAY CHANGE
    # fName = .cub filepath
    # cubitPath = local cubit \bin filepath
    # elmsList = list of elms
    # E_eigs, B_eigs = parameters to calculate error
    # 
    # outputs same as readAbq (or flag if no need to iterate?)
    #
    # ERROR CODES:  AMR1 # idklol
    
    # CHECK IF MATRICES ARE OPERABLE
    if len(E_eigs) != len(B_eigs): # PLACEHOLDER
        raise Exception("code: AMR1")

    elmErr = []

    # perform error calculations
    for elm in range(len(elmsList)):
        errCalc = errorGen(elm) # PLACEHOLDER
        elmErr.append(errCalc) # ASSUMES order of elmsList is identical to cubit list, hence implicit indexing

    tolerance = 0.5 # threshold to refine: USER-DEFINED
    print(f"total error: {sum(elmErr)}")
    if sum(elmErr) < tolerance:
        print("below tolerance, refinement cancelled")
        return
    
    tolerance = 0.5 # elements with more than 50% max error are refined: USER-DEFINED
    badElms = [elm + 1 for elm in range(len(elmErr)) if elmErr[elm] > tolerance * max(elmErr)]
        # obtain list of bad elements according to cubit indices
    print(badElms)

    # begin cubit operations
    sys.path.append(f'{cubitPath}')
    import cubit
    cubit.init
    cubit.cmd(f'open "{fName}"')
    # tolerance = 0.5 # if above 50% elements need refining: LIKELY UNNEEDED
    # if len(badElms) > tolerance * len(elmsList):
    #     # cubit.cmd(f"refine volume {} numsplit 1 smooth")
    #     # HOW TO QUERY VOLUMES TO REFINE: refine all is not a command
    #     pass
    # else:
    t = time.time()
    for idx in range(len(badElms)):
        cubit.cmd(f"refine element {badElms[idx]} depth 1 smooth")
        # tet/tri if needed
        # querying existence of elm may not be necessary: cubit just throws warning
    print(f"elapsed time: {time.time() - t}s")
    cubit.cmd("compress all")
    cubit.cmd("export abaqus 'C:/Users/dongg/Documents/Work/VIP279/temp.inp' mesh_only overwrite dimension 3")
    cubit.cmd("save cub5 'C:/Users/dongg/Documents/Work/VIP279/amr2af2.cub5' overwrite") #### testing
    [n, e, b] = readAbq('C:/Users/dongg/Documents/Work/VIP279/temp.inp')

    # p-adaptivity?

    return[n, e, b]

def errorGen(elm):
    return 5 if elm == 10331 or elm == 10340 or elm == 10349 else 0 #### testing

# testing
if __name__ == '__main__':
    # f = r"C:\Users\dongg\Documents\Work\VIP279\test3.inp"
    # f = 'cubicMesh.inp'
    # [n, e, b] = readAbq(f)
    # [n, e, b] = readAbaqusMesh2D(f)
    # fl = open('output2', 'w', encoding='utf-8')
    # for idx in range(len(b)):
    # out1 = str(n) + "\n\n" + str(e) + "\n\n" + str(b)
    # fl.write(out1)
    # fl.close
    f = r"C:\Users\dongg\Documents\Work\VIP279\amr2af1.cub5" #### testing
    meshRefine(f, r'D:\Coreform Cubit 2023.8\bin', [0] * 10486, [0], [0]) #### testing