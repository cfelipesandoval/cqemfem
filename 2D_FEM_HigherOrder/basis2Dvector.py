from basis2Dscalar import *
def basis2Dvector(n1, n2, n3, nEval1,nEval2, points):
    # inputs: 
    # n1, n2, n3 are the points of the triangle nodes 
    # nEval1 and nEval2 are the nodes (1, 2, or 3) to evaluate at, 1 < 2
    # points is a list of points to evaluate the function for 

    # output: 
    vectN = []

    scal1, grad1 = basis2Dscalar(n1, n2, n3, nEval1, points)
    scal2, grad2 = basis2Dscalar(n1, n2, n3, nEval2, points)

    import math 

    # if nodes are 1 and 2
    if ((nEval1 == 1 and nEval2 == 2) or (nEval1 == 2 and nEval2 == 1)):
        len = math.sqrt((n1[0] - n2[0])**2+(n1[1]-n2[1])**2)

    # if nodes are 1 and 3
    elif ((nEval1 == 1 and nEval2 == 3) or (nEval1 == 3 and nEval2 == 1)):
        len = math.sqrt((n1[0] - n3[0])**2+(n1[1]-n3[1])**2)

    # if nodes are 2 and 3
    elif ((nEval1 == 2 and nEval2 == 3) or (nEval1 == 3 and nEval2 == 2)):
        len = math.sqrt((n2[0] - n3[0])**2+(n2[1]-n3[1])**2)

    for i, val in enumerate(points):
        vectX = len*(scal1[i]* grad2[0] - scal2[i] * grad1[0])
        vectY = len*(scal1[i]* grad2[1] - scal2[i] * grad1[1])
        vectN.append([vectX, vectY])
        
    return vectN
