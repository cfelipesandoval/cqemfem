# Function evaluates 2D scalar basis functions and their derivatives 
#   over a single triangle

def basis2Dscalar(n1, n2, n3, nodeEval, points):
    # inputs:
    # n1, n2, and n3 are the points of the triangle vertices
    # nodeEval is 1, 2, or 3, which node the function is evaluated at
    # points is an array of the points (x,y) to evaluate the function at

    # outputs:
    fvals = [] # array of function values at given points
    gradvals = [] # array of derivative (gradient) function value at given node
    
    import numpy as np # numpy used for matrix determinant 

    # get area of triangle with determinant 
    triMatrix = np.array([[1, 1, 1],
                         [n1[0], n2[0], n3[0]],
                         [n1[1], n2[1], n3[1]]])
    area = np.linalg.det(triMatrix)

    if (nodeEval == 1): # evaluating at node 1
        # find the derivative (gradient x and y) of the function
        grad_x = np.linalg.det([[0,1,1], [1, n2[0], n3[0]], [0, n2[1], n3[1]]]) / area
        grad_y = np.linalg.det([[0,1,1], [0, n2[0], n3[0]], [1, n2[1], n3[1]]]) / area
        gradvals = [grad_x, grad_y]

        for row in points: # find the scalar value at each point
            ptMatrix = np.array([[1, 1, 1],
                                [row[0], n2[0], n3[0]],
                                [row[1], n2[1], n3[1]]])
            val = abs(np.linalg.det(ptMatrix) / area)
            fvals.append(val)

    elif (nodeEval == 2): # evaluating at node 2
        # find gradient x and y
        grad_x = np.linalg.det([[1, 0, 1], [n1[0], 1, n3[0]], [n1[1], 0, n3[1]]]) / area
        grad_y = np.linalg.det([[1, 0, 1], [n1[0], 0, n3[0]], [n1[1], 1, n3[1]]]) / area
        gradvals = [grad_x, grad_y]
        
        for row in points: # find scalar values at points
            ptMatrix = np.array([[1, 1, 1],
                                [n1[0], row[0], n3[0]],
                                [n1[1], row[1], n3[1]]])
            val = abs(np.linalg.det(ptMatrix) / area)
            fvals.append(val)

    else: # evaluating at node 3
        #find gradient x and y
        grad_x = np.linalg.det([[1, 1, 0], [n1[0], n2[0], 1], [n1[1], n2[1], 0]]) / area
        grad_y = np.linalg.det([[1, 1, 0], [n1[0], n2[0], 0], [n1[1], n2[1], 1]]) / area
        gradvals = [grad_x, grad_y] 

        for row in points: # find scalar values at points 
            ptMatrix = np.array([[1, 1, 1],
                                [n1[0], n2[0], row[0]],
                                [n1[1], n2[1], row[1]]])
            val = abs(np.linalg.det(ptMatrix) / area)
            fvals.append(val)

    return fvals, gradvals
