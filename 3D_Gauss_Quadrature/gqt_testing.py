import numpy as np
import scipy as sp
from gauss_quadrature_tetra import *
import matplotlib.pyplot as plt

"""Testing File

Tests quadrature over a linear basis function using Georgia's basis3Dscalar function.

"""


# Define the kernel function
# this file tests the basis function for 3d scalar 
# uses test points within tetrahedral, outputs a 3d plot of function
#   values at the given points

# Evaluates 3D scalar basis functions and their derivative over a tetrahedron
def basis3Dscalar(n1, n2, n3, n4, points, node_eval):
    # inputs:
    # n1, n2, n3, and n4 are the points of the tetrahedron
    # points is a list of points to evaluate at 
    # node_eval is the node (1, 2, 3, or 4) to evaluate from

    # outputs:
    # fvals is a list of the function value at the given points
    # gradvals is the gradient of the function in form [x,y,z]
    fvals = []
    gradvals = []

    import numpy as np

    # calculate determinant of tetrahedron 
    tetra = np.array([[n1[0], n1[1], n1[2], 1],
             [n2[0], n2[1], n2[2], 1],
             [n3[0], n3[1], n3[2], 1],
             [n4[0], n4[1], n4[2], 1]])
    # v = det(matrix)
    v = abs(np.linalg.det(tetra))

    if (node_eval == 1): 
        # make gradient matrices 
        grad_x_matrix = tetra
        grad_x_matrix[0] = [1,0,0,0]
        grad_y_matrix = tetra
        grad_y_matrix[0] = [0,1,0,0]
        grad_z_matrix = tetra
        grad_z_matrix[0] = [0,0,1,0]

        # evaluate gradient values
        grad_x = np.linalg.det(grad_x_matrix) / v
        grad_y = np.linalg.det(grad_y_matrix) / v
        grad_z = np.linalg.det(grad_z_matrix) / v
        gradvals = [grad_x, grad_y, grad_z]

        # evaluate function at points
        for row in points:
            tetra1 = np.array([[1, row[0], row[1], row[2]],
                    [1, n2[0], n2[1], n2[2]],
                    [1, n3[0], n3[1], n3[2]],
                    [1, n4[0], n4[1], n4[2]]]) 
            # val = det(matrix) / (6 * v)
            val = abs(np.linalg.det(tetra1)) / v
            fvals.append(val)

    elif (node_eval == 2):
        # make gradient matrices 
        grad_x_matrix = tetra
        grad_x_matrix[1] = [1,0,0,0]
        grad_y_matrix = tetra
        grad_y_matrix[1] = [0,1,0,0]
        grad_z_matrix = tetra
        grad_z_matrix[1] = [0,0,1,0]

        # evaluate gradient values
        grad_x = np.linalg.det(grad_x_matrix) / v
        grad_y = np.linalg.det(grad_y_matrix) / v
        grad_z = np.linalg.det(grad_z_matrix) / v
        gradvals = [grad_x, grad_y, grad_z]

        # evaluate function at points
        for row in points:
            tetra2 = np.array([[1, n1[0], n1[1], n1[2]],
                    [1, row[0], row[1], row[2]],
                    [1, n3[0], n3[1], n3[2]],
                    [1, n4[0], n4[1], n4[2]]]) 
            val = abs(np.linalg.det(tetra2)) / v
            fvals.append(val)

    elif (node_eval == 3):
        # make gradient matrices 
        grad_x_matrix = tetra
        grad_x_matrix[2] = [1,0,0,0]
        grad_y_matrix = tetra
        grad_y_matrix[2] = [0,1,0,0]
        grad_z_matrix = tetra
        grad_z_matrix[2] = [0,0,1,0]

        # evaluate gradient values
        grad_x = np.linalg.det(grad_x_matrix) / v
        grad_y = np.linalg.det(grad_y_matrix) / v
        grad_z = np.linalg.det(grad_z_matrix) / v
        gradvals = [grad_x, grad_y, grad_z]

        # evaluate function at points
        for row in points:
            tetra3 = np.array([[1, n1[0], n1[1], n1[2]],
                    [1, n2[0], n2[1], n2[2]],
                    [1, row[0], row[1], row[2]],
                    [1, n4[0], n4[1], n4[2]]]) 
            val = abs(np.linalg.det(tetra3)) / v
            fvals.append(val)
    
    else: # node_eval == 4
        # make gradient matrices 
        grad_x_matrix = tetra
        grad_x_matrix[3] = [1,0,0,0]
        grad_y_matrix = tetra
        grad_y_matrix[3] = [0,1,0,0]
        grad_z_matrix = tetra
        grad_z_matrix[3] = [0,0,1,0]

        # evaluate gradient values
        grad_x = np.linalg.det(grad_x_matrix) / v
        grad_y = np.linalg.det(grad_y_matrix) / v
        grad_z = np.linalg.det(grad_z_matrix) / v
        gradvals = [grad_x, grad_y, grad_z]

        # evaluate function at points
        for row in points:
            tetra4 = np.array([[1, n1[0], n1[1], n1[2]],
                    [1, n2[0], n2[1], n2[2]],
                    [1, n3[0], n3[1], n3[2]],
                    [1, row[0], row[1], row[2]]]) 
            val = abs(np.linalg.det(tetra4)) / v
            fvals.append(val)

    return fvals, gradvals

# Main function
def main():
    nodes_tetra1 = np.array([[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]])
    
    #returns a list of points inside tetrahedron
    loc = gauss_quadrature_tetra(nodes_tetra1, 4)

    #Plots our tetrahedron
    fig = plt.figure()
    ax = fig.add_subplot(projection="3d")
    ax.scatter(loc[:, 0], loc[:, 1], loc[:, 2])

    #Returns list of function values
    flist, grad= basis3Dscalar([0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], loc, 1)

    #Performs quadrature integration
    val = gauss_quadrature_tetra(nodes_tetra1, 4, flist)
    print(val)
    plt.show()



if __name__ == "__main__":
    main()
