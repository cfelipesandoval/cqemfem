import numpy as np
from basis2Dscalar import basis2Dscalar as b2s
from basis2Dvector import basis2Dvector as b2v

def basis2DSquad(nodes, points):
    # nodes is a list of tri vertices
    # points is a list of points to evaluate
    
    [n1, n2, n3] = nodes                                        # readability
    fvals = []                                                  # list of function values for given points: [N1, N2, N3, ...] where Nn is interpolation for each point
    gvalsx = []                                                 # list of function x-gradients: [gN1, gN2, gN3, ...] where gNn is x-gradient associated with Nn
    gvalsy = []                                                 # list of function y-gradients: [gN1, gN2, gN3, ...] where gNn is y-gradient associated with Nn
    aMat = np.array([[1,     1,     1    ],
                     [n1[0], n2[0], n3[0]],
                     [n1[1], n2[1], n3[1]]])
    area = np.linalg.det(aMat)                             # area = 2A
    # a = [n1[0] * n2[1] - n2[0] * n1[1],
    #      n3[0] * n1[1] - n1[0] * n3[1],
    #      n1[0] * n3[1] - n3[0] * n1[1]]
    a = [n2[0] * n3[1] - n3[0] * n2[1],
         n3[0] * n1[1] - n1[0] * n3[1],
         n1[0] * n2[1] - n2[0] * n1[1]]
    b = [n2[1] - n3[1], n3[1] - n1[1], n1[1] - n2[1]]
    c = [n3[0] - n2[0], n1[0] - n3[0], n2[0] - n1[0]]

    for pt in points:
        fval = []
        gvx = []
        gvy = []
        lvals = []

        # EVALUATE N1, N2, N3
        for i in range(3):
            lval, _ = b2s(n1, n2, n3, i + 1, [pt])
            lval = lval.pop()
            fval.append(lval * (2 * lval - 1))                  # Nn = Ln(2Ln - 1)
            gvx.append(4 * b[i] / area * lval - b[i] / area)    # gxNn = 2b/A * Ln - b/2A: b is respective node coefficient
            gvy.append(4 * c[i] / area * lval - c[i] / area)    # gyNn = 2c/A * Ln - c/2A: c is respective node coefficient
            lvals.append(lval)

        # EVALUATE N4 = 4L1L2, N5 = 4L2L3, N6 = 4L1L3
        # gxNn = b2(a1 + c1y) + b1(a2 + 2b2x + c2y) / A^2
        # gyNn = c2(a1 + b1x) + c1(a2 + b2x + 2c2y) / A^2
        for i in range(3):
            fval.append(4 * lvals[i] * lvals[i - 2])
            gvx.append(4*(b[i-2] * (a[i] + c[i] * pt[1]) + b[i] * (a[i-2] + 2*b[i-2] * pt[0] + c[i-2] * pt[1])) / area ** 2)
            gvy.append(4*(c[i-2] * (a[i] + b[i] * pt[0]) + c[i] * (a[i-2] + b[i-2] * pt[0] + 2*c[i-2] * pt[1])) / area ** 2)
        
        fvals.append(fval)
        gvalsx.append(gvx)
        gvalsy.append(gvy)

    return np.array(fvals), np.array(gvalsx), np.array(gvalsy)


def basis2DSqbic(nodes, points):
    # nodes is a list of tri vertices
    # points is a list of points to evaluate
    
    [n1, n2, n3] = nodes                                        # readability
    fvals = []                                                  # list of function values for given points: [N1, N2, N3, ...] where Nn is interpolation for each point
    gvalsx = []                                                 # list of function x-gradients: [gN1, gN2, gN3, ...] where gNn is x-gradient associated with Nn
    gvalsy = []                                                 # list of function y-gradients: [gN1, gN2, gN3, ...] where gNn is y-gradient associated with Nn
    aMat = np.array([[1,     1,     1    ],
                     [n1[0], n2[0], n3[0]],
                     [n1[1], n2[1], n3[1]]])
    area = np.linalg.det(aMat)                                  # area = 2A
    a = [n2[0] * n3[1] - n3[0] * n2[1],
         n3[0] * n1[1] - n1[0] * n3[1],
         n1[0] * n2[1] - n2[0] * n1[1]]
    b = [n2[1] - n3[1], n3[1] - n1[1], n1[1] - n2[1]]
    c = [n3[0] - n2[0], n1[0] - n3[0], n2[0] - n1[0]]

    for pt in points:
        fval = []
        gvx = []
        gvy = []
        lvals = []

        # EVALUATE N1, N2, N3
        for i in range(3):
            lval, _ = b2s(n1, n2, n3, i + 1, [pt])
            lval = lval.pop()
            fval.append(.5 * lval * (3 * lval - 1) * (3 * lval - 2))                # Nn = 0.5Ln(3Ln - 1)(3Ln - 2), Ln = f(x,y), hence chain rule:
            gvx.append(b[i]/area * (27/2*lval**2 - 9*lval + 1))                     # gxNn = b/area * (27/2Ln^2 - 9Ln + 1)
            gvy.append(c[i]/area * (27/2*lval**2 - 9*lval + 1))                     # gyNn = c/area * (27/2Ln^2 - 9Ln + 1)
            lvals.append(lval)

        # EVALUATE N4-N9
        for i in range(3):
            fval.append(9/2 * lvals[i] * lvals[i - 2] * (3*lvals[i] - 1))            # each edge uses vertices' corresponding Lm and Ln to make
            fval.append(9/2 * lvals[i] * lvals[i - 2] * (3*lvals[i - 2] - 1))        # Nn = 4.5Lm*Ln*(3Lm/n - 1), chain and product rules are not trivial:
            # 27/2*Lm'*Lm*Ln + 9/2*Lm'(3Lm-1)3Ln + 9/2Ln'Lm(3Lm-1)
            # where Lm appears in the (3Lm-1) term and Lm'/Ln' = b/area (d/dx) or c/area (d/dy)
            gvx.append((27/2*b[i]*lvals[i]*lvals[i-2] + 9/2*b[i]*(3*lvals[i]-1)*lvals[i-2] + 9/2*b[i-2]*lvals[i]*(3*lvals[i]-1))/area)
            gvx.append((27/2*b[i-2]*lvals[i]*lvals[i-2] + 9/2*b[i-2]*(3*lvals[i-2]-1)*lvals[i] + 9/2*b[i]*lvals[i-2]*(3*lvals[i-2]-1))/area)
            gvy.append((27/2*c[i]*lvals[i]*lvals[i-2] + 9/2*c[i]*(3*lvals[i]-1)*lvals[i-2] + 9/2*c[i-2]*lvals[i]*(3*lvals[i]-1))/area)
            gvy.append((27/2*c[i-2]*lvals[i]*lvals[i-2] + 9/2*c[i-2]*(3*lvals[i-2]-1)*lvals[i] + 9/2*c[i]*lvals[i-2]*(3*lvals[i-2]-1))/area)
        
        # EVALUATE N10
        fval.append(27 * np.prod(lvals))                                                # 27 * L1 * L2 * L3, chain/product rules gives:
        gvx.append(27/(area**3) * sum([b[i]*lvals[i-1]*lvals[i-2] for i in range(3)]))  # 27(L1'L2L3 + L1L2'L3 + L1L2L3')
        gvy.append(27/(area**3) * sum([c[i]*lvals[i-1]*lvals[i-2] for i in range(3)]))  # where Ln' = b/area (d/dx) or c/area (d/dy)
        
        fvals.append(fval)
        gvalsx.append(gvx)
        gvalsy.append(gvy)

    return np.array(fvals), np.array(gvalsx), np.array(gvalsy)


def basis2DVquad(nodes, points):
    # nodes is a list of tri vertices
    # points is a list of points to evaluate

    [n1, n2, n3] = nodes                                        # readability
    fvals = []

    for pt in points:
        fval = []

        W1 = np.array(b2v(n1, n2, n3, 1, 2, [pt]))              # Wn here are actually ln * Wn
        W2 = np.array(b2v(n1, n2, n3, 2, 3, [pt]))
        W3 = np.array(b2v(n1, n2, n3, 3, 1, [pt]))
        N1, _ = b2s(n1, n2, n3, 1, [pt])
        N2, _ = b2s(n1, n2, n3, 2, [pt])
        N3, _ = b2s(n1, n2, n3, 3, [pt])
        
        # EVALUATE EDGE BASES
        fval.append(W1[0] * (3*N1[0] - 1))
        fval.append(W1[0] * (3*N2[0] - 1))
        fval.append(W2[0] * (3*N2[0] - 1))
        fval.append(W2[0] * (3*N3[0] - 1))
        fval.append(W3[0] * (3*N3[0] - 1))
        fval.append(W3[0] * (3*N1[0] - 1))

        # EVALUATE INTERIOR BASES: one can be discarded
        fval.append(3/2 * W1[0] * 3*N3[0])
        fval.append(3/2 * W2[0] * 3*N1[0])
        fval.append(3/2 * W3[0] * 3*N2[0])

        fvals.append(fval)

    return fvals


def basis2DVqbic(nodes, points):
    # nodes is a list of tri vertices
    # points is a list of points to evaluate

    [n1, n2, n3] = nodes                                        # readability
    fvals = []

    for pt in points:
        fval = []

        W1 = np.array(b2v(n1, n2, n3, 1, 2, [pt]))              # Wn here are actually ln * Wn
        W2 = np.array(b2v(n1, n2, n3, 2, 3, [pt]))
        W3 = np.array(b2v(n1, n2, n3, 3, 1, [pt]))
        N1, _ = b2s(n1, n2, n3, 1, [pt])
        N2, _ = b2s(n1, n2, n3, 2, [pt])
        N3, _ = b2s(n1, n2, n3, 3, [pt])
        
        # EVALUATE EDGE BASES
        fval.append(W1[0] * (4*N1[0] - 1) * (4*N1[0] - 2))
        fval.append(W1[0] * (3*N1[0] - 1) * (4*N2[0] - 1))
        fval.append(W1[0] * (3*N2[0] - 1) * (4*N2[0] - 2))
        fval.append(W2[0] * (3*N2[0] - 1) * (4*N2[0] - 2))
        fval.append(W2[0] * (3*N2[0] - 1) * (4*N3[0] - 1))
        fval.append(W2[0] * (3*N3[0] - 1) * (4*N3[0] - 2))
        fval.append(W3[0] * (3*N3[0] - 1) * (4*N3[0] - 2))
        fval.append(W3[0] * (3*N3[0] - 1) * (4*N1[0] - 1))
        fval.append(W3[0] * (3*N1[0] - 1) * (4*N1[0] - 2))

        # EVALUATE INTERIOR BASES: one of each three can be discarded
        fval.append(4/3 * W1[0] * (4*N1[0] - 1) * 4*N3[0])
        fval.append(  2 * W2[0] * 4*N1[0] * (4*N1[0] - 1))
        fval.append(4/3 * W3[0] * (4*N1[0] - 1) * 4*N2[0])
        fval.append(4/3 * W1[0] * (4*N2[0] - 1) * 4*N3[0])
        fval.append(4/3 * W2[0] * 4*N1[0] * (4*N2[0] - 1))
        fval.append(  2 * W3[0] * 4*N2[0] * (4*N2[0] - 1))
        fval.append(  2 * W1[0] * 4*N3[0] * (4*N3[0] - 1))
        fval.append(4/3 * W2[0] * 4*N1[0] * (4*N3[0] - 1))
        fval.append(4/3 * W3[0] * 4*N2[0] * (4*N3[0] - 1))

        fvals.append(fval)

    return fvals


def basis2Dlinr(x1, y1, x2, y2, x3, y3, loc_node_num, loc):
    """
    nodal_int_eval calculates the values of the linear interpolant and the
    gradient of the linear interpolant for a triangle  
    %   Inputs:
    %   xi, the x value at the ith node
    %   yi, the y value at the ith node
        x1, y1, x2, y2, x3, y3: Coordinates of the three nodes.
    %   loc_node_num, the local node number that the evaluation is being done
    %   for
    %   loc, the [x,y] values that the evaluation is to take place at
    %
    %   Outputs:
    %   N, the evaluated linear interpolant
    %   grad_N, the evaluated gradient of the linear interpolant, formated as
    %   a 2D array of x_hat and y_hat values in columns and different
    %   evaluation locations as the rows
    """
    # Calculate the length of the first column of loc
    length = len(loc[:, 0])
    # Array w/ same number of rows as loc first column and has one column
    N = np.zeros((length, 1))

    # Array w/ same number of rows as loc first column and has two columns
    grad_N = np.zeros((length, 2))

    # Create area_mat as a 3x3 NumPy array
    area_mat = np.array([[1, 1, 1], [x1, x2, x3], [y1, y2, y3]])

    # Calculate the absolute determinant of area_mat
    area = abs(np.linalg.det(area_mat))

    if loc_node_num == 0:
        # Create area_mat as a 3x3 NumPy array
        Nx_mat = np.array([[0, 1, 1], 
                           [1, x2, x3], 
                           [0, y2, y3]])

        # Calculate grad_Nx
        grad_Nx = np.linalg.det(Nx_mat) / area

        # Create another area_mat for grad_Ny
        Ny_mat = np.array([[0, 1, 1], 
                           [0, x2, x3], 
                           [1, y2, y3]])

        # Calculate grad_Ny
        grad_Ny = np.linalg.det(Ny_mat) / area

        for m in range(len(loc[:, 1])):
            loc_mat = np.array([[1, 1, 1], [loc[m, 0], x2, x3], [loc[m, 1], y2, y3]])
            N[m] = np.linalg.det(loc_mat) / area
            grad_N[m, 0] = grad_Nx
            grad_N[m, 1] = grad_Ny

    elif loc_node_num == 1:
        # Create area_mat as a 3x3 NumPy array
        Nx_mat = np.array([[1, 0, 1], [x1, 1, x3], [y1, 0, y3]])

        grad_Nx = np.linalg.det(Nx_mat) / area

        # Create another area_mat for grad_Ny
        Ny_mat = np.array([[1, 0, 1], [x1, 0, x3], [y1, 1, y3]])

        grad_Ny = np.linalg.det(Ny_mat) / area

        for m in range(len(loc[:, 1])):
            loc_mat = np.array([[1, 1, 1], [x1, loc[m, 0], x3], [y1, loc[m, 1], y3]])
            N[m] = np.linalg.det(loc_mat) / area
            grad_N[m, 0] = grad_Nx
            grad_N[m, 1] = grad_Ny

    elif loc_node_num == 2:
        # Create area_mat as a 3x3 NumPy array
        Nx_mat = np.array([[1, 1, 0], [x1, x2, 1], [y1, y2, 0]])

        grad_Nx = np.linalg.det(Nx_mat) / area

        # Create another area_mat for grad_Ny
        Ny_mat = np.array([[1, 1, 0], [x1, x2, 0], [y1, y2, 1]])

        grad_Ny = np.linalg.det(Ny_mat) / area

        for m in range(len(loc[:, 1])):
            loc_mat = np.array([[1, 1, 1], [x1, x2, loc[m, 0]], [y1, y2, loc[m, 1]]])
            N[m] = np.linalg.det(loc_mat) / area
            grad_N[m, 0] = grad_Nx
            grad_N[m, 1] = grad_Ny

    else:
        # error statement
        print(
            f"This function only works for a 2D triangle, loc_node_num should be ranging from 1 to 3"
        )
        return

    return N, grad_N


if __name__ == '__main__':
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.tri as tri
    import matplotlib.animation as anim

    nodes = np.array([[-.5, 1/6], [.5, 0], [0, 0.5]])
    # ff is a real Cubit trimesh of the above tri
    pts = np.array([[0.5,0],
                    [-0.5,0.1666667],
                    [0,0.5],
                    [0.4565217,0.007246377],
                    [0.4130435,0.01449275],
                    [0.3695652,0.02173913],
                    [0.326087,0.02898551],
                    [0.2826087,0.03623188],
                    [0.2391304,0.04347826],
                    [0.1956522,0.05072464],
                    [0.1521739,0.05797101],
                    [0.1086957,0.06521739],
                    [0.06521739,0.07246377],
                    [0.02173913,0.07971014],
                    [-0.02173913,0.08695652],
                    [-0.06521739,0.0942029],
                    [-0.1086957,0.1014493],
                    [-0.1521739,0.1086957],
                    [-0.1956522,0.115942],
                    [-0.2391304,0.1231884],
                    [-0.2826087,0.1304348],
                    [-0.326087,0.1376812],
                    [-0.3695652,0.1449275],
                    [-0.4130435,0.1521739],
                    [-0.4565217,0.1594203],
                    [-0.4642857,0.1904762],
                    [-0.4285714,0.2142857],
                    [-0.3928571,0.2380952],
                    [-0.3571429,0.2619048],
                    [-0.3214286,0.2857143],
                    [-0.2857143,0.3095238],
                    [-0.25,0.3333333],
                    [-0.2142857,0.3571429],
                    [-0.1785714,0.3809524],
                    [-0.1428571,0.4047619],
                    [-0.1071429,0.4285714],
                    [-0.07142857,0.452381],
                    [-0.03571429,0.4761905],
                    [0.03125,0.46875],
                    [0.0625,0.4375],
                    [0.09375,0.40625],
                    [0.125,0.375],
                    [0.15625,0.34375],
                    [0.1875,0.3125],
                    [0.21875,0.28125],
                    [0.25,0.25],
                    [0.28125,0.21875],
                    [0.3125,0.1875],
                    [0.34375,0.15625],
                    [0.375,0.125],
                    [0.40625,0.09375],
                    [0.4375,0.0625],
                    [0.46875,0.03125],
                    [0.2690909,0.173912],
                    [0.2073117,0.2385617],
                    [-0.001311541,0.4575276],
                    [-0.03048676,0.431748],
                    [-0.06420741,0.4079583],
                    [-0.1039348,0.3850688],
                    [-0.139649,0.3612593],
                    [-0.1753633,0.3374497],
                    [-0.2467919,0.2898307],
                    [-0.3177308,0.2423115],
                    [-0.3539348,0.2184021],
                    [-0.3881371,0.1912344],
                    [-0.2547006,0.1647879],
                    [-0.1676375,0.1499721],
                    [-0.3418752,0.178879],
                    [-0.2980723,0.1717112],
                    [-0.2818735,0.2630129],
                    [0.1801886,0.0920011],
                    [0.2668981,0.07715098],
                    [-0.2111158,0.1572185],
                    [0.3106234,0.07026197],
                    [-0.2110776,0.3136402],
                    [0.3574139,0.07060067],
                    [-0.1241592,0.1427257],
                    [0.3970088,0.05493521],
                    [-0.08068098,0.1354794],
                    [-0.03720271,0.128233],
                    [0.006275546,0.1209866],
                    [0.04975381,0.1137402],
                    [0.09323207,0.1064939],
                    [0.01933202,0.425356],
                    [0.04688581,0.3945881],
                    [0.0811416,0.363692],
                    [0.1135617,0.3323117],
                    [0.1448117,0.3010617],
                    [0.1760617,0.2698117],
                    [0.1367103,0.09924748],
                    [0.2385617,0.2073117],
                    [0.2236669,0.08475473],
                    [0.3001809,0.1447728],
                    [0.3326117,0.1118796],
                    [0.2902008,0.1078025],
                    [-0.3118428,0.2080264],
                    [-0.01028763,0.3908268],
                    [0.1333734,0.2583734],
                    [0.1021234,0.2896234],
                    [0.06960087,0.3213596],
                    [0.03250679,0.3546508],
                    [0.1219088,0.1414474],
                    [0.07776848,0.1477703],
                    [0.03429022,0.1550167],
                    [-0.009188038,0.1622631],
                    [-0.0526663,0.1695095],
                    [-0.09614456,0.1767558],
                    [0.2530253,0.1242132],
                    [0.1639066,0.1332624],
                    [0.2082033,0.1260312],
                    [-0.1831011,0.1912486],
                    [-0.1396228,0.1840022],
                    [-0.272465,0.2113576],
                    [-0.227348,0.1993218],
                    [-0.2080795,0.2709526],
                    [-0.2433367,0.2448002],
                    [-0.1721552,0.2939471],
                    [-0.1364409,0.3177566],
                    [-0.1007266,0.3415662],
                    [-0.06003955,0.3668251],
                    [0.1961597,0.1967674],
                    [0.1646234,0.2271234],
                    [0.2273542,0.1643207],
                    [-0.2010536,0.2312957],
                    [0.185432,0.1614951],
                    [0.1511417,0.1820135],
                    [-0.0166039,0.3479352],
                    [-0.06040536,0.323415],
                    [-0.09751854,0.2980635],
                    [-0.1321051,0.2693605],
                    [-0.1699543,0.2550455],
                    [-0.154291,0.2247578],
                    [-0.1125741,0.2202847],
                    [-0.06866464,0.2116377],
                    [-0.02465162,0.2035395],
                    [0.01901658,0.1965497],
                    [0.06069192,0.1878161],
                    [0.105668,0.1810949],
                    [0.02329089,0.3111772],
                    [0.05943512,0.2781851],
                    [0.1225704,0.2166644],
                    [0.09122206,0.248917],
                    [-0.08914947,0.2544678],
                    [0.08107055,0.2149739],
                    [0.04727559,0.2350265],
                    [0.02127777,0.2688635],
                    [0.002403749,0.2364429],
                    [-0.04113978,0.2442434],
                    [-0.05756381,0.2806459],
                    [-0.02086516,0.3101452],
                    [-0.01427875,0.2775174]])
    # test = False # T if vector, F if scalar
    # dim = False # T if cubic, F if quadratic

    for test in [False]:
        for dim in [True]:
            # TESTING QUADRATIC SCALAR
            if not test and not dim:
                fvals, gvx, gvy = basis2DSquad(nodes, pts)
                def rotator(angle): # from MatPlotLib docs
                    # Normalize the angle to the range [-180, 180] for display
                    azim = (angle + 180) % 360 - 180
                    ax.view_init(None, azim)
                for i in range(6):
                    # PLOTTING BASIS
                    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
                    ax.plot_trisurf(pts[:, 0], pts[:, 1], fvals[:, i], cmap=cm.inferno)
                    ax.set_xlabel("x")
                    ax.set_xlim(-.5, .5)
                    ax.set_ylabel("y")
                    ax.set_ylim(-.5, .5)
                    ax.set_zlabel("interpolation")
                    plot = anim.FuncAnimation(fig, rotator, range(0, 361, 5))
                    plot.save(f'graphs/2DQuadScal/BasisN{i + 1}.gif', anim.PillowWriter())
                    # PLOTTING x-GRADIENT
                    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
                    ax.plot_trisurf(pts[:, 0], pts[:, 1], gvx[:, i], cmap=cm.inferno)
                    ax.set_xlabel("x")
                    ax.set_xlim(-.5, .5)
                    ax.set_ylabel("y")
                    ax.set_ylim(-.5, .5)
                    ax.set_zlabel("interpolation")
                    plot = anim.FuncAnimation(fig, rotator, range(0, 361, 5))
                    plot.save(f'graphs/2DQuadScal/GradxN{i + 1}.gif', anim.PillowWriter())
                    # PLOTTING y-GRADIENT
                    fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
                    ax.plot_trisurf(pts[:, 0], pts[:, 1], gvy[:, i], cmap=cm.inferno)
                    ax.set_xlabel("x")
                    ax.set_xlim(-.5, .5)
                    ax.set_ylabel("y")
                    ax.set_ylim(-.5, .5)
                    ax.set_zlabel("interpolation")
                    plot = anim.FuncAnimation(fig, rotator, range(0, 361, 5))
                    plot.save(f'graphs/2DQuadScal/GradyN{i + 1}.gif', anim.PillowWriter())
            
            # TESTING QUADRATIC VECTOR
            if test and not dim:
                fvals = np.array(basis2DVquad(nodes, pts))
                # print(fvals)

                bases = ['E1-N210', 'E1-N120', 'E2-N021', 'E2-N012', 'E3-N102', 'E3-N201', 'E1-N111', 'E2-N111', 'E3-N111']

                for vec in range(len(fvals[0])):
                    fig, ax = plt.subplots()
                    ax.scatter(nodes[:, 0], nodes[:, 1])
                    ax.set_xlabel("x")
                    ax.set_ylabel("y")
                    for i, txt in enumerate(['n1', 'n2', 'n3']):
                        ax.text(nodes[i, 0], nodes[i, 1], txt)
                    for idx, pt in enumerate(pts):
                        x = pt[0]
                        y = pt[1]
                        q = plt.quiver(pt[0], pt[1], fvals[idx, vec, 0], fvals[idx, vec, 1], width = .001, scale = 40)
                        plt.title(bases[vec])
                    plt.savefig(f'graphs/2DQuadVec/{bases[vec]}.png', dpi=300)
                    plt.close()

            # TESTING CUBIC SCALAR
            if not test and dim:
                fvals, gvx, gvy = basis2DSqbic(nodes, pts)
                print(fvals)
                print("\n\n")
                print(gvx)
                print("\n\n")
                print(gvy)
                # def rotator(angle): # from MatPlotLib docs
                #     # Normalize the angle to the range [-180, 180] for display
                #     azim = (angle + 180) % 360 - 180
                #     ax.view_init(None, azim)
                # for i in range(10):
                #     # PLOTTING BASIS
                #     fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
                #     ax.plot_trisurf(pts[:, 0], pts[:, 1], fvals[:, i], cmap=cm.inferno)
                #     ax.set_xlabel("x")
                #     ax.set_xlim(-.5, .5)
                #     ax.set_ylabel("y")
                #     ax.set_ylim(-.5, .5)
                #     ax.set_zlabel("interpolation")
                #     plot = anim.FuncAnimation(fig, rotator, range(0, 361, 5))
                #     plot.save(f'graphs/2DQbicScal/BasisN{i + 1}.gif', anim.PillowWriter())
                #     plt.close()
                #     # PLOTTING x-GRADIENT
                #     fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
                #     ax.plot_trisurf(pts[:, 0], pts[:, 1], gvx[:, i], cmap=cm.inferno)
                #     ax.set_xlabel("x")
                #     ax.set_xlim(-.5, .5)
                #     ax.set_ylabel("y")
                #     ax.set_ylim(-.5, .5)
                #     ax.set_zlabel("interpolation")
                #     plot = anim.FuncAnimation(fig, rotator, range(0, 361, 5))
                #     plot.save(f'graphs/2DQbicScal/GradxN{i + 1}.gif', anim.PillowWriter())
                #     plt.close()
                #     # PLOTTING y-GRADIENT
                #     fig, ax = plt.subplots(subplot_kw={'projection': '3d'})
                #     ax.plot_trisurf(pts[:, 0], pts[:, 1], gvy[:, i], cmap=cm.inferno)
                #     ax.set_xlabel("x")
                #     ax.set_xlim(-.5, .5)
                #     ax.set_ylabel("y")
                #     ax.set_ylim(-.5, .5)
                #     ax.set_zlabel("interpolation")
                #     plot = anim.FuncAnimation(fig, rotator, range(0, 361, 5))
                #     plot.save(f'graphs/2DQbicScal/GradyN{i + 1}.gif', anim.PillowWriter())
                #     plt.close()

            # TESTING CUBIC VECTOR
            if test and dim:
                fvals = np.array(basis2DVqbic(nodes, pts))
                # print(fvals)

                bases = ['E1-N310', 'E1-N220', 'E1-N130', 'E2-N031', 'E2-N022', 'E2-N013', 'E3-N103', 'E3-N202', 'E3-N301',
                        'E1-N211', 'E2-N211', 'E3-N211', 'E1-N121', 'E2-N121', 'E3-N121', 'E1-N112', 'E2-N112', 'E3-N112']

                for vec in range(len(fvals[0])):
                    fig, ax = plt.subplots()
                    ax.scatter(nodes[:, 0], nodes[:, 1])
                    ax.set_xlabel("x")
                    ax.set_ylabel("y")
                    for i, txt in enumerate(['n1', 'n2', 'n3']):
                        ax.text(nodes[i, 0], nodes[i, 1], txt)
                    for idx, pt in enumerate(pts):
                        x = pt[0]
                        y = pt[1]
                        q = plt.quiver(pt[0], pt[1], fvals[idx, vec, 0], fvals[idx, vec, 1], width = .001, scale = 40)
                        plt.title(bases[vec])
                    plt.savefig(f'graphs/2DQbicVec/{bases[vec]}.png', dpi=300)
                    plt.close()