from mesh import *
import matplotlib.pyplot as plt
import gauss_quadrature_tetra as q

nodes = [[0,0,0],
         [1,0,0],
         [.5,1,0],
         [0.5,0.5,-1]]

nodes = np.array(nodes)

loc = q.gauss_quadrature_tetra(nodes, 4)

thing = basis3Dscalar(nodes, 3, loc)

print("yay")