import ReadAbaqusMesh as abq
import quadrature as q

import numpy as np
import matplotlib.pyplot as plt
import sys


"""
initMesh initializes the mesh elements from abaqus file  
    Inputs:
    FileName, the name of the abaqus file
    quadPoints, the number of points to evaluate quadrature
    
    Outputs:
    elements, the list of elements
    boundaryList, list of nodes on the boundary
"""
def initMesh(FileName, quadPoints):
  nodeList, elementList, boundaryList = abq.readAbq(FileName)
  
  nodeList = np.array(nodeList)
  elementList = np.array(elementList)
  boundaryList = np.array(boundaryList) - 1 # -1 to take matlab offset into account

  nodeList[:,0] = nodeList[:, 0] - 1 # -1 to take matlab offset into account

  elementList = elementList[:, (1, 2, 3)] # Not taking permitivity into consideration
  elementList = elementList.astype(np.float64) # turning into floats from str
  elementList = elementList - 1 # -1 to take matlab offset into account
  
  elements = [] # Initializing Output element list
  
  for el in range(0,len(elementList)):
    coords = []
    nodeNums = []
    
    for i in elementList[el]:
      coords.append(nodeList[(nodeList[:,0] == i), (1,2,3)]) # Creating list with coordinates for each node
      nodeNums.append(nodeList[(nodeList[:,0] == i), 0]) # Creating list with node number for each node
      
    elements.append(Triangle(el, coords, nodeNums, quadPoints, boundaryList)) # Making triangle element object

  return elements, nodeList, boundaryList


'''
Triangle element for mesh
    Inputs:
      elNum, the number of the element in the mesh
      nodes, list of node coordinates of the element
      nodeNums, the list of each node number of global mesh
      points, the points per division of mesh
      boundaryList, list of nodes on the boundary
    
    Attributes:
      self.elNum, the mesh element number
      self.loc, discretized points array of element
      self.points, points per divsion
      self.nodes[i], contains attributes for each node
        self.nodes[i].loc, the coordinate location of the node
        self.nodes[i].globalNum, the node number in the global mesh
        self.nodes[i].num, the node number in non-boundary mesh 
        self.nodes[i].inBoundary, boolean for if in boundary
        self.nodes[i].N, the basis function for the node
        self.nodes[i].grad_N, the gradient of basis function of node
'''
class Triangle():
  def __init__(self, elNum, nodes: list, nodeNums, points, boundaryList):
    # Input Error Catching
    if(type is not np.ndarray):
      nodes = np.array(nodes)
    if(nodes.shape[0] != 3):
      print(f'Error: Wrong number of nodes for triangular element! Should be 3')
      sys.exit()
    if(nodes.shape != np.unique(nodes, axis=0).shape):
      print("Error: Nodes cannot be repeated")
      sys.exit()
    
    # Declarations
    self.elNum = elNum
    self.nodes = []
    
    # Initializing node properties
    for i in range(len(nodes)):
      self.nodes.append(type("Node", (), {})()) # Making a "node" object that will contain properties for each node
      
      self.nodes[i].loc = nodes[i] # Coordinate of Node
      self.nodes[i].globalNum = int(nodeNums[i]) # Global Node number
      self.nodes[i].inBoundary = (self.nodes[i].globalNum in boundaryList) # In Boundary Boolean
      self.nodes[i].num = -1 if self.nodes[i].inBoundary else self.nodes[i].globalNum # Reduced Node Number or -1 if in boundary meant for TM mode
      
    # Initializing
    self.points = points # These are the number of points for the quadrature

    self.loc = q.gauss_quadrature_trisimp_init(self.points, self.nodes[0].loc, self.nodes[1].loc, self.nodes[2].loc) # Initializing location of discretized points within element
    
    for i in range(len(nodes)):
      # Calculating and initializing the basis and gradient of the basis function for each node
      N, grad_N = nodal_int_eval(self.nodes[0].loc[0], self.nodes[0].loc[1], self.nodes[1].loc[0], self.nodes[1].loc[1], self.nodes[2].loc[0], self.nodes[2].loc[1], i, self.loc)
      self.nodes[i].N = np.array(N, dtype=np.complex_)
      self.nodes[i].grad_N = np.array(grad_N, dtype=np.complex_)

  '''Find if node number in element
  Inputs:
    nodeNum, node number to look for
    reduced = False, boolean if to look for in the global or reduced set of nodes
    
    Outputs:
    index, index of desired node
    -1 if not in element
  '''
  def findNodeIndex(self, nodeNum, reduced = False):
    index = 0
  
    for node in self.nodes:
      if((node.globalNum == nodeNum and not reduced) or (node.num == nodeNum and reduced)):
        return index
      index += 1
    
    return -1  


'''
THIS SHOULD BE RE-WRITTEN (SUPER SLOW)
reIndex finds the indices of the nodes not on the boundary
and writes them into [element].nodes[i].num
    Inputs:
    Elements, the list of mesh elements
    nodeList, the list of all nodes
    
    Outputs:
    elements, the element list with new indices 
'''
def reIndex(elements: list, nodeList: list):
  index = np.zeros(len(nodeList))
  indexShift = 0
  
  for i in range(len(nodeList)):
    for el in elements:
      for node in el.nodes:
        
        if(node.inBoundary):
          index[node.globalNum] = -1
        else:
          index[node.globalNum] = node.globalNum
  
  newIndex = []
  
  for val in index:
    if(val == -1):
      indexShift -= 1
      newIndex.append(-1)
    else:
      newIndex.append(val + indexShift)

  for i in range(len(nodeList)):
    for el in elements:
      for node in el.nodes:
        currIndex = el.findNodeIndex(i)
        if(currIndex == -1):
          continue
        
        if(newIndex[i] == -1):
          continue
        else:
          el.nodes[currIndex].num = int(newIndex[i])

  return elements


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
def nodal_int_eval(x1, y1, x2, y2, x3, y3, loc_node_num, loc):
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