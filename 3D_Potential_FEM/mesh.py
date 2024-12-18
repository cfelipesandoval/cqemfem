import ReadAbaqusMesh as abq
import gauss_quadrature_tetra as q

import numpy as np
import matplotlib.pyplot as plt
import sys
import time


ep0 = 1 * 10 ** -9 / (36 * np.pi)
mu0 = 4 * np.pi * 10**-7


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
  
  nodeList, elementList, boundaryList, ports = abq.readAbq_ports(FileName)

  nodeList = np.array(nodeList)
  elementList = np.array(elementList)
  elementList = elementList[:,1:8]
  elementList[:,0:4] = np.array(elementList[:,0:4]) - 1
  boundaryList = np.array(boundaryList) - 1 # -1 to take matlab offset into account
  portList = np.array(ports) - 1

  nodeList[:,0] = nodeList[:,0] - 1 # -1 to take matlab offset into account
 
  elements = [] # Initializing Output element list
  
  nodesOnBoundary = np.unique(boundaryList)
  
  index = np.zeros(len(nodeList))
  indexShift = 0
  
  for i in range(len(nodeList)):
    if(nodeList[i,0] in nodesOnBoundary):
      index[i] = -1
    else:
      index[i] = int(nodeList[i,0])
  
  newIndex = []
  
  for val in index:
    if(val == -1):
      indexShift -= 1
      newIndex.append(-1)
    else:
      newIndex.append(val + indexShift)

  for el in range(0, len(elementList)):
    coords = []
    nodeNums = []
    num = []
    
    for i in elementList[el, 0:4].astype(int):
      coords.append(nodeList[(nodeList[:,0] == i), (1,2,3)]) # Creating list with coordinates for each node
      nodeNums.append(i)
      num.append(newIndex[i])
      
    coords = np.array(coords)
    
    elements.append(Tetrahedron(el, coords, nodeNums, quadPoints, num, elementList[el,4], elementList[el,5])) # Making triangle element object
  
  numberOfEdges = len(Tetrahedron.edgeList)
  
  boundaryEdgeNums = edgesOnBoudnary(boundaryList)
  
  Tetrahedron.nodeList = np.array(Tetrahedron.nodeList)
  
  elements = reIndex(elements, boundaryEdgeNums)
  
  for el in elements:
    el.getBasisVectors()
    
  unknownNodes = int(len(nodeList) - len(np.unique(boundaryList)))
  unknownEdges = numberOfEdges - len(boundaryEdgeNums)

  return elements, nodeList, boundaryList, unknownNodes, unknownEdges

class Tetrahedron:
  edgeList = np.empty((0,2))
  nodeList = []
  
  def __init__(self, elNum, nodes: list, nodeNums, points, num, epr = 1, mur = 1):
    # Input Error Catching
    if(type is not np.ndarray):
      nodes = np.array(nodes)
    if(nodes.shape[0] != 4):
      print(f'Error: Wrong number of nodes for triangular element! Should be 4')
      sys.exit()
    if(nodes.shape[0] != len(nodeNums)):
      print(f'Error: Only {len(nodeNums)} global node numbers were given, Should be 4!')
      sys.exit()
    if(nodes.shape != np.unique(nodes, axis=0).shape):
      print("Error: Nodes cannot be repeated")
      sys.exit()
    
    # Declarations
    self.num = elNum
    
    self.ep = epr * ep0
    self.mu = mur * mu0
    
    self.edges = []
    self.nodes = []
    
    # Initializing
    self.gaussPoints = points # These are the number of points for the quadrature

    self.loc = q.gauss_quadrature_tetra(nodes, self.gaussPoints) # Initializing location of discretized points within element
    self.loc = np.array(self.loc)
    
    nodeNums = np.array(nodeNums)
    num = np.array(num)
    inds = nodeNums.argsort()
    nodes = nodes[inds,:]
    nodeNums = nodeNums[inds]
    num = num[inds]
    
    self.nodeList = np.copy(nodes)
    
    # Node Elements
    for node in range(0,4):
      self.nodes.append(type("Node", (), {})()) # Making a "node" object that will contain properties for each node
      
      self.nodes[node].loc = nodes[node] # Coordinate of Node
      self.nodes[node].globalNum = int(nodeNums[node]) # Global Node number
      self.nodes[node].num = int(num[node]) # Reduced Node Number or -1 if in boundary meant for TM mode
      self.nodes[node].GradPhi = 0
      
      if(self.nodes[node].num != -1):
        scal, grad = basis3Dscalar(nodes, node, self.loc)
          
        self.nodes[node].N = scal
        self.nodes[node].gradN = grad
        
      else:
        self.nodes[node].N = 0
        self.nodes[node].gradN = 0
    
    # Edge Elements
    edgeNums = np.array([[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]])
    
    for edge in range(edgeNums.shape[0]):
      self.edges.append(type("Edge", (), {})())

      self.edges[edge].nodesLoc = nodes[edgeNums[edge]]
      self.edges[edge].globalNodesNum = nodeNums[edgeNums[edge]].reshape((1,2)) #Probably change the reshape      
      self.edges[edge].num = num[edgeNums[edge]]

      # Add edges to the list, and check whether they're already there to determine numbering
      listInd = np.where(np.all(self.edges[edge].globalNodesNum[0] == Tetrahedron.edgeList, axis=1))[0]
      if np.size(listInd) == 0:
        listInd = np.where(np.all(self.edges[edge].globalNodesNum[0][::-1] == Tetrahedron.edgeList, axis=1))[0] # Check the reverse order
        if np.size(listInd) == 0:
          Tetrahedron.edgeList = np.append(Tetrahedron.edgeList, self.edges[edge].globalNodesNum, axis = 0)
          self.edges[edge].globalNum = len(Tetrahedron.edgeList)-1
      else:
          self.edges[edge].globalNum = listInd[0]
      
      self.edges[edge].vectN = []
      self.edges[edge].curlN = []

      self.edges[edge].num = self.edges[edge].globalNum # Reduced Node Number or -1 if in boundary

      Tetrahedron.nodeList.append([int(self.edges[edge].globalNum), elNum, edge])
    
      self.edges[edge].A = np.zeros((self.loc.shape[0], 3), dtype = np.complex_)
    
    self.field = np.zeros((self.loc.shape[0], 3), dtype = np.complex_)
  
  def getBasisVectors(self):
    edgeNums = np.array([[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]])
    
    for edge in range(edgeNums.shape[0]):
      if(self.edges[edge].num != -1):
        self.edges[edge].vectN, self.edges[edge].curlN = basis3Dvector(self.nodeList, edgeNums[edge][0], edgeNums[edge][1], self.loc)

  def getA(self, a):
    for edge in self.edges:
      if(edge.num != -1):
        edge.A = edge.vectN * a[edge.num]
  
  def getGradPhi(self, d, w):
    for node in self.nodes:
      if(node.num != -1):     
        node.GradPhi += (d[node.num] * node.gradN / (1j * w))
  
  def getE(self, w):
    for edge in self.edges:
      if(edge.num != -1):
        self.field += 1j * w * edge.A
      
    for node in self.nodes:
      if(node.num != -1):
        self.field -= node.GradPhi
    
  '''Find if node number in element
  Inputs:
    nodeNum, node number to look for
    reduced = False, boolean if to look for in the global or reduced set of nodes
    
    Outputs:
    index, index of desired node
    -1 if not in element
  '''
  def findEdgeIndex(self, edgeNum, reduced = False):
    index = 0
  
    for edge in self.edges:
      if((edge.globalNum == edgeNum)):
        return index
      
      index += 1
    
    return -1  


'''
reIndex finds the indices of the indeces not on the boundary
and writes them into [element].edges[i].num
    Inputs:
    Elements, the list of mesh elements
    edgeList, the list of all edges on boundary
    
    Outputs:
    elements, the element list with new indices 
'''
def reIndex(elements: list, edgeList: list):
  index = np.zeros(len(Tetrahedron.edgeList))
  indeces = []
  indexShift = 0
  
  Tetrahedron.nodeList = Tetrahedron.nodeList[Tetrahedron.nodeList[:,0].argsort()]
  for i in Tetrahedron.nodeList:
    if(i[0] in edgeList):
      elements[i[1]].edges[i[2]].num = -1
      index[i[0]] = -1
    else:
      index[i[0]] = elements[i[1]].edges[i[2]].num

  newIndex = []
  
  for val in index:
    if(val == -1):
      indexShift -= 1
      newIndex.append(-1)
    else:
      newIndex.append(val + indexShift)
  newIndex = np.array(newIndex)

  for i in Tetrahedron.nodeList:
    elements[i[1]].edges[i[2]].num = int(newIndex[i[0]])
  
  return elements

def edgesOnBoudnary(boundaryList):
  boundaryEdges = np.empty((0,2))
  # elNums = np.empty((0,1))
  edgeNums = np.array([[0,1], [0,2], [1,2]])
    
  for el in boundaryList:
    for edge in range(edgeNums.shape[0]):
      
      boundaryEdges = np.append(boundaryEdges, el[edgeNums[edge]].reshape((1,2)), axis = 0)

  boundaryEdgesNums = []
  
  for edge in boundaryEdges:
    temp = np.where(np.all(Tetrahedron.edgeList == edge, axis = 1))[0]
    
    if(temp.size == 0):
      boundaryEdgesNums.append(np.where(np.all(Tetrahedron.edgeList == edge[::-1], axis = 1))[0][0])
    else:
      boundaryEdgesNums.append(np.where(np.all(Tetrahedron.edgeList == edge, axis = 1))[0][0])
  
  return np.unique(boundaryEdgesNums)


# Evaluates 3D scalar basis functions and their derivative over a tetrahedron
def basis3Dscalar22222222222222222222222222222(nodes: list, node_eval, points):
    # inputs:
    # n1, n2, n3, and n4 are the points of the tetrahedron
    # points is a list of points to evaluate at 
    # node_eval is the node (1, 2, 3, or 4) to evaluate from

    # outputs:
    # fvals is a list of the function value at the given points
    # gradvals is the gradient of the function in form [x,y,z]
    
    n1 = nodes[0]
    n2 = nodes[1]
    n3 = nodes[2]
    n4 = nodes[3]

    fvals = []
    gradvals = []

    # calculate determinant of tetrahedron 
    # print(nodes)
    tetra = np.array([[1, n1[0], n1[1], n1[2]],
                      [1, n2[0], n2[1], n2[2]],
                      [1, n3[0], n3[1], n3[2]],
                      [1, n4[0], n4[1], n4[2]]])
    # v = det(matrix)
    v = np.linalg.det(tetra)
    
    grad_x_matrix = np.copy(tetra)
    grad_x_matrix[node_eval] = [0,1,0,0]
    grad_y_matrix = np.copy(tetra)
    grad_y_matrix[node_eval] = [0,0,1,0]
    grad_z_matrix = np.copy(tetra)
    grad_z_matrix[node_eval] = [0,0,0,1]
    
    grad_x = np.linalg.det(grad_x_matrix) / v
    grad_y = np.linalg.det(grad_y_matrix) / v
    grad_z = np.linalg.det(grad_z_matrix) / v

    gradvals = [grad_x, grad_y, grad_z]
    
    for row in points:
      tetra1 = np.array([[1, n1[0], n1[1], n1[2]],
                          [1, n2[0], n2[1], n2[2]],
                          [1, n3[0], n3[1], n3[2]],
                          [1, n4[0], n4[1], n4[2]]]) 
      
      tetra1[node_eval,:] = [1, row[0], row[1], row[2]];
      
      val = (np.linalg.det(tetra1)) / v
      fvals.append(val)
    
    # return abs(np.array(fvals)), np.array(gradvals) * np.sign(fvals[0])
    return np.array(fvals), np.array(gradvals)
  

# Evaluates 3D scalar basis functions and their derivative over a tetrahedron
def basis3Dscalar(nodes: list, node_eval, points):
    # inputs:
    # n1, n2, n3, and n4 are the points of the tetrahedron
    # points is a list of points to evaluate at 
    # node_eval is the node (1, 2, 3, or 4) to evaluate from

    # outputs:
    # fvals is a list of the function value at the given points
    # gradvals is the gradient of the function in form [x,y,z]
    
    n1 = nodes[0]
    n2 = nodes[1]
    n3 = nodes[2]
    n4 = nodes[3]

    fvals = []
    gradvals = []

    # calculate determinant of tetrahedron 
    # print(nodes)
    tetra = np.array([[1, n1[0], n1[1], n1[2]],
                      [1, n2[0], n2[1], n2[2]],
                      [1, n3[0], n3[1], n3[2]],
                      [1, n4[0], n4[1], n4[2]]])
    # v = det(matrix)
    v = np.linalg.det(tetra)

    if (node_eval == 0): 
        # make gradient matrices 
        grad_x_matrix = np.copy(tetra)
        grad_x_matrix[0] = [0,1,0,0]
        grad_y_matrix = np.copy(tetra)
        grad_y_matrix[0] = [0,0,1,0]
        grad_z_matrix = np.copy(tetra)
        grad_z_matrix[0] = [0,0,0,1]

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
            val = (np.linalg.det(tetra1)) / v
            fvals.append(val)

    elif (node_eval == 1):
        # make gradient matrices 
        grad_x_matrix = np.copy(tetra)
        grad_x_matrix[1] = [0,1,0,0]
        grad_y_matrix = np.copy(tetra)
        grad_y_matrix[1] = [0,0,1,0]
        grad_z_matrix = np.copy(tetra)
        grad_z_matrix[1] = [0,0,0,1]

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
            val = (np.linalg.det(tetra2)) / v
            fvals.append(val)

    elif (node_eval == 2):
        # make gradient matrices 
        grad_x_matrix = np.copy(tetra)
        grad_x_matrix[2] = [0,1,0,0]
        grad_y_matrix = np.copy(tetra)
        grad_y_matrix[2] = [0,0,1,0]
        grad_z_matrix = np.copy(tetra)
        grad_z_matrix[2] = [0,0,0,1]

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
            val = (np.linalg.det(tetra3)) / v
            fvals.append(val)
    
    else: # node_eval == 4
        # make gradient matrices 
        grad_x_matrix = np.copy(tetra)
        grad_x_matrix[3] = [0,1,0,0]
        grad_y_matrix = np.copy(tetra)
        grad_y_matrix[3] = [0,0,1,0]
        grad_z_matrix = np.copy(tetra)
        grad_z_matrix[3] = [0,0,0,1]

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
            val = (np.linalg.det(tetra4)) / v
            fvals.append(val)

    # return abs(np.array(fvals)), np.array(gradvals) * np.sign(fvals[0])
    return np.array(fvals), np.array(gradvals)
  
def basis3Dvector(nodes: list, nEval1, nEval2, points):
    # inputs: 
    # n1, n2, n3 are the points of the triangle nodes 
    # nEval1 and nEval2 are the nodes (0, 1, or 2) to evaluate at, 0 < 1
    # points is a list of points to evaluate the function for 

    # output: 
    
    scal1, grad1 = basis3Dscalar(nodes, nEval1, points)
    scal2, grad2 = basis3Dscalar(nodes, nEval2, points)

    length = np.linalg.norm(nodes[nEval2] - nodes[nEval1])

    vectN = np.empty((scal1.shape[0], 3))
    
    vectN[:,0] = (scal1 * grad2[0] - scal2 * grad1[0])
    vectN[:,1] = (scal1 * grad2[1] - scal2 * grad1[1])
    vectN[:,2] = (scal1 * grad2[2] - scal2 * grad1[2])
    vectN *= length
    curlN = 2 * length * np.cross(grad1, grad2)
    curlN = np.tile(curlN, (vectN.shape[0],1))

    return vectN, curlN
