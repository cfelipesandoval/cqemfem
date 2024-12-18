import ReadAbaqusMesh as abq
import gauss_quadrature_tetra as q

import numpy as np
import matplotlib.pyplot as plt
import sys
import time


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
  elementList = np.array(elementList) - 1
  boundaryList = np.array(boundaryList) - 1 # -1 to take matlab offset into account

  nodeList[:,0] = nodeList[:,0] - 1 # -1 to take matlab offset into account
 
  elements = [] # Initializing Output element list
  
  
  start_time = time.time()
  for el in range(0,len(elementList)):
    coords = []
    nodeNums = []
    
    for i in elementList[el]:
      coords.append(nodeList[(nodeList[:,0] == i), (1,2,3)]) # Creating list with coordinates for each node
      nodeNums.append(nodeList[(nodeList[:,0] == i), 0]) # Creating list with node number for each node
    coords = np.array(coords)
    
    elements.append(Tetrahedron(el, coords, nodeNums, quadPoints, boundaryList)) # Making triangle element object
  print(f'time: {time.time()-start_time}')
  return elements, nodeList, boundaryList

class Tetrahedron:
  edgeList = np.empty((0,2))
  nodeList = []
  
  def __init__(self, elNum, nodes: list, nodeNums, points, boundaryList):
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
    self.elNum = elNum
    self.nodeList = np.copy(nodes)
    nodeNums = np.array(nodeNums)
    self.edges = []
    tempNodes = []
    
    # Initializing
    self.gaussPoints = points # These are the number of points for the quadrature

    self.loc = q.gauss_quadrature_tetra(nodes, self.gaussPoints) # Initializing location of discretized points within element
    self.loc = np.array(self.loc)
    
    inds = nodeNums[:,0].argsort()
    nodes = nodes[inds,:]
    nodeNums = nodeNums[inds]
    
    # Edge elements
    edgeNums = np.array([[0,1], [0,2], [0,3], [1,2], [1,3], [2,3]])

    for edge in range(edgeNums.shape[0]):
      self.edges.append(type("Edge", (), {})())

      self.edges[edge].nodesLoc = nodes[edgeNums[edge]]
      self.edges[edge].globalNodesNum = nodeNums[edgeNums[edge]].reshape((1,2)) #Probably change the reshape      
      
      self.edges[edge].edgeDir = (self.edges[edge].nodesLoc[1] - self.edges[edge].nodesLoc[0]) / (np.linalg.norm(self.edges[edge].nodesLoc[1] - self.edges[edge].nodesLoc[0]))
      
      # Initializing node properties
      for i in range(0,2):
        tempNodes.append(type("Node", (), {})()) # Making a "node" object that will contain properties for each node
        tempNodes[i].loc = []
        tempNodes[i].globalNum = []
        
        tempNodes[i].loc = (self.edges[edge].nodesLoc[i]) # Coordinate of Node
        tempNodes[i].globalNum = (self.edges[edge].globalNodesNum[0][i]) # Global Node number
        
        scal, grad = basis3Dscalar(nodes, i, self.loc)
        tempNodes[i].N = scal
        tempNodes[i].gradN = grad
      
      self.edges[edge].nodes = tempNodes

      if(self.edges[edge].globalNodesNum[0].tolist() not in Tetrahedron.edgeList.tolist() and self.edges[edge].globalNodesNum[0][::-1].tolist() not in Tetrahedron.edgeList.tolist()): 
        
        self.edges[edge].vectN, self.edges[edge].curlN = basis3Dvector(nodes, edgeNums[edge][0], edgeNums[edge][1], self.loc)
        
        Tetrahedron.edgeList = np.append(Tetrahedron.edgeList, self.edges[edge].globalNodesNum, axis = 0) 
        
        self.edges[edge].globalNum = np.where(np.all(Tetrahedron.edgeList == self.edges[edge].globalNodesNum, axis = 1))[0]
        # Tetrahedron.edgeThing.append([int(self.edges[edge].globalNum[0]), elNum, edge])
      else:

        self.edges[edge].globalNum = np.where(np.all(Tetrahedron.edgeList == self.edges[edge].globalNodesNum, axis = 1))[0]

        if(not self.edges[edge].globalNum.size):
          self.edges[edge].globalNum = np.where(np.all(Tetrahedron.edgeList == self.edges[edge].globalNodesNum[0][::-1], axis = 1))[0]
        
        self.edges[edge].vectN, self.edges[edge].curlN = basis3Dvector(nodes, edgeNums[edge][0], edgeNums[edge][1], self.loc)
        
      self.edges[edge].num = self.edges[edge].globalNum # Reduced Node Number or -1 if in boundary

      Tetrahedron.nodeList.append([int(self.edges[edge].globalNum[0]), elNum, edge])
      
    self.field = np.array([0] * self.loc, dtype=np.complex_)
    self.avg = np.array([0] * 3)

          
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

def getVects(elements, eigVects, mode):
  for el in elements:
    for edge in el.edges:
      if(edge.num == -1):
        continue
      
      el.field += edge.vectN * eigVects[edge.num,mode]
  
  return elements
'''
THIS SHOULD BE RE-WRITTEN (SUPER SLOW)
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

  # for ind in Tetrahedron.edgeThing: 
  #   if(elements[ind[1]].edges[ind[2]].globalNum in edgeList):
  #     indeces.append(ind[0])

  for i in Tetrahedron.nodeList:
    if(i[0] in edgeList):
      elements[i[1]].edges[i[2]].num = -1
      index[i[0]] = -1
    else:
      index[i[0]] = elements[i[1]].edges[i[2]].num

  # # THIS IS SUPER SLOW
  # for i in range(len(Tetrahedron.edgeList)):
  #   for el in elements:
  #     for edge in el.edges:
  #       index[edge.globalNum] = edge.num
  
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
  elNums = np.empty((0,1))
  edgeNums = np.array([[0,1], [0,2], [1,2]])
    
  for el in boundaryList:
    for edge in range(edgeNums.shape[0]):
      
      # if(el[edgeNums[edge]] not in boundaryEdges):
      #   print(el[edgeNums[edge]])
      boundaryEdges = np.append(boundaryEdges, el[edgeNums[edge]].reshape((1,2)), axis = 0)
  # print(boundaryEdges)

  boundaryEdgesNums = []
  
  for edge in boundaryEdges:
    temp = np.where(np.all(Tetrahedron.edgeList == edge, axis = 1))[0]
    
    if(temp.size == 0):
      boundaryEdgesNums.append(np.where(np.all(Tetrahedron.edgeList == edge[::-1], axis = 1))[0][0])
    else:
      boundaryEdgesNums.append(np.where(np.all(Tetrahedron.edgeList == edge, axis = 1))[0][0])
  
  return np.unique(boundaryEdgesNums)
    
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
    v = abs(np.linalg.det(tetra))

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
    
    # n1 = nodes[0]
    # n2 = nodes[1]
    # n3 = nodes[2]
    # n4 = nodes[3]
    
    # Eval1 = nEval1
    # Eval2 = nEval2

    scal1, grad1 = basis3Dscalar(nodes, nEval1, points)
    scal2, grad2 = basis3Dscalar(nodes, nEval2, points)
    
    # grad1 = grad1 / np.linalg.norm(grad1)
    # grad2 = grad2 / np.linalg.norm(grad2)

    length = np.linalg.norm(nodes[nEval2] - nodes[nEval1])
    
    # if nodes are 0 and 1
    # if ((Eval1 == 0 and Eval2 == 1) or (Eval1 == 1 and Eval2 == 0)):
    #   length = np.sqrt((n1[0] - n2[0])**2 + (n1[1] - n2[1])**2 + (n1[2] - n2[2])**2)

    # # if nodes are 0 and 2
    # elif ((Eval1 == 0 and Eval2 == 2) or (Eval1 == 2 and Eval2 == 0)):
    #   length = np.sqrt((n1[0] - n3[0])**2 + (n1[1] - n3[1])**2 + (n1[2] - n3[2])**2)

    # # if nodes are 0 and 3
    # elif ((Eval1 == 0 and Eval2 == 3) or (Eval1 == 3 and Eval2 == 0)):
    #   length = np.sqrt((n1[0] - n4[0])**2 + (n1[1] - n4[1])**2 + (n1[2] - n4[2])**2)
        
    # # if nodes are 1 and 2
    # elif ((Eval1 == 1 and Eval2 == 2) or (Eval1 == 2 and Eval2 == 1)):
    #   length = np.sqrt((n2[0] - n3[0])**2 + (n2[1] - n3[1])**2 + (n2[2] - n3[2])**2)
    
    # # if nodes are 1 and 3
    # elif ((Eval1 == 1 and Eval2 == 3) or (Eval1 == 3 and Eval2 == 1)):
    #   length = np.sqrt((n2[0] - n4[0])**2 + (n2[1] - n4[1])**2 + (n2[2] - n4[2])**2)
    
    # # if nodes are 2 and 3
    # elif ((Eval1 == 2 and Eval2 == 3) or (Eval1 == 3 and Eval2 == 2)):
    #   length = np.sqrt((n3[0] - n4[0])**2 + (n3[1] - n4[1])**2 + (n3[2] - n4[2])**2)

    vectN = np.empty((scal1.shape[0], 3))
    
    vectN[:,0] = (scal1 * grad2[0] - scal2 * grad1[0])
    vectN[:,1] = (scal1 * grad2[1] - scal2 * grad1[1])
    vectN[:,2] = (scal1 * grad2[2] - scal2 * grad1[2])
    vectN *= length
    curlN = 2 * length * np.cross(grad1, grad2)
    # print(curlN)
    curlN = np.tile(curlN, (vectN.shape[0],1))
    # print(curlN)
    return vectN, curlN
