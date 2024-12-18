import numpy as np
from mesh import *
from plotModes import *
import gauss_quadrature_tetra as q
import scipy as sc
import time
import copy


def FEM(fileName): 
  begin_time = time.time()
  #%% Initialization
  print("Beginning FEM")
  
  start_time = time.time() 
  print("...Initializing elements...")
  elements, nodeList, boundaryList, unknownNodes, unknownEdges = initMesh(fileName, 4)
  
  print(f"Done, time elapsed = {time.time() - start_time} s")
  
  print(f"Number of Unknown Nodes: {unknownNodes}")
  print(f"Number of Unknown Edges: {unknownEdges}")

  #%% FEM
  K1 = np.zeros((unknownEdges, unknownEdges), dtype=np.complex_)
  K2 = np.zeros((unknownEdges, unknownEdges), dtype=np.complex_)
  Ken = np.zeros((unknownEdges, unknownNodes), dtype=np.complex_)
  Gnn = np.zeros((unknownNodes, unknownNodes), dtype=np.complex_)
  Kne = np.zeros((unknownNodes, unknownEdges), dtype=np.complex_)

  start_time = time.time()  
  print("Beginning FEM")
  
  for el in elements:
    for i in range(0, 6):
      if(el.edges[i].num == -1):
        continue
      for j in range(0, 6):
        if(el.edges[j].num == -1):
          continue
        
        K1[el.edges[i].num, el.edges[j].num] += q.gauss_quadrature_tetra(el.nodeList, el.gaussPoints, 
        (1/el.mu) * np.sum(el.edges[i].curlN * el.edges[j].curlN, axis = 1))
        
        K2[el.edges[i].num, el.edges[j].num] += q.gauss_quadrature_tetra(el.nodeList, el.gaussPoints, 
        el.ep * np.sum(el.edges[i].vectN * el.edges[j].vectN, axis = 1))

  for el in elements:
    for i in range(0,6):
      if(el.edges[i].num == -1):
        continue
      for j in range(0,4):
          if(el.nodes[j].num == -1):
            continue          
          
          Ken[el.edges[i].num, el.nodes[j].num] += q.gauss_quadrature_tetra(el.nodeList, el.gaussPoints, 
          el.ep * np.sum(el.edges[i].vectN * el.nodes[j].gradN, axis = 1))
          
          Kne[el.nodes[j].num, el.edges[i].num] += q.gauss_quadrature_tetra(el.nodeList, el.gaussPoints, 
          -el.ep / (mu0 * ep0 ** 2) * np.sum(el.edges[i].vectN * el.nodes[j].gradN, axis = 1))

  for el in elements:
    for i in range(0,4):
      if(el.nodes[i].num == -1):
        continue
      for j in range(0,4):
        if(el.nodes[j].num == -1):
          continue
                
        Gnn[el.nodes[i].num, el.nodes[j].num] += q.gauss_quadrature_tetra(el.nodeList, el.gaussPoints, 
        el.nodes[i].N * el.nodes[j].N) 
  
  print(f"FEM Done, time elapsed = {time.time() - start_time} s")
  
  #%% Linear Algebra
  start_time = time.time()
  print("Finding eigen values and vectors")

  A = K1 - (Ken @ np.linalg.inv(Gnn) @ Kne)
  B = K2
  
  w2, aArray = sc.sparse.linalg.eigs(A, M=B, which = 'SM', k = 100)
  # w2, aArray = sc.linalg.eig(A, b=B)

  print(f"Eigen Done, time elapsed = {time.time() - start_time} s")
  print(f'Run Time: {time.time() - begin_time} s')
  
  dArray = np.linalg.inv(Gnn) @ Kne @ aArray

  w2 = np.array(w2)
  aArray = np.array(aArray)
  dArray = np.array(dArray)
  
  wArray = np.sqrt(w2)

  w2_inds = np.argsort(w2)
  aArray = aArray[:,w2_inds]
  dArray = dArray[:,w2_inds]
  wArray = wArray[w2_inds]
  
  return elements, wArray, aArray, dArray
  

if __name__ == "__main__":
  ## modeInfo is of form [[mode, axis]]
  # with open('test.npy', 'rb') as f:
  # with open('test.txt', 'r') as f:
  #   print((f.readlines))
  #   print(np.load(f))
  
  filename = r'1x0.5x0.75cm.inp'
  # filename = r"C:\Users\selkin\OneDrive - purdue.edu\Documents\Research\Code\A-phi-Testing\Meshes\microstrip.inp"
  elements, wArray, aArray, dArray = FEM(filename)
    
  modeInfo = np.array([[0, 1]])

  for info in modeInfo:
    elms = copy.deepcopy(elements)
    w = wArray[info[0]]
    a = aArray[:,info[0]]
    d = dArray[:,info[0]]
    
    getFields(elms, w, a, d)
    plotMode(elms, info)