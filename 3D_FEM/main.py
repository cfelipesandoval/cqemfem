import numpy as np
from mesh import *
import gauss_quadrature_tetra as q
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib import cm
import scipy as sc
import time

def main(fileName, mode): 
  #%% Initialization
  print("Beginning FEM")
  start_time = time.time() 
  print("...Initializing elements...")
  elements, nodeList, boundaryList = initMesh(fileName, 4)
  print(f"Done, time elapsed = {time.time() - start_time} s")
  
  numberOfEdges = len(Tetrahedron.edgeList)
  
  boundaryEdgeNums = edgesOnBoudnary(boundaryList)
  
  Tetrahedron.nodeList = np.array(Tetrahedron.nodeList)

  unknowns = numberOfEdges - len(boundaryEdgeNums)
  
  print(f"Number of Unknown Edges: {unknowns}")
  
  print("...Reindexing Elements...")
  start_time = time.time()
  elements = reIndex(elements, boundaryEdgeNums)
  print(f"Done, time elapsed = {time.time() - start_time} s")
  #%%
  #%% FEM
  A = np.zeros((unknowns, unknowns), dtype=np.complex_)
  B = np.zeros((unknowns, unknowns), dtype=np.complex_)
  
  
  start_time = time.time()
  print("Beginning FEM")

  for el in elements:
    for i in range(0, 6):
      for j in range(0, 6):
        if(el.edges[i].num == -1 or el.edges[j].num == -1):
          continue
        
        A[el.edges[i].num, el.edges[j].num] += q.gauss_quadrature_tetra(el.nodeList, el.gaussPoints, 
        np.sum(el.edges[i].curlN * el.edges[j].curlN, axis = 1))

        B[el.edges[i].num, el.edges[j].num] += q.gauss_quadrature_tetra(el.nodeList, el.gaussPoints, 
        np.sum(el.edges[i].vectN * el.edges[j].vectN, axis = 1))
  
  print(f"FEM Done, time elapsed = {time.time() - start_time} s")
  #%%
  
  #%% Linear Algebra
  start_time = time.time()
  print("Finding eigen values and vectors")
  
  Ek, Ev = sc.sparse.linalg.eigs(A, M=B, which='SM', k=50)#unknowns-2)
  E_eigVals = np.array(Ek)
  E_eigVecs = np.array(Ev)
  
  # Sorting Eigenvalues
  E_eigVecs = E_eigVecs[:,E_eigVals > 0.5]
  E_eigVals = E_eigVals[E_eigVals > 0.5]
  
  # plot_mode = 0
  
  # elements = getVects(elements, E_eigVecs, plot_mode) # Adding vector basis in each element
  
  print(f"Eigen Done, time elapsed = {time.time() - start_time} s")
  #%% Extra
  
  analytical = np.array([5.236, 7.025, 7.531, 7.531, 8.179, 8.179, 8.886, 8.947])
  #%% Plotting 
  index = 1

  fig = plt.figure()
  fig.suptitle(f"Mode: {mode}", fontsize = 15)
  
  
  for i in [0,1,2]:
    ax = fig.add_subplot(1,3, i + 1, projection="3d")
    defElements = np.copy(elements)
    newElements = getVects(defElements, E_eigVecs, mode) # Adding vector basis in each element
    
    maxX = max(np.real(newElements[0].field[:,0]))
    minX = min(np.real(newElements[0].field[:,0]))
    
    for el in newElements:
      if(max(np.real(el.field[:,0])) > maxX):
        maxX = max(np.real(el.field[:,0]))
      if(min(np.real(el.field[:,0])) < minX):
        minX = min(np.real(el.field[:,0]))
    
    maxY = max(np.real(newElements[0].field[:,1]))
    minY = min(np.real(newElements[0].field[:,1]))
    
    for el in newElements:
      if(max(np.real(el.field[:,1])) > maxY):
        maxY = max(np.real(el.field[:,1]))
      if(min(np.real(el.field[:,1])) < minY):
        minY = min(np.real(el.field[:,1]))
        
    maxZ = max(np.real(newElements[0].field[:,2]))
    minZ = min(np.real(newElements[0].field[:,2]))
    
    for el in newElements:
      if(max(np.real(el.field[:,2])) > maxZ):
        maxZ = max(np.real(el.field[:,2]))
      if(min(np.real(el.field[:,2])) < minZ):
        minZ = min(np.real(el.field[:,2]))
    
    nX = colors.Normalize(vmin=minX, vmax=maxX)
    nY = colors.Normalize(vmin=minY, vmax=maxY)
    nZ = colors.Normalize(vmin=minZ, vmax=maxZ)
    
    normList = [nX, nY, nZ]
    
    # print(maxX, maxY, maxZ)
    # print(minX, minY, minZ)
    
    # rangeList = [maxX - minX, maxY - minY, maxZ - minZ]
    
    # component = rangeList.index(max(rangeList))
  
    # ax = fig.add_subplot(projection='3d')
    
    index += 1
    # ax.view_init(azim=90, elev=0)
    # ax.w_yaxis.line.set_lw(0.)
    # ax.set_yticks([])
    # ax.tick_params(axis='both', labelsize=8)
    # ax.grid(False)
    
    for el in newElements:
      ax.scatter(el.loc[:,0], el.loc[:,1], el.loc[:,2], c=np.real((el.field[:,i])), cmap='coolwarm', norm=normList[i])

  plt.show()
  #%%
  

if __name__ == "__main__":
  # main("3D_FEM/cubicMesh.inp")
  main("3D_FEM/1x0.5x0.75Mesh.inp", 14)