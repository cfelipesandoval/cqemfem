from mesh import *
import scipy as sc
from matplotlib import colors
import time


"""TM modes FEM Program
Program returns desired TM mode from
ordered 
"""
def main_TM(modes, fileName):
  print("Beginning TM mode FEM")
  
  # Element Initialization
  start_time = time.time()
  print("...Initializing elements...")
  elements, nodeList, boundaryList = initMesh(fileName, 9)
  boundaryList = np.unique(boundaryList)
  print(f"Done, time elapsed = {time.time() - start_time} s")
  
  if(max(modes) + 1 > (len(nodeList) - len(boundaryList))):
    print(f"Error, not enough elements to calculate modes above {(len(nodeList) - len(boundaryList))}")
    sys.exit()
  else:
    kMax = max(modes) + 1 # Max eigevalue/vector to calculate
    
  # Re indexing elements to only take into account those not on the boundary
  print("...Reindexing Elements...")
  start_time = time.time()
  elements = reIndex(elements, nodeList)
  print(f"Done, time elapsed = {time.time() - start_time} s")
  
  # Initializing matrices from equation AHz = (Kc^2)BHz
  A = np.zeros((len(nodeList) - len(boundaryList), len(nodeList) - len(boundaryList)), dtype=np.complex_)
  B = np.zeros((len(nodeList) - len(boundaryList), len(nodeList) - len(boundaryList)), dtype=np.complex_)
  
  # Calculating FEM
  start_time = time.time()
  print("Beginning FEM")
  for el in elements:
    for i in range(0, 3):
      for j in range(0, 3):
        if(el.nodes[i].num == -1 or el.nodes[j].num == -1):
          continue
        
        # Calculating matrix elements for A = ∫∫(∇Ni⋅∇Nj)dΩ and B = ∫∫(Ni⋅Nj)dΩ
        A[el.nodes[i].num, el.nodes[j].num] += q.gauss_quadrature_trisimp(el.nodes[0].loc, el.nodes[1].loc, el.nodes[2].loc, 
        (np.sum(el.nodes[i].grad_N * el.nodes[j].grad_N, axis=1)))
        B[el.nodes[i].num, el.nodes[j].num] += q.gauss_quadrature_trisimp(el.nodes[0].loc, el.nodes[1].loc, el.nodes[2].loc, 
        (el.nodes[i].N * el.nodes[j].N))
  print(f"FEM Done, time elapsed = {time.time() - start_time} s")

  # Finding eigen values and vectors (6 first by default)
  start_time = time.time()
  print("Finding eigen values and vectors")
  Ek, Ev = sc.sparse.linalg.eigs(A, M=B, which = "SM", k = kMax)
  E_eigVals = np.array(Ek)
  E_eigVecs = np.array(Ev)
  print(f"Eigen Done, time elapsed = {time.time() - start_time} s")

  fig = plt.figure(figsize=(3.5,5))
  fig.suptitle("TM Modes", fontsize = 15)
  index = 1
  # Plotting
  for mode in modes:
    start_time = time.time()
    print("Plotting")
    
    E_eigVecs = np.real(E_eigVecs)

    
    # ax = fig.add_subplot(projection="3d")
    ax = fig.add_subplot(2, 1, index, projection="3d")
    index += 1
    ax.view_init(azim=-0.01, elev=89.99)
    
    x = []
    y = []
    z = []
    for el in elements:
      for node in el.nodes:
        if(node.num == -1):
          x.append(node.loc[0])
          y.append(node.loc[1])
          z.append(0)
          continue

        x.append(node.loc[0])
        y.append(node.loc[1])
        z.append(E_eigVecs[:,mode][node.num])

    ax.set_xlim(-1,1)
    ax.set_ylim(-1,1)

    ax.w_zaxis.line.set_lw(0.)
    ax.set_zticks([])
    ax.tick_params(axis='both', labelsize=8)
    ax.grid(False)
    
    # ax.set_title(f"TM mode[{mode}]\nKc^2 = {np.round(np.real(E_eigVals[mode]), 3)}")
    # ax.set_title(f"TM mode[{mode}]")
    ax.set_ylabel("Length")
    ax.set_xlabel("Width")

    p = ax.plot_trisurf(x, y, np.real(z), cmap="coolwarm")  
    
    cb = fig.colorbar(p, extend="both")
    # cb.set_label("Re{Ez}")
        
    print(f"Plotting Done, time elapsed = {time.time() - start_time} s") 
    print(f"Current mode kc^2 = {np.real(E_eigVals[mode])}")
    
  # Analytical cutoff frequencies plotting
  # fig2 = plt.figure(figsize=(10, 6))
  # ax2 = fig2.add_subplot()
  
  # analytical = np.array([[(i * np.pi / 2) ** 2 + (j * np.pi / 1) ** 2 for i in range(0,5)] for j in range(0,5)])
  
  # x = np.array([i for i in range(0,8)])
  
  # y1 = np.real(E_eigVals[0:8])
  # y2 = [analytical[1,1], analytical[1,2], analytical[1,3], analytical[2,1], analytical[1,4], analytical[2,2], analytical[2,3], analytical[2,4]]
  # print(analytical)
  # print(np.real(E_eigVals))
  # ax2.set_title("Comparison between FEM calculated and Analytical Results")
  
  # ax2.plot(x, y1,'ro', label='FEM')
  # ax2.plot(x, y2,'bo', label='Analytical')
  # ax2.set_xlabel("Mode")
  # ax2.set_ylabel("Kc^2")
  # ax2.legend()
  
  # plt.show()
  
  print("TM mode FEM Done")  
       
# if __name__ == "__main__":
#   main_TM([0, 4], "rectangularMesh.inp")  
#   plt.show()