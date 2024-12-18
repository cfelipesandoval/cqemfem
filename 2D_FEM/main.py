from main_TM import *
from main_TE import *

def main(modes, fileName):
  # main_TM(modes, fileName)
  main_TE(modes, fileName)

if __name__ == "__main__":
  main([0,9], "2D_FEM/rectangularMesh.inp") # Currently have rectangular, circular, and elliptical mesh examples
  plt.show()