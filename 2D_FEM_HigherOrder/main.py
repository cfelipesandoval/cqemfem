from main_TM import *
from main_TE import *

def main(modes, fileName, order):
  # main_TM(modes, fileName)
  main_TE(modes, fileName, order)

if __name__ == "__main__":
  main([1,9], "./rectangularMesh.inp", 3) # Currently have rectangular, circular, and elliptical mesh examples
  plt.show()