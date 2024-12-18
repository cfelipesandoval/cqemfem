import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors

def getFields(elements, w, a, d):
  
  
  #aArray = aArray[:,np.real(kArray) > 5]
  #wArray = wArray[np.real(kArray) > 5]
  #kArray = kArray[np.real(kArray) > 5]
  
  for el in elements:
    el.getA(a)
    el.getGradPhi(d, w)
    el.getE(w)
  

def plotMode(elements, modeInfo):
  modes = modeInfo[0]
  axis = modeInfo[1]
  
  done = []
  out = np.empty((0,2))
  
  for el in elements:
    if(el.num in done):
      continue
    
    done.append(el.num)
    
    out = np.append(out, np.array((el.num, np.imag(el.field[0,modeInfo[1]]))).reshape(1,2), axis = 0)
    # print(out)
  with open('test.txt', 'w') as f:
    # np.save(f, out)
    for l in out:
      f.write(str(l))
  # with open('test.txt', 'r') as f:
  #   print(np.load(f))
  
  # print(f'w: {np.real(w)}')
  # print(f'k: {np.real(k)}')
  
  
  fig = plt.figure()
  # fig.suptitle(f"Mode: {modes[mode]}", fontsize = 15)

  ax = fig.add_subplot(projection="3d")
  
  # ax.w_zaxis.line.set_lw(0.)
  # ax.set_zticks([])
  # ax.tick_params(axis='both', labelsize=8)
  # ax.grid(False)

  x1 = []
  y1 = []
  z1 = []
  c1 = []
  
  x2 = []
  y2 = []
  z2 = []
  c2 = []
  
  sliceRange = 0.2
  
  for el in elements:
    loc = np.average(el.loc, axis = 0)
    if(loc[1] < sliceRange and loc[1] > -sliceRange):
      x1.extend(el.loc[:,0].tolist())
      # y.append(np.zeros(el.loc[:,0].shape))
      y1.extend(el.loc[:,1].tolist())
      z1.extend(el.loc[:,2].tolist())
      c1.extend(np.imag((el.field[:, axis])).tolist())
    
    if(loc[0] < sliceRange and loc[0] > -sliceRange):
      x2.extend(el.loc[:,0].tolist())
      # y.append(np.zeros(el.loc[:,0].shape))
      y2.extend(el.loc[:,1].tolist())
      z2.extend(el.loc[:,2].tolist())
      c2.extend(np.imag((el.field[:, axis])).tolist())

  x1 = np.array(x1)
  y1 = np.array(y1)
  z1 = np.array(z1)
  c1 = np.array(c1)

  x2 = np.array(x2)
  y2 = np.array(y2)
  z2 = np.array(z2)
  c2 = np.array(c2)    
  
  # c = c / max(abs(c))
  
  match axis:
    case 0:
      x1 *= 0
      y2 *= 0
    
    case 1:
      y1 *= 0
      x2 *= 0
      # ax.view_init(azim=-0.01, elev=89.99)
    
    case 2:
      z1 *= 0
      y2 *= 0

      
  x = np.append(x1, x2)
  y = np.append(y1, y2)
  z = np.append(z1, z2)
  c = np.append(c1, c2)
  
  maxX = max(np.imag(elements[0].field[:,0]))
  minX = min(np.imag(elements[0].field[:,0]))
  
  maxY = max(np.imag(elements[0].field[:,1]))
  minY = min(np.imag(elements[0].field[:,1]))
  
  maxZ = max(np.imag(elements[0].field[:,2]))
  minZ = min(np.imag(elements[0].field[:,2]))
  
  for el in elements:
    if(max(np.imag(el.field[:,0])) > maxX):
      maxX = max(np.imag(el.field[:,0]))
    if(min(np.imag(el.field[:,0])) < minX):
      minX = min(np.imag(el.field[:,0]))
      
    if(max(np.imag(el.field[:,1])) > maxY):
      maxY = max(np.imag(el.field[:,1]))
    if(min(np.imag(el.field[:,1])) < minY):
      minY = min(np.imag(el.field[:,1]))
      
    if(max(np.imag(el.field[:,2])) > maxZ):
      maxZ = max(np.imag(el.field[:,2]))
    if(min(np.imag(el.field[:,2])) < minZ):
      minZ = min(np.imag(el.field[:,2]))
  
  maxs = [maxX, maxY, maxZ]
  mins = [minX, minY, minZ]
  
  n = colors.Normalize(vmin=mins[axis], vmax=maxs[axis])
  
  p = ax.scatter(x, y, z, c=c, cmap ='Spectral_r', norm=n)
  # for el in elements:
  #   p = ax.scatter(loc[:,0], loc[:,1], loc[:,2], c=el.field[:,axis], cmap ='Spectral_r', norm=n)
      
  fig.colorbar(p, extend="min")
  ax.grid(False)
  ax.axis(False)

# Hide axes ticks
  ax.set_xticks([])
  ax.set_yticks([])
  ax.set_zticks([])
  plt.show()


def errorPlot(kArray):
  anModes = np.array([[1,0,1], [1,1,0], [0,1,1], [2,0,1], [1,1,1], [1,1,1], [2,1,0], [1,0,2]])
  a = 0.01
  b = 0.005
  c = 0.0075
  
  ks = np.sqrt((anModes[:,0] * np.pi / a) ** 2 + (anModes[:,1] * np.pi / b) ** 2 + (anModes[:,2] * np.pi / c) ** 2)
  
  sorted = np.sort(kArray)
  
  x = range(0,len(ks))
  
  plt.plot(x, ks, "rx", markersize = 10)
  plt.plot(x, sorted[0:8], "bo")
  plt.legend(["Analytical Solutions", "FEM Solutions"], fontsize = 20)
  plt.title("Analytical and FEM Results", fontsize = 20)
  plt.xlabel("Mode", fontsize = 20)
  plt.ylabel("Cutoff Wavenumber ($m^{-1}$)", fontsize = 20)
  
  # # error = np.abs(ks-sorted[0:8])/ks * 100
  # # plt.plot(x, error, "bo")
  
  plt.show()