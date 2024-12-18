import matplotlib.pyplot as plt
import numpy as np
error = [1.9056680394078117e-06j, 3.603267287969239e-07j, 6.630187057153376e-08j, 1.4121327632695936e-08j, 3.185536836229494e-09j, 1.0728143071664621e-09]
error = np.imag(error)
error = error / max(error)

plt.plot(range(len(error)),(error),'bo')
# plt.grid()
plt.ylabel("Normalized Error", fontsize = 20)
plt.xlabel("Iteration", fontsize = 20)
plt.title("Adaptive Mesh Refinement Iteration Error", fontsize = 20)
plt.show()