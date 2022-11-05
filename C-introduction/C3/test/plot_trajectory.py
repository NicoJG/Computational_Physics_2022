import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt('trajectory.dat')
plt.plot(data[:, 1], data[:, 2]) # Plot x versus y
plt.xlabel('x')
plt.xlabel('y')
plt.savefig("trajectory.pdf")