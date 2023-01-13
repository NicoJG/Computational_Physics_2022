import numpy as np
import matplotlib.pyplot as plt

N, E_T = np.genfromtxt("walkers.csv", delimiter=",", unpack=True)

fig,axs = plt.subplots(2,1,figsize=(5,8))

plt.sca(axs[0])
plt.plot(N)

plt.sca(axs[1])
plt.plot(E_T)

plt.savefig("H3.pdf")