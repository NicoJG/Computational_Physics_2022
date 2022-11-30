# %%
import numpy as np
import matplotlib.pyplot as plt

i_step, E, P, r = np.genfromtxt("../data/H2a_sim.csv", delimiter=",", unpack=True)


# %%
plt.plot(i_step, E)
plt.show()
plt.plot(i_step, P)
plt.show()
plt.plot(i_step, r)
plt.show()
# %%
