import numpy as np
import matplotlib.pyplot as plt

i_step, tau, N, E_T = np.genfromtxt("data/task1.csv", delimiter=",", unpack=True)

E0_exact = 3/8

plt.figure(figsize=(8,4))
plt.plot(i_step, N)
plt.xlabel("iteration step")
plt.ylabel("$N$ (number of walkers)")
plt.tight_layout()
plt.savefig("plots/task1_N.pdf")

plt.figure(figsize=(8,4))
plt.axhline(E0_exact, color="k", alpha=0.5, linestyle="--")
plt.plot(i_step, E_T)
plt.xlabel("iteration step")
plt.ylabel("$E_T \: / \: $ (a.u.)")
plt.tight_layout()
plt.savefig("plots/task1_E_T.pdf")