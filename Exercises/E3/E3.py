import numpy as np
import matplotlib.pyplot as plt
import json

##################################
# Task 1
##################################
N, I, var = np.genfromtxt("data/E3_1.csv", delimiter=",", unpack=True)

plt.figure(figsize=(5,4))
plt.errorbar(N, I, np.sqrt(var), fmt="x", label="MC integration (uniform)")
plt.axhline(1./6., linestyle="--", color="k", label="Exact value")
plt.xscale("log")
plt.xlabel("N")
plt.ylabel(r"$\int_0^1 x(1-x) \;\mathrm{d}x$")
plt.legend()
plt.tight_layout()
plt.savefig("plots/E3_1.png")

##################################
# Task 2
##################################
N2, I2, var2 = np.genfromtxt("data/E3_2.csv", delimiter=",", unpack=True)

plt.figure(figsize=(5,4))
plt.errorbar(N, I, np.sqrt(var), fmt="x", label="MC integration (uniform)")
plt.errorbar(N2, I2, np.sqrt(var2), fmt="x", label="MC integration (importance)")
plt.axhline(1./6., linestyle="--", color="k", label="Exact value")
plt.xscale("log")
plt.xlabel("N")
plt.ylabel(r"$\int_0^1 x(1-x) \;\mathrm{d}x$")
plt.legend()
plt.tight_layout()
plt.savefig("plots/E3_2.png")

# compare errors
plt.figure(figsize=(5,4))
plt.plot(N, np.sqrt(var), "o", label="MC integration (uniform)")
plt.plot(N2, np.sqrt(var2), "o", label="MC integration (importance)")
plt.xscale("log")
plt.yscale("log")
plt.xlabel("N")
plt.ylabel(r"$\sigma$")
plt.legend()
plt.tight_layout()
plt.savefig("plots/E3_2_sigma.png")

# check the x points
x, f, p = np.genfromtxt("data/E3_2_x.csv", delimiter=",", unpack=True)

plt.figure(figsize=(5,4))
plt.hist(x, bins=50, density=True, label="histogram (density)")
sort_idxs = np.argsort(x)
plt.plot(x[sort_idxs],p[sort_idxs], "-", label="p(x)")
plt.plot(x[sort_idxs],f[sort_idxs], "-", label="f(x)")
#plt.axhline(np.pi/2, linestyle="--", color="k", alpha=0.5, label="$\\pi/2$")
plt.xlabel("x")
plt.legend()
plt.tight_layout()
plt.savefig("plots/E3_2_x.png")


##################################
# Task 3
##################################
# get header metadata
with open("data/E3_3_samples.csv", "r") as file:
    metadata_str = "".join([file.readline() for i in range(1)])
    metadata_str = metadata_str.replace("# ","")
    metadata = json.loads(metadata_str)

N_burn = metadata["N_burn"]
N_steps = metadata["N_steps"]

samples = np.genfromtxt("data/E3_3_samples.csv", delimiter=",")

import corner

plt.figure(figsize=(5,4))
corner.corner(samples[N_burn:], bins=50, labels=["$x_0$", "$x_1$", "$x_2$"])
#plt.hist2d(samples[:,0],samples[:,1], bins=100)
plt.savefig("plots/E3_3_samples.png")