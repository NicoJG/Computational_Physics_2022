# %%
import numpy as np
import matplotlib.pyplot as plt
import json

i_step, tau, N, E_T = np.genfromtxt("data/task1.csv", delimiter=",", unpack=True)

# get header metadata
with open("data/task1.csv", "r") as file:
    metadata_str = "".join([file.readline() for i in range(1)])
    metadata_str = metadata_str.replace("# ","")
    metadata = json.loads(metadata_str)
    
n_eq_steps = metadata["n_eq_steps"]
dtau = metadata["dtau"]

E0_exact = 3/8

E_mean = E_T[n_eq_steps:].mean()
E_std = E_T[n_eq_steps:].std()

N_mean = N[n_eq_steps:].mean()
N_std = N[n_eq_steps:].std()

fig, axs = plt.subplots(2,1,figsize=(8,8))

plt.sca(axs[0])
plt.axvline(tau[n_eq_steps], color="k", alpha=0.5, linestyle="--")
plt.axhline(N_mean, color="C0", linestyle=":", label=f"$<N> = {N_mean:.0f} \\pm {N_std:.1f}$")
plt.plot(tau, N, color="C0")
plt.xlabel(r"$\tau \: / \:$ a.u.")
plt.ylabel(r"$N \: / \:$ number of walkers")
plt.legend()

plt.sca(axs[1])
plt.axvline(tau[n_eq_steps], color="k", linestyle="--")
plt.axhline(E_mean, color="C1", linestyle=":", label=f"$<E_T> = {E_mean:.4f} \\pm {E_std:.4f}$")
plt.axhline(E0_exact, color="k", alpha=0.5, linestyle=":", label=f"$E_0 = {E0_exact:.4f}$")
plt.plot(tau, E_T, color="C1")
plt.xlabel(r"$\tau \: / \:$ a.u.")
plt.ylabel(r"$E_T \: / \: $ a.u.")
plt.legend()

plt.tight_layout()
plt.savefig("plots/task1.pdf")

# %%
# plot the histogram

x_left, x_right, x_center, bin_density = np.genfromtxt("data/task1_x_hist.csv", delimiter=",", unpack=True)

x_min = x_left[0]
x_max = x_right[-1]

x_lin = np.linspace(x_min, x_max, 1000)
phi0 = np.sqrt(2) * np.exp(-np.exp(-x_lin)-x_lin/2)

plt.figure(figsize=(5,4))
plt.plot(x_center, bin_density, label=f"histogram (integral={np.trapz(bin_density,x_center):.4f})")
plt.plot(x_lin, phi0**2, label=f"exact (integral={np.trapz(phi0**2,x_lin):.4f})")
plt.ylabel("$\Phi_0^2$")
plt.xlabel(r"$x \: / \:$ a.u.")
plt.legend()
plt.tight_layout()
plt.savefig("plots/task1_x_hist.pdf")
