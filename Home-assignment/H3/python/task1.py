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
# plot the distribution during the timesteps
x_data = np.genfromtxt("data/task1_x.csv", delimiter=",")

i_step = x_data[:,0]
x_2d = x_data[:,1:]

i_step_2d = np.tile(i_step,(x_2d.shape[1],1)).T

mask = ~np.isnan(x_2d.flatten())
x = x_2d.flatten()[mask]
i_step = i_step_2d.flatten()[mask]

tau = i_step*dtau
# %%
fig, axs = plt.subplots(1,2, gridspec_kw={"width_ratios":[0.8,0.2]}, figsize=(10,4), sharey=True)
# plot the evolution of the distribution
plt.sca(axs[0])
plt.hist2d(tau, x, bins=(x_2d.shape[0], 40), rasterized=True)
plt.axvline(n_eq_steps*dtau, color="white", alpha=0.7, linestyle="--")
plt.xlabel(r"$\tau \: / \:$ a.u.")
plt.ylabel(r"$x \: / \:$ a.u.")

x_min = x.min()
x_max = x.max()

x_lin = np.linspace(x_min, x_max, 1000)
phi0 = np.sqrt(2) * np.exp(-np.exp(-x_lin)-x_lin/2)

mask2 = i_step>=n_eq_steps
weights = 1./N[i_step.astype(int)]

mean_hist, bin_edges = np.histogram(x[mask2], bins=100, weights=weights[mask2], density=True)
bin_centers = (bin_edges[:-1]+bin_edges[1:])/2

#mean_hist = 
plt.sca(axs[1])
plt.plot(phi0,x_lin)
plt.hist(bin_centers, bins=bin_edges, weights=mean_hist, orientation="horizontal")
plt.xlabel("$\Phi_0$")
plt.ylabel(r"$x \: / \:$ a.u.")
axs[1].yaxis.set_tick_params(labelright=True)
axs[1].yaxis.tick_right()
axs[1].yaxis.set_label_position("right")

plt.tight_layout()
plt.savefig("plots/task1_dist.pdf")


# %%
