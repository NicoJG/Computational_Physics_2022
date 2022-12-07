import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json


# Plot a few full simulations
full_simulations_dir = Path("data/full_simulations/")
for full_simulation_file in full_simulations_dir.iterdir():
    # get header metadata
    with open(full_simulation_file, "r") as file:
        metadata_str = "".join([file.readline() for i in range(2)])
        metadata_str = metadata_str.replace("# ","")
        metadata = json.loads(metadata_str)
        
    T = metadata["T[K]"]
    n_eq_steps = metadata["n_eq_steps"]
    
    E_avg = metadata["E[eV]"]
    P_avg = metadata["P"]
    r_avg = metadata["r"]
    C = metadata["C[eV/K]"]
    
    print(f"Plot full simulation T={T:.2f} K ...")
    
    i_step, E, P, r = np.genfromtxt(full_simulation_file, delimiter=",", unpack=True)
    fig, axs = plt.subplots(1,3, figsize=(15,4))
    
    fig.suptitle(f"$T = {T:2f} \,\mathrm{{K}}$, $N_\mathrm{{eq}} = {n_eq_steps:.1e}$, \n$<P> = {P_avg:.4f} \mathrm{{eV}}$, $<E> = {E_avg:.4f} \mathrm{{eV}}$, $<C> = {C:.4f} \mathrm{{eV/K}}$, $<r> = {r_avg:.4f}$")
    
    plt.sca(axs[0])
    plt.plot(i_step, P)
    plt.axvline(n_eq_steps, linestyle="--", color="k", alpha=0.5)
    plt.xlabel("simulation step")
    plt.ylabel("P")
    
    plt.sca(axs[1])
    plt.plot(i_step, E, rasterized=True)
    plt.axvline(n_eq_steps, linestyle="--", color="k", alpha=0.5)
    plt.xlabel("simulation step")
    plt.ylabel("energy / eV")
    
    plt.sca(axs[2])
    plt.plot(i_step, r)
    plt.axvline(n_eq_steps, linestyle="--", color="k", alpha=0.5, label=rf"$N_\mathrm{{eq}} = {n_eq_steps:.1e}$")
    plt.xlabel("simulation step")
    plt.ylabel("r")
    
    plt.tight_layout()
    plt.savefig(f"plots/H2a_simsteps_T{T:.2f}K.pdf")
    
    
# plot the results
print(f"Plot final results...")
T, E, P, r, C = np.genfromtxt("data/H2a.csv", delimiter=",", unpack=True)

# get header metadata
with open("data/H2a.csv", "r") as file:
    metadata_str = "".join([file.readline() for i in range(1)])
    metadata_str = metadata_str.replace("# ","")
    metadata = json.loads(metadata_str)
    
n_eq_steps = metadata["n_eq_steps"]
n_steps = metadata["n_steps"]

# calculate T_c
k_B = 8.617333262e-5 # eV/K
E_AA = -436e-3 # eV
E_BB = -113e-3 # eV
E_AB = -294e-3 # eV
dE = E_AA + E_BB - 2*E_AB
T_c = 2*dE/k_B

fig, axs = plt.subplots(1,3, figsize=(15,4))

fig.suptitle(f"$T_c = {T_c:.3f} \, \mathrm{{K}}$, $N_\mathrm{{eq}} = {n_eq_steps:.1e}$, $N = {n_steps:.1e}$")

plt.sca(axs[0])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.plot(T,P)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$P$")

plt.sca(axs[1])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.plot(T,E)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$E \:/\: \mathrm{eV}$")

plt.sca(axs[2])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.plot(T,C)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$C \:/\: \mathrm{eV \, K^{-1}}$")

plt.tight_layout()
plt.savefig("plots/task2.pdf")

# plot r
plt.figure(figsize=(5,4))
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.plot(T,r)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$r$")
plt.tight_layout()
plt.savefig("plots/task2_r.pdf")
