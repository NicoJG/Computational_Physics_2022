# %%
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json


# plot the results
print(f"Plot final results...")
T, E, E_std, s_E, P, P_std, s_P, r, r_std, s_r, C = np.genfromtxt("data/H2a.csv", delimiter=",", unpack=True)

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
T_c_exp = 468 + 273.15

E_err = np.sqrt(s_E*E_std**2/n_steps)
P_err = np.sqrt(s_P*P_std**2/n_steps)
r_err = np.sqrt(s_r*r_std**2/n_steps)

fig, axs = plt.subplots(4,1, figsize=(10,10))

#fig.suptitle(f"$T_{{c,task1}} = {T_c:.3f} \, \mathrm{{K}}$, $N_\mathrm{{eq}} = {n_eq_steps:.1e}$, $N = {n_steps:.1e}$")


plt.sca(axs[0])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.9)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.9)
plt.fill_between(T, E-E_err, E+E_err, color="C0", alpha=0.3)
plt.plot(T,E, color="C0")
plt.ylim(-2360,-2280)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$E \:/\: \mathrm{eV}$")
plt.grid()

plt.sca(axs[1])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.9)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.9)
plt.plot(T,C, color="C1")
plt.ylim(0,0.5)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$C \:/\: \mathrm{eV \, K^{-1}}$")
plt.grid()

plt.sca(axs[2])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.9)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.9)
plt.fill_between(T, np.abs(P)-P_err, np.abs(P)+P_err, color="C2", alpha=0.3)
plt.plot(T,P, color="green", linestyle=":", alpha=0.5, label="$P$")
plt.plot(T,np.abs(P), color="C2", label="$|P|$")
plt.ylim(-0.5,1)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$P$")
plt.grid()
plt.legend()

plt.sca(axs[3])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.9)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.9)
plt.fill_between(T, r-r_err, r+r_err, color="C3", alpha=0.3)
plt.plot(T,r, color="C3")
plt.ylim(0,1)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$r$")
plt.grid()

plt.tight_layout()
plt.savefig("plots/task2.pdf")

# plot the uncertainties
fig, axs = plt.subplots(3,1, figsize=(10,7.5))

plt.sca(axs[0])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.9)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.9)
plt.plot(T,E_err, color="C0")
plt.yscale("log")
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$\sigma_{\langle E \rangle} \:/\: \mathrm{eV}$")
plt.grid()

plt.sca(axs[1])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.9)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.9)
plt.plot(T,P_err, color="C2")
plt.yscale("log")
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$\sigma_{\langle P \rangle}$")
plt.grid()

plt.sca(axs[2])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.9)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.9)
plt.plot(T,r_err, color="C3")
plt.yscale("log")
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$\sigma_{\langle r \rangle}$")
plt.grid()

plt.tight_layout()
plt.savefig("plots/task2_uncertainties.pdf")


# plot the statistical inefficiency
plt.figure(figsize=(4,3))
plt.plot(T,s_E, "C0", label="$s(E)$ and \n$s(r)$")
plt.plot(T,s_P, "C2", label="$s(P)$")
#plt.plot(T,s_r, "C3", linestyle=":", label="$s(r)$")
plt.yscale("log")
plt.xlabel(r"$T \: / \: \mathrm{K}$")
plt.ylabel(r"$s$ (statistical inefficiency)")
plt.legend()
plt.tight_layout()
plt.savefig("plots/H2a_stat_ineff.pdf")


# Plot a few full simulations
full_simulations_dir = Path("data/full_simulations/")
fig, axs = plt.subplots(3,1, figsize=(10,8))

#fig.suptitle(f"$T = {T:2f} \,\mathrm{{K}}$, $N_\mathrm{{eq}} = {n_eq_steps:.1e}$, \n$<P> = {P_avg:.4f} \mathrm{{eV}}$, $<E> = {E_avg:.4f} \mathrm{{eV}}$, $<C> = {C:.4f} \mathrm{{eV/K}}$, $<r> = {r_avg:.4f}$")
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
    
    print(f"Plot full simulation T={T:.0f} K ...")
    
    i_step, E, P, r = np.genfromtxt(full_simulation_file, delimiter=",", unpack=True)
    
    plt.sca(axs[0])
    plt.plot(i_step, P, label=f"$T = {T:.0f} \\mathrm{{K}}$", rasterized=True)
    
    plt.sca(axs[1])
    plt.plot(i_step, E, label=f"$T = {T:.0f} \\mathrm{{K}}$", rasterized=True)
    
    plt.sca(axs[2])
    plt.plot(i_step, r, label=f"$T = {T:.0f} \\mathrm{{K}}$", rasterized=True)
    

plt.sca(axs[0])
plt.axvline(n_eq_steps, linestyle="--", color="k", alpha=0.5, label=rf"$N_\mathrm{{eq}} = {n_eq_steps:.1e}$")
plt.xlabel("simulation step")
plt.ylabel("P")
plt.legend(loc="center right")

plt.sca(axs[1])
plt.axvline(n_eq_steps, linestyle="--", color="k", alpha=0.5, label=rf"$N_\mathrm{{eq}} = {n_eq_steps:.1e}$")
plt.xlabel("simulation step")
plt.ylabel("energy / eV")
plt.legend(loc="center right")

plt.sca(axs[2])
plt.axvline(n_eq_steps, linestyle="--", color="k", alpha=0.5, label=rf"$N_\mathrm{{eq}} = {n_eq_steps:.1e}$")
plt.xlabel("simulation step")
plt.ylabel("r")
plt.legend(loc="center right")
    
plt.tight_layout()
plt.savefig(f"plots/H2a_simsteps.pdf")
    

# Plot a few full simulation correlation functions
corr_simulations_dir = Path("data/corr_simulations/")
plt.figure(figsize=(4,3))

for i,corr_simulation_file in enumerate(corr_simulations_dir.iterdir()):
    # get header metadata
    with open(corr_simulation_file, "r") as file:
        metadata_str = "".join([file.readline() for i in range(1)])
        metadata_str = metadata_str.replace("# ","")
        metadata = json.loads(metadata_str)
        
    T = metadata["T[K]"]
    s_E = metadata["s_E"]
    
    print(f"Plot correlation function T={T:.0f} K ...")
    
    k, Phi = np.genfromtxt(corr_simulation_file, delimiter=",", unpack=True)
    
    plt.plot(k, Phi, color=f"C{i}", label=f"$T = {T:.0f} \\mathrm{{K}}$", rasterized=True)
    plt.axvline(s_E, linestyle="--", color=f"C{i}", alpha=0.5)
    
plt.axhline(np.exp(-2), linestyle=":", color="k", alpha=0.5, label=f"$\exp(-2)$")
plt.xscale("log")
#plt.yscale("log")
plt.xlabel("$k$")
plt.ylabel(r"$\Phi_E(k)$")
plt.legend(loc="lower left")
    
plt.tight_layout()
plt.savefig(f"plots/H2a_corr_func.pdf")

# %%
