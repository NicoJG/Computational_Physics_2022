import numpy as np
import matplotlib.pyplot as plt
import json


def plot_dmc(data_file_path, plot_file_path, E0_exact):
    i_step, tau, N, E_T = np.genfromtxt(data_file_path, delimiter=",", unpack=True)
    
    # get header metadata
    with open(data_file_path, "r") as file:
        metadata_str = "".join([file.readline() for i in range(1)])
        metadata_str = metadata_str.replace("# ","")
        metadata = json.loads(metadata_str)
        
    n_eq_steps = metadata["n_eq_steps"]
    dtau = metadata["dtau"]
    
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
    plt.axvline(tau[n_eq_steps], color="k", alpha=0.5, linestyle="--")
    plt.axhline(E_mean, color="C1", linestyle=":", label=f"$<E_T> = {E_mean:.4f} \\pm {E_std:.4f}$")
    plt.axhline(E0_exact, color="k", alpha=0.5, linestyle=":", label=f"$E_0 = {E0_exact:.4f}$")
    plt.plot(tau, E_T, color="C1")
    plt.xlabel(r"$\tau \: / \:$ a.u.")
    plt.ylabel(r"$E_T \: / \: $ a.u.")
    plt.legend()

    plt.tight_layout()
    plt.savefig(plot_file_path)
    
    
    
    
################################################
# Task 2
################################################
E0_exact = -2.90372 # https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Quantum_Mechanics/10%3A_Multi-electron_Atoms/8%3A_The_Helium_Atom
print("Plotting Task 2...")
plot_dmc("data/task2.csv", "plots/task2.pdf", E0_exact)

################################################
# Task 3
################################################
print("Plotting Task 3...")
plot_dmc("data/task3_1st_order.csv", "plots/task3_1st_order.pdf", E0_exact)