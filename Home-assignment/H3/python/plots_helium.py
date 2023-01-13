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
E0_exact = -2.903385
print("Plotting Task 2...")
plot_dmc("data/task2.csv", "plots/task2.pdf", E0_exact)

################################################
# Task 3
################################################
print("Plotting Task 3...")
print("-First Order...")
plot_dmc("data/task3_1st_order.csv", "plots/task3_1st_order.pdf", E0_exact)
print("-Second Order...")
plot_dmc("data/task3_2nd_order.csv", "plots/task3_2nd_order.pdf", E0_exact)




################################################
# Task 4
################################################
print("Plotting Task 4...")
dtau, E_T_1st_order, E_T_std_1st_order, E_T_2nd_order, E_T_std_2nd_order, N0_1st_order, N0_2nd_order = np.genfromtxt("data/task4.csv", delimiter=",", unpack=True)

from scipy.optimize import curve_fit

def E_linear(dtau,a,b):
    return a*dtau + b

def E_quadratic(dtau, a, b):
    return a*dtau**2 + b 

(a_1st, b_1st), pcov = curve_fit(E_linear, dtau, E_T_1st_order)
(a_2nd, b_2nd), pcov = curve_fit(E_quadratic, dtau, E_T_2nd_order)

dtau_lin = np.linspace(0,dtau.max(), 1000)


fig, axs = plt.subplots(2,1,figsize=(5,8))
plt.sca(axs[0])
plt.axhline(E0_exact, color="k", linestyle="--", alpha=0.5, label=f"$E_0 = {E0_exact:.5f}$")

plt.plot(dtau_lin, E_linear(dtau_lin, a_1st, b_1st), "C0-", label=f"1st order fit\nE(0)={b_1st:.5f}")
plt.plot(dtau_lin, E_quadratic(dtau_lin, a_2nd, b_2nd), "C1-", label=f"2nd order fit\nE(0)={b_2nd:.5f}")

plt.errorbar(dtau, E_T_1st_order, yerr=E_T_std_1st_order, fmt="C0x", label="1st order")
plt.errorbar(dtau, E_T_2nd_order, yerr=E_T_std_2nd_order, fmt="C1x", label="2nd order")
plt.xlabel(r"$\Delta \tau \: / \:$(a.u.)")
plt.ylabel(r"$E_T \: / \:$(a.u.)")
plt.xlim(0,dtau.max())
plt.legend()

plt.sca(axs[1])
plt.plot(dtau, N0_1st_order, "x", label="1st order")
plt.plot(dtau, N0_2nd_order, "x", label="2nd order")
plt.xlabel(r"$\Delta \tau \: / \:$(a.u.)")
plt.ylabel(r"$N_0$")
plt.xlim(0,dtau.max())
plt.legend()

plt.tight_layout()
plt.savefig("plots/task4.pdf")
