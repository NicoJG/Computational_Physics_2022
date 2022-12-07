import numpy as np
import matplotlib.pyplot as plt
import json

def plot_energies_task2(input_file, output_file):
    # get header metadata
    with open(input_file, "r") as file:
        metadata_str = "".join([file.readline() for i in range(4)])
        metadata_str = metadata_str.replace("# ","")
        metadata = json.loads(metadata_str)
        
    dt = metadata["dt"]

    t, E_pot, E_kin, T, P, a0 = np.genfromtxt(input_file, delimiter=",", unpack=True)

    T_mean = np.mean(T)

    fig, axs = plt.subplots(1,3,figsize=(12,3))
    fig.suptitle(f"T = {T_mean:.5f} K ; dt = {dt*1e3:.1f} fs")
    
    plt.sca(axs[0])
    plt.plot(t, E_pot, "C0", label=r"$E_\mathrm{pot}$")
    plt.xlabel("time [ps]")
    plt.ylabel("potential energy [eV]")
    plt.legend()
    
    plt.sca(axs[1])
    plt.plot(t, E_kin, "C1", label=r"$E_\mathrm{kin}$")
    plt.xlabel("time [ps]")
    plt.ylabel("kinetic energy [eV]")
    plt.legend()
    
    plt.sca(axs[2])
    plt.plot(t, E_pot+E_kin, "C2", label=r"$E_\mathrm{tot}$")
    plt.xlabel("time [ps]")
    plt.ylabel("total energy [eV]")
    plt.legend()
    
    
    plt.tight_layout()
    plt.savefig(output_file)

plot_energies_task2("data/H1_2_small_enough.csv", "plots/H1_2_small_enough.pdf")
plot_energies_task2("data/H1_2_far_too_large.csv", "plots/H1_2_far_too_large.pdf")
plot_energies_task2("data/H1_2_little_too_large.csv", "plots/H1_2_little_too_large.pdf")
