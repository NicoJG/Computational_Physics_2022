import numpy as np
import matplotlib.pyplot as plt
import json

def plot_energies_task2(input_file, output_file):
    # get header metadata
    with open(input_file, "r") as file:
        metadata_str = file.readline()[1:].strip()
        metadata = json.loads(metadata_str)
        
    dt = metadata["dt"]

    t, E_pot, E_kin = np.genfromtxt(input_file, delimiter=",", unpack=True)

    k_B = 8.617333262e-5 # eV/K (https://en.wikipedia.org/wiki/Boltzmann_constant#Value_in_different_units)
    N = 256

    T = 2/(3*N*k_B) * np.mean(E_kin)

    plt.figure(figsize=(5,4))
    plt.title(f"T = {T:.5f} K ; dt = {dt*1e3:.4f} fs")
    plt.plot(t, E_pot, label=r"$E_\mathrm{pot}$")
    plt.plot(t, E_kin, label=r"$E_\mathrm{kin}$")
    plt.plot(t, E_pot+E_kin, label=r"$E_\mathrm{tot}$")

    plt.xlabel("time [ps]")
    plt.ylabel("energy [eV]")

    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file)

plot_energies_task2("data/H1_2_small_enough.csv", "plots/H1_2_small_enough.pdf")
plot_energies_task2("data/H1_2_far_too_large.csv", "plots/H1_2_far_too_large.pdf")
plot_energies_task2("data/H1_2_little_too_large.csv", "plots/H1_2_little_too_large.pdf")
