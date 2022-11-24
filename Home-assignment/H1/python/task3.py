import numpy as np
import matplotlib.pyplot as plt
import json
from task2 import plot_energies_task2


t, E_pot, E_kin, T, P, a0 = np.genfromtxt("data/H1_3_temp_scaling.csv", delimiter=",", unpack=True)
t2, E_pot2, E_kin2, T2, P2, a02 = np.genfromtxt("data/H1_3_pressure_scaling.csv", delimiter=",", unpack=True)


# get header metadata
with open("data/H1_3_pressure_scaling.csv", "r") as file:
    metadata_str = "".join([file.readline() for i in range(3)])
    metadata_str = metadata_str.replace("# ","")
    metadata = json.loads(metadata_str)
    
T_desired = metadata["T_desired"]
P_desired = metadata["P_desired"]

t_switch_on_P_scaling = (t[-1]+t2[0])/2

fig, axs = plt.subplots(2,1,figsize=(5,6))
plt.sca(axs[0])
plt.axvline(t_switch_on_P_scaling, linestyle="--", color="k", alpha=0.5)
plt.plot(t, T, "C0")
plt.plot(t2, T2, "C0")
plt.axhline(T_desired, linestyle=":", color="k")
plt.xlabel("time [ps]")
plt.ylabel("temperature [K]")

plt.sca(axs[1])
plt.axvline(t_switch_on_P_scaling, linestyle="--", color="k", alpha=0.5)
plt.plot(t, P, "C1")
plt.plot(t2, P2, "C1")
plt.axhline(P_desired, linestyle=":", color="k")
plt.xlabel("time [ps]")
plt.ylabel("pressure [eV Å^-3]")

plt.tight_layout()
plt.savefig("plots/H1_3_temp_scaling.pdf")
