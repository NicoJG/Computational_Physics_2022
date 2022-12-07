import numpy as np
import matplotlib.pyplot as plt
import json
from task2 import plot_energies_task2


t, E_pot, E_kin, T, P, a0 = np.genfromtxt("data/H1_3_temp_scaling.csv", delimiter=",", unpack=True)
t2, E_pot2, E_kin2, T2, P2, a02 = np.genfromtxt("data/H1_3_pressure_scaling.csv", delimiter=",", unpack=True)

E_tot = E_kin + E_pot
E_tot2 = E_kin2 + E_pot2

# get header metadata
with open("data/H1_3_pressure_scaling.csv", "r") as file:
    metadata_str = "".join([file.readline() for i in range(4)])
    metadata_str = metadata_str.replace("# ","")
    metadata = json.loads(metadata_str)
    
T_desired = metadata["T_desired"]
P_desired = metadata["P_desired"]

t_switch_on_P_scaling = (t[-1]+t2[0])/2

fig, axs = plt.subplots(4,1,figsize=(8,8))

# plot the temperature
plt.sca(axs[0])
plt.axvline(t_switch_on_P_scaling, linestyle="--", color="k", alpha=0.5)#, label=f"start pressure \nequilibration \nat $t={t_switch_on_P_scaling:.2f} \\:\mathrm{{ps}}$")
plt.plot(t, T, "C0")
plt.plot(t2, T2, "C0")
plt.axhline(T_desired, linestyle=":", color="k", label=f"desired $T={T_desired:.2f} \\: \\mathrm{{K}}$")
plt.xlabel("time [ps]")
plt.ylabel("temperature [K]")
plt.legend(loc="upper right")

# plot the pressure
plt.sca(axs[1])
plt.axvline(t_switch_on_P_scaling, linestyle="--", color="k", alpha=0.5)#, label=f"start pressure \nequilibration \nat $t={t_switch_on_P_scaling:.2f} \\:\mathrm{{ps}}$")
plt.plot(t, P, "C1")
plt.plot(t2, P2, "C1")
plt.axhline(P_desired, linestyle=":", color="k", label=f"desired $P={P_desired:.2f} \\: \\mathrm{{bar}}$")
plt.xlabel("time [ps]")
plt.ylabel("pressure [bar]")
plt.legend(loc="upper right")

# plot the total energy
plt.sca(axs[2])
plt.axvline(t_switch_on_P_scaling, linestyle="--", color="k", alpha=0.5)#, label=f"start pressure \nequilibration \nat $t={t_switch_on_P_scaling:.2f} \\:\mathrm{{ps}}$")
plt.plot(t, E_tot, "C2")
plt.plot(t2, E_tot2, "C2")
plt.axhline(E_tot2[-1], linestyle=":", color="k", label=f"end $E_\mathrm{{tot}}={E_tot2[-1]:.2f} \\: \\mathrm{{eV}}$")
plt.xlabel("time [ps]")
plt.ylabel("total energy [eV]")
plt.legend(loc="lower right")

# plot the lattice constant
plt.sca(axs[3])
plt.axvline(t_switch_on_P_scaling, linestyle="--", color="k", alpha=0.5)#, label=f"start pressure \nequilibration \nat $t={t_switch_on_P_scaling:.2f} \\:\mathrm{{ps}}$")
plt.plot(t, a0, "C3")
plt.plot(t2, a02, "C3")
plt.axhline(a02[-1], linestyle=":", color="k", label=f"end $a_0={a02[-1]:.4f} \\: \\mathrm{{Å}}$")
plt.xlabel("time [ps]")
plt.ylabel("lattice constant [Å]")
plt.legend(loc="lower right")

plt.tight_layout()
plt.savefig("plots/H1_3_equilibration.pdf")



##########################################################
# Plot the simulation after the equilibration


# get header metadata
with open("data/H1_3_after_scaling.csv", "r") as file:
    metadata_str = "".join([file.readline() for i in range(4)])
    metadata_str = metadata_str.replace("# ","")
    metadata = json.loads(metadata_str)
    
data = np.genfromtxt("data/H1_3_after_scaling.csv", delimiter=",",unpack=True)

t, E_pot, E_kin, T, P, a0 = data[:6]
E_tot = E_kin + E_pot
pos = data[6:]
x = pos[::3]
y = pos[1::3]
z = pos[2::3]

# plot the positions
fig, axs = plt.subplots(3,1, figsize=(8,6))

plt.sca(axs[0])
for i,x_i in enumerate(x):
    plt.plot(t, x_i[0]-x_i, label=f"x_{i}")
    
plt.xlabel("time [ps]")
plt.ylabel("$x(t=0) - x(t)$ [Å]")
plt.legend(loc="upper right")

plt.sca(axs[1])
for i,y_i in enumerate(y):
    plt.plot(t, y_i[0]-y_i, label=f"y_{i}")
    
plt.xlabel("time [ps]")
plt.ylabel("$y(t=0) - y(t)$ [Å]")
plt.legend(loc="upper right")

plt.sca(axs[2])
for i,z_i in enumerate(z):
    plt.plot(t, z_i[0]-z_i, label=f"z_{i}")
    
plt.xlabel("time [ps]")
plt.ylabel("$z(t=0) - z(t)$ [Å]")
plt.legend(loc="upper right")

plt.tight_layout()
plt.savefig("plots/H1_3_positions.pdf")



############################################################
# plot the equilibration plots for the simulation without scaling
fig, axs = plt.subplots(4,1,figsize=(8,8))

# plot the temperature
plt.sca(axs[0])
plt.plot(t, T, "C0")
plt.axhline(T_desired, linestyle=":", alpha=0.5, color="k", label=f"desired $T={T_desired:.2f} \\: \\mathrm{{K}}$")
plt.axhline(np.mean(T), linestyle="--", alpha=0.5, color="k", label=f"mean $T={np.mean(T):.2f} \\: \\mathrm{{K}}$")
plt.plot(t[0],T[0],alpha=0., label=f"std $T={np.std(T):.2f} \\: \\mathrm{{K}}$")
plt.xlabel("time [ps]")
plt.ylabel("temperature [K]")
plt.legend(loc="upper right")

# plot the pressure
plt.sca(axs[1])
plt.plot(t, P, "C1")
plt.axhline(P_desired, linestyle=":", alpha=0.5, color="k", label=f"desired $P={P_desired:.2f} \\: \\mathrm{{bar}}$")
plt.axhline(np.mean(P), linestyle="--", alpha=0.5, color="k", label=f"mean $P={np.mean(P):.2f} \\: \\mathrm{{bar}}$")
plt.plot(t[0],P[0],alpha=0., label=f"std $P={np.std(P):.2f} \\: \\mathrm{{bar}}$")
plt.xlabel("time [ps]")
plt.ylabel("pressure [bar]")
plt.legend(loc="upper right")

# plot the total energy
plt.sca(axs[2])
plt.plot(t, E_tot, "C2")
plt.axhline(E_tot2[-1], linestyle=":", alpha=0.5, color="k", label=f"end $E_\mathrm{{tot}}={E_tot2[-1]:.2f} \\: \\mathrm{{eV}}$")
plt.axhline(np.mean(E_tot), linestyle="--", alpha=0.5, color="k", label=f"mean $E_\mathrm{{tot}}={np.mean(E_tot):.2f} \\: \\mathrm{{eV}}$")
plt.plot(t[0],E_tot[0],alpha=0., label=f"std $E_\mathrm{{tot}}={np.std(E_tot):.4f} \\: \\mathrm{{eV}}$")
plt.xlabel("time [ps]")
plt.ylabel("total energy [eV]")
plt.legend(loc="lower right")

# plot the lattice constant
plt.sca(axs[3])
plt.plot(t, a0, "C3")
plt.axhline(a02[-1], linestyle=":", alpha=0.5, color="k", label=f"end $a_0={a02[-1]:.4f} \\: \\mathrm{{Å}}$")
plt.axhline(np.mean(a0), linestyle="--", alpha=0.5, color="k", label=f"mean $a_0={np.mean(a0):.4f} \\: \\mathrm{{Å}}$")
plt.xlabel("time [ps]")
plt.ylabel("lattice constant [Å]")
plt.legend(loc="lower right")

plt.tight_layout()
plt.savefig("plots/H1_3_after_equilibration.pdf")