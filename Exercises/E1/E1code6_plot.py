#!/usr/bin/env python
###############################################################################
# E1code6_plot
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
data = np.genfromtxt('E1_6.csv', delimiter=',', comments="#", unpack=True)
t, q_1, q_2, q_3, v_1, v_2, v_3, E_kin, E_pot = data
E_tot = E_kin + E_pot

# powerspectrum
data2 = np.genfromtxt('E1_6_powerspectrum.csv', delimiter=',', comments="#", unpack=True)
f, p_1, p_2, p_3 = data2

fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(8,12))

# plot the displacement over time
ax1.plot(t, q_1, "-", label="q_1")
ax1.plot(t, q_2, "-", label="q_2")
ax1.plot(t, q_3, "-", label="q_3")
ax1.set_xlim(t.min(), t.max())
ax1.set_xlabel("time [ps]")
ax1.set_ylabel("displacement [Å]")
ax1.grid()
ax1.legend()

ax2.plot(t, E_kin, "-", label="E_kin")
ax2.plot(t, E_pot, "-", label="E_pot")
ax2.plot(t, E_tot, "-", label="E_tot")
ax2.set_xlim(t.min(), t.max())
ax2.set_xlabel("time [ps]")
ax2.set_ylabel("energy [eV]")
ax2.grid()
ax2.legend()

ax3.plot(f, p_1, "-", label="p_1")
ax3.plot(f, p_2, "-", label="p_2")
ax3.plot(f, p_3, "-", label="p_3")
ax3.set_xlim(-100, 100)
ax3.set_xlabel("frequency [1/ps]")
ax3.set_ylabel("powerspectrum [$Å^{-2}$]")
ax3.grid()
ax3.legend()

print(f"resonance frequencies p1: {np.unique(f[p_1>0.005])} THz")
print(f"resonance frequencies p2: {np.unique(f[p_2>0.005])} THz")
print(f"resonance frequencies p3: {np.unique(f[p_3>0.005])} THz")

c = 299792458 # m/s
for nu in [1388, 2349, 667]:
    print(f"experimental nu={nu:.4e} 1/cm ; f={nu*10**2 * c / 10**12:.5e} THz")

plt.tight_layout()
fig.savefig('E1_6.pdf')
fig.savefig('E1_6.png')
