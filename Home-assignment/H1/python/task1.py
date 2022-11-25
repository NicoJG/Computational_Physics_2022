import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def E_pot_fit(a0, c1, c2, c3):
    return c1*a0**2 + c2*a0 + c3

a0, E_pot = np.genfromtxt("data/H1_1.csv", delimiter=",", unpack=True)


# find index of nearest value to the cut point
cut_idx = len(a0)-1 #(np.abs(a0 - 4.03425425)).argmin()

params, pcov = curve_fit(E_pot_fit, a0[:cut_idx], E_pot[:cut_idx])

E_pot_fitted = E_pot_fit(a0, *params)

# minimum point:
a0_min = -0.5*(params[1]/params[0])

plt.figure(figsize=(5,4))
plt.axvline(a0_min, linestyle="--", color="k", alpha=0.5, label=rf"$\mathrm{{argmin}}_{{a_0}}(E_\mathrm{{pot}}) = {a0_min:.5f} \: \mathrm{{Å}}$")
plt.plot(a0, E_pot_fitted, "--", label="fit")
plt.plot(a0, E_pot, "-", label="simulation data")
plt.xlabel(r"$a_0 \: [\mathrm{Å}]$")
plt.ylabel(r"$E_\mathrm{pot} \: [\mathrm{eV}]$")
plt.legend()
plt.tight_layout()
plt.savefig("plots/H1_1.pdf")