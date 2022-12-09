import numpy as np
import matplotlib.pyplot as plt
import json

# get header metadata
with open("data/H1_6.csv", "r") as file:
    metadata_str = "".join([file.readline() for i in range(1)])
    metadata_str = metadata_str.replace("# ","")
    metadata = json.loads(metadata_str)
    
n_atoms = metadata["n_atoms"]
L_box = metadata["L_box[Å]"]
V = L_box**3

r, N_r = np.genfromtxt("data/H1_6.csv", delimiter=",", unpack=True)

n_bins = len(r)
dr = np.diff(r)[0]
r_max = r[-1] + dr/2

print(f"n_bins = {n_bins} ; r_max = {r_max:.4f} Å ; dr = {dr:.4f} Å")

N_ideal = [(n_atoms-1)*4*np.pi*(3*k**2-3*k+1)*dr**3/(3*V) for k in range(1,n_bins+1)]
N_ideal = np.array(N_ideal)

g_r = N_r/N_ideal

# find the first local minimum
r_m_idx = np.argmin(g_r[np.argmax(g_r):])+np.argmax(g_r)
r_m = r[r_m_idx]

I_r_m = 4*np.pi*(n_atoms/V)*np.trapz(g_r[:r_m_idx]*r[:r_m_idx]**2 ,r[:r_m_idx])

plt.figure(figsize=(5,4))
plt.axvline(r_m, linestyle="--", color="k", alpha=0.5, label=f"$r_m = {r_m:.4f} \\: \\mathrm{{Å}}$")
plt.fill_between(r[:r_m_idx],0,g_r[:r_m_idx], color="C1", alpha=0., label=f"$\\mathcal{{I}}(r_m) = {I_r_m:.4f}$")
plt.plot(r,g_r)
plt.xlabel(r"$r \: / \: \mathrm{Å}$")
plt.ylabel(r"$g(r)$")
plt.legend()
plt.tight_layout()
plt.savefig("plots/H1_6.pdf")