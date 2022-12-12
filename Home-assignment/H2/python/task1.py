# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root,minimize
from tqdm.auto import tqdm

E_AA = -436e-3 # eV
E_BB = -113e-3 # eV
E_AB = -294e-3 # eV

k_B = 8.617333262e-5 # eV/K

N = 100000

E_0_per_N = 2*(E_AA+E_BB+2*E_AB)
dE = E_AA + E_BB - 2*E_AB

T_c = 2*dE/k_B

T = np.linspace(0,1200,500) # K
#T = np.linspace(906,1200,1000) # K

def F_per_N(P,T):
    if P==1:
        return - 2*k_B*T*np.log(2)
    return E_0_per_N - 2*P**2*dE - 2*k_B*T*np.log(2) + k_B*T*((1+P)*np.log(1+P)+(1-P)*np.log(1-P))

P_0 = 0.5
P = []
for T_i in tqdm(T):
    P_0 = minimize(F_per_N, P_0, args=T_i, bounds=((0,1),)).x[0]
    P.append(P_0)
P = np.array(P)

N_AA_per_N = 2*(1-P**2) 
N_BB_per_N = 2*(1-P**2)
N_AB_per_N = 4*(1+P**2)

E_per_N = N_AA_per_N*E_AA + N_BB_per_N*E_BB + N_AB_per_N*E_AB
E_MFA_per_N = E_0_per_N - 2*P**2*dE

C_per_N = np.gradient(E_per_N, T)

fig, axs = plt.subplots(1,3, figsize=(15,4))

fig.suptitle(f"$T_c = {T_c:.3f} \, \mathrm{{K}}$")

plt.sca(axs[0])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.plot(T,P)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$P$")

plt.sca(axs[1])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.plot(T,E_per_N)
#plt.plot(T,E_MFA_per_N)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$E/N \:/\: \mathrm{eV}$")

plt.sca(axs[2])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.plot(T,C_per_N)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$C/N \:/\: \mathrm{eV \, K^{-1}}$")

plt.tight_layout()
plt.savefig("plots/task1.pdf")
