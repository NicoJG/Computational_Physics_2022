# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root,minimize
from tqdm.auto import tqdm

E_AA = -436e-3 # eV
E_BB = -113e-3 # eV
E_AB = -294e-3 # eV

k_B = 8.617333262e-5 # eV/K

N = 1

E_0 = 2*N*(E_AA+E_BB+2*E_AB)
dE = E_AA + E_BB - 2*E_AB

T_c = 2*dE/k_B

T = np.linspace(0,1.5*T_c,1000) # K
#T = np.linspace(906,1200,1000) # K

def F(P,T):
    if P==1:
        return - 2*N*k_B*T*np.log(2)
    return E_0 - 2*N*P**2*dE - 2*N*k_B*T*np.log(2) + N*k_B*T*((1+P)*np.log(1+P)+(1-P)*np.log(1-P))

P = [minimize(F, 0.5, args=T_i, bounds=((0,1),)).x[0] for T_i in tqdm(T)]
P = np.array(P)

plt.figure(figsize=(5,4))
plt.axvline(T_c, linestyle=":", color="k", alpha=0.001, label=fr"$T_c = {T_c:.3f} \, \mathrm{{K}}$")
plt.plot(T,P)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$P$")
plt.legend()

# %%
N_AA = 2*(1-P**2)*N 
N_BB = 2*(1-P**2)*N
N_AB = 4*(1+P**2)*N

E = N_AA*E_AA + N_BB*E_BB + N_AB*E_AB

plt.figure(figsize=(5,4))
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5, label=fr"$T_c = {T_c:.3f} \, \mathrm{{K}}$")
plt.plot(T,E)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$E \:/\: \mathrm{eV}$")
plt.legend()
# %%
