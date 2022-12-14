# %%
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from tqdm.auto import tqdm

E_AA = -436e-3 # eV
E_BB = -113e-3 # eV
E_AB = -294e-3 # eV

k_B = 8.617333262e-5 # eV/K

N = 100000

E_0_per_N = 2*(E_AA+E_BB+2*E_AB)
dE = E_AA + E_BB - 2*E_AB

T_c = 2*dE/k_B
T_c_exp = 468 + 273.15

print(f"T_(c,mfa) = {T_c:.2f} K ; T_(c,exp) = {T_c_exp:.2f} K")

T = np.linspace(0,1200,1200) # K
#T = np.linspace(906,1200,1000) # K

def F_per_N(P,T):
    if P==1 or P==-1:
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

fig, axs = plt.subplots(2,3, figsize=(12,6))

#fig.suptitle(f"$T_c = {T_c:.3f} \, \mathrm{{K}}$")

plt.sca(axs[0][0])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.5)
plt.plot(T,P, color="C2")
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$P$")

plt.sca(axs[0][1])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.5)
plt.plot(T,E_per_N, color="C0")
#plt.plot(T,E_MFA_per_N)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$E/N \:/\: \mathrm{eV}$")

plt.sca(axs[0][2])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.5)
plt.plot(T,C_per_N, color="C1")
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$C/N \:/\: \mathrm{eV \, K^{-1}}$")

# log scale plots below
plt.sca(axs[1][0])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.5)
plt.plot(T,P, color="C2")
plt.yscale("log")
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$P$")

plt.sca(axs[1][1])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.5)
plt.plot(T,E_per_N[-1]-E_per_N, color="C0")
plt.yscale("log")
#plt.plot(T,E_MFA_per_N)
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(rf"$E(T={T[-1]:.0f}\mathrm{{K}})/N - E/N \:/\: \mathrm{{eV}}$")

plt.sca(axs[1][2])
plt.axvline(T_c, linestyle=":", color="k", alpha=0.5)
plt.axvline(T_c_exp, linestyle="--", color="k", alpha=0.5)
plt.plot(T,C_per_N, color="C1")
plt.yscale("log")
plt.xlabel(r"$T \:/\: \mathrm{K}$")
plt.ylabel(r"$C/N \:/\: \mathrm{eV \, K^{-1}}$")

plt.tight_layout()
plt.savefig("plots/task1.pdf")



# plot 3 different F(P) to visualize the minimization of F

def F_per_N(P,T):
    res = np.zeros_like(P)
    mask = (P==1)|(P==-1)
    res[mask] = - 2*k_B*T*np.log(2)
    P = P[~mask]
    res[~mask] = E_0_per_N - 2*P**2*dE - 2*k_B*T*np.log(2) + k_B*T*((1+P)*np.log(1+P)+(1-P)*np.log(1-P))
    return res

Ts = [100,500,750,10000] # K

P_lin = np.linspace(-1,1,1000)

fig, axs = plt.subplots(1,len(Ts),figsize=(len(Ts)*3,2.5))
for i,T_i in enumerate(Ts):
    F_ = F_per_N(P_lin, T_i)
    plt.sca(axs[i])
    plt.plot(P_lin, F_per_N(P_lin, T_i))
    plt.ylim(np.min(F_[F_<-1.5])*1.001, np.max(F_[F_<-1.5])*0.999)
    plt.xlabel(r"$P$")
    plt.ylabel(rf"$F (P, T={T_i:.0f}\mathrm{{K}}) \: / \: \mathrm{{eV}}$")
    
plt.tight_layout()
plt.savefig("plots/task1_F.pdf")
