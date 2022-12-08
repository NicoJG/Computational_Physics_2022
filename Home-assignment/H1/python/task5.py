import numpy as np
import matplotlib.pyplot as plt
import json

eV = 1.602176634e-19 # J
celsius = 273.15 # K
k_B = 8.617333262e-5 # eV/K
m_u = 1.660539066e-27 # kg

m_Al = 26.98153853 * m_u * 10**3 # g

def C_V(E,N,T):
    E_var = np.var(E)
    return 3*N*k_B/2 / (1 - 2*E_var/(3*N*k_B**2*T**2))

def C_V_specific(E,N,T):
    C_V_ = C_V(E,N,T) # eV/K
    return C_V_ * eV / (N*m_Al) # J/(g K)
    

# solid state
data = np.genfromtxt("data/H1_3_after_scaling.csv", delimiter=",",unpack=True)

with open("data/H1_3_pressure_scaling.csv", "r") as file:
    metadata_str = "".join([file.readline() for i in range(4)])
    metadata_str = metadata_str.replace("# ","")
    metadata = json.loads(metadata_str)
    
N = metadata["n_atoms"]
T_desired = metadata["T_desired"]

t, E_pot, E_kin, T, P, a0 = data[:6]
E_tot = E_kin + E_pot
pos = data[6:]
x = pos[::3]
y = pos[1::3]
z = pos[2::3]

T_mean = np.mean(T)
T_std = np.std(T)

print("Task 5 - Solid state:")
print(f"T_desired = {T_desired:.2f} K ; T_mean = {T_mean:.2f} +- {T_std:.2f} K")
print(f"T_desired = {T_desired-celsius:.2f} °C ; T_mean = {T_mean-celsius:.2f} +- {T_std:.2f} °C")
print(f"C_V (E_kin, T_mean)    = {C_V(E_kin,N,T_mean):.6f} eV/K")
print(f"C_V (E_pot, T_mean)    = {C_V(E_pot,N,T_mean):.6f} eV/K")
print(f"c_V (E_kin, T_mean)    = {C_V_specific(E_kin,N,T_mean):.6f} J/(g K)")
print(f"c_V (E_pot, T_mean)    = {C_V_specific(E_pot,N,T_mean):.6f} J/(g K)")
print(f"C_V (E_kin, T_desired) = {C_V(E_kin,N,T_desired):.6f} eV/K")
print(f"C_V (E_pot, T_desired) = {C_V(E_pot,N,T_desired):.6f} eV/K")
print(f"c_V (E_kin, T_desired) = {C_V_specific(E_kin,N,T_desired):.6f} J/(g K)")
print(f"c_V (E_pot, T_desired) = {C_V_specific(E_pot,N,T_desired):.6f} J/(g K)")


# fluid state
data = np.genfromtxt("data/H1_4_after_scaling.csv", delimiter=",",unpack=True)

with open("data/H1_4_temp_decreasing.csv", "r") as file:
    metadata_str = "".join([file.readline() for i in range(4)])
    metadata_str = metadata_str.replace("# ","")
    metadata = json.loads(metadata_str)
    
N = metadata["n_atoms"]
T_desired = metadata["T_desired"]

t, E_pot, E_kin, T, P, a0 = data[:6]
E_tot = E_kin + E_pot
pos = data[6:]
x = pos[::3]
y = pos[1::3]
z = pos[2::3]

T_mean = np.mean(T)
T_std = np.std(T)

print("Task 5 - Liquid state:")
print(f"T_desired = {T_desired:.2f} K ; T_mean = {T_mean:.2f} +- {T_std:.2f} K")
print(f"T_desired = {T_desired-celsius:.2f} °C ; T_mean = {T_mean-celsius:.2f} +- {T_std:.2f} °C")
print(f"C_V (E_kin, T_mean)    = {C_V(E_kin,N,T_mean):.6f} eV/K")
print(f"C_V (E_pot, T_mean)    = {C_V(E_pot,N,T_mean):.6f} eV/K")
print(f"c_V (E_kin, T_mean)    = {C_V_specific(E_kin,N,T_mean):.6f} J/(g K)")
print(f"c_V (E_pot, T_mean)    = {C_V_specific(E_pot,N,T_mean):.6f} J/(g K)")
print(f"C_V (E_kin, T_desired) = {C_V(E_kin,N,T_desired):.6f} eV/K")
print(f"C_V (E_pot, T_desired) = {C_V(E_pot,N,T_desired):.6f} eV/K")
print(f"c_V (E_kin, T_desired) = {C_V_specific(E_kin,N,T_desired):.6f} J/(g K)")
print(f"c_V (E_pot, T_desired) = {C_V_specific(E_pot,N,T_desired):.6f} J/(g K)")

