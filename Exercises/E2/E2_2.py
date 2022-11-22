# %%
import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt("E2_2.csv", delimiter=",", comments="#")
print(data.shape)

vars = ["u","p","Q","P","E"]
t = data[:,0]
n_particles = (data.shape[1]-1)//len(vars)
n_timesteps = data.shape[0]

print(f"{n_particles = }")
print(f"{n_timesteps = }")

u = data[:,1:n_particles+1]
p = data[:,1+n_particles:1+2*n_particles]
Q = data[:,1+2*n_particles:1+3*n_particles]
P = data[:,1+3*n_particles:1+4*n_particles]
E = data[:,1+4*n_particles:1+5*n_particles]

print(p.shape)

plt.plot(t,Q)
plt.show()
plt.plot(t,P)
plt.show()
plt.plot(t,E)
plt.show()
# %%
