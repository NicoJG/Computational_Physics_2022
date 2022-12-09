import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import binned_statistic

q_x, q_y, q_z, q, S_q = np.genfromtxt("data/H1_7.csv", delimiter=",", unpack=True)

df = pd.DataFrame(np.column_stack([q,S_q]), columns=["q","S_q"])
temp = df.groupby('q').mean().reset_index().values
q_mean = temp[:,0]
S_q_mean = temp[:,1]
temp = df.groupby('q').std().reset_index().values
S_q_std = temp[:,1]

q_min = 0.1
q_max = np.max(q)
n_bins = 300
bin_edges = np.linspace(q_min, q_max, n_bins+1)
bin_width = np.diff(bin_edges)[0]
bin_centers = bin_edges[:-1] + bin_width/2

S_q_binned_mean = binned_statistic(q,S_q,"mean",bin_edges)[0]
S_q_binned_std = binned_statistic(q,S_q,"std",bin_edges)[0]

plt.figure(figsize=(5,4))
#plt.plot(q[1:], S_q[1:], "x",alpha=0.3, label=r"$S(q)$")
plt.plot(q_mean[1:], S_q_mean[1:], "C1.", markersize=1., label="$<S(q)>$")
#plt.fill_between(q_mean,S_q_mean-S_q_std, S_q_mean+S_q_std, color="C0", alpha=0.2, label="uncertainty")
plt.fill_between(bin_edges[:-1],S_q_binned_mean-S_q_binned_std, S_q_binned_mean+S_q_binned_std, step="post", color="C0", alpha=0.2, label="binned std")
plt.step(bin_edges[:-1], S_q_binned_mean, where="post", label="binned mean")
plt.xlabel("$q \: / \: \mathrm{Å^-1}$")
plt.ylabel("$S(q)$")
plt.legend()
plt.tight_layout()
plt.savefig("plots/H1_7.pdf")
