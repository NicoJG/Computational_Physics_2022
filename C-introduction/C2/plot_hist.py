# %%
import numpy as np
import matplotlib.pyplot as plt

x = np.genfromtxt("random_numbers.csv", delimiter=",")
print("Read numpy shape: ",x.shape)
plt.hist(x, bins=100, histtype="stepfilled")
plt.ylabel("counts")
plt.xlabel("random number")
plt.tight_layout()
plt.savefig("random_numbers_hist.pdf")
# %%
