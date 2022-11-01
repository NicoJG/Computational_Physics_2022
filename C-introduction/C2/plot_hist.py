# %%
import numpy as np
import matplotlib.pyplot as plt

x = np.genfromtxt("random_numbers.csv", delimiter=",")
print("Read numpy shape: ",x.shape)
plt.hist(x, bins=100, histtype="stepfilled")
plt.savefig("random_numbers_hist.pdf")
plt.show()
# %%
