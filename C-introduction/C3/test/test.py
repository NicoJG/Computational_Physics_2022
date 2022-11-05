# %%
a1 = -27.0014581862
v1 = [-93.3060103996,
        -34.0071584722,
        38.1271408583,
        -15.5026635693,
        -58.7469722884]
norm1 = 122.5059859961

import numpy as np
v1 = np.array(v1)
print(f"|v| = {np.linalg.norm(v1)}")
print(f"|a*v| = {np.linalg.norm(a1*v1)}")
print(f"|a|*|v| = {np.abs(a1)*np.linalg.norm(v1)}")
# %%
