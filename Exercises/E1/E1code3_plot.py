#!/usr/bin/env python
###############################################################################
# E1code3_plot
###############################################################################

import numpy as np
import matplotlib.pyplot as plt

# skip_header skips the first
# row in data.csv
array = np.genfromtxt('powerspectrum.csv', delimiter=',', skip_header=1)

fig, ax = plt.subplots()
#ax.plot(array[:, 0], array[:, 1])
ax.plot(array[:, 1], array[:, 0])

ax.set_xlabel('frequency (arb.unit)')
ax.set_ylabel('power spectrum (arb.unit)')
ax.grid()

fig.savefig('powerspectrum.pdf')
