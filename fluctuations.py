import matplotlib.pyplot as plt
import numpy as np

# Read data from file
data = np.loadtxt('fluctuations.dat')

plt.rc('font', family='Times')
plt.rc('mathtext', fontset='cm')

# Extract columns from data
m = data[:, 0]
r = data[:, 1]

# Set up figure and subplots
fig = plt.figure(figsize=(6,4))
ax = fig.add_axes([0.2,0.1,0.75,0.8])

ax.tick_params(direction='in', top=True, right=True)

# Plot
ax.plot(m, r, linestyle='-', marker='', color='red')
ax.set_xlabel('molecule', fontfamily='Times')
ax.set_ylabel('$r$ $(\\AA)$', fontfamily='Times')
#ax.set_xlim(0, 16.00)
ax.set_ylim(2.9999,3.0001)

plt.margins(0, 0)

# Show plot
plt.savefig('fluctuations.pdf')
