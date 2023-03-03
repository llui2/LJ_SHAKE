import matplotlib.pyplot as plt
import numpy as np

# Read data from file
data = np.loadtxt('gdr.dat')

plt.rc('font', family='Times')
plt.rc('mathtext', fontset='cm')

# Extract columns from data
r = data[:, 0]
g = data[:, 1]

# Set up figure and subplots
fig, ax = plt.subplots(1, 1, figsize=(8, 3))

ax.tick_params(direction='in', top=True, right=True)

# Plot
plt.title(f"No Shake")
ax.plot(r, g, linestyle='-', marker='', color='red')
ax.set_xlabel('$r$', fontfamily='Times')
ax.set_ylabel('$g(r)$', fontfamily='Times')
ax.set_xlim(0, 16.00)
ax.set_ylim(0, 6)

plt.margins(0, 0)

# Show plot
plt.savefig('gdr.pdf')
