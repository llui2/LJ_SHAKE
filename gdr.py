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

fig = plt.figure(figsize=(8, 3))

ax = fig.add_axes([0.1,0.15,0.8,0.75])
#ax1 = fig.add_axes([0.4,0.4,0.48,0.45])

ax.tick_params(direction='in', top=True, right=True)

# Plot
ax.set_title(f"Shake")
ax.plot(r, g, linestyle='-', marker='', color='red')
ax.set_xlabel('$r$ $(\\AA)$', fontfamily='Times')
ax.set_ylabel('$g(r)$', fontfamily='Times')
ax.set_xlim(0, 16.00)
ax.set_ylim(0,)

# ax1.plot(r, g, linestyle='-', marker='', color='red', linewidth=0.5)
# ax1.set_xlabel('$r$', fontfamily='Times')
# ax1.set_ylabel('$g(r)$', fontfamily='Times')
# ax1.set_xlim(0, 10.00)
# ax1.set_ylim(0,6)

plt.margins(0, 0)

# Show plot
plt.savefig('gdr.pdf')
