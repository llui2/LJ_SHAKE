import matplotlib.pyplot as plt
import numpy as np

# Read data from file
data = np.loadtxt('thermodynamics.dat')

plt.rc('font', family='Times')
plt.rc('mathtext', fontset='cm')


# Extract columns from data
step = data[:, 1]
E = data[:, 2]
T = data[:, 3]


# Set up figure and subplots
fig, axs = plt.subplots(2, 1, figsize=(10,4))

plt.subplots_adjust(wspace=0.35)
axs[0].tick_params(direction='in', top=True, right=True)
axs[1].tick_params(direction='in', top=True, right=True)

# Plot magnetization data with error bars
axs[0].plot(step, E, linestyle='-', marker='', color='blue')
#axs[0].set_xlabel('tiempo (ps)', fontfamily='Times')
axs[0].set_xlabel('')
axs[0].set_ylabel('Energia (kcal/mol)', fontfamily='Times')
axs[0].set_xlim(0, 20)
axs[0].tick_params(
    axis='x',          # changes apply to the x-axis
    which='both',      # both major and minor ticks are affected
    bottom=False,      # ticks along the bottom edge are off
    top=False,         # ticks along the top edge are off
    labelbottom=False) # labels along the bottom edge are off
axs[0].set_ylim(-6000000, 0)

# Plot energy data with y error bars
axs[1].plot(step, T, linestyle='-', marker='', color='red')
axs[1].set_xlabel('tiempo (ps)', fontfamily='Times')
axs[1].set_ylabel('Temperatura', fontfamily='Times')
axs[1].set_xlim(0, 20)
axs[1].set_ylim(0, 400)

plt.margins(0, 0)

# Show plot
plt.savefig('thermo.pdf')
