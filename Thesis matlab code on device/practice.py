import numpy as np
import matplotlib.pyplot as plt

# Define Delta values for each phase
delta_topo = np.linspace(0, 145, 100)
Eg_topo = np.linspace(0.1, 0, 100)

delta_trivial = np.linspace(145, 250, 100)
Eg_trivial = np.linspace(0, 0.11, 100)

# Plot
plt.figure(figsize=(8, 6))

# Colors
color_topo = "#269e1b"  # Teal
color_trivial = "#d91f02"  # Gold

# Topological phase
plt.plot(delta_topo, Eg_topo, color=color_topo)
plt.fill_between(delta_topo, Eg_topo, color=color_topo, alpha=1.0,
                 label='Topological Insulator\nphase ($Z_2 = 1$)')

# Trivial phase
plt.plot(delta_trivial, Eg_trivial, color=color_trivial)
plt.fill_between(delta_trivial, Eg_trivial, color=color_trivial,
                 alpha=1.0, label='Trivial Insulator\nphase ($Z_2 = 0$)')

# Phase transition line
plt.axvline(x=0.3, color='black', linestyle='--', linewidth=1)

# Labels and formatting
plt.xlabel(r'$Electric field,E$ (mV/$\mathrm{\AA}$)', fontsize=20)
plt.ylabel(r'$Bandgap$ (eV)', fontsize=20)
plt.xlim(0, 250)
plt.ylim(0, 0.12)
plt.title('Topological to Trivial Phase Transition', fontsize=20)
# Increase legend font size
plt.legend(loc='upper right', fontsize=20)
plt.grid(True)
plt.tight_layout()
# Save the figure with high resolution (e.g., 300 dpi)
plt.savefig('D:/1.a Journal paper 2025/Figures/figure4_bandgap_E.png', dpi=300)
plt.show()
