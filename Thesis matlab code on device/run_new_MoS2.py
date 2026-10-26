import tbmodels
import numpy as np
import matplotlib.pyplot as plt
import os
# use only *_hr.dat
model = tbmodels.Model.from_wannier_files(
    hr_file='D:/1.Wte2_HSE_Jan25_part3/MoS2_hse_hr.dat')

# use all files
model = tbmodels.Model.from_wannier_files(
    hr_file='D:/1.Wte2_HSE_Jan25_part3/MoS2_hse_hr.dat',
    wsvec_file='D:/1.Wte2_HSE_Jan25_part3/MoS2_wsvec.dat',
    xyz_file='D:/1.Wte2_HSE_Jan25_part3/MoS2_centres.xyz',
    win_file='D:/1.Wte2_HSE_Jan25_part3/MoS2.win'
)


def kpgen(kpath, nopot):
    kpathp = len(kpath)
    kpoints_cryst = []

    for k in range(3):
        kxf = []
        for j in range(kpathp - 1):
            if j == 0:
                if kpath[j, k] == kpath[j + 1, k]:
                    kx = np.full(nopot[j], kpath[j, k])
                else:
                    kx = np.linspace(kpath[j, k], kpath[j + 1, k], nopot[j])
                kxf.extend(kx)
            else:
                if kpath[j, k] == kpath[j + 1, k]:
                    kx = np.full(nopot[j], kpath[j, k])
                else:
                    kx = np.linspace(kpath[j, k], kpath[j + 1, k], nopot[j])
                kxf.extend(kx[1:])

        kpoints_cryst.append(kxf)

    kpoints_cryst = np.array(kpoints_cryst).T
    # kpoints_cart = 2 * np.pi * kpoints_cryst
    kpoints_cart = kpoints_cryst

    return kpoints_cart, kpoints_cryst


# Example usage:
kpath = np.array([[-0.5, 0.0, 0], [0.0, 0, 0], [0.5, 0.0, 0]])
nopot = [40, 40]

kpoints_cart, kpoints_cryst = kpgen(kpath, nopot)

print(kpoints_cart)
# print(model.eigenval(kpath))


# bands_tb = np.array(model.eigenval(kpoints_cart))


# Initialize variables
# Given values
# Ef = 3.3958
dMo1 = 17.738845
dMo2 = 17.903203
dS1 = 19.558748
dS2 = 19.139499
dS3 = 16.083334
dS4 = 16.486848


voltages = np.arange(0, 0.6, 0.1)
# voltages = [0, 0]
bandgap = np.zeros(len(voltages))  # Initialize bandgap array
ii = 0
for Voltage in voltages:
    # Calculate E
    E = -(2 * Voltage) / (dS3 - dS1)
    # Initialize the Unit_cell_H matrix
    Unit_cell_H = np.zeros((48, 48))

    # Populate the diagonal elements
    for i in range(44):
        if i < 12:
            Unit_cell_H[i, i] = E * (dS3 - dMo1) + Voltage
        elif i < 24:
            Unit_cell_H[i, i] = E * (dS3 - dMo2) + Voltage
        elif i < 30:
            Unit_cell_H[i, i] = E * (dS3 - dS1) + Voltage
        elif i < 36:
            Unit_cell_H[i, i] = E * (dS3 - dS2) + Voltage
        elif i < 42:
            Unit_cell_H[i, i] = E * (dS3 - dS3) + Voltage
        elif i < 48:
            Unit_cell_H[i, i] = E * (dS3 - dS4) + Voltage

    # Calculate eigenvalues for each k-point
    bands_tb_E = np.array([
        np.real(np.linalg.eigvalsh(model.hamilton(
            kpoints_cart[x, :]) + Unit_cell_H))
        for x in range(len(kpoints_cart))
    ])

    Ef = max(bands_tb_E[:, 26])  # Fermi energy

    # Separate spin-up and spin-down bands
    column = bands_tb_E.shape[1]
    Esup_new = np.zeros((bands_tb_E.shape[0], (column + 1) // 2))
    Esdown_new = np.zeros((bands_tb_E.shape[0], column // 2))

    xx, yy = 0, 0
    for i in range(column):
        if i % 2 == 0:
            Esup_new[:, xx] = bands_tb_E[:, i]
            xx += 1
        else:
            Esdown_new[:, yy] = bands_tb_E[:, i]
            yy += 1
    bandgap[ii] = (np.min(bands_tb_E[:, 28]) -
                   np.max(bands_tb_E[:, 27]))  # Bandgap calculation
    ii += 1
    # Plot the results
    # plt.figure()
    # plt.plot(Esup_new, color="C1")
    # plt.plot(Esdown_new, color="C0")
    # plt.text(2, Ef, f'E={E:.4f} V/Å', fontsize=12,
    #        bbox=dict(facecolor='black', alpha=0.5))
    # plt.ylim([Ef - 1.5, Ef + 1.5])
    # plt.title(f'Band Structure at Voltage = {Voltage:.4f} V')
    # plt.xlabel('k-point index')
    # plt.ylabel('Energy (eV)')
    # plt.legend()
    # plt.show()

# bandgap = (np.min(bands_tb_E[:, 28]) -
 #          np.max(bands_tb_E[:, 27]))  # Bandgap calculation
# print(bands_tb_E.shape)


# print(f"Bandgap: {bandgap}")


# Plot the data
# Blue line with line width of 1.5
plt.plot(voltages, bandgap, color='b', linewidth=1.5)
plt.xlabel('Δ (eV)', fontsize=14)
plt.ylabel('Eg (eV)', fontsize=14)
# Custom annotation
plt.text(
    0.3, 0.4, 'Topological Insulator\nphase ($\\mathbb{Z}_2 = 1$)', fontsize=12, color='black')
plt.text(
    2.1, 0.4, 'Trivial Insulator\nphase ($\\mathbb{Z}_2 = 0$)', fontsize=12, color='black')

plt.ylim(0, 0.12)
plt.xlim(0, 0.5)
plt.tight_layout()
plt.show()
