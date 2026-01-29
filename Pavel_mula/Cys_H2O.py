import numpy as np
import matplotlib.pyplot as plt
from funct import load_data, dihedral_histogram
from funct import data_KDE, dihedral_KDE

"""
# DIHEDRAL HISTOGRAMS
angles_H2O_N = load_data("H2O", "N")
histo_H2O_N = dihedral_histogram(angles_H2O_N, 300, 3.6, "H2O", "N")

angles_H2O_O = load_data("H2O", "O")
angles_H2O_OXT = load_data("H2O", "OXT")
angles_H2O_O_OXT = np.hstack((angles_H2O_O, angles_H2O_OXT))
histo_H2O_O_OXT = dihedral_histogram(angles_H2O_O_OXT, 300, 3.6, "H2O", "O+OXT")


# DIHEDRAL KDE
angles_H2O_N_ext = data_KDE(angles_H2O_N)
kde_H2O_N = dihedral_KDE(angles_H2O_N_ext, 300, "H2O", "N", 0.02)

angles_H2O_O_OXT_ext = data_KDE(angles_H2O_O_OXT)
kde_H2O_O_OXT = dihedral_KDE(angles_H2O_O_OXT_ext, 300, "H2O", "O+OXT", 0.02)


# RADIAL DISTRIBUTION FUNCTION:
data_rdf_H = np.loadtxt("Data_Analysis_H2O/rdfH2O_HG_WH.xvg", comments=["@", "#"])
r_rdf_H = data_rdf_H[:, 0]
g_rdf_H = data_rdf_H[:, 1]

plt.plot(r_rdf_H, g_rdf_H, lw=2)
plt.xlabel("Distance r (nm)")
plt.ylabel("g(r)")
plt.title("(H2O) RDF: CYS HG - Water H")
plt.tight_layout()
plt.savefig("images_H2O/rdf_HG_WH.png", dpi=300)
plt.close()

data_rdf_S = np.loadtxt("Data_Analysis_H2O/rdfH2O_SG_WO.xvg", comments=["@", "#"])
r_rdf_S = data_rdf_S[:, 0]
g_rdf_S = data_rdf_S[:, 1]

plt.plot(r_rdf_S, g_rdf_S, lw=2)
plt.xlabel("Distance r (nm)")
plt.ylabel("g(r)")
plt.title("(H2O) RDF: CYS SG - Water O")
plt.tight_layout()
plt.savefig("images_H2O/rdf_SG_WO.png", dpi=300)
plt.close()
"""

# POTENTIAL ENERGY TIME EVOLUTION:
data_H = np.loadtxt("Data_Analysis_H2O/potential_H2O.xvg", comments=["@", "#"])
data_D = np.loadtxt("Data_Analysis_D2O/potential_D2O.xvg", comments=["@", "#"])

pot_en_H = data_H[:, 1] 
pot_en_D = data_D[:, 1] 


bin_width = abs((min(pot_en_H) - max(pot_en_H) ) / 100)
bins_H = np.arange(min(pot_en_H), max(pot_en_H) + bin_width, bin_width)
bin_centers_H = 0.5 * (bins_H[:-1] + bins_H[1:])
bins_D = np.arange(min(pot_en_D), max(pot_en_D) + bin_width, bin_width)
bin_centers_D = 0.5 * (bins_D[:-1] + bins_D[1:])

Nframes = len(pot_en_H)

# ============================
# HISTOGRAM PLOT
# ============================
counts_H, _ = np.histogram(pot_en_H, bins=bins_H)
prob_H = counts_H / Nframes
counts_D, _ = np.histogram(pot_en_D, bins=bins_D)
prob_D = counts_D / Nframes
plt.figure(figsize=(8, 5))
plt.bar(
    bin_centers_H,
    prob_H,
    width=bin_width,
    align="center",
    alpha=0.35,
    edgecolor="k",
    label="Pot. energy distribution (H2O)"
)
plt.bar(
    bin_centers_D,
    prob_D,
    width=bin_width,
    align="center",
    alpha=0.35,
    edgecolor="k",
    label="Pot. energy distribution (D2O)"
)
plt.xlabel("Potential energy (KJ/mol)")
plt.ylabel("Probability")
plt.xlim(max(pot_en_H), min(pot_en_D))
plt.legend()
plt.tight_layout()
plt.savefig("images_H2O/pot_energy.png", dpi=300)
plt.close()

# ============================
# FREE ENERGY PROFILE
# ============================
kb = 0.008314462618  # kJ/(mol*K)
temp = 300      # K
prob_safe_H = prob_H + 1e-12
prob_safe_D = prob_D + 1e-12
free_en_H = - kb * temp * np.log(prob_safe_H)
free_en_D = - kb * temp * np.log(prob_safe_D)

plt.figure(figsize=(8, 5))
plt.bar(
    bin_centers_H,
    free_en_H,
    width=bin_width,
    align="center",
    alpha=0.35,
    edgecolor="k",
    label=f"Gibbs free energy (H2O)"
)
plt.bar(
    bin_centers_D,
    free_en_D,
    width=bin_width,
    align="center",
    alpha=0.35,
    edgecolor="k",
    label=f"Gibbs free energy (D2O)"
)
plt.xlabel(f"Potential energy (kJ/mol)")
plt.ylabel("Free Energy (kJ/mol)")
plt.xlim(max(pot_en_H), min(pot_en_D))
plt.legend()
plt.tight_layout()
plt.savefig("images_H2O/free_energy.png", dpi=300)
plt.close()

