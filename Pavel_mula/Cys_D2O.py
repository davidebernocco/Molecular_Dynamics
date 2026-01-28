import numpy as np
import matplotlib.pyplot as plt
from funct import load_data, dihedral_histogram
from funct import data_KDE, dihedral_KDE


# DIHEDRAL HISTOGRAMS
angles_D2O_N = load_data("D2O", "N")
histo_D2O_N = dihedral_histogram(angles_D2O_N, 300, 3.6, "D2O", "N")

angles_D2O_O = load_data("D2O", "O")
angles_D2O_OXT = load_data("D2O", "OXT")
angles_D2O_O_OXT = np.hstack((angles_D2O_O, angles_D2O_OXT))
histo_D2O_O_OXT = dihedral_histogram(angles_D2O_O_OXT, 300, 3.6, "D2O", "O+OXT")

# DIHEDRAL KDE
angles_D2O_N_ext = data_KDE(angles_D2O_N)
kde_D2O_N = dihedral_KDE(angles_D2O_N_ext, 300, "D2O", "N", 0.02)

angles_D2O_O_OXT_ext = data_KDE(angles_D2O_O_OXT)
kde_D2O_O_OXT = dihedral_KDE(angles_D2O_O_OXT_ext, 300, "D2O", "O+OXT", 0.02)


"""
# RADIAL DISTRIBUTION FUNCTION:
data_rdf_H = np.loadtxt("Data_Analysis_D2O/rdfD2O_HG_WD.xvg", comments=["@", "#"])
r_rdf_H = data_rdf_H[:, 0]
g_rdf_H = data_rdf_H[:, 1]

plt.plot(r_rdf_H, g_rdf_H, lw=2)
plt.xlabel("Distance r (nm)")
plt.ylabel("g(r)")
plt.title("RDF: CYS HG - Water H")
plt.tight_layout()
plt.savefig("images_D2O/rdf_HG_WD.png", dpi=300)
plt.close()

data_rdf_S = np.loadtxt("Data_Analysis_D2O/rdfD2O_SG_WO.xvg", comments=["@", "#"])
r_rdf_S = data_rdf_S[:, 0]
g_rdf_S = data_rdf_S[:, 1]

plt.plot(r_rdf_S, g_rdf_S, lw=2)
plt.xlabel("Distance r (nm)")
plt.ylabel("g(r)")
plt.title("RDF: CYS SG - Water O")
plt.tight_layout()
plt.savefig("images_D2O/rdf_SG_WO.png", dpi=300)
plt.close()

"""