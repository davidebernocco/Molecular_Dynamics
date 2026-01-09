import numpy as np
import matplotlib.pyplot as plt


# DIHEDRAL TIME EVOLUTION:
data = np.loadtxt("5_Data_Analysis/dih_time.xvg", comments=["@", "#"])  # Load XVG, skipping comments

time = data[:, 0]   # ps
angle = data[:, 1]  # degrees

plt.plot(time, angle)
plt.xlabel("Time (ps)")
plt.ylabel("Chi1 dihedral (degrees)")
plt.title("Cysteine χ₁ dihedral evolution")
plt.tight_layout()

plt.savefig("images/dihedral_time.png", dpi=300)
plt.close()



# DIHEDRAL HISTOGRAM:
data_histo = np.loadtxt("5_Data_Analysis/dih_histo.xvg", comments=["@", "#"])

angle_histo = data_histo[:, 0]    # degrees
prob_histo  = data_histo[:, 1]    # fraction of frames in each bin

plt.bar(angle_histo, prob_histo, width=1.0, edgecolor="black")
plt.xlabel("Dihedral angle (deg)")
plt.ylabel("Probability")
plt.tight_layout()
plt.savefig("images/dihedral_histogram.png", dpi=300)



# RADIAL DISTRIBUTION FUNCTION:
data_rdf = np.loadtxt("5_Data_Analysis/rdf_HG_WO.xvg", comments=["@", "#"])
r_rdf = data_rdf[:, 0]
g_rdf = data_rdf[:, 1]

plt.plot(r_rdf, g_rdf, lw=2)
plt.xlabel("Distance r (nm)")
plt.ylabel("g(r)")
plt.title("RDF: CYS HG - Water O")
plt.tight_layout()
plt.savefig("images/rdf_HG_OW.png", dpi=300)