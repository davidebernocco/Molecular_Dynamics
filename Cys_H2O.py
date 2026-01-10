import numpy as np
import matplotlib.pyplot as plt

"""
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
"""


# DIHEDRAL HISTOGRAM:
data_histo = np.loadtxt("5_Data_Analysis/dih_histo.xvg", comments=["@", "#"])

angle_histo = data_histo[:, 0]    # degrees
prob_histo  = data_histo[:, 1]    # fraction of frames in each bin

bin_width = 10.0
bins = np.arange(-180, 180 + bin_width, bin_width)
bin_centers = 0.5 * (bins[:-1] + bins[1:])

hist, _ = np.histogram(angle_histo, bins=bins, weights=prob_histo)

N = 125
counts = hist * N
count_err = np.sqrt(counts)
prob_err = count_err/N


"""
plt.figure(figsize=(8, 5))

plt.bar(
    bin_centers,
    hist,
    width=bin_width,
    align='center',
    alpha=0.7,
    edgecolor='k',
    label="Dihedral distribution"
)

plt.errorbar(
    bin_centers,
    hist,
    yerr=prob_err,
    fmt='o',
    color='black',
    capsize=3,
    label="Poisson error (N=125)"
)

plt.xlabel("Dihedral angle (degrees)")
plt.ylabel("Probability")
plt.xlim(-180, 180)
plt.legend()
plt.tight_layout()
plt.savefig("images/dihedral_histogram.png", dpi=300)
"""


k_b = 1.380649e-23
temp = 300
free_en = - k_b * temp * np.log(hist+1e-10)
free_en_err = k_b * temp * (1/(hist+1e-10)) * prob_err

plt.figure(figsize=(8, 5))

plt.bar(
    bin_centers,
    free_en,
    width=bin_width,
    align='center',
    alpha=0.7,
    edgecolor='k',
    label="Gibbs free energy"
)

plt.errorbar(
    bin_centers,
    free_en,
    yerr=free_en_err,
    fmt='o',
    color='black',
    capsize=3,
    label="Poisson error (N=125)"
)

plt.xlabel("Dihedral angle (degrees)")
plt.ylabel("Free Energy (Joule)")
plt.xlim(-180, 180)
plt.legend()
plt.tight_layout()
plt.savefig("images/free_energy.png", dpi=300)



"""
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

"""