import numpy as np
import matplotlib.pyplot as plt
from funct_Cys_H2O import corr_plot



suffixes = ["N", "O", "OXT"]

kb = 0.0019872  # kcal/mol/K
temp = 300      # K

bin_width = 10.0
bins = np.arange(-180, 180 + bin_width, bin_width)
bin_centers = 0.5 * (bins[:-1] + bins[1:])

Nblocks = 10   # adjust if you have longer trajectories

for suf in suffixes:

    # ============================
    # LOAD TIME SERIES
    # ============================
    data = np.loadtxt(f"5_Data_Analysis/dih_time_{suf}.xvg", comments=["@", "#"])
    angles = data[:, 1]

    Nframes = len(angles)


    # ============================
    # CORRELATION TIME ESTIMATION
    # ============================
    corr_plot(suf, Nframes, 150, angles)


    # ============================
    # SPLIT INTO BLOCKS
    # ============================
    block_size = Nframes // Nblocks
    angles = angles[:block_size * Nblocks]  # trim remainder

    block_probs = []

    for i in range(Nblocks):
        block = angles[i * block_size : (i + 1) * block_size]
        counts, _ = np.histogram(block, bins=bins)
        prob = counts / len(block)
        block_probs.append(prob)

    block_probs = np.array(block_probs)

    # ============================
    # MEAN & ERROR
    # ============================
    prob_mean = block_probs.mean(axis=0)
    prob_std = block_probs.std(axis=0, ddof=1)   # sample std
    prob_err = prob_std / np.sqrt(Nblocks)


    # ============================
    # HISTOGRAM PLOT
    # ============================
    plt.figure(figsize=(8, 5))

    plt.bar(
        bin_centers,
        prob_mean,
        width=bin_width,
        align="center",
        alpha=0.7,
        edgecolor="k",
        label="Dihedral distribution"
    )

    plt.errorbar(
        bin_centers,
        prob_mean,
        yerr=prob_err,
        fmt="o",
        color="black",
        capsize=3,
        label=f"Block-avg error ({Nblocks} blocks)"
    )

    plt.xlabel(f"Dihedral χ_{suf} angle (degrees)")
    plt.ylabel("Probability")
    plt.xlim(-180, 180)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"images/dihedral_histogram_{suf}.png", dpi=300)
    plt.close()


    # ============================
    # FREE ENERGY PROFILE
    # ============================
    prob_safe = prob_mean + 1e-12
    free_en = - kb * temp * np.log(prob_safe)
    free_en_err = kb * temp * (prob_err / prob_safe)

    plt.figure(figsize=(8, 5))

    plt.bar(
        bin_centers,
        free_en,
        width=bin_width,
        align="center",
        alpha=0.7,
        edgecolor="k",
        label="Gibbs free energy"
    )

    plt.errorbar(
        bin_centers,
        free_en,
        yerr=free_en_err,
        fmt="o",
        color="black",
        capsize=3,
        label="Block-avg error"
    )

    plt.xlabel(f"Dihedral χ_{suf} angle (degrees)")
    plt.ylabel("Free Energy (kcal/mol)")
    plt.xlim(-180, 180)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"images/free_energy_{suf}.png", dpi=300)
    plt.close()





# RADIAL DISTRIBUTION FUNCTION:
data_rdf_H = np.loadtxt("5_Data_Analysis/rdf_HG_WH.xvg", comments=["@", "#"])
r_rdf_H = data_rdf_H[:, 0]
g_rdf_H = data_rdf_H[:, 1]

plt.plot(r_rdf_H, g_rdf_H, lw=2)
plt.xlabel("Distance r (nm)")
plt.ylabel("g(r)")
plt.title("RDF: CYS HG - Water H")
plt.tight_layout()
plt.savefig("images/rdf_HG_WH.png", dpi=300)
plt.close()

data_rdf_S = np.loadtxt("5_Data_Analysis/rdf_SG_WO.xvg", comments=["@", "#"])
r_rdf_S = data_rdf_S[:, 0]
g_rdf_S = data_rdf_S[:, 1]

plt.plot(r_rdf_S, g_rdf_S, lw=2)
plt.xlabel("Distance r (nm)")
plt.ylabel("g(r)")
plt.title("RDF: CYS SG - Water O")
plt.tight_layout()
plt.savefig("images/rdf_SG_WO.png", dpi=300)
plt.close()

