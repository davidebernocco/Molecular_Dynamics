import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit




# Calculates correlations over the sampled points
def corr(n, N_max, lst):
    
    corr_arr = np.zeros(N_max, dtype = np.float32)
        
    for j in range(N_max):
        
        xi = 0
        x2i = 0
        xi_xij = 0
        
        for i in range(n - j):
            xi += lst[i]
            x2i += lst[i] ** 2 
            xi_xij += lst[i] * lst[i + j] 
            
        xi = xi / (n - j)
        x2i = x2i / (n - j)
        xi_xij = xi_xij / (n - j)
        
        corr_arr[j] = (xi_xij - xi ** 2) / (x2i - xi ** 2)
        
    return corr_arr


# Exponential function for the correlation fit
def exp_func(x, A, tau):
    return A * np.exp(-x / tau)



# Displays the correlation for the selected dihedral data
# and print the expoential fit parameters to 
# estimate the correlation length (=tau)
def corr_plot(s, Nf, j_max, data_arr):
    correlation = corr(Nf, j_max, data_arr)
    j_arr = np.arange(j_max)
    
    popt, pcov = curve_fit(exp_func, j_arr, correlation)

    A, tau = popt
    err = np.sqrt(np.diag(pcov))

    print(f"A_{s}   = {A:.4f} ± {err[0]:.4f}")
    print(f"tau_{s} = {tau:.4f} ± {err[1]:.4f}")

    xfit = np.linspace(j_arr.min(), j_arr.max(), 300)
    yfit = exp_func(xfit, *popt)

    plt.scatter(j_arr, correlation, label=f'C(j) for dihedral = {s} ', marker='o', s=50)
    plt.plot(xfit, yfit, label="Exponential fit", linewidth=2)
    plt.xlabel('j', fontsize=12)
    plt.ylabel(r'$ C(j) $', fontsize=12)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"images/dihedral_correlation_{s}.png", dpi=300)
    plt.close()

"""
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
    plt.savefig(f"images_bis/dihedral_histogram_{suf}.png", dpi=300)
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
    plt.savefig(f"images_bis/free_energy_{suf}.png", dpi=300)
    plt.close()

"""




