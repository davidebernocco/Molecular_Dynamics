import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import gaussian_kde




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



# ============================
# CORRELATION TIME ESTIMATION
# ============================
# ATTENTION: since the dihedral data are periodic, the autocorrelation
# should be computed using GRAM MATRICES (the actual algorithm has to be mod.)
# corr_plot(suf, Nframes, 150, angles)


# LOAD DIHEDRAL TIME SERIES
def load_data(w_type, s):
    data = np.loadtxt(f"Data_Analysis_{w_type}/dih_time_{w_type}_{s}.xvg", comments=["@", "#"])
    angles_180 = data[:, 1]
    return np.mod(angles_180, 360)


# PLOT DIHEDRAL HISTOGRAM AND FREE ENERGY
def dihedral_histogram(angles, temp, bin_w, water_type, suf):

    # temp: temperature of the system in K
    kb = 0.0019872  # kcal/mol/K

    bin_width = bin_w
    bins = np.arange(0, 360 + bin_width, bin_width)
    bin_centers = 0.5 * (bins[:-1] + bins[1:])

    Nframes = len(angles)

    # ============================
    # HISTOGRAM PLOT
    # ============================
    counts, _ = np.histogram(angles, bins=bins)
    prob = counts / Nframes
    plt.figure(figsize=(8, 5))
    plt.bar(
        bin_centers,
        prob,
        width=bin_width,
        align="center",
        alpha=0.7,
        edgecolor="k",
        label=f"Dihedral distribution {water_type}"
    )
    plt.xlabel(f"Dihedral χ {suf} angle (degrees)")
    plt.ylabel("Probability")
    plt.xlim(0, 360)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"images_{water_type}/dihedral_histogram_{water_type}_{suf}.png", dpi=300)
    plt.close()

    # ============================
    # FREE ENERGY PROFILE
    # ============================
    prob_safe = prob + 1e-12
    free_en = - kb * temp * np.log(prob_safe)

    plt.figure(figsize=(8, 5))
    plt.bar(
        bin_centers,
        free_en,
        width=bin_width,
        align="center",
        alpha=0.7,
        edgecolor="k",
        label=f"Gibbs free energy {water_type}"
    )
    plt.xlabel(f"Dihedral χ {suf} angle (degrees)")
    plt.ylabel("Free Energy (kcal/mol)")
    plt.xlim(0, 360)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f"images_{water_type}/free_energy_{water_type}_{suf}.png", dpi=300)
    plt.close()



# LOAD DIHEDRAL TIME SERIES FOR KDE TREATMENT
def data_KDE(ang_360):
    angles_ext = np.concatenate([ang_360 - 360.0, ang_360, ang_360 + 360.0])
    return angles_ext



# RETURN DIHEDRAL KDE DISTRIBUTION AND FREE ENERGY
def dihedral_KDE(ang, temp, water_type, suf, kde_band):
    kb = 0.0019872  # kcal/mol/K
    kde = gaussian_kde(ang, bw_method=kde_band)
    phi = np.linspace(0, 360, 720)   # 0.5° resolution
    P = kde(phi)
    P /= np.trapz(P, phi)

    plt.figure(figsize=(7,5))
    plt.plot(phi, P, lw=2)
    plt.xlabel(f"Dihedral χ {suf} angle (deg)")
    plt.ylabel(f"KDE dihedral distribution {water_type} (kJ/mol)")
    plt.xlim(0, 360)
    plt.tight_layout()
    plt.savefig(f"images_{water_type}/KDE_dihedral_{water_type}_{suf}.png", dpi=300)
    plt.close()


    F = -kb * temp * np.log(P)
    F -= F.min()

    plt.figure(figsize=(7,5))
    plt.plot(phi, F, lw=2)
    plt.xlabel(f"Dihedral χ {suf} angle (deg)")
    plt.ylabel(f"Free energy {water_type} (kJ/mol)")
    plt.xlim(0, 360)
    plt.tight_layout()
    plt.savefig(f"images_{water_type}/free_energyKDE_{water_type}_{suf}.png", dpi=300)
    plt.close()
