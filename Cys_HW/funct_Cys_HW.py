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
    plt.savefig(f"images/dihedral_HW_correlation_{s}.png", dpi=300)
    plt.close()











