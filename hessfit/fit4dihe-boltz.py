#!/usr/bin/env python3

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import argparse


def build_parser():
    """
    Build options for parser
    """
    parser = argparse.ArgumentParser()
    txt = "file containg phi angle, QM, and MM energies in 1st, 2nd, and 3rd columns"
    parser.add_argument('file', type=str, help=txt)
    # txt = 'Plot the PES with matplotlib'
    # parser.add_argument('--plot', action="store_true", help=txt)
    # txt = 'Name of output png file'
    # parser.add_argument('--name', default='scan', help=txt)
    return parser


def oplsa_matrix(xval_deg):
    xrad = np.deg2rad(xval_deg)
    A = np.column_stack([
        np.ones_like(xrad),                      # constant term
        1 + np.cos(xrad),                        # V1
        1 - np.cos(2 * xrad),                    # V2
        1 + np.cos(3 * xrad),                    # V3
        1 - np.cos(4 * xrad)                     # V4
    ])
    return A

def boltzmann_weights(E_qm, T):
    kB_kcal = 1.9872036e-3  # kcal/mol·K
    return np.exp(-E_qm / (kB_kcal * T))

def weighted_error(coeffs, A, E_qm, T):
    E_mm_model = A @ coeffs
    residual = E_mm_model - E_qm
    w = boltzmann_weights(E_qm, T)
    err = np.sqrt(np.sum((residual ** 2) * w) / len(E_qm))
    return err

def plot_fit(xval, qm, mm, yfit, T):
    fig = plt.subplots(1, figsize=(10, 8), dpi=96)
    plt.plot(xval, qm, label='QM', marker='o', markersize=8,
             fillstyle='none', c='black', linestyle='dashed', dashes=(5, 10))
    plt.plot(xval, mm, label='MM', marker='o', markersize=8,
             fillstyle='none', c='orangered', linestyle='dashed', dashes=(5, 10))
    plt.plot(xval, yfit, label=f'Fitted$_{{BW}}$ T={T}K', marker='o',
             markersize=10, c='blue', alpha=0.6)
    plt.xlabel('Angle (°)', fontsize=18)
    plt.ylabel('$\Delta$E (kcal/mol)', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(prop={'size': 20})
    plt.tight_layout()
    plt.savefig(f'fit_oplsa_T{T}.png')
    plt.show()

# === Load input data ===

# Load scan data (angle, EQM, EMM)
parser = build_parser()
args = parser.parse_args()
fname = args.file
data = pd.read_csv(fname, delim_whitespace=True, header=None)
angles = data.iloc[:, 0].to_numpy()
E_qm_raw = data.iloc[:, 1].to_numpy()
E_mm_raw = data.iloc[:, 2].to_numpy()

# Convert to relative energies in kcal/mol
E_qm = (E_qm_raw - np.min(E_qm_raw)) * 627.503
E_mm = (E_mm_raw - np.min(E_mm_raw)) * 627.503

# Build OPLS design matrix
A = oplsa_matrix(angles)

# Optimize using Boltzmann-weighted least squares at different temperatures
for T in [0, 500, 1000, 2000]:
    result = minimize(
        weighted_error,
        x0=np.zeros(5),        # 5 OPLS-AA parameters
        args=(A, E_qm, T),
        method='BFGS'
    )

    coeffs = result.x
    yfit = A @ coeffs

    print(f"\nFitted Coefficients at T = {T} K:")
    for i, c in enumerate(coeffs):
        print(f"  V{i} = {c:.4f} kcal/mol")
    print(f"Boltzmann-weighted RMSD: {result.fun:.4f} kcal/mol")

    # Plot results
    plot_fit(angles, E_qm, E_mm, yfit, T)
